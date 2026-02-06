//! FFT Tree for FALCON.
//!
//! This module implements the FFT tree structure used in FALCON's
//! Fast Fourier Sampling (ffSampling) algorithm. The tree stores the
//! Gram-Schmidt orthogonalized secret key basis in FFT form.
//!
//! The tree has log2(n) + 1 levels, where each level stores the
//! Gram-Schmidt data for that recursion level.

use crate::fft::{split_fft, Complex};

/// A node in the FFT tree storing LDL* decomposition data.
#[derive(Clone, Debug)]
pub struct FftNode {
    /// The sigma values: sqrt of the LDL* diagonal at this level.
    pub sigma: Vec<f64>,
    /// The l10 values: off-diagonal factor from LDL* decomposition.
    /// At the base case (n=1), l10 has one element.
    pub l10: Vec<Complex>,
}

impl FftNode {
    /// Creates an empty node.
    pub fn empty() -> Self {
        FftNode {
            sigma: Vec::new(),
            l10: Vec::new(),
        }
    }

    /// Creates a node with the given data.
    pub fn new(sigma: Vec<f64>, l10: Vec<Complex>) -> Self {
        FftNode { sigma, l10 }
    }
}

/// The FFT tree structure for FALCON's sampler.
///
/// The tree stores the Gram-Schmidt orthogonalization of the secret key
/// in a form suitable for the recursive FFT sampling algorithm.
#[derive(Clone, Debug)]
pub struct FftTree {
    /// The polynomial degree n.
    pub n: usize,
    /// Tree depth (log2(n) + 1).
    pub depth: usize,
    /// Tree nodes, organized by level.
    /// Level 0 has 1 node (root), level k has 2^k nodes.
    nodes: Vec<FftNode>,
}

impl FftTree {
    /// Creates a new FFT tree for the given polynomial degree.
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two(), "n must be a power of 2");
        let depth = n.trailing_zeros() as usize + 1;

        // Total nodes: 1 + 2 + 4 + ... + n = 2n - 1
        let total_nodes = 2 * n - 1;
        let nodes = vec![FftNode::empty(); total_nodes];

        FftTree { n, depth, nodes }
    }

    /// Returns the index of a node at a given level and position.
    fn node_index(&self, level: usize, pos: usize) -> usize {
        // Level 0 starts at index 0
        // Level k starts at index 2^k - 1
        (1 << level) - 1 + pos
    }

    /// Gets a reference to a node.
    pub fn get_node(&self, level: usize, pos: usize) -> &FftNode {
        &self.nodes[self.node_index(level, pos)]
    }

    /// Gets a mutable reference to a node.
    pub fn get_node_mut(&mut self, level: usize, pos: usize) -> &mut FftNode {
        let idx = self.node_index(level, pos);
        &mut self.nodes[idx]
    }

    /// Sets a node's data.
    pub fn set_node(&mut self, level: usize, pos: usize, node: FftNode) {
        let idx = self.node_index(level, pos);
        self.nodes[idx] = node;
    }

    /// Returns the number of nodes at a given level.
    pub fn nodes_at_level(&self, level: usize) -> usize {
        1 << level
    }

    /// Builds the FFT tree from the secret key polynomials using ffLDL*.
    ///
    /// Given (f, g, F, G) in FFT form, computes the initial Gram matrix and
    /// builds the LDL* tree recursively per FIPS 206.
    pub fn build_from_basis(
        f_fft: &[Complex],
        g_fft: &[Complex],
        big_f_fft: &[Complex],
        big_g_fft: &[Complex],
    ) -> Self {
        let n = f_fft.len();
        let mut tree = FftTree::new(n);

        // Compute the initial Gram matrix G = B * adj(B)
        // where B = [[f, g], [F, G]] in FFT form.
        // g00[i] = |f[i]|^2 + |g[i]|^2  (real)
        // g01[i] = f[i]*conj(F[i]) + g[i]*conj(G[i])  (complex)
        // g11[i] = |F[i]|^2 + |G[i]|^2  (real)
        let g00: Vec<f64> = (0..n)
            .map(|i| f_fft[i].norm_sq() + g_fft[i].norm_sq())
            .collect();
        let g01: Vec<Complex> = (0..n)
            .map(|i| f_fft[i] * big_f_fft[i].conj() + g_fft[i] * big_g_fft[i].conj())
            .collect();
        let g11: Vec<f64> = (0..n)
            .map(|i| big_f_fft[i].norm_sq() + big_g_fft[i].norm_sq())
            .collect();

        tree.ffldl_recursive(0, 0, &g00, &g01, &g11);

        tree
    }

    /// Recursively builds the ffLDL* tree from a 2x2 Hermitian Gram matrix.
    ///
    /// At each level, computes the LDL* decomposition:
    ///   G = [[g00, g01], [conj(g01), g11]] = L * D * L*
    /// where L = [[1, 0], [l10, 1]], D = diag(d0, d1)
    ///
    /// Then stores sigma = sqrt(d0) and l10 at this node, and recurses
    /// on the children using the FFT-split of d0 and d1.
    fn ffldl_recursive(
        &mut self,
        level: usize,
        pos: usize,
        g00: &[f64],
        g01: &[Complex],
        g11: &[f64],
    ) {
        let n = g00.len();

        if n == 1 {
            // Base case: scalar LDL* decomposition
            let d0 = g00[0];
            let l10 = g01[0].conj().scale(1.0 / d0);
            let d1 = g11[0] - g01[0].norm_sq() / d0;

            let sigma = vec![d0.sqrt(), d1.max(0.0).sqrt()];
            self.set_node(level, pos, FftNode::new(sigma, vec![l10]));
            return;
        }

        // Pointwise LDL* decomposition
        let mut d0 = vec![0.0f64; n];
        let mut l10 = vec![Complex::ZERO; n];
        let mut d1 = vec![0.0f64; n];

        for i in 0..n {
            d0[i] = g00[i];
            l10[i] = g01[i].conj().scale(1.0 / g00[i]);
            d1[i] = g11[i] - g01[i].norm_sq() / g00[i];
            // Clamp d1 to avoid negative values from floating-point error
            if d1[i] < 0.0 {
                d1[i] = 0.0;
            }
        }

        // Store sigma = sqrt(d0) and l10 at this node
        let sigma: Vec<f64> = d0.iter().map(|&x| x.sqrt()).collect();
        self.set_node(level, pos, FftNode::new(sigma, l10));

        // Split d0 and d1 into even/odd halves via FFT splitting
        let d0_complex: Vec<Complex> = d0.iter().map(|&x| Complex::from_real(x)).collect();
        let d1_complex: Vec<Complex> = d1.iter().map(|&x| Complex::from_real(x)).collect();

        let (d0_even, d0_odd) = split_fft(&d0_complex);
        let (d1_even, d1_odd) = split_fft(&d1_complex);

        // Left child: Gram matrix from split of d0
        // [[Re(d0_even), d0_odd], [conj(d0_odd), Re(d0_even)]]
        let g00_left: Vec<f64> = d0_even.iter().map(|c| c.re).collect();
        let g11_left = g00_left.clone();

        // Right child: Gram matrix from split of d1
        let g00_right: Vec<f64> = d1_even.iter().map(|c| c.re).collect();
        let g11_right = g00_right.clone();

        self.ffldl_recursive(level + 1, 2 * pos, &g00_left, &d0_odd, &g11_left);
        self.ffldl_recursive(level + 1, 2 * pos + 1, &g00_right, &d1_odd, &g11_right);
    }

    /// Gets the sigma value at a given level and position.
    pub fn get_sigma(&self, level: usize, pos: usize, index: usize) -> f64 {
        let node = self.get_node(level, pos);
        node.sigma.get(index).copied().unwrap_or(1.0)
    }
}

/// Gram-Schmidt data for the full basis.
///
/// Stores the Gram-Schmidt orthogonalization of [[f, g], [F, G]]
/// in a form suitable for sampling.
#[derive(Clone, Debug)]
pub struct GramSchmidt {
    /// The [[f, g], [F, G]] basis in FFT form.
    pub f_fft: Vec<Complex>,
    pub g_fft: Vec<Complex>,
    pub big_f_fft: Vec<Complex>,
    pub big_g_fft: Vec<Complex>,

    /// The FFT tree for sampling.
    pub tree: FftTree,

    /// Precomputed sigma values for the full basis.
    pub sigma_fg: Vec<f64>,
    pub sigma_FG: Vec<f64>,
}

impl GramSchmidt {
    /// Creates a new Gram-Schmidt structure from the basis polynomials.
    pub fn new(
        f_fft: Vec<Complex>,
        g_fft: Vec<Complex>,
        big_f_fft: Vec<Complex>,
        big_g_fft: Vec<Complex>,
    ) -> Self {
        // Build the FFT tree
        let tree = FftTree::build_from_basis(&f_fft, &g_fft, &big_f_fft, &big_g_fft);

        // Compute sigma values
        let n = f_fft.len();
        let sigma_fg: Vec<f64> = (0..n)
            .map(|i| (f_fft[i].norm_sq() + g_fft[i].norm_sq()).sqrt())
            .collect();

        let sigma_FG: Vec<f64> = (0..n)
            .map(|i| (big_f_fft[i].norm_sq() + big_g_fft[i].norm_sq()).sqrt())
            .collect();

        GramSchmidt {
            f_fft,
            g_fft,
            big_f_fft,
            big_g_fft,
            tree,
            sigma_fg,
            sigma_FG,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fft::fft;

    #[test]
    fn test_fft_tree_structure() {
        let tree = FftTree::new(8);
        assert_eq!(tree.n, 8);
        assert_eq!(tree.depth, 4); // log2(8) + 1 = 4
        assert_eq!(tree.nodes.len(), 15); // 2*8 - 1 = 15
    }

    #[test]
    fn test_node_indexing() {
        let tree = FftTree::new(8);

        // Level 0: 1 node at index 0
        assert_eq!(tree.node_index(0, 0), 0);

        // Level 1: 2 nodes at indices 1, 2
        assert_eq!(tree.node_index(1, 0), 1);
        assert_eq!(tree.node_index(1, 1), 2);

        // Level 2: 4 nodes at indices 3, 4, 5, 6
        assert_eq!(tree.node_index(2, 0), 3);
        assert_eq!(tree.node_index(2, 3), 6);

        // Level 3: 8 nodes at indices 7-14
        assert_eq!(tree.node_index(3, 0), 7);
        assert_eq!(tree.node_index(3, 7), 14);
    }

    #[test]
    fn test_build_from_basis() {
        // Create simple test polynomials in FFT form
        let n = 4;
        let mut f: Vec<Complex> = vec![1.0, 0.0, 0.0, 0.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();
        let mut g: Vec<Complex> = vec![0.0, 1.0, 0.0, 0.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();
        let mut big_f: Vec<Complex> = vec![0.0, 0.0, 1.0, 0.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();
        let mut big_g: Vec<Complex> = vec![0.0, 0.0, 0.0, 1.0]
            .into_iter()
            .map(|x| Complex::from_real(x))
            .collect();

        fft(&mut f);
        fft(&mut g);
        fft(&mut big_f);
        fft(&mut big_g);

        let tree = FftTree::build_from_basis(&f, &g, &big_f, &big_g);
        assert_eq!(tree.n, 4);
        assert_eq!(tree.depth, 3);
    }

    #[test]
    fn test_gram_schmidt_creation() {
        let n = 4;
        let f_fft: Vec<Complex> = (0..n).map(|i| Complex::new(1.0 + i as f64, 0.0)).collect();
        let g_fft: Vec<Complex> = (0..n).map(|i| Complex::new(0.0, 1.0 + i as f64)).collect();
        let big_f_fft: Vec<Complex> = (0..n).map(|i| Complex::new(2.0 * i as f64, 0.0)).collect();
        let big_g_fft: Vec<Complex> = (0..n).map(|i| Complex::new(0.0, 2.0 * i as f64)).collect();

        let gs = GramSchmidt::new(f_fft, g_fft, big_f_fft, big_g_fft);

        assert_eq!(gs.sigma_fg.len(), n);
        assert_eq!(gs.sigma_FG.len(), n);
        assert!(gs.sigma_fg.iter().all(|&s| s > 0.0));
    }
}
