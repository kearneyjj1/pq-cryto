//! FFT Tree for FALCON.
//!
//! This module implements the FFT tree structure used in FALCON's
//! Fast Fourier Sampling (ffSampling) algorithm. The tree stores the
//! Gram-Schmidt orthogonalized secret key basis in FFT form.
//!
//! The tree has log2(n) + 1 levels, where each level stores the
//! Gram-Schmidt data for that recursion level.

use crate::fft::{merge_fft, split_fft, Complex};

/// A node in the FFT tree storing Gram-Schmidt data.
#[derive(Clone, Debug)]
pub struct FftNode {
    /// The sigma values (standard deviations) for this node.
    pub sigma: Vec<f64>,
    /// The Gram-Schmidt coefficients for this node.
    pub gs_coeffs: Vec<Complex>,
}

impl FftNode {
    /// Creates an empty node.
    pub fn empty() -> Self {
        FftNode {
            sigma: Vec::new(),
            gs_coeffs: Vec::new(),
        }
    }

    /// Creates a node with the given data.
    pub fn new(sigma: Vec<f64>, gs_coeffs: Vec<Complex>) -> Self {
        FftNode { sigma, gs_coeffs }
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

    /// Builds the FFT tree from the secret key polynomials.
    ///
    /// Given (f, g, F, G) in FFT form, builds the tree for sampling.
    pub fn build_from_basis(
        f_fft: &[Complex],
        g_fft: &[Complex],
        big_f_fft: &[Complex],
        big_g_fft: &[Complex],
    ) -> Self {
        let n = f_fft.len();
        let mut tree = FftTree::new(n);

        // Build recursively
        tree.build_level(0, 0, f_fft, g_fft, big_f_fft, big_g_fft);

        tree
    }

    /// Recursively builds a level of the tree.
    fn build_level(
        &mut self,
        level: usize,
        pos: usize,
        f_fft: &[Complex],
        g_fft: &[Complex],
        big_f_fft: &[Complex],
        big_g_fft: &[Complex],
    ) {
        let n = f_fft.len();

        if n == 1 {
            // Base case: compute the Gram-Schmidt data for a 2x2 matrix
            // [[f, g], [F, G]] forms the basis
            // Gram-Schmidt: b1 = (f, g), b2 = (F, G) - proj_{b1}(b2) * b1

            let f = f_fft[0];
            let g = g_fft[0];
            let big_f = big_f_fft[0];
            let big_g = big_g_fft[0];

            // ||(f, g)||^2
            let norm_b1_sq = f.norm_sq() + g.norm_sq();

            // <(F, G), (f, g)>
            let inner = big_f * f.conj() + big_g * g.conj();

            // mu = inner / norm_b1_sq
            let mu = inner.scale(1.0 / norm_b1_sq);

            // b2* = (F, G) - mu * (f, g)
            let b2_star_f = big_f - mu * f;
            let b2_star_g = big_g - mu * g;
            let norm_b2_star_sq = b2_star_f.norm_sq() + b2_star_g.norm_sq();

            // Sigma values: sqrt of the squared norms
            let sigma = vec![norm_b1_sq.sqrt(), norm_b2_star_sq.sqrt()];
            let gs_coeffs = vec![mu];

            self.set_node(level, pos, FftNode::new(sigma, gs_coeffs));
            return;
        }

        // Recursive case: split into two half-size problems
        let (f0, f1) = split_fft(f_fft);
        let (g0, g1) = split_fft(g_fft);
        let (big_f0, big_f1) = split_fft(big_f_fft);
        let (big_g0, big_g1) = split_fft(big_g_fft);

        // Compute the Gram-Schmidt data for this level
        // The Gram matrix at this level involves the field norms

        // For the FFT tree, we store:
        // - sigma: the Gram-Schmidt sigma values
        // - gs_coeffs: the projection coefficients

        // Compute Gram matrix entries
        // G = [[<(f,g), (f,g)>, <(f,g), (F,G)>],
        //      [<(F,G), (f,g)>, <(F,G), (F,G)>]]

        let mut norm_fg_sq = vec![Complex::ZERO; n];
        let mut inner_fgFG = vec![Complex::ZERO; n];
        let mut norm_FG_sq = vec![Complex::ZERO; n];

        for i in 0..n {
            norm_fg_sq[i] = Complex::from_real(f_fft[i].norm_sq() + g_fft[i].norm_sq());
            inner_fgFG[i] = f_fft[i] * big_f_fft[i].conj() + g_fft[i] * big_g_fft[i].conj();
            norm_FG_sq[i] = Complex::from_real(big_f_fft[i].norm_sq() + big_g_fft[i].norm_sq());
        }

        // Compute sigma from the diagonal of the Gram matrix
        let sigma: Vec<f64> = norm_fg_sq.iter().map(|c| c.re.sqrt()).collect();

        // Store this level's data
        self.set_node(level, pos, FftNode::new(sigma, inner_fgFG.clone()));

        // Recurse on children
        // Left child uses f0, g0
        // Right child uses f1, g1
        self.build_level(level + 1, 2 * pos, &f0, &g0, &big_f0, &big_g0);
        self.build_level(level + 1, 2 * pos + 1, &f1, &g1, &big_f1, &big_g1);
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
