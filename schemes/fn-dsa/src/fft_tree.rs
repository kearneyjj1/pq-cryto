//! LDL* Tree for FALCON (FIPS 206).
//!
//! This module implements the LDL* tree structure used in FALCON's
//! Fast Fourier Sampling (ffSampling) algorithm. The tree stores the
//! LDL* decomposition of the Gram matrix G = B · B^* where B is the
//! secret key basis.
//!
//! # FIPS 206 Compliance
//!
//! This implementation follows the FALCON specification for the ffLDL
//! algorithm, which recursively computes the LDL* decomposition in the
//! FFT domain. The tree has log2(n) + 1 levels, where each level stores
//! the L factors (off-diagonal) and sigma values (sqrt of D diagonal).
//!
//! # Algorithm Overview
//!
//! For a 2×2 Gram matrix G = [[G00, G01], [G10, G11]]:
//! - d0 = G00 (first diagonal of D)
//! - l10 = G10 / G00 (off-diagonal of L)
//! - d1 = G11 - |l10|² · G00 (second diagonal of D)
//!
//! The recursion splits each FFT representation in half, creating
//! two independent sub-problems at each level.

use crate::fft::{merge_fft, split_fft, Complex};
use crate::gaussian::{SIGMA_MIN_512, SIGMA_MIN_1024};

// ============================================================================
// LDL* Tree Node
// ============================================================================

/// A node in the LDL* tree storing decomposition data.
///
/// Each node stores the L factor (off-diagonal element) and the sigma
/// values (square roots of the D diagonal elements) for its level.
#[derive(Clone, Debug)]
pub struct LdlNode {
    /// The sigma values (standard deviations) for this node.
    /// These are sqrt(d0) and sqrt(d1) from the LDL* decomposition.
    /// For leaf nodes: [sigma0, sigma1] for the 2×2 base case.
    /// For internal nodes: sigma values for the polynomial coefficients.
    pub sigma: Vec<f64>,

    /// The L factor (off-diagonal element l10 = G10/G00).
    /// This is a complex polynomial in FFT form.
    pub l10: Vec<Complex>,
}

impl LdlNode {
    /// Creates an empty node.
    pub fn empty() -> Self {
        LdlNode {
            sigma: Vec::new(),
            l10: Vec::new(),
        }
    }

    /// Creates a node with the given data.
    pub fn new(sigma: Vec<f64>, l10: Vec<Complex>) -> Self {
        LdlNode { sigma, l10 }
    }

    /// Creates a leaf node for the 2×2 base case.
    pub fn leaf(sigma0: f64, sigma1: f64, l10: Complex) -> Self {
        LdlNode {
            sigma: vec![sigma0, sigma1],
            l10: vec![l10],
        }
    }
}

// ============================================================================
// LDL* Tree Structure
// ============================================================================

/// The LDL* tree structure for FALCON's sampler.
///
/// The tree stores the LDL* decomposition of the Gram matrix G = B · B^*
/// in a form suitable for the recursive ffSampling algorithm.
///
/// Tree organization:
/// - Level 0 (root): 1 node covering all n FFT positions
/// - Level k: 2^k nodes, each covering n/2^k positions
/// - Level log2(n): n leaf nodes, each a 2×2 base case
#[derive(Clone, Debug)]
pub struct LdlTree {
    /// The polynomial degree n.
    pub n: usize,
    /// Tree depth (log2(n) + 1).
    pub depth: usize,
    /// Tree nodes, organized by level.
    /// Level k starts at index 2^k - 1.
    nodes: Vec<LdlNode>,
    /// Minimum sigma value (for clamping during sampling).
    pub sigma_min: f64,
}

impl LdlTree {
    /// Creates a new empty LDL* tree for the given polynomial degree.
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two(), "n must be a power of 2");
        let depth = n.trailing_zeros() as usize + 1;

        // Total nodes: 1 + 2 + 4 + ... + n = 2n - 1
        let total_nodes = 2 * n - 1;
        let nodes = vec![LdlNode::empty(); total_nodes];

        // Select sigma_min based on n
        let sigma_min = if n <= 512 {
            SIGMA_MIN_512
        } else {
            SIGMA_MIN_1024
        };

        LdlTree {
            n,
            depth,
            nodes,
            sigma_min,
        }
    }

    /// Returns the index of a node at a given level and position.
    #[inline]
    fn node_index(&self, level: usize, pos: usize) -> usize {
        // Level k starts at index 2^k - 1
        (1 << level) - 1 + pos
    }

    /// Gets a reference to a node.
    #[inline]
    pub fn get_node(&self, level: usize, pos: usize) -> &LdlNode {
        &self.nodes[self.node_index(level, pos)]
    }

    /// Gets a mutable reference to a node.
    #[inline]
    pub fn get_node_mut(&mut self, level: usize, pos: usize) -> &mut LdlNode {
        let idx = self.node_index(level, pos);
        &mut self.nodes[idx]
    }

    /// Sets a node's data.
    #[inline]
    pub fn set_node(&mut self, level: usize, pos: usize, node: LdlNode) {
        let idx = self.node_index(level, pos);
        self.nodes[idx] = node;
    }

    /// Returns the number of nodes at a given level.
    #[inline]
    pub fn nodes_at_level(&self, level: usize) -> usize {
        1 << level
    }

    /// Gets the sigma value at a given level, position, and index.
    /// Returns sigma_min if the value is too small or not found.
    pub fn get_sigma(&self, level: usize, pos: usize, index: usize) -> f64 {
        let node = self.get_node(level, pos);
        node.sigma
            .get(index)
            .copied()
            .unwrap_or(self.sigma_min)
            .max(self.sigma_min)
    }

    /// Gets the l10 value at a given level, position, and index.
    pub fn get_l10(&self, level: usize, pos: usize, index: usize) -> Complex {
        let node = self.get_node(level, pos);
        node.l10.get(index).copied().unwrap_or(Complex::ZERO)
    }

    // ========================================================================
    // ffLDL Algorithm (FIPS 206)
    // ========================================================================

    /// Builds the LDL* tree from the secret key polynomials.
    ///
    /// Given the Gram matrix elements in FFT form, computes the LDL*
    /// decomposition recursively following the ffLDL algorithm.
    ///
    /// # Arguments
    ///
    /// * `g00` - G[0,0] = ||(f, g)||² in FFT form
    /// * `g01` - G[0,1] = <(f,g), (F,G)> in FFT form
    /// * `g11` - G[1,1] = ||(F, G)||² in FFT form
    ///
    /// Note: G[1,0] = conj(G[0,1]) for Hermitian matrices.
    pub fn build_from_gram(g00: &[Complex], g01: &[Complex], g11: &[Complex]) -> Self {
        let n = g00.len();
        let mut tree = LdlTree::new(n);

        // Start the recursive LDL decomposition
        tree.ffldl_recursive(0, 0, g00, g01, g11);

        tree
    }

    /// Builds the LDL* tree from the basis polynomials (f, g, F, G).
    ///
    /// First computes the Gram matrix G = [[<b0,b0>, <b0,b1>], [<b1,b0>, <b1,b1>]]
    /// where b0 = (f, g) and b1 = (F, G), then builds the LDL* tree.
    pub fn build_from_basis(
        f_fft: &[Complex],
        g_fft: &[Complex],
        big_f_fft: &[Complex],
        big_g_fft: &[Complex],
    ) -> Self {
        let n = f_fft.len();

        // Compute Gram matrix elements in FFT form
        // G[0,0] = |f|² + |g|² (norm squared of first basis vector)
        // G[0,1] = f·conj(F) + g·conj(G) (inner product)
        // G[1,1] = |F|² + |G|² (norm squared of second basis vector)
        let mut g00 = Vec::with_capacity(n);
        let mut g01 = Vec::with_capacity(n);
        let mut g11 = Vec::with_capacity(n);

        for i in 0..n {
            // G[0,0] = ||b0||² = |f_i|² + |g_i|²
            g00.push(Complex::from_real(
                f_fft[i].norm_sq() + g_fft[i].norm_sq(),
            ));

            // G[0,1] = <b0, b1> = f_i · conj(F_i) + g_i · conj(G_i)
            g01.push(f_fft[i] * big_f_fft[i].conj() + g_fft[i] * big_g_fft[i].conj());

            // G[1,1] = ||b1||² = |F_i|² + |G_i|²
            g11.push(Complex::from_real(
                big_f_fft[i].norm_sq() + big_g_fft[i].norm_sq(),
            ));
        }

        Self::build_from_gram(&g00, &g01, &g11)
    }

    /// Recursive ffLDL algorithm (FIPS 206).
    ///
    /// Computes the LDL* decomposition at the given tree position.
    fn ffldl_recursive(
        &mut self,
        level: usize,
        pos: usize,
        g00: &[Complex],
        g01: &[Complex],
        g11: &[Complex],
    ) {
        let n = g00.len();

        if n == 1 {
            // Base case: 2×2 LDL* decomposition
            self.ffldl_base_case(level, pos, g00[0], g01[0], g11[0]);
            return;
        }

        // Recursive case: split and recurse

        // Compute L and D factors for this level
        // d0 = G[0,0]
        // l10 = G[1,0] / G[0,0] = conj(G[0,1]) / G[0,0]
        // d1 = G[1,1] - |l10|² · G[0,0] = G[1,1] - G[1,0]·conj(G[1,0])/G[0,0]

        let mut l10 = Vec::with_capacity(n);
        let mut d0 = Vec::with_capacity(n);
        let mut d1 = Vec::with_capacity(n);

        for i in 0..n {
            // d0[i] = G00[i] (already a squared norm, so real)
            let d0_i = g00[i].re;
            d0.push(d0_i);

            // l10[i] = G10[i] / G00[i] = conj(G01[i]) / G00[i]
            // Since G00 is real (squared norm), division is simple
            let g10_i = g01[i].conj(); // G[1,0] = conj(G[0,1])
            let l10_i = if d0_i > 1e-10 {
                Complex::new(g10_i.re / d0_i, g10_i.im / d0_i)
            } else {
                Complex::ZERO
            };
            l10.push(l10_i);

            // d1[i] = G11[i] - |l10[i]|² · G00[i]
            //       = G11[i] - l10[i] · conj(l10[i]) · G00[i]
            //       = G11[i] - |G10[i]|² / G00[i]
            let d1_i = g11[i].re - l10_i.norm_sq() * d0_i;
            d1.push(d1_i.max(0.0)); // Ensure non-negative
        }

        // Compute sigma values (sqrt of D diagonal)
        let sigma: Vec<f64> = d0
            .iter()
            .zip(d1.iter())
            .flat_map(|(&d0_i, &d1_i)| vec![d0_i.sqrt(), d1_i.sqrt()])
            .collect();

        // Store this level's data
        self.set_node(level, pos, LdlNode::new(sigma, l10));

        // Split Gram matrix elements for recursion
        let (g00_0, g00_1) = split_fft(g00);
        let (g01_0, g01_1) = split_fft(g01);
        let (g11_0, g11_1) = split_fft(g11);

        // Recurse on left child (even indices)
        self.ffldl_recursive(level + 1, 2 * pos, &g00_0, &g01_0, &g11_0);

        // Recurse on right child (odd indices)
        self.ffldl_recursive(level + 1, 2 * pos + 1, &g00_1, &g01_1, &g11_1);
    }

    /// Base case LDL* decomposition for a single 2×2 Hermitian matrix.
    ///
    /// G = [[g00, g01], [g10, g11]] where g10 = conj(g01)
    ///
    /// LDL* gives:
    /// L = [[1, 0], [l10, 1]]
    /// D = [[d0, 0], [0, d1]]
    ///
    /// Where:
    /// - d0 = g00
    /// - l10 = g10 / g00
    /// - d1 = g11 - |l10|² · g00
    fn ffldl_base_case(
        &mut self,
        level: usize,
        pos: usize,
        g00: Complex,
        g01: Complex,
        g11: Complex,
    ) {
        // g00 and g11 should be real (squared norms)
        let d0 = g00.re;

        // l10 = g10 / g00 = conj(g01) / g00
        let g10 = g01.conj();
        let l10 = if d0 > 1e-10 {
            Complex::new(g10.re / d0, g10.im / d0)
        } else {
            Complex::ZERO
        };

        // d1 = g11 - |l10|² · g00
        let d1 = (g11.re - l10.norm_sq() * d0).max(0.0);

        // Sigma values
        let sigma0 = d0.sqrt().max(self.sigma_min);
        let sigma1 = d1.sqrt().max(self.sigma_min);

        self.set_node(level, pos, LdlNode::leaf(sigma0, sigma1, l10));
    }
}

// ============================================================================
// Legacy Compatibility (FftTree alias)
// ============================================================================

/// Alias for backward compatibility.
pub type FftTree = LdlTree;

/// Alias for backward compatibility.
pub type FftNode = LdlNode;

// ============================================================================
// Gram-Schmidt Data Structure
// ============================================================================

/// Gram-Schmidt data for the full basis.
///
/// Stores the secret key basis and its LDL* tree for sampling.
#[derive(Clone, Debug)]
pub struct GramSchmidt {
    /// The [[f, g], [F, G]] basis in FFT form.
    pub f_fft: Vec<Complex>,
    /// g polynomial in FFT form.
    pub g_fft: Vec<Complex>,
    /// F polynomial in FFT form.
    pub big_f_fft: Vec<Complex>,
    /// G polynomial in FFT form.
    pub big_g_fft: Vec<Complex>,

    /// The LDL* tree for sampling.
    pub tree: LdlTree,

    /// Precomputed sigma values for the (f, g) vector.
    pub sigma_fg: Vec<f64>,
    /// Precomputed sigma values for the (F, G) vector.
    pub sigma_fg_star: Vec<f64>,
}

impl GramSchmidt {
    /// Creates a new Gram-Schmidt structure from the basis polynomials.
    pub fn new(
        f_fft: Vec<Complex>,
        g_fft: Vec<Complex>,
        big_f_fft: Vec<Complex>,
        big_g_fft: Vec<Complex>,
    ) -> Self {
        let n = f_fft.len();

        // Build the LDL* tree
        let tree = LdlTree::build_from_basis(&f_fft, &g_fft, &big_f_fft, &big_g_fft);

        // Compute sigma values for (f, g) - these are sqrt(||b0||²)
        let sigma_fg: Vec<f64> = (0..n)
            .map(|i| (f_fft[i].norm_sq() + g_fft[i].norm_sq()).sqrt())
            .collect();

        // Compute orthogonalized sigma for (F, G) - sqrt(d1) from LDL*
        // This is ||b1*||² = ||b1||² - |<b0,b1>|²/||b0||²
        let sigma_fg_star: Vec<f64> = (0..n)
            .map(|i| {
                let norm_b0_sq = f_fft[i].norm_sq() + g_fft[i].norm_sq();
                let norm_b1_sq = big_f_fft[i].norm_sq() + big_g_fft[i].norm_sq();
                let inner = f_fft[i] * big_f_fft[i].conj() + g_fft[i] * big_g_fft[i].conj();
                let inner_sq = inner.norm_sq();

                if norm_b0_sq > 1e-10 {
                    (norm_b1_sq - inner_sq / norm_b0_sq).max(0.0).sqrt()
                } else {
                    norm_b1_sq.sqrt()
                }
            })
            .collect();

        GramSchmidt {
            f_fft,
            g_fft,
            big_f_fft,
            big_g_fft,
            tree,
            sigma_fg,
            sigma_fg_star,
        }
    }

    /// Gets the sigma value for the first basis vector at position i.
    pub fn sigma0(&self, i: usize) -> f64 {
        self.sigma_fg.get(i).copied().unwrap_or(1.0)
    }

    /// Gets the sigma value for the orthogonalized second basis vector at position i.
    pub fn sigma1(&self, i: usize) -> f64 {
        self.sigma_fg_star.get(i).copied().unwrap_or(1.0)
    }

    /// Gets the l10 coefficient for the Gram-Schmidt projection at position i.
    /// This is <b1, b0> / ||b0||²
    pub fn l10(&self, i: usize) -> Complex {
        let norm_b0_sq = self.f_fft[i].norm_sq() + self.g_fft[i].norm_sq();
        if norm_b0_sq > 1e-10 {
            let inner =
                self.big_f_fft[i] * self.f_fft[i].conj() + self.big_g_fft[i] * self.g_fft[i].conj();
            Complex::new(inner.re / norm_b0_sq, inner.im / norm_b0_sq)
        } else {
            Complex::ZERO
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fft::fft;

    #[test]
    fn test_ldl_tree_structure() {
        let tree = LdlTree::new(8);
        assert_eq!(tree.n, 8);
        assert_eq!(tree.depth, 4); // log2(8) + 1 = 4
        assert_eq!(tree.nodes.len(), 15); // 2*8 - 1 = 15
    }

    #[test]
    fn test_node_indexing() {
        let tree = LdlTree::new(8);

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
    fn test_ldl_base_case() {
        // Test 2×2 LDL* decomposition
        // G = [[4, 2+i], [2-i, 5]]
        let g00 = Complex::from_real(4.0);
        let g01 = Complex::new(2.0, 1.0);
        let g11 = Complex::from_real(5.0);

        let mut tree = LdlTree::new(1);
        tree.ffldl_base_case(0, 0, g00, g01, g11);

        let node = tree.get_node(0, 0);

        // d0 = 4
        // l10 = (2-i)/4 = 0.5 - 0.25i
        // d1 = 5 - |l10|² * 4 = 5 - (0.25 + 0.0625) * 4 = 5 - 1.25 = 3.75
        assert!((node.sigma[0] - 2.0).abs() < 1e-10); // sqrt(4)
        assert!((node.sigma[1] - 3.75_f64.sqrt()).abs() < 1e-10);
        assert!((node.l10[0].re - 0.5).abs() < 1e-10);
        assert!((node.l10[0].im - (-0.25)).abs() < 1e-10);
    }

    #[test]
    fn test_build_from_basis() {
        // Create simple test polynomials in FFT form
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

        let tree = LdlTree::build_from_basis(&f, &g, &big_f, &big_g);
        assert_eq!(tree.n, 4);
        assert_eq!(tree.depth, 3);

        // Check that all nodes have been populated
        for level in 0..tree.depth {
            for pos in 0..tree.nodes_at_level(level) {
                let node = tree.get_node(level, pos);
                assert!(!node.sigma.is_empty(), "Node at level {} pos {} should have sigma values", level, pos);
            }
        }
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
        assert_eq!(gs.sigma_fg_star.len(), n);
        assert!(gs.sigma_fg.iter().all(|&s| s > 0.0));
    }

    #[test]
    fn test_ldl_hermitian_property() {
        // Verify that the LDL* decomposition preserves the Hermitian property
        // G = L · D · L^* should equal the original Gram matrix

        let g00 = Complex::from_real(9.0);
        let g01 = Complex::new(3.0, 4.0);
        let g11 = Complex::from_real(10.0);

        let mut tree = LdlTree::new(1);
        tree.ffldl_base_case(0, 0, g00, g01, g11);

        let node = tree.get_node(0, 0);
        let d0 = node.sigma[0] * node.sigma[0];
        let d1 = node.sigma[1] * node.sigma[1];
        let l10 = node.l10[0];

        // Reconstruct G from LDL*
        // G[0,0] = d0
        // G[1,0] = l10 · d0
        // G[0,1] = conj(l10) · d0
        // G[1,1] = |l10|² · d0 + d1

        let g00_reconstructed = d0;
        let g10_reconstructed = l10.scale(d0);
        let g11_reconstructed = l10.norm_sq() * d0 + d1;

        assert!((g00_reconstructed - g00.re).abs() < 1e-10);
        assert!((g10_reconstructed.re - g01.conj().re).abs() < 1e-10);
        assert!((g10_reconstructed.im - g01.conj().im).abs() < 1e-10);
        assert!((g11_reconstructed - g11.re).abs() < 1e-10);
    }

    #[test]
    fn test_sigma_positivity() {
        // All sigma values should be positive (or at least sigma_min)
        let n = 8;
        let f_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0 + i as f64, 0.5))
            .collect();
        let g_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(0.5, 1.0 + i as f64))
            .collect();
        let big_f_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(2.0 + i as f64 * 0.5, 0.1))
            .collect();
        let big_g_fft: Vec<Complex> = (0..n)
            .map(|i| Complex::new(0.1, 2.0 + i as f64 * 0.5))
            .collect();

        let gs = GramSchmidt::new(f_fft, g_fft, big_f_fft, big_g_fft);

        // Check all sigma values are positive
        for &s in &gs.sigma_fg {
            assert!(s > 0.0, "sigma_fg should be positive");
        }
        for &s in &gs.sigma_fg_star {
            assert!(s >= 0.0, "sigma_fg_star should be non-negative");
        }
    }

    // Legacy compatibility test
    #[test]
    fn test_fft_tree_alias() {
        let tree: FftTree = LdlTree::new(4);
        assert_eq!(tree.n, 4);
    }
}
