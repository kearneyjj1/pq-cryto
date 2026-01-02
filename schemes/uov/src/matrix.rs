//! Matrix operations over GF(2^8).
//!
//! This module provides matrix utilities for the UOV signature scheme,
//! including transposition, multiplication, rank computation, and
//! conversions between upper triangular and full matrix representations.

use rand::{CryptoRng, RngCore};

use crate::field::F;

/// Computes the index into an upper triangular matrix stored as a flat vector.
///
/// For an n×n upper triangular matrix, element (i,j) where i <= j is stored at
/// position i*n - i*(i+1)/2 + j.
#[inline]
pub fn idx_ut(n: usize, i: usize, j: usize) -> usize {
    debug_assert!(i <= j, "idx_ut requires i <= j");
    i * n - (i * (i + 1)) / 2 + j
}

/// Transposes a square matrix.
pub fn transpose(a: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = a.len();
    let mut t = vec![vec![F::ZERO; n]; n];
    for i in 0..n {
        for j in 0..n {
            t[j][i] = a[i][j];
        }
    }
    t
}

/// Multiplies two square matrices over GF(256).
pub fn mat_mul(a: &[Vec<F>], b: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = a.len();
    let mut c = vec![vec![F::ZERO; n]; n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i][k];
            if !aik.is_zero() {
                for j in 0..n {
                    c[i][j] += aik * b[k][j];
                }
            }
        }
    }
    c
}

/// Converts an upper triangular matrix (packed as a vector) to a symmetric full matrix.
///
/// The upper triangular coefficients are reflected to create a symmetric matrix.
pub fn ut_to_full(n: usize, ut: &[F]) -> Vec<Vec<F>> {
    let mut m = vec![vec![F::ZERO; n]; n];
    for i in 0..n {
        for j in i..n {
            let c = ut[idx_ut(n, i, j)];
            m[i][j] = c;
            m[j][i] = c;
        }
    }
    m
}

/// Converts a symmetric full matrix to upper triangular packed representation.
pub fn full_to_ut(n: usize, m: &[Vec<F>]) -> Vec<F> {
    let mut ut = vec![F::ZERO; n * (n + 1) / 2];
    for i in 0..n {
        for j in i..n {
            ut[idx_ut(n, i, j)] = m[i][j];
        }
    }
    ut
}

/// Computes the rank of a matrix using Gaussian elimination.
///
/// Note: This function modifies the input matrix in place.
pub fn rank(a: &mut [Vec<F>]) -> usize {
    let n = a.len();
    if n == 0 {
        return 0;
    }
    let cols = a[0].len();
    let mut r = 0;

    for col in 0..cols {
        // Find pivot
        let mut piv = None;
        for row in r..n {
            if !a[row][col].is_zero() {
                piv = Some(row);
                break;
            }
        }

        if let Some(p) = piv {
            a.swap(r, p);
            let inv = a[r][col].inverse();

            // Scale pivot row
            for j in col..cols {
                a[r][j] *= inv;
            }

            // Eliminate column
            for i in 0..n {
                if i == r {
                    continue;
                }
                let f = a[i][col];
                if !f.is_zero() {
                    for j in col..cols {
                        let pivot_val = a[r][j];
                        a[i][j] -= f * pivot_val;
                    }
                }
            }
            r += 1;
        }
    }

    r
}

/// Samples a random invertible n×n matrix over GF(256).
///
/// Uses rejection sampling: generates random matrices until one with full rank is found.
pub fn sample_invertible<R: RngCore + CryptoRng>(rng: &mut R, n: usize) -> Vec<Vec<F>> {
    loop {
        let mut a = vec![vec![F::ZERO; n]; n];
        for row in a.iter_mut() {
            for elem in row.iter_mut() {
                *elem = F(rng.next_u32() as u8);
            }
        }
        let mut tmp = a.clone();
        if rank(&mut tmp) == n {
            return a;
        }
    }
}

/// Solves the linear system Ax = b over GF(256) using Gaussian elimination.
///
/// Returns `None` if the system is singular (no unique solution).
pub fn solve(mut a: Vec<Vec<F>>, mut b: Vec<F>) -> Option<Vec<F>> {
    let n = b.len();
    if n == 0 {
        return Some(vec![]);
    }

    // Forward elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let mut piv = None;
        for r in col..n {
            if !a[r][col].is_zero() {
                piv = Some(r);
                break;
            }
        }
        let p = piv?;

        // Swap rows
        if p != col {
            a.swap(p, col);
            b.swap(p, col);
        }

        // Scale pivot row
        let inv = a[col][col].inverse();
        for j in col..n {
            a[col][j] *= inv;
        }
        b[col] *= inv;

        // Eliminate column
        for i in 0..n {
            if i == col {
                continue;
            }
            let f = a[i][col];
            if !f.is_zero() {
                for j in col..n {
                    let pivot_val = a[col][j];
                    a[i][j] -= f * pivot_val;
                }
                let b_pivot = b[col];
                b[i] -= f * b_pivot;
            }
        }
    }

    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_idx_ut() {
        // For n=3, upper triangular indices should be:
        // (0,0)=0, (0,1)=1, (0,2)=2, (1,1)=3, (1,2)=4, (2,2)=5
        assert_eq!(idx_ut(3, 0, 0), 0);
        assert_eq!(idx_ut(3, 0, 1), 1);
        assert_eq!(idx_ut(3, 0, 2), 2);
        assert_eq!(idx_ut(3, 1, 1), 3);
        assert_eq!(idx_ut(3, 1, 2), 4);
        assert_eq!(idx_ut(3, 2, 2), 5);
    }

    #[test]
    fn test_transpose() {
        let a = vec![
            vec![F(1), F(2), F(3)],
            vec![F(4), F(5), F(6)],
            vec![F(7), F(8), F(9)],
        ];
        let t = transpose(&a);
        assert_eq!(t[0][1], F(4));
        assert_eq!(t[1][0], F(2));
        assert_eq!(t[2][0], F(3));
    }

    #[test]
    fn test_ut_full_roundtrip() {
        let n = 3;
        let ut = vec![F(1), F(2), F(3), F(4), F(5), F(6)];
        let full = ut_to_full(n, &ut);
        let ut2 = full_to_ut(n, &full);
        assert_eq!(ut, ut2);
    }

    #[test]
    fn test_solve_simple() {
        // Simple 2x2 system: x + y = 3, x = 1 => y = 2
        let a = vec![vec![F(1), F(1)], vec![F(1), F(0)]];
        let b = vec![F(3), F(1)];
        let x = solve(a, b).unwrap();
        assert_eq!(x[0], F(1));
        assert_eq!(x[1], F(2));
    }

    #[test]
    fn test_solve_singular() {
        // Singular system: two identical rows
        let a = vec![vec![F(1), F(1)], vec![F(1), F(1)]];
        let b = vec![F(1), F(2)];
        assert!(solve(a, b).is_none());
    }
}
