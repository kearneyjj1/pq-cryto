//! Classic UOV signing (demo): pick random vinegar, solve linear system for oil.

use rand::{RngCore, CryptoRng};
use super::keygen::{Params, F, PublicKey, SecretKey, Signature, hash_to_field};

#[inline] fn fadd(a:F,b:F)->F{F(a.0^b.0)} #[inline] fn fsub(a:F,b:F)->F{F(a.0^b.0)}
#[inline] fn fmul(mut a:F, mut b:F)->F{
    let (mut aa,mut bb,mut r)=(a.0,b.0,0u8);
    while bb!=0 { if (bb&1)!=0 { r^=aa; } let hi=(aa&0x80)!=0; aa<<=1; if hi{aa^=0x1b;} bb>>=1; } F(r)
}
#[inline] fn fsq(x:F)->F{ fmul(x,x) }
#[inline] fn finv(a:F)->F{
    let x2=fsq(a); let x4=fsq(x2); let x8=fsq(x4); let x16=fsq(x8);
    let x32=fsq(x16); let x64=fsq(x32); let x128=fsq(x64);
    let mut r=x128; r=fmul(r,x64); r=fmul(r,x32); r=fmul(r,x16); r=fmul(r,x8); r=fmul(r,x4); r=fmul(r,x2); r
}
#[inline] fn idx_ut(n:usize,i:usize,j:usize)->usize{ i*n - (i*(i+1))/2 + j }

/// Solve A x = b for x over GF(256) (dense Gaussian, returns None if singular).
fn solve(mut a: Vec<Vec<F>>, mut b: Vec<F>) -> Option<Vec<F>> {
    let n = b.len();
    // forward
    for col in 0..n {
        // pivot
        let mut piv = None;
        for r in col..n { if a[r][col].0 != 0 { piv = Some(r); break; } }
        let r = piv?;
        if r != col { a.swap(r,col); b.swap(r,col); }
        // scale
        let inv = finv(a[col][col]);
        for j in col..n { a[col][j] = fmul(a[col][j], inv); }
        b[col] = fmul(b[col], inv);
        // eliminate
        for i in 0..n {
            if i==col { continue; }
            let f = a[i][col];
            if f.0 != 0 {
                for j in col..n { a[i][j] = fsub(a[i][j], fmul(f, a[col][j])); }
                b[i] = fsub(b[i], fmul(f, b[col]));
            }
        }
    }
    Some(b)
}

/// Sign by inverting central OV map:
/// 1) Choose random vinegar v
/// 2) Build linear system in oil variables so that F(Tx)=t
/// 3) Solve; output x plus salt
pub fn sign<R: RngCore + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    sk: &SecretKey,
    msg: &[u8],
) -> Option<Signature> {
    let Params{v,m} = sk.params; let n = v + m;
    // salt
    let mut salt = vec![0u8; m];
    rng.fill_bytes(&mut salt);
    let t = super::keygen::hash_to_field(msg, &salt, m);

    // Try a few vinegar attempts
    for _attempt in 0..64 {
        // 1) choose vinegar
        let mut x = vec![F(0); n];
        for i in 0..v { x[i] = F(rng.next_u32() as u8); }

        // 2) Work in y = T x  (we will solve for y to satisfy central F(y)=t, then x = T^{-1} y)
        // Build A_oil * y_oil = b  (m eqs, m unknowns)
        // Central forms are in sk.f_quads (UT). Evaluate vinegar-vinegar and vinegar-oil contributions.
        let mut a_mat = vec![vec![F(0); m]; m];
        let mut b_vec = t.clone();

        for eq in 0..m {
            let ut = &sk.f_quads[eq];
            // subtract vv and v*o cross-terms dependent on vinegar (move to RHS)
            // Compute constant part c = sum_{i<=j < v} a_ij x_i x_j  +  sum_{i<v, j>=v} a_ij x_i * y_j  (the latter is linear in y_oil; move to LHS)
            // Build LHS coefficients for y_oil and RHS for constants.
            // First vv:
            let mut const_term = F(0);
            for i in 0..v {
                for j in i..v {
                    let aij = ut[idx_ut(n,i,j)];
                    const_term = fadd(const_term, fmul(aij, fmul(x[i], x[j])));
                }
            }
            // Move vv to RHS: t := t - const_term
            b_vec[eq] = fsub(b_vec[eq], const_term);

            // Now vinegar-oil terms: sum_{i<v, j>=v} a_ij x_i y_j  (linear in y_oil)
            for i in 0..v {
                for j in v..n {
                    let aij = ut[idx_ut(n,i,j)];
                    let coeff = fmul(aij, x[i]);
                    // add to column (j-v)
                    a_mat[eq][j - v] = fadd(a_mat[eq][j - v], coeff);
                }
            }

            // Oil-oil block is zero by construction in central map (classic UOV), so no quadratic y_oil*y_oil terms.
        }

        // Solve for y_oil
        if let Some(y_oil) = solve(a_mat, b_vec) {
            // Compose y = [vinegar | oil]
            let mut y = vec![F(0); n];
            for i in 0..v { y[i] = x[i]; }
            for j in 0..m { y[v + j] = y_oil[j]; }

            // Compute x = T^{-1} y
            let x = solve_t_inv(&sk.t, &y)?;
            return Some(Signature { x, salt });
        }
    }
    None
}

/// Solve T x = y for x (dense Gauss with RHS vector)
fn solve_t_inv(t: &[Vec<F>], y: &[F]) -> Option<Vec<F>> {
    let n = t.len();
    let mut a = t.to_vec();
    let mut b = y.to_vec();
    // standard solve
    super::sign::solve(a, b)
}
