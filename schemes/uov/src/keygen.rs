use rand::{CryptoRng, RngCore};
use sha3::{digest::{Update, ExtendableOutput, XofReader}, Shake256};

/// GF(256) with AES polynomial x^8 + x^4 + x^3 + x + 1 (0x11B)
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct F(pub u8);
#[inline] fn fadd(a: F, b: F) -> F { F(a.0 ^ b.0) }
#[inline] fn fsub(a: F, b: F) -> F { fadd(a,b) }
#[inline] fn fmul(a: F, b: F) -> F {                // <- removed `mut` params
    let mut r = 0u8; let mut aa=a.0; let mut bb=b.0;
    while bb != 0 { if (bb & 1)!=0 { r ^= aa; } let hi=(aa&0x80)!=0; aa<<=1; if hi{aa^=0x1b;} bb>>=1; }
    F(r)
}
#[inline] fn fsq(x: F) -> F { fmul(x,x) }
#[inline] fn finv(a: F) -> F {
    let x2=fsq(a); let x4=fsq(x2); let x8=fsq(x4); let x16=fsq(x8);
    let x32=fsq(x16); let x64=fsq(x32); let x128=fsq(x64);
    let mut r = x128; r = fmul(r,x64); r = fmul(r,x32); r = fmul(r,x16); r = fmul(r,x8); r = fmul(r,x4); r = fmul(r,x2); r
}

#[derive(Clone, Copy)]                                 // <- make Params Copy
pub struct Params { pub v: usize, pub m: usize }  // n = v + m
impl Params { pub fn n(&self)->usize { self.v + self.m } }

#[derive(Clone)]
pub struct QuadForm { pub coeffs: Vec<F> } // length n*(n+1)/2

#[derive(Clone)]
pub struct PublicKey { pub params: Params, pub polys: Vec<QuadForm> }

pub struct SecretKey {
    pub params: Params,
    pub t: Vec<Vec<F>>,            // invertible linear map T (n x n)
    pub f_quads: Vec<Vec<F>>,      // central OV quadratic forms (UT packed)
}

pub struct Signature { pub x: Vec<F>, pub salt: Vec<u8> }

#[inline] fn idx_ut(n: usize, i: usize, j: usize) -> usize { i*n - (i*(i+1))/2 + j }

pub fn hash_to_field(msg: &[u8], salt: &[u8], m: usize) -> Vec<F> {
    let mut h = Shake256::default();
    h.update(b"pqsigs-uov-v1");
    h.update(msg);
    h.update(salt);
    let mut xof = h.finalize_xof();
    let mut out = vec![0u8; m];
    xof.read(&mut out);
    out.into_iter().map(F).collect()
}

fn transpose(a: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = a.len(); let mut t = vec![vec![F(0); n]; n];
    for i in 0..n { for j in 0..n { t[j][i] = a[i][j]; } } t
}
fn mat_mul(a: &[Vec<F>], b: &[Vec<F>]) -> Vec<Vec<F>> {
    let n=a.len(); let mut c=vec![vec![F(0);n];n];
    for i in 0..n { for k in 0..n { let aik=a[i][k]; for j in 0..n { c[i][j]=fadd(c[i][j], fmul(aik,b[k][j])); }}}
    c
}
fn ut_to_full(n: usize, ut: &[F]) -> Vec<Vec<F>> {
    let mut m = vec![vec![F(0); n]; n];
    for i in 0..n { for j in i..n { let c = ut[idx_ut(n,i,j)]; m[i][j]=c; m[j][i]=c; } } m
}
fn full_to_ut(n: usize, m: &[Vec<F>]) -> Vec<F> {
    let mut ut = vec![F(0); n*(n+1)/2];
    for i in 0..n { for j in i..n { ut[idx_ut(n,i,j)] = m[i][j]; } } ut
}
fn rank(a: &mut [Vec<F>]) -> usize {
    let n=a.len(); let mut r=0usize;
    for col in 0..n {
        let mut piv=None;
        for row in r..n { if a[row][col].0!=0 { piv=Some(row); break; } }
        if let Some(p) = piv {
            a.swap(r,p);
            let inv = finv(a[r][col]);
            for j in col..n { a[r][j] = fmul(a[r][j], inv); }
            for i in 0..n {
                if i==r {continue;}
                let f = a[i][col];
                if f.0!=0 { for j in col..n { a[i][j] = fsub(a[i][j], fmul(f, a[r][j])); } }
            }
            r+=1;
        }
    }
    r
}
fn sample_gl<R: RngCore + CryptoRng>(rng: &mut R, n: usize) -> Vec<Vec<F>> {
    loop {
        let mut a = vec![vec![F(0); n]; n];
        for i in 0..n { for j in 0..n { a[i][j] = F(rng.next_u32() as u8); } }
        let mut tmp = a.clone();
        if rank(&mut tmp) == n { return a; }
    }
}

/// Demo Level-1ish params (we’ll swap to the exact Round-2 set next)
pub const PARAMS_L1: Params = Params { v: 16, m: 8 };

/// Build central OV map F (no oil*oil terms), sample T, and publish P(x) = F(Tx)
pub fn keygen<R: RngCore + CryptoRng>(rng: &mut R, params: Params) -> (PublicKey, SecretKey) {
    let n = params.n(); let ut_len = n*(n+1)/2; let v = params.v;
    // central OV forms
    let mut f_quads = Vec::with_capacity(params.m);
    for _ in 0..params.m {
        let mut a = vec![F(0); ut_len];
        for i in 0..n {
            for j in i..n {
                let oil_i = i >= v; let oil_j = j >= v;
                if oil_i && oil_j { continue; }
                a[idx_ut(n,i,j)] = F(rng.next_u32() as u8);
            }
        }
        f_quads.push(a);
    }
    let t  = sample_gl(rng, n);
    let tT = transpose(&t);
    // compose P(x) = F(Tx) via A' = T^T A T
    let mut polys = Vec::with_capacity(params.m);
    for k in 0..params.m {
        let a_full = ut_to_full(n, &f_quads[k]);
        let tmp = mat_mul(&tT, &a_full);
        let a_p = mat_mul(&tmp, &t);
        polys.push(QuadForm { coeffs: full_to_ut(n, &a_p) });
    }
    (PublicKey { params, polys }, SecretKey { params, t, f_quads })
}
