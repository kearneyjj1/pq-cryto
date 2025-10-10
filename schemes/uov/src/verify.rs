use super::keygen::{F, QuadForm, PublicKey, Signature, hash_to_field};

#[inline] fn idx_ut(n: usize, i: usize, j: usize) -> usize { i*n - (i*(i+1))/2 + j }
#[inline] fn fadd(a:F,b:F)->F{F(a.0^b.0)} #[inline] fn fmul(a:F,b:F)->F{
    let (mut aa,mut bb,mut r)=(a.0,b.0,0u8);
    while bb!=0 { if (bb&1)!=0 { r^=aa; } let hi=(aa&0x80)!=0; aa<<=1; if hi{aa^=0x1b;} bb>>=1; }
    F(r)
}

fn eval_quad(q:&QuadForm, x:&[F]) -> F {
    let n=x.len(); let mut acc=F(0);
    for i in 0..n { for j in i..n {
        let c = q.coeffs[idx_ut(n,i,j)];
        acc = fadd(acc, fmul(c, fmul(x[i], x[j])));
    }}
    acc
}
fn eval_public(pk:&PublicKey, x:&[F]) -> Vec<F> {
    pk.polys.iter().map(|q| eval_quad(q,x)).collect()
}

pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    if sig.x.len() != pk.params.n() { return false; }
    let t  = hash_to_field(msg, &sig.salt, pk.params.m);
    let t2 = eval_public(pk, &sig.x);
    t == t2
}
