//! Modular inversion via divsteps (Bernstein-Yang safegcd).
//!
//! Port of libsecp256k1 modinv64. Uses signed 62-bit limbs; constant-time.

/// Signed 62-bit limb representation: value = sum(v[i] * 2^(62*i)), i=0..4.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Signed62 {
    pub v: [i64; 5],
}

/// Modulus info for secp256k1 field: p = 2^256 - 2^32 - 977.
/// modulus in signed62: {{-0x1000003D1, 0, 0, 0, 256}}
/// modulus^{-1} mod 2^62
pub struct ModInv64ModInfo {
    pub modulus: Signed62,
    pub modulus_inv62: u64,
}

/// secp256k1 field modulus: p = 2^256 - 2^32 - 977
pub const SECP256K1_FE_MODINV_MODINFO: ModInv64ModInfo = ModInv64ModInfo {
    modulus: Signed62 {
        v: [-0x1000003D1, 0, 0, 0, 256],
    },
    modulus_inv62: 0x27C7F6E22DDACACF,
};

/// secp256k1 scalar modulus (group order n): n = 2^256 - 0x14551231950B75FC4402DA1722FC9BAEE
pub const SECP256K1_SCALAR_MODINV_MODINFO: ModInv64ModInfo = ModInv64ModInfo {
    modulus: Signed62 {
        v: [0x3FD25E8CD0364141, 0x2ABB739ABD2280EE, -0x15, 0, 256],
    },
    modulus_inv62: 0x34F20099AA774EC1,
};

/// 2x2 transition matrix: t = [u v; q r]
#[derive(Clone, Copy)]
struct Trans2x2 {
    u: i64,
    v: i64,
    q: i64,
    r: i64,
}

const M62: u64 = u64::MAX >> 2;

/// 59 divsteps (constant-time). Returns new zeta.
#[inline(always)]
fn divsteps_59(zeta: i64, f0: u64, g0: u64, t: &mut Trans2x2) -> i64 {
    let mut u: u64 = 8;
    let mut v: u64 = 0;
    let mut q: u64 = 0;
    let mut r: u64 = 8;
    let mut f = f0;
    let mut g = g0;
    let mut zeta = zeta;

    for _ in 3..62 {
        let c1 = (zeta >> 63) as u64;
        let mask1 = c1;
        let c2 = g & 1;
        let mask2 = c2.wrapping_neg();

        let x = (f ^ mask1).wrapping_sub(mask1);
        let y = (u ^ mask1).wrapping_sub(mask1);
        let z = (v ^ mask1).wrapping_sub(mask1);

        g = g.wrapping_add(x & mask2);
        q = q.wrapping_add(y & mask2);
        r = r.wrapping_add(z & mask2);

        let mask1 = mask1 & mask2;
        zeta = (zeta ^ mask1 as i64).wrapping_sub(1);

        f = f.wrapping_add(g & mask1);
        u = u.wrapping_add(q & mask1);
        v = v.wrapping_add(r & mask1);

        g >>= 1;
        u <<= 1;
        v <<= 1;
    }

    t.u = u as i64;
    t.v = v as i64;
    t.q = q as i64;
    t.r = r as i64;
    zeta
}

/// Normalize r from (-2*modulus, modulus) to [0, modulus).
/// sign < 0 means negate. Limbs in (-2^62, 2^62) -> [0, 2^62).
#[inline(always)]
fn normalize_62(r: &mut Signed62, sign: i64, modinfo: &ModInv64ModInfo) {
    let m = &modinfo.modulus;
    let mut r0 = r.v[0];
    let mut r1 = r.v[1];
    let mut r2 = r.v[2];
    let mut r3 = r.v[3];
    let mut r4 = r.v[4];

    let m62 = M62 as i64;

    let cond_add = r4 >> 63;
    r0 += m.v[0] & cond_add;
    r1 += m.v[1] & cond_add;
    r2 += m.v[2] & cond_add;
    r3 += m.v[3] & cond_add;
    r4 += m.v[4] & cond_add;

    let cond_negate = sign >> 63;
    r0 = (r0 ^ cond_negate) - cond_negate;
    r1 = (r1 ^ cond_negate) - cond_negate;
    r2 = (r2 ^ cond_negate) - cond_negate;
    r3 = (r3 ^ cond_negate) - cond_negate;
    r4 = (r4 ^ cond_negate) - cond_negate;

    r1 += r0 >> 62;
    r0 &= m62;
    r2 += r1 >> 62;
    r1 &= m62;
    r3 += r2 >> 62;
    r2 &= m62;
    r4 += r3 >> 62;
    r3 &= m62;

    let cond_add = r4 >> 63;
    r0 += m.v[0] & cond_add;
    r1 += m.v[1] & cond_add;
    r2 += m.v[2] & cond_add;
    r3 += m.v[3] & cond_add;
    r4 += m.v[4] & cond_add;

    r1 += r0 >> 62;
    r0 &= m62;
    r2 += r1 >> 62;
    r1 &= m62;
    r3 += r2 >> 62;
    r2 &= m62;
    r4 += r3 >> 62;
    r3 &= m62;

    r.v[0] = r0;
    r.v[1] = r1;
    r.v[2] = r2;
    r.v[3] = r3;
    r.v[4] = r4;
}

/// Update [d,e] = (t/2^62) * [d,e] mod modulus.
#[inline(always)]
fn update_de_62(d: &mut Signed62, e: &mut Signed62, t: &Trans2x2, modinfo: &ModInv64ModInfo) {
    let d0 = d.v[0];
    let d1 = d.v[1];
    let d2 = d.v[2];
    let d3 = d.v[3];
    let d4 = d.v[4];
    let e0 = e.v[0];
    let e1 = e.v[1];
    let e2 = e.v[2];
    let e3 = e.v[3];
    let e4 = e.v[4];
    let u = t.u;
    let v = t.v;
    let q = t.q;
    let r = t.r;
    let m = &modinfo.modulus;
    let mi = modinfo.modulus_inv62;

    let sd = d4 >> 63;
    let se = e4 >> 63;
    let mut md = (u & sd) + (v & se);
    let mut me = (q & sd) + (r & se);

    let mut cd: i128 = (u as i128) * (d0 as i128) + (v as i128) * (e0 as i128);
    let mut ce: i128 = (q as i128) * (d0 as i128) + (r as i128) * (e0 as i128);

    md -= (mi.wrapping_mul(cd as u64).wrapping_add(md as u64) & M62) as i64;
    me -= (mi.wrapping_mul(ce as u64).wrapping_add(me as u64) & M62) as i64;

    cd += (m.v[0] as i128) * (md as i128);
    ce += (m.v[0] as i128) * (me as i128);

    cd >>= 62;
    ce >>= 62;

    cd += (u as i128) * (d1 as i128) + (v as i128) * (e1 as i128);
    ce += (q as i128) * (d1 as i128) + (r as i128) * (e1 as i128);
    if m.v[1] != 0 {
        cd += (m.v[1] as i128) * (md as i128);
        ce += (m.v[1] as i128) * (me as i128);
    }
    d.v[0] = ((cd as u64) & M62) as i64;
    e.v[0] = ((ce as u64) & M62) as i64;
    cd >>= 62;
    ce >>= 62;

    cd += (u as i128) * (d2 as i128) + (v as i128) * (e2 as i128);
    ce += (q as i128) * (d2 as i128) + (r as i128) * (e2 as i128);
    if m.v[2] != 0 {
        cd += (m.v[2] as i128) * (md as i128);
        ce += (m.v[2] as i128) * (me as i128);
    }
    d.v[1] = ((cd as u64) & M62) as i64;
    e.v[1] = ((ce as u64) & M62) as i64;
    cd >>= 62;
    ce >>= 62;

    cd += (u as i128) * (d3 as i128) + (v as i128) * (e3 as i128);
    ce += (q as i128) * (d3 as i128) + (r as i128) * (e3 as i128);
    if m.v[3] != 0 {
        cd += (m.v[3] as i128) * (md as i128);
        ce += (m.v[3] as i128) * (me as i128);
    }
    d.v[2] = ((cd as u64) & M62) as i64;
    e.v[2] = ((ce as u64) & M62) as i64;
    cd >>= 62;
    ce >>= 62;

    cd += (u as i128) * (d4 as i128) + (v as i128) * (e4 as i128) + (m.v[4] as i128) * (md as i128);
    ce += (q as i128) * (d4 as i128) + (r as i128) * (e4 as i128) + (m.v[4] as i128) * (me as i128);
    d.v[3] = ((cd as u64) & M62) as i64;
    e.v[3] = ((ce as u64) & M62) as i64;
    cd >>= 62;
    ce >>= 62;

    d.v[4] = cd as i64;
    e.v[4] = ce as i64;
}

/// Update [f,g] = (t/2^62) * [f,g].
#[inline(always)]
fn update_fg_62(f: &mut Signed62, g: &mut Signed62, t: &Trans2x2) {
    let f0 = f.v[0];
    let f1 = f.v[1];
    let f2 = f.v[2];
    let f3 = f.v[3];
    let f4 = f.v[4];
    let g0 = g.v[0];
    let g1 = g.v[1];
    let g2 = g.v[2];
    let g3 = g.v[3];
    let g4 = g.v[4];
    let u = t.u;
    let v = t.v;
    let q = t.q;
    let r = t.r;

    let mut cf: i128 = (u as i128) * (f0 as i128) + (v as i128) * (g0 as i128);
    let mut cg: i128 = (q as i128) * (f0 as i128) + (r as i128) * (g0 as i128);
    cf >>= 62;
    cg >>= 62;

    cf += (u as i128) * (f1 as i128) + (v as i128) * (g1 as i128);
    cg += (q as i128) * (f1 as i128) + (r as i128) * (g1 as i128);
    f.v[0] = ((cf as u64) & M62) as i64;
    g.v[0] = ((cg as u64) & M62) as i64;
    cf >>= 62;
    cg >>= 62;

    cf += (u as i128) * (f2 as i128) + (v as i128) * (g2 as i128);
    cg += (q as i128) * (f2 as i128) + (r as i128) * (g2 as i128);
    f.v[1] = ((cf as u64) & M62) as i64;
    g.v[1] = ((cg as u64) & M62) as i64;
    cf >>= 62;
    cg >>= 62;

    cf += (u as i128) * (f3 as i128) + (v as i128) * (g3 as i128);
    cg += (q as i128) * (f3 as i128) + (r as i128) * (g3 as i128);
    f.v[2] = ((cf as u64) & M62) as i64;
    g.v[2] = ((cg as u64) & M62) as i64;
    cf >>= 62;
    cg >>= 62;

    cf += (u as i128) * (f4 as i128) + (v as i128) * (g4 as i128);
    cg += (q as i128) * (f4 as i128) + (r as i128) * (g4 as i128);
    f.v[3] = ((cf as u64) & M62) as i64;
    g.v[3] = ((cg as u64) & M62) as i64;
    cf >>= 62;
    cg >>= 62;

    f.v[4] = cf as i64;
    g.v[4] = cg as i64;
}

/// Compute modular inverse of x. Constant-time in x.
/// x must be in [0, modulus). If x=0, result is 0.
#[inline(always)]
pub fn modinv64(x: &mut Signed62, modinfo: &ModInv64ModInfo) {
    let mut d = Signed62 { v: [0, 0, 0, 0, 0] };
    let mut e = Signed62 { v: [1, 0, 0, 0, 0] };
    let mut f = modinfo.modulus;
    let mut g = *x;
    let mut zeta: i64 = -1;

    for _ in 0..10 {
        let mut t = Trans2x2 {
            u: 0,
            v: 0,
            q: 0,
            r: 0,
        };
        zeta = divsteps_59(zeta, f.v[0] as u64, g.v[0] as u64, &mut t);
        update_de_62(&mut d, &mut e, &t, modinfo);
        update_fg_62(&mut f, &mut g, &t);
    }

    normalize_62(&mut d, f.v[4], modinfo);
    *x = d;
}
