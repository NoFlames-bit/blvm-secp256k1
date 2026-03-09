//! BIP 340 Schnorr half-aggregation.
//!
//! Non-interactive aggregation of n BIP 340 signatures into a single aggregate
//! of 32*(n+1) bytes (~half combined size). Third party can aggregate.
//!
//! Reference: BlockstreamResearch/cross-input-aggregation (half-aggregation.mediawiki)

use crate::ecmult;
use crate::group::Gej;
use crate::scalar::Scalar;
use crate::schnorr::{lift_x, tagged_hash, tagged_hash_challenge};

const TAG_HALFAGG_RANDOMIZER: &[u8] = b"HalfAgg/randomizer";

/// (pubkey_x32, msg, sig64)
pub type PmS = ([u8; 32], Vec<u8>, [u8; 64]);
/// (pubkey_x32, msg)
pub type Pm = ([u8; 32], Vec<u8>);

/// Compute randomizer z_i for index i. z_0 = 1; else z_i = hash_halfagg(r_0||pk_0||m_0||...||r_i||pk_i||m_i) mod n.
fn randomizer(pmr: &[(/*r*/ [u8; 32], /*pk*/ [u8; 32], /*m*/ &[u8])], index: usize) -> Scalar {
    if index == 0 {
        let mut one = Scalar::zero();
        one.set_int(1);
        return one;
    }
    let mut data = Vec::new();
    for p in pmr.iter().take(index + 1) {
        data.extend_from_slice(&p.0);
        data.extend_from_slice(&p.1);
        data.extend_from_slice(p.2);
    }
    let h = tagged_hash(TAG_HALFAGG_RANDOMIZER, &data);
    let mut z = Scalar::zero();
    z.set_b32(&h);
    z
}

/// Aggregate BIP 340 signatures. Returns aggregate as Vec<u8> or None on failure.
///
/// Input: slice of (pubkey_x32, msg, sig64). Messages can be arbitrary length.
/// Output: r_0 || r_1 || ... || r_u || bytes(s), length 32*(u+1).
///
/// Caller should verify individual signatures before aggregating if invalid-input
/// acceptance is undesired (aggregate can accept some invalid triples).
pub fn aggregate(pms: &[PmS]) -> Option<Vec<u8>> {
    let empty: Vec<u8> = vec![0u8; 32];
    inc_aggregate(&empty, &[], pms)
}

/// Incrementally aggregate additional signatures into an existing aggregate.
///
/// - `aggsig`: existing aggregate (32*(v+1) bytes for v already-aggregated (pk,msg) pairs)
/// - `pm_aggd`: (pk, msg) pairs corresponding to aggsig
/// - `pms_to_agg`: additional (pk, msg, sig) triples to aggregate
pub fn inc_aggregate(aggsig: &[u8], pm_aggd: &[Pm], pms_to_agg: &[PmS]) -> Option<Vec<u8>> {
    let v = pm_aggd.len();
    let u = pms_to_agg.len();
    if v.saturating_add(u) >= 0x10000 {
        return None;
    }
    if aggsig.len() != 32 * (v + 1) {
        return None;
    }

    let mut pmr: Vec<([u8; 32], [u8; 32], Vec<u8>)> = Vec::with_capacity(v + u);

    for i in 0..v {
        let (pk, msg) = &pm_aggd[i];
        let r: [u8; 32] = aggsig[32 * i..32 * (i + 1)].try_into().ok()?;
        pmr.push((r, *pk, msg.clone()));
    }

    let mut s_bytes = [0u8; 32];
    s_bytes.copy_from_slice(&aggsig[32 * v..32 * (v + 1)]);
    let mut s = Scalar::zero();
    s.set_b32(&s_bytes);

    for i in v..(v + u) {
        let (pk, msg, sig) = &pms_to_agg[i - v];
        let r: [u8; 32] = sig[0..32].try_into().unwrap();
        pmr.push((r, *pk, msg.clone()));

        let z = randomizer(
            &pmr.iter()
                .map(|(r, pk, m)| (*r, *pk, m.as_slice()))
                .collect::<Vec<_>>(),
            i,
        );

        let si_bytes: [u8; 32] = sig[32..64].try_into().unwrap();
        let mut si = Scalar::zero();
        if si.set_b32(&si_bytes) {
            return None; // s >= n
        }

        let mut z_si = Scalar::zero();
        z_si.mul(&z, &si);
        let mut s_new = Scalar::zero();
        s_new.add(&s, &z_si); // reduces mod n on overflow
        s = s_new;
    }

    let mut out = Vec::with_capacity(32 * (v + u + 1));
    for (r, _, _) in &pmr {
        out.extend_from_slice(r);
    }
    let mut s_out = [0u8; 32];
    s.get_b32(&mut s_out);
    out.extend_from_slice(&s_out);
    Some(out)
}

/// Verify a half-aggregate signature against (pk, msg) pairs.
///
/// Returns true iff verification passes. aggsig must be 32*(u+1) bytes for u pairs.
pub fn verify_aggregate(aggsig: &[u8], pm_aggd: &[Pm]) -> bool {
    let u = pm_aggd.len();
    if u >= 0x10000 {
        return false;
    }
    if aggsig.len() != 32 * (u + 1) {
        return false;
    }

    let mut scalars = Vec::with_capacity(2 * u);
    let mut points = Vec::with_capacity(2 * u);

    for i in 0..u {
        let (pk, msg) = &pm_aggd[i];
        let p = match lift_x(pk) {
            Some(p) => p,
            None => return false,
        };

        let r_bytes: [u8; 32] = aggsig[32 * i..32 * (i + 1)].try_into().unwrap();
        let r = match lift_x(&r_bytes) {
            Some(r) => r,
            None => return false,
        };

        let mut challenge_data = Vec::with_capacity(32 + 32 + msg.len());
        challenge_data.extend_from_slice(&r_bytes);
        challenge_data.extend_from_slice(pk);
        challenge_data.extend_from_slice(msg);
        let e_hash = tagged_hash_challenge(&challenge_data);
        let mut e = Scalar::zero();
        e.set_b32(&e_hash);

        let pmr_ref: Vec<([u8; 32], [u8; 32], &[u8])> = (0..=i)
            .map(|j| {
                let (pk_j, msg_j) = &pm_aggd[j];
                let r_j: [u8; 32] = aggsig[32 * j..32 * (j + 1)].try_into().unwrap();
                (r_j, *pk_j, msg_j.as_slice())
            })
            .collect();
        let z = randomizer(&pmr_ref, i);

        let mut z_e = Scalar::zero();
        z_e.mul(&z, &e);

        scalars.push(z);
        points.push(r);
        scalars.push(z_e);
        points.push(p);
    }

    let s_bytes: [u8; 32] = aggsig[32 * u..32 * (u + 1)].try_into().unwrap();
    let mut s = Scalar::zero();
    if s.set_b32(&s_bytes) {
        return false; // s >= n
    }

    let zero = Scalar::zero();
    let mut sum_left = Gej::default();
    ecmult::ecmult_multi(&mut sum_left, &zero, &scalars, &points);

    // s*G
    let mut s_g = Gej::default();
    ecmult::ecmult_gen(&mut s_g, &s);

    // Check sum_left == s*G
    let mut neg_s_g = Gej::default();
    neg_s_g.neg(&s_g);
    let mut diff = Gej::default();
    diff.add_var(&sum_left, &neg_s_g);
    diff.is_infinity()
}
