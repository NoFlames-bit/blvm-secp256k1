//! Tests for BIP 340 Schnorr half-aggregation.
//! Vectors from BlockstreamResearch/cross-input-aggregation hacspec-halfagg.

use blvm_secp256k1::schnorr;
use blvm_secp256k1::schnorr_half_agg;

fn hex32(s: &str) -> [u8; 32] {
    let v = hex::decode(s).unwrap();
    v.try_into().unwrap()
}

#[test]
fn test_verify_vectors() {
    // Empty aggregate
    let aggsig =
        hex::decode("0000000000000000000000000000000000000000000000000000000000000000").unwrap();
    let pm: Vec<([u8; 32], Vec<u8>)> = vec![];
    assert!(schnorr_half_agg::verify_aggregate(&aggsig, &pm));

    // 1 signature
    let pm1 = vec![(
        hex32("1b84c5567b126440995d3ed5aaba0565d71e1834604819ff9c17f5e9d5dd078f"),
        hex32("0202020202020202020202020202020202020202020202020202020202020202").to_vec(),
    )];
    let aggsig1 = hex::decode("b070aafcea439a4f6f1bbfc2eb66d29d24b0cab74d6b745c3cfb009cc8fe4aa80e066c34819936549ff49b6fd4d41edfc401a367b87ddd59fee38177961c225f").unwrap();
    assert!(schnorr_half_agg::verify_aggregate(&aggsig1, &pm1));

    // 2 signatures
    let pm2 = vec![
        (
            hex32("1b84c5567b126440995d3ed5aaba0565d71e1834604819ff9c17f5e9d5dd078f"),
            hex32("0202020202020202020202020202020202020202020202020202020202020202").to_vec(),
        ),
        (
            hex32("462779ad4aad39514614751a71085f2f10e1c7a593e4e030efb5b8721ce55b0b"),
            hex32("0505050505050505050505050505050505050505050505050505050505050505").to_vec(),
        ),
    ];
    let aggsig2 = hex::decode("b070aafcea439a4f6f1bbfc2eb66d29d24b0cab74d6b745c3cfb009cc8fe4aa8a3afbdb45a6a34bf7c8c00f1b6d7e7d375b54540f13716c87b62e51e2f4f22ffbf8913ec53226a34892d60252a7052614ca79ae939986828d81d2311957371ad").unwrap();
    assert!(schnorr_half_agg::verify_aggregate(&aggsig2, &pm2));
}

#[test]
fn test_single_sig_aggregate_matches_individual_verify() {
    // Use key that's not d=1 (avoids P=G edge case in ecmult)
    let mut seckey = [0u8; 32];
    seckey[31] = 4; // d = 4
    let msg = [2u8; 32];
    let aux = [3u8; 32];
    let sig = schnorr::schnorr_sign(&seckey, &msg, &aux).expect("sign");
    assert!(schnorr::schnorr_verify(
        &sig,
        &msg,
        &schnorr::xonly_pubkey_from_secret(&seckey).unwrap()
    ));
    let pms = vec![(
        schnorr::xonly_pubkey_from_secret(&seckey).unwrap(),
        msg.to_vec(),
        sig,
    )];
    let aggsig = schnorr_half_agg::aggregate(&pms).expect("aggregate");
    let pm = vec![(pms[0].0, pms[0].1.clone())];
    assert!(schnorr_half_agg::verify_aggregate(&aggsig, &pm));
}

#[test]
fn test_aggregate_verify_roundtrip() {
    let mut pms: Vec<([u8; 32], Vec<u8>, [u8; 64])> = vec![];
    for i in 0..3 {
        let seckey = [i as u8 + 1; 32];
        let msg = [i as u8 + 2; 32];
        let aux = [i as u8 + 3; 32];
        let sig = schnorr::schnorr_sign(&seckey, &msg, &aux).expect("sign");
        let pk = schnorr::xonly_pubkey_from_secret(&seckey).expect("pubkey");
        pms.push((pk, msg.to_vec(), sig));
    }

    let aggsig = schnorr_half_agg::aggregate(&pms).expect("aggregate");
    let pm_tuples: Vec<([u8; 32], Vec<u8>)> =
        pms.iter().map(|(pk, m, _)| (*pk, m.clone())).collect();
    assert!(schnorr_half_agg::verify_aggregate(&aggsig, &pm_tuples));
}

#[test]
fn test_inc_aggregate() {
    let seckey0 = [1u8; 32];
    let msg0 = [2u8; 32];
    let aux0 = [3u8; 32];
    let sig0 = schnorr::schnorr_sign(&seckey0, &msg0, &aux0).unwrap();
    let pk0 = schnorr::xonly_pubkey_from_secret(&seckey0).unwrap();

    let pms0 = vec![(pk0, msg0.to_vec(), sig0)];
    let aggsig0 = schnorr_half_agg::aggregate(&pms0).unwrap();
    let pm0 = vec![(pk0, msg0.to_vec())];
    assert!(schnorr_half_agg::verify_aggregate(&aggsig0, &pm0));

    let seckey1 = [4u8; 32];
    let msg1 = [5u8; 32];
    let aux1 = [6u8; 32];
    let sig1 = schnorr::schnorr_sign(&seckey1, &msg1, &aux1).unwrap();
    let pk1 = schnorr::xonly_pubkey_from_secret(&seckey1).unwrap();

    let pms1 = vec![(pk1, msg1.to_vec(), sig1)];
    let aggsig1 = schnorr_half_agg::inc_aggregate(&aggsig0, &pm0, &pms1).unwrap();
    let pm01 = vec![(pk0, msg0.to_vec()), (pk1, msg1.to_vec())];
    assert!(schnorr_half_agg::verify_aggregate(&aggsig1, &pm01));
}
