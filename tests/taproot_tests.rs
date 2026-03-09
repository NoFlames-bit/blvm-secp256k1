//! Tests for BIP 341 Taproot utilities.

use blvm_secp256k1::ecdsa::{ge_from_compressed, ge_to_compressed, pubkey_from_secret};
use blvm_secp256k1::scalar::Scalar;
use blvm_secp256k1::schnorr::xonly_pubkey_from_secret;
use blvm_secp256k1::taproot;

#[test]
fn test_taproot_output_key() {
    // BIP 341 uses tagged hash TapTweak, not plain SHA256. We verify our impl is consistent.
    let mut sk = [0u8; 32];
    sk[31] = 1;
    let internal_key = xonly_pubkey_from_secret(&sk).unwrap();
    let merkle_root = [0u8; 32];

    let output = taproot::taproot_output_key(&internal_key, &merkle_root).unwrap();
    assert_ne!(output, internal_key);
    assert_ne!(output, [0u8; 32]);
    // Output must be valid x-coordinate (lift_x succeeds)
    assert!(blvm_secp256k1::schnorr::lift_x(&output).is_some());
}

#[test]
fn test_tap_tweak_hash_deterministic() {
    let internal = [1u8; 32];
    let root = [2u8; 32];
    let h1 = taproot::tap_tweak_hash(&internal, &root);
    let h2 = taproot::tap_tweak_hash(&internal, &root);
    assert_eq!(h1, h2);
}

#[test]
fn test_tap_leaf_hash() {
    let script = vec![0x51]; // OP_1
    let h = taproot::tap_leaf_hash(0xc0, &script);
    assert_ne!(h, [0u8; 32]);
}

#[test]
fn test_tap_branch_hash() {
    let left = [1u8; 32];
    let right = [2u8; 32];
    let h = taproot::tap_branch_hash(&left, &right);
    assert_ne!(h, [0u8; 32]);
    let h_swap = taproot::tap_branch_hash(&right, &left);
    assert_ne!(h, h_swap, "order matters for TapBranch");
}

#[test]
fn test_xonly_pubkey_tweak_add() {
    let mut sk = [4u8; 32];
    sk[31] = 1;
    let pk = xonly_pubkey_from_secret(&sk).unwrap();
    let tweak = [1u8; 32];
    let (out, parity) = taproot::xonly_pubkey_tweak_add(&pk, &tweak).unwrap();
    assert_ne!(out, pk);
    assert!(parity == 0 || parity == 1);
}

#[test]
fn test_xonly_from_point() {
    let sk = [1u8; 32];
    let pk = xonly_pubkey_from_secret(&sk).unwrap();
    let ge = blvm_secp256k1::schnorr::lift_x(&pk).unwrap();
    let (bytes, parity) = taproot::xonly_from_point(&ge);
    assert_eq!(bytes, pk);
    assert_eq!(parity, 0, "lift_x returns even y");
}

#[test]
fn test_pubkey_combine() {
    let mut sk1_bytes = [0u8; 32];
    sk1_bytes[31] = 1;
    let mut sk2_bytes = [0u8; 32];
    sk2_bytes[31] = 2;
    let mut sk1 = Scalar::zero();
    sk1.set_b32(&sk1_bytes);
    let mut sk2 = Scalar::zero();
    sk2.set_b32(&sk2_bytes);
    let pk1 = ge_to_compressed(&pubkey_from_secret(&sk1));
    let pk2 = ge_to_compressed(&pubkey_from_secret(&sk2));

    let combined = taproot::pubkey_combine(&[pk1, pk2]).unwrap();
    assert_ne!(combined, pk1);
    assert_ne!(combined, pk2);
    assert!(ge_from_compressed(&combined).is_some());
}
