//! Scalar arithmetic tests and cross-validation.

use blvm_secp256k1::scalar::Scalar;
use subtle::ConstantTimeEq;

#[test]
fn test_scalar_zero_one() {
    let z = Scalar::zero();
    assert!(z.is_zero());
    assert!(!z.is_one());

    let o = Scalar::one();
    assert!(!o.is_zero());
    assert!(o.is_one());

    let mut buf = [0u8; 32];
    o.get_b32(&mut buf);
    assert_eq!(buf[31], 1);
    assert!(buf[..31].iter().all(|&b| b == 0));
}

#[test]
fn test_scalar_set_get_roundtrip() {
    let bytes =
        hex::decode("0000000000000000000000000000000000000000000000000000000000001234").unwrap();
    let mut arr = [0u8; 32];
    arr.copy_from_slice(&bytes);

    let mut s = Scalar::zero();
    s.set_b32(&arr);
    let mut out = [0u8; 32];
    s.get_b32(&mut out);
    assert_eq!(arr, out, "roundtrip failed");
}

#[test]
fn test_scalar_add() {
    let a = Scalar::one();
    let b = Scalar::one();
    let mut sum = Scalar::zero();
    sum.add(&a, &b);

    let mut expected = Scalar::zero();
    expected.set_int(2);
    assert!(bool::from(sum.ct_eq(&expected)));
}

#[test]
fn test_scalar_mul() {
    let mut a = Scalar::zero();
    a.set_int(2);
    let mut b = Scalar::zero();
    b.set_int(3);

    let mut ab = Scalar::zero();
    ab.mul(&a, &b);

    let mut expected = Scalar::zero();
    expected.set_int(6);
    assert!(bool::from(ab.ct_eq(&expected)));
}

#[test]
fn test_scalar_negate() {
    let mut a = Scalar::zero();
    a.set_int(5);

    let mut neg = Scalar::zero();
    neg.negate(&a);

    let mut sum = Scalar::zero();
    sum.add(&a, &neg);
    assert!(sum.is_zero(), "a + (-a) should be 0");
}

#[test]
fn test_scalar_split_128() {
    let mut k = Scalar::zero();
    k.set_int(0x1234);
    k.d[2] = 0x5678;
    k.d[3] = 0x9abc;

    let mut r1 = Scalar::zero();
    let mut r2 = Scalar::zero();
    Scalar::split_128(&mut r1, &mut r2, &k);

    assert_eq!(r1.d[0], 0x1234);
    assert_eq!(r1.d[1], 0);
    assert_eq!(r1.d[2], 0);
    assert_eq!(r1.d[3], 0);

    assert_eq!(r2.d[0], 0x5678);
    assert_eq!(r2.d[1], 0x9abc);
    assert_eq!(r2.d[2], 0);
    assert_eq!(r2.d[3], 0);
}

#[test]
fn test_scalar_inv_var() {
    let one = Scalar::one();
    let mut inv1 = Scalar::zero();
    inv1.inv_var(&one);
    assert!(bool::from(inv1.ct_eq(&one)), "inv(1) = 1");

    // 3^2 = 9 mod n
    let mut three = Scalar::zero();
    three.set_int(3);
    let mut nine = Scalar::zero();
    nine.set_int(9);
    let mut three_sq = Scalar::zero();
    three_sq.mul(&three, &three);
    assert!(bool::from(three_sq.ct_eq(&nine)), "3*3 = 9 mod n");

    // n * 1 = 0 mod n (curve order)
    const N_0: u64 = 0xBFD25E8CD0364141;
    const N_1: u64 = 0xBAAEDCE6AF48A03B;
    const N_2: u64 = 0xFFFFFFFFFFFFFFFE;
    const N_3: u64 = 0xFFFFFFFFFFFFFFFF;
    let n_scalar = Scalar {
        d: [N_0, N_1, N_2, N_3],
    };
    let mut n_times_one = Scalar::zero();
    n_times_one.mul(&n_scalar, &Scalar::one());
    assert!(n_times_one.is_zero(), "n*1 = 0 mod n");

    // 3^4 = 81 mod n
    let mut exp4 = Scalar::zero();
    exp4.set_int(4);
    let mut res4 = Scalar::one();
    for i in (0..256).rev() {
        let r2 = res4;
        res4.mul(&r2, &r2);
        let limb_idx = (i / 64) as usize;
        let bit_idx = i % 64;
        if (exp4.d[limb_idx] >> bit_idx) & 1 == 1 {
            let r3 = res4;
            res4.mul(&r3, &three);
        }
    }
    let mut eighty_one = Scalar::zero();
    eighty_one.set_int(81);
    assert!(bool::from(res4.ct_eq(&eighty_one)), "3^4 = 81");

    // Manual 3^2 via square-and-multiply with exp=2 (bit 1 set)
    let exp2 = Scalar { d: [2, 0, 0, 0] };
    let mut res = Scalar::one();
    for i in (0..256).rev() {
        let r2 = res;
        res.mul(&r2, &r2);
        let limb_idx = (i / 64) as usize;
        let bit_idx = i % 64;
        if (exp2.d[limb_idx] >> bit_idx) & 1 == 1 {
            let r3 = res;
            res.mul(&r3, &three);
        }
    }
    assert!(bool::from(res.ct_eq(&nine)), "3^2 via pow = 9");

    // 2^2 = 4 via full 256-iter exponentiation (sanity check for inv loop)
    let exp_two = Scalar { d: [2, 0, 0, 0] };
    let mut res2 = Scalar::one();
    let base_two = Scalar { d: [2, 0, 0, 0] };
    for i in (0..256).rev() {
        let r2 = res2;
        res2.mul(&r2, &r2);
        let limb_idx = (i / 64) as usize;
        let bit_idx = i % 64;
        if (exp_two.d[limb_idx] >> bit_idx) & 1 == 1 {
            let r3 = res2;
            res2.mul(&r3, &base_two);
        }
    }
    let mut four = Scalar::zero();
    four.set_int(4);
    assert!(bool::from(res2.ct_eq(&four)), "2^2 = 4 via pow");

    // inv(2) = (n+1)/2 per libsecp256k1
    let inv2_hex = "7fffffffffffffffffffffffffffffff5d576e7357a4501ddfe92f46681b20a1";
    let inv2_bytes = hex::decode(inv2_hex).unwrap();
    assert_eq!(inv2_bytes.len(), 32, "inv2 hex must be 32 bytes");
    let mut buf = [0u8; 32];
    buf.copy_from_slice(&inv2_bytes);
    let mut inv2_known = Scalar::zero();
    inv2_known.set_b32(&buf);
    let mut two = Scalar::zero();
    two.set_int(2);
    let mut prod_known = Scalar::zero();
    prod_known.mul(&inv2_known, &two);
    assert!(
        bool::from(prod_known.ct_eq(&Scalar::one())),
        "known inv(2)*2 = 1"
    );

    // inv(2) via our inv_var
    let mut inv2 = Scalar::zero();
    inv2.inv_var(&two);
    assert!(
        bool::from(inv2.ct_eq(&inv2_known)),
        "inv_var(2) matches known (n+1)/2"
    );

    let mut inv3 = Scalar::zero();
    inv3.inv_var(&three);
    let mut prod = Scalar::zero();
    prod.mul(&inv3, &three);
    assert!(bool::from(prod.ct_eq(&Scalar::one())), "inv(3)*3 = 1 mod n");
}

#[test]
fn test_scalar_split_lambda_smoke() {
    let mut k = Scalar::zero();
    k.set_int(42);

    let mut r1 = Scalar::zero();
    let mut r2 = Scalar::zero();
    Scalar::split_lambda(&mut r1, &mut r2, &k);

    // r1 + lambda * r2 == k (mod n). Full verification via secp256k1 cross-check later.
    assert!(!r1.is_zero() || !r2.is_zero());
}
