//! 10x26 field element layout for 32-bit ARM.
//!
//! Uses vendored ASM for mul/sqr; all other ops in Rust.
//! Field modulus: p = 2^256 - 2^32 - 977 (SEC2 secp256k1)

use subtle::{Choice, ConstantTimeEq};

extern "C" {
    fn blvm_secp256k1_fe_mul_inner(r: *mut u32, a: *const u32, b: *const u32);
    fn blvm_secp256k1_fe_sqr_inner(r: *mut u32, a: *const u32);
}

/// Field element in 10x26 limb representation.
/// Each limb holds 26 bits (except limb 9 which holds 22).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct FieldElement {
    pub n: [u32; 10],
}

const M: u32 = 0x3FFFFFF;
const M9: u32 = 0x03FFFFF;

impl FieldElement {
    pub fn zero() -> Self {
        Self {
            n: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        }
    }

    pub fn one() -> Self {
        Self {
            n: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        }
    }

    pub fn set_b32_limit(&mut self, bytes: &[u8; 32]) -> bool {
        self.set_b32_mod(bytes);
        !((self.n[9] == M9)
            && (self.n[8] & self.n[7] & self.n[6] & self.n[5] & self.n[4] & self.n[3] & self.n[2])
                == M
            && (self.n[1]
                .wrapping_add(0x40)
                .wrapping_add((self.n[0].wrapping_add(0x3D1) >> 26))
                > M))
    }

    pub fn set_b32_mod(&mut self, a: &[u8; 32]) {
        self.n[0] = u32::from(a[31])
            | (u32::from(a[30]) << 8)
            | (u32::from(a[29]) << 16)
            | ((u32::from(a[28]) & 0x3) << 24);
        self.n[1] = ((a[28] >> 2) & 0x3f) as u32
            | (u32::from(a[27]) << 6)
            | (u32::from(a[26]) << 14)
            | ((u32::from(a[25]) & 0xF) << 22);
        self.n[2] = ((a[25] >> 4) & 0xF) as u32
            | (u32::from(a[24]) << 4)
            | (u32::from(a[23]) << 12)
            | ((u32::from(a[22]) & 0x3F) << 20);
        self.n[3] = ((a[22] >> 6) & 0x3) as u32
            | (u32::from(a[21]) << 2)
            | (u32::from(a[20]) << 10)
            | (u32::from(a[19]) << 18);
        self.n[4] = u32::from(a[18])
            | (u32::from(a[17]) << 8)
            | (u32::from(a[16]) << 16)
            | ((u32::from(a[15]) & 0x3) << 24);
        self.n[5] = ((a[15] >> 2) & 0x3f) as u32
            | (u32::from(a[14]) << 6)
            | (u32::from(a[13]) << 14)
            | ((u32::from(a[12]) & 0xF) << 22);
        self.n[6] = ((a[12] >> 4) & 0xF) as u32
            | (u32::from(a[11]) << 4)
            | (u32::from(a[10]) << 12)
            | ((u32::from(a[9]) & 0x3F) << 20);
        self.n[7] = ((a[9] >> 6) & 0x3) as u32
            | (u32::from(a[8]) << 2)
            | (u32::from(a[7]) << 10)
            | (u32::from(a[6]) << 18);
        self.n[8] = u32::from(a[5])
            | (u32::from(a[4]) << 8)
            | (u32::from(a[3]) << 16)
            | ((u32::from(a[2]) & 0x3) << 24);
        self.n[9] = ((a[2] >> 2) & 0x3f) as u32 | (u32::from(a[1]) << 6) | (u32::from(a[0]) << 14);
    }

    pub fn get_b32(&self, r: &mut [u8; 32]) {
        let a = &self.n;
        r[0] = (a[9] >> 14) as u8;
        r[1] = (a[9] >> 6) as u8;
        r[2] = (((a[9] & 0x3F) << 2) | ((a[8] >> 24) & 0x3)) as u8;
        r[3] = (a[8] >> 16) as u8;
        r[4] = (a[8] >> 8) as u8;
        r[5] = a[8] as u8;
        r[6] = (a[7] >> 18) as u8;
        r[7] = (a[7] >> 10) as u8;
        r[8] = (a[7] >> 2) as u8;
        r[9] = (((a[7] & 0x3) << 6) | ((a[6] >> 20) & 0x3f)) as u8;
        r[10] = (a[6] >> 12) as u8;
        r[11] = (a[6] >> 4) as u8;
        r[12] = (((a[6] & 0xF) << 4) | ((a[5] >> 22) & 0xF)) as u8;
        r[13] = (a[5] >> 14) as u8;
        r[14] = (a[5] >> 6) as u8;
        r[15] = (((a[5] & 0x3f) << 2) | ((a[4] >> 24) & 0x3)) as u8;
        r[16] = (a[4] >> 16) as u8;
        r[17] = (a[4] >> 8) as u8;
        r[18] = a[4] as u8;
        r[19] = (a[3] >> 18) as u8;
        r[20] = (a[3] >> 10) as u8;
        r[21] = (a[3] >> 2) as u8;
        r[22] = (((a[3] & 0x3) << 6) | ((a[2] >> 20) & 0x3f)) as u8;
        r[23] = (a[2] >> 12) as u8;
        r[24] = (a[2] >> 4) as u8;
        r[25] = (((a[2] & 0xF) << 4) | ((a[1] >> 22) & 0xF)) as u8;
        r[26] = (a[1] >> 14) as u8;
        r[27] = (a[1] >> 6) as u8;
        r[28] = (((a[1] & 0x3f) << 2) | ((a[0] >> 24) & 0x3)) as u8;
        r[29] = (a[0] >> 16) as u8;
        r[30] = (a[0] >> 8) as u8;
        r[31] = a[0] as u8;
    }

    pub fn normalize(&mut self) {
        let mut t0 = self.n[0];
        let mut t1 = self.n[1];
        let mut t2 = self.n[2];
        let mut t3 = self.n[3];
        let mut t4 = self.n[4];
        let mut t5 = self.n[5];
        let mut t6 = self.n[6];
        let mut t7 = self.n[7];
        let mut t8 = self.n[8];
        let mut t9 = self.n[9];

        let mut x = t9 >> 22;
        t9 &= M9;

        t0 += x * 0x3D1;
        t1 += x << 6;
        t1 += t0 >> 26;
        t0 &= M;
        let mut m = t2;
        t2 += t1 >> 26;
        t1 &= M;
        m &= t3;
        t3 += t2 >> 26;
        t2 &= M;
        m &= t4;
        t4 += t3 >> 26;
        t3 &= M;
        m &= t5;
        t5 += t4 >> 26;
        t4 &= M;
        m &= t6;
        t6 += t5 >> 26;
        t5 &= M;
        m &= t7;
        t7 += t6 >> 26;
        t6 &= M;
        m &= t8;
        t8 += t7 >> 26;
        t7 &= M;
        m &= t9;
        t9 += t8 >> 26;
        t8 &= M;

        x = (t9 >> 22)
            | (((t9 == M9) as u32)
                & ((m == M) as u32)
                & ((t1
                    .wrapping_add(0x40)
                    .wrapping_add((t0.wrapping_add(0x3D1) >> 26))
                    > M) as u32));

        t0 += x * 0x3D1;
        t1 += x << 6;
        t1 += t0 >> 26;
        t0 &= M;
        t2 += t1 >> 26;
        t1 &= M;
        t3 += t2 >> 26;
        t2 &= M;
        t4 += t3 >> 26;
        t3 &= M;
        t5 += t4 >> 26;
        t4 &= M;
        t6 += t5 >> 26;
        t5 &= M;
        t7 += t6 >> 26;
        t6 &= M;
        t8 += t7 >> 26;
        t7 &= M;
        t9 += t8 >> 26;
        t8 &= M;
        t9 &= M9;

        self.n[0] = t0;
        self.n[1] = t1;
        self.n[2] = t2;
        self.n[3] = t3;
        self.n[4] = t4;
        self.n[5] = t5;
        self.n[6] = t6;
        self.n[7] = t7;
        self.n[8] = t8;
        self.n[9] = t9;
    }

    pub fn is_zero(&self) -> bool {
        (self.n[0]
            | self.n[1]
            | self.n[2]
            | self.n[3]
            | self.n[4]
            | self.n[5]
            | self.n[6]
            | self.n[7]
            | self.n[8]
            | self.n[9])
            == 0
    }

    pub fn add_assign(&mut self, a: &FieldElement) {
        for i in 0..10 {
            self.n[i] += a.n[i];
        }
    }

    pub fn negate(&mut self, a: &FieldElement, m: u32) {
        let m = m as u32;
        self.n[0] = 0x3FFFC2F * 2 * (m + 1) - a.n[0];
        self.n[1] = 0x3FFFFBF * 2 * (m + 1) - a.n[1];
        self.n[2] = M * 2 * (m + 1) - a.n[2];
        self.n[3] = M * 2 * (m + 1) - a.n[3];
        self.n[4] = M * 2 * (m + 1) - a.n[4];
        self.n[5] = M * 2 * (m + 1) - a.n[5];
        self.n[6] = M * 2 * (m + 1) - a.n[6];
        self.n[7] = M * 2 * (m + 1) - a.n[7];
        self.n[8] = M * 2 * (m + 1) - a.n[8];
        self.n[9] = M9 * 2 * (m + 1) - a.n[9];
    }

    pub fn mul(&mut self, a: &FieldElement, b: &FieldElement) {
        unsafe {
            blvm_secp256k1_fe_mul_inner(self.n.as_mut_ptr(), a.n.as_ptr(), b.n.as_ptr());
        }
    }

    pub fn sqr(&mut self, a: &FieldElement) {
        unsafe {
            blvm_secp256k1_fe_sqr_inner(self.n.as_mut_ptr(), a.n.as_ptr());
        }
    }

    pub fn inv(&mut self, a: &FieldElement) {
        let mut t = *a;
        t.normalize();
        let exp: [u64; 4] = [
            0xFFFFFFFEFFFFFC2D,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
        ];
        let mut res = FieldElement::one();
        let mut base = t;
        for i in 0..256 {
            let limb_idx = i / 64;
            let bit_idx = i % 64;
            if (exp[limb_idx] >> bit_idx) & 1 == 1 {
                let mut tmp = FieldElement::zero();
                tmp.mul(&res, &base);
                res = tmp;
            }
            if i < 255 {
                let mut tmp = FieldElement::zero();
                tmp.sqr(&base);
                base = tmp;
            }
        }
        *self = res;
    }
}

impl ConstantTimeEq for FieldElement {
    fn ct_eq(&self, other: &Self) -> Choice {
        let mut c = Choice::from(1);
        for i in 0..10 {
            c = c & self.n[i].ct_eq(&other.n[i]);
        }
        c
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl Eq for FieldElement {}
