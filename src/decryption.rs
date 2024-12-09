use ark_ec::{pairing::Pairing, PrimeGroup, VariableBaseMSM};
use ark_ff::Field;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::Zero;
use std::collections::BTreeMap;

use crate::{
    dealer::CRS,
    encryption::Ciphertext,
    utils::{hash_to_bytes, lagrange_interp_eval, open_all_values, xor},
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey<E: Pairing> {
    sk_share: E::ScalarField,
}

impl<E: Pairing> SecretKey<E> {
    pub fn new(sk_share: E::ScalarField) -> Self {
        SecretKey { sk_share }
    }

    pub fn get_pk(&self) -> E::G2 {
        E::G2::generator() * self.sk_share
    }

    /// each party in the committee computes a partial decryption
    pub fn partial_decrypt(
        &self,
        ct: &Vec<Ciphertext<E>>,
        hid: E::G1,
        pk: E::G2,
        crs: &CRS<E>,
    ) -> E::G1 {
        let batch_size = crs.powers_of_g.len();
        for i in 0..batch_size {
            ct[i].verify(crs.htau, pk);
        }

        let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size).unwrap();

        let mut fevals = vec![E::ScalarField::zero(); batch_size];
        for i in 0..batch_size {
            let tg_bytes = hash_to_bytes(ct[i].gs);
            fevals[i] = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
        }
        let fcoeffs = tx_domain.ifft(&fevals);
        let com = <E::G1 as VariableBaseMSM>::msm(&crs.powers_of_g, &fcoeffs).unwrap();
        let delta = hid - com;

        let pd = delta * self.sk_share;

        pd
    }
}

/// decrypts all the ciphertexts in a batch
pub fn decrypt_all<E: Pairing>(
    public_keys: &Vec<E::G2>,
    partial_decryptions: &BTreeMap<usize, E::G1>,
    ct: &Vec<Ciphertext<E>>,
    hid: E::G1,
    crs: &CRS<E>,
) -> Vec<[u8; 32]> {
    let batch_size = ct.len();

    let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size).unwrap();

    // compute fevals by hashing gs of the ciphertexts to get fevals
    let mut fevals = vec![E::ScalarField::zero(); batch_size];
    for i in 0..batch_size {
        let tg_bytes = hash_to_bytes(ct[i].gs);
        fevals[i] = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
    }

    let fcoeffs = tx_domain.ifft(&fevals);

    let com = <E::G1 as VariableBaseMSM>::msm(&crs.powers_of_g, &fcoeffs).unwrap();
    let delta = hid - com;

    // check that all partial_decryptions are valid
    let h_inv = -E::G2::generator();
    for (&i, &pd) in partial_decryptions.iter() {
        let should_be_zero = E::multi_miller_loop([delta, pd], [public_keys[i - 1], h_inv]);
        let should_be_zero = E::final_exponentiation(should_be_zero).unwrap();
        assert!(should_be_zero.is_zero());
    }

    // interpolate partial decryptions to get the signature
    // gather partial decryptions into a vec by iterating over the BTreeMap
    let mut evals = Vec::new();
    let mut eval_points = Vec::new();
    // Iterate over the map and collect keys and values
    for (&key, &value) in partial_decryptions.iter() {
        evals.push(value);
        eval_points.push(E::ScalarField::from(key as u64));
    }

    let sigma = lagrange_interp_eval(&eval_points, &vec![E::ScalarField::zero()], &evals)[0];

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    let pi = open_all_values::<E>(&crs.y, &fcoeffs, &tx_domain);

    // now decrypt each of the ciphertexts as m = ct(1) xor H(e(delta,ct2).e(pi,ct3)e(-sigma,ct4))
    let mut m = vec![[0u8; 32]; batch_size];
    for i in 0..batch_size {
        let mask = E::multi_miller_loop([pi[i], delta, -sigma], [ct[i].ct2, ct[i].ct3, ct[i].ct4]);
        let mask = E::final_exponentiation(mask).unwrap();

        let hmask = hash_to_bytes(mask);
        m[i] = xor(&ct[i].ct1, &hmask).as_slice().try_into().unwrap();
    }

    m
}
