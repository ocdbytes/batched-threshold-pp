use ark_ec::{pairing::Pairing, scalar_mul::ScalarMul, PrimeGroup};
use ark_poly::{domain::EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{rand::RngCore, One, UniformRand, Zero};
use rand::thread_rng;
use std::{iter, vec};

use crate::utils::lagrange_interp_eval;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct CRS<E: Pairing> {
    pub powers_of_g: Vec<E::G1Affine>,
    pub htau: E::G2,

    pub y: Vec<E::G1Affine>, // Preprocessed Toeplitz matrix to compute opening proofs at all points
}

/// Dealer sets up the CRS and secret shares sk. Assumes the shares are over (1..n) and the secret key is stored at 0
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Dealer<E: Pairing> {
    batch_size: usize,
    n: usize,
    t: usize, // t+1 parties need to agree to decrypt
    sk: E::ScalarField,
}

impl<E> Dealer<E>
where
    E: Pairing,
{
    pub fn new(batch_size: usize, n: usize, t: usize) -> Self {
        let rng = &mut thread_rng();
        Self {
            batch_size,
            n,
            t,
            sk: E::ScalarField::rand(rng),
        }
    }

    pub fn get_pk(&self) -> E::G2 {
        E::G2::generator() * self.sk
    }

    pub fn setup<R: RngCore>(&mut self, rng: &mut R) -> (CRS<E>, Vec<E::ScalarField>) {
        // Sample tau and compute its powers ==========================================================
        let tau = E::ScalarField::rand(rng);
        let powers_of_tau: Vec<E::ScalarField> =
            iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(self.batch_size)
                .collect();

        // Generators
        let g = E::G1::generator();
        let h = E::G2::generator();

        // Compute powers of g
        let powers_of_g = g.batch_mul(&powers_of_tau);

        // Compute the Toeplitz matrix preprocessing ==================================================
        let mut top_tau = powers_of_tau.clone();
        top_tau.truncate(self.batch_size);
        top_tau.reverse();
        top_tau.resize(2 * self.batch_size, E::ScalarField::zero());

        let top_domain =
            Radix2EvaluationDomain::<E::ScalarField>::new(2 * self.batch_size).unwrap();
        let top_tau = top_domain.fft(&top_tau);

        // Compute powers of top_tau
        let y = g.batch_mul(&top_tau);

        let mut sk_poly = vec![E::ScalarField::zero(); self.t + 1];
        sk_poly[0] = self.sk;
        for i in 1..self.t {
            sk_poly[i] = E::ScalarField::rand(rng);
        }

        let share_domain = (1..=self.n)
            .map(|i| E::ScalarField::from(i as u64))
            .collect::<Vec<_>>();

        let eval_domain = (0..=self.t)
            .map(|i| -E::ScalarField::from(i as u64))
            .collect::<Vec<_>>();

        let sk_shares = lagrange_interp_eval(&eval_domain, &share_domain, &sk_poly);

        // let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        // share_domain.fft_in_place(&mut sk_shares);

        let crs = CRS::<E> {
            powers_of_g,
            htau: h * tau,
            y,
        };

        (crs, sk_shares)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    type E = Bls12_381;
    type Fr = <E as Pairing>::ScalarField;
    type G2 = <E as Pairing>::G2;

    #[test]
    fn test_dealer() {
        let mut rng = ark_std::test_rng();
        let batch_size = 1 << 5;
        let n = 1 << 4;
        let t = n / 2 - 1;

        let mut dealer = Dealer::<E>::new(batch_size, n, t);
        let (crs, sk_shares) = dealer.setup(&mut rng);

        let share_domain = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
        let should_be_sk = lagrange_interp_eval(&share_domain, &vec![Fr::zero()], &sk_shares)[0];
        assert_eq!(dealer.sk, should_be_sk);

        let pk = dealer.get_pk();
        let should_be_pk = G2::generator() * should_be_sk;
        assert_eq!(pk, should_be_pk);

        let g_sk_shares = sk_shares
            .iter()
            .map(|ski| G2::generator() * ski)
            .collect::<Vec<_>>();

        let interp_pk = lagrange_interp_eval(&share_domain, &vec![Fr::zero()], &g_sk_shares)[0];
        assert_eq!(pk, interp_pk);

        assert_eq!(crs.powers_of_g.len(), batch_size);
        assert_eq!(sk_shares.len(), n);
    }
}
