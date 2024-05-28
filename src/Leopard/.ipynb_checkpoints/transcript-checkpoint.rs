use blstrs::*;
use merlin::Transcript;
use digest::Digest;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use sha3::Sha3_256;
use group::ff::Field;
use group::Group;

pub trait TranscriptProtocol {
    fn innerproduct_domain_sep(&mut self, n: u64);

    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar);

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar;
    
    fn validate_and_append_point(&mut self, label: &'static [u8], point: &Gt);
    
    fn append_point(&mut self, label: &'static [u8], point: &Gt);
}

impl TranscriptProtocol for Transcript{
    fn innerproduct_domain_sep(&mut self, n: u64) {
        self.append_message(b"dom-sep", b"Asiacrypt24_Cougar");
        self.append_u64(b"n", n);
    }
    
    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar) {
        self.append_message(label, &scalar.to_bytes_le());
    }
    
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);

        let mut sha3 = Sha3_256::new();
        sha3.update(b"Asiacrypt24_Cougar");
        sha3.update(buf);

        let mut rng = ChaCha20Rng::from_seed(sha3.finalize().into());
        Scalar::random(&mut rng)
    }
    
    fn validate_and_append_point(
        &mut self,
        label: &'static [u8],
        point: &Gt) -> () {
            if Gt::is_identity(point).into() {
                ()
            } else {
                self.append_point(label, point)
            }
        
    }
    
    fn append_point(&mut self, label: &'static [u8], point: &Gt) {
        let compressed = point.compress().unwrap();
        let mut buf: Vec<u8> = Vec::new();
        point.write_compressed(&mut buf).unwrap();
        self.append_message(label, &buf[..]);
    }
}