"""A module to multiply polynomials using the Fast Fourier Transform (FFT), Number Theoretic Transform (NTT),
and Fermat Theoretic Transform (FTT). See https://rijndael.ece.vt.edu/schaum/pdf/papers/2013hostb.pdf.
"""
from math import log, pi, cos, sin
import util.number_theory as nbtheory
from util.bit_operations import bit_reverse_vec, reverse_bits
from util.montgomery import MontgomeryReducer

class NTTContext:
    """An instance of Number/Fermat Theoretic Transform parameters.
    Here, R is the quotient ring Z_a[x]/f(x), where f(x) = x^d + 1.
    The NTTContext keeps track of the ring degree d, the coefficient
    modulus a, a root of unity w so that w^2d = 1 (mod a), and
    precomputations to perform the NTT/FTT and the inverse NTT/FTT.
    Attributes:
        coeff_modulus (int): Modulus for coefficients of the polynomial.
        degree (int): Degree of the polynomial ring.
        roots_of_unity (list): The ith member of the list is w^i, where w
            is a root of unity.
        roots_of_unity_inv (list): The ith member of the list is w^(-i),
            where w is a root of unity.
        scaled_rou_inv (list): The ith member of the list is 1/n * w^(-i),
            where w is a root of unity.
        reversed_bits (list): The ith member of the list is the bits of i
            reversed, used in the iterative implementation of NTT.
    """
    def __init__(self, poly_degree, coeff_modulus, root_of_unity=None):
        assert (poly_degree & (poly_degree - 1)) == 0, \
            "Polynomial degree must be a power of 2. d = " + str(poly_degree) + " is not."
        self.coeff_modulus = coeff_modulus
        self.degree = poly_degree
        if not root_of_unity:
            root_of_unity = nbtheory.root_of_unity(order=2 * poly_degree, modulus=coeff_modulus)
        self.precompute_ntt(root_of_unity)

    def precompute_ntt(self, root_of_unity):
        self.roots_of_unity = [1] * self.degree
        for i in range(1, self.degree):
            self.roots_of_unity[i] = (self.roots_of_unity[i - 1] * root_of_unity) % self.coeff_modulus
        root_of_unity_inv = nbtheory.mod_inv(root_of_unity, self.coeff_modulus)
        self.roots_of_unity_inv = [1] * self.degree
        for i in range(1, self.degree):
            self.roots_of_unity_inv[i] = (self.roots_of_unity_inv[i - 1] * root_of_unity_inv) % self.coeff_modulus
        self.reversed_bits = [0] * self.degree
        width = int(log(self.degree, 2))
        for i in range(self.degree):
            self.reversed_bits[i] = reverse_bits(i, width) % self.degree

    def ntt(self, coeffs, rou):
        num_coeffs = len(coeffs)
        assert len(rou) == num_coeffs, \
            "Length of the roots of unity is too small. Length is {}".format(len(rou))
        reducer = MontgomeryReducer(self.coeff_modulus)
        result = bit_reverse_vec(coeffs)
        log_num_coeffs = int(log(num_coeffs, 2))
        for logm in range(1, log_num_coeffs + 1):
            for j in range(0, num_coeffs, (1 << logm)):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))
                    rou_idx = (i << (1 + log_num_coeffs - logm))
                    omega = rou[rou_idx]
                    # Montgomery multiply for omega * result[index_odd]
                    omega_mont = reducer.to_montgomery(omega)
                    odd_mont = reducer.to_montgomery(result[index_odd])
                    omega_factor = reducer.montgomery_mul(omega_mont, odd_mont)
                    omega_factor = reducer.from_montgomery(omega_factor)
                    butterfly_plus = (result[index_even] + omega_factor) % self.coeff_modulus
                    butterfly_minus = (result[index_even] - omega_factor) % self.coeff_modulus
                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus
        return result

    def ftt_fwd(self, coeffs):
        num_coeffs = len(coeffs)
        assert num_coeffs == self.degree, "ftt_fwd: input length does not match context degree"
        reducer = MontgomeryReducer(self.coeff_modulus)
        # FTT input multiplied in Montgomery domain
        ftt_input = [reducer.montgomery_mul(reducer.to_montgomery(int(coeffs[i])), 
                                            reducer.to_montgomery(self.roots_of_unity[i]))
                     for i in range(num_coeffs)]
        # Convert back from Montgomery domain
        ftt_input = [reducer.from_montgomery(x) for x in ftt_input]
        return self.ntt(coeffs=ftt_input, rou=self.roots_of_unity)

    def ftt_inv(self, coeffs):
        num_coeffs = len(coeffs)
        assert num_coeffs == self.degree, "ntt_inv: input length does not match context degree"
        reducer = MontgomeryReducer(self.coeff_modulus)
        to_scale_down = self.ntt(coeffs=coeffs, rou=self.roots_of_unity_inv)
        poly_degree_inv = nbtheory.mod_inv(self.degree, self.coeff_modulus)
        result = [
            reducer.montgomery_mul(
                reducer.montgomery_mul(
                    reducer.to_montgomery(int(to_scale_down[i])),
                    reducer.to_montgomery(self.roots_of_unity_inv[i])
                ),
                reducer.to_montgomery(poly_degree_inv)
            ) % self.coeff_modulus for i in range(num_coeffs)
        ]
        # Convert back from Montgomery domain
        result = [reducer.from_montgomery(x) for x in result]
        return result
