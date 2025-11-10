"""A module to handle polynomial arithmetic in the quotient ring Z_a[x]/f(x)."""

from util.ntt import NTTContext, FFTContext
from util.montgomery import MontgomeryReducer  # Added for Montgomery reduction

class Polynomial:
    """A polynomial in the ring R_a.
    Here, R is the quotient ring Z[x]/f(x), where f(x) = x^d + 1.
    The polynomial keeps track of the ring degree d, the coefficient
    modulus a, and the coefficients in an array.
    Attributes:
        ring_degree (int): Degree d of polynomial that determines the
            quotient ring R.
        coeffs (array): Array of coefficients of polynomial, where coeffs[i]
            is the coefficient for x^i.
    """
    def __init__(self, degree, coeffs):
        """Inits Polynomial in the ring R_a with the given coefficients.
        Args:
            degree (int): Degree of quotient polynomial for ring R_a.
            coeffs (array): Array of integers of size degree, representing
                coefficients of polynomial.
        """
        self.ring_degree = degree
        assert len(coeffs) == degree, 'Size of polynomial array %d is not equal to degree %d of ring' % (len(coeffs), degree)
        self.coeffs = coeffs

    def add(self, poly, coeff_modulus=None):
        assert isinstance(poly, Polynomial)
        poly_sum = Polynomial(self.ring_degree, [self.coeffs[i] + poly.coeffs[i] for i in range(self.ring_degree)])
        if coeff_modulus:
            reducer = MontgomeryReducer(coeff_modulus)
            poly_sum.coeffs = [reducer.reduce(c) for c in poly_sum.coeffs]
        return poly_sum

    def subtract(self, poly, coeff_modulus=None):
        assert isinstance(poly, Polynomial)
        poly_diff = Polynomial(self.ring_degree, [self.coeffs[i] - poly.coeffs[i] for i in range(self.ring_degree)])
        if coeff_modulus:
            reducer = MontgomeryReducer(coeff_modulus)
            poly_diff.coeffs = [reducer.reduce(c) for c in poly_diff.coeffs]
        return poly_diff

    def multiply(self, poly, coeff_modulus, ntt=None, crt=None):
        if crt:
            return self.multiply_crt(poly, crt)
        if ntt:
            a = ntt.ftt_fwd(self.coeffs)
            b = ntt.ftt_fwd(poly.coeffs)
            ab = [a[i] * b[i] for i in range(self.ring_degree)]
            prod = ntt.ftt_inv(ab)
            reducer = MontgomeryReducer(coeff_modulus)
            prod = [reducer.reduce(x) for x in prod]
            return Polynomial(self.ring_degree, prod)
        return self.multiply_naive(poly, coeff_modulus)

    def multiply_naive(self, poly, coeff_modulus=None):
        assert isinstance(poly, Polynomial)
        reducer = MontgomeryReducer(coeff_modulus) if coeff_modulus else None
        poly_prod = Polynomial(self.ring_degree, [0] * self.ring_degree)
        for d in range(2 * self.ring_degree - 1):
            index = d % self.ring_degree
            sign = int(d < self.ring_degree) * 2 - 1
            coeff = 0
            for i in range(self.ring_degree):
                if 0 <= d - i < self.ring_degree:
                    a = self.coeffs[i]
                    b = poly.coeffs[d - i]
                    if reducer:
                        a = reducer.to_montgomery(a)
                        b = reducer.to_montgomery(b)
                        prod = reducer.montgomery_mul(a, b)
                        prod = reducer.from_montgomery(prod)
                        coeff += prod
                    else:
                        coeff += a * b
            poly_prod.coeffs[index] += sign * coeff
            if coeff_modulus:
                poly_prod.coeffs[index] = reducer.reduce(poly_prod.coeffs[index])
        return poly_prod

    def scalar_multiply(self, scalar, coeff_modulus=None):
        if coeff_modulus:
            reducer = MontgomeryReducer(coeff_modulus)
            new_coeffs = []
            for c in self.coeffs:
                a = reducer.to_montgomery(scalar)
                b = reducer.to_montgomery(c)
                prod = reducer.montgomery_mul(a, b)
                prod = reducer.from_montgomery(prod)
                new_coeffs.append(prod)
        else:
            new_coeffs = [scalar * c for c in self.coeffs]
        return Polynomial(self.ring_degree, new_coeffs)

    def scalar_integer_divide(self, scalar, coeff_modulus=None):
        if coeff_modulus:
            reducer = MontgomeryReducer(coeff_modulus)
            new_coeffs = [(c // scalar) % coeff_modulus for c in self.coeffs]
            new_coeffs = [reducer.reduce(c) for c in new_coeffs]
        else:
            new_coeffs = [(c // scalar) for c in self.coeffs]
        return Polynomial(self.ring_degree, new_coeffs)

    # All other methods (rotate, conjugate, round, floor, mod, mod_small, etc.) remain unchanged,
    # except modular reductions should use MontgomeryReducer where applicable.
    def mod(self, coeff_modulus):
        reducer = MontgomeryReducer(coeff_modulus)
        new_coeffs = [reducer.reduce(c) for c in self.coeffs]
        return Polynomial(self.ring_degree, new_coeffs)

    def mod_small(self, coeff_modulus):
        reducer = MontgomeryReducer(coeff_modulus)
        new_coeffs = [reducer.reduce(c) for c in self.coeffs]
        new_coeffs = [c - coeff_modulus if c > coeff_modulus // 2 else c for c in new_coeffs]
        return Polynomial(self.ring_degree, new_coeffs)

    # Other routines (rotate, conjugate, round, etc.) do not involve modular reduction and can stay as they are.
