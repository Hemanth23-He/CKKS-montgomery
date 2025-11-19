class MontgomeryReducer:
    def __init__(self, modulus):
        self.modulus = modulus
        self.r_bits = modulus.bit_length() + 1
        self.r = 1 << self.r_bits             # R = 2^(bit length of modulus + 1)
        self.r_inv = pow(self.r, -1, modulus) # R^-1 mod modulus
        self.np = -pow(modulus, -1, self.r)   # −N^(-1) mod R

    def to_montgomery(self, x):
        """Convert x to Montgomery domain."""
        return (x * self.r) % self.modulus

    def from_montgomery(self, x):
        """Convert x from Montgomery domain."""
        return (x * self.r_inv) % self.modulus

    def reduce(self, t):
        """Montgomery reduction: computes t * R^(-1) mod modulus."""
        m = ((t & (self.r - 1)) * self.np) & (self.r - 1)
        u = (t + m * self.modulus) >> self.r_bits
        if u >= self.modulus:
            u -= self.modulus
    def _mod_inverse(self, a, m):
    """Efficiently finds the modular inverse of a mod m using the Extended Euclidean algorithm."""
    # Ensure a and m are coprime
    g, x, y = self._extended_gcd(a, m)
    if g != 1:
        raise ValueError(f"No modular inverse for {a} mod {m}")
    return x % m

    def _extended_gcd(self, a, b):
    """Returns (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
      if a == 0:
          return (b, 0, 1)
       else:
          g, y, x = self._extended_gcd(b % a, a)
          return (g, x - (b // a) * y, y)

        return u

