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
        return u

