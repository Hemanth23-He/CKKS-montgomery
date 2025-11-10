"""A module to perform Montgomery reduction for efficient modular arithmetic."""

class MontgomeryReducer:
    """A class to perform Montgomery reduction for modular arithmetic.
    
    Montgomery reduction is an efficient method for computing modular multiplication
    without expensive division operations.
    
    Attributes:
        modulus (int): The modulus for arithmetic operations.
        r (int): Montgomery parameter, a power of 2 greater than modulus.
        r_inv (int): Modular inverse of r modulo modulus.
        n_inv (int): Negative inverse of modulus modulo r.
    """
    
    def __init__(self, modulus):
        """Initializes MontgomeryReducer with the given modulus.
        
        Args:
            modulus (int): The modulus for Montgomery reduction.
        """
        self.modulus = modulus
        # Choose r as the smallest power of 2 greater than modulus
        self.r = 1
        bit_length = modulus.bit_length()
        self.r = 1 << ((bit_length + 63) // 64 * 64)  # Round up to nearest 64-bit boundary
        
        # Compute r_inv = r^(-1) mod modulus
        self.r_inv = self._mod_inverse(self.r % modulus, modulus)
        
        # Compute n_inv = -modulus^(-1) mod r
        mod_inv = self._mod_inverse(modulus % self.r, self.r)
        self.n_inv = (self.r - mod_inv) % self.r
    
    def _mod_inverse(self, a, m):
        """Computes modular inverse using extended Euclidean algorithm.
        
        Args:
            a (int): Value to find inverse of.
            m (int): Modulus.
            
        Returns:
            Modular inverse of a modulo m.
        """
        if m == 1:
            return 0
        
        m0, x0, x1 = m, 0, 1
        
        while a > 1:
            q = a // m
            m, a = a % m, m
            x0, x1 = x1 - q * x0, x0
        
        if x1 < 0:
            x1 += m0
        
        return x1
    
    def to_montgomery(self, x):
        """Converts a value to Montgomery form.
        
        Args:
            x (int): Value to convert.
            
        Returns:
            Value in Montgomery form: x * r mod modulus.
        """
        return (x * self.r) % self.modulus
    
    def from_montgomery(self, x):
        """Converts a value from Montgomery form to standard form.
        
        Args:
            x (int): Value in Montgomery form.
            
        Returns:
            Value in standard form: x * r^(-1) mod modulus.
        """
        return self.reduce(x)
    
    def montgomery_mul(self, a, b):
        """Multiplies two Montgomery-form numbers.
        
        Computes (a * b * r^(-1)) mod modulus efficiently.
        
        Args:
            a (int): First value in Montgomery form.
            b (int): Second value in Montgomery form.
            
        Returns:
            Product in Montgomery form: (a * b * r^(-1)) mod modulus.
        """
        t = a * b
        return self.reduce(t)
    
    def reduce(self, t):
        """Performs Montgomery reduction.
        
        Computes t * r^(-1) mod modulus efficiently.
        
        Args:
            t (int): Value to reduce.
            
        Returns:
            Reduced value: t * r^(-1) mod modulus.
        """
        # Montgomery reduction algorithm
        m = ((t % self.r) * self.n_inv) % self.r
        u = (t + m * self.modulus) // self.r
        
        if u >= self.modulus:
            u -= self.modulus
        
        return u
