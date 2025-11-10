"""A module to perform computations on ciphertexts in CKKS."""
import math
from ckks.ckks_bootstrapping_context import CKKSBootstrappingContext
from util.ciphertext import Ciphertext
from util.plaintext import Plaintext
from util.polynomial import Polynomial

# If you have any explicit modular arithmetic in *this* file, import the reducer
from util.montgomery import MontgomeryReducer

class CKKSEvaluator:
    # ... (all class and method definitions exactly as before) ...

    def __init__(self, params):
        self.degree = params.poly_degree
        self.big_modulus = params.big_modulus
        self.scaling_factor = params.scaling_factor
        self.boot_context = CKKSBootstrappingContext(params)
        # If you have any CRT context in use, remove/ignore for pure Montgomery

    # All the arithmetic, such as:
    #   add, subtract, multiply, multiply_plain, relinearize,
    #   rescale, lower_modulus, switch_key, rotate, conjugate, etc.
    # can call the updated Polynomial methods, which already use Montgomery.

    # Example: (unchanged logic; Polynomial handles the reduction)
    def add(self, ciph1, ciph2):
        assert isinstance(ciph1, Ciphertext)
        assert isinstance(ciph2, Ciphertext)
        assert ciph1.scaling_factor == ciph2.scaling_factor
        assert ciph1.modulus == ciph2.modulus
        modulus = ciph1.modulus
        c0 = ciph1.c0.add(ciph2.c0, modulus)
        c0 = c0.mod_small(modulus)
        c1 = ciph1.c1.add(ciph2.c1, modulus)
        c1 = c1.mod_small(modulus)
        return Ciphertext(c0, c1, ciph1.scaling_factor, modulus)

    # For all multiplication, addition, modular operations, use only Poly methods
    # e.g. c0 = ciph1.c0.multiply(ciph2.c0, modulus)
    # Since Poly.multiply now uses MontgomeryReducer internally, you do not need to modify this

    # If there’s any explicit manual modular reduction in this file,
    # replace usages like:
    #     x = x % modulus
    # with:
    #     reducer = MontgomeryReducer(modulus)
    #     x = reducer.reduce(x)

    # Otherwise, keep all your business logic unchanged.

# The remaining business logic (multiply, multiply_plain, relinearize, etc.)
# can remain unmodified—so long as your Poly objects use the correct arithmetic!
