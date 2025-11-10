"""A module to perform computations on ciphertexts in CKKS."""
import math
from ckks.ckks_bootstrapping_context import CKKSBootstrappingContext
from util.ciphertext import Ciphertext
from util.plaintext import Plaintext
from util.polynomial import Polynomial

class CKKSEvaluator:
    """An evaluator to perform operations on CKKS ciphertexts.
    Attributes:
        degree (int): Polynomial degree for ciphertext space.
        big_modulus (int): Large modulus for multiplications.
        scaling_factor (int): Scaling factor for encoding.
        boot_context (CKKSBootstrappingContext): Bootstrapping context.
    """
    def __init__(self, params):
        """Inits evaluator with parameters.
        Args:
            params (Parameters): Parameters including polynomial degree,
                plaintext, and ciphertext modulus.
        """
        self.degree = params.poly_degree
        self.big_modulus = params.big_modulus
        self.scaling_factor = params.scaling_factor
        self.boot_context = CKKSBootstrappingContext(params)

    def add(self, ciph1, ciph2):
        """Adds ciphertexts.
        Adds two ciphertexts and returns the result.
        Args:
            ciph1 (Ciphertext): Ciphertext to be added.
            ciph2 (Ciphertext): Ciphertext to be added.
        Returns:
            A ciphertext encrypting the sum of the two ciphertext's
            plaintexts.
        """
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

    def add_plain(self, ciph, plain):
        """Adds a ciphertext with a plaintext.
        Adds a ciphertext with a plaintext and returns the result.
        Args:
            ciph (Ciphertext): Ciphertext to be added.
            plain (Plaintext): Plaintext to be added.
        Returns:
            A ciphertext encrypting the sum of the ciphertext and plaintext's
            plaintexts.
        """
        assert isinstance(ciph, Ciphertext)
        assert isinstance(plain, Plaintext)
        assert ciph.scaling_factor == plain.scaling_factor
        modulus = ciph.modulus
        c0 = ciph.c0.add(plain.poly, modulus)
        c0 = c0.mod_small(modulus)
        return Ciphertext(c0, ciph.c1, ciph.scaling_factor, modulus)

    def subtract(self, ciph1, ciph2):
        """Subtracts ciphertexts.
        Subtracts the second ciphertext from the first and returns the result.
        Args:
            ciph1 (Ciphertext): Ciphertext to be subtracted from.
            ciph2 (Ciphertext): Ciphertext to be subtracted.
        Returns:
            A ciphertext encrypting the difference of the two ciphertext's
            plaintexts.
        """
        assert isinstance(ciph1, Ciphertext)
        assert isinstance(ciph2, Ciphertext)
        assert ciph1.scaling_factor == ciph2.scaling_factor
        assert ciph1.modulus == ciph2.modulus
        modulus = ciph1.modulus
        c0 = ciph1.c0.subtract(ciph2.c0, modulus)
        c0 = c0.mod_small(modulus)
        c1 = ciph1.c1.subtract(ciph2.c1, modulus)
        c1 = c1.mod_small(modulus)
        return Ciphertext(c0, c1, ciph1.scaling_factor, modulus)

    def multiply(self, ciph1, ciph2, relin_key):
        """Multiplies ciphertexts.
        Multiplies two ciphertexts and returns the result, which is relinearized
        using the relinearization key.
        Args:
            ciph1 (Ciphertext): Ciphertext to be multiplied.
            ciph2 (Ciphertext): Ciphertext to be multiplied.
            relin_key (RelinKey): Relinearization key used for relinearization.
        Returns:
            A ciphertext encrypting the product of the two ciphertext's
            plaintexts.
        """
        assert isinstance(ciph1, Ciphertext)
        assert isinstance(ciph2, Ciphertext)
        modulus = self.big_modulus
        scaling_factor = ciph1.scaling_factor * ciph2.scaling_factor
        # Perform polynomial multiplication
        c0 = ciph1.c0.multiply(ciph2.c0, modulus)
        c0 = c0.mod_small(modulus)
        c1 = ciph1.c0.multiply(ciph2.c1, modulus)
        c1 = c1.add(ciph1.c1.multiply(ciph2.c0, modulus), modulus)
        c1 = c1.mod_small(modulus)
        c2 = ciph1.c1.multiply(ciph2.c1, modulus)
        c2 = c2.mod_small(modulus)
        # Relinearize
        return self.relinearize(c0, c1, scaling_factor, relin_key, c2)

    def multiply_plain(self, ciph, plain):
        """Multiplies a ciphertext with a plaintext.
        Multiplies a ciphertext with a plaintext and returns the result.
        Args:
            ciph (Ciphertext): Ciphertext to be multiplied.
            plain (Plaintext): Plaintext to be multiplied.
        Returns:
            A ciphertext encrypting the product of the ciphertext and plaintext's
            plaintexts.
        """
        assert isinstance(ciph, Ciphertext)
        assert isinstance(plain, Plaintext)
        modulus = ciph.modulus
        scaling_factor = ciph.scaling_factor * plain.scaling_factor
        c0 = ciph.c0.multiply(plain.poly, modulus)
        c0 = c0.mod_small(modulus)
        c1 = ciph.c1.multiply(plain.poly, modulus)
        c1 = c1.mod_small(modulus)
        return Ciphertext(c0, c1, scaling_factor, modulus)

    def relinearize(self, c0, c1, scaling_factor, relin_key, c2):
        """Relinearizes a ciphertext.
        Relinearizes a ciphertext using the relinearization key.
        Args:
            c0 (Polynomial): First polynomial in ciphertext.
            c1 (Polynomial): Second polynomial in ciphertext.
            scaling_factor (float): Scaling factor of ciphertext.
            relin_key (RelinKey): Relinearization key used for relinearization.
            c2 (Polynomial): Third polynomial in ciphertext (before relinearization).
        Returns:
            A relinearized ciphertext encrypting the same plaintext.
        """
        modulus = self.big_modulus
        modulus_squared = modulus ** 2
        # Perform key switching
        p0 = relin_key.p0.multiply(c2, modulus_squared)
        p0 = p0.mod_small(modulus_squared)
        p1 = relin_key.p1.multiply(c2, modulus_squared)
        p1 = p1.mod_small(modulus_squared)
        # Scale down
        p0 = p0.scalar_integer_divide(modulus, modulus_squared)
        p1 = p1.scalar_integer_divide(modulus, modulus_squared)
        # Add to c0, c1
        new_c0 = c0.add(p0, modulus)
        new_c0 = new_c0.mod_small(modulus)
        new_c1 = c1.add(p1, modulus)
        new_c1 = new_c1.mod_small(modulus)
        return Ciphertext(new_c0, new_c1, scaling_factor, modulus)

    def rescale(self, ciph, divisions):
        """Rescales a ciphertext.
        Rescales a ciphertext by dividing the ciphertext modulus by some value.
        Args:
            ciph (Ciphertext): Ciphertext to be rescaled.
            divisions (int): Number to divide the scaling factor by.
        Returns:
            A rescaled ciphertext.
        """
        assert isinstance(ciph, Ciphertext)
        modulus = ciph.modulus
        new_scaling_factor = ciph.scaling_factor / divisions
        c0 = ciph.c0.scalar_integer_divide(divisions, modulus)
        c1 = ciph.c1.scalar_integer_divide(divisions, modulus)
        return Ciphertext(c0, c1, new_scaling_factor, modulus)

    def lower_modulus(self, ciph, new_modulus):
        """Lowers the ciphertext modulus.
        Lowers the ciphertext modulus to a smaller value.
        Args:
            ciph (Ciphertext): Ciphertext whose modulus is to be lowered.
            new_modulus (int): New modulus value.
        Returns:
            A ciphertext with a lower modulus.
        """
        assert isinstance(ciph, Ciphertext)
        c0 = ciph.c0.mod_small(new_modulus)
        c1 = ciph.c1.mod_small(new_modulus)
        return Ciphertext(c0, c1, ciph.scaling_factor, new_modulus)

    def switch_key(self, ciph, switch_key):
        """Performs key switching on a ciphertext.
        Args:
            ciph (Ciphertext): Ciphertext to perform key switching on.
            switch_key (PublicKey): Key switching key.
        Returns:
            A key-switched ciphertext.
        """
        assert isinstance(ciph, Ciphertext)
        modulus = self.big_modulus
        modulus_squared = modulus ** 2
        # Perform key switching
        p0 = switch_key.p0.multiply(ciph.c1, modulus_squared)
        p0 = p0.mod_small(modulus_squared)
        p1 = switch_key.p1.multiply(ciph.c1, modulus_squared)
        p1 = p1.mod_small(modulus_squared)
        # Scale down
        p0 = p0.scalar_integer_divide(modulus, modulus_squared)
        p1 = p1.scalar_integer_divide(modulus, modulus_squared)
        # Add to c0, replace c1
        new_c0 = ciph.c0.add(p0, modulus)
        new_c0 = new_c0.mod_small(modulus)
        new_c1 = p1.mod_small(modulus)
        return Ciphertext(new_c0, new_c1, ciph.scaling_factor, modulus)

    def rotate(self, ciph, rotation, rotation_key):
        """Rotates a ciphertext.
        Rotates the plaintext slots of a ciphertext by a given amount.
        Args:
            ciph (Ciphertext): Ciphertext to be rotated.
            rotation (int): Number of slots to rotate by.
            rotation_key (RotationKey): Rotation key for the given rotation.
        Returns:
            A rotated ciphertext.
        """
        assert isinstance(ciph, Ciphertext)
        # Rotate c0
        c0_rot = ciph.c0.rotate(rotation)
        # Rotate c1 and perform key switching
        c1_rot = ciph.c1.rotate(rotation)
        ciph_rot = Ciphertext(c0_rot, c1_rot, ciph.scaling_factor, ciph.modulus)
        return self.switch_key(ciph_rot, rotation_key.key)

    def conjugate(self, ciph, conj_key):
        """Conjugates a ciphertext.
        Conjugates the plaintext slots of a ciphertext.
        Args:
            ciph (Ciphertext): Ciphertext to be conjugated.
            conj_key (PublicKey): Conjugation key.
        Returns:
            A conjugated ciphertext.
        """
        assert isinstance(ciph, Ciphertext)
        # Conjugate c0
        c0_conj = ciph.c0.conjugate()
        # Conjugate c1 and perform key switching
        c1_conj = ciph.c1.conjugate()
        ciph_conj = Ciphertext(c0_conj, c1_conj, ciph.scaling_factor, ciph.modulus)
        return self.switch_key(ciph_conj, conj_key)

    def raise_modulus(self, ciph, new_modulus):
        """Raises the ciphertext modulus.
        Raises the ciphertext modulus to a larger value.
        Args:
            ciph (Ciphertext): Ciphertext whose modulus is to be raised.
            new_modulus (int): New modulus value.
        Returns:
            A ciphertext with a higher modulus.
        """
        assert isinstance(ciph, Ciphertext)
        return Ciphertext(ciph.c0, ciph.c1, ciph.scaling_factor, new_modulus)

    def bootstrap(self, ciph, rotation_keys):
        """Bootstraps a ciphertext.
        Bootstraps a ciphertext to refresh the noise and modulus.
        Args:
            ciph (Ciphertext): Ciphertext to be bootstrapped.
            rotation_keys (list): List of rotation keys for bootstrapping.
        Returns:
            A bootstrapped ciphertext.
        """
        return self.boot_context.bootstrap(ciph, self, rotation_keys)
