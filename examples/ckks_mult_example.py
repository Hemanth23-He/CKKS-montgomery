# This CKKS benchmark is configured for Montgomery modular reduction.
# All polynomial and arithmetic logic uses MontgomeryReducer internally.

import time
import random
from memory_profiler import memory_usage
from ckks.ckks_decryptor import CKKSDecryptor
from ckks.ckks_encoder import CKKSEncoder
from ckks.ckks_encryptor import CKKSEncryptor
from ckks.ckks_evaluator import CKKSEvaluator
from ckks.ckks_key_generator import CKKSKeyGenerator
from ckks.ckks_parameters import CKKSParameters

def random_complex_vector(poly_degree):
    length = poly_degree // 2
    # Generate complex numbers with real/imag parts rounded to 1 decimal place
    return [complex(round(random.uniform(0, 1), 1), round(random.uniform(0, 1), 1)) for _ in range(length)]

def ckks_benchmark():
    poly_degree = 1 << 12
    ciph_modulus = 1 << 600
    big_modulus = 1 << 40
    scaling_factor = 1 << 30

    start_time = time.time()

    params = CKKSParameters(poly_degree=poly_degree,
                            ciph_modulus=ciph_modulus,
                            big_modulus=big_modulus,
                            scaling_factor=scaling_factor)
    key_generator = CKKSKeyGenerator(params)
    public_key = key_generator.public_key
    secret_key = key_generator.secret_key
    relin_key = key_generator.relin_key
    encoder = CKKSEncoder(params)
    encryptor = CKKSEncryptor(params, public_key, secret_key)
    decryptor = CKKSDecryptor(params, secret_key)
    evaluator = CKKSEvaluator(params)

    message1 = random_complex_vector(poly_degree)
    message2 = random_complex_vector(poly_degree)

    plain1 = encoder.encode(message1, scaling_factor)
    plain2 = encoder.encode(message2, scaling_factor)

    ciph1 = encryptor.encrypt(plain1)
    ciph2 = encryptor.encrypt(plain2)

    ciph_prod = evaluator.multiply(ciph1, ciph2, relin_key)

    decrypted_prod = decryptor.decrypt(ciph_prod)

    decoded_prod = encoder.decode(decrypted_prod)

    end_time = time.time()
    time_taken = end_time - start_time

    round_complex = lambda x: round(x.real, 1) + round(x.imag, 1) * 1j

    results = {
        "Length of message1": len(message1),
        "Message1": message1,
        "Length of message2": len(message2),
        "Message2": message2,
        "Time taken": time_taken,
        "plaintext1": plain1,
        "plaintext2": plain2,
        "ciphertext1": ciph1,
        "ciphertext2": ciph2,
        "ciphertext_prod": ciph_prod,
        "decrypted_prod": decrypted_prod,
        "Decoded product": [round_complex(x) for x in decoded_prod],
    }
    return results

if __name__ == '__main__':
    mem, results = memory_usage((ckks_benchmark, (), {}), retval=True, max_usage=True)
    for key, val in results.items():
        print(f"{key}: {val}")
    print(f"Max Memory Usage: {mem * 1.04858:.2f} MB")
