# This code is based on PyMatching document: https://pymatching.readthedocs.io/en/stable/toric-code-example.html
# Modified for the C++ UnionFind Project <https://github.com/chaeyeunpark/UnionFind>, and https://github.com/nbi-hyq/speed_decoder

import numpy as np
import matplotlib.pyplot as plt
import UnionFindPy
import pymatching
import sys
import time
from run_from_py import toric_code_x_logicals, toric_code_x_stabilisers, TannerGraphDecoder
from scipy.sparse import csr_matrix


def num_decoding_failures(DecoderClass, H, logicals, p, num_trials):
    decoder = DecoderClass(H)
    num_errors = 0
    total_time = 0.0
    for i in range(num_trials):
        noise = np.random.binomial(1, p, 2 * L * L)
        syndrome = H @ noise % 2

        start_time = time.perf_counter()
        correction = decoder.decode(syndrome)
        end_time = time.perf_counter()
        
        total_time += (end_time - start_time)
        error = (noise + correction) % 2
        if np.any(error @ logicals.T % 2):
            num_errors += 1
    return total_time, num_errors


def num_decoding_failures_b(DecoderClass, H, logicals, p, num_trials):
    decoder = DecoderClass(H)
    num_errors = 0
    total_time = 0.0
    for i in range(num_trials):
        noise = np.random.binomial(1, p, 2 * L * L)
        syndrome = np.zeros(decoder.nnode, dtype=np.uint8)
        syndrome[:decoder.nsyndromes] = H @ noise % 2

        start_time = time.perf_counter()
        decoder.decode(syndrome)
        end_time = time.perf_counter()
        
        total_time += (end_time - start_time)
        error = (noise + decoder.correction[decoder.nsyndromes:]) % 2
        if np.any(error @ logicals.T % 2):
            num_errors += 1
    return total_time, num_errors


if __name__ == "__main__":
    DEFAULT_NUM_TRIALS = 10000

    if len(sys.argv) not in [4, 5]:
        print("Usage: {} [uf | matching] L p [num_trials]".format(sys.argv[0]))
        print(f'   Default value for num_trials is {DEFAULT_NUM_TRIALS}')

    if sys.argv[1].strip().lower() == 'uf':
        DecoderClass = UnionFindPy.Decoder
    elif sys.argv[1].strip().lower() == 'matching':
        DecoderClass = pymatching.Matching
    elif sys.argv[1].strip().lower() == 'uf2':
        DecoderClass = TannerGraphDecoder
    else:
        print("The first argument must be 'uf', 'matching', or 'uf2'")
        sys.exit(1)

    L = int(sys.argv[2])
    p = float(sys.argv[3])

    if len(sys.argv) == 5:
        num_trials = int(sys.argv[4])
    else:
        num_trials = DEFAULT_NUM_TRIALS
    
    Hx = toric_code_x_stabilisers(L)
    logX = toric_code_x_logicals(L)
    if sys.argv[1].strip().lower() == 'uf2':
        total_time, num_errors = num_decoding_failures_b(DecoderClass, Hx, logX, p, num_trials)
    elif sys.argv[1].strip().lower() == 'matching':
        total_time, num_errors = num_decoding_failures(DecoderClass, Hx, logX, p, num_trials)
    elif sys.argv[1].strip().lower() == 'uf':  # https://github.com/chaeyeunpark/UnionFind needs csr_matrix
        total_time, num_errors = num_decoding_failures(DecoderClass, csr_matrix(Hx), logX, p, num_trials)

    print('{}\t{}'.format(total_time / num_trials, num_errors / num_trials))
