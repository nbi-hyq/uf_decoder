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
    num_trials = 5000
    l_decoder = ['uf2', 'matching', 'uf']
    l_L = [10, 20, 40, 80]
    l_p = [0.006*i for i in range(20)]
    l_col1 = ['-r', '-g', '-b']
    l_col2 = [':r', ':g', ':b']
    fig1, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    for i_dec, dec in enumerate(l_decoder):
        print(dec)
        if dec == 'uf':
            DecoderClass = UnionFindPy.Decoder
        elif dec == 'matching':
            DecoderClass = pymatching.Matching
        elif dec == 'uf2':
            DecoderClass = TannerGraphDecoder
        for L in l_L:
            print(L)
            Hx = toric_code_x_stabilisers(L)
            logX = toric_code_x_logicals(L)
            l_time = []
            l_error = []
            for p in l_p:
                if dec == 'uf2':
                    total_time, num_errors = num_decoding_failures_b(DecoderClass, Hx, logX, p, num_trials)
                elif dec == 'matching':
                    total_time, num_errors = num_decoding_failures(DecoderClass, Hx, logX, p, num_trials)
                elif dec == 'uf':  # https://github.com/chaeyeunpark/UnionFind needs csr_matrix
                    total_time, num_errors = num_decoding_failures(DecoderClass, csr_matrix(Hx), logX, p, num_trials)
                l_error.append(num_errors / num_trials)
                l_time.append(total_time / num_trials)
            ax1.semilogy(l_p, l_error, l_col1[i_dec])
            ax2.semilogy(l_p, l_time, l_col2[i_dec])
    ax1.set_xlabel('physical error rate')
    ax1.set_ylabel('logical error rate')
    ax2.set_ylabel('time (s)')
    plt.savefig('error_and_time.pdf')
            