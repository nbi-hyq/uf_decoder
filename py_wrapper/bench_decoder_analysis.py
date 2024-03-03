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


# call library only once for block of many repetitions
def num_decoding_failures_batch(DecoderClass, H, logicals, p, num_trials):
    decoder = DecoderClass(H)
    decoder.correction = np.zeros(decoder.nnode * num_trials, dtype=np.uint8)
    syndrome = np.zeros(decoder.nnode * num_trials, dtype=np.uint8)
    l_noise = []
    for i in range(num_trials):
        l_noise.append(np.random.binomial(1, p, H.shape[1]))
        syndrome[i*decoder.nnode:i*decoder.nnode+decoder.nsyndromes] = H @ l_noise[i] % 2

    start_time = time.perf_counter()
    decoder.decode_batch(syndrome, num_trials)
    end_time = time.perf_counter()
    total_time = end_time - start_time

    num_errors = 0
    for i in range(num_trials):
        error = (l_noise[i] + decoder.correction[i*decoder.nnode+decoder.nsyndromes:(i+1)*decoder.nnode]) % 2
        if np.any(error @ logicals.T % 2):
            num_errors += 1
    return total_time, num_errors


# call library for every repetition
def num_decoding_failures(DecoderClass, H, logicals, p, num_trials):
    decoder = DecoderClass(H)
    num_errors = 0
    total_time = 0.0
    for i in range(num_trials):
        noise = np.random.binomial(1, p, H.shape[1])
        syndrome = H @ noise % 2

        start_time = time.perf_counter()
        correction = decoder.decode(syndrome)
        end_time = time.perf_counter()
        
        total_time += (end_time - start_time)
        error = (noise + correction) % 2
        if np.any(error @ logicals.T % 2):
            num_errors += 1
    return total_time, num_errors


# call library for every repetition
def num_decoding_failures_b(DecoderClass, H, logicals, p, num_trials):
    decoder = DecoderClass(H)
    num_errors = 0
    total_time = 0.0
    for i in range(num_trials):
        noise = np.random.binomial(1, p, H.shape[1])
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
    l_col1 = ['-r', '-g', '-b']
    l_col2 = [':r', ':g', ':b']
    l_decoder = ['uf2', 'matching', 'uf']
    batch_evaluate = True  # do error sampling in C-progrm and do all repetitions of same parameter in one block (only used for uf2-decoder)

    # speed scaling with size
    num_trials = 3000
    l_L = [10*(i+1) for i in range(14)]
    l_p = [0.02, 0.04, 0.06, 0.08]
    fig1, ax1 = plt.subplots()
    a_time = np.zeros((len(l_decoder), len(l_L), len(l_p)), dtype=np.double)
    for i_dec, dec in enumerate(l_decoder):
        print(dec)
        if dec == 'uf':
            DecoderClass = UnionFindPy.Decoder
        elif dec == 'matching':
            DecoderClass = pymatching.Matching
        elif dec == 'uf2':
            DecoderClass = TannerGraphDecoder
        for i_L, L in enumerate(l_L):
            print(L)
            Hx = toric_code_x_stabilisers(L)
            logX = toric_code_x_logicals(L)
            for i_p, p in enumerate(l_p):
                if dec == 'uf2':
                    if batch_evaluate:
                        total_time, num_errors = num_decoding_failures_batch(DecoderClass, Hx, logX, p, num_trials)
                    else:
                        total_time, num_errors = num_decoding_failures_b(DecoderClass, Hx, logX, p, num_trials)
                elif dec == 'matching':
                    total_time, num_errors = num_decoding_failures(DecoderClass, Hx, logX, p, num_trials)
                elif dec == 'uf':  # https://github.com/chaeyeunpark/UnionFind needs csr_matrix
                    total_time, num_errors = num_decoding_failures(DecoderClass, csr_matrix(Hx), logX, p, num_trials)
                a_time[i_dec, i_L, i_p] = total_time / num_trials
        for i_p, _ in enumerate(l_p):
            ax1.loglog(np.array(l_L)**2, a_time[i_dec, :, i_p], l_col1[i_dec])
    ax1.loglog(np.array([l_L[0]**2, l_L[-1]**2]), 10**-7 * np.array([l_L[0]**2, l_L[-1]**2]), '--k')  # linear reference
    ax1.set_xlabel('size')
    ax1.set_ylabel('time')
    plt.savefig('time_scaling.pdf')

    # comparison of threhsold and speed (different sizes, different error rates)
    num_trials = 5000
    l_L = [10, 20, 30, 40]
    l_p = [0.006*i for i in range(20)]
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
                    if batch_evaluate:
                        total_time, num_errors = num_decoding_failures_batch(DecoderClass, Hx, logX, p, num_trials)
                    else:
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
            
