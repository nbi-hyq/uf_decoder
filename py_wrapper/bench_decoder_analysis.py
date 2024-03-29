import numpy as np
import matplotlib.pyplot as plt
import UnionFindPy
import pymatching
import sys
import time
from run_from_py import toric_code_x_logicals, toric_code_x_stabilisers, UFDecoder
from scipy.sparse import csr_matrix


# call library only once for block of many repetitions
def num_decoding_failures_batch(decoder, logicals, p_err, p_erase, num_rep):
    # apply noise and erasures
    H = decoder.h
    nnode = H.shape[0] + H.shape[1]
    syndrome = np.zeros(nnode * num_rep, dtype=np.uint8)
    erasure = np.zeros(nnode * num_rep, dtype=np.uint8)
    l_noise = []
    for i in range(num_rep):
        error_pauli = np.random.binomial(1, p_err, H.shape[1])
        erasure[i*nnode+H.shape[0]:(i+1)*nnode] = np.random.binomial(1, p_erase, H.shape[1])
        l_noise.append(np.logical_or(error_pauli, np.logical_and(erasure[i*nnode+H.shape[0]:(i+1)*nnode], np.random.binomial(1, 0.5, H.shape[1]))))
        syndrome[i*nnode:i*nnode+H.shape[0]] = H @ l_noise[i] % 2

    # decode batch
    decoder.correction = np.zeros(nnode * num_rep, dtype=np.uint8)  # create space of batch
    start_time = time.perf_counter()
    decoder.decode_batch(syndrome, erasure, num_rep)
    end_time = time.perf_counter()
    time_decode = end_time - start_time

    # evaluate decoding
    n_err = 0
    for i in range(num_rep):
        error = (l_noise[i] + decoder.correction[i*nnode+H.shape[0]:(i+1)*nnode]) % 2
        if np.any(error @ logicals.T % 2):
            n_err += 1
    return time_decode, n_err


# call library for every repetition (correction is written into existing array decoder.correction)
def num_decoding_failures(decoder, logicals, p_err, p_erase, num_rep):
    H = decoder.h
    n_err = 0
    time_decode = 0.0
    for i in range(num_rep):
        error_pauli = np.random.binomial(1, p_err, H.shape[1])
        erasure = np.random.binomial(1, p_erase, decoder.nnode)
        noise = np.logical_or(error_pauli, np.logical_and(erasure[H.shape[0]:], np.random.binomial(1, 0.5, H.shape[1])))
        syndrome = np.zeros(decoder.nnode, dtype=np.uint8)
        syndrome[:decoder.nsyndromes] = H @ noise % 2

        start_time = time.perf_counter()
        decoder.decode(syndrome, erasure)
        end_time = time.perf_counter()
        time_decode += (end_time - start_time)

        error = (noise + decoder.correction[decoder.nsyndromes:]) % 2
        if np.any(error @ logicals.T % 2):
            n_err += 1
    return time_decode, n_err


# call library for every repetition (for other decoders that return correction as new array)
def num_decoding_failures_others(decoder, H, logicals, p_err, num_rep):
    n_err = 0
    time_decode = 0.0
    for i in range(num_rep):
        noise = np.random.binomial(1, p_err, H.shape[1])
        syndrome = H @ noise % 2

        start_time = time.perf_counter()
        correction = decoder.decode(syndrome)
        end_time = time.perf_counter()
        time_decode += (end_time - start_time)

        error = (noise + correction) % 2
        if np.any(error @ logicals.T % 2):
            n_err += 1
    return time_decode, n_err


# for pymatching: call library only once for block of many repetitions
def num_decoding_failures_pymatching_batch(decoder, H, logicals, p_err, num_rep):
    noise = (np.random.random((num_rep, H.shape[1])) < p_err).astype(np.uint8)
    shots = (noise @ H.T) % 2

    start_time = time.perf_counter()
    correction = decoder.decode_batch(shots)
    end_time = time.perf_counter()
    time_decode = end_time - start_time

    error = (noise + correction) % 2
    n_err = np.sum(np.any(error @ logicals.T % 2, axis=1))
    return time_decode, n_err


# for pymatching: call library only once for block of many repetitions (trick with faults_matrix=logicals)
def num_decoding_failures_pymatching_batch_no_correction(decoder, H, logicals, p_err, num_rep):
    noise = (np.random.random((num_rep, H.shape[1])) < p_err).astype(np.uint8)
    shots = (noise @ H.T) % 2
    actual_observables = (noise @ logicals.T) % 2
    start_time = time.perf_counter()
    predicted_observables = decoder.decode_batch(shots)
    end_time = time.perf_counter()
    time_decode = end_time - start_time
    n_err = np.sum(np.any(predicted_observables != actual_observables, axis=1))
    return time_decode, n_err


if __name__ == "__main__":
    l_col1 = ['-r', '-g', '-b']
    l_col2 = [':r', ':g', ':b']
    l_decoder = ['uf', 'matching', 'uf_reference']
    batch_evaluate = True  # do all repetitions of same parameter in one block

    # speed scaling with size
    num_trials = 300
    l_L = [10*(i+1) for i in range(16)]
    l_p = [0.01, 0.05, 0.09]  # probability of Pauli error
    fig1, ax1 = plt.subplots()
    a_time = np.zeros((len(l_decoder), len(l_L), len(l_p)), dtype=np.double)
    for i_dec, dec in enumerate(l_decoder):
        print(dec)
        for i_L, L in enumerate(l_L):
            print(L)
            Hx = toric_code_x_stabilisers(L)
            logX = toric_code_x_logicals(L)
            if dec == 'uf':
                decoder = UFDecoder(Hx)
            elif dec == 'matching':
                decoder = pymatching.Matching.from_check_matrix(Hx)
                # decoder = pymatching.Matching.from_check_matrix(Hx, faults_matrix=logicals)
            elif dec == 'uf_reference':
                decoder = UnionFindPy.Decoder(Hx)
            for i_p, p in enumerate(l_p):
                if dec == 'uf':
                    if batch_evaluate:
                        total_time, num_errors = num_decoding_failures_batch(decoder, logX, p, 0.0, num_trials)
                    else:
                        total_time, num_errors = num_decoding_failures(decoder, logX, p, 0.0, num_trials)
                elif dec == 'matching':
                    if batch_evaluate:
                        total_time, num_errors = num_decoding_failures_pymatching_batch(decoder, Hx, logX, p, num_trials)
                    else:
                        total_time, num_errors = num_decoding_failures_others(decoder, Hx, logX, p, num_trials)
                elif dec == 'uf_reference':  # https://github.com/chaeyeunpark/UnionFind needs csr_matrix
                    total_time, num_errors = num_decoding_failures_others(decoder, csr_matrix(Hx), logX, p, num_trials)
                a_time[i_dec, i_L, i_p] = total_time / num_trials
        for i_p, _ in enumerate(l_p):
            ax1.loglog(np.array(l_L)**2, a_time[i_dec, :, i_p], l_col1[i_dec])
    ax1.loglog(np.array([l_L[0]**2, l_L[-1]**2]), 10**-7 * np.array([l_L[0]**2, l_L[-1]**2]), '--k')  # linear reference
    ax1.set_xlabel('size')
    ax1.set_ylabel('time')
    plt.savefig('time_scaling.pdf')

    # comparison of threshold and speed (different sizes, different error rates)
    num_trials = 500
    l_L = [10, 20, 30, 40]
    l_p = [0.006*i for i in range(20)]  # probability of Pauli error
    fig1, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    for i_dec, dec in enumerate(l_decoder):
        print(dec)
        for L in l_L:
            print(L)
            Hx = toric_code_x_stabilisers(L)
            logX = toric_code_x_logicals(L)
            if dec == 'uf':
                decoder = UFDecoder(Hx)
            elif dec == 'matching':
                decoder = pymatching.Matching.from_check_matrix(Hx)
            elif dec == 'uf_reference':
                decoder = UnionFindPy.Decoder(Hx)
            l_time = []
            l_error = []
            for p in l_p:
                if dec == 'uf':
                    if batch_evaluate:
                        total_time, num_errors = num_decoding_failures_batch(decoder, logX, p, 0.0, num_trials)
                    else:
                        total_time, num_errors = num_decoding_failures(decoder, logX, p, 0.0, num_trials)
                elif dec == 'matching':
                    if batch_evaluate:
                        total_time, num_errors = num_decoding_failures_pymatching_batch(decoder, Hx, logX, p, num_trials)
                    else:
                        total_time, num_errors = num_decoding_failures_others(decoder, Hx, logX, p, num_trials)
                elif dec == 'uf_reference':  # https://github.com/chaeyeunpark/UnionFind needs csr_matrix
                    total_time, num_errors = num_decoding_failures_others(decoder, csr_matrix(Hx), logX, p, num_trials)
                l_error.append(num_errors / num_trials)
                l_time.append(total_time / num_trials)
            ax1.semilogy(l_p, l_error, l_col1[i_dec])
            ax2.semilogy(l_p, l_time, l_col2[i_dec])
    ax1.set_xlabel('physical error rate')
    ax1.set_ylabel('logical error rate')
    ax2.set_ylabel('time (s)')
    plt.savefig('error_and_time.pdf')
