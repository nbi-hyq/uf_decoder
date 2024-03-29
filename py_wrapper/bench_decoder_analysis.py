import numpy as np
import matplotlib.pyplot as plt
import time
from run_from_py import toric_code_x_logicals, toric_code_x_stabilisers, UFDecoder


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


if __name__ == "__main__":
    batch_evaluate = True  # do all repetitions of same parameter in one block

    # speed scaling with size
    num_trials = 3000
    l_L = [10*(i+1) for i in range(16)]
    l_p = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11]  # probability of Pauli error
    fig1, ax1 = plt.subplots()
    a_time = np.zeros((len(l_L), len(l_p)), dtype=np.double)
    for i_L, L in enumerate(l_L):
        print(L)
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        decoder = UFDecoder(Hx)
        for i_p, p in enumerate(l_p):
            if batch_evaluate:
                total_time, num_errors = num_decoding_failures_batch(decoder, logX, p, 0.0, num_trials)
            else:
                total_time, num_errors = num_decoding_failures(decoder, logX, p, 0.0, num_trials)
            a_time[i_L, i_p] = total_time / num_trials
    for i_p, _ in enumerate(l_p):
        ax1.loglog(np.array(l_L)**2, a_time[:, i_p])
    ax1.loglog(np.array([l_L[0]**2, l_L[-1]**2]), 10**-7 * np.array([l_L[0]**2, l_L[-1]**2]), '--k')  # linear reference
    ax1.set_xlabel('size')
    ax1.set_ylabel('time')
    plt.savefig('time_scaling.pdf')
    np.savetxt('time.txt', a_time)

    # threshold and speed (different sizes, sweep error rates)
    num_trials = 5000
    l_L = [10, 20, 30, 40]
    l_p = [0.006*i for i in range(20)]  # probability of Pauli error
    fig1, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    a_time_2 = np.zeros((len(l_L), len(l_p)), dtype=np.double)
    a_error_2 = np.zeros((len(l_L), len(l_p)), dtype=np.double)
    for i_L, L in enumerate(l_L):
        print(L)
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        decoder = UFDecoder(Hx)
        for i_p, p in enumerate(l_p):
            if batch_evaluate:
                total_time, num_errors = num_decoding_failures_batch(decoder, logX, p, 0.0, num_trials)
            else:
                total_time, num_errors = num_decoding_failures(decoder, logX, p, 0.0, num_trials)
            a_error_2[i_L, i_p] = num_errors / num_trials
            a_time_2[i_L, i_p] = total_time / num_trials
        ax1.semilogy(l_p, a_error_2[i_L, :])
        ax2.semilogy(l_p, a_time_2[i_L, :], linestyle='dashed')
    ax1.set_xlabel('physical error rate')
    ax1.set_ylabel('logical error rate')
    ax2.set_ylabel('time (s)')
    plt.savefig('error_and_time.pdf')
    np.savetxt('time2.txt', a_time_2)
    np.savetxt('error2.txt', a_error_2)

