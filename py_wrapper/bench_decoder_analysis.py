import numpy as np
import matplotlib.pyplot as plt
import time
from py_decoder import UFDecoder
from some_codes import toric_code_x_logicals, toric_code_x_stabilisers


# call the decoding library only once for block of many repetitions
def num_decoding_failures_batch(decoder, logicals, p_err, p_erase, num_rep):
    # apply noise and erasures
    H = decoder.h
    n_syndr = H.shape[0]
    n_qbt = H.shape[1]
    syndrome = np.zeros(n_syndr * num_rep, dtype=np.uint8)
    erasure = np.zeros(n_qbt * num_rep, dtype=np.uint8)
    l_noise = []
    for i in range(num_rep):
        error_pauli = np.random.binomial(1, p_err, n_qbt).astype(np.uint8)
        erasure[i*n_qbt:(i+1)*n_qbt] = np.random.binomial(1, p_erase, n_qbt).astype(np.uint8)
        l_noise.append(np.logical_or(np.logical_and(np.logical_not(erasure[i*n_qbt:(i+1)*n_qbt]), error_pauli), np.logical_and(erasure[i*n_qbt:(i+1)*n_qbt], np.random.binomial(1, 0.5, n_qbt))))
        syndrome[i*n_syndr:(i+1)*n_syndr] = (H @ l_noise[i] % 2).astype(np.uint8)

    # decode batch
    decoder.correction = np.zeros(n_qbt * num_rep, dtype=np.uint8)  # create space of batch
    start_time = time.perf_counter()
    decoder.decode_batch(syndrome, erasure, num_rep)
    end_time = time.perf_counter()
    time_decode = end_time - start_time

    # evaluate decoding
    n_err = 0
    for i in range(num_rep):
        error = (l_noise[i] + decoder.correction[i*n_qbt:(i+1)*n_qbt]) % 2
        if np.any(error @ logicals.T % 2):
            n_err += 1
    return time_decode, n_err


# call library for every repetition (correction is written into existing array decoder.correction)
def num_decoding_failures(decoder, logicals, p_err, p_erase, num_rep):
    H = decoder.h
    n_err = 0
    time_decode = 0.0
    for i in range(num_rep):
        error_pauli = np.random.binomial(1, p_err, n_qbt).astype(np.uint8)
        erasure = np.random.binomial(1, p_erase, n_qbt).astype(np.uint8)
        noise = np.logical_or(np.logical_and(np.logical_not(erasure), error_pauli), np.logical_and(erasure, np.random.binomial(1, 0.5, n_qbt)))
        syndrome = (H @ noise % 2).astype(np.uint8)

        start_time = time.perf_counter()
        decoder.decode(syndrome, erasure)
        end_time = time.perf_counter()
        time_decode += (end_time - start_time)

        error = (noise + decoder.correction) % 2
        if np.any(error @ logicals.T % 2):
            n_err += 1
    return time_decode, n_err


if __name__ == "__main__":
    batch_evaluate = True  # do all repetitions of the same parameter in one block

    # speed scaling with size
    print('------ speed benchmark -------')
    num_trials = 100
    l_L = [10*(i+1) for i in range(16)]
    l_p = [0.001, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11]  # probability of Pauli error
    fig1, ax1 = plt.subplots()
    a_time = np.zeros((len(l_L), len(l_p)), dtype=np.double)
    for i_L, L in enumerate(l_L):
        print('L =', L)
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        decoder = UFDecoder(Hx)
        for i_p, p in enumerate(l_p):
            if batch_evaluate:
                total_time, _ = num_decoding_failures_batch(decoder, logX, p, 0.0, num_trials)
            else:
                total_time, _ = num_decoding_failures(decoder, logX, p, 0.0, num_trials)
            a_time[i_L, i_p] = total_time / num_trials
    for i_p, _ in enumerate(l_p):
        ax1.loglog(2*np.array(l_L)**2, a_time[:, i_p])
    ax1.loglog(np.array([2*l_L[0]**2, 2*l_L[-1]**2]), 3*10**-8 * np.array([l_L[0]**2, l_L[-1]**2]), '--k')  # linear reference
    ax1.loglog(np.array([2*l_L[0]**2, 5000, 2*l_L[-1]**2]), np.array([4.4*10**-4*2*l_L[0]**2/5000, 4.4*10**-4, 4.4*10**-4*2*l_L[-1]**2/5000]), '--b')  # linear reference arxiv1709.06218 (p=0.01)
    ax1.set_xlabel('number of qubits (2L^2)')
    ax1.set_ylabel('time')
    plt.savefig('time_scaling.pdf')
    np.savetxt('time.txt', a_time)

    # threshold and speed (different sizes, sweep error rates)
    print('------ benchmark threshold and speed -------')
    num_trials = 1000
    l_L = [10, 20, 40]
    l_p = [0.006*i for i in range(20)]  # probability of Pauli error
    fig1, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    a_time = np.zeros((len(l_L), len(l_p)), dtype=np.double)
    a_error = np.zeros((len(l_L), len(l_p)), dtype=np.double)
    for i_L, L in enumerate(l_L):
        print('L =', L)
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        decoder = UFDecoder(Hx)
        for i_p, p in enumerate(l_p):
            if batch_evaluate:
                total_time, num_errors = num_decoding_failures_batch(decoder, logX, p, 0.0, num_trials)
            else:
                total_time, num_errors = num_decoding_failures(decoder, logX, p, 0.0, num_trials)
            a_error[i_L, i_p] = num_errors / num_trials
            a_time[i_L, i_p] = total_time / num_trials
        ax1.semilogy(l_p, a_error[i_L, :])
        ax2.semilogy(l_p, a_time[i_L, :], linestyle='dashed')
    ax1.set_xlabel('physical error rate')
    ax1.set_ylabel('logical error rate')
    ax2.set_ylabel('time (s)')
    plt.savefig('error_and_time.pdf')
    np.savetxt('time2.txt', a_time)
    np.savetxt('error2.txt', a_error)

    # threshold as a function of Pauli errors and erasures
    print('------ Pauli error vs erasure thresholds -------')
    num_trials = 500
    l_L = [20, 40]
    l_p = [0.003*i for i in range(40)]  # probability of Pauli error
    l_ers = [0.05*i for i in range(10)]  # probability of erasure
    a_error = np.zeros((len(l_L), len(l_p), len(l_ers)), dtype=np.double)
    for i_L, L in enumerate(l_L):
        print('L =', L)
        Hx = toric_code_x_stabilisers(L)
        logX = toric_code_x_logicals(L)
        decoder = UFDecoder(Hx)
        for i_p, p in enumerate(l_p):
            for i_e, e in enumerate(l_ers):
                if batch_evaluate:
                    total_time, num_errors = num_decoding_failures_batch(decoder, logX, p, e, num_trials)
                else:
                    total_time, num_errors = num_decoding_failures(decoder, logX, p, e, num_trials)
                a_error[i_L, i_p, i_e] = num_errors / num_trials
    l_threshold = []
    for i_e, e in enumerate(l_ers):
        i_thresh = np.argmax(a_error[0, :, i_e] < a_error[1, :, i_e])
        l_threshold.append(l_p[i_thresh])
        fig1, ax1 = plt.subplots()
        for i_L, _ in enumerate(l_L):
            ax1.semilogy(l_p, np.squeeze(a_error[i_L, :, i_e]))
        ax1.set_xlabel('physical error rate')
        ax1.set_ylabel('logical error rate')
        ax1.set_title('erasure: ' + str(e))
        plt.savefig('error_erasure' + str(i_e) + '.pdf')
        np.savetxt('error_erasure' + str(i_e) + '.txt', np.squeeze(a_error[:, :, i_e]))
    fig1, ax1 = plt.subplots()
    ax1.plot(l_ers, l_threshold)
    ax1.set_xlabel('erasure rate')
    ax1.set_ylabel('Pauli error rate')
    ax1.set(xlim=(0, 0.5), ylim=(0, 0.1))
    plt.savefig('error_vs_erasure.pdf')
    np.savetxt('error_vs_erasure.txt', l_threshold)

