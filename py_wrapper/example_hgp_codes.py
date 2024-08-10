import numpy as np
import matplotlib.pyplot as plt
from qldpc import codes
from py_decoder import UFDecoder

'''
union-find decoding examples for hypergraph product codes
codes are generated using https://github.com/Infleqtion/qLDPC/tree/main
'''


def num_decoding_failures_batch(decoder, logicals, p_err, p_erase, num_rep, topological=False):
    # apply random Pauli errors and erasures
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
    decoder.correction = np.zeros(n_qbt * num_rep, dtype=np.uint8)  # create space for batch of decodings
    if topological:
        decoder.decode_batch(syndrome, erasure, num_rep)  # use faster decoder for topological codes
    else:
        decoder.ldpc_decode_batch(syndrome, erasure, num_rep)  # use more general, yet slower decoder for non-topological codes

    # evaluate decoding
    n_err = 0
    for i in range(num_rep):       
        error = (l_noise[i] + decoder.correction[i*n_qbt:(i+1)*n_qbt]) % 2
        if (H @ error % 2).any():
            print('decoding invalid')
        if np.any(error @ logicals.T % 2):
            n_err += 1
    return n_err


##### create 2d toric codes via hypergraph product construction #####
num_trials = 1000
a_pauli_error_rate = [0.006*i for i in range(20)]
erasure_rate = 0.00

l_ring_code = [codes.RingCode(10), codes.RingCode(20)]
for ring_code in l_ring_code:
    toric_code = codes.HGPCode(ring_code, ring_code)
    Hz = np.array(toric_code.matrix_z).astype(np.uint8)
    lz = np.array(toric_code.get_logical_ops()[1, :, :]).astype(np.uint8)
    uf_decoder = UFDecoder(Hz)
    l_logical_error_rate = []
    for p in a_pauli_error_rate:
        num_err = num_decoding_failures_batch(uf_decoder, lz, p, erasure_rate, num_trials, topological=True)
        l_logical_error_rate.append(num_err / num_trials)
    plt.semilogy(a_pauli_error_rate, l_logical_error_rate, 'o-')
plt.xlabel('Pauli error rate')
plt.ylabel('logical error rate')
plt.show()

##### simulate non-topological hypergraph product code from two classical codes #####
num_trials = 1000
a_pauli_error_rate = [0.006*i for i in range(30)]
erasure_rate = 0.02

l_bch_code = [codes.BCHCode(7, 4), codes.BCHCode(7, 4), codes.BCHCode(7, 4)]
l_ring_code = [codes.RingCode(5), codes.RingCode(10), codes.RingCode(15)]
for i, _ in enumerate(l_bch_code):
    hgp_code = codes.HGPCode(l_bch_code[i], l_ring_code[i])
    Hz = np.array(hgp_code.matrix_z).astype(np.uint8)
    lz = np.array(hgp_code.get_logical_ops()[1, :, :]).astype(np.uint8)
    uf_decoder = UFDecoder(Hz)
    l_logical_error_rate = []
    for p in a_pauli_error_rate:
        num_err = num_decoding_failures_batch(uf_decoder, lz, p, erasure_rate, num_trials, topological=False)
        l_logical_error_rate.append(num_err / num_trials)
    plt.semilogy(a_pauli_error_rate, l_logical_error_rate, 'o-')
plt.xlabel('Pauli error rate')
plt.ylabel('logical error rate')
plt.show()

