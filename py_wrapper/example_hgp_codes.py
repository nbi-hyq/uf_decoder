import numpy as np
import matplotlib.pyplot as plt
from qldpc import codes
from py_decoder import UFDecoder
from simulation_run import num_decoding_failures_batch


'''
union-find decoding examples for hypergraph product codes
codes are generated using https://github.com/Infleqtion/qLDPC/tree/main
'''

##### create 2d toric codes via hypergraph product construction #####
num_trials = 1000
a_pauli_error_rate = [0.006*i for i in range(20)]
erasure_rate = 0.00

d_code = [10, 20]
l_ring_code = [codes.RingCode(d) for d in d_code]
for i, ring_code in enumerate(l_ring_code):
    toric_code = codes.HGPCode(ring_code, ring_code)
    Hz = np.array(toric_code.matrix_z).astype(np.uint8)
    lAll = toric_code.get_logical_ops()
    lz = np.array(lAll[np.logical_not(np.any(lAll[:, :Hz.shape[1]], axis=1)), Hz.shape[1]:]).astype(np.uint8)  # only logicals of Z-type (CSS)
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
    lAll = hgp_code.get_logical_ops()
    lz = np.array(lAll[np.logical_not(np.any(lAll[:, :Hz.shape[1]], axis=1)), Hz.shape[1]:]).astype(np.uint8)  # only logicals of Z-type (CSS)
    uf_decoder = UFDecoder(Hz)
    l_logical_error_rate = []
    for p in a_pauli_error_rate:
        num_err = num_decoding_failures_batch(uf_decoder, lz, p, erasure_rate, num_trials, topological=False)
        l_logical_error_rate.append(num_err / num_trials)
    plt.semilogy(a_pauli_error_rate, l_logical_error_rate, 'o-')
plt.xlabel('Pauli error rate')
plt.ylabel('logical error rate')
plt.show()

