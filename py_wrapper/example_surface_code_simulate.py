import numpy as np
import matplotlib.pyplot as plt
from py_decoder import UFDecoder
from some_codes import toric_code, surface_code_non_periodic
from simulation_run import num_decoding_failures_batch


########### simulate 2d surface code (with and without periodic boundaries): ###########
l_len = [5, 10, 20]  # surface code distances
num_trials = 1000  # repetitions for averaging
p_erasure = 0.0  # erasure rates
a_pauli_error_rate = np.linspace(0.03, 0.12, num=20)  # Pauli error rates

# 2d surface code with non-periodic boundaries
for i, l in enumerate(l_len):
    H, logical = surface_code_non_periodic(l)  # create parity-check matrix and logicals
    uf_decoder = UFDecoder(H)  # setup decoder
    l_logical_error_rate = []
    for p_err in a_pauli_error_rate:
        num_err = num_decoding_failures_batch(uf_decoder, logical, p_err, p_erasure, num_trials, topological=True)
        l_logical_error_rate.append(num_err / num_trials)
    plt.semilogy(a_pauli_error_rate, l_logical_error_rate, 'ro-', alpha=0.4+0.6*((i+1)/len(l_len)))

# 2d toric code with periodic boundaries
for i, l in enumerate(l_len):
    H, logical = toric_code(l)  # create parity-check matrix and logicals
    uf_decoder = UFDecoder(H)  # setup decoder
    l_logical_error_rate = []
    for p_err in a_pauli_error_rate:
        num_err = num_decoding_failures_batch(uf_decoder, logical, p_err, p_erasure, num_trials, topological=True)
        l_logical_error_rate.append(num_err / num_trials)
    plt.semilogy(a_pauli_error_rate, l_logical_error_rate, 'bo:', alpha=0.4+0.6*((i+1)/len(l_len)))

plt.show()

