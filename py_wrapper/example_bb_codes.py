import numpy as np
import matplotlib.pyplot as plt
from qldpc import codes
from py_decoder import UFDecoder
from simulation_run import num_decoding_failures_batch


'''
union-find decoding example for bivariate bicycle (BB) codes from https://arxiv.org/pdf/2308.07915 Table 3
codes are generated using https://github.com/Infleqtion/qLDPC/tree/main
'''

# code parameters from https://arxiv.org/pdf/2308.07915 Table 3:
from sympy.abc import x, y 
l_dims = [{x: 6, y: 6}, {x: 15, y: 3}, {x: 9, y: 6}, {x: 12, y: 6}, {x: 12, y: 12}]
l_terms_a = [x**3+y+y**2, x**9+y+y**2, x**3+y+y**2, x**3+y+y**2, x**3+y**2+y**7]
l_terms_b = [y**3+x+x**2, 1+x**2+x**7, y**3+x+x**2, y**3+x+x**2, y**3+x+x**2]
n_list = [72, 90, 108, 144, 288]
k_list = [12, 8, 8, 12, 12]

# simulation parameters
num_trials = 10000
a_pauli_error_rate = np.logspace(-3, -1, num=10)
erasure_rate = 0.02

# run a simulation using UF-decoder and batch evaluation
for i, n_qbt in enumerate(n_list):
    print('n = ', n_qbt)
    bicycle_code = codes.BBCode(l_dims[i], l_terms_a[i], l_terms_b[i])
    Hz = np.array(bicycle_code.matrix_z).astype(np.uint8)
    lz = np.array(bicycle_code.get_logical_ops()[1, :, :]).astype(np.uint8)

    uf_decoder = UFDecoder(Hz)

    l_logical_error_rate = []
    for p in a_pauli_error_rate:
        num_err = num_decoding_failures_batch(uf_decoder, lz, p, erasure_rate, num_trials, topological=False)
        l_logical_error_rate.append(num_err / num_trials)
    plt.loglog(a_pauli_error_rate, l_logical_error_rate, 'o-')
    plt.loglog(a_pauli_error_rate, 1 - (1 - a_pauli_error_rate)**k_list[i], color='gray')
plt.xlabel('Pauli error rate')
plt.ylabel('logical error rate')
plt.title('erasure rate = ' + str(erasure_rate))
plt.show()

