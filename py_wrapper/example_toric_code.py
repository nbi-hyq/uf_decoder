import numpy as np
from py_decoder import UFDecoder
from some_codes import toric_code_x_stabilisers, plt_2d_square_toric_code


########### decode and plot 2d toric code (periodic boundaries): ###########
# build H-matrix/Tanner graph
L = 40
H = toric_code_x_stabilisers(L)
g = UFDecoder(H)

# apply noise and get syndromes
p_err = 0.05
p_erasure = 0.10
erasure = np.random.binomial(1, p_erasure, g.n_qbt).astype(np.uint8)
pauli_err = np.random.binomial(1, p_err, g.n_qbt).astype(np.uint8)
error = np.logical_or(np.logical_and(np.logical_not(erasure), pauli_err), np.logical_and(erasure, np.random.binomial(1, 0.5, g.n_qbt))).astype(np.uint8)
syndrome = (g.h @ error % 2).astype(np.uint8)
plt_2d_square_toric_code(L, error, g.correction, syndrome, g.n_syndr)

# decode
g.decode(syndrome, erasure)

# recompute syndrome to check decoding
err_plus_correction = np.logical_xor(error, g.correction)
syndrome = g.h @ err_plus_correction % 2
print('sum syndrome: ', np.sum(syndrome))
plt_2d_square_toric_code(L, error, g.correction, syndrome, g.n_syndr)

