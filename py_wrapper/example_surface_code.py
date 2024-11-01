import numpy as np
from py_decoder import UFDecoder
from some_codes import surface_code_non_periodic, plt_surface_code_non_periodic


########### decode and plot 2d surface code (non-periodic boundaries): ###########
# build H-matrix/Tanner graph
L = 40
H = surface_code_non_periodic(L)
g = UFDecoder(H)  # look at H-matrix with rough boundaries

# apply noise and get syndromes
p_err = 0.05
p_erasure = 0.1
erasure = np.random.binomial(1, p_erasure, g.n_qbt).astype(np.uint8)
pauli_err = np.random.binomial(1, p_err, g.n_qbt).astype(np.uint8)
error = np.logical_or(np.logical_and(np.logical_not(erasure), pauli_err), np.logical_and(erasure, np.random.binomial(1, 0.5, g.n_qbt))).astype(np.uint8)
syndrome = (g.h @ error % 2).astype(np.uint8)
plt_surface_code_non_periodic(L, error, g.correction, syndrome)

# decode
g.decode(syndrome, erasure)

# recompute syndrome to check decoding
err_plus_correction = np.logical_xor(error, g.correction)
syndrome = g.h @ err_plus_correction % 2
print('sum syndrome: ', np.sum(syndrome))
plt_surface_code_non_periodic(L, error, g.correction, syndrome)

