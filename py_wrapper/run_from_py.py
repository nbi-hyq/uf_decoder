import numpy as np
import matplotlib.pyplot as plt
import ctypes
import time
import scipy
from scipy.sparse import hstack, kron, eye, csc_matrix, block_diag, csr_matrix


def repetition_code(n):
    """
    Parity check matrix of a repetition code with length n.
    from: https://pymatching.readthedocs.io/en/stable/toric-code-example.html
    """
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i + 1) % n)))
    data = np.ones(2 * n, dtype=np.uint8)
    return csc_matrix((data, (row_ind, col_ind)))


def toric_code_x_stabilisers(L):
    """
    Sparse check matrix for the X stabilisers of a toric code with
    lattice size L, constructed as the hypergraph product of
    two repetition codes. from: https://pymatching.readthedocs.io/en/stable/toric-code-example.html
    """
    Hr = repetition_code(L)
    H = hstack(
        [kron(Hr, eye(Hr.shape[1])), kron(eye(Hr.shape[0]), Hr.T)],
        dtype=np.uint8
    )
    H.data = H.data % 2
    H.eliminate_zeros()
    return H


def toric_code_x_logicals(L):
    """
    Sparse binary matrix with each row corresponding to an X logical operator
    of a toric code with lattice size L. Constructed from the
    homology groups of the repetition codes using the Kunneth
    theorem. from: https://pymatching.readthedocs.io/en/stable/toric-code-example.html
    """
    H1 = csc_matrix(([1], ([0], [0])), shape=(1, L), dtype=np.uint8)
    H0 = csc_matrix(np.ones((1, L), dtype=np.uint8))
    x_logicals = block_diag([kron(H1, H0), kron(H0, H1)])
    x_logicals.data = x_logicals.data % 2
    x_logicals.eliminate_zeros()
    return x_logicals


def plt_2d_square_toric_code(size, error, correction, syndrome, nsyndromes):
    for i in range(nsyndromes):
        if syndrome[i]:
            plt.plot([i%size], [np.floor(i/size)], 'or')
        else:
            plt.plot([i%size], [np.floor(i/size)], 'ok')
    for i in range(nsyndromes):
        if error[i] and not correction[i]:
            plt.plot([i%size], [np.floor(i/size) - 0.5], '.r')
        elif not error[i] and correction[i]:
            plt.plot([i%size], [np.floor(i/size) - 0.5], '.m')
        elif error[i] and correction[i]:
            plt.plot([i%size], [np.floor(i/size) - 0.5], '.g')
    for i in range(nsyndromes):
        if error[i+nsyndromes] and not correction[i+nsyndromes]:
            plt.plot([i%size + 0.5], [np.floor(i/size)], '.r')
        elif not error[i+nsyndromes] and correction[i+nsyndromes]:
            plt.plot([i%size + 0.5], [np.floor(i/size)], '.m')
        elif error[i+nsyndromes] and correction[i+nsyndromes]:
            plt.plot([i%size + 0.5], [np.floor(i/size)], '.g')
    plt.show()


def s_power(m, power):
    return np.roll(np.eye(m), -power, axis=0) # shifts colummns the right by power (cyclic shift matrix)


def x_power(l, m, power):
    s_p = s_power(l, power)
    return np.kron(s_p, np.eye(m))


def y_power(l, m, power):
    s_p = s_power(m, power)
    return np.kron(np.eye(l), s_p)


def get_AB(l, m, A_powers, B_powers):
    """
    Constructs matrices A and B where A_powers and B_powers are lists.
    from https://arxiv.org/pdf/2308.07915.pdf Table 3 1st example
    """
    A1 = x_power(l, m, A_powers[0])
    A2 = y_power(l, m, A_powers[1])
    A3 = y_power(l, m, A_powers[2])
    B1 = y_power(l, m, B_powers[0])
    B2 = x_power(l, m, B_powers[1])
    B3 = x_power(l, m, B_powers[2])
    A = (A1 + A2 + A3) % 2
    B = (B1 + B2 + B3) % 2
    return A, B


def get_H(l, m, A_powers, B_powers):
    A, B = get_AB(l, m, A_powers, B_powers)
    Hx = np.hstack((A, B))
    Hz = np.hstack((B.T, A.T))
    return Hx, Hz


def toric_code_3d_stabilizers(L):
    H = np.zeros((L**3, 3*L**3), dtype=np.uint8)
    for x in range(L):
        for y in range(L):
            for z in range(L):
                H[x*L**2 + y*L + z, 3*(x*L**2 + y*L + z):3*(x*L**2 + y*L + z)+3] = np.array([1, 1, 1], dtype=np.uint8)
                H[x*L**2 + y*L + z, 3*(((x+1) % L)*L**2 + y*L + z)] = 1
                H[x*L**2 + y*L + z, 3*(x*L**2 + ((y+1) % L)*L + z) + 1] = 1
                H[x*L**2 + y*L + z, 3*(x*L**2 + y*L + ((z+1) % L)) + 2] = 1
    return H


def toric_code_3d_logicals(L):
    log = np.zeros((3, 3*L**3), dtype=np.uint8)
    for i in range(L):
        for j in range(L):
            log[0, 3*(i*L**2 + j*L) + 2] = 1
            log[1, 3*(i*L**2 + j) + 1] = 1
            log[2, 3*(j*L + i) + 0] = 1
    return log


class UFDecoder:
    def __init__(self, h):  # h can be scipy coo_matrix, csr_matrix, or numpy array
        self.h = h
        self.nsyndromes = self.h.shape[0]
        self.nnode = self.nsyndromes + self.h.shape[1]  # number of syndromes and qubits
        if type(h) != np.ndarray and h.getformat() == 'coo':
            cnt = np.zeros(self.nsyndromes, dtype=np.uint8)  # count number of qubits per syndrome
            for i in self.h.row:
                cnt[i] += 1
        elif type(h) != np.ndarray and h.getformat() == 'csr':
            cnt = np.zeros(self.nsyndromes, dtype=np.uint8)  # count number of qubits per syndrome
            for row in range(self.h.shape[0]):
                cnt[row] = len(h.getrow(row).indices)
        elif type(h) == np.ndarray:
            cnt = np.sum(h, axis=1, dtype=np.uint8)
        else:
            print('invalid parity check matrix')
        self.num_nb_max = max(2, cnt.max())  # maximum vertex degree
        self.nn = np.zeros(self.nnode * self.num_nb_max, dtype=np.int32)
        self.len_nb = np.zeros(self.nnode, dtype=np.uint8)
        self.is_qbt = np.zeros(self.nnode, dtype=np.uint8)
        self.correction = np.zeros(self.nnode, dtype=np.uint8)
        self.h_matrix_to_tanner_graph()
        self.decode_lib = ctypes.cdll.LoadLibrary('../build/libSpeedDecoder.so')

    def add_from_h_row_and_col(self, r, c):
        self.nn[self.num_nb_max * r + self.len_nb[r]] = self.nsyndromes + c  # syndromes come 1st in indexing
        self.nn[self.num_nb_max * (self.nsyndromes + c) + self.len_nb[self.nsyndromes + c]] = r
        self.len_nb[r] += 1
        self.len_nb[self.nsyndromes + c] += 1

    def h_matrix_to_tanner_graph(self):
        self.is_qbt[self.nsyndromes:] += 1
        if type(self.h) != np.ndarray and self.h.getformat() == 'coo':
            for i in range(len(self.h.col)):
                c = self.h.col[i]
                r = self.h.row[i]
                self.add_from_h_row_and_col(r, c)
        elif type(self.h) != np.ndarray and self.h.getformat() == 'csr':
            for r in range(self.h.shape[0]):
                for c in self.h.getrow(r).indices:
                    self.add_from_h_row_and_col(r, c)
        elif type(self.h) == np.ndarray:
            for r in range(self.h.shape[0]):
                for c in range(self.h.shape[1]):
                    if self.h[r, c]:
                        self.add_from_h_row_and_col(r, c)

    def decode(self, a_syndrome, a_erasure):
        self.decode_lib.collect_graph_and_decode(ctypes.c_int(self.nnode), ctypes.c_uint8(self.num_nb_max),
                                           ctypes.c_void_p(self.nn.ctypes.data), ctypes.c_void_p(self.len_nb.ctypes.data),
                                           ctypes.c_void_p(self.is_qbt.ctypes.data), ctypes.c_void_p(a_syndrome.ctypes.data),
                                           ctypes.c_void_p(a_erasure.ctypes.data), ctypes.c_void_p(self.correction.ctypes.data))

    def decode_batch(self, a_syndrome, a_erasure, nrep):
        self.decode_lib.collect_graph_and_decode_batch(ctypes.c_int(self.nnode), ctypes.c_uint8(self.num_nb_max),
                                           ctypes.c_void_p(self.nn.ctypes.data), ctypes.c_void_p(self.len_nb.ctypes.data),
                                           ctypes.c_void_p(self.is_qbt.ctypes.data), ctypes.c_void_p(a_syndrome.ctypes.data),
                                           ctypes.c_void_p(a_erasure.ctypes.data), ctypes.c_void_p(self.correction.ctypes.data), ctypes.c_int(nrep))

    def ldpc_decode(self, a_syndrome, a_erasure):
        self.decode_lib.ldpc_collect_graph_and_decode(ctypes.c_int(self.nnode), ctypes.c_uint8(self.num_nb_max),
                                           ctypes.c_void_p(self.nn.ctypes.data), ctypes.c_void_p(self.len_nb.ctypes.data),
                                           ctypes.c_void_p(self.is_qbt.ctypes.data), ctypes.c_void_p(a_syndrome.ctypes.data),
                                           ctypes.c_void_p(a_erasure.ctypes.data), ctypes.c_void_p(self.correction.ctypes.data))

    def ldpc_decode_batch(self, a_syndrome, a_erasure, nrep):
        self.decode_lib.ldpc_collect_graph_and_decode_batch(ctypes.c_int(self.nnode), ctypes.c_uint8(self.num_nb_max),
                                           ctypes.c_void_p(self.nn.ctypes.data), ctypes.c_void_p(self.len_nb.ctypes.data),
                                           ctypes.c_void_p(self.is_qbt.ctypes.data), ctypes.c_void_p(a_syndrome.ctypes.data),
                                           ctypes.c_void_p(a_erasure.ctypes.data), ctypes.c_void_p(self.correction.ctypes.data), ctypes.c_int(nrep))

if __name__ == '__main__':
    ########### 2d surface code: ###########
    # build H-matrix/Tanner graph
    L = 40
    H = toric_code_x_stabilisers(L)
    g = UFDecoder(H)

    # apply noise and get sydromes
    p_err = 0.05
    p_erasure = 0.10
    erasure = np.zeros(g.nnode, dtype=np.uint8)
    erasure[g.nsyndromes:] = np.random.binomial(1, p_erasure, g.nnode - g.nsyndromes)
    pauli_err = np.random.binomial(1, p_err, g.nnode - g.nsyndromes)
    error = np.logical_or(pauli_err, np.logical_and(erasure[g.nsyndromes:], np.random.binomial(1, 0.5, g.nnode - g.nsyndromes)))
    syndrome = np.zeros(g.nnode, dtype=np.uint8)
    syndrome[:g.nsyndromes] = g.h @ error % 2
    plt_2d_square_toric_code(L, error, g.correction[g.nsyndromes:], syndrome, g.nsyndromes)

    # decode
    t0 = time.perf_counter()
    g.decode(syndrome, erasure)
    t1 = time.perf_counter()
    print('time (s): ', t1 - t0)

    # recompute syndrome to check decoding
    err_plus_correction = np.logical_xor(error, g.correction[g.nsyndromes:])
    syndrome[:g.nsyndromes] = g.h @ err_plus_correction % 2
    print('sum syndrome: ', np.sum(syndrome[:g.nsyndromes]))
    plt_2d_square_toric_code(L, error, g.correction[g.nsyndromes:], syndrome, g.nsyndromes)

    ########### bivariate bicycle LDPC code: ###########
    # build H-matrix/Tanner graph
    H, _ = get_H(6, 6, [3, 1, 2], [3, 1, 2])
    g = UFDecoder(H)

    # apply noise and get sydromes
    for p_err in [0.01*i for i in range(10)]:
        for p_erasure in [0.05*i for i in range(10)]:
            erasure = np.zeros(g.nnode, dtype=np.uint8)
            erasure[g.nsyndromes:] = np.random.binomial(1, p_erasure, g.nnode - g.nsyndromes)
            pauli_err = np.random.binomial(1, p_err, g.nnode - g.nsyndromes)
            error = np.logical_or(pauli_err, np.logical_and(erasure[g.nsyndromes:], np.random.binomial(1, 0.5, g.nnode - g.nsyndromes)))
            syndrome = np.zeros(g.nnode, dtype=np.uint8)
            syndrome[:g.nsyndromes] = g.h @ error % 2

            # decode
            t0 = time.perf_counter()
            g.ldpc_decode(syndrome, erasure)
            t1 = time.perf_counter()
            print('time (s): ', t1 - t0)

            # recompute syndrome to check decoding
            err_plus_correction = np.logical_xor(error, g.correction[g.nsyndromes:])
            syndrome[:g.nsyndromes] = g.h @ err_plus_correction % 2
            print('sum syndrome: ', np.sum(syndrome[:g.nsyndromes]))
 
