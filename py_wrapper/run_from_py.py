import numpy as np
import matplotlib.pyplot as plt
import ctypes
import numpy as np
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


class TannerGraphDecoder:
    def __init__(self, h):  # h can be coo_matrix, csr_matrix, or dense array
        self.h = h
        self.nsyndromes = self.h.shape[0]
        self.nnode = self.nsyndromes + self.h.shape[1]  # number of syndromes and qubits
        if type(h) == scipy.sparse._coo.coo_matrix:
            cnt = np.zeros(self.nsyndromes, dtype=np.uint8)  # count number of qubits per syndrome
            for i in self.h.row:
                cnt[i] += 1
        elif type(h) == scipy.sparse._csr.csr_matrix:
            cnt = np.zeros(self.nsyndromes, dtype=np.uint8)  # count number of qubits per syndrome
            for row in self.h.shape[0]:
                cnt[row] = len(h.getrow(row).indices)
        elif type(h) == numpy.ndarray:
            cnt = np.sum(h, axis=1)
        else:
            print('invalid parity check matrix')
        self.num_nb_max = max(2, cnt.max())  # maximum vertex degree
        self.nn = np.zeros(self.nnode * self.num_nb_max, dtype=np.int32)
        self.len_nb = np.zeros(self.nnode, dtype=np.uint8)
        self.is_qbt = np.zeros(self.nnode, dtype=np.uint8)
        self.erasure = np.zeros(self.nnode, dtype=np.uint8)
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
        if type(self.h) == scipy.sparse._coo.coo_matrix:
            for i in range(len(self.h.col)):
                c = self.h.col[i]
                r = self.h.row[i]
                self.add_from_h_row_and_col(r, c)
        elif type(self.h) == scipy.sparse._csr.csr_matrix:
            for r in self.h.shape[0]:
                for c in self.h.getrow(r).indices:
                    self.add_from_h_row_and_col(r, c)
        elif type(self.h) == numpy.ndarray:
            for r in self.h.shape[0]:
                for c in self.h.shape[1]:
                    if self.h[r, c]:
                        self.add_from_h_row_and_col(r, c)

    def decode(self, a_syndrome):
        self.decode_lib.collect_graph_and_decode(ctypes.c_int(self.nnode), ctypes.c_uint8(self.num_nb_max),
                                           ctypes.c_void_p(self.nn.ctypes.data), ctypes.c_void_p(self.len_nb.ctypes.data),
                                           ctypes.c_void_p(self.is_qbt.ctypes.data), ctypes.c_void_p(a_syndrome.ctypes.data),
                                           ctypes.c_void_p(self.erasure.ctypes.data), ctypes.c_void_p(self.correction.ctypes.data))


if __name__ == '__main__':
    # build H-matrix/Tanner graph
    L = 40
    H = toric_code_x_stabilisers(L)
    g = TannerGraphDecoder(H)

    # apply noise and get sydromes
    p_err = 0.05
    p_erasure = 0.10
    g.erasure[g.nsyndromes:] = np.random.binomial(np.ones(g.nnode - g.nsyndromes, dtype=int), p_erasure)
    pauli_err = np.random.binomial(np.ones(g.nnode - g.nsyndromes, dtype=int), p_err)
    error = np.logical_or(pauli_err, np.logical_and(g.erasure[g.nsyndromes:], np.random.binomial(np.ones(g.nnode - g.nsyndromes, dtype=int), 0.5)))
    syndrome = np.zeros(g.nnode, dtype=np.uint8)
    syndrome[:g.nsyndromes] = g.h @ error % 2
    plt_2d_square_toric_code(L, error, g.correction[g.nsyndromes:], syndrome, g.nsyndromes)

    # decode
    t0 = time.perf_counter()
    g.decode(syndrome)
    t1 = time.perf_counter()
    print('time (s): ', t1 - t0)

    # recompute syndrome to check decoding
    err_plus_correction = np.logical_xor(error, g.correction[g.nsyndromes:])
    syndrome[:g.nsyndromes] = g.h @ err_plus_correction % 2
    print('sum syndrome: ', np.sum(syndrome[:g.nsyndromes]))
    plt_2d_square_toric_code(L, error, g.correction[g.nsyndromes:], syndrome, g.nsyndromes)
