import numpy as np
import matplotlib.pyplot as plt
import ctypes
import numpy as np
import time
from scipy.sparse import hstack, kron, eye, csc_matrix, block_diag

# code is partly from: https://pymatching.readthedocs.io/en/stable/toric-code-example.html

def repetition_code(n):
    """
    Parity check matrix of a repetition code with length n.
    """
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i + 1) % n)))
    data = np.ones(2 * n, dtype=np.uint8)
    return csc_matrix((data, (row_ind, col_ind)))


def toric_code_x_stabilisers(L):
    """
    Sparse check matrix for the X stabilisers of a toric code with
    lattice size L, constructed as the hypergraph product of
    two repetition codes.
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
    theorem.
    """
    H1 = csc_matrix(([1], ([0], [0])), shape=(1, L), dtype=np.uint8)
    H0 = csc_matrix(np.ones((1, L), dtype=np.uint8))
    x_logicals = block_diag([kron(H1, H0), kron(H0, H1)])
    x_logicals.data = x_logicals.data % 2
    x_logicals.eliminate_zeros()
    return x_logicals


class TannerGraph:
    def __init__(self):
        self.nnode = 0
        self.nsyndromes = 0
        self.nn = None  # adjacency list representation of the Tannder graph
        self.len_nb = None
        self.is_qbt = None
        self.num_nb_max = None  # maximum vertex degree
        self.syndrome = None
        self.erasure = None
        self.decode = None
        self.h = None  # matrix mapping from qubit to syndrome in sparse form

    def h_matrix_to_tanner_graph(self):
        self.nsyndromes = self.h.shape[0]
        for i in range(len(self.h.col)):
            c = self.h.col[i]
            r = self.h.row[i]
            self.nn[self.num_nb_max * r + self.len_nb[r]] = self.nsyndromes + c  # syndromes come 1st in indexing
            self.nn[self.num_nb_max * (self.nsyndromes + c) + self.len_nb[self.nsyndromes + c]] = r
            self.len_nb[r] += 1
            self.len_nb[self.nsyndromes + c] += 1
        self.is_qbt[self.nsyndromes:] += 1

    def build_2d_surface_code(self, size):
        self.nnode = size * size * 3  # size * size syndromes + 2 * size * size qubits
        self.num_nb_max = 4  # maximum vertex degree per node
        self.nn = np.zeros(self.nnode * self.num_nb_max, dtype=np.int32)
        self.len_nb = np.zeros(self.nnode, dtype=np.uint8)
        self.is_qbt = np.zeros(self.nnode, dtype=np.uint8)
        self.syndrome = np.zeros(self.nnode, dtype=np.uint8)
        self.erasure = np.zeros(self.nnode, dtype=np.uint8)
        self.decode = np.zeros(self.nnode, dtype=np.uint8)
        self.h = toric_code_x_stabilisers(size)
        self.h_matrix_to_tanner_graph()

    def plt_2d_surface_code(self, size, error, decode):
        for i in range(g.nsyndromes):
            if g.syndrome[i]:
                plt.plot([i%size], [np.floor(i/size)], 'or')
            else:
                plt.plot([i%size], [np.floor(i/size)], 'ok')
        for i in range(g.nsyndromes):
            if error[i] and not decode[i]:
                plt.plot([i%size], [np.floor(i/size) - 0.5], '.r')
            elif not error[i] and decode[i]:
                plt.plot([i%size], [np.floor(i/size) - 0.5], '.m')
            elif error[i] and decode[i]:
                plt.plot([i%size], [np.floor(i/size) - 0.5], '.g')
        for i in range(g.nsyndromes):
            if error[i+g.nsyndromes] and not decode[i+g.nsyndromes]:
                plt.plot([i%size + 0.5], [np.floor(i/size)], '.r')
            elif not error[i+g.nsyndromes] and decode[i+g.nsyndromes]:
                plt.plot([i%size + 0.5], [np.floor(i/size)], '.m')
            elif error[i+g.nsyndromes] and decode[i+g.nsyndromes]:
                plt.plot([i%size + 0.5], [np.floor(i/size)], '.g')
        plt.show()

if __name__ == '__main__':
    # build H-matrix/Tanner graph
    size = 10
    g = TannerGraph()
    g.build_2d_surface_code(size)

    # apply noise and get sydromes
    p_err = 0.01
    p_erasure = 0.05
    g.erasure[g.nsyndromes:] = np.random.binomial(np.ones(g.nnode - g.nsyndromes, dtype=int), p_erasure)
    pauli_err = np.random.binomial(np.ones(g.nnode - g.nsyndromes, dtype=int), p_err)
    error = np.logical_or(pauli_err, np.logical_and(g.erasure[g.nsyndromes:], np.random.binomial(np.ones(g.nnode - g.nsyndromes, dtype=int), 0.5)))
    g.syndrome[:g.nsyndromes] = g.h @ error % 2
    g.plt_2d_surface_code(size, error, g.decode[g.nsyndromes:])

    # decode
    lib_graph = ctypes.cdll.LoadLibrary('../build/libSpeedDecoder.so')
    t0 = time.time()
    num_syndrome = np.sum(g.syndrome)
    lib_graph.collect_graph_and_decode(ctypes.c_int(g.nnode), ctypes.c_int(num_syndrome), ctypes.c_uint8(g.num_nb_max),
                                       ctypes.c_void_p(g.nn.ctypes.data), ctypes.c_void_p(g.len_nb.ctypes.data),
                                       ctypes.c_void_p(g.is_qbt.ctypes.data), ctypes.c_void_p(g.syndrome.ctypes.data),
                                       ctypes.c_void_p(g.erasure.ctypes.data), ctypes.c_void_p(g.decode.ctypes.data))
    t1 = time.time()
    print('time (s): ', t1 - t0)

    # do some checks
    err_plus_decode = np.logical_xor(error, g.decode[g.nsyndromes:])
    g.syndrome[:g.nsyndromes] = g.h @ err_plus_decode % 2  # recompute syndrome to check decoding
    print('syndromes: ', g.syndrome[:g.nsyndromes])
    print('error+decode: ', err_plus_decode)
    print('sum syndrome, sum E+C: ', np.sum(g.syndrome[:g.nsyndromes]), np.sum(err_plus_decode))
    g.plt_2d_surface_code(size, error, g.decode[g.nsyndromes:])

