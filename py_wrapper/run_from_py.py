import numpy as np
import matplotlib.pyplot as plt
import ctypes
import numpy as np
from scipy.sparse import hstack, kron, eye, csc_matrix, block_diag


def repetition_code(n):
    """
    Parity check matrix of a repetition code with length n.
    """
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i+1)%n)))
    data = np.ones(2*n, dtype=np.uint8)
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
    H1 = csc_matrix(([1], ([0],[0])), shape=(1,L), dtype=np.uint8)
    H0 = csc_matrix(np.ones((1, L), dtype=np.uint8))
    x_logicals = block_diag([kron(H1, H0), kron(H0, H1)])
    x_logicals.data = x_logicals.data % 2
    x_logicals.eliminate_zeros()
    return x_logicals


class TannerGraph:
	def __init__(self, nnode, num_nb_max):
    	self.nnode = 0
		self.nn = None  # adjacency list representation of the Tannder graph
		self.len_nb = None
		self.is_qbt = None
		self.num_nb_max = None  # maximum vertex degree
		self.syndrome = None
		self.erasure = None
		self.decode = None

	def build_2d_surface_code(self, size):
    	self.nnode = size*size
		self.nn = np.zeros(nnode*num_nb_max, dtype=np.int32)
		self.len_nb = np.zeros(nnode, dtype=np.uint8)
		self.is_qbt = np.zeros(nnode, dtype=np.uint8)
		self.num_nb_max = num_nb_max  # maximum vertex degree
		self.syndrome = np.zeros(nnode, dtype=np.uint8)
		self.erasure = np.zeros(nnode, dtype=np.uint8)
		self.decode = np.zeros(nnode, dtype=np.uint8)


# gcc -Wall -g -shared -rdynamic percolation_main.c -o percolation_main.so,  M.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)) C++?
if __name__ == '__main__':
    size = 10
    num_nb_max = 4
	nnode = size*size*3
    g = TannerGraph(nnode, num_nb_max)
    g.build_2d_surface_code()

    # apply noise and get sydromes

    num_syndrome = np.sum(g.syndrome)

	lib_graph = ctypes.cdll.LoadLibrary('./percolation_main.so')
	lib_graph.collect_graph_and_decode(ctypes.c_int(nnode), ctypes.c_int(num_syndrome), ctypes.c_uint8(num_nb_max), ctypes.c_void_p(g.nn.ctypes.data), ctypes.c_void_p(g.len_nb.ctypes.data), ctypes.c_void_p(g.syndrome.ctypes.data), ctypes.c_void_p(g.erasure.ctypes.data), ctypes.c_void_p(g.decode.ctypes.data))

    print(g.decode)
    print(g.syndrome)
