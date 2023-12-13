import numpy as np
import matplotlib.pyplot as plt
import ctypes


class TannerGraph:
	def __init__(self, nnode, num_nb_max, lxy):
		self.ptr = np.zeros(nnode, dtype=np.int32)
		self.nn = np.zeros(nnode*num_nb_max, dtype=np.int32)
		self.len_nb = np.zeros(nnode, dtype=np.uint8)
		self.is_qbt = np.zeros(nnode, dtype=np.uint8)
		self.num_nb_max = num_nb_max
		self.nnode = nnode
		self.bfs_list = np.zeros(nnode, dtype=np.int32)
		self.visited = np.zeros(nnode, dtype=np.uint8)
		self.syndrome = np.zeros(nnode, dtype=np.uint8)
		self.erasure = np.zeros(nnode, dtype=np.uint8)
		self.error = np.zeros(nnode, dtype=np.uint8)
		self.parity = np.zeros(nnode, dtype=np.uint8)
		self.decode = np.zeros(nnode, dtype=np.uint8)
		self.crr_surf_x = np.zeros(lxy, dtype=np.uint8)
		self.crr_surf_y = np.zeros(lxy, dtype=np.uint8)
		self.num_parity = 0
		self.big = 0

	def build_2d_surface_code(self):
		# implement in python


# gcc -Wall -g -shared -rdynamic percolation_main.c -o percolation_main.so
def run_surface_code(size=5):
	num_nb_max = 4
	nnode = size*size*3
	ptr = np.zeros(nnode, dtype=np.int32)
	nn = np.zeros(nnode*num_nb_max, dtype=np.int32)
	len_nb = np.zeros(nnode, dtype=np.uint8)
	is_qbt = np.zeros(nnode, dtype=np.uint8)

	# do full simulation?
	lib_graph = ctypes.cdll.LoadLibrary('./graph_type.so')
	lib_graph.collect_graph_and_simulate(ctypes.c_void_p(ptr.ctypes.data), ctypes.c_void_p(nn.ctypes.data), ctypes.c_void_p(len_nb.ctypes.data), ctypes.c_void_p(is_qbt.ctypes.data), ctypes.c_int(nnode), ctypes.c_int(num_nb_max))

	# or apply_erasure_and_error, get_even_clusters_bfs, get_forest, peel_forest in separate steps?


if __name__ == '__main__':
	run_surface_code()

