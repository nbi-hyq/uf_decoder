import numpy as np
import matplotlib.pyplot as plt
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


# Constructs BB codes from https://arxiv.org/pdf/2308.07915.pdf Table 3, where A_powers and B_powers are lists.
def get_bb(l, m, A_powers, B_powers):
    A, B = get_AB(l, m, A_powers, B_powers)
    Hx = np.hstack((A, B))
    Hz = np.hstack((B.T, A.T))
    return Hx, Hz


def get_AB(l, m, A_powers, B_powers):
    A1 = x_power(l, m, A_powers[0])
    A2 = y_power(l, m, A_powers[1])
    A3 = y_power(l, m, A_powers[2])
    B1 = y_power(l, m, B_powers[0])
    B2 = x_power(l, m, B_powers[1])
    B3 = x_power(l, m, B_powers[2])
    A = (A1 + A2 + A3) % 2
    B = (B1 + B2 + B3) % 2
    return A, B


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

