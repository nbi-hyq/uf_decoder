import numpy as np
from scipy.sparse import hstack, kron, eye, csc_matrix, block_diag, csr_matrix


def repetition_code(n):
    """
    Parity check matrix of a repetition code with length n.
    from: https://pymatching.readthedocs.io/en/stable/toric-code-example.html
    """
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i + 1) % n)))
    data = np.ones(2 * n, dtype=np.uint8)
    return csc_matrix((data, (row_ind, col_ind)))


def toric_code(L):
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

    return H, x_logicals


def plt_2d_square_toric_code(size, error, correction, syndrome, nsyndromes):
    import matplotlib.pyplot as plt
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
    plt.title('or = non-zero syndrome | .r = error | .m = correction | .g = error+correction')
    plt.show()


# square lattice surface code with rough boundaries at x=0 and x=L-1, L: code distance (same for X- and Z-type errors)
# (primal and dual have the same parity check matrix, but rotated by 90°)
def surface_code_non_periodic(L):
    H = np.zeros((L*(L-1), 2*L*L*2 - 2*L + 1), dtype=np.uint8)  # primal or dual
    logical = np.zeros(2*L*L*2 - 2*L + 1, dtype=np.uint8)
    for y in range(L):
        logical[y*L] = 1
    offset = L**2  # offset of vertical edges
    for y in range(L):
        for x in range(L-1):
            H[y*(L-1)+x, y*L + x] = 1
            H[y*(L-1)+x, y*L + x + 1] = 1
            if y < L-1:
                H[y*(L-1)+x, offset + y*(L-1) + x] = 1
            if y > 0:
                H[y*(L-1)+x, offset + (y-1)*(L-1) + x] = 1
    return H, logical


# plot square lattice surface code with rough boundaries
def plt_surface_code_non_periodic(L, error, correction, syndrome):
    import matplotlib.pyplot as plt
    for y in range(L):
        for x in range(L-1):
            if syndrome[y*(L-1)+x]:
                plt.plot([x], [y], 'or')
            else:
                plt.plot([x], [y], 'ok')
    for y in range(L):
        for x in range(L):
            if error[y*L+x] and not correction[y*L+x]:
                plt.plot([x - 0.5], [y], '.r')
            elif not error[y*L+x] and correction[y*L+x]:
                plt.plot([x - 0.5], [y], '.m')
            elif error[y*L+x] and correction[y*L+x]:
                plt.plot([x - 0.5], [y], '.g')
    offset = L**2  # offset of vertical edges
    for y in range(L-1):
        for x in range(L-1):
            if error[offset + y*(L-1) + x] and not correction[offset + y*(L-1) + x]:
                plt.plot([x], [y + 0.5], '.r')
            elif not error[offset + y*(L-1) + x] and correction[offset + y*(L-1) + x]:
                plt.plot([x], [y + 0.5], '.m')
            elif error[offset + y*(L-1) + x] and correction[offset + y*(L-1) + x]:
                plt.plot([x], [y + 0.5], '.g')
    plt.title('or = non-zero syndrome | .r = error | .m = correction | .g = error+correction')
    plt.show()


def s_power(m, power):
    return np.roll(np.eye(m), -power, axis=0) # shifts columns the right by power (cyclic shift matrix)


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


def toric_code_2plus1d_stabilizers(L):
    H = np.zeros((L**3, 3*L**3), dtype=np.uint8)
    for x in range(L):
        for y in range(L):
            for z in range(L):
                H[x*L**2 + y*L + z, 3*(x*L**2 + y*L + z):3*(x*L**2 + y*L + z)+3] = np.array([1, 1, 1], dtype=np.uint8)
                H[x*L**2 + y*L + z, 3*(((x+1) % L)*L**2 + y*L + z)] = 1
                H[x*L**2 + y*L + z, 3*(x*L**2 + ((y+1) % L)*L + z) + 1] = 1
                H[x*L**2 + y*L + z, 3*(x*L**2 + y*L + ((z+1) % L)) + 2] = 1
    return H


def toric_code_2plus1d_logicals(L):
    log = np.zeros((3, 3*L**3), dtype=np.uint8)
    for i in range(L):
        for j in range(L):
            log[0, 3*(i*L**2 + j*L) + 2] = 1
            log[1, 3*(i*L**2 + j) + 1] = 1
            log[2, 3*(j*L + i) + 0] = 1
    return log

