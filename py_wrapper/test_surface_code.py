import numpy as np
import unittest
from py_decoder import UFDecoder
from some_codes import surface_code_non_periodic


# test decoding of square-lattice surface code with non-periodic boundaries
class TestDecodeSurfaceCode(unittest.TestCase):
    # test that error + correction gives zeros syndromes
    def test_decoding(self, l_l, l_p_err, l_p_ers):
        for L in l_l:
            H, _ = surface_code_non_periodic(L)  # create parity-check matrix
            g = UFDecoder(H)  # set up decoder
            for p_erasure in l_p_ers:
                for p_err in l_p_err:
                    # apply noise and get syndromes
                    erasure = np.random.binomial(1, p_erasure, g.n_qbt).astype(np.uint8)
                    pauli_err = np.random.binomial(1, p_err, g.n_qbt).astype(np.uint8)
                    error = np.logical_or(np.logical_and(np.logical_not(erasure), pauli_err), np.logical_and(erasure, np.random.binomial(1, 0.5, g.n_qbt))).astype(np.uint8)
                    syndrome = (g.h @ error % 2).astype(np.uint8)

                    # decode
                    g.decode(syndrome, erasure)

                    # recompute syndrome to check decoding
                    err_plus_correction = np.logical_xor(error, g.correction)
                    syndrome = g.h @ err_plus_correction % 2
                    self.assertEqual(0, np.sum(syndrome))


if __name__ == '__main__':
    t = TestDecodeSurfaceCode()
    l_p_error = [0.005 * (i+1) for i in range(25)]
    l_p_erasure = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    l_len = [10, 20, 30, 40]
    t.test_decoding(l_len, l_p_error, l_p_erasure)

