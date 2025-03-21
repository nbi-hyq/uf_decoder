import numpy as np
import unittest
from some_codes import get_bb
from py_decoder import UFDecoder


'''
union-find decoding example for bivariate bicycle (BB) codes from https://arxiv.org/pdf/2308.07915 Table 3
codes are generated using https://github.com/Infleqtion/qLDPC/tree/main
'''

class TestDecodeBBCode(unittest.TestCase):
    # test that error + correction gives zeros syndromes
    def test_decoding(self, l_p_err, l_p_erasure):
        # code parameters from https://arxiv.org/pdf/2308.07915 Table 3:
        l_dims = [[6, 6], [15, 3], [9, 6], [12, 6], [12, 12]]
        pow_a = [[3, 1, 2], [9, 1, 2], [3, 1, 2], [3, 1, 2], [3, 2, 7]]
        pow_b = [[3, 1, 2], [0, 2, 7], [3, 1, 2], [3, 1, 2], [3, 1, 2]]

        for i in range(len(l_dims)):
            Hx, Hz = get_bb(l_dims[i][0], l_dims[i][1], pow_a[i], pow_b[i])  #.astype(np.uint8)
            for h_css in [Hz, Hx]:
                g = UFDecoder(h_css.astype(np.uint8))
                for p_erasure in l_p_erasure:
                    for p_err  in l_p_err:
                        # apply noise and get syndromes
                        erasure = np.random.binomial(1, p_erasure, g.n_qbt).astype(np.uint8)
                        pauli_err = np.random.binomial(1, p_err, g.n_qbt).astype(np.uint8)
                        error = np.logical_or(np.logical_and(np.logical_not(erasure), pauli_err), np.logical_and(erasure, np.random.binomial(1, 0.5, g.n_qbt))).astype(np.uint8)
                        syndrome = (g.h @ error % 2).astype(np.uint8)

                        # decode
                        g.ldpc_decode(syndrome, erasure)

                        # recompute syndrome to check decoding
                        err_plus_correction = np.logical_xor(error, g.correction)
                        syndrome = g.h @ err_plus_correction % 2
                        self.assertEqual(0, np.sum(syndrome))


if __name__ == '__main__':
    t = TestDecodeBBCode()
    l_error = [0.01 * (i+1) for i in range(12)]
    l_erasure = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    t.test_decoding(l_error, l_erasure)

