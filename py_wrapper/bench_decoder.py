import numpy as np
import matplotlib.pyplot as plt
import UnionFindPy
import pymatching
import sys
import time
from run_from_py import toric_code_x_logicals, toric_code_x_stabilisers, TannerGraphDecoder
from bench_decoder_analysis import num_decoding_failures_batch, num_decoding_failures_pymatching_batch, num_decoding_failures_others
from scipy.sparse import csr_matrix


DEFAULT_NUM_TRIALS = 10000

if len(sys.argv) not in [4, 5]:
    print("Usage: {} [uf | matching] L p [num_trials]".format(sys.argv[0]))
    print(f'   Default value for num_trials is {DEFAULT_NUM_TRIALS}')

if sys.argv[1].strip().lower() not in ['uf', 'matching', 'uf_reference']:
    print("The first argument must be 'uf', 'matching', 'uf_reference'")
    sys.exit(1)

L = int(sys.argv[2])
p = float(sys.argv[3])

if len(sys.argv) == 5:
    num_trials = int(sys.argv[4])
else:
    num_trials = DEFAULT_NUM_TRIALS

Hx = toric_code_x_stabilisers(L)
logX = toric_code_x_logicals(L)
if sys.argv[1].strip().lower() == 'uf':
    total_time, num_errors = num_decoding_failures_batch(Hx, logX, p, num_trials)
elif sys.argv[1].strip().lower() == 'matching':
    total_time, num_errors = num_decoding_failures_pymatching_batch(Hx, logX, p, num_trials)
elif sys.argv[1].strip().lower() == 'uf_reference':  # https://github.com/chaeyeunpark/UnionFind needs csr_matrix
    total_time, num_errors = num_decoding_failures_others(UnionFindPy.Decoder, csr_matrix(Hx), logX, p, num_trials)

print('{}\t{}'.format(total_time / num_trials, num_errors / num_trials))
