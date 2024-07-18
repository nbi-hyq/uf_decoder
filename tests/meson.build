tests = [
  'toric2d_square_alg1',
  'toric2d_square_alg2',
  'toric2d_square_alg4',
  'toric2d_square_ldpc',
  'toric2d_triangular_alg1',
  'toric2d_triangular_alg2',
  'toric2d_triangular_alg4',
  'toric2d_triangular_ldpc',
  'logical_errors',
]

foreach t : tests
  test_exe = executable(
    'test_' + t,
    [t + '.c', '../src/decoder_main.c', '../src/decoder_ldpc.c', '../src/stabilizer_main.c', '../src/graph_type.c', '../src/graph_construction.c'],
    include_directories : ['../inc/'],
  )
  test(t, test_exe)
endforeach
