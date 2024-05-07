tests = [
  'toric_2d_square_alg1',
  'toric_2d_square_alg2',
  'toric_2d_square_alg3',
  'toric_2d_triangular_alg1',
  'toric_2d_triangular_alg2',
  'toric_2d_triangular_alg3',
  'logical_errors',
  'toric_2d_square_ldpc',
]

foreach t : tests
  test_exe = executable(
    'test_' + t,
    [t + '.c', '../src/decoder_main.c', '../src/decoder_ldpc.c', '../src/stabilizer_main.c', '../src/graph_type.c', '../src/graph_construction.c'],
    include_directories : ['../inc/'],
  )
  test(t, test_exe)
endforeach
