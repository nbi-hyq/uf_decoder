tests = [
  'toric_2d_square_alg1',
  'toric_2d_square_alg2',
  'toric_2d_square_alg3',
  'toric_2d_triangular_alg1',
  'toric_2d_triangular_alg2',
  'toric_2d_triangular_alg3',
  'logical_errors',
]

foreach t : tests
  test_exe = executable(
    'test_' + t,
    [t + '.c', '../src/decoder_main.c', '../src/graph_type.c', '../src/graph_construction.c'],
    include_directories : ['../inc/'],
  )
  test(t, test_exe)
endforeach
