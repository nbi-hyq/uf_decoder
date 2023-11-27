tests = [
  'toric_2d_square',
  'toric_2d_triangular',
]

foreach t : tests
  test_exe = executable(
    'test_' + t,
    [t + '.c', '../src/percolation_main.c', '../src/graph_type.c', '../src/graph_construction.c'],
    include_directories : ['../inc/'],
  )
  test(t, test_exe)
endforeach
