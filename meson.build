project('percolation_c', 'c',
  version : '0.1',
  default_options : ['warning_level=3'])

m_dep = meson.get_compiler('c').find_library('m', required : true)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

examples = [
  'example_2d_toric_decode',
  'example_2d_toric_decode_speed',
  'example_2d_toric_decode_ldpc',
  'example_2d_toric_decode_and_analyze',
  'example_2d_toric_decode_and_analyze_ldpc',
  'example_2d_toric_decode_erasure',
  'example_2d_toric_triangular_decode_erasure',
]

foreach t : examples
  example_exe = executable(t, ['examples/' + t + '.c', 'src/decoder_main.c', 'src/decoder_ldpc.c', 'src/stabilizer_main.c', 'src/graph_type.c', 'src/graph_construction.c'])
endforeach

if get_option('tests') and not meson.is_cross_build()
  subdir('tests')
endif

lib = library('SpeedDecoder', ['src/decoder_main.c', 'src/decoder_ldpc.c', 'src/stabilizer_main.c', 'src/graph_type.c', 'src/graph_construction.c'])
