conex_sources = [
  'file.cpp',
  'file.iterator.cpp'
]

libconex = shared_library(
  'conex', conex_sources,
  include_directories : inc_src,
  cpp_args: cpp_args,
  link_args: ld_args,
  install: true
)

subdir('extensions')