conex_extension_sources = [
  'file.cpp',
  'event.cpp',
  'interaction_tree.cpp'
]

libconexextensions = shared_library(
  'conex_extensions', conex_extension_sources,
  include_directories: inc_src,
  cpp_args: cpp_args,
  link_args: ld_args,
  install: true
)