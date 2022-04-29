conex_extension_sources = [
  'file.cpp',
  'event.cpp',
  'interaction.cpp',
  'particle.cpp',
  'projectile.cpp'
]

libconexextensions = library(
  'conex_extensions', conex_extension_sources,
  include_directories: inc_src,
  cpp_args: cpp_args,
  link_args: ld_args,
  link_with: [libutil],
  install: true
)