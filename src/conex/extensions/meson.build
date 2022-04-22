conex_extension_sources = [
  'file.cpp',
  'event.cpp',
  'interaction.cpp',
  'particle.cpp',
  'projectile.cpp'
]

libconexextensions = library(
  'conex_extensions', conex_extension_sources,
  cpp_args: cpp_args,
  link_args: ld_args,
  link_with: [libutil],
  include_directories: include_directories('../../')
)