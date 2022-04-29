corsika_sources = [
  'word_t.cpp',
  'binarystream.cpp',
  'binarystream.iterator.cpp',
  'subblock.cpp',
  'binaryfile.cpp',
  'shower.cpp',
  'particle.cpp',
  'particle_iterator.cpp',
  'longfile.cpp',
  'longfile.iterator.cpp',
  'profile.cpp'
]

libcorsika = library(
  'corsika',
  corsika_sources,
  include_directories : inc_src,
  link_with: [libutil],
  install: true
)