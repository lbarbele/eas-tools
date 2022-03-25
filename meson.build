project('eas-sim-readers', 'cpp')

sources = [
  'src/main.cpp',
  'src/corsika-word_t.cpp',
  'src/corsika-fstream.cpp',
  'src/corsika-fstream-iterator.cpp',
  'src/corsika-subblock.cpp',
  'src/corsika-file.cpp',
  'src/corsika-shower.cpp'
]

sources += [
  'src/conex-file.cpp'
]

cpp_args = run_command('root-config', '--cflags').stdout().strip().split(' ')
ld_args = run_command('root-config', '--libs').stdout().strip().split(' ')

includes = include_directories('src')

executable('main', sources, include_directories: includes, cpp_args: cpp_args, link_args: ld_args)