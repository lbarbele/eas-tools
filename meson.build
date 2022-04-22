project('eas-sim-readers', 'cpp')

cpp_args = run_command('root-config', '--cflags').stdout().strip().split(' ')
ld_args = run_command('root-config', '--libs').stdout().strip().split(' ')

inc_src = include_directories('src')

subdir('src')