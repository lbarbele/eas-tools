project('eas-tools', 'cpp')

cpp_args = run_command('root-config', '--cflags').stdout().strip().split(' ')
ld_args = run_command('root-config', '--libs').stdout().strip().split(' ')

inc_src = include_directories('include')

subdir('src')
subdir('examples')

# install headers
install_subdir('include', install_dir: 'include', strip_directory: true)