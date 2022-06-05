example_dirs = [
  'double-bump',
]

foreach dir : example_dirs
  subdir(dir)
endforeach

example_progs = [
  'make_stack',
  'check_stack',
  'get_alpha',
  'fit_conex_profiles',
  'plot_conex_profiles',
  'plot_corsika_profiles',
  'plot_coreas_data'
]

foreach name : example_progs
  executable(
    name, name+'.cpp',
    include_directories: inc_src,
    cpp_args: cpp_args,
    link_args: ld_args,
    link_with: [libcorsika, libconex, libconexextensions, libutil, libmodels],
    install: true
  )
endforeach