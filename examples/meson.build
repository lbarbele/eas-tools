example_progs = [
  'make_stack',
  'check_stack',
  'get_alpha',
  'fit_conex_profiles',
  'fit_conex_anomalous',
  'select_anomalous',
  'plot_conex_profiles'
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