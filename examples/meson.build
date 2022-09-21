example_progs = [
  'get_alpha',
  'parametrize_conex_profiles',
  'plot_conex_profiles',
  'plot_corsika_profiles',
  'plot_coreas_data',
  'plot_lateral_distribution'
]

foreach name : example_progs
  executable(
    name, name+'.cpp',
    include_directories: inc_src,
    cpp_args: cpp_args,
    link_args: ld_args,
    link_with: [libcorsika, libconex, libconexextensions, libmodels],
    install: true
  )
endforeach
