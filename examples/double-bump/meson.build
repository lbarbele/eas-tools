double_bump_progs = [
  'fit_conex_anomalous',
  'select_anomalous',
  'plot_anomalous_profile'
]

foreach name : double_bump_progs
  executable(
    name, name+'.cpp',
    include_directories: inc_src,
    cpp_args: cpp_args,
    link_args: ld_args,
    link_with: [libcorsika, libconex, libconexextensions, libutil, libmodels],
    install: true
  )
endforeach