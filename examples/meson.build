example_progs = {
  'make_stack' : ['make_stack.cpp'],
  'get_alpha' : ['get_alpha.cpp'],
  'fit_conex_profiles' : ['fit_conex_profiles.cpp']
}

foreach name, sources : example_progs
  executable(
    name, sources,
    include_directories: inc_src,
    cpp_args: cpp_args,
    link_args: ld_args,
    link_with: [libcorsika, libconex, libconexextensions, libutil],
    install: true
  )
endforeach