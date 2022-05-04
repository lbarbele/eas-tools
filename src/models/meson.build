models_sources = [
  'dedx_profile.cpp'
]

libmodels = library(
  'models',
  models_sources,
  include_directories : inc_src,
  link_with: [libutil],
  install: true
)