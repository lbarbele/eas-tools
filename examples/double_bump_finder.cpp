#include <iostream>
#include <iomanip>

#include <conex/file.h>
#include <conex/shower.h>


int
main(
  int argc,
  char** argv
)
{
  for (int ifile = 1; ifile < argc; ++ifile) {
    conex::file cxFile(argv[ifile]);

  }

  return 0;
}