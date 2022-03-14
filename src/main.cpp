#include <iostream>
#include <iomanip>

#include "corsika-file.h"

int
main(int argc, char const** argv)
{
  corsika::file file(argv[1]);

  return 0;
}
