#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <vector>

std::vector<double> masses = { -1,
  0.e0        ,.51099893e-3,.51099893e-3,  0.e0     ,.105658372e0, 
  .105658372e0,  .1349766e0, .13957018e0,.13957018e0, 0.497611e0 , //10
  0.493677e0  , 0.493677e0 ,.93956538e0 ,.93827205e0,.93827205e0 ,
  0.497611e0  , 0.547862e0 , 1.115683e0 , 1.18937e0 , 1.192642e0 , //20
  1.197449e0  , 1.31486e0  , 1.32171e0  , 1.67245e0 ,.93956538e0 ,
  1.115683e0  , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31486e0  , //30
  1.32171e0   , 1.67245e0  , 0.e0       , 0.e0      , 0.e0       , 
  0.e0        , 0.e0       , 0.e0       , 0.e0      , 0.e0       , //40
  1.e9        , 580.e0     , 1.e5       , 0.e0      , 0.e0       ,
  0.e0        , 0.e0       , 0.95778e0  , 1.019461e0, 0.78265e0  , //50
  0.7690e0    , 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
  1.2331e0    , 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   , //60
  1.2349e0    , 0.89581e0  , 0.89166e0  , 0.89166e0 , 0.89581e0  ,
  0.e0        , 0.e0       , 0.e0       , 0.e0      , 0.e0       , //70
  0.547862e0  , 0.547862e0 , 0.547862e0 , 0.547862e0, 0.e0       ,
             0,           0,           0,          0,           0, //80
             0,           0,           0,          0,           0,
             0,           0,           0,          0,           0, //90
             0,           0,           0,          0,           0,
             0,           0,           0,          0,           0, //100
             0,           0,           0,          0,           0,
             0,           0,           0,          0,           0, //110
             0,           0,           0,          0,           0,
   1.86484e0  , 1.86961e0  , 1.86961e0  , 1.86484e0  , 1.9683e0   , //120
   1.9683e0   , 2.9836e0   , 2.00697e0  , 2.01027e0  , 2.01027e0  ,
   2.00697e0  , 2.1121e0   , 2.1121e0   , 0.0e0      , 3.096916e0 , //130
   1.77686e0  , 1.77686e0  , 0.e0       , 0.e0       , 0.e0       ,
   0.e0       , 2.28646e0  , 2.46793e0  , 2.47085e0  , 2.45397e0  , //140
   2.4529e0   , 2.45375e0  , 2.5757e0   , 2.5779e0   , 2.6952e0   ,
   0.e0       , 0.e0       , 0.e0       , 2.28646e0  , 2.46793e0  , //150
   2.47085e0  , 2.45397e0  , 2.4529e0   , 2.45375e0  , 2.5757e0   ,
   2.5779e0   , 2.6952e0   , 0.e0       , 0.e0       , 0.e0       , //160
   2.51841e0  , 2.5175e0   , 2.51848e0  , 0.e0       , 0.e0       ,
           0  ,            0,           0,          0,           0, //170
   2.51841e0  , 2.5175e0   , 2.51848e0  , 0.e0       , 0.e0       ,
   5.27961e0  , 5.27929e0  , 5.27929e0  , 5.27961e0  , 5.36679e0  , //180
   5.36679e0  , 6.2751e0   , 6.2751e0   , 5.61951e0  , 5.8155e0   ,
   5.8113e0   , 5.7918e0   , 5.7944e0   , 6.0480e0   , 5.61951e0  , //190
   5.8155e0   , 5.8113e0   , 5.7918e0   , 5.7944e0   , 6.0480e0 
};

int
main(
  int argc,
  char** argv
)
{
  std::ifstream stackFile(argv[1]);

  if (!stackFile.is_open()) {
    std::cout << "stack file not open" << std::endl;
    return 1;
  }

  int npart;
  double totalEnergy;
  double fixHeight;
  double showerZenith;
  double showerAzimuth;

  stackFile >> npart >> totalEnergy >> fixHeight >> showerZenith >> showerAzimuth;

  std::array<double, 17> secpar;

  // accumulators
  double esum = 0;
  double pxsum = 0;
  double pysum = 0;
  double pzsum = 0;

  while (stackFile.good()) {
    for (auto& x : secpar) {
      stackFile >> x;
    }

    if (!stackFile.good()) {
      break;
    }

    // check id
    if (secpar[0] < 0) {
      std::cout << "check_stack: bad id in " << argv[1] << std::endl;
      return 1;
    }

    // check gamma factor
    if (secpar[1] < 1 && secpar[0] != 1) {
      std::cout << "check_stack: bad gamma in " << argv[1] << std::endl;
      return 1;
    }

    // check altitude
    if (!(secpar[5] > 0)) {
      std::cout << "check_stack: bad altitude in " << argv[1] << std::endl;
      return 1;
    }

    // check time
    if (!(secpar[6] >= 0)) {
      std::cout << "check_stack: bad time in " << argv[1] << std::endl;
      return 1;
    }

    // check polarization
    if (secpar[11] != 0 || secpar[12] != 0) {
      std::cout << "check_stack: bad polarization in " << argv[1] << std::endl;
      return 1;
    }

    // check weight
    if (secpar[13] != 1) {
      std::cout << "check_stack: bad weight in " << argv[1] << std::endl;
      return 1;
    }

    // check apparent height
    if (!(secpar[14] > 0)) {
      std::cout << "check_stack: bad happ in " << argv[1] << std::endl;
      return 1;
    }

    // check costhetaearth
    if (!(secpar[16] <= 1 && secpar[16] > 0)) {
      std::cout << "check_stack: bad costhetaEarth in " << argv[1] << std::endl;
      return 1;
    }

    // accumulators
    const double& m = secpar[0] > 200?
      secpar[0]/100 * masses.at(13) : 
      masses.at(int(secpar[0]));

    const double& gamma = secpar[1];
    const double e = m > 0? m*gamma : gamma;
    const double p = std::sqrt((e+m)*(e-m));
    const double pz = p * secpar[2];
    const double px = p * secpar[3];
    const double py = p * secpar[4];

    pzsum += pz;
    pxsum += px;
    pysum += py;

    esum += e;
  }

  const double ptot = std::sqrt(pxsum*pxsum+pysum*pysum+pzsum*pzsum);
  const double dirx = pxsum/ptot;
  const double diry = pysum/ptot;
  const double dirz = pzsum/ptot;

  const double theta = showerZenith*std::acos(-1)/180.0;
  const double phi = showerAzimuth*std::acos(-1)/180.0;

  const double devEnergy = std::fabs(esum/totalEnergy - 1);
  const double devPx = std::fabs(dirx/(std::sin(theta)*std::cos(phi)) - 1);
  const double devPy = std::fabs(diry/(std::sin(theta)*std::sin(phi)) - 1);
  const double devPz = std::fabs(dirz/std::cos(theta) - 1);

  if (devEnergy > 0.1) {
    std::cout << "check_stack: bad energy sum in " << argv[1] << std::endl;
    std::cout << devEnergy << std::endl;
    return 1;
  }

  if (devPx > 0.1) {
    std::cout << "check_stack: bad px sum in " << argv[1] << std::endl;
    return 1;
  }

  if (devPy > 0.1) {
    std::cout << "check_stack: bad py sum in " << argv[1] << std::endl;
    return 1;
  }

  if (devPz > 0.1) {
    std::cout << "check_stack: bad pz sum in " << argv[1] << std::endl;
    return 1;
  }

  return 0;
}