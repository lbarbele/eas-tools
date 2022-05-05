#include <iostream>
#include <iomanip>
#include <string>
#include <array>
#include <list>
#include <algorithm>
#include <cctype>
#include <regex>

#include <conex/file.h>
#include <conex/extensions/file.h>
#include <util/vector.h>

// use aliases to avoid long namespace chains
namespace ext {using namespace conex::extensions;}
namespace cx  {using namespace conex;}

// a struct encapsulating the SECPAR data
struct secpar_t : public std::array<double, 17> {
  // create a secpar object given a mother projectile and its child secondary particle (also a t0 value in seconds)
  static secpar_t create(const ext::projectile& prj, const ext::particle& prt, const double t);
};

// a stack is a list of secpar_t's
struct stack_t : public std::list<secpar_t> {};

// secpar printer
template<class CharT> std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& stream, const secpar_t& p);
// conex rotate subroutine
util::vector rotate(const util::vector& v, const double st, const double ct, const double sp, const double cp);
// function that actually creates a stack
stack_t makeStack(cx::shower& shower, ext::event& event, const double threshold);
// compare a particle and a projectile
bool compare(const ext::projectile& mother, const ext::particle& child, const ext::projectile& proj);
// convert CONEX (NEXUS) particle ID to CORSIKA ID
int idToCorsika(const int id);

int
main(
  int argc,
  char** argv
)
{
  // parse command line
  // syntax is: ./make_stack conex_standard_file.root conex_extensions_file.root eventNumber
  if (argc != 4) {
    std::cerr << "syntax error!" << std::endl;
    return 1;
  }

  // open the conex file and check
  cx::file cxFile(argv[1]);
  if (!cxFile.is_open()) {
    std::cerr << "bad conex file!" << std::endl;
    return 1;
  }

  // conex extensions file
  ext::file extFile(argv[2]);
  if (!extFile.is_open()) {
    std::cerr << "bad conex extensions file!" << std::endl;
    return 1;
  }

  // get the event number
  const std::string str = argv[3];
  if (!std::regex_match(str, std::regex("[-+]?[0-9]+"))) {
    std::cerr << "bad event nubmer " << str << std::endl;
    return 1;
  }

  const int offset = std::stoi(argv[3]);

  if (offset >= cxFile.get_n_showers() || offset+cxFile.get_n_showers() < 0) {
    std::cerr << "event number is larger than the number of events" << std::endl;
    return 1;
  } 

  const uint ievent = offset > 0? offset : cxFile.get_n_showers() + offset;
  
  if (ievent >= extFile.get_n_events()) {
    std::cerr << "event number is larger than the number of events" << std::endl;
    return 1;
  }

  // get the event
  auto shower = cxFile.get_shower(ievent);
  auto event = extFile.get_event(ievent);

  // create the stack
  auto stack = makeStack(shower, event, 0.005);

  // check the stack
  if (stack.empty()) {
    std::cerr << "failed to create the stack. it is empty!" << std::endl;
    return 1;
  }

  // print the stack header containing
  // - number of particles (required by CORSIKA)
  // - total energy (required by CORSIKA) [GeV]
  // - first interaction altitude [cm]
  // - zenith angle of shower axis [deg]
  // - azimuth angle of shower axis [deg]
  // the first three fields are intented to be used in the corsika input data cards
  std::cout
    << std::setw(13) << stack.size()
    << std::setw(13) << shower.get_energy_gev()
    << std::setw(13) << event.get_interaction(0).get_projectile().get_height()*100
    << std::setw(13) << shower.get_zenith_deg()
    << std::setw(13) << std::fmod(shower.get_azimuth_deg() + 90, 360)
    << std::endl;

  // print all particles in the stack
  for (auto& secpar : stack) {
    std::cout << secpar << std::endl;
  }

  return 0;
}



template<class CharT>
std::basic_ostream<CharT>&
operator<<(
  std::basic_ostream<CharT>& stream,
  const secpar_t& p
)
{
  for (const auto& x : p) {stream << std::setw(13) << x;}
  return stream;
}



util::vector
rotate(
  const util::vector& v,
  const double st,
  const double ct,
  const double sp,
  const double cp
)
{
  return {
     cp*v[0] + ct*sp*v[1] + st*sp*v[2],
    -sp*v[0] + ct*cp*v[1] + st*cp*v[2],
           0 -    st*v[1] +    ct*v[2]
  };
}



stack_t
makeStack(
  cx::shower& shower,
  ext::event& event,
  const double threshold
)
{
  // the stack to be filled
  stack_t stack;

  // consistency checks
  double esum = 0;
  util::vector psum = {0,0,0};

  // time of the first interaction
  const double t0 = event.get_interaction(0).get_projectile().get_time() / util::constants::c;

  // create a list of all interactions except the first one
  std::list<ext::interaction> remaining;
  for (int i = 1; i < event.get_n_interactions(); ++i) {
    remaining.push_back(event.get_interaction(i));
  }

  // start a queue of interesting interacions, whose secondaries will either:
  // - be searched as a projectile of the remaining interactions, making such
  //   interaction enter into this queue, or
  // - be sent to the stack, if not participating of any other interaction as
  //   a projectile
  // we start this queue with the first interaction only
  std::list<ext::interaction> queue{event.get_interaction(0)};

  // analyze each interaction in the queue until it is empty
  while (!queue.empty()) {
    // get a reference to the current mother interaction
    const auto& interaction = queue.front();

    // get the mother projectile
    const auto& mother = interaction.get_projectile();

    // now, loop over all secondaries emerging from this interaction and check
    // if they participate in some other interaction
    for (const auto& child : interaction) {

      // some particles have zero momenta. skip them
      if (child.get_momentum().norm() <= 0) {
        continue;
      }

      // this flag indicates if we have found the current child, secondary particle
      // as a projectile in another interaction
      bool found = false;

      // we only look for particles whose energy is above a certain threshold, given
      // as a fraction of the shower energy
      if (child.get_energy() > threshold * shower.get_energy_gev()) {

        // iterate over the remaining interactions to search for a projectile that matches
        // the current child particle
        auto it = remaining.begin();
        while (it != remaining.end()) {
          // get the projectile of this interaction
          const auto& proj = it->get_projectile();

          // check if the projectile matches the characteristics of the child particle
          const bool isSame = compare(mother, child, proj);

          // if there is a match between the particle and the projectile, then
          // - add the child interaction to the queue of interactions to be analyzed
          // - erase the child interaction from the remaining list
          // - set the "found" flag to true, indicating we have found the child interaction
          // - break the loop over the remaining interactions
          // otherwise, simply continue the loop and look into the next interaction
          if (isSame) {
            queue.push_back(*it);
            it = remaining.erase(it);
            found = true;
            break;
          } else {
            ++it;
          }
        } // loop over remaining interactions

      } // if (energy cut)

      // if the particle hasn't been found in any other interaction, create a secpar_t
      // object from it, corresponding to CORSIKA's secpar fields and put it into the stack
      if (!found) {
        esum += child.get_energy();
        psum += rotate(child.get_momentum(), mother.s0s, mother.c0s, mother.s0xs, mother.c0xs).to_obs(mother.get_phip(), mother.get_thetap());
        stack.push_back(secpar_t::create(mother, child, t0));
      }

    } // loop over child particles


    // pop the current mother interaction from the queue
    queue.pop_front();
  }

  //
  // consistency checks of energy/momentum
  //
  const double energyDev = esum/shower.get_energy_gev() - 1;
  if (std::fabs(energyDev) > 1e-5) {
    std::cerr << "warning: bad energy sum. deviation is of " << 100*energyDev <<  "%" << std::endl;
    return stack_t();
  }

  const auto pprim = event.get_interaction(0).get_projectile().get_momentum();
  if ((psum-pprim).norm()/pprim.norm() > 1e-5) {
    std::cerr << "warning: bad p sum " << std::endl;
    std::cerr << (psum-pprim)/pprim.norm() << std::endl;
    return stack_t();
  }
  
  return stack;
}



bool
compare(
  const ext::projectile& mother,
  const ext::particle& child,
  const ext::projectile& proj
)
{
  static const double tolerance = 1e-5;

  // first thing is to check IDs
  if (proj.get_id() != child.get_id()) {
    return false;
  }

  // also, check energy
  const double dev_energy = proj.get_energy()/child.get_energy() - 1;
  if (std::fabs(dev_energy) > tolerance) {
    return false;
  }

  // check momenta (in observer frame)
  const auto ppart = rotate(child.get_momentum(), mother.s0s, mother.c0s, mother.s0xs, mother.c0xs).to_obs(mother.get_phip(), mother.get_thetap());
  const auto pproj = proj.get_momentum();

  const double pdev = (pproj-ppart).norm()/ppart.norm();
  if (std::fabs(pdev) > tolerance) {
    return false;
  }

  // check position (in observer frame)
  auto posproj = proj.get_position();
  auto pospart = mother.get_position() + (ppart/child.get_energy()) * (proj.get_time() - mother.get_time());

  const double rdev = (posproj - pospart).norm()/pospart.norm();
  if (std::fabs(rdev) > tolerance) {
    return false;
  }

  return true;
}



secpar_t
secpar_t::create(
  const ext::projectile& proj,
  const ext::particle& part,
  const double tStart
)
{
  // the secpar object that will be filled, then returned
  secpar_t secpar;

  // some constants we use below
  using util::constants::earth_radius;
  using util::constants::c;

  // * compute everyting in CONEX units (distance in m, id from NEXUS)
  // * then convert to CORSIKA when copying to secpar below

  // particle's CONEX id
  const double id = part.get_id();

  // particle's generation
  const double gen = proj.get_generation() + 1;

  // energy and mass (both in GeV)
  const double mass = part.get_mass();
  const double energy = part.get_energy();

  // the Lorentz factor
  const double gammaFactor = mass > 0? energy/mass : energy;

  // time since tStart
  const double t = proj.get_time() / c - tStart;

  // particle's position in the observer frame
  const double x = proj.get_x();      // x coordinate in cartesian frame with origin in the shower core
  const double y = proj.get_y();      // y coordinate in cartesian frame with origin in the shower core
  const double r = proj.get_r();      // sqrt(x*x + y*y)
  const double h = proj.get_height(); // height a.s.l. (measured in the line that connects the particle position and the Earth center)

  // sine/cosine of the projection of the particle position on the ground
  const double sinp = y/r;
  const double cosp = x/r;
  
  // particle momentum in its own frame. this frame follows from the definition of
  // the curved atmosphere coordinates both in corsika and conex. this frame is specified by:
  // - z axis points from the particle position towards the Earth center
  // - y axis is contained in the plane that contains the z axis and the shower core position
  // - x axis is given assuming the frame is right-handed
  // the values of s0s, c0s, s0xs, c0xs are stored by CONEX in the extensions file to
  // allow for this rotation. originally, the particle momentum is given in the interaction frame,
  // in which the z axis coincides with the projectile moving direction
  const auto p = rotate(part.get_momentum(), proj.s0s, proj.c0s, proj.s0xs, proj.c0xs);

  // particle direction in CONEX
  const auto pNorm = p / p.norm();
  // particle direction in CORSIKA (simple rotation around the z axis)
  const auto pNormRot = util::vector::to_obs(pNorm, sinp, cosp, 0.0, 1.0);

  // particle position in CORSIKA's curved atmosphere:
  // - distance between particle and Earth center
  const double distCenter = h + earth_radius;
  // - apparent height, measured in the line the connects the shower core position and the earth center
  const double happ = std::sqrt((distCenter-r)*(distCenter+r)) - util::constants::earth_radius;
  // - zenith of particle position in a frame at the Earth center
  const double costea = (happ + earth_radius) / distCenter;
  // - distance between shower core and projection of the particle position at Earth, considering Earth's curvature
  const double rCurved = earth_radius * std::acos(costea);


  secpar[0]  = idToCorsika(id);
  secpar[1]  = gammaFactor;
  secpar[2]  = pNorm[2];     // the z axis in CORSIKA coincides with the z axis in CONEX
  secpar[3]  =  pNormRot[1]; // a rotation around the z axis of (pi/2) is applied to obtain
  secpar[4]  = -pNormRot[0]; // CORSIKA's coordinates as a function of CONEX's ones
  secpar[5]  = h*100;        // altitude.
  secpar[6]  = t;            // time since tStart
  secpar[7]  =  rCurved * sinp * 100; // x and y coordinates in CORSIKA is measured along the Earth surface
  secpar[8]  = -rCurved * cosp * 100; // instead of along a plane observation level as in CONEX
  secpar[9]  = gen;          // particle's generation is gen_proj + 1
  secpar[10] = h*100;        // level of last interaction, same as particle's height
  secpar[11] = 0;            // polarization, not used
  secpar[12] = 0;            // polarization, not used
  secpar[13] = 1;            // thinning weight is simply one. we don't use thinning in CONEX
  secpar[14] = happ * 100;   // apparent height
  secpar[15] = -1;           // apparent zenith angle of particle position (takes into account obs. level, computed in CORSIKA)
  secpar[16] = costea;       // cosine of zenith of particle position for a frame at the earth center

  return secpar;
}



int
idToCorsika(
  const int idConex
)
{
  std::map<int, int> codevt = {
    {   10,   1}, // gamma
    {  -12,   2}, // e+
    {   12,   3}, // e-
    
    {  -14,   5}, // mu+
    {   14,   6}, // mu-
    {  110,   7}, // pi0
    {  120,   8}, // pi+
    { -120,   9}, // pi-
    {  -20,  10}, // klong
    {  130,  11}, // k+
    { -130,  12}, // k-
    { 1220,  13}, // n
    { 1120,  14}, // p
    {-1120,  15}, // p bar
    {   20,  16}, // kshort
    {  220,  17}, // eta
    { 2130,  18}, // lambda
    { 1130,  19}, // sigma+
    { 1230,  20}, // sigma0
    { 2230,  21}, // sigma-
    { 1330,  22}, // xi0
    { 2330,  23}, // xi-
    { 3331,  24}, // omega-
    {-1220,  25}, // n bar
    {-2130,  26}, // lambda bar
    {-1130,  27}, // sigma- bar
    {-1230,  28}, // sigma0 bar
    {-2230,  29}, // sigma+ bar
    {-1330,  30}, // xi0 bar
    {-2330,  31}, // xi+ bar
    {-3331,  32}, // omega+ bar
    {  230,  33}, // k0
    { -230,  34}, // k0b
    
    {   41,  41}, // QBall
    
    {   43,  43}, // Monopole
    
    {  221,  50}, // omega
    {  111,  51}, // rho0  (=-10 in QII (See below))
    {  121,  52}, // rho+
    { -121,  53}, // rho-
    { 1111,  54}, // delta++
    { 1121,  55}, // delta+
    { 1221,  56}, // delta0
    { 2221,  57}, // delta-
    {-1111,  58}, // delta--
    {-1121,  59}, // delta-
    {-1221,  60}, // delta0 bar
    {-2221,  61}, // delta+
    {  231,  62}, // k*0
    {  131,  63}, // k*+
    { -131,  64}, // k*-
    { -231,  65}, // k*0b
    {   11,  66}, // nu_e-
    {  -11,  67}, // nu_e+
    {   13,  68}, // nu_mu-
    {  -13,  69}, // nu_mu+
    
    { -140, 116}, // D0(1.864)
    { -240, 117}, // D(1.869)+
    {  240, 118}, // Db(1.869)-
    {  140, 119}, // D0b(1.864)
    { -340, 120}, // Ds+
    {  340, 121}, // Ds-
    {  440, 122}, // etac
    { -141, 123}, // D*0
    { -241, 124}, // D*+
    {  241, 125}, // D*-
    {  141, 126}, // D*0b
    { -341, 127}, // Ds*+
    {  341, 128}, // Ds*-
    
    {  441, 130}, // J/psi
    {  -16, 131}, // tau+
    {   16, 132}, // tau-
    {   15, 133}, // nu_tau-
    {  -15, 134}, // nu_tau+
    
    { 2140, 137}, // LambdaC(2.285)+
    { 3140, 138}, // Xic+
    { 3240, 139}, // Xic0
    { 1140, 140}, // sigmac++
    { 1240, 141}, // sigmac+
    { 2240, 142}, // sigmac0
    { 1340, 143}, // Xi'c+
    { 2340, 144}, // Xi'c0
    { 3340, 145}, // omegac0

    {-2140, 149}, // LambdaC(2.285)-
    {-3140, 150}, // Xic-
    {-3240, 151}, // Xic0 bar
    {-1140, 152}, // sigmac--
    {-1240, 153}, // sigmac-
    {-2240, 154}, // sigmac0 bar
    {-1340, 155}, // Xi'c-
    {-2340, 156}, // Xi'c0 bar
    {-3340, 157}, // omegac0 bar
    
    { 1141, 161}, // sigma*c++
    { 1241, 162}, // sigma*c+
    { 2241, 163}, // sigma*c0

    {-1141, 171}, // sigma*c--
    {-1241, 172}, // sigma*c-
    {-2241, 173}, // sigma*c0 bar
    
    {   17, 201}, // Deuteron
    
    {   18, 301}, // Triton
      
    {   19, 402} // Alpha
  };

  int idCorsika = -1;

  if (idConex == 41) {
    // qball
    idCorsika = idConex;
  } else if (idConex == 43) {
    // monopole
    idCorsika = idConex;
  } else {
    if (idConex%100 != 0) {
      idCorsika = codevt.at(idConex);
    } else if (idConex < 0) {
      // strangelet
      idCorsika = 42;
    } else {
      // nuclei
      int nChrg = idConex/100;
      if (idConex == 3) {
        nChrg = 1;
      } else {
        nChrg = int(double(nChrg) / 2.15 + 0.7);
      }
      idCorsika = idConex + nChrg;
    }
  }

  if (idCorsika <= 0) {
    throw("bad code conversion!");
  }

  return idCorsika;
}