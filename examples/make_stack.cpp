#include <iostream>
#include <iomanip>
#include <string>
#include <array>
#include <list>
#include <functional>
#include <stdexcept>

#include <tclap/CmdLine.h>

#include <conex/file.h>
#include <conex/extensions/file.h>
#include <conex/extensions/interaction_tree.h>
#include <models/atmosphere.h>
#include <util/math.h>
#include <util/point.h>
#include <util/vector.h>

int idToCorsika(const int idConex);

int main(int argc, char** argv) {

  using namespace units::literals;

  // - aliases

  namespace cx = conex;
  namespace ce = conex::extensions;

  using secpar_t = std::array<double, 17>;

  // - parse the command line

  // * define the command line parameters

  TCLAP::CmdLine cmdLine("make_stack");

  TCLAP::ValueArg<std::string> conexFilePath("c", "conex-file",
    "CONEX default root file containing shower data, as wrote by the CxRoot interface.",
    true, "", "conex_file.root", cmdLine);

  TCLAP::ValueArg<std::string> extensionsFilePath("e", "extensions-file",
    "CONEX extensions file containing information about the leading interactions"
    ", activated by the CONEX_EXTENSIONS macro in CONEX",
    true, "", "extensions_file.root", cmdLine);

  TCLAP::ValueArg<double> thresholdRatio("t", "threshold",
    "Minimum fraction of the primary energy required for an interaction to be"
    "analyzed. Interactions with energy below this threshold will be dismissed"
    " and particles below this threshold will go directly to the stack",
    false, 0.05, "threshold", cmdLine
  );

  TCLAP::ValueArg<double> observationLevelInput("o", "obslev",
    "Altitude above sea level of the detector in CORSIKA.",
    false, 0.0, "obslev",  cmdLine
  );

  TCLAP::ValueArg<int> eventNumber("n", "event-number",
    "Number (index) of event in the CONEX files to be analyzed. Only one event can"
    " be declared. Negative values are interpreted as values starting from the end."
    " That is, 0 is the first event in the file and -1 is the last.",
    true, 0, "index", cmdLine
  );

  TCLAP::ValueArg<unsigned int> printPrecision("p", "precision",
    "Numeric precision used to print numbers. Value is limited to 5 <= prec <= 30.",
    false, 12, "prec", cmdLine
  );

  TCLAP::SwitchArg doConsistencyChecks("x", "no-checks",
    "Disable consistency checks when parsing the CONEX extension file.",
    cmdLine, true
  );

  // * actually parse the command line

  cmdLine.parse(argc, argv);

  if (thresholdRatio < 0) {
    std::cout << "threshold ratio must be >= 0. instead it was " << thresholdRatio << std::endl;
    return 1;
  }

  if (printPrecision.getValue() > 30) {
    printPrecision.getValue() = 30;
  } else if (printPrecision.getValue() < 5) {
    printPrecision.getValue() = 5;
  }

  const units::meter_t<double> observationLevel(observationLevelInput);

  // -  open input files and check

  cx::file cxFile(conexFilePath.getValue());
  if (!cxFile.is_open()) {
    std::cerr << "unable to open CONEX file " << conexFilePath.getValue() << std::endl;
    return 1;
  }

  ce::file cxPartFile(extensionsFilePath.getValue());
  if (!cxPartFile.is_open()) {
    std::cerr << "unable to open CONEX extensions file " << extensionsFilePath.getValue() << std::endl;
    return 1;
  }

  // -  check event index and read the event

  const unsigned int trueEvtIndex = eventNumber >= 0? eventNumber : eventNumber + cxPartFile.get_n_events();

  if (trueEvtIndex >= cxPartFile.get_n_events()) {
    std::cerr << "bad event index " << eventNumber << std::endl;
    std::cerr << "file has " << cxPartFile.get_n_events() << " events" << std::endl;
    return 1;
  }
  
  const cx::shower& shower = cxFile.get_shower(trueEvtIndex);
  const ce::event& evt = cxPartFile.get_event(trueEvtIndex, thresholdRatio, doConsistencyChecks);

  // -  get primary interaction data

  ce::interaction_ptr primaryInteraction = evt.get_interaction(0);
  ce::projectile_ptr primaryParticle = primaryInteraction->get_projectile();

  const units::energy_t primaryEnergy = primaryParticle->get_energy();
  const units::energy_t energyThreshold = primaryEnergy * double(thresholdRatio);
  const units::time_t t0 = primaryParticle->get_time();

  // - definition of the detector frame (considering the obs. level)

  // * by default, the detector frame coincides with the standard frame

  util::frame_ptr detectorFrameCnx = util::frame::conex_observer;

  units::angle_t showerThetaObs = shower.get_zenith();

  // * if the observation level is not 0, a change of frame applies

  if (observationLevel > 1_cm) {

    // shower angles and shower axis
    const units::angle_t showerTheta = shower.get_zenith();
    const units::angle_t showerPhi = shower.get_azimuth();
    const util::vector_d showerAxis = shower.get_axis();

    const double sinShowerTheta = util::math::sin(showerTheta);
    const double cosShowerTheta = util::math::cos(showerTheta);

    // auxiliar quantities
    const double aux_a = observationLevel / util::constants::earth_radius;
    const double aux_b = util::math::sqrt((1 + aux_a + sinShowerTheta) * (1 + aux_a - sinShowerTheta)) - cosShowerTheta;

    // the zenithal angle changes due to the altitude of the observation level
    const units::angle_t thetaDet = util::math::asin(sinShowerTheta / (1 + aux_a));

    // the vertical at the detector is rotated with respect to the conex observer frame by the following angle
    const units::angle_t detectorRotAngle = showerTheta - thetaDet;

    // build a rotation to go from the conex observer frame to the detector frame
    const util::square_matrix_d detectorRotationMatrix =
      (util::axis::z, showerPhi) *        // align x axis with projection of the shower axis at ground
      (util::axis::y, detectorRotAngle) * // rotate the z axis by an angle theta_obs_rot
      (util::axis::z, -showerPhi);        // recover original x and y axes

    // the detector is also displaced with respect to the conex observer origin
    const auto detectorPosition = util::point_t<units::length_t>(util::frame::conex_observer) +
      (aux_b * showerAxis - util::vector_d(0, 0, aux_a, util::frame::conex_observer)) *
      util::constants::earth_radius / (1 + aux_a);

    // create the detector frame in conex and define the observation angle
    detectorFrameCnx = util::frame::create(detectorRotationMatrix, detectorPosition, util::frame::conex_observer);

    showerThetaObs = thetaDet;
  }

  // -  compute the interaction tree

  ce::interaction_tree_ptr tree = ce::interaction_tree::create(evt, energyThreshold);

  // -  the stack is a list that will hold secpar_t objects

  std::list<secpar_t> stack;

  // - prepare the consistency checks

  units::energy_t totalEnergySum(0);
  util::vector_t<units::momentum_t> totalMomentumSum(util::frame::corsika_observer);

  // -  loop over the interaction tree to compute the CORSIKA stack

  tree->apply_recursive([&](const ce::interaction_tree& tree){

    // get the projectile of the current interaction
    ce::interaction_ptr interaction = tree.get_interaction();
    ce::projectile_ptr proj = interaction->get_projectile();

    // * particle position in the detector frame

    const util::point_t<units::length_t> pos_cnx = proj->get_position().on_frame(detectorFrameCnx);

    const units::length_t x_cnx = pos_cnx[0];
    const units::length_t y_cnx = pos_cnx[1];
    const units::length_t z_cnx = pos_cnx[2];
    const units::length_t r_cnx = util::math::hypot(x_cnx, y_cnx);

    // * definition of the conex particle frame and relevant angles

    const units::length_t distance_to_earth_center = proj->get_height() + util::constants::earth_radius;

    const units::angle_t phi_cnx = r_cnx > 1_nm? util::math::atan2<units::radian>(y_cnx, x_cnx) : 0_rad;

    const double sin_theta_ea = r_cnx / distance_to_earth_center;
    const double cos_theta_ea = (z_cnx + util::constants::earth_radius) / distance_to_earth_center;
    const units::angle_t theta_ea = util::math::asin(sin_theta_ea);

    const util::square_matrix_d particle_frame_rot_cnx =
      (util::axis::x, 1_pi - theta_ea) *
      (util::axis::z, phi_cnx - 0.5_pi);
    
    const util::frame_ptr particle_frame_cnx = util::frame::create(particle_frame_rot_cnx, detectorFrameCnx);

    // * definition of the corsika particle frame

    // rotation matrix connecting CONEX particle frame to CORSIKA particle frame
    const util::square_matrix_d to_corsika =
      (util::axis::z, -phi_cnx) *
      (util::axis::y, 1_pi);

    // build the corsika particle frame based on the conex particle frame
    const util::frame_ptr particle_frame_csk = util::frame::create(to_corsika, particle_frame_cnx);

    // * compute parameters common to all secondary particles

    secpar_t secpar;

    // true particle height [cm]
    secpar[5] = proj->get_height() / 1_cm;

    // accumulated time since the first interaction [s]
    secpar[6] = (proj->get_time() - t0) / 1_s;

    // particle position in CORSIKA curved coordinates (following earth's curvature) [cm]
    const units::angle_t phi_csk = phi_cnx - 0.5_pi;
    const units::length_t r_csk = util::constants::earth_radius * (theta_ea / 1_rad);
    secpar[7] = util::math::cos(phi_csk) * r_csk / 1_cm;
    secpar[8] = util::math::sin(phi_csk) * r_csk / 1_cm;

    // particle generation
    secpar[9] = 1 + tree.get_generation();

    // level of the last interaction [cm]
    secpar[10] = proj->get_height() / 1_cm;

    // polarization (not used)
    secpar[11] = 0;
    secpar[12] = 0;

    // particle weight
    secpar[13] = proj->get_weight();

    // apparent height (above sea level) [cm]
    secpar[14] = z_cnx / 1_cm;

    // cosine of apparent polar angle of particle
    const units::length_t z_above_obs = z_cnx - observationLevel;
    secpar[15] = z_above_obs / util::math::hypot(x_cnx, y_cnx, z_above_obs);

    // cosine of polar angle measured at earth center
    secpar[16] = cos_theta_ea;

    // * loop over final products of the current interaction

    for (const ce::particle_ptr& particle : tree.get_final_products()) {

      // particle ID
      secpar[0] = idToCorsika(particle->get_id());

      // gamma factor (or energy in GeV, if particle is massless)
      secpar[1] = particle->get_mass() > units::mass_t(0)?
        particle->get_energy() / (particle->get_mass() * 1_c * 1_c) :
        particle->get_energy() / 1_GeV;

      // particle direction in CORSIKA particle frame
      const util::vector_d direction_corsika = particle->get_momentum().on_frame(particle_frame_csk).get_normalized(1);
      secpar[2]  = -direction_corsika[2];
      secpar[3]  = direction_corsika[0];
      secpar[4]  = direction_corsika[1];

      // add to the stack
      stack.push_back(secpar);

      // * consistency checks

      totalEnergySum += particle->get_energy();
      totalMomentumSum += particle->get_momentum();
    }
  });

  // - perform the consistency checks

  if (std::fabs(totalEnergySum/primaryEnergy - 1) > 1e-5) {
    std::cerr << "bad energy sum" << std::endl;
    return 1;
  }

  const util::vector_t<units::momentum_t> pprim = primaryParticle->get_momentum();

  if ((totalMomentumSum - pprim).norm()/pprim.norm() > 1e-5) {
    std::cerr << "bad momentum sum" << std::endl;
    return 1;
  }

  // * width used to print secpar fields
  unsigned int fieldWidth = printPrecision + 6;
  std::cout.precision(printPrecision);

  // * print the stackin file header
  std::cout
    << std::setw(fieldWidth) << stack.size() << ' '
    << std::setw(fieldWidth) << shower.get_energy()/1_GeV << ' '
    << std::setw(fieldWidth) << primaryParticle->get_height()/1_cm << ' '
    << std::setw(fieldWidth) << showerThetaObs/1_deg << ' '
    << std::setw(fieldWidth) << util::math::fmod(shower.get_azimuth() + 90_deg, 360_deg)/1_deg << '\n';

  // * print the stack
  for (const auto& secpar : stack) {
    for (size_t i = 0; i < secpar.size(); ++i) {
      std::cout << std::setw(fieldWidth) << secpar[i] << (i == secpar.size()-1 ? '\n' : ' ');
    }
  }

  return 0;
};

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