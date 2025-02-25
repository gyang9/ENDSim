#include <HitmanProc.hh>

#if TENSORFLOW_Enabled && NLOPT_Enabled
#include <cppflow/cppflow.h>
#include <cppflow/model.h>

#include <End.hh>
#include <RAT/DB.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/Log.hh>
#include <RAT/Processor.hh>
#include <RAT/Rat.hh>
#include <filesystem>
#include <iostream>
#include <nlopt.hpp>
#include <numeric>

namespace END {

double dispatch_function(const std::vector<double>& hypothesis, std::vector<double>& grad, void* hitman) {
  // Used to pass the function to be evaluated (NLLH) to NLopt as it requires a global function
  return ((HitmanProc*)hitman)->NLLH(hypothesis);
};

HitmanProc::HitmanProc() : RAT::Processor("hitman") {}

RAT::Processor::Result HitmanProc::Event(RAT::DS::Root* ds, RAT::DS::EV* ev) {
  // Run HITMAN reconstruction in line with simulation on an event by event basis
  RAT::DS::FitResult* fit = new RAT::DS::FitResult("hitman");

  // First we must acquire the observational information from our event
  // The observation is comprised of each PMT hit
  // the x,y,z postion of the PMT hit, the measured photon time, and p.e. pulse integral
  RAT::DS::MC* mc = ds->GetMC();

  int totalhits = 0;

  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
    RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
    int npe = mcpmt->GetMCPhotonCount();

    for (int iPE = 0; iPE < npe; iPE++) {
      double time = mcpmt->GetMCPhoton(iPE)->GetFrontEndTime();
      int pmtID = mcpmt->GetID();

      hitman_hits.push_back(static_cast<double>(x_pos_pmt[pmtID]));
      hitman_hits.push_back(static_cast<double>(y_pos_pmt[pmtID]));
      hitman_hits.push_back(static_cast<double>(z_pos_pmt[pmtID]));
      hitman_hits.push_back(static_cast<double>(time));
      hitman_hits.push_back(
          static_cast<double>(1.0));  // Placeholder for charge in a single p.e. pulse (may be removed)

      totalhits++;
    }
  }
  hitman_charge.push_back(static_cast<double>(totalhits));
  hitman_charge.push_back(static_cast<double>(totalhits));  // one day this will be total detector charge

  // Set the particle kinematics hypothesis seed for the optimization routine
  // Right now this is hard coded, Ideally we would replace these lines with a fast algorithm like quadfit
  double x_seed = 0.0;
  double y_seed = 0.0;
  double z_seed = 0.0;
  double ze_seed = 0.7;
  double az_seed = 3.14;
  double t_seed = 0.0;
  double E_seed = 1.7;

  std::vector<double> hyp_seed{x_seed, y_seed, z_seed, ze_seed, az_seed, t_seed, E_seed};

  // Run the optimization routine
  Optimize(hyp_seed, fit);

  ev->AddFitResult(fit);
  hitman_hits.clear();
  hitman_charge.clear();
  return RAT::Processor::OK;
}

void HitmanProc::BeginOfRun(RAT::DS::Run* run) {
  RAT::DB* db = RAT::DB::Get();
  // Hitman variables stored in ratdb
  RAT::DBLinkPtr table = db->GetLink("Fitter", "HITMAN");
  std::string directory = table->GetS("directory");

  // Set up relative reference EndSimulations/models
  std::string path = "";
  for (auto dir : End::model_directories) {
    std::string fullpath = dir + directory;
    if (std::filesystem::exists(fullpath)) {
      path = fullpath;
      break;
    }
  }

  // Check to see if model exists as an 'approved' EndSimulations model, otherwise use absolute path
  if (!std::filesystem::exists(path)) {
    RAT::warn << "Could not find HITMAN model at EndSimulations/models" + directory << newline;
    RAT::warn << "Attempting to use 'unapproved' network at absolute path " + directory << newline;
    path = directory;
    if (!std::filesystem::exists(path)) {
      RAT::Log::Die("Could not find HITMAN model at: " + path);
    }
  }

  // load HITMAN models in cppflow, note HITMAN requires both a hitnet and a chargenet to function
  this->hitnet = new cppflow::model(path + "/hitnet/");
  this->chargenet = new cppflow::model(path + "/chargenet/");

  // get map from PMT ID to PMT positions for input to HITMAN
  run_pmtinfo = run->GetPMTInfo();
  int npmt = run_pmtinfo->GetPMTCount();

  for (int apmt = 0; apmt < npmt; apmt++) {
    TVector3 pos = run_pmtinfo->GetPosition(apmt);
    x_pos_pmt.push_back(pos.X());
    y_pos_pmt.push_back(pos.Y());
    z_pos_pmt.push_back(pos.Z());
  }
}

void HitmanProc::Optimize(std::vector<double> seed, RAT::DS::FitResult* fit) {
  TVector3 fitPosition(0.0, 0.0, 0.0);
  TVector3 fitDirection(1.0, 0.0, 0.0);
  double fitEnergy = -1.0;
  double fitTime = 0.0;
  double nllh = 0.0;
  bool fitsuccess = false;
  int dimensions = 7;
  std::vector<double> lowerbounds{-950, -950, -950, 0.006, 0, -5, 0.5};  // avoid poles
  std::vector<double> upperbounds{950, 950, 950, 3.138, 6.2832, 5, 6};   // avoid poles

  nlopt::opt opt(nlopt::LN_SBPLX, dimensions);
  opt.set_lower_bounds(lowerbounds);
  opt.set_upper_bounds(upperbounds);
  opt.set_min_objective(dispatch_function, this);
  opt.set_xtol_rel(1e-4);
  double minf;

  try {
    nlopt::result result = opt.optimize(seed, minf);
    fitPosition.SetXYZ(seed[0], seed[1], seed[2]);
    fitDirection.SetTheta(seed[3]);
    fitDirection.SetPhi(seed[4]);
    fitDirection.SetMag(1.);
    fitTime = seed[5];
    fitEnergy = seed[6];
    fitsuccess = true;
  } catch (std::exception& e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
  fit->SetBoolFigureOfMerit("FitValid", fitsuccess);
  fit->SetDoubleFigureOfMerit("NLLH", minf);
  fit->SetPosition(fitPosition);
  fit->SetDirection(fitDirection);
  fit->SetEnergy(fitEnergy);
  fit->SetTime(fitTime);
}

double HitmanProc::NLLH(const std::vector<double>& hypothesis) {
  /* calculate a quantity proportiaonal to the
   * Negative Logarithm of the Likelihood (NLLH) by summing hitnet and chargenet terms
   * with a predefined:
   * hitman_hits flattened std::vector dim (number of p.e. * 5) [x, y, z, t, q, x, y, z, ...])
   * hitman_charge std::vector dim (2) with [nhit,total_q])
   * User specified:
   * hypothesis which contains particle inference parameters {x,y,z,ze,az,t,E} of particle
   * ouptut:
   * scalar, the negative logarithm of the likelihood to evidence ratio
   */

  double sum = 0;

  // copy the hypothesis so its first dimension matches the first dimension of the hit vector
  int totalhits = static_cast<int>(hitman_charge[0]);
  std::vector<double> temp_hyp(totalhits * 7, 0.0);

  for (int ahit = 0; ahit < totalhits; ahit++) {
    int index = ahit * 7;

    temp_hyp[index] = hypothesis[0];
    temp_hyp[index + 1] = hypothesis[1];
    temp_hyp[index + 2] = hypothesis[2];
    temp_hyp[index + 3] = hypothesis[3];
    temp_hyp[index + 4] = hypothesis[4];
    temp_hyp[index + 5] = hypothesis[5];
    temp_hyp[index + 6] = hypothesis[6];
  }

  // prepare input tensors for hitnet by converting them to a cppflow tensors
  cppflow::tensor tensor_hyp(temp_hyp, {totalhits, 7});
  cppflow::tensor tensor_hits(hitman_hits, {totalhits, 5});

  // evaluate hitnet and sum the values produced over each hit
  tensor_hits = cppflow::cast(tensor_hits, TF_DOUBLE, TF_FLOAT);
  tensor_hyp = cppflow::cast(tensor_hyp, TF_DOUBLE, TF_FLOAT);
  cppflow::tensor output =
      (*hitnet)({{"serving_default_input_1:0", tensor_hits}, {"serving_default_input_2:0", tensor_hyp}},
                {"StatefulPartitionedCall:0"})[0];
  for (int index = 0; index < output.get_data<float>().size(); index++) {
    sum -= static_cast<double>(output.get_data<float>()[index]);
  }
  // prepare input tensors for chargenet
  cppflow::tensor tensor_chyp(hypothesis, {1, 7});
  cppflow::tensor tensor_charge(hitman_charge, {1, 2});

  // evaluate chargenet and add it's value to the total llh
  tensor_charge = cppflow::cast(tensor_charge, TF_DOUBLE, TF_FLOAT);
  tensor_chyp = cppflow::cast(tensor_chyp, TF_DOUBLE, TF_FLOAT);
  cppflow::tensor coutput =
      (*chargenet)({{"serving_default_input_3:0", tensor_charge}, {"serving_default_input_4:0", tensor_chyp}},
                   {"StatefulPartitionedCall:0"})[0];
  sum -= static_cast<double>(coutput.get_data<float>()[0]);

  return sum;
}

}  // namespace END
#endif
