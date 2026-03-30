#include <iostream>
#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "Pythia8/Pythia.h"

using namespace fastjet;
using namespace Pythia8;
using namespace std;

// e+e- eec changed in .cmnd file 

// Function to compute Delta R
float deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  float dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  float deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}

// Function to calulate opening angle

float costheta(const PseudoJet& p1, const PseudoJet& p2) {
    float dotprod = (p1.px() * p2.px() +  p1.py() * p2.py() +  p1.pz() * p2.pz());
    float normp1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz()); 
    float normp2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz()); 
    
    if (normp1 == 0 || normp2 == 0) // division by 0 check 
        return 0; 
    
    return dotprod / (normp1 * normp2);
}

const Int_t nDeltaRBinsEECw = 25;
const Float_t minDeltaREECw = 1e-4; // Avoid log(0) issue
const Float_t maxDeltaREECw = 0.5;
const Float_t binnerShiftw = 0.00001;
const Float_t deltaRlogBinWidthw =
  (TMath::Log(maxDeltaREECw + binnerShiftw) - TMath::Log(minDeltaREECw + binnerShiftw)) / nDeltaRBinsEECw;
Float_t deltaRBinsEECw[nDeltaRBinsEECw + 1];

const Int_t nDeltaRBinsEEC = 13;
const Float_t minDeltaREEC = .025; // Avoid log(0) issue                             
const Float_t maxDeltaREEC = 1;
const Float_t binnerShift = 0.1;
const Float_t deltaRlogBinWidth =
  (TMath::Log(maxDeltaREEC + binnerShift) - TMath::Log(minDeltaREEC + binnerShift)) \
  / nDeltaRBinsEEC;
Float_t deltaRBinsEEC[nDeltaRBinsEEC + 1];


//allows for the .cmnd file to be accepted as an argument 
int main(int argc, char* argv[]) {
  // Check if a .cmnd file is provided
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <config.cmnd>" << endl;
    return 1;
  } 

    string configFile = argv[1];
    string outputFileName = argv[2];
    
    // Create output ROOT file
    TFile *out = new TFile(outputFileName.c_str(), "RECREATE");
    
    //TFile *out = new TFile("testthermal.root", "RECREATE");
    
    if (!out || out->IsZombie()) {
        cerr << "Error: Could not open output file: " << outputFileName << endl;
        return 1;
    }

    out->cd();
    cout << "Created output file: " << outputFileName << endl;
  
  Double_t eecbounds[41] = {1.00000000e-05, 1.71770469e-05, 2.95050939e-05, 5.06810379e-05,
    8.70550563e-05, 1.49534878e-04, 2.56856761e-04, 4.41204061e-04,
    7.57858283e-04, 1.30177672e-03, 2.23606798e-03, 3.84090444e-03,
    6.59753955e-03, 1.13326246e-02, 1.94661024e-02, 3.34370152e-02,
    5.74349177e-02, 9.86562273e-02, 1.69462264e-01, 2.91086125e-01,
    5.00000000e-01,
    0.7089138754019768, 0.8305377361330085, 0.9013437726906997, 0.9425650822501482, 
    0.9665629847511789, 0.9805338976261914, 0.9886673753979593, 0.9934024604461356, 
    0.9961590955587668, 0.9977639320225002, 0.9986982232761837, 0.9992421417167447, 
    0.9995587959386655, 0.9997431432392584, 0.9998504651218779, 0.9999129449436703, 
    0.9999493189620527, 0.9999704949061466, 0.9999828229531487, 1.0};
  // Initialize histogram
  TH1::SetDefaultSumw2();
 //TH2::SetDefaultSumw2();
  int bins = 40;
  
  TH1F EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
  //TH1F EEC_w_p("EEC_w_p", "Whole Event EEC", bins, 0, bins);
  
  
  TH1F Aj_spectrum("Aj_spectrum", "Jet Assymetry spectrum", 60, 0.0, 1.0);
  
  // Initialize Pythia for 200 GeV pp collision
  Pythia pythia;
 
  if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
  }
  
  pythia.init();
  

  //float jet_radius = 0.4;
  fastjet::JetDefinition jet_def(ee_kt_algorithm);
  int dijet_event_counter = 0;
  Sphericity sph;
   
  for (int iEvent = 0; iEvent < 100000; ++iEvent) { // 10000 events for now
    if (!pythia.next()) continue;
    
    // implementing event sphericity cut from lunas paper 
    sph.analyze(pythia.event);
    if (sph.sphericity() < (7*M_PI)/36 || sph.sphericity() < (29*M_PI)/36) continue;
    
    std::vector<fastjet::PseudoJet> event;

    // pT cut > .5 GeV
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      
      if (p.isFinal() && p.pT() > .2 && p.isCharged() == 1) {
	    event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
	    
      }
    double totale = 0;
    
    for (int i = 0; i < event.size(); ++i) {
        totale = totale + event.at(i).e();
     }
    cout  << "Total energy = "<< totale << endl;
    
    }//! pseudojet of all particles in event loop close
    
    
    //jet asymmetry  
  //float aj = (jets[0].e()-jets[1].e()) / (jets[0].e() + jets[1].e());
  //Aj_spectrum.Fill(aj);
 
  for (size_t i = 0; i < event.size(); ++i) {
    for (size_t j = i + 1; j < event.size(); ++j) {
      float eec = event.at(i).e() * event.at(j).e();
      //float delr = deltaR(event.at(i), event.at(j));
      float q2 = 91.2*91.2;
      float ctheta = costheta(event[i], event[j]);
      float z = (1 - ctheta)/2 ; 
      EEC_w.Fill(z, eec/(q2)); // weighted by the energy scale q2^2
      // filling balanced and unbalanced whole EECs
      
    } 
  }//!Whole event EEC loop close 

}//! event loop 
  
  out->Write();
  
  // Pythia cleanup
  pythia.stat();
  out->Close();
  return 0;
}

