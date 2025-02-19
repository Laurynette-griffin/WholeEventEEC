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



// Function to compute Delta R
float deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  float dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  float deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
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

  // Initialize logarithmic binning for whole event (make this a callable function please)
  for (int iDeltaRw = 0; iDeltaRw <= nDeltaRBinsEECw; iDeltaRw++) {
    deltaRBinsEECw[iDeltaRw] =
      (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw;
  }

std::vector<double> binedges;

    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 0; iDeltaRw <= 25; iDeltaRw++) {
        binedges.push_back((minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 24; iDeltaRw >= 0; iDeltaRw--) {
        binedges.push_back(1. - (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    for(int i = 0; i < binedges.size(); ++i)
        cout<<binedges.at(i)<<", ";
        
    cout<<endl;

  for (int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++) {
    deltaRBinsEEC[iDeltaR] =
      (minDeltaREEC + binnerShift) * TMath::Exp(iDeltaR * deltaRlogBinWidth) - binnerShift;
  }

  TFile *out = new TFile("../offline/pythia_pp_rhic_10pthat60_eec_star+fullevent_feb5.root", "RECREATE");
    
    if (!out || out->IsZombie()) {
        std::cerr << "Error: Could not open output" << std::endl;
        return 1;
    }
    
  out->cd();

  std::cout << "created main511.root!" << std::endl;
  
  Double_t eecbounds[51] = {0.0001, 0.000144062, 0.000205774, 0.000292206, 0.000413259, 
                            0.000582802, 0.000820258, 0.00115283, 0.00161862, 0.00227099,
                            0.00318468, 0.00446436, 0.00625663, 0.00876683, 0.0122825,
                            0.0172065, 0.0241028, 0.0337616, 0.0472893, 0.0662357
                            , 0.0927715, 0.129937, 0.181989, 0.254891, 0.356996,
                            0.5, 0.642984, 0.745089, 0.817991, 0.870043, 0.907208,
                            0.933744, 0.952691, 0.966218, 0.975877, 0.982774,
                            0.987697, 0.991213, 0.993723, 0.995516, 0.996795,
                            0.997709, 0.998361, 0.998827, 0.99916, 0.999397,
                            0.999567, 0.999688, 0.999774, 0.999836, 0.99988};
  // Initialize histogram
  TH1::SetDefaultSumw2();
 //TH2::SetDefaultSumw2();

  
  TH1F EEC_w("EEC_w", "Energy Energy Correlator", 50, eecbounds);
  TH1F EEC_w_p("EEC_w_p", "Energy Energy Correlator", 50, 0, 50);
  
  TH1F EEC_w_b("EEC_w_b", "Energy Energy Correlator", 100, eecbounds);
  TH1F EEC_w_pb("EEC_w_pb", "Energy Energy Correlator", 50, 0, 50);
  
  TH1F EEC_w_ub2("EEC_w_ub2", "Energy Energy Correlator", 50, eecbounds);
  TH1F EEC_w_pub2("EEC_w_pub2", "Energy Energy Correlator", 50, 0, 50);
  
  TH1F EEC_w_b2("EEC_w_b2", "Energy Energy Correlator", 50, eecbounds);
  TH1F EEC_w_pb2("EEC_w_pb2", "Energy Energy Correlator", 50, 0, 50);
  
  TH1F EEC_w_ub("EEC_w_ub", "Energy Energy Correlator", 50, eecbounds);
  TH1F EEC_w_pub("EEC_w_pub", "Energy Energy Correlator", 50, 0, 50);
  
  
  TH1F EEC_star("EEC_star", "Energy Energy Correlator",nDeltaRBinsEEC, deltaRBinsEEC);
  TH1F Aj_spectrum("Aj_spectrum", "Jet Assymetry spectrum", 60, 0.0, 1.0);
  
  // Initialize Pythia for 200 GeV pp collision
  Pythia pythia;
 
  if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
  }
  
  pythia.init();

  float jet_radius = 0.4;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_radius);
  int dijet_event_counter = 0;

  for (int iEvent = 0; iEvent < 100000; ++iEvent) { // 10000 events for now
    if (!pythia.next()) continue;

    std::vector<fastjet::PseudoJet> event;

    // pT cut > .5 GeV
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      if (p.isFinal() && p.pT() > .5) {
	    event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
      }
    }//! pseudojet of all particles in event loop close
    
    // Perform jet clustering
    fastjet::ClusterSequence cs(event, jet_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

  // done here so there is no dijet cut 
  for (size_t i = 0; i < jets.size(); ++i) {
  // star cuts from paper need to add charged particle
  if (jets[i].eta() < -0.6 || jets[i].eta() > 0.6) continue;
  if (jets[i].pt() < 15.0 || jets[i].pt() > 20.0) continue;
    
  vector<PseudoJet> constituents = jets[i].constituents();
  for (size_t j = 0; j < constituents.size(); ++j) {
    for (size_t k = j + 1; k < constituents.size(); ++k) {
      //float eec = constituents.at(j).pt() * constituents.at(k).pt();
      if (constituents.at(j).pt() < 2 || constituents.at(k).pt() < 2) continue; // constituent pt cut of 2 GeV 
      float eec = (constituents.at(j).pt() * constituents.at(k).pt()) / (jets[i].pt() * jets[i].pt()); 
      float weight = 1 / eec ; // normalization of pairs 
      float delr = deltaR( constituents.at(j),  constituents.at(k));
      EEC_star.Fill(delr, (weight * eec));
      }	
    }   
  } //!EEC star outer loop close 

    // Require atleast two jets
  if (jets.size() < 2) continue;
    
    //! select di-jets on jet pT 
  if(jets[0].pt() < 15 || jets[1].pt() < 15) continue;
    
    // Check for dijet condition
  float dphi = std::abs(jets[0].phi() - jets[1].phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  if (dphi < (7.0 * M_PI / 8.0)) continue;
    
  dijet_event_counter++;
    
    //jet asymmetry  
  float aj = (jets[0].pt()-jets[1].pt()) / (jets[0].pt() + jets[1].pt());
  Aj_spectrum.Fill(aj);
    
  for (size_t i = 0; i < event.size(); ++i) {
    for (size_t j = i + 1; j < event.size(); ++j) {
      float eec = event.at(i).pt() * event.at(j).pt();
      float delr = deltaR(event.at(i), event.at(j));
      float z = (1 - TMath::Cos(delr))/2 ; 
      EEC_w.Fill(z, eec);
      // filling balanced and unbalanced whole EECs
      if (aj < 0.667) EEC_w_b.Fill(z, eec);
      if (aj > 0.667) EEC_w_ub.Fill(z, eec);
      if (aj < 0.1) EEC_w_b2.Fill(z, eec);
      if (aj > 0.1) EEC_w_ub2.Fill(z, eec);
    } 
  }//!Whole event EEC loop close 

}//! event loop 
  
  out->Write();
  
  // Pythia cleanup
  pythia.stat();
  out->Close();
  return 0;
}

