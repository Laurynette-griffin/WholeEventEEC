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

const Int_t nDeltaRBinsEECw = 50;
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
    for (int iDeltaRw = 0; iDeltaRw <= 50; iDeltaRw++) {
        binedges.push_back((minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 49; iDeltaRw >= 0; iDeltaRw--) {
        binedges.push_back(1. - (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    for(int i = 0; i < binedges.size(); ++i)
        cout<<binedges.at(i)<<", ";
        
    cout<<endl;

  for (int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++) {
    deltaRBinsEEC[iDeltaR] =
      (minDeltaREEC + binnerShift) * TMath::Exp(iDeltaR * deltaRlogBinWidth) - binnerShift;
  }

 TFile *out = new TFile("pythia_wholeeventEEC_starEEC_Feb17.root", "RECREATE");
 
  out->cd();

  std::cout << "created main5111.root!" << std::endl;
  
  Double_t eecbounds[101] = {9.99998e-05, 0.000924117, 0.00181548, 0.00277957, 0.00382233, 
                              0.00495017, 0.00617004, 0.00748945, 0.00891652, 0.01046,
                              0.0121295, 0.0139351, 0.0158882, 0.0180005, 0.0202852, 
                              0.0227564, 0.0254292, 0.02832, 0.0314468, 0.0348287,
                              0.0384865, 0.0424428, 0.0467219, 0.0513502, 0.0563561,
                              0.0617705, 0.0676266, 0.0739606, 0.0808115, 0.0882213, 
                              0.0962357, 0.104904, 0.11428, 0.124421, 0.135389, 
                              0.147252, 0.160083, 0.173961, 0.188971, 0.205207,
                              0.222766, 0.241759, 0.262302, 0.28452, 0.308552,
                              0.334545, 0.362658, 0.393065, 0.425954, 0.461526, 
                              0.5, 0.518474, 0.554046, 0.586935, 0.617342, 0.645455, // this line has 6 
                              0.671448, 0.69548, 0.717698, 0.738241, 0.757234, 
                              0.774793, 0.791029, 0.806039, 0.819917, 0.832748, 
                              0.844611, 0.855579, 0.86572, 0.875096, 0.883764, 
                              0.891779, 0.899189, 0.906039, 0.912373, 0.91823, 
                              0.923644, 0.92865, 0.933278, 0.937557, 0.941513,
                              0.945171, 0.948553, 0.95168, 0.954571, 0.957244, 
                              0.959715, 0.961999, 0.964112, 0.966065, 0.967871,
                              0.96954, 0.971083, 0.972511, 0.97383, 0.97505,
                              0.976178, 0.97722, 0.978185, 0.979076, 0.9799};
  // Initialize histogram
  TH1::SetDefaultSumw2();
 //TH2::SetDefaultSumw2();

  
  TH1F EEC_w("EEC_w", "Energy Energy Correlator", 100, eecbounds);
  TH1F EEC_w_p("EEC_w_p", "Energy Energy Correlator", 100, 0, 100);
  
  TH1F EEC_w_b("EEC_w_b", "Energy Energy Correlator", 100, eecbounds);
  TH1F EEC_w_pb("EEC_w_pb", "Energy Energy Correlator", 100, 0, 100);
  
  TH1F EEC_w_ub2("EEC_w_ub2", "Energy Energy Correlator", 100, eecbounds);
  TH1F EEC_w_pub2("EEC_w_pub2", "Energy Energy Correlator", 100, 0, 100);
  
  TH1F EEC_w_b2("EEC_w_b2", "Energy Energy Correlator", 100, eecbounds);
  TH1F EEC_w_pb2("EEC_w_pb2", "Energy Energy Correlator", 100, 0, 100);
  
  TH1F EEC_w_ub("EEC_w_ub", "Energy Energy Correlator", 100, eecbounds);
  TH1F EEC_w_pub("EEC_w_pub", "Energy Energy Correlator", 100, 0, 100);
  
  
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
      if (aj < 0.167) EEC_w_b.Fill(z, eec);
      if (aj > 0.167) EEC_w_ub.Fill(z, eec);
      if (aj < 0.333) EEC_w_b2.Fill(z, eec);
      if (aj > 0.333) EEC_w_ub2.Fill(z, eec);
    } 
  }//!Whole event EEC loop close 

}//! event loop 
  
  out->Write();
  
  // Pythia cleanup
  pythia.stat();
  out->Close();
  return 0;
}

