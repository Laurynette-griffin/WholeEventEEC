#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"
#include "TH2.h"
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
    
    return (dotprod / (normp1 * normp2));
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
    deltaRBinsEECw[iDeltaRw] = (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw;
  }

    std::vector<double> binedges;

    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 0; iDeltaRw <= 25; iDeltaRw++) {
        binedges.push_back((minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 26; iDeltaRw >= 50; iDeltaRw--) {
        binedges.push_back(1. - (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
   for(int i = 0; i < binedges.size(); ++i){ 
        cout<<binedges.at(i)<<", " <<endl;
    }

    for (int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++) {
        deltaRBinsEEC[iDeltaR] = (minDeltaREEC + binnerShift) * TMath::Exp(iDeltaR * deltaRlogBinWidth) - binnerShift;
    }

    TFile *out = new TFile("../offline/pythia_wholeeventEEC_ee_Feb21.root", "RECREATE");
    
    if (!out || out->IsZombie()) {
        std::cerr << "Error: Could not open output" << std::endl;
        return 1;
    }
    
    out->cd();

    std::cout << "created main5112.root!" << std::endl;
  
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
    TH2::SetDefaultSumw2();
    
    TH1F EEC_w("EEC_w", "Energy Energy Correlator", 50, eecbounds);
    TH1F EEC_w_p("EEC_w_p", "Whole Event EEC", 50, 0, 50);
    
    TH1F Aj_spectrum("Aj_spectrum", "Jet Assymetry spectrum", 60, 0.0, 1.0);
    TH1F Costheta_spectrum("Costheta_spectrum", "Cosine theta", 100, -1.2, 1.2);
    TH1F JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 60, 0.0, 60);
    TH2F etaphi_spectrum("etaphi_spectrum", "Eta Phi Distribution", 60 , 0.0, 6.283, 60, -5.0, 5.0);
    
    // Initialize Pythia for 91.2 GeV e+e- collision
    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
  
    pythia.init();
    
    std::ofstream outFile("particle_momentum.txt"); // Open file
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    //float jet_radius = 0.4;
    fastjet::JetDefinition jet_def(ee_genkt_algorithm, .4, -1);
    
    int dijet_event_counter = 0;
    Sphericity sph;

    for (int iEvent = 0; iEvent < 1000000; ++iEvent) { // 10000 events for now
        if (!pythia.next()) continue;
    
     
    // implementing event sphericity cut from lunas paper 
    sph.analyze(pythia.event);
    
    Vec4 sphericityaxis = sph.eventAxis(3);
    
    double thetaS = std::acos(sphericityaxis.pz() / sphericityaxis.pAbs());
    
    
    if (thetaS <= (7*M_PI)/36 ||  (29*M_PI)/36 <= thetaS) continue; // added cut from luna's paper 
    
    std::vector<fastjet::PseudoJet> event; // pseudojets that will hold all particles passing the cuts
    
    float totale = 0;
    
    for (size_t i = 0; i < pythia.event.size(); ++i) {
        Particle& p = pythia.event[i];
         if (p.pT() < 0 || !p.isCharged() || !p.isVisible() || !p.isFinal()) continue;
         totale =  p.e() + totale;
         event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        }
        if (totale < 15.0) continue;
    
    float q2 = 91.2*91.2;
 
    for (size_t i = 0; i < event.size(); ++i) {
        float eta = event.at(i).eta();
        float phi = event.at(i).phi();
        
        etaphi_spectrum.Fill(phi,eta);

        for (size_t j = i + 1; j < event.size(); ++j) {
            float eec = event.at(i).pt() * event.at(j).pt();
            float ctheta = costheta(event.at(i), event.at(j));
            float z = (1 - ctheta) / 2;
            
            Costheta_spectrum.Fill(z);
            EEC_w.Fill(z, eec / q2);
        }
    }
        
    
    
    //cout << "Total e for " << iEvent << "is" << totale << endl;
    //Perform jet clustering
    fastjet::ClusterSequence cs(event, jet_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    
    for (size_t i = 0; i < jets.size(); ++i) {
        JetSpectrum.Fill(jets[i].pt());
    }
    
    //jet asymmetry  
    float aj = (jets[0].e()-jets[1].e()) / (jets[0].e() + jets[1].e());
    Aj_spectrum.Fill(aj);
   

    }//! event loop 
  out->Write();
  
  // Pythia cleanup
  pythia.stat();
  out->Close();
  return 0;
}

