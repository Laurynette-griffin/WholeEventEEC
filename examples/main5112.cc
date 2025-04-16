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

const Int_t nDeltaRBinsEECw = 30;
const Float_t minDeltaREECw = 1e-5; // Avoid log(0) issue
const Float_t maxDeltaREECw = 0.5;
const Float_t binnerShiftw = 1e-7;
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
    for (int iDeltaRw = 0; iDeltaRw <= 30; iDeltaRw++) {
        binedges.push_back((minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 29; iDeltaRw >= 0; iDeltaRw--) {
        binedges.push_back(1. - (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
   for(int i = 0; i < binedges.size(); ++i){ 
        cout<<binedges.at(i)<<", " <<endl;
    }

    for (int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++) {
        deltaRBinsEEC[iDeltaR] = (minDeltaREEC + binnerShift) * TMath::Exp(iDeltaR * deltaRlogBinWidth) - binnerShift;
    }

    TFile *out = new TFile("../offline/pythia_wholeeventEEC_ee_March5.root", "RECREATE");
    
    if (!out || out->IsZombie()) {
        std::cerr << "Error: Could not open output" << std::endl;
        return 1;
    }
    
    out->cd();

    std::cout << "created main5112.root!" << std::endl;
  
    Double_t eecbounds[61] = {1e-05, 1.43814e-05, 2.06634e-05, 2.96705e-05, 4.25849e-05, 
                                6.11016e-05, 8.76508e-05, 0.000125717, 0.000180296, 0.000258552, 
                                0.000370755, 0.000531632, 0.000762296, 0.00109302, 0.00156722, 
                                0.00224712, 0.00322196, 0.00461969, 0.00662375, 0.00949717, 
                                0.0136171, 0.0195242, 0.0279938, 0.0401376, 0.0575492, 
                                0.0825141, 0.118309, 0.169631, 0.243217, 0.348724, 
                                0.5, 
                                0.651276, 0.756783, 0.830369, 0.881691, 0.917486, 
                                0.942451,  0.959862, 0.972006, 0.980476, 0.986383, 
                                0.990503, 0.993376, 0.99538, 0.996778, 0.997753, 
                                0.998433, 0.998907, 0.999238, 0.999468, 0.999629, 
                                0.999741, 0.99982, 0.999874, 0.999912, 0.999939, 
                                0.999957, 0.99997, 0.999979, 0.999985, 0.99999};
    // Initialize histogram
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    TH1F EEC_w("EEC_w", "Energy Energy Correlator", 60, eecbounds);
    TH1F EEC_w_p("EEC_w_p", "Whole Event EEC", 60, 0, 60);
    
    TH1F Aj_spectrum("Aj_spectrum", "Jet Assymetry spectrum", 60, 0.0, 1.0);
    TH1F Costheta_spectrum("Costheta_spectrum", "Cosine theta", 100, -1., 1.);
    TH1F JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 60, 0.0, 60);
    TH2F etaphi_spectrum("etaphi_spectrum", "Eta Phi Distribution", 60 , 0.0, 6.283, 60, -5.0, 5.0);
    
    // Initialize Pythia for 91.2 GeV e+e- collision
    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
  
    pythia.init();


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
     float costheta_p = std::cos(p.theta());
     if (!p.isFinal() || p.pT() < .2 || !p.isCharged()|| std::abs(costheta_p >= .94)) continue;
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
            float dr = deltaR(event.at(i), event.at(j));
            Costheta_spectrum.Fill(ctheta);
            EEC_w.Fill(dr, eec / q2);
        }
    }
        
    
    
    //cout << "Total e for " << iEvent << "is" << totale << endl;
    //Perform jet clustering
    fastjet::ClusterSequence cs(event, jet_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    
    for (size_t i = 0; i < jets.size(); ++i) {
        if (jets[i].pt() < 5) continue;
        JetSpectrum.Fill(jets[i].pt());
    }
    
    //jet asymmetry 
    float aj = (jets[0].e()-jets[1].e()) / (jets[0].e() + jets[1].e());
    Aj_spectrum.Fill(aj);
   

    }//! event loop 
cout << " total # of dijet events = " << dijet_event_counter << endl;

  out->Write();
 
  // Pythia cleanup
  pythia.stat();
  out->Close();
  return 0;
}

