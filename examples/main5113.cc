#include <iostream>
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

// 8 GeV pt cut pthatmin = 5 in .cmnd

// Function to compute Delta R
float deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  float dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  float deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}


float costheta(const PseudoJet& p1, const PseudoJet& p2) {
    float dotprod = (p1.px() * p2.px() +  p1.py() * p2.py() +  p1.pz() * p2.pz());
    float normp1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz()); 
    float normp2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz()); 
    
    if (normp1 == 0 || normp2 == 0) // division by 0 check 
        return 0; 
    
    return (dotprod / (normp1 * normp2));
}


float deltaphitheta(const PseudoJet& p1, const PseudoJet& p2) {
    float eta1 = p1.eta();
    float eta2 = p2.eta();
    float theta1 = 2 * std::atan(std::exp(-1*eta1));
    float theta2 = 2 * std::atan(std::exp(-1*eta2));
    // to make sure the value of delta theta is not greater than 180 
    float dtheta = std::acos(std::cos(theta1 - theta2));
    
    float dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    
    return std::sqrt(dtheta * dtheta + dphi * dphi);
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
    (TMath::Log(maxDeltaREEC + binnerShift) - TMath::Log(minDeltaREEC + binnerShift)) / nDeltaRBinsEEC;
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
    for (int iDeltaRw = 0; iDeltaRw <= 30; iDeltaRw++) {
        binedges.push_back((minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    //!   print the values of the bin edges for the left side (1e-4 to 0.5) 
    for (int iDeltaRw = 29; iDeltaRw >= 0; iDeltaRw--) {
        binedges.push_back(1. - (minDeltaREECw + binnerShiftw) * TMath::Exp(iDeltaRw * deltaRlogBinWidthw) - binnerShiftw);
    }    
    
    for(int i = 0; i < binedges.size(); ++i){
        cout<<binedges.at(i)<<", ";
        cout<<endl;}

    for (int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++) {
        deltaRBinsEEC[iDeltaR] =
        (minDeltaREEC + binnerShift) * TMath::Exp(iDeltaR * deltaRlogBinWidth) - binnerShift;
    }

    TFile *out = new TFile("../offline/pythia_pp_rhic_15GeV_fullevent_5pthat60_Apr15.root", "RECREATE");
    
    if (!out || out->IsZombie()) {
        std::cerr << "Error: Could not open output" << std::endl;
        return 1;
    }
    
    out->cd();

    std::cout << "created output.root!" << std::endl;

                                
    Double_t tryit[101] =   {1e-05, 0.0444099, 0.0888098, 0.1332097, 0.1776096,
                            0.22200950000000003, 0.2664094000000001, 0.31080930000000007, 0.35520920000000006, 0.39960910000000005,
                            0.44400900000000004, 0.4884089000000001, 0.5328088000000001, 0.5772087, 0.6216086000000001,
                            0.6660085, 0.7104084, 0.7548083000000001, 0.7992082, 0.8436081000000001,
                            0.888008, 0.9324079000000001, 0.9768078000000001, 1.0212077000000002, 1.0656076000000003,
                            1.1100075000000003, 1.1544074000000002, 1.1988073000000001, 1.2432072000000003, 1.2876071000000002,
                            1.3320070000000002, 1.3764069000000003, 1.4208068000000003, 1.4652067000000002, 1.5096066000000004,
                            1.5540065000000003, 1.5984064000000002, 1.6428063000000004, 1.6872062000000003, 1.7316061000000003,
                            1.7760060000000002, 1.8204059000000004, 1.8648058000000003, 1.9092057000000002, 1.9536056000000004,
                            1.9980055000000003, 2.0424054000000003, 2.0868053000000004, 2.1312052000000006, 2.1756051000000003,
                            2.2200050000000005, 2.2644049, 2.3088048000000003, 2.3532047000000005, 2.3976046,
                            2.4420045000000004, 2.4864044000000005, 2.5308043000000002, 2.5752042000000004, 2.6196041000000005,
                            2.6640040000000003, 2.7084039000000004, 2.7528038000000006, 2.7972037000000003, 2.8416036000000005,
                            2.8860035000000006, 2.9304034000000003, 2.9748033000000005, 3.0192032000000006, 3.0636031000000004,
                            3.1080030000000005, 3.1524029000000007, 3.1968028000000004, 3.2412027000000005, 3.2856026000000007,
                            3.3300025000000004, 3.3744024000000006, 3.4188023000000007, 3.4632022000000005, 3.5076021000000006,
                            3.5520020000000003, 3.5964019000000005, 3.6408018000000006, 3.6852017000000004, 3.7296016000000005,
                            3.7740015000000007, 3.8184014000000004, 3.8628013000000005, 3.9072012000000007, 3.9516011000000004,
                            3.9960010000000006, 4.0404009, 4.0848008, 4.1292007, 4.1736006,
                            4.2180005000000005, 4.262400400000001, 4.3068003, 4.3512002, 4.3956001,
                            4.44};
    
    Double_t eecbounds[29] = {1e-05, 9.829499e-05, 9.661906e-04, 0.00949717, 
                                0.0136171, 0.0195242, 0.0279938, 0.0401376, 0.0575492, 
                                0.0825141, 0.118309, 0.169631, 0.243217, 0.348724, 
                                0.5, 
                                0.651276, 0.756783, 0.830369, 0.881691, 0.9174859, 
                                0.9424508, 0.9598624, 0.9720062, 0.9804758, 0.9863829, 
                                0.99050283, 0.9990338094, 0.9999017, 0.99999};

    // Initialize histogram'
    
    int bins = 28;
    TH1::SetDefaultSumw2();
    //TH2::SetDefaultSumw2();
    
    
    TH1F EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_w_p("EEC_w_p", "Energy Energy Correlator", bins, 0, bins);
    
    TH1F EEC_w_low("EEC_w_low", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_w_l("EEC_w_l", "Energy Energy Correlator", bins, 0, bins);

    TH1F EEC_w_mid("EEC_w_mid", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_w_m("EEC_w_m", "Energy Energy Correlator", bins, 0, bins);

    TH1F EEC_w_high("EEC_w_high", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_w_h("EEC_w_h", "Energy Energy Correlator", bins, 0, bins);

    TH1F JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 60, 0, 70);
    TH1F LeadingJetSpectrum("LeadingJetspectrum", "Leading Jet Spectrum", 60 , 0, 70);
    TH1F SubleadingJetSpectrum("SubleadingJetSpectrum", "Subleading Jet Spectrum ", 60, 0, 70);
    
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

    for (int iEvent = 0; iEvent < 3000000; ++iEvent) { // 2M events (ran March 4th @ 3pm)
        if (!pythia.next()) continue;
    
        
        std::vector<fastjet::PseudoJet> event;
        std::vector<fastjet::PseudoJet> charged_event;

       // pT cut > .2 GeV
        for (int i = 0; i < pythia.event.size(); ++i) {
        Particle& p = pythia.event[i];
    
        // charged particles for eec
        if (!p.isFinal() || p.pT() < .2 || !p.isCharged() || abs(p.eta()) > 1.1 ) continue;
        charged_event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
    
        
        // all particles for jet finding 
        if (!p.isFinal() || p.pT() < .2  || abs(p.eta()) > 1.1 ) continue;
        event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        }
        
        // Perform jet clustering
        fastjet::ClusterSequence cs(event, jet_def);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

        // Require atleast two jets
        if (jets.size() < 2) continue;
        
        //erases (remove if statement) detector acceptance cut jets 
        jets.erase( std::remove_if(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet& jet) { 
                return std::abs(jet.eta()) > 0.7; }),  
            jets.end());
        
        //! select di-jets on jet pT 
        if(jets[0].pt() < 20.9 || jets[0].pt() >= 60.8) continue;
        if(jets[1].pt() < 9.4) continue;
        
        
        // Check for dijet condition
        float dphi = std::abs(jets[0].phi() - jets[1].phi());
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        if (dphi < (3.0 * M_PI / 4.0)) continue;
    
        dijet_event_counter++;
        
        // Jet spectra of dijet events that pass the cut 
         for (size_t i = 0; i < jets.size(); ++i) {
            if (jets[i].pt() < 5) continue;
            JetSpectrum.Fill(jets[i].pt());
        }
        
        //jet spectra leading and subleading 
        LeadingJetSpectrum.Fill(jets[0].pt());
        SubleadingJetSpectrum.Fill(jets[1].pt());
    
        for (size_t i = 0; i < charged_event.size(); ++i) {
            for (size_t j = i + 1; j < charged_event.size(); ++j) {

                float eec = charged_event.at(i).pt() * charged_event.at(j).pt();  
                float ctheta = costheta(charged_event.at(i), charged_event.at(j));
                float pmq2 = ((jets[0].pt() + jets[1].pt())/2) *  ((jets[0].pt() + jets[1].pt())/2); // average of leading jets pt squared 
                float z = (1 - ctheta)/2; 
                EEC_w.Fill(z, eec/pmq2);
                
                // work PLEASE! 
                if (jets[0].pt() < 31.2) EEC_w_low.Fill(z, eec);
                if (31.2<= jets[0].pt() < 40.7) EEC_w_mid.Fill(z, eec);
                if (40.7<= jets[0].pt() < 60.8) EEC_w_high.Fill(z, eec);

                } 
            }//!Whole event EEC loop close 

    }//! event loop 
    
    cout << " total # of dijet events = " << dijet_event_counter << endl;
  
    out->Write();
  
    // Pythia cleanup
    pythia.stat();
    out->Close();
    return 0;
    }

