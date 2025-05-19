#include <iostream>
#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "Pythia8/Pythia.h"
#include "TRandom3.h"
#include "TLine.h"

using namespace fastjet;
using namespace Pythia8;
using namespace std;

// 8 GeV pt cut pthatmin = 5 in .cmnd

// Function to compute Delta R
double deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  double deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}


double costheta(const PseudoJet& p1, const PseudoJet& p2) {
    double dotprod = (p1.px() * p2.px() +  p1.py() * p2.py() +  p1.pz() * p2.pz());
    double normp1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz()); 
    double normp2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz()); 
    
    if (normp1 == 0 || normp2 == 0) // division by 0 check 
        return 0; 
    
    return (dotprod / (normp1 * normp2));
}

double deltaphi(const PseudoJet& p1, const PseudoJet& p2) {
    double dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    return dphi;
}


const Int_t nDeltaRBinsEECw = 30;
const double_t minDeltaREECw = 1e-5; // Avoid log(0) issue
const double_t maxDeltaREECw = 0.5;
const double_t binnerShiftw = 1e-7;
const double_t deltaRlogBinWidthw =
  (TMath::Log(maxDeltaREECw + binnerShiftw) - TMath::Log(minDeltaREECw + binnerShiftw)) / nDeltaRBinsEECw;
double_t deltaRBinsEECw[nDeltaRBinsEECw + 1];

const Int_t nDeltaRBinsEEC = 13;
const double_t minDeltaREEC = .025; // Avoid log(0) issue                             
const double_t maxDeltaREEC = 1;
const double_t binnerShift = 0.1;
const double_t deltaRlogBinWidth = 
    (TMath::Log(maxDeltaREEC + binnerShift) - TMath::Log(minDeltaREEC + binnerShift)) / nDeltaRBinsEEC;
double_t deltaRBinsEEC[nDeltaRBinsEEC + 1];

    //lol lets try thermal particle production :)
    
vector<PseudoJet> thermalGen(int nParticles = 1500){
    
    TF1 *mb = new TF1("mb","exp(-[0]*x/[1])",0.3,10); // pt range from .3 to 10 Gev
    mb->SetParameters(1,0.26); // shape and temperature (260 MeV)
    
    double mass = 0.13957; //charged pions 
    
    TRandom3 *r = new TRandom3(0); // be sure to set seed to 0 which is a time dependant seed if you pass nothing it is a default so each one is the same
    
    vector<PseudoJet> thermalParticles;
    for(int i=0; i<nParticles; i++)
        {
        	double pT = mb->GetRandom();
        	
        	double eta = 2.2*r->Rndm() - 1.1; 
        	double phi = 2.0*TMath::Pi()*r->Rndm() - TMath::Pi();
            
        	double p = pT*cosh(eta);
        	double pz = pT*sinh(eta);
        
        	double px = pT*cos(phi);
        	double py = pT*sin(phi);
        
        	double E = sqrt(mass*mass + p*p);
        
        	PseudoJet tmpPart(px, py, pz, E);
        	tmpPart.set_user_index(1); // to keep track that this is a thermal particle 
        	thermalParticles.push_back(tmpPart);
        }

return thermalParticles;
}


    

//allows for the .cmnd file to be accepted as an argument 
int main(int argc, char* argv[]) {
    // Check for required arguments
    if (argc < 1) {
        cerr << "Usage: " << argv[0] << " <config.cmnd> <output.root>" << endl;
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
    
    Double_t topi[61] = {1.00000000e-05, 1.49006082e-05, 2.22028125e-05, 3.30835411e-05,
                            4.92964884e-05, 7.34547660e-05, 1.09452069e-04, 1.63090240e-04,
                            2.43014377e-04, 3.62106202e-04, 5.39560265e-04, 8.03977611e-04,
                            1.19797554e-03, 1.78505642e-03, 2.65984263e-03, 3.96332730e-03,
                            5.90559873e-03, 8.79970130e-03, 1.31120901e-02, 1.95378118e-02,
                            2.91125279e-02, 4.33794373e-02, 6.46379999e-02, 9.63145513e-02,
                            1.43514539e-01, 2.13845393e-01, 3.18642641e-01, 4.74796916e-01,
                            7.07476283e-01, 1.05418269e+00, 1.57079633e+00, 
                            2.0874099625242373, 2.4341163709009024, 2.6667957376383757, 2.8229500122560345, 
                            2.927747261033487, 2.9980781141479063, 3.0452781022917392, 3.0769546536504087, 
                            3.098213216320398, 3.1124801256716323, 3.122054841779965, 3.1284805634480732, 
                            3.132792952294784, 3.1356870548603424, 3.137629326292187, 3.1389328109578147, 
                            3.139807597172383, 3.1403946780498013, 3.1407886759785417, 3.14105309332515, 
                            3.1412305473879445, 3.141349639213069, 3.1414295633499876, 3.1414832015208316, 
                            3.1415191988238047, 3.1415433571014066, 3.141559570048725, 3.1415704507772717, 
                            3.141577752981577, 3.141582653589793};
                            
                            
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

    // Initialize histogram'
    
    int bins = 60;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // whole event with thermal and mixed seperated out 
    TH1F EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_td("EEC_t", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_m("EEC_m", "Energy Energy Correlator", bins, eecbounds);
    TH1F EEC_s("EEC_s", "Energy Energy Correlator", bins, eecbounds);
    
    //phi + thermal components
    TH1F EEC_w_phi("EEC_w_phi", "Energy Energy Correlator", bins, topi);
    TH1F EEC_tphi("EEC_w_lowphi", "Energy Energy Correlator", bins, topi);
    TH1F EEC_mphi("EEC_w_midphi", "Energy Energy Correlator", bins, topi);
    TH1F EEC_sphi("EEC_w_hiphi", "Energy Energy Correlator", bins, topi);
  

    TH1F JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 70, 0, 70);
    TH1F LeadingJetSpectrum("LeadingJetspectrum", "Leading Jet Spectrum", 70 , 0, 70);
    TH1F SubleadingJetSpectrum("SubleadingJetSpectrum", "Subleading Jet Spectrum ", 70, 0, 70);
    //TH1F ThermalPhi("ThermalPhi", "phi distrubution", 31, 0, 2*TMath::Pi());
    TH1F Deltaphi("Deltaphi", "linear del phi", 31, 0, TMath::Pi());
    //TH2F eectvdelphit("ptvphi", "eect vs delPhi",10, 0, 10, 32, -2*TMath::Pi(), 2*TMath::Pi()); 
    
    TH3F ptvetaphit("ptvetaphi", "pt vs delphi", 32, -2*TMath::Pi(), 2*TMath::Pi(), 22, -1.1, 1.1, 10, 0, 10); 
    
    // reco vs truth q2
    TH2F qvspmq("comparisonpthat", "pthat^2 vs Average leading and subleading jet pt^2", 12, 0, 60, 10, 0, 60);
    
    for (int iEvent = 0; iEvent < 1000; ++iEvent) { // Simulate 100 events
    
    vector<PseudoJet> thermalParticles = thermalGen(1500);
    
    for (size_t i = 0; i < thermalParticles.size(); ++i) {
        for (size_t j = i + 1; j < thermalParticles.size(); ++j) {
            double eect = thermalParticles[i].pt()*thermalParticles[j].pt();
            double cthetat = costheta(thermalParticles[i], thermalParticles[j]);
            double delphit = deltaphi(thermalParticles[i], thermalParticles[j]);
            double zt = (1.0 - cthetat) / 2.0;
            double bigdp = TMath::Pi() - delphit;
            EEC_td.Fill(delphit, eect);
            // Histogram filling logic
            if (delphit < TMath::Pi() / 2) {
                EEC_td.Fill(delphit, eect *delphit);
                EEC_m.Fill(delphit, eect );
            } else {
                EEC_td.Fill(delphit, eect * bigdp);
               EEC_m.Fill(delphit, eect *(1/bigdp));
            }
        }
     }
        /*
        //cout << "Initial Charged event size = " << charged_event.size() << endl;
        
        charged_event.insert(charged_event.end(), thermalParticles.begin(), thermalParticles.end());
        
        for (size_t i = 0; i < charged_event.size(); ++i) {
            for (size_t j = i + 1; j < charged_event.size(); ++j) {

                double eec = charged_event.at(i).pt() * charged_event.at(j).pt();  
                double ctheta = costheta(charged_event.at(i), charged_event.at(j));
                double avejet = jets[0].pt() + jets[1].pt();
                double pmq2 = (avejet*avejet) / 4; // average of leading jets pt squared 
                double z = (1 - ctheta)/2; 
                double delphi = deltaphi(charged_event.at(i), charged_event.at(j));
                EEC_w.Fill(z, eec/pmq2);
                int ui = charged_event.at(i).user_index();
                int uj = charged_event.at(j).user_index();
                if ((ui == 1 || uj == 1) && (ui == 1 || uj == 1))  EEC_t.Fill(z, eec/pmq2);
                
                if ((ui == -1 || uj == -1) && (ui == 1 || uj == 1)) EEC_m.Fill(z, eec/pmq2);
                
                if ((ui == -1 || uj == -1) && (ui == -1 || uj == -1)) EEC_s.Fill(z, eec/pmq2);
                
                EEC_w_phi.Fill(delphi, eec/pmq2);
                
                if ((ui == 1 || uj == 1) && (ui == 1 || uj == 1))  EEC_tphi.Fill(delphi, eec/pmq2);
                
                if ((ui == -1 || uj == -1) && (ui == 1 || uj == 1)) EEC_mphi.Fill(delphi, eec/pmq2);
                
                if ((ui == -1 || uj == -1) && (ui == -1 || uj == -1)) EEC_sphi.Fill(delphi, eec/pmq2);

            }
        }//!Whole event EEC loop close 
        */
    }//! event loop 
    
  
    out->Write();
  
    // Pythia cleanup
    out->Close();
    return 0;
    }

