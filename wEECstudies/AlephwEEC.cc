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
double deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  double deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}

// Function to calulate opening angle

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

//allows for the .cmnd file to be accepted as an argument 
int main(int argc, char* argv[]) {
  // Check if a .cmnd file is provided
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <config.cmnd>" << endl;
    return 1;
  } 
    string configFile = argv[1];
    string outputFileName = argv[2]; 
    
    TFile *out = new TFile(outputFileName.c_str(), "RECREATE");
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

    Double_t topi[61] =   {1.00000000e-05, 1.49006082e-05, 2.22028125e-05, 3.30835411e-05,
                                4.92964884e-05, 7.34547660e-05, 1.09452069e-04, 1.63090240e-04,
                                2.43014377e-04, 3.62106202e-04, 5.39560265e-04, 8.03977611e-04,
                                1.19797554e-03, 1.78505642e-03, 2.65984263e-03, 3.96332730e-03,
                                5.90559873e-03, 8.79970130e-03, 1.31120901e-02, 1.95378118e-02,
                                2.91125279e-02, 4.33794373e-02, 6.46379999e-02, 9.63145513e-02,
                                1.43514539e-01, 2.13845393e-01, 3.18642641e-01, 4.74796916e-01,
                                7.07476283e-01, 1.05418269e+00, 1.57079633e+00, 2.08761365,
                                2.43432005, 2.66799942, 2.82315369, 2.92795094, 
                                2.99828179, 3.04548178, 3.07715833, 3.09841689, 
                                3.11268380, 3.12225852, 3.12868424, 3.13299663, 
                                3.13589073, 3.13783300, 3.13813649, 3.13981128, 
                                3.14059876, 3.14099235, 3.14115677, 3.14123423, 
                                3.14135331, 3.14143323, 3.14148688, 3.14155555, 
                                3.14156300, 3.14157005, 3.14157413, 3.14157740, 
                                3.14158265};
    // Initialize histogram'
    
    int bins = 60;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // whole event with thermal and mixed seperated out 
    TH1D EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_e("EEC_e", "wEEC in linear costheta", bins, 0.0, 1.0);
    TH1D EEC_d("EEC_d","wEEC in delphi", bins, topi);
    
    TH1D Aj_spectrum("Aj_spectrum", "Jet Assymetry spectrum", 60, 0.0, 1.0);
    TH1D Costheta_spectrum("Costheta_spectrum", "Cosine theta", 100, -1., 1.);
    TH1D JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 60, 0.0, 60);
    TH2D etaphi_spectrum("etaphi_spectrum", "Eta Phi Distribution", 60 , 0.0, 6.283, 60, -5.0, 5.0);
    
    // Initialize Pythia for 91.2 GeV e+e- collision
    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
  
    pythia.init();
    
    Sphericity sph;

    for (int iEvent = 0; iEvent < 1000000; ++iEvent) { // 10000 events for now
        if (!pythia.next()) continue;
    
     
    // implementing event sphericity cut from lunas paper 
    sph.analyze(pythia.event);
    
    Vec4 sphericityaxis = sph.eventAxis(3);
    
    double thetaS = std::acos(sphericityaxis.pz() / sphericityaxis.pAbs());
    
    
    if (thetaS <= (7*M_PI)/36 ||  (29*M_PI)/36 >= thetaS) continue; // added cut from luna's paper 
    
    std::vector<fastjet::PseudoJet> event; // pseudojets that will hold all particles passing the cuts
    
    double totale = 0;
    
    for (size_t i = 0; i < pythia.event.size(); ++i) {
     Particle& p = pythia.event[i];
     double costheta_p = std::cos(p.theta());
     if (!p.isFinal() || p.pT() < .2 || !p.isCharged()|| std::abs(costheta_p) >= .94) continue;
     totale =  p.e() + totale;
     event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
    }
    
    if (totale < 15.0) continue;

    double q2 = 91.2*91.2;
 
    for (size_t i = 0; i < event.size(); ++i) {
        double eta = event.at(i).eta();
        double phi = event.at(i).phi();
        
        etaphi_spectrum.Fill(phi,eta);

        for (size_t j = i + 1; j < event.size(); ++j) {
            double eec = event.at(i).pt() * event.at(j).pt();
            double ctheta = costheta(event.at(i), event.at(j));
            double z = (1 - ctheta) / 2;
            double delphi = deltaphi(event.at(i), event.at(j));
            Costheta_spectrum.Fill(ctheta);
            EEC_w.Fill(z, eec / q2);
            EEC_e.Fill(ctheta, eec / q2);
            EEC_d.Fill(delphi, eec / q2);
            
        }
    }
        
    
  
    }//! event loop 

  out->Write();
 
  // Pythia cleanup
  pythia.stat();
  out->Close();
  return 0;
}

