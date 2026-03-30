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
TRandom3 *r = new TRandom3(0); // be sure to set seed to 0 which is a time dependant seed if you pass nothing it is a default so each one is the same


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
 
vector<PseudoJet> thermalGen(int nParticles = 1500){

    TF1 *mb = new TF1("mb","exp(-[0]*x/[1])",0.3,10); // pt range from .3 to 10 Gev
    mb->SetParameters(1,0.26); // shape and temperature (260 MeV)
    
    TF1 *flowPhi = new TF1("flowPhi","1+2*[0]*cos(2*x)",-TMath::Pi(),TMath::Pi());
    flowPhi->SetParameter(0,0.015); 
    double phi = flowPhi->GetRandom();
    
    double mass = 0.13957; //charged pions 
    
    vector<PseudoJet> thermalParticles;
    for(int i=0; i<nParticles; i++){
        double pT = mb->GetRandom();
        
        double eta = 5.0*r->Rndm() - 2.5; 
        double phi = flowPhi->GetRandom();
        
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

vector<PseudoJet> thermalGenne(int nParticles = 1500){

    TF1 *mb = new TF1("mb","exp(-[0]*x/[1])",0.3,10); // pt range from .3 to 10 Gev
    mb->SetParameters(1,0.26); // shape and temperature (260 MeV)
    
    TF1 *flowPhi = new TF1("flowPhi","1+2*[0]*cos(2*x)",-TMath::Pi(),TMath::Pi());
    flowPhi->SetParameter(0,0.015); 
    double phi = flowPhi->GetRandom();
    
    double mass = 0.13957; //charged pions 
    
    vector<PseudoJet> thermalParticles;
    for(int i=0; i<nParticles; i++){
        double pT = mb->GetRandom();
        
        double eta = 5.0*r->Rndm() - 2.5; 
        double phi = flowPhi->GetRandom();
        
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
    if (argc < 3) {
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
    
    Double_t topi[41] =   {1.00000000e-05, 1.81888815e-05, 3.30835411e-05, 6.01752609e-05,
    1.09452069e-04, 1.99081071e-04, 3.62106202e-04, 6.58630680e-04,
    1.19797554e-03, 2.17898352e-03, 3.96332730e-03, 7.20884906e-03,
    1.31120901e-02, 2.38494254e-02, 4.33794373e-02, 7.89023445e-02,
    1.43514539e-01, 2.61036895e-01, 4.74796916e-01, 8.63602485e-01,
    1.57079633e+00,  2.2779901689071522, 2.6667957376383757, 2.880555758264255, 
    2.9980781141479063, 3.0626903091318147, 3.098213216320398, 3.1177432281926607,
    3.1284805634480732, 3.1343838045285803, 3.137629326292187, 3.139413670074593,
    3.1403946780498013, 3.14093402290975, 3.1412305473879445, 3.1413935725184277,
    3.1414832015208316, 3.141532478328942, 3.141559570048725, 3.1415744647082806,
    M_PI};
    
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
    /*
    Double_t topi[31] = {1.00000000e-05, 2.22028125e-05, 4.92964884e-05, 1.09452069e-04, 2.43014377e-04,
                        5.39560265e-04, 1.19797554e-03, 2.65984263e-03, 5.90559874e-03, 1.31120902e-02,
                        2.91125280e-02, 6.46380000e-02, 1.43514540e-01, 3.18642642e-01, 7.07476284e-01,
                        1.57079633e+00,
                        2.43411637e+00, 2.82295001e+00, 2.99807811e+00, 3.07695465e+00, 3.11248013e+00,
                        3.12848056e+00, 3.13568705e+00, 3.13893281e+00, 3.14039468e+00, 3.14105309e+00,
                        3.14134964e+00, 3.14148320e+00, 3.14154336e+00, 3.14157045e+00, 3.14158265e+00
        
    };
    
    
    
    
    
    Double_t eecbounds[31] =   {1.00000000e-05, 2.05714387e-05, 4.23184092e-05, 8.70550563e-05, 1.79084776e-04,
                            3.68403150e-04, 7.57858283e-04, 1.55902353e-03, 3.20713570e-03, 6.59753955e-03,
                            1.35720881e-02, 2.79197379e-02, 5.74349177e-02, 1.18151889e-01, 2.43055435e-01,
                            5.00000000e-01,
                            7.56944565e-01, 8.81848111e-01, 9.42565082e-01, 9.72080262e-01, 9.86427912e-01,
                            9.93402460e-01, 9.96792864e-01, 9.98440976e-01, 9.99242142e-01, 9.99631597e-01,
                            9.99820915e-01, 9.99912945e-01, 9.99957682e-01, 9.99979429e-01, 1.00000000};
                            
      */
    // Initialize histogram'
    
    int bins = 40;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // whole event with thermal and mixed seperated out   
    /*
    TH1D EEC_w("EEC_w", "Energy Energy Correlator", bins, 0,1);
    TH1D EEC_t("EEC_t", "Energy Energy Correlator", bins, 0,1);
    TH1D EEC_m("EEC_m", "Energy Energy Correlator", bins, 0,1);
    TH1D EEC_s("EEC_s", "Energy Energy Correlator", bins, 0,1);
    
    TH1D EEC_mt("EEC_mt", "Energy Energy Correlator", bins, 0,1);
    TH1D EEC_ms("EEC_ms", "Energy Energy Correlator", bins, 0,1); // mix signal
    TH1D EEC_mm("EEC_mm", "Energy Energy Correlator", bins, 0,1); // mix1 mix1 
    TH1D EEC_m2("EEC_m2", "Energy Energy Correlator", bins, 0,1); // mix1 mix2
    
    //phi + thermal components
    TH1D EEC_wp("EEC_wp", "Energy Energy Correlator", bins, topi);
    TH1D EEC_tp("EEC_tp", "Energy Energy Correlator", bins, topi);
    TH1D EEC_mp("EEC_mp", "Energy Energy Correlator", bins, topi);
    TH1D EEC_sp("EEC_sp", "Energy Energy Correlator", bins, topi);
    
    TH1D EEC_msp("EEC_msp", "Energy Energy Correlator", bins, topi); // mix signal phi
    TH1D EEC_mmp("EEC_mmp", "Energy Energy Correlator", bins, topi); // mix1 mix1 phi 
    TH1D EEC_m2p("EEC_m2p", "Energy Energy Correlator", bins, topi); // mix1 mix2 phi
    
    // phi linear bins 
    TH1D EEC_wpl("EEC_wpl", "Energy Energy Correlator", bins, 0,  M_PI);
    TH1D EEC_tpl("EEC_tpl", "Energy Energy Correlator", bins, 0,  M_PI);
    TH1D EEC_mpl("EEC_mpl", "Energy Energy Correlator", bins, 0,  M_PI);
    TH1D EEC_spl("EEC_spl", "Energy Energy Correlator", bins, 0,  M_PI);
    
    TH1D EEC_mspl("EEC_mspl", "Energy Energy Correlator", bins, 0, M_PI); // mix signal phi
    TH1D EEC_mmpl("EEC_mmpl", "Energy Energy Correlator", bins, 0,  M_PI); // mix1 mix1 phi 
    TH1D EEC_m2pl("EEC_m2pl", "Energy Energy Correlator", bins, 0,  M_PI); // mix1 mix2 phi
*/
    
    TH1D EEC_eta50("EEC_eta50", "Energy Energy Correlator", bins, 0,7);
    TH1D EEC_eta20("EEC_eta20", "Energy Energy Correlator", bins, 0,1);
    TH1D EEC_eta25("EEC_eta25", "Energy Energy Correlator", bins, 0,1); 
    TH1D EEC_eta15("EEC_eta15", "Energy Energy Correlator", bins, 0,1); 
    TH1D EEC_eta10("EEC_eta10", "Energy Energy Correlator", bins, 0,1); 
    TH1D EEC_eta05("EEC_eta05", "Energy Energy Correlator", bins, 0,1); 
    //TH1D EEC_t3t3("EEC_t3t3", "Energy Energy Correlator", bins, eecbounds);
    
    TH2D etavsphi("etavsphi", " eta vs phi", 100, -5., 5., 64, 0, 2*M_PI);
    
    TH1D EEC_t1t1p("EEC_t1t1p", "Energy Energy Correlator", bins, 0, M_PI);
    TH1D EEC_t2t2p("EEC_t2t2p", "Energy Energy Correlator", bins, 0, M_PI); 
    //TH1D EEC_t3t3p("EEC_t3t3p", "Energy Energy Correlator", bins, 0, M_PI);
    
    TH1D EEC_t1t2p("EEC_t1t2p", "Energy Energy Correlator", bins, 0, M_PI); 
    //TH1D EEC_t1t3p("EEC_t1t3p", "Energy Energy Correlator", bins, 0, M_PI); 
    TH1D EEC_t2t3p("EEC_t2t3p", "Energy Energy Correlator", bins, 0, M_PI); 
    //TH1D EEC_t2t3p("EEC_t2t3p", "Energy Energy Correlator", bins, 0, M_PI); 
    
    
    TH1D JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 70, 0, 70);
    
    // Initialize Pythia for 200 GeV pp collision
    //Pythia pythia;
    /*
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
    */
    //pythia.init();
    
    double jet_radius = 0.4;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_radius);
    int dijet_event_counter = 0;
    
    long nEvt = 1000;
     for (int iEvent = 0; iEvent < nEvt; ++iEvent) {
       // if (!pythia.next()) continue;
        cout << iEvent << endl;
        vector<PseudoJet> m1 = thermalGen(1500);
        vector<PseudoJet> m2 = thermalGenne(1500);

        for (size_t j = 0; j < m1.size(); ++j) {
            auto &p1 = m1[j];
            double e1 = p1.eta();
            
            for (size_t k = j + 1; k < m1.size(); ++k) {
                auto &p2 = m1[k];
                double e2 = p2.eta();
                double mm_eec = p1.pt() * p2.pt();
                double mm_ctheta = costheta(p1, p2);
                double mm_z = (1 - mm_ctheta) / 2;
                double mm_delphi = deltaphi(p1, p2);
                double mm_delr = deltaR(p1,p2);
                
                EEC_eta50.Fill(mm_delr, mm_eec);
                if (e1 < fabs(2.5) && e2 < fabs(2.5)) EEC_eta25.Fill(mm_z, mm_eec);
                if (e1 < fabs(2.0) && e2 < fabs(2.0)) EEC_eta20.Fill(mm_z, mm_eec);
                if (e1 < fabs(1.5) && e2 < fabs(1.5)) EEC_eta15.Fill(mm_z, mm_eec);
                if (e1 < fabs(1.0) && e2 < fabs(1.0)) EEC_eta10.Fill(mm_z, mm_eec);
                if (e1 < fabs(0.5) && e2 < fabs(0.5)) EEC_eta05.Fill(mm_z, mm_eec);
                
            }
        }
        
    }//! event loop 
    
    out->Write();
    
    //pythia.stat();
    out->Close();
    return 0;
}
