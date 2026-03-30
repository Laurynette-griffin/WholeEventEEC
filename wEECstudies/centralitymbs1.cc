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

TRandom3 *r = new TRandom3(0);
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
    // be sure to set seed to 0 which is a time dependant seed if you pass nothing it is a default so each one is the same
    
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

    // Initialize histogram'
    
    int bins = 40;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // whole event with thermal and mixed seperated out   
    TH1D EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_t("EEC_t", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_m("EEC_m", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_s("EEC_s", "Energy Energy Correlator", bins, eecbounds);
    
    TH1D EEC_ms("EEC_ms", "Energy Energy Correlator", bins, eecbounds); // mixed + whole event 
    TH1D EEC_mm("EEC_mm", "Energy Energy Correlator", bins, eecbounds);
    //every pythia + thermal
    TH1D EEC_pt2("EEC_pt2","Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_pt3("EEC_pt3","Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_pt4("EEC_pt4","Energy Energy Correlator", bins, eecbounds);
 
    // every self 
    TH1D EEC_t1t1("EEC_t1t1", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_t2t2("EEC_t2t2", "Energy Energy Correlator", bins, eecbounds); // mix signal
    TH1D EEC_t3t3("EEC_t3t3", "Energy Energy Correlator", bins, eecbounds); // mix1 mix1 
    TH1D EEC_t4t4("EEC_t4t4", "Energy Energy Correlator", bins, eecbounds); // mix1 mix2
    
    // every mixed mixed
    TH1D EEC_t1t2("EEC_t1t2", "Energy Energy Correlator", bins, eecbounds); // mix signal
    TH1D EEC_t1t3("EEC_t1t3", "Energy Energy Correlator", bins, eecbounds); // mix1 mix1 
    TH1D EEC_t1t4("EEC_t1t4", "Energy Energy Correlator", bins, eecbounds); // mix1 mix1 
    TH1D EEC_t2t3("EEC_t2t3", "Energy Energy Correlator", bins, eecbounds); // mix1 mix2
    TH1D EEC_t2t4("EEC_t2t4", "Energy Energy Correlator", bins, eecbounds); // mix1 mix2
    TH1D EEC_t3t4("EEC_t3t4", "Energy Energy Correlator", bins, eecbounds); // mix1 mix2
    
    //phi + thermal components
    TH1D EEC_wp("EEC_wp", "Energy Energy Correlator", bins, topi);
    TH1D EEC_tp("EEC_tp", "Energy Energy Correlator", bins, topi);
    TH1D EEC_mp("EEC_mp", "Energy Energy Correlator", bins, topi);
    TH1D EEC_sp("EEC_sp", "Energy Energy Correlator", bins, topi);

    TH1D JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 70, 0, 70);
    
    TH1D q2vspmq2("q2vspmq2", "Q^2 vs pmq^2", 20, 0, 90);
    // Initialize Pythia for 200 GeV pp collision
    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
    
    pythia.init();
    
    double jet_radius = 0.4;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_radius);
    int dijet_event_counter = 0;

    for (int iEvent = 0; iEvent < 100000; ++iEvent) { // 1M events
        if (!pythia.next()) continue;
    
        
        std::vector<fastjet::PseudoJet> event;
        std::vector<fastjet::PseudoJet> charged_event;

       // pT cut > .2 GeV
        for (int i = 0; i < pythia.event.size(); ++i) {
        Particle& p = pythia.event[i];
    
        // charged particles for eec
        if (!p.isFinal() || p.pT() < 0.2 || !p.isCharged() || abs(p.eta()) > 1.1 ) continue;
        charged_event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
    
        
        // all particles for jet finding 
        if (!p.isFinal() || p.pT() < .2  || abs(p.eta()) > 1.1 ) continue;
        event.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        }
        
        // Perform jet clustering
        fastjet::ClusterSequence cs(event, jet_def);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

        
        //erases (remove if statement) detector acceptance cut jets 
        jets.erase( std::remove_if(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet& jet) { 
                return std::abs(jet.eta()) > 0.7; }),  
            jets.end());
        
        if (jets.size() < 2) continue;
        
        double avjet = (jets[0].pt() + jets[1].pt())/2;
        double pmq2 = pow(avjet, 2); // average of leading jets pt squared
        
        double q2 = pythia.info.Q2Fac();
        
        q2vspmq2.Fill(q2, pmq2);
        
        //! select di-jets on jet pT 
        if(jets[0].pt() < 31.2  || jets[0].pt() >= 40.7) continue;
        if(jets[1].pt() < 20.9 || jets[1].pt() >= 31.2) continue;
        
        // Check for dijet condition
        double dphi = std::abs(jets[0].phi() - jets[1].phi());
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        if (dphi < (3.0 * M_PI / 4.0)) continue;
      
        dijet_event_counter++;
        
        // Jet spectra of dijet events that pass the cut 
         for (size_t i = 0; i < jets.size(); ++i) {
            if (jets[i].pt() < 5) continue;
            JetSpectrum.Fill(jets[i].pt());
        }
        
        //putting thermal particle pushback here 
        vector<PseudoJet> thermalParticles = thermalGen(1500);

        
        vector<PseudoJet> m6 = thermalGen(1500);
        vector<PseudoJet> m7 = thermalGen(1500);
        vector<PseudoJet> m5 = thermalGen(1500);
        vector<PseudoJet> m3 = thermalGen(1500);
        vector<PseudoJet> m4 = thermalGen(1500);
        vector<PseudoJet> m1 = thermalGen(1500);
        vector<PseudoJet> m2 = thermalGen(1500);
 
  
        charged_event.insert(charged_event.end(), thermalParticles.begin(), thermalParticles.end());
            for (size_t i = 0; i < charged_event.size(); ++i) {
                if (charged_event[i].pt() < 0.2) continue; // redundant but safe
                int ui = charged_event[i].user_index();
                for (size_t j = i + 1; j < charged_event.size(); ++j) {
                    if (charged_event[j].pt() < 0.2) continue;
                    double eec = charged_event[i].pt() * charged_event[j].pt();
                    double ctheta = costheta(charged_event[i], charged_event[j]);
                    double z = (1.0 - ctheta) / 2.0;
                    int uj = charged_event[j].user_index();
    
                    EEC_w.Fill(z, eec / pmq2);
                    if (ui == 1 && uj == 1) EEC_t.Fill(z, eec / pmq2);
                    if (ui == -1 && uj == -1) EEC_s.Fill(z, eec / pmq2);
                    if ((ui == -1 && uj == 1) || (ui == 1 && uj == -1)) EEC_m.Fill(z, eec / pmq2);
                }
    
                // signal-background mixes with m1..m4
                for (size_t k = 0; k < m1.size(); ++k) {
                    double sm_eec = charged_event[i].pt() * m1[k].pt();
                    double sm_ctheta = costheta(charged_event[i], m1[k]);
                    double sm_z = (1.0 - sm_ctheta) / 2.0;
                    EEC_ms.Fill(sm_z, sm_eec / pmq2);
                }
                for (size_t k = 0; k < m2.size(); ++k) {
                    double sm2_eec = charged_event[i].pt() * m2[k].pt();
                    double sm2_ctheta = costheta(charged_event[i], m2[k]);
                    double sm2_z = (1.0 - sm2_ctheta) / 2.0;
                    EEC_pt2.Fill(sm2_z, sm2_eec / pmq2);
                }
                for (size_t k = 0; k < m3.size(); ++k) {
                    double sm3_eec = charged_event[i].pt() * m3[k].pt();
                    double sm3_ctheta = costheta(charged_event[i], m3[k]);
                    double sm3_z = (1.0 - sm3_ctheta) / 2.0;
                    EEC_pt3.Fill(sm3_z, sm3_eec / pmq2);
                }
                for (size_t k = 0; k < m4.size(); ++k) {
                    double sm4_eec = charged_event[i].pt() * m4[k].pt();
                    double sm4_ctheta = costheta(charged_event[i], m4[k]);
                    double sm4_z = (1.0 - sm4_ctheta) / 2.0;
                    EEC_pt4.Fill(sm4_z, sm4_eec / pmq2);
                }
            } // end charged_event pair loops
    
        // Now mixed-mixed and other pairings (avoid double counting using k>i)
        // m1-m1
        for (size_t i = 0; i < m1.size(); ++i) {
            for (size_t k = i + 1; k < m1.size(); ++k) {
                double mm11_eec = m1[i].pt() * m1[k].pt();
                double mm_ctheta = costheta(m1[i], m1[k]);
                double mm_z = (1.0 - mm_ctheta) / 2.0;
                EEC_mm.Fill(mm_z, mm11_eec / pmq2);
                EEC_t1t1.Fill(mm_z, mm11_eec / pmq2);
            }
            // m1-m2
            for (size_t l = 0; l < m2.size(); ++l) {
                double mm12_eec = m1[i].pt() * m2[l].pt();
                double mm12_ctheta = costheta(m1[i], m2[l]);
                double mm12_z = (1.0 - mm12_ctheta) / 2.0;
                EEC_t1t2.Fill(mm12_z, mm12_eec / pmq2);
            }
            // m1-m3
            for (size_t m = 0; m < m3.size(); ++m) {
                double mm13_eec = m1[i].pt() * m3[m].pt();
                double mm13_ctheta = costheta(m1[i], m3[m]);
                double mm13_z = (1.0 - mm13_ctheta) / 2.0;
                EEC_t1t3.Fill(mm13_z, mm13_eec / pmq2);
            }
            // m1-m4
            for (size_t n = 0; n < m4.size(); ++n) {
                double mm14_eec = m1[i].pt() * m4[n].pt();
                double mm14_ctheta = costheta(m1[i], m4[n]);
                double mm14_z = (1.0 - mm14_ctheta) / 2.0;
                EEC_t1t4.Fill(mm14_z, mm14_eec / pmq2);
            }
        }

        // m2 group
        for (size_t i = 0; i < m2.size(); ++i) {
            for (size_t k = i + 1; k < m2.size(); ++k) {
                double mm22_eec = m2[i].pt() * m2[k].pt();
                double mm22_ctheta = costheta(m2[i], m2[k]);
                double mm22_z = (1.0 - mm22_ctheta) / 2.0;
                EEC_t2t2.Fill(mm22_z, mm22_eec / pmq2);
            }
            // m2-m3
            for (size_t l = 0; l < m3.size(); ++l) {
                double mm23_eec = m2[i].pt() * m3[l].pt();
                double mm23_ctheta = costheta(m2[i], m3[l]);
                double mm23_z = (1.0 - mm23_ctheta) / 2.0;
                EEC_t2t3.Fill(mm23_z, mm23_eec / pmq2);
            }
            // m2-m4
            for (size_t m = 0; m < m4.size(); ++m) {
                double mm24_eec = m2[i].pt() * m4[m].pt();
                double mm24_ctheta = costheta(m2[i], m4[m]);
                double mm24_z = (1.0 - mm24_ctheta) / 2.0;
                EEC_t2t4.Fill(mm24_z, mm24_eec / pmq2);
            }
        }

        // m3-m3 and m3-m4
        for (size_t i = 0; i < m3.size(); ++i) {
            for (size_t k = i + 1; k < m3.size(); ++k) {
                double mm33_eec = m3[i].pt() * m3[k].pt();
                double mm33_ctheta = costheta(m3[i], m3[k]);
                double mm33_z = (1.0 - mm33_ctheta) / 2.0;
                EEC_t3t3.Fill(mm33_z, mm33_eec / pmq2);
            }
            for (size_t k = 0; k < m4.size(); ++k) {
                double mm34_eec = m3[i].pt() * m4[k].pt();
                double mm34_ctheta = costheta(m3[i], m4[k]);
                double mm34_z = (1.0 - mm34_ctheta) / 2.0;
                EEC_t3t4.Fill(mm34_z, mm34_eec / pmq2);
            }
        }

        // m4-m4 (avoid double count)
        for (size_t i = 0; i < m4.size(); ++i) {
            for (size_t k = i + 1; k < m4.size(); ++k) {
                double mm44_eec = m4[i].pt() * m4[k].pt();
                double mm44_ctheta = costheta(m4[i], m4[k]);
                double mm44_z = (1.0 - mm44_ctheta) / 2.0;
                EEC_t4t4.Fill(mm44_z, mm44_eec / pmq2);
            }
        }

        
    }//! event loop 
    
  
    out->Write();
  
    // Pythia cleanup
    pythia.stat();
    out->Close();
    return 0;
    }
