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

TRandom3 *r = new TRandom3(0); // be sure to set seed to 0 which is a time dependant seed once before thermal particle gen 

vector<PseudoJet> thermalGen(int nParticles = 1500){
    
    TF1 *mb = new TF1("mb","exp(-[0]*x/[1])",0.3,10); // pt range from .3 to 10 Gev
    mb->SetParameters(1,0.26); // shape and temperature (260 MeV)
    
    double mass = 0.13957; //charged pions 
    
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

vector<double> createDoubleLogBins(int nBins, double minVal, double maxVal) {
    vector<double> bins;
    
    if (minVal <= 0) minVal = 1e-6; 

    int nHalf = nBins / 2;

    double midPoint = maxVal / 2;
    double logMin = log10(minVal);
    double logMax = log10(midPoint);

    //! Generates lower half
    for (int i = 0; i <= nHalf; i++) {
        double val = pow(10, logMin + i * (logMax - logMin) / nHalf);
        bins.push_back(val);
    }

    //! generates upper half
    for (int i = bins.size() - 2; i >= 0; i--) {
        double mirroredVal = 1.0 - bins[i];
        bins.push_back(mirroredVal);
    }
    return bins;
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
    int bins = 40;
    //vector<double> topi;
    //topi = createDoubleLogBins(bins, 1e-3, M_PI);
    //cout << topi.data() << endl;

    
    
     Double_t topi[41] =   {1.000000000000e-03,1.444794214505e-03,2.087430322267e-03,3.015907252794e-03,4.357365350321e-03,
                            6.295496248628e-03,9.095696557457e-03,1.314140976311e-02,1.898663279618e-02,2.743177721685e-02,
                            3.963327301649e-02,5.726192355613e-02,8.273169586533e-02,1.195302755424e-01,1.726966505619e-01,
                            2.495111215962e-01,3.604922249369e-01,5.208370809628e-01,7.525024012748e-01,1.087211115763e+00,
                            1.570796330000e+00,
                            2.054381537827e+00,2.389090252315e+00,2.620755572627e+00,2.781100428653e+00,2.892081531994e+00,
                            2.968896003028e+00,3.022062378047e+00,3.058860957724e+00,3.084330730034e+00,3.101959380573e+00,
                            3.114160876373e+00,3.122606020794e+00,3.128451243827e+00,3.132496957032e+00,3.135297157341e+00,
                            3.137235288239e+00,3.138576746337e+00,3.139505223268e+00,3.140147859375e+00,M_PI};
    
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
    
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    // whole event with thermal and mixed seperated out   
    TH1D EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_t("EEC_t", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_m("EEC_m", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_s("EEC_s", "Energy Energy Correlator", bins, eecbounds);
    
    TH1D EEC_ms("EEC_ms", "Energy Energy Correlator", bins, eecbounds); // mix signal
    TH1D EEC_mm("EEC_mm", "Energy Energy Correlator", bins, eecbounds); // mix1 mix1 
    TH1D EEC_m2("EEC_m2", "Energy Energy Correlator", bins, eecbounds); // mix1 mix2
 
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


    TH1D JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 70, 0, 70);
    
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
        if (!p.isFinal() || p.pT() < .2 || !p.isCharged() || abs(p.eta()) > 1.1 ) continue;
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
        
        // Require atleast two jets
        if (jets.size() < 2) continue;
        
        //! select di-jets on jet pT 
        if(jets[0].pt() < 31.2  || jets[0].pt() >= 40.7) continue;
        if(jets[1].pt() < 20.9 || jets[1].pt() >= 31.2) continue;
        
        
        // Check for dijet condition
        double dphi = std::abs(jets[0].phi() - jets[1].phi());
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        if (dphi < (3.0 * M_PI / 4.0)) continue;
      
        dijet_event_counter++;
        double avjet = (jets[0].pt() + jets[1].pt())/2;
        double pmq2 = pow(avjet, 2); // average of leading jets pt squared
        
        // Jet spectra of dijet events that pass the cut 
         for (size_t i = 0; i < jets.size(); ++i) {
            if (jets[i].pt() < 5) continue;
            JetSpectrum.Fill(jets[i].pt());
        }
        
        double avthermalpt = 0;
        //putting thermal particle pushback here 
        vector<PseudoJet> thermalParticles = thermalGen(1500);
        vector<PseudoJet> m3 = thermalGen(1500);
        vector<PseudoJet> m4 = thermalGen(1500);
        
        vector<PseudoJet> m1 = thermalGen(1500);
        vector<PseudoJet> m2 = thermalGen(1500);
        charged_event.insert(charged_event.end(), thermalParticles.begin(), thermalParticles.end());
        for (size_t i = 0; i < charged_event.size(); ++i) {
            // mixed event loop? I just need 1 signal particle so this should work and it doesnt increase the loops
            for (size_t j = i + 1; j < charged_event.size(); ++j) {

                double eec = charged_event.at(i).pt() * charged_event.at(j).pt();  
                double ctheta = costheta(charged_event.at(i), charged_event.at(j));
                double z = (1 - ctheta)/2; 
                double delphi = deltaphi(charged_event.at(i), charged_event.at(j));
                int ui = charged_event.at(i).user_index();
                int uj = charged_event.at(j).user_index();
                
                EEC_w.Fill(z, eec/pmq2);
                EEC_wp.Fill(delphi, eec/pmq2);
                EEC_wpl.Fill(delphi, eec/pmq2);
                
                if (ui == 1 && uj == 1) { //thermal
                    EEC_t.Fill(z, eec/pmq2);
                    EEC_tp.Fill(delphi, eec/pmq2);
                    EEC_tpl.Fill(delphi, eec/pmq2);
                }
                
                if (ui == -1 &&  uj == -1) { //signal
                    EEC_s.Fill(z, eec/pmq2);
                    EEC_sp.Fill(delphi, eec/pmq2);
                    EEC_spl.Fill(delphi, eec/pmq2);                
                    
                }
                    
                if ((ui == -1 &&  uj == 1) || (ui == 1 &&  uj == -1 ))  { //mixed
                    EEC_m.Fill(z, eec/pmq2);
                    EEC_mp.Fill(delphi, eec/pmq2);
                    EEC_mpl.Fill(delphi, eec/pmq2);
                }
             
            }
      
            for (size_t k = 0; k < m1.size(); ++k) { // every k with the one i and this loops over i 
        
            // Signal Background 
            double sm_eec = charged_event.at(i).pt() * m1.at(k).pt();
            
            double sm_ctheta = costheta(charged_event.at(i), m1.at(k));
            double sm_z = (1 - sm_ctheta)/2; 
            double sm_delphi = deltaphi(charged_event.at(i),m1.at(k));
            
            EEC_ms.Fill(sm_z, sm_eec/pmq2);
            EEC_msp.Fill(sm_delphi,sm_eec/pmq2);
            EEC_mspl.Fill(sm_delphi,sm_eec/pmq2);
    
            }
        }//!Whole event EEC loop close 
        
        // mix1 mix1 
        for (size_t i = 0; i < m1.size(); ++i) {
                for (size_t k = i + 1; k < m1.size(); ++k) {
                    double mm11_eec = m1[i].pt() * m1[k].pt();
                    double mm_ctheta = costheta(m1[i], m1[k]);
                    double mm_z = (1 - mm_ctheta) / 2;
                    double mm_delphi = deltaphi(m1[i], m1[k]);
            
                    EEC_mm.Fill(mm_z, mm11_eec / pmq2);
                    EEC_mmp.Fill(mm_delphi, mm11_eec / pmq2);
                    EEC_mmpl.Fill(mm_delphi, mm11_eec / pmq2);
                }
            }
        
        //mix1 mix2
        for (size_t i = 0; i < m1.size(); ++i) {
            for (size_t k = 0; k < m2.size(); ++k) {
                double mm12_eec = m1[i].pt() * m2[k].pt();
                double mm12_ctheta = costheta(m1[i], m2[k]);
                double mm12_z = (1 - mm12_ctheta) / 2;
                double mm12_delphi = deltaphi(m1[i], m2[k]);
        
                EEC_m2.Fill(mm12_z, mm12_eec / pmq2);
                EEC_m2p.Fill(mm12_delphi, mm12_eec / pmq2);
                EEC_m2pl.Fill(mm12_delphi, mm12_eec / pmq2);
            }
        }
        
    }//! event loop 
    
  
    out->Write();
  
    // Pythia cleanup
    pythia.stat();
    out->Close();
    return 0;
    }
