#include <iostream>
#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
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
    
vector<PseudoJet> thermalGen(int nParticles = 570){
    
    TF1 *mb = new TF1("mb","exp(-[0]*x/[1])",0.3,10); // pt range from .3 to 10 Gev
    mb->SetParameters(1,0.26); // shape and temperature (300 MeV)
    
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
    TH1F EEC_t("EEC_t", "Energy Energy Correlator", bins, eecbounds);
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
    //TH1F Deltaphi("Deltaphi", "linear del phi", 31, 0, TMath::Pi());
    //TH2F etavphit("hEtaPhi", "Eta vs Phi;#eta;#phi",22, -1.1, 1.1, 32, -2*TMath::Pi(), 2*TMath::Pi()); 
    
    // reco vs truth q2
    TH2F qvspmq("comparisonpthat", "pthat^2 vs Average leading and subleading jet pt^2", 12, 0, 60, 10, 0, 60);
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

    for (int iEvent = 0; iEvent < 1000000; ++iEvent) { // 1M events
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
        if(jets[0].pt() < 31.2  || jets[0].pt() >= 40.7) continue;
        if(jets[1].pt() < 27.3 || jets[1].pt() >= 31.2) continue;
        
        
        // Check for dijet condition
        double dphi = std::abs(jets[0].phi() - jets[1].phi());
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        if (dphi < (3.0 * M_PI / 4.0)) continue;
       
        // filling q 2d hist
        double avjet = (jets[0].pt() + jets[1].pt())/2;
        double pmq = avjet;
        double pthat = pythia.info.pTHat();
        qvspmq.Fill(avjet, pthat);
        
        if (iEvent < 1000) cout << "pmq = " << pmq << " pTHat = " << pthat <<  endl;
        
        //! select di-jets on jet pT 
        //if(jets[0].pt() < 31.2  || jets[0].pt() >= 40.7) continue;
        //if(jets[1].pt() < 27.3 || jets[1].pt() >= 31.2) continue;
        
        dijet_event_counter++;
        
        
        // Jet spectra of dijet events that pass the cut 
         for (size_t i = 0; i < jets.size(); ++i) {
            if (jets[i].pt() < 5) continue;
            JetSpectrum.Fill(jets[i].pt());
        }
        
        //jet spectra
        LeadingJetSpectrum.Fill(jets[0].pt());
        SubleadingJetSpectrum.Fill(jets[1].pt());
        
        //putting thermal particle pushback here 
        vector<PseudoJet> thermalParticles = thermalGen(570);
        
        /*
        for (size_t i = 0; i < thermalParticles.size(); ++i) {
            //etavphit.Fill(thermalParticles.at(i).eta(), thermalParticles.at(i).phi());
            ThermalPhi.Fill(thermalParticles.at(i).phi());
            for (size_t j = i + 1; j < thermalParticles.size(); ++j) {
                double eect = thermalParticles.at(i).pt() * thermalParticles.at(j).pt();  
                double cthetat = costheta(thermalParticles.at(i), thermalParticles.at(j));
                double delphit = deltaphi(thermalParticles.at(i), thermalParticles.at(j));
                double zt = (1 - cthetat)/2;
                //EEC_t.Fill(zt, eect);
                //EEC_td.Fill(delphit, eect);
                //Deltaphi.Fill(delphit, eect);
                
            }
        }
        */
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
    
    }//! event loop 
    
    cout << " total # of dijet events = " << dijet_event_counter << endl;
  
    out->Write();
  
    // Pythia cleanup
    pythia.stat();
    out->Close();
    return 0;
    }

