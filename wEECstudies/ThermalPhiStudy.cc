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

using namespace std;
using namespace fastjet;

// deltaphi 
double deltaphi(const PseudoJet& p1, const PseudoJet& p2) {
 double dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi; // ensuring delta phi is between 0 and pi 
    return dphi;
}

// 1500 particles in most central 200 GeV AuAu collisions 
vector <PseudoJet> thermalbackground(int nParticles = 1500){ 

    double mass = .13957; // charged pion mass 
    
    TRandom3 *r = new TRandom3(0); // random number generator with time dependent seed which is truly random
    
    vector <PseudoJet> thermalparticles;
    for (int i=0; i<nParticles; i++){
        // im lazy so im just gonna make a uniform pt distribution between .3 and 10 
        double pt = r->Uniform(.3, 10);
        
        // r->Rndm() is what actually gives you a number and its between 0 and 1 
        double eta =2.2*r->Rndm() - 1.1;
        double phi = 2.0*TMath::Pi()*r->Rndm() - TMath::Pi();
        
        double p = pt*cosh(eta);
        double px = pt*cos(phi);
        double py = pt*sin(phi);
        double pz = pt*sinh(eta);
        
        double E = sqrt(mass*mass + p*p);
        
        PseudoJet part(px, py, pz, E);
        thermalparticles.push_back(part);
        }
    return thermalparticles;
    }    

int main(int argc, char* argv[]){
     if (argc < 2){
          cerr << "No output file :(" << endl;
        return 1;
    }
    
    string outputFileName = argv[1];
    TFile *out = new TFile(outputFileName.c_str(), "RECREATE");
    
    out->cd();
    cout << "Created output file: " << outputFileName << endl;
   
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    
    // initialize histogram stuff double log bins from 0 to pi
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

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    int bins = 60;
    int morebins = 3600;
    
    // hist initialization
    TH1F EEC_linear1("EEC_l1", "Linear Delta Phi EEC 1h", bins, 0, TMath::Pi());
    TH1F EEC_linear2("EEC_l2", "Linear Delta Phi EEC 1k", bins, 0, TMath::Pi());
    TH1F EEC_linear3("EEC_l3", "Linear Delta Phi EEC 10k", bins, 0, TMath::Pi());
    
    TH1F EEC_biglinear1("EEC_bl1", "Linear Delta Phi EEC 1h Granular", morebins, 0,TMath::Pi());
    TH1F EEC_biglinear2("EEC_bl2", "Linear Delta Phi EEC 1k Granular", morebins, 0,TMath::Pi());
    TH1F EEC_biglinear3("EEC_bl3", "Linear Delta Phi EEC 10k Granular", morebins, 0,TMath::Pi());
    
    TH1F EEC_doublelog1("EEC_double1", "Double log Delta Phi EEC 1h", bins, topi);
    TH1F EEC_doublelog2("EEC_double2", "Double log Delta Phi EEC 1k", bins, topi);
    TH1F EEC_doublelog3("EEC_double3", "Double log Delta Phi EEC 10k", bins, topi);
    
    TH1F particlept("particlept", "p_{T} of particles", 10, 0, 10);
    TH2F etaphi("etaphi", "eta phi map", 22, -1.1, 1.1, 63, 0, 2*TMath::Pi());
    
    
    // funxtion :)
    for (int iEvent = 0; iEvent < 10000; iEvent ++){
        
        vector<PseudoJet> thermalparticles = thermalbackground(1500);
        
        for(size_t i = 0; i < thermalparticles.size(); i++){
            etaphi.Fill(thermalparticles.at(i).eta(),thermalparticles.at(i).phi());
            particlept.Fill(thermalparticles.at(i).pt());
            
            for(size_t j = i + 1; j < thermalparticles.size(); j++){
                double eec = (thermalparticles.at(i).pt()*thermalparticles.at(j).pt());
                double delphi = deltaphi(thermalparticles.at(i),thermalparticles.at(j));
                
                if (iEvent < 100){
                    EEC_linear1.Fill(delphi, eec);
                    EEC_biglinear1.Fill(delphi, eec);
                    EEC_doublelog1.Fill(delphi, eec);
                }
                
                 if (iEvent < 1000){
                    EEC_linear2.Fill(delphi, eec);
                    EEC_biglinear2.Fill(delphi, eec);
                    EEC_doublelog2.Fill(delphi, eec);
                }
                
                EEC_linear3.Fill(delphi, eec);
                EEC_biglinear3.Fill(delphi, eec);
                EEC_doublelog3.Fill(delphi, eec);
            }
        }
    }  
out->Write();
out->Close();
return 0;
}
    
    