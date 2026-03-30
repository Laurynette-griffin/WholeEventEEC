//! Authored: Laurynette Griffin Feb 3, 2025
//! Creates a TTree with just particle info for dijet events between 
//! specified pt ranges and back-to-back dijet requirement
//! As well as saving jet information according to different recombination schemes

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace fastjet;
using namespace Pythia8;
using namespace std;

bool dijet_selection(vector<PseudoJet> jets) { 
    
    jets.erase( std::remove_if(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet& jet) { 
                return std::abs(jet.eta()) > 0.7; }),  
            jets.end());
    
    //! For some reason if I ever want to actually look into the jets but I do not need all the soft ones 
    jets.erase( std::remove_if(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet& jet) { 
                return jet.pt() < 5.0; }),  
            jets.end());
    
    if (jets.size() < 2) return false;
    
    if(jets[0].pt() < 31.2  || jets[0].pt() >= 40.7) return false;
    if(jets[1].pt() < 20.9 || jets[1].pt() >= 31.2) return false;
     
    double back2backcut = (3.0 * M_PI / 4.0);
    double dphi = std::abs(jets[0].phi() - jets[1].phi());
    
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    if (dphi < back2backcut) return false;
    
    return true;
}

bool passes_selection(vector<PseudoJet> jets) { 
    
    jets.erase( std::remove_if(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet& jet) { 
                return std::abs(jet.eta()) > 0.7; }),  
            jets.end());
    
    jets.erase( std::remove_if(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet& jet) { 
                return jet.pt() < 20; }),  
            jets.end());
            
    
    jets.erase( std::remove_if(jets.begin(), jets.end(),
        [](const fastjet::PseudoJet& jet) { 
            return jet.pt() > 30; }),  
        jets.end());
    
    if (jets.size() < 1) return false; 
    
    return true;
}

double deltaphi(const PseudoJet& p1, const PseudoJet& p2) {
    float dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    return dphi;
}

double deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  double dphi = fabs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  double deta = p1.eta() - p2.eta();
  return sqrt(deta * deta + dphi * dphi);
}


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
    
    
    if (!out || out->IsZombie()) {
    cerr << "Error: Could not open output file: " << outputFileName << endl;
    return 1;
    }
    
    out->cd();
    cout << "Created output file: " << outputFileName << endl;

    //! Weec tree
    TTree *tree = new TTree("events", "Pythia dijet event data  ");
    tree->SetDirectory(out);
    //! Particle info
    int event;
    vector<double> px, py, pz, e; //! 4-momentum of the particles 
    vector<int> id; 
    
    tree->Branch("event", &event);
    tree->Branch("px", &px); // "branch name" & variable associated with the branch 
    tree->Branch("py", &py);
    tree->Branch("pz", &pz);
    tree->Branch("e", &e);
    tree->Branch("id", &id);
    
    vector<double> *jet_pt  = new vector<double>();
    vector<double> *jet_eta = new vector<double>();
    vector<double> *jet_phi = new vector<double>();

    vector<vector<double>> *c_px = new vector<vector<double>>();
    vector<vector<double>> *c_py = new vector<vector<double>>();
    vector<vector<double>> *c_pz = new vector<vector<double>>();
    vector<vector<double>> *c_e  = new vector<vector<double>>();
    vector<vector<int>>    *c_pid = new vector<vector<int>>();

    //! Jet tree setup 
    auto setupJetTree = [&](TTree* t) {
        t->Branch("event", &event);
        //t->Branch("jet_index", &jet_index);
        t->Branch("jet_pt", &jet_pt);
        t->Branch("jet_eta", &jet_eta);
        t->Branch("jet_phi", &jet_phi);
        t->Branch("c_px", &c_px);
        t->Branch("c_py", &c_py);
        t->Branch("c_pz", &c_pz);
        t->Branch("c_e", &c_e);
        t->Branch("c_pid", &c_pid);
    };
    
    auto pushbackBranches = [&](const vector<PseudoJet>& jets_to_save) {
        
        jet_pt->clear(); jet_eta->clear(); jet_phi->clear();
        c_px->clear();   c_py->clear();    c_pz->clear(); 
        c_e->clear();    c_pid->clear();

        //! jets loop
        for (const auto &j : jets_to_save) {
            jet_pt->push_back(j.pt());
            jet_eta->push_back(j.eta());
            jet_phi->push_back(j.phi());

            //! temp vectors for jet constituents
            vector<double> t_px, t_py, t_pz, t_e;
            vector<int> t_pid;

            for (const auto& c : j.constituents()) {
                t_px.push_back(c.px());
                t_py.push_back(c.py());
                t_pz.push_back(c.pz());
                t_e.push_back(c.e());
                t_pid.push_back(c.user_index());
            }
            //! Push constituent vectors into  vector of vectors
            c_px->push_back(t_px);
            c_py->push_back(t_py);
            c_pz->push_back(t_pz);
            c_e->push_back(t_e);
            c_pid->push_back(t_pid);
        }
    };
    
    TTree *dijet04Tree_s = new TTree("dijets_R04_antikt_s", "E scheme anti-kt R=0.4 jets pt cut 31.2 <= jets[0].pt() < 40.7 & 20.9 <= jets[0].pt() < 31.2");
    setupJetTree(dijet04Tree_s);
    dijet04Tree_s->SetDirectory(out);
    
    TTree *jet04Tree_s = new TTree("jets_R04_antikt_s", "E scheme anti-kt R=0.4 jets");
    setupJetTree(jet04Tree_s);
    jet04Tree_s->SetDirectory(out);
    
    TTree *jet04WTA_s = new TTree("jets_R04_wta_s", "WTA-scheme jets");
    setupJetTree(jet04WTA_s);
    jet04WTA_s->SetDirectory(out);

    TTree *jet04HardCone_s = new TTree("jets_R04_HardCone_s","R=0.4 jets using constituent pt > 2.0 axis + all particles > 0.5");
    setupJetTree(jet04HardCone_s);
    jet04HardCone_s->SetDirectory(out);
    
    //! Jet clustering definitions 
    double r = 0.4;
    JetDefinition jetDef_hc(antikt_algorithm, r, E_scheme); //! Redundant but it makes me feel better :)
    JetDefinition jetDef_es(antikt_algorithm, r, E_scheme);
    JetDefinition jetDef_wta(antikt_algorithm, r, WTA_pt_scheme);
    Recluster reclust(jetDef_wta);

    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
    
    pythia.init();
    
    //! Event loop 
    int jetevent_cnt = 0;
    int dijetevent_cnt = 0;
    int es_cnt = 0;
    int wta_cnt = 0;
    int hc_cnt = 0;
    long int nEvents = 100000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        event = iEvent;
        if (!pythia.next()) continue;
        
        px.clear(); py.clear(); pz.clear(); e.clear();
        id.clear(); 

        std::vector<PseudoJet> particles;
        std::vector<PseudoJet> charged_particles;
        std::vector<PseudoJet> hc_particles;
        for (int i = 0; i < pythia.event.size(); ++i) {
            Particle& p = pythia.event[i];
            
            if (!p.isFinal() || fabs(p.eta()) > 1.1 || p.pT() < 0.5) continue;
            double px_i = pythia.event[i].px();
            double py_i = pythia.event[i].py();
            double pz_i = pythia.event[i].pz();
            double e_i  = pythia.event[i].e();
            double id_i = pythia.event[i].id();
            
            px.push_back(px_i);
            py.push_back(py_i);
            pz.push_back(pz_i);
            e.push_back(e_i);
            id.push_back(id_i);
            
            PseudoJet pj(px_i, py_i, pz_i, e_i);
            pj.set_user_index(id_i);
            particles.push_back(pj);
            if(pythia.event[i].isCharged()) {
                charged_particles.push_back(pj); //! ch
                
                //! def should be a better way to redefine 
            }
            if(pj.pt() > 2.0) hc_particles.push_back(pj);
 
        }//! particle loop close

        vector<PseudoJet> axis_particles;
        
        //! Cluster jets
        ClusterSequence clustSeq_es(particles, jetDef_es);
        ClusterSequence clustSeq(particles, jetDef_es);
        ClusterSequence clustSeq_hc_pt2(hc_particles, jetDef_hc);
        ClusterSequence clustSeq_wta(particles, jetDef_wta);
        //! will be reclusted into wta jets? 
        
        
        vector<PseudoJet> jets = sorted_by_pt(clustSeq.inclusive_jets()); 
        //! are the e scheme jets
        vector<PseudoJet> es_jets = sorted_by_pt(clustSeq_es.inclusive_jets()); 
        //! base for hc jet finding 
        vector<PseudoJet> hc_jets_pt2 = sorted_by_pt(clustSeq_hc_pt2.inclusive_jets()); 
        //! idk but maybe wta jets? idk 
        //vector<PseudoJet> wta_jets = sorted_by_pt(clustSeq_wta.inclusive_jets()); 
        vector<PseudoJet> wta_jets;
        
        
        //! Dijet selection for weec
        if (dijet_selection(jets)){
            //! fills tree & dijet tree if cuts are met
            tree->Fill();
            
            pushbackBranches(jets);
            dijet04Tree_s->Fill();
            dijetevent_cnt ++;
        }
        
        if (!passes_selection(es_jets)) continue;
        
        for(auto& jet : es_jets){
            PseudoJet wta_jet = reclust(jet);
            wta_jets.push_back(wta_jet);
        }
        
        
        if (!passes_selection(wta_jets)) continue;
        
        //! Need to use the hc jet axis and do a recombination using escheme for the purpose of applying cuts 
        vector<PseudoJet> hc_jets;
        
        //! vector to keep our ClusterSequences alive
        vector<shared_ptr<ClusterSequence>> kept_cs; 

        for (const auto& hcjet : hc_jets_pt2) {
            //! Makes sure the particles in the new jet are only from this hc jet 
            axis_particles.clear();
            
            double eta0 = hcjet.eta();
            double phi0 = hcjet.phi();
            for (const auto& p : particles) {
                double dr = deltaR(hcjet, p);
                if (dr < r) {
                    axis_particles.push_back(p);
                }
            }
        
            if (axis_particles.empty()) continue;
            
            //! Dynamically allocate the ClusterSequence using shared_ptr
            auto cs = make_shared<ClusterSequence>(axis_particles, jetDef_es);
            
            auto jets_hc_rc = sorted_by_pt(cs->inclusive_jets());
            
            if (passes_selection(jets_hc_rc)) {
                hc_jets.push_back(jets_hc_rc[0]);
                //! NEW: Push the shared_ptr to our keeping vector so it doesn't get destroyed
                kept_cs.push_back(cs); 
            }
        }
        
        if (hc_jets.empty()) continue;
        //! After checking all 3 jets exists fill trees
        es_cnt++ ;
        wta_cnt++ ;
        hc_cnt++ ;
        jetevent_cnt++ ;
        
        pushbackBranches(hc_jets);
        jet04HardCone_s->Fill();
        
        pushbackBranches(es_jets);
        jet04Tree_s->Fill();
        
        pushbackBranches(wta_jets);
        jet04WTA_s->Fill();

        jetevent_cnt+= 1;
        
    }//! event loop close
    cout<< "dijet Event count = " << dijetevent_cnt << endl;
    cout<< "Event count = " << jetevent_cnt << endl;
    cout<< "e scheme jet count = " << es_cnt << endl;
    cout<< "wta jet count = " << wta_cnt << endl;
    cout<< "hard cone count = " << hc_cnt << endl;
    pythia.stat();
    //! Writes the information to our root file :) 
    out->Write();
    out->Close();
}
