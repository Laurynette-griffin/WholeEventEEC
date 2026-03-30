#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include "fastjet/ClusterSequence.hh"
#include "TLorentzVector.h"
#include "TComplex.h"
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
#include "TTree.h"

using namespace fastjet;
using namespace Pythia8;
using namespace std;

//! Delta R calculation (uses eta)
double deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  double deta = p1.eta() - p2.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}

//! costheta calculation
double costheta(const PseudoJet& p1, const PseudoJet& p2) {
    double dotprod = (p1.px() * p2.px() +  p1.py() * p2.py() +  p1.pz() * p2.pz());
    double normp1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz()); 
    double normp2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz()); 
    
    if (normp1 == 0 || normp2 == 0) // division by 0 check 
        return 0; 
    return (dotprod / (normp1 * normp2));
}

//! delphi calculation 
double deltaphi(const PseudoJet& p1, const PseudoJet& p2) {
    double dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    return dphi;
}

int getCentralityBin(int mult, const vector<double>& edges) {
    for (int i = 0; i < 5; ++i) {
        //! edges[i] is the upper bound, edges[i+1] is the lower bound
        if (mult <= edges[i] && mult > edges[i+1]) return i;
    }
    return -1;
}

struct ThermalEvent {
    vector<PseudoJet> particles;
    double header_EPangle;    //! True event plane angle from file
    TComplex Qn[3];           //! Q-vectors for n=1, 2, 3
    double psi_n[3];          //! Calculated event plane angles
    double vn[3];  //! Event-averaged v1, v2, v3
    int mult;
    
};

//! Hydro reader: parses file and pre-calculates Q-vectors/v2
vector<ThermalEvent> readSingleThermalFile(const string& filename) {
    vector<ThermalEvent> bkgevents;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open thermal file: " << filename << endl;
        return bkgevents;
    }

    string line;
    while (getline(file, line)) {
        if (line.rfind("#\tEvent\t", 0) == 0) {
            int event_number, n_hadrons;
            double EPangle;
            if (sscanf(line.c_str(), "#\tEvent\t%d\tweight\t%*f\tEPangle\t%lf\tN_hadrons\t%d",
                       &event_number, &EPangle, &n_hadrons) != 3) continue;

            ThermalEvent ev;
            ev.header_EPangle = EPangle;
            ev.particles.reserve(n_hadrons);

            for (int i = 0; i < n_hadrons; ++i) {
                if (!getline(file, line)) break;
                int index, pid, status;
                double E, px, py, pz;
                stringstream ss(line);
                if (!(ss >> index >> pid >> status >> E >> px >> py >> pz)) continue;

                //! Charged particle selection
                if (!(abs(pid) == 211 || abs(pid) == 321 || abs(pid) == 2212)) continue;

                ev.particles.emplace_back(px, py, pz, E);
                auto& pj = ev.particles.back();
                pj.set_user_index(1); //! Mark as Thermal Background

                //! Kinematic cuts
                if(pj.pt() < 0.5 || fabs(pj.eta()) > 1.1) {
                    ev.particles.pop_back();
                }
            }

            if (!ev.particles.empty()) {
                
                ev.mult = ev.particles.size();
                
                //! Calculate Q-vectors and vn for the loaded event
                for(int iH=0; iH<3; iH++) {
                    int n = iH + 1;
                    ev.Qn[iH] = TComplex(0, 0);
                    for(const auto& p : ev.particles) {
                        ev.Qn[iH] += TComplex(cos(n * p.phi()), sin(n * p.phi()));
                    }
                    ev.psi_n[iH] = TMath::ATan2(ev.Qn[iH].Im(), ev.Qn[iH].Re()) / n;
                    
                    double sum_vn = 0;
                    for(const auto& p : ev.particles) sum_vn += cos(n * (p.phi() - ev.psi_n[iH]));
                    ev.vn[iH] = sum_vn / ev.particles.size();
                }
                bkgevents.push_back(move(ev));
            }
        }
    }
    return bkgevents;
}
/*

//!Hydro reader
vector<vector<fastjet::PseudoJet>>
readSingleThermalFile(const string& filename)
{
    vector<vector<PseudoJet>> bkgevents;

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open thermal file: " << filename << endl;
        return bkgevents;
    }

    cout << "Reading thermal file: " << filename << endl;

    string line;
    int event_number = -1;
    int n_hadrons = 0;
    double EPangle = 0.0;

    while (getline(file, line)) {

        // Only trigger on REAL event headers
        if (line.rfind("#\tEvent\t", 0) == 0) {

            if (sscanf(line.c_str(),
                "#\tEvent\t%d\tweight\t%*f\tEPangle\t%lf\tN_hadrons\t%d",
                &event_number, &EPangle, &n_hadrons) != 3)
            {
                cerr << "WARNING: Could not parse event header:\n"
                     << line << endl;
                continue;
            }

            vector<PseudoJet> bkgevent;
            bkgevent.reserve(n_hadrons);

            for (int i = 0; i < n_hadrons; ++i) {
                if (!getline(file, line)) {
                    cerr << "ERROR: Unexpected EOF in " << filename << endl;
                    break;
                }

                int index, pid, status;
                double E, px, py, pz;

                stringstream ss(line);
                if (!(ss >> index >> pid >> status >> E >> px >> py >> pz)) {
                    cerr << "ERROR parsing particle line:\n"
                         << line << endl;
                    continue;
                }
                //! Charged Mult Sim (sPHENIX) does not include e/mu
                bool isChargedFinal = (abs(pid) ==  211  ||
                                        abs(pid) == 321  || 
                                        abs(pid) == 2212 );
                                        //||
                                        //abs(pid) == 11   ||
                                        //abs(pid) == 13);
                                        
                                
                if (!isChargedFinal) continue;
                
                PseudoJet pj(px, py, pz, E);
                pj.set_user_index(1);
                
                if(pj.pt() < 0.5) continue; //! This pt cut is not based on the paper based on the mult i get out. 
                if(fabs(pj.eta()) > 1.1) continue;
                bkgevent.push_back(pj);
            }

            if (!bkgevent.empty()) {
                bkgevents.push_back(move(bkgevent));
            }
        }
    }

    cout << "  -> Loaded " << bkgevents.size()
         << " events from " << filename << endl;

    return bkgevents;
}
*/
//! Thermal struct
struct ThermalLibrary {
    vector<ThermalEvent> events;
    mutable mt19937 gen;
    mutable uniform_int_distribution<> dist;

    ThermalLibrary() : gen(random_device{}()), dist(0, 0) {}

    void initialize(vector<ThermalEvent>&& input) {
        events = move(input);
        if (!events.empty()) {
            dist = uniform_int_distribution<>(0, events.size() - 1);
        }
    }

    const ThermalEvent& randomEvent() const {
        if (events.empty()) {
            throw runtime_error("Thermal library empty");
        }
        return events[dist(gen)];
    }

    size_t size() const { return events.size(); }
};

ThermalLibrary loadThermalEvents(const vector<string>& filenames){
    vector<ThermalEvent> allEvents;

    for (const auto& f : filenames) {
        auto evts = readSingleThermalFile(f);
        allEvents.insert(allEvents.end(),
                         make_move_iterator(evts.begin()),
                         make_move_iterator(evts.end()));
    }

    ThermalLibrary lib;
    lib.initialize(move(allEvents));
    return lib;
}

class PythiaTreeReader {
private:
    TFile* file = nullptr;
    TTree* tree = nullptr;

    vector<float>* px = nullptr;
    vector<float>* py = nullptr;
    vector<float>* pz = nullptr;
    vector<float>* e  = nullptr;
    vector<int>*   id = nullptr;
    vector<int>*   status = nullptr;

    long nEntries = 0;

public:
    PythiaTreeReader(const string& filename)
    {
        file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie())
            throw runtime_error("Cannot open Pythia ROOT file");

        tree = dynamic_cast<TTree*>(file->Get("events"));
        if (!tree)
            throw runtime_error("Cannot find TTree 'events'");

        tree->SetBranchAddress("px", &px);
        tree->SetBranchAddress("py", &py);
        tree->SetBranchAddress("pz", &pz);
        tree->SetBranchAddress("e",  &e);
        tree->SetBranchAddress("id", &id);
        tree->SetBranchAddress("status", &status);

        nEntries = tree->GetEntries();
        cout << "Opened Pythia tree with " << nEntries << " events" << endl;
    }

    ~PythiaTreeReader() {
        if (file) {
            file->Close();
            delete file;
            file = nullptr;
        }
    }


    long entries() const { return nEntries; }

    vector<PseudoJet> getEvent(long i)
    {
        if (i < 0 || i >= nEntries)
            throw out_of_range("Event index out of range");

        tree->GetEntry(i);

        vector<PseudoJet> signal_events;
        signal_events.reserve(px->size());

        for (size_t j = 0; j < px->size(); ++j) {

            PseudoJet pj(px->at(j), py->at(j),
                          pz->at(j), e->at(j));
            pj.set_user_index(id->at(j));
            signal_events.push_back(pj);
        }
        return signal_events;
    }
};



//allows for the .cmnd file to be accepted as an argument 
int main(int argc, char* argv[]) {
    // Check for required arguments
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <config.cmnd> <output.root> <PythiaROOTFile>  <thermalFiles>" << endl;
        return 1;
    }

    string configFile = argv[1];
    string outputFileName = argv[2];
    string pythiaFile = argv[3];
    // ARGUMENTS FOR PARALLELIZATION
    long int iEventStart = stol(argv[4]); // Start processing from this event index
    long int iEventEnd = stol(argv[5]);   // Stop processing at this event index (exclusive)
    
    vector<string> thermalFiles;
    for (int i = 6; i < argc; ++i) {
        thermalFiles.push_back(argv[i]);
    }
    
    TFile *in = TFile::Open(pythiaFile.c_str(), "READ");
    if (!in || in->IsZombie()) {
        cerr << "Error: Could not open Pythia input file: " << pythiaFile << endl;
        return 1;
    }
    
    TTree *tree = (TTree*)in->Get("events");  // our tree was named "events"
    if (!tree) {
        cerr << "Error: could not find TTree named 'events'" << std::endl;
        return 1;
    }
    TTree *dijets_R04_antikt_s = (TTree*)in->Get("dijets_R04_antikt_s");  // our tree was named "events"
    if (!dijets_R04_antikt_s) {
        cerr << "Error: could not find TTree named 'dijets_R04_antikt_s'" << std::endl;
        return 1;
    }
    

    TFile *out = new TFile(outputFileName.c_str(), "RECREATE");
    if (!out || out->IsZombie()) {
        cerr << "Error: Could not open output file: " << outputFileName << endl;
        return 1;
    }
    out->cd();
    cout << "Created output file: " << outputFileName << endl;
    
    vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *e = nullptr;
    vector<int> *id = nullptr;
    
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("id", &id);
    
    vector<double> *jet_pt = nullptr;
    dijets_R04_antikt_s->SetBranchAddress("jet_pt", &jet_pt);
    
    //PythiaTreeReader* pythiaReader = nullptr; //! just will access the tree directly 
    ThermalLibrary thermalLibrary;

    try {
        //pythiaReader = new PythiaTreeReader(pythiaFile);
        thermalLibrary = loadThermalEvents(thermalFiles);
    } 
    catch (const runtime_error& e) {
        cerr << "Initialization Error: " << e.what() << endl;
        return 1;
    }
    
    int nThermalPool = (int)thermalLibrary.size();
    
       if (nThermalPool == 0) {
        cerr << "ERROR: Thermal library is empty. Cannot perform background subtraction." << endl;
        return 1;
    }
    //! Multiplicity finder for these 1k hydro events
    vector<int> allMults;
    for (const auto& ev : thermalLibrary.events) {
        allMults.push_back(ev.mult);
    }
    
    sort(allMults.begin(), allMults.end(), greater<int>());
    int min_mult = allMults[0]; //! 10% centrality 
    int max_mult = allMults[allMults.size() - 1]; //! 0% centrality 
    
    TH1D hMult("hMult", "Charged Multiplicity;N_{ch};Counts", 25, min_mult, max_mult); 
    for(const auto& ev : thermalLibrary.events) hMult.Fill(ev.mult);
    
    hMult.Scale(1.0/hMult.Integral());

    int nBins = 20; 
    vector<double> centEdges;
    
    //! 0% edge is the high-multiplicity boundary (far right of hMult)
    centEdges.push_back(hMult.GetBinLowEdge(hMult.GetNbinsX() + 1)); 
    
    double runningArea = 0.0;
    int currentTargetbin = 1; //! current target centrality bin
    
    //! going from left to right because 0% centrality is far left 
    for (int i = hMult.GetNbinsX(); i >= 1; --i) { 
        double binContent = hMult.GetBinContent(i);
        double areaBeforeBin = runningArea;
        runningArea += binContent;
    
        //! Target area for the current bin boundary
        double targetArea = (double)currentTargetbin / nBins;
    
        //! 3. Check if the current bin has crossed the target threshold
        if (runningArea >= targetArea) {
            //! Calculate distance to target from current bin edge vs previous bin edge
            double distCurrent = fabs(runningArea - targetArea); //! distance of current bin edge to centrality target 
            double distPrev = fabs(areaBeforeBin - targetArea);
    
            double chosenEdge;
            if (distCurrent < distPrev) {
                //! Current bin low edge is closer
                chosenEdge = hMult.GetBinLowEdge(i);
                cout << "Chosen bin edge cooresponds to centrality " <<  areaBeforeBin << endl;
            } else {
                //! Previous bin low edge (which is this bin's high edge) is closer
                chosenEdge = hMult.GetBinLowEdge(i + 1);
                cout << "Chosen bin edge cooresponds to centrality " << runningArea<< endl;
            }
    
            centEdges.push_back(chosenEdge);
            currentTargetbin++;
           
    
            //! Stop if we have found all requested bin edges
            if (currentTargetbin > nBins) break;
        }
    }
    
    
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
    
    TH1D v1("v1", "avg linear flow", 20 , 0.0, 0.8);
    TH1D v2("v2", "avg elliptical flow", 20 , 0.0, 0.8);
    TH1D v3("v3", "avg triangular flow", 20 , 0.0, 0.8);
    //int nEntries = tree->GetEntries();
    
    long int nEvents = tree->GetEntries();
    if (iEventEnd > nEvents) iEventEnd = nEvents;
    if (iEventStart >= nEvents) return 1;
    for (long int iEvent = iEventStart; iEvent < iEventEnd; iEvent++) {
        if(iEvent % 100 == 0) cout << "Processing Event " << iEvent << "/" << nEvents << endl;
        //! Init thermal rng here 
        vector<fastjet::PseudoJet> charged_event;
        tree->GetEntry(iEvent);
        
        for (size_t k = 0; k < px->size(); ++k) {
            int pdg = id->at(k);
            int absPdg = std::abs(pdg);
            
            bool isCharged = (absPdg == 211 ||  // pions
                  absPdg == 321 ||  // kaons
                  absPdg == 2212 || // protons
                  absPdg == 11 ||   // electrons
                  absPdg == 13);    // muons
            
            if (!isCharged) continue;
            PseudoJet p(px->at(k), py->at(k), pz->at(k), e->at(k));
            p.set_user_index(-1); //! Mark as Signal
            charged_event.push_back(p);
        }
        
        dijets_R04_antikt_s->GetEntry(iEvent); 
        
        if (!jet_pt || jet_pt->size() < 2) continue;
        
        double jet0_pt = jet_pt->at(0);
        double jet1_pt = jet_pt->at(1);
        
        double avjet = (jet0_pt + jet1_pt) / 2.0;
        
        double pmq2 = pow(avjet, 2); //! average of leading jets pt squared
    
        //!putting thermal particle pushback here 

        const ThermalEvent& thermalEvent = thermalLibrary.randomEvent();
        
        int targetBin = getCentralityBin(thermalEvent.mult, centEdges);
        const vector<PseudoJet>& thermalParticles = thermalEvent.particles;
        
        double v2_thermal = thermalEvent.vn[1];
        /*
        const ThermalEvent *m1Ptr = nullptr, *m2Ptr = nullptr;
        int safety = 0;
        while ((!m1Ptr || !m2Ptr) && safety < 1000) {
            const ThermalEvent& cand = thermalLibrary.randomEvent();
            int candBin = getCentralityBin(cand.mult, centEdges);
            //! the + .001 is for the v2 = 0 case 
            if (candBin == targetBin && fabs(cand.vn[1] - v2_thermal) < 0.05 * (v2_thermal + .0001)) {
                if (!m1Ptr) m1Ptr = &cand;
                else if (!m2Ptr && &cand != m1Ptr) m2Ptr = &cand;
            }
            safety++;
        }
        if (!m1Ptr || !m2Ptr) continue;
        */
        
        //const ThermalEvent& m1Event = *m1Ptr;
        //const ThermalEvent& m2Event = *m2Ptr;
        
        v1.Fill(thermalEvent.vn[0]);
        v2.Fill(thermalEvent.vn[1]);
        v3.Fill(thermalEvent.vn[2]);
        
        //const vector<PseudoJet>& m1 = m1Event.particles;
        //const vector<PseudoJet>& m2 = m2Event.particles;
        const vector<PseudoJet>& m1 = thermalEvent.particles;
        const vector<PseudoJet>& m2 = thermalEvent.particles;
        
        charged_event.insert(charged_event.end(), thermalParticles.begin(), thermalParticles.end());
    
        for (size_t i = 0; i < charged_event.size(); ++i) {
            if (charged_event.at(i).pt() < 1.0) continue;
            // mixed event loop? I just need 1 signal particle so this should work and it doesnt increase the loops
            for (size_t j = i + 1; j < charged_event.size(); ++j) {
                if (charged_event.at(j).pt() < 1.0) continue;
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
            if (m1.at(k).pt() < 1.0) continue;
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
            if (m1.at(i).pt() < 1.0) continue;
                for (size_t k = i + 1; k < m1.size(); ++k) {
                    if (m1.at(k).pt() < 1.0) continue;
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
            if (m1.at(i).pt() < 1.0) continue;
            for (size_t k = 0; k < m2.size(); ++k) {
                if (m2.at(k).pt() < 1.0) continue;
                double mm12_eec = m1[i].pt() * m2[k].pt();
                double mm12_ctheta = costheta(m1[i], m2[k]);
                double mm12_z = (1 - mm12_ctheta) / 2;
                double mm12_delphi = deltaphi(m1[i], m2[k]);
        
                EEC_m2.Fill(mm12_z, mm12_eec / pmq2);
                EEC_m2p.Fill(mm12_delphi, mm12_eec / pmq2);
                EEC_m2pl.Fill(mm12_delphi, mm12_eec / pmq2);
            }
        }
        
        if (iEvent % 10 == 0){
            cout << "event number = " << iEvent << endl;
        }
    }
    
    out->cd();
    
    EEC_w.Write(); EEC_t.Write(); EEC_m.Write(); EEC_s.Write();
    EEC_ms.Write(); EEC_mm.Write(); EEC_m2.Write();
    
    EEC_wp.Write(); EEC_tp.Write(); EEC_mp.Write(); EEC_sp.Write();
    EEC_msp.Write(); EEC_mmp.Write(); EEC_m2p.Write();
    
    EEC_wpl.Write(); EEC_tpl.Write(); EEC_mpl.Write(); EEC_spl.Write();
    EEC_mspl.Write(); EEC_mmpl.Write(); EEC_m2pl.Write();
    
    v1.Write(); v2.Write(); v3.Write();

    hMult.Write();
    
    cout << "Analysis complete. Histograms saved to " << outputFileName << endl;

    
    out->Close();
    in->Close();
}


