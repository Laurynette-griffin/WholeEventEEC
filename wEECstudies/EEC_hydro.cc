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

double trackptcut = 0.5;

//! Delta R calculation (uses eta)
double deltaR(const PseudoJet& p1, const PseudoJet& p2) {
    double dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    double deta = p1.eta() - p2.eta();
    return std::sqrt(deta * deta + dphi * dphi);
}

//! costheta calculation
double costheta(const PseudoJet& p1, const PseudoJet& p2) {
    double dotprod = (p1.px() * p2.px() + p1.py() * p2.py() + p1.pz() * p2.pz());
    double normp1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz());
    double normp2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz());

    if (normp1 == 0 || normp2 == 0) return 0;
    return (dotprod / (normp1 * normp2));
}

//! delphi calculation 
double deltaphi(const PseudoJet& p1, const PseudoJet& p2) {
    double dphi = std::abs(p1.phi() - p2.phi());
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    return dphi;
}

int getCentralityBin(int mult, const vector<double>& edges) {
    for (int i = 0; i < edges.size(); ++i) {
        //! edges[i] is the upper bound, edges[i+1] is the lower bound
        if (mult <= edges[i] && mult > edges[i + 1]) return i;
    }
    return -1;
}

struct ThermalEvent {
    vector<PseudoJet> particles;
    double header_EPangle; //! True event plane angle from file
    TComplex Qn[3];        //! Q-vectors for n=1, 2, 3
    double psi_n[3];       //! Calculated event plane angles
    double vn[3];          //! Event-averaged v1, v2, v3
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
                if (pj.pt() < trackptcut || fabs(pj.eta()) > 1.1) {
                    ev.particles.pop_back();
                }
            }

            if (!ev.particles.empty()) {
                ev.mult = ev.particles.size();

                //! Calculate Q-vectors and vn for the loaded event
                for (int iH = 0; iH < 3; iH++) {
                    int n = iH + 1;
                    ev.Qn[iH] = TComplex(0, 0);
                    for (const auto& p : ev.particles) {
                        ev.Qn[iH] += TComplex(cos(n * p.phi()), sin(n * p.phi()));
                    }

                    ev.psi_n[iH] = TMath::ATan2(ev.Qn[iH].Im(), ev.Qn[iH].Re()) / n;
                    if (ev.psi_n[iH] < 0) ev.psi_n[iH] += 2 * TMath::Pi() / n;
                    double sum_vn = 0;
                    for (const auto& p : ev.particles) sum_vn += cos(n * (p.phi() - ev.psi_n[iH]));
                    ev.vn[iH] = sum_vn / ev.particles.size();
                }
                bkgevents.push_back(move(ev));
            }
        }
    }
    return bkgevents;
}

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

ThermalLibrary loadThermalEvents(const vector<string>& filenames) {
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
    vector<float>* e = nullptr;
    vector<int>* id = nullptr;
    vector<int>* status = nullptr;
    long nEntries = 0;

public:
    PythiaTreeReader(const string& filename) {
        file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie())
            throw runtime_error("Cannot open Pythia ROOT file");

        tree = dynamic_cast<TTree*>(file->Get("events"));
        if (!tree)
            throw runtime_error("Cannot find TTree 'events'");

        tree->SetBranchAddress("px", &px);
        tree->SetBranchAddress("py", &py);
        tree->SetBranchAddress("pz", &pz);
        tree->SetBranchAddress("e", &e);
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

    vector<PseudoJet> getEvent(long i) {
        if (i < 0 || i >= nEntries)
            throw out_of_range("Event index out of range");

        tree->GetEntry(i);
        vector<PseudoJet> signal_events;
        signal_events.reserve(px->size());

        for (size_t j = 0; j < px->size(); ++j) {
            PseudoJet pj(px->at(j), py->at(j), pz->at(j), e->at(j));
            pj.set_user_index(id->at(j));
            signal_events.push_back(pj);
        }
        return signal_events;
    }
};

//! Just copied it and made a jet tree reader class :) 
class PythiaJetTreeReader {
private:
    TFile* file = nullptr;
    TTree* tree = nullptr;
    vector<float>* jet_px = nullptr;
    vector<float>* jet_py = nullptr;
    vector<float>* jet_pz = nullptr;
    vector<float>* jet_e = nullptr;
    vector<float>* jet_pt = nullptr;
    vector<int>* jet_id = nullptr;
    long nEntries = 0;

public:
    PythiaJetTreeReader(const string& filename, const string& ttreeName) {
        file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie())
            throw runtime_error("Cannot open Pythia ROOT file");

        tree = dynamic_cast<TTree*>(file->Get(ttreeName.c_str()));
        if (!tree)
            throw runtime_error("Cannot find TTree: " + ttreeName);

        tree->SetBranchAddress("jet_px", &jet_px);
        tree->SetBranchAddress("jet_py", &jet_py);
        tree->SetBranchAddress("jet_pz", &jet_pz);
        tree->SetBranchAddress("jet_e", &jet_e);
        tree->SetBranchAddress("jet_pt", &jet_pt);
        tree->SetBranchAddress("jet_id", &jet_id);

        nEntries = tree->GetEntries();
        cout << "Opened jet tree '" << ttreeName << "' with " << nEntries << " events" << endl;
    }

    ~PythiaJetTreeReader() {
        if (file) {
            file->Close();
            delete file;
            file = nullptr;
        }
    }

    long entries() const { return nEntries; }

    //! Returns vector of jets with pt in [ptMin, ptMax]
    vector<PseudoJet> getJets(long iEvent, double ptMin = 20.0, double ptMax = 30.0) {
        if (iEvent < 0 || iEvent >= nEntries)
            throw out_of_range("Event index out of range");

        tree->GetEntry(iEvent);
        vector<PseudoJet> selectedJets;
        for (size_t i = 0; i < jet_px->size(); ++i) {
            double pt = jet_pt->at(i);
            if (pt < ptMin || pt > ptMax) continue;

            PseudoJet jet(jet_px->at(i), jet_py->at(i), jet_pz->at(i), jet_e->at(i));
            jet.set_user_index(jet_id->at(i));
            selectedJets.push_back(jet);
        }
        return selectedJets;
    }
};

//allows for the .cmnd file to be accepted as an argument 
int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <config.cmnd> <output.root> <PythiaROOTFile> <iStart> <iEnd> <thermalFiles>" << endl;
        return 1;
    }

    string configFile = argv[1];
    string outputFileName = argv[2];
    string pythiaFile = argv[3];
    long int iEventStart = stol(argv[4]);
    long int iEventEnd = stol(argv[5]);

    vector<string> thermalFiles;
    for (int i = 6; i < argc; ++i) {
        thermalFiles.push_back(argv[i]);
    }

    string treename = "jets_R04_antikt_s";
    PythiaJetTreeReader jetReader(pythiaFile, treename);

    TFile *out = new TFile(outputFileName.c_str(), "RECREATE");
    out->cd();

    ThermalLibrary thermalLibrary;
    try {
        thermalLibrary = loadThermalEvents(thermalFiles);
    } catch (const runtime_error& e) {
        cerr << "Initialization Error: " << e.what() << endl;
        return 1;
    }

    int nThermalPool = (int)thermalLibrary.size();
    if (nThermalPool == 0) {
        cerr << "ERROR: Thermal library is empty." << endl;
        return 1;
    }

    //! Multiplicity finder for these hydro events
    vector<int> allMults;
    for (const auto& ev : thermalLibrary.events) {
        allMults.push_back(ev.mult);
    }

    sort(allMults.begin(), allMults.end(), greater<int>());
    int min_mult = allMults[allMults.size() - 1];
    int max_mult = allMults[0];

    TH1D hMult("hMult", "Charged Multiplicity;N_{ch};Counts", 25, min_mult, max_mult);
    for (const auto& ev : thermalLibrary.events) hMult.Fill(ev.mult);
    hMult.Scale(1.0 / hMult.Integral());

    int nBins = 20;
    vector<double> centEdges;
    centEdges.push_back(hMult.GetBinLowEdge(hMult.GetNbinsX() + 1));

    double runningArea = 0.0;
    int currentTargetbin = 1;

    for (int i = hMult.GetNbinsX(); i >= 1; --i) {
        double binContent = hMult.GetBinContent(i);
        double areaBeforeBin = runningArea;
        runningArea += binContent;
        double targetArea = (double)currentTargetbin / nBins;

        if (runningArea >= targetArea) {
            double distCurrent = fabs(runningArea - targetArea);
            double distPrev = fabs(areaBeforeBin - targetArea);
            double chosenEdge = (distCurrent < distPrev) ? hMult.GetBinLowEdge(i) : hMult.GetBinLowEdge(i + 1);
            centEdges.push_back(chosenEdge);
            currentTargetbin++;
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

    Double_t delredge[61] = {1.000000000000e-05, 1.238977277499e-05, 1.535064694159e-05, 1.901910275555e-05, 2.356423615255e-05,
                            2.919555315464e-05, 3.617262696261e-05, 4.481706287414e-05, 5.552732254531e-05, 6.879709091401e-05,
                            8.523803240051e-05, 1.056079853230e-04, 1.308458941376e-04, 1.621150896906e-04, 2.008569124664e-04,
                            2.488571505746e-04, 3.083283549051e-04, 3.820118257361e-04, 4.733039718231e-04, 5.864128664390e-04,
                            7.265522167511e-04, 9.001816874713e-04, 1.115304656398e-03, 1.381837126766e-03, 1.712064801268e-03,
                            2.121209386377e-03, 2.628130230540e-03, 3.256193637948e-03, 4.034349928555e-03, 4.998467890960e-03,
                            6.192988139209e-03,
                            7.672971584303e-03, 9.506637443848e-03, 1.177850777835e-02, 1.459330350023e-02, 1.808077144043e-02,
                            2.240166497435e-02, 2.775515388137e-02, 3.438800499252e-02, 4.260595680426e-02, 5.278781236659e-02,
                            6.540290005110e-02, 8.103270704587e-02, 1.003976827641e-01, 1.243904476583e-01, 1.541169381866e-01,
                            1.909473844909e-01, 2.365794705822e-01, 2.931165883741e-01, 3.631647926536e-01, 4.499529260856e-01,
                            5.574814513643e-01, 6.907068508677e-01, 8.557700936381e-01, 1.060279700781e+00, 1.313662457061e+00,
                            1.627597934603e+00, 2.016556857878e+00, 2.498468125696e+00, 3.095545236293e+00, 3.835310209239e+00};


    int bins = 60;
    TH1::SetDefaultSumw2();

    TH1D EEC_w("EEC_w", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_t("EEC_t", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_m("EEC_m", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_s("EEC_s", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_ms("EEC_ms", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_mm("EEC_mm", "Energy Energy Correlator", bins, eecbounds);
    TH1D EEC_m2("EEC_m2", "Energy Energy Correlator", bins, eecbounds);

    TH1D EEC_wr("EEC_wr", "EEC DeltaR", bins, delredge);
    TH1D EEC_tr("EEC_tr", "EEC DeltaR", bins, delredge);
    TH1D EEC_mr("EEC_mr", "EEC DeltaR", bins, delredge);
    TH1D EEC_sr("EEC_sr", "EEC DeltaR", bins, delredge);
    TH1D EEC_msr("EEC_msr", "EEC DeltaR", bins, delredge);
    TH1D EEC_mmr("EEC_mmr", "EEC DeltaR", bins, delredge);
    TH1D EEC_m2r("EEC_m2r", "EEC DeltaR", bins, delredge);

    TH1D EEC_wrl("EEC_wrl", "EEC DeltaR Linear", bins, 0.0, 3.8354);
    TH1D EEC_trl("EEC_trl", "EEC DeltaR Linear", bins, 0.0, 3.8354);
    TH1D EEC_mrl("EEC_mrl", "EEC DeltaR Linear", bins, 0.0, 3.8354);
    TH1D EEC_srl("EEC_srl", "EEC DeltaR Linear", bins, 0.0, 3.8354);
    TH1D EEC_msrl("EEC_msrl", "EEC DeltaR Linear", bins, 0.0, 3.8354);
    TH1D EEC_mmrl("EEC_mmrl", "EEC DeltaR Linear", bins, 0.0, 3.8354);
    TH1D EEC_m2rl("EEC_m2rl", "EEC DeltaR Linear", bins, 0.0, 3.8354);

    TH1D EEC_wp("EEC_wp", "EEC Phi", bins, topi);
    TH1D EEC_tp("EEC_tp", "EEC Phi", bins, topi);
    TH1D EEC_mp("EEC_mp", "EEC Phi", bins, topi);
    TH1D EEC_sp("EEC_sp", "EEC Phi", bins, topi);
    TH1D EEC_msp("EEC_msp", "EEC Phi", bins, topi);
    TH1D EEC_mmp("EEC_mmp", "EEC Phi", bins, topi);
    TH1D EEC_m2p("EEC_m2p", "EEC Phi", bins, topi);

    TH1D EEC_wpl("EEC_wpl", "EEC Phi Linear", bins, 0, M_PI);
    TH1D EEC_tpl("EEC_tpl", "EEC Phi Linear", bins, 0, M_PI);
    TH1D EEC_mpl("EEC_mpl", "EEC Phi Linear", bins, 0, M_PI);
    TH1D EEC_spl("EEC_spl", "EEC Phi Linear", bins, 0, M_PI);
    TH1D EEC_mspl("EEC_mspl", "EEC Phi Linear", bins, 0, M_PI);
    TH1D EEC_mmpl("EEC_mmpl", "EEC Phi Linear", bins, 0, M_PI);
    TH1D EEC_m2pl("EEC_m2pl", "EEC Phi Linear", bins, 0, M_PI);

    TH1D psit("psit", "#Psi_{2} hydro 1", 63, 0.0, 2 * M_PI);
    TH1D psim1("psim1", "#Psi_{2} hydro 2", 63, 0.0, 2 * M_PI);
    TH1D psim2("psim2", "#Psi_{2} hydro 3", 63, 0.0, 2 * M_PI);

    long int nEvents = jetReader.entries();
    if (iEventEnd > nEvents) iEventEnd = nEvents;
    
    //! these are the min in the pythiajettreereader function anyway 
    double minJetPt = 20.0;
    double maxJetPt = 30.0;

    for (long int iEvent = iEventStart; iEvent < iEventEnd; ++iEvent) {
        if (iEvent % 100 == 0) cout << "Processing Event " << iEvent << "/" << iEventEnd << endl;

        vector<PseudoJet> selectedJets = jetReader.getJets(iEvent, minJetPt, maxJetPt);
        if (selectedJets.empty()) continue;

        vector<fastjet::PseudoJet> charged_event;
        for (const auto& jet : selectedJets) {
            for (auto& p : jet.constituents()) {
                fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
                pj.set_user_index(-1); 
                charged_event.push_back(pj);
            }
        }

        const ThermalEvent& thermalEvent = thermalLibrary.randomEvent();
        int targetBin = getCentralityBin(thermalEvent.mult, centEdges);
        const vector<PseudoJet>& thermalParticles = thermalEvent.particles;
        double psin_thermal = thermalEvent.psi_n[1];

        const ThermalEvent *m1Ptr = nullptr, *m2Ptr = nullptr;
        int safety = 0;
        while ((!m1Ptr || !m2Ptr) && safety < 1000) {
            const ThermalEvent& cand = thermalLibrary.randomEvent();
            int candBin = getCentralityBin(cand.mult, centEdges);
            double dpsi = fabs(cand.psi_n[1] - psin_thermal);
            if (dpsi > M_PI / 2) dpsi = M_PI - dpsi;

            if (candBin == targetBin && dpsi < 0.05 * (M_PI / 2)) {
                if (!m1Ptr) m1Ptr = &cand;
                else if (!m2Ptr && &cand != m1Ptr) m2Ptr = &cand;
            }
            safety++;
        }
        if (!m1Ptr || !m2Ptr) continue;

        const ThermalEvent& m1Event = *m1Ptr;
        const ThermalEvent& m2Event = *m2Ptr;
        psit.Fill(thermalEvent.psi_n[0]);
        psim1.Fill(m1Event.psi_n[1]);
        psim2.Fill(m2Event.psi_n[2]);

        const vector<PseudoJet>& m1 = m1Event.particles;
        const vector<PseudoJet>& m2 = m2Event.particles;

        charged_event.insert(charged_event.end(), thermalParticles.begin(), thermalParticles.end());

        double pmq2 = 1.0;
        for (const auto& jet : selectedJets) {
            vector<PseudoJet> jetConeParticles;
            for (const auto& p : charged_event) {
                if (deltaR(jet, p) < 0.4) {
                    jetConeParticles.push_back(p);
                }
            }

            for (size_t i = 0; i < jetConeParticles.size(); ++i) {
                if (jetConeParticles[i].pt() < trackptcut) continue;
                for (size_t j = i + 1; j < jetConeParticles.size(); ++j) {
                    if (jetConeParticles[j].pt() < trackptcut) continue;

                    double eec = jetConeParticles[i].pt() * jetConeParticles[j].pt();
                    double ctheta = costheta(jetConeParticles[i], jetConeParticles[j]);
                    double z = (1 - ctheta) / 2;
                    double delphi = deltaphi(jetConeParticles[i], jetConeParticles[j]);
                    double delr = deltaR(jetConeParticles[i], jetConeParticles[j]);

                    int ui = jetConeParticles[i].user_index();
                    int uj = jetConeParticles[j].user_index();

                    EEC_w.Fill(z, eec / pmq2);
                    EEC_wp.Fill(delphi, eec / pmq2);
                    EEC_wpl.Fill(delphi, eec / pmq2);
                    EEC_wr.Fill(delr, eec / pmq2);
                    EEC_wrl.Fill(delr, eec / pmq2);

                    if (ui == 1 && uj == 1) { //! thermal
                        EEC_t.Fill(z, eec / pmq2);
                        EEC_tp.Fill(delphi, eec / pmq2);
                        EEC_tpl.Fill(delphi, eec / pmq2);
                        EEC_tr.Fill(delr, eec / pmq2);
                        EEC_trl.Fill(delr, eec / pmq2);
                    }

                    if (ui == -1 && uj == -1) { //! signal
                        EEC_s.Fill(z, eec / pmq2);
                        EEC_sp.Fill(delphi, eec / pmq2);
                        EEC_spl.Fill(delphi, eec / pmq2);
                        EEC_sr.Fill(delr, eec / pmq2);
                        EEC_srl.Fill(delr, eec / pmq2);
                    }

                    if ((ui == -1 && uj == 1) || (ui == 1 && uj == -1)) { //! mixed
                        EEC_m.Fill(z, eec / pmq2);
                        EEC_mp.Fill(delphi, eec / pmq2);
                        EEC_mpl.Fill(delphi, eec / pmq2);
                        EEC_mr.Fill(delr, eec / pmq2);
                        EEC_mrl.Fill(delr, eec / pmq2);
                    }
                }

                for (size_t k = 0; k < m1.size(); ++k) {
                    if (m1.at(k).pt() < trackptcut) continue;
                    if (deltaR(m1.at(k), jet) > 0.4) continue;
                    double sm_eec = jetConeParticles.at(i).pt() * m1.at(k).pt();
                    double sm_ctheta = costheta(jetConeParticles.at(i), m1.at(k));
                    double sm_z = (1 - sm_ctheta) / 2;
                    double sm_delphi = deltaphi(jetConeParticles.at(i), m1.at(k));
                    double sm_delr = deltaR(jetConeParticles.at(i), m1.at(k));

                    EEC_ms.Fill(sm_z, sm_eec / pmq2);
                    EEC_msp.Fill(sm_delphi, sm_eec / pmq2);
                    EEC_mspl.Fill(sm_delphi, sm_eec / pmq2);
                    EEC_msr.Fill(sm_delr, sm_eec / pmq2);
                    EEC_msrl.Fill(sm_delr, sm_eec / pmq2);
                }
            }//! closes the jet + hydro loop 

            for (size_t i = 0; i < m1.size(); ++i) {
                if (m1.at(i).pt() < trackptcut) continue;
                if (deltaR(m1.at(i), jet) > 0.4) continue;
                for (size_t k = i + 1; k < m1.size(); ++k) {
                    if (m1.at(k).pt() < trackptcut) continue;
                    if (deltaR(m1.at(k), jet) > 0.4) continue;
                    double mm_eec = m1[i].pt() * m1[k].pt();
                    double mm_ctheta = costheta(m1[i], m1[k]);
                    double mm_z = (1 - mm_ctheta) / 2;
                    double mm_delphi = deltaphi(m1[i], m1[k]);
                    double mm_delr = deltaR(m1[i], m1[k]);

                    EEC_mm.Fill(mm_z, mm_eec / pmq2);
                    EEC_mmp.Fill(mm_delphi, mm_eec / pmq2);
                    EEC_mmpl.Fill(mm_delphi, mm_eec / pmq2);
                    EEC_mmr.Fill(mm_delr, mm_eec / pmq2);
                    EEC_mmrl.Fill(mm_delr, mm_eec / pmq2);
                }
            }

            for (size_t i = 0; i < m1.size(); ++i) {
                if (m1.at(i).pt() < trackptcut) continue;
                if (deltaR(m1.at(i), jet) > 0.4) continue;
                for (size_t k = 0; k < m2.size(); ++k) {
                    if (m2.at(k).pt() < trackptcut) continue;
                    if (deltaR(m2.at(k), jet) > 0.4) continue;
                    double m2_eec = m1[i].pt() * m2[k].pt();
                    double m2_ctheta = costheta(m1[i], m2[k]);
                    double m2_z = (1 - m2_ctheta) / 2;
                    double m2_delphi = deltaphi(m1[i], m2[k]);
                    double m2_delr = deltaR(m1[i], m2[k]);

                    EEC_m2.Fill(m2_z, m2_eec / pmq2);
                    EEC_m2p.Fill(m2_delphi, m2_eec / pmq2);
                    EEC_m2pl.Fill(m2_delphi, m2_eec / pmq2);
                    EEC_m2r.Fill(m2_delr, m2_eec / pmq2);
                    EEC_m2rl.Fill(m2_delr, m2_eec / pmq2);
                }
            }
        }//! close the jet loop :)
    }

    out->cd();
    
    EEC_w.Write(); EEC_t.Write(); EEC_m.Write(); EEC_s.Write();
    EEC_ms.Write(); EEC_mm.Write(); EEC_m2.Write();
    
    EEC_wp.Write(); EEC_tp.Write(); EEC_mp.Write(); EEC_sp.Write();
    EEC_msp.Write(); EEC_mmp.Write(); EEC_m2p.Write();
    
    EEC_wr.Write(); EEC_tr.Write(); EEC_mr.Write(); EEC_sr.Write();
    EEC_msr.Write(); EEC_mmr.Write(); EEC_m2r.Write();
    
    EEC_wrl.Write(); EEC_trl.Write(); EEC_mrl.Write(); EEC_srl.Write();
    EEC_msrl.Write(); EEC_mmrl.Write(); EEC_m2rl.Write();
    
    EEC_wpl.Write(); EEC_tpl.Write(); EEC_mpl.Write(); EEC_spl.Write();
    EEC_mspl.Write(); EEC_mmpl.Write(); EEC_m2pl.Write();
    
    psit.Write(); psim1.Write(); psim2.Write();
    
    hMult.Write();

    cout << "Analysis complete. Histograms saved to " << outputFileName << endl;
    out->Close();
    return 0;
}