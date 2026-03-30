//! Author: Laurynette Griffin 
//! Date: Dec 9 2025

//! Takes ASCII file and makes the particle info into a TTREE
//! Assumes the input file format for each hadron is: Px Py Pz E PID

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;
};

vector<vector<PseudoJet>>
readSingleThermalFile(const string& filename)
{
    vector<vector<PseudoJet>> events;

    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "ERROR: Cannot open thermal file: " << filename << endl;
        return events;
    }

    cout << "Reading thermal file: " << filename << endl;

    string line;
    vector<PseudoJet> currentEvent;

    while (getline(inputFile, line)) {
        stringstream ss(line);
        string token;
        ss >> token;

        if (token == "EVENT") {
            currentEvent.clear();
        }
        else if (token == "N_HADRONS") {
            int nHadrons = 0;
            ss >> nHadrons;

            currentEvent.reserve(nHadrons);

            for (int i = 0; i < nHadrons; ++i) {
                if (!getline(inputFile, line)) {
                    cerr << "ERROR: Unexpected EOF in thermal file" << endl;
                    break;
                }

                double px, py, pz, e;
                int pid;

                stringstream hs(line);
                if (!(hs >> px >> py >> pz >> e >> pid)) {
                    cerr << "ERROR parsing thermal particle: " << line << endl;
                    continue;
                }

                PseudoJet pj(px, py, pz, e);
                pj.set_user_index(pid);
                currentEvent.push_back(pj);
            }

            if (!currentEvent.empty()) {
                events.push_back(currentEvent);
            }
        }
    }

    cout << "  -> Loaded " << events.size()
         << " thermal events from " << filename << endl;

    return events;
}

//! Thermal struct
struct ThermalLibrary {
    vector<vector<PseudoJet>> events;
    mutable mt19937 gen;
    mutable uniform_int_distribution<> dist;

    ThermalLibrary() : gen(random_device{}()), dist(0, 0) {}

    void initialize(vector<vector<PseudoJet>>&& input) {
        events = move(input);
        if (!events.empty()) {
            dist = uniform_int_distribution<>(0, events.size() - 1);
        }
    }

    const vector<PseudoJet>& randomEvent() const {
        if (events.empty()) {
            throw runtime_error("Thermal library empty");
        }
        return events[dist(gen)];
    }

    size_t size() const { return events.size(); }
};

ThermalLibrary loadThermalEvents(const vector<string>& filenames)
{
    vector<vector<PseudoJet>> allEvents;

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

//! Pythia reader
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
        if (file) file->Close();
        delete file;
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
            if (status->at(j) != 1) continue; // final-state only

            PseudoJet pj(px->at(j), py->at(j),
                          pz->at(j), e->at(j));
            pj.set_user_index(id->at(j));
            signal_events.push_back(pj);
        }
        return signal_events;
    }
};

int main(int argc, char* argv[])
{
    if (argc < 7) {
        cerr << "Usage: " << argv[0]
             << " <pythia.root> <N> <out> <startEvt> <endEvt> <thermal1> [thermal2 ...]"
             << endl;
        return 1;
    }

    string pythiaFile = argv[1];
    int N = stoi(argv[2]);
    string outName = argv[3];
    long evtStart = stol(argv[4]);
    long evtEnd   = stol(argv[5]);

    vector<string> thermalFiles;
    for (int i = 6; i < argc; ++i)
        thermalFiles.push_back(argv[i]);

    try {
        PythiaTreeReader pythia(pythiaFile);
        ThermalLibrary thermal = loadThermalEvents(thermalFiles);

        if (thermal.size() == 0) {
            cerr << "ERROR: Thermal library empty" << endl;
            return 1;
        }

        cout << "Thermal pool size: " << thermal.size() << endl;

        for (long i = evtStart; i < evtEnd; ++i) {
            auto signal = pythia.getEvent(i);
            auto background = thermal.randomEvent();

            cout << "Event " << i
                 << " | signal = " << signal.size()
                 << " | thermal = " << background.size()
                 << endl;

            // >>> N-point correlator / subtraction happens here <<<
        }
    }
    catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return 1;
    }
    cout << "yay" << endl;
    return 0;
}
