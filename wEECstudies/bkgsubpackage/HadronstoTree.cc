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

vector<vector<fastjet::PseudoJet>>
readSingleThermalFile(const string& filename)
{
    vector<vector<PseudoJet>> events;

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open thermal file: " << filename << endl;
        return events;
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

            vector<PseudoJet> event;
            event.reserve(n_hadrons);

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

        
                //! take out 22 if want true em 
                bool isChargedFinal = (abs(pid) ==  211  ||
                                        abs(pid) == 321  || 
                                        abs(pid) == 2212 );
                                        /*||
                                        abs(pid) == 11   ||
                                        abs(pid) == 13);
                                        */
                                
                if (!isChargedFinal) continue;
                
                PseudoJet pj(px, py, pz, E);
                if (pj.pt() <  0.4) continue; 
                if (abs(pj.eta()) > 1.1 ) continue;
                pj.set_user_index(pid);

                event.push_back(pj);
            }

            if (!event.empty()) {
                events.push_back(move(event));
            }
        }
    }

    cout << "  -> Loaded " << events.size()
         << " events from " << filename << endl;

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

        double bkgavg = 0;
        for (long i = evtStart; i < evtEnd; ++i) {
            auto signal = pythia.getEvent(i);
            auto background = thermal.randomEvent();

            bkgavg += background.size();
            cout << "Event " << i
                 //<< " | signal = " << signal.size()
                 << " | thermal = " << background.size()
                 << endl;

            if (i == 1) {
                cout << "Background PIDs (event 1):" << endl;
        
                for (size_t j = 0; j < background.size(); ++j) {
                    cout << " j=" << j
                         << " pid=" << background[j].user_index()
                         << endl;
                }
                
            }
            
            // >>> N-point correlator / subtraction happens here <<<
        }
        cout << "total bkg " << bkgavg << endl;
    }
    catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return 1;
    }
    cout << "yay" << endl;
    return 0;
}
