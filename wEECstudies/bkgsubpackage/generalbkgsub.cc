//! Authors: Laurynette Griffin & Benjamin Kimmelman

//! Full EEC using TermData partition generator
//! Generates term histograms from partitions 
//! Fills each term histogram for every event and for all needed thermal-block combinations
//! Auto-binning including double log Z bins, linear bins for DR/DPHI


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <cmath>
#include <random>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

//! Helper functions for observables

//! Delta r with eta
double deltaR(const PseudoJet& p1, const PseudoJet& p2) {
  double dphi = fabs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  double deta = p1.eta() - p2.eta();
  return sqrt(deta * deta + dphi * dphi);
}

double costheta(const PseudoJet& p1, const PseudoJet& p2) {
  double dotprod = (p1.px() * p2.px() +  p1.py() * p2.py() +  p1.pz() * p2.pz());
  double normp1 = sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz()); 
  double normp2 = sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz()); 

  if (normp1 == 0 || normp2 == 0) // division by 0 check 
    return 0; 

  return (dotprod / (normp1 * normp2));
}

double deltaphi(const PseudoJet& p1, const PseudoJet& p2) {
  double dphi = fabs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return dphi;
}

//! Loading in Thermalevents in ASCII format

struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;
};

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
                                        /*||
                                        abs(pid) == 11   ||
                                        abs(pid) == 13);
                                        */
                                
                if (!isChargedFinal) continue;
                
                PseudoJet pj(px, py, pz, E);
                pj.set_user_index(pid);
                
                if(pj.pt() < 0.4) continue; //! This pt cut is not based on the paper based on the mult i get out. 
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

//! Thermal struct
struct ThermalLibrary {
    vector<vector<PseudoJet>> bkgevents;
    mutable mt19937 gen;
    mutable uniform_int_distribution<> dist;

    ThermalLibrary() : gen(random_device{}()), dist(0, 0) {}

    void initialize(vector<vector<PseudoJet>>&& input) {
        bkgevents = move(input);
        if (!bkgevents.empty()) {
            dist = uniform_int_distribution<>(0, bkgevents.size() - 1);
        }
    }

    const vector<PseudoJet>& randomEvent() const {
        if (bkgevents.empty()) {
            throw runtime_error("Thermal library empty");
        }
        return bkgevents[dist(gen)];
    }

    size_t size() const { return bkgevents.size(); }
};

ThermalLibrary loadThermalEvents(const vector<string>& filenames)
{
    vector<vector<PseudoJet>> allbkgEvents;

    for (const auto& f : filenames) {
        auto bkgevts = readSingleThermalFile(f);
        allbkgEvents.insert(allbkgEvents.end(),
                         make_move_iterator(bkgevts.begin()),
                         make_move_iterator(bkgevts.end()));
    }

    ThermalLibrary lib;
    lib.initialize(move(allbkgEvents));
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


//! Partition holding backround sub sequence info 
struct TermData {
    string termString;       //! readable strings 
    vector<int> partition;   //! event type combinations 
    int jPower;              //! number of data/Pythia vectors = n - alpha
    int coefficient;         //! coefficient for eec bkg sub
};

int calculate_coefficient(const vector<int>& partition)
{
    if(partition.empty()) return 1;

    vector<int> powerCounts;
    if(!partition.empty()) {
        powerCounts.push_back(1);
        for(size_t i = 1; i < partition.size(); i++) {
            if(partition[i] == partition[i - 1]) powerCounts.back()++;
            else powerCounts.push_back(1);
        }
    }

    int num = 1;
    int den = 1;
    for(int i = 1; i <= (int)partition.size(); i++) num *= i;
    for(int count : powerCounts) for(int i = 1; i <= count; i++) den *= i;

    int sign = (partition.size() % 2 == 0) ? 1 : -1;
    int coeff = sign * (num/den);
    return coeff;
}

string partition_to_string(const vector<int>& partition)
{
    if(partition.empty()) return "";
    string s = "";
    char current_bk = 'A';
    for(int power : partition) {
        for(int i=0;i<power;i++) s.push_back(current_bk);
        current_bk++;
    }
    return s;
}

void generatePartition(int alpha,
                       int remaining_sum,
                       int max_partition_size,
                       vector<int>& current_partition,
                       vector<TermData>& results,
                       const int n,
                       int max_nModels)
{
    if(current_partition.size() > (size_t)max_nModels) return;

    if(remaining_sum == 0) {
        int nJ = n - alpha;
        string j_part(nJ, 'J');
        string b_part = partition_to_string(current_partition);
        int coeff = calculate_coefficient(current_partition);
        results.push_back({j_part + b_part, current_partition, nJ, coeff});
        return;
    }

    for(int i=1; i <= remaining_sum && i <= max_partition_size; i++) {
        current_partition.push_back(i);
        generatePartition(alpha, remaining_sum - i, i, current_partition, results, n, max_nModels);
        current_partition.pop_back();
    }
}

vector<TermData> generate_all_terms(const int n)
{
    vector<TermData> all_terms;
    const int max_nModels = n;

    for(int alpha = 0; alpha <= n; alpha++) {
        if(alpha == 0) {
            all_terms.push_back({string(n,'J'), {}, n, 1});
            continue;
        }

        int remaining_sum = alpha;
        int max_partition_size = alpha;
        vector<int> current_partition;
        generatePartition(alpha, remaining_sum, max_partition_size, current_partition, all_terms, n, max_nModels);
    }

    return all_terms;
}

//! Binning types
enum class EECVariable { Z, DELTA_R, DELTA_PHI };

vector<double> linearBins(int nBins, double min, double max){
    vector<double> bins(nBins+1);
    double step = (max - min)/nBins;
    for(int i=0;i<=nBins;i++) bins[i] = min + i*step;
    return bins;
}

vector<double> logBins(int nBins, double min, double max){
    vector<double> bins(nBins+1);
    double logMin = std::log10(min);
    double logMax = std::log10(max);
    double step = (logMax - logMin)/nBins;
    for(int i=0;i<=nBins;i++) bins[i] = pow(10, logMin + i*step);
    return bins;
}

vector<double> createZLogBins(int nBins, double minVal) {
    //! double log bins for z = (1-cos(theta))/2
    vector<double> bins;
    
    if (minVal <= 0) minVal = 1e-6; 

    int nHalf = nBins / 2;

    double midPoint = 0.5;
    double logMin = log10(minVal);
    double logMax = log10(midPoint);

    //! Generates z <= 0.5 
    for (int i = 0; i <= nHalf; i++) {
        double val = pow(10, logMin + i * (logMax - logMin) / nHalf);
        bins.push_back(val);
    }

    //! Generate z > 0.5 (Mirrored)
    for (int i = bins.size() - 2; i >= 0; i--) {
        double mirroredVal = 1.0 - bins[i];
        bins.push_back(mirroredVal);
    }
    return bins;
}

struct EECResult {
    double ptProduct;
    vector<double> obsValues; //! Contains all N(N-1)/2 distances (Z, Delta_R, or Delta_Phi)
};


EECResult calculateAllData(const vector<fastjet::PseudoJet>& comb, EECVariable var) 
{
    size_t N = comb.size();
    EECResult result = {1.0, {}}; //! Initialize ptProduct to 1.0, maxObsValue empty

    for (const auto& p : comb) {
        result.ptProduct *= p.pt();
    }

    //! Max of the observable between pairs
    for (size_t a = 0; a < N; a++) {
        for (size_t b = a + 1; b < N; b++) {
            
            double obsValue = 0;
            const fastjet::PseudoJet& pA = comb[a];
            const fastjet::PseudoJet& pB = comb[b];

            
            if (var == EECVariable::DELTA_R) {
                obsValue = deltaR(pA, pB);
            }
            else if (var == EECVariable::Z) {
                // Apply Z transformation (1 - cos(angle))/2
                obsValue = (1.0 - costheta(pA,pB)) / 2.0;
            }
            else { // DELTA_PHI
                obsValue = deltaphi(pA,pB);
            }
            result.obsValues.push_back(obsValue); 
        }
        
    }

    if (N == 3) {
        sort(result.obsValues.begin(), result.obsValues.end());
    }
     return result;
}


//! Fill a single  histogram by recursively selecting one particle from each slot
void recursiveFillTerm(int level,
                       int N_totalSlots,
                       const vector<const vector<fastjet::PseudoJet>*>& vectors,
                       vector<int>& indices,
                       TH1D* hist,
                       EECVariable var)
    {
    if(level == N_totalSlots){
        vector<fastjet::PseudoJet> comb;
        comb.reserve(N_totalSlots);
        for(int i=0;i<N_totalSlots;i++) comb.push_back((*vectors[i])[indices[i]]);
        
        EECResult data = calculateAllData(comb, var);
        double weight = data.ptProduct; // Use raw pT product as the weight
        //! add a normalization by leading and subleading jet pt later 
        
        if (N_totalSlots == 3) { //!!!! this should be a vector of vectors for the distance with short med long as 1 and 0 is ptproduct
            hist->Fill(data.obsValues[2], weight); // Long
        }
        else {
            
            for (double distance : data.obsValues) {
                    hist->Fill(distance, weight);
                }
        }

    return;
    }
    
    int start_index = 0;
    //! Unique starting index for loops over same event
    if (level > 0 && vectors[level] == vectors[level-1]) {
        start_index = indices[level-1] + 1;
    }

    auto &vec = *vectors[level];
    for(size_t i = start_index; i < vec.size(); i++){
        indices[level] = (int)i;
        recursiveFillTerm(level + 1, N_totalSlots, vectors, indices, hist, var); 
    }
}

//! Helper to generate block-index combinations
void generateBlockIndexCombos(int blocks, 
                              int nThermalEvents, 
                              vector<vector<int>>& outCombos, 
                              vector<int>& current, 
                              int level=0) 
{ 
    if(level == blocks){ 
        outCombos.push_back(current); 
        return; 
    } 
    int start = (level==0) ? 0 : current[level-1]; // non-decreasing 
    for(int idx = start; idx < nThermalEvents; ++idx){ 
        current[level] = idx; 
        generateBlockIndexCombos(blocks, nThermalEvents, outCombos, current, level+1); 
    } 
}

 int main(int argc, char* argv[]) {
   //! <PythiaROOTFile> <N> <OutputName> <StartEvent number in pythia> <EndEvent number in pythia> <ThermalFile1> ...
   //! arg 3 & 4 allow for all signal(J) events to be run over when submiting to the grid :)
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <PythiaROOTFile> <N-point> <OutputName> <StartEvent> <EndEvent> <ThermalFile1> [ThermalFile2] ..." << endl;
        return 1;
    }

    // Command line arguments
    string pythiaFile = argv[1];
    int N = stoi(argv[2]);
    string outputBaseName = argv[3];
    
    // ARGUMENTS FOR PARALLELIZATION
    long int iEventStart = stol(argv[4]); // Start processing from this event index
    long int iEventEnd = stol(argv[5]);   // Stop processing at this event index (exclusive)
    
    // COLLECT ALL THERMAL FILES from argv[6] onwards (Indices shifted)
    vector<string> thermalFiles;
    for (int i = 6; i < argc; ++i) {
        thermalFiles.push_back(argv[i]);
    }
    EECVariable var = EECVariable::Z;
    int bins = 50;

    //! binning Setup 
    vector<double> binEdges;
    string varName = "";
    if (var == EECVariable::Z) {
        varName = "Z";
        binEdges = createZLogBins(bins, 1e-4);
    } 
    else if (var == EECVariable::DELTA_R) {
        varName = "DR";
        binEdges = linearBins(bins, 0.0, 6.5);
    } 
    else {
        varName = "DPHI";
        binEdges = linearBins(bins, 0.0, M_PI);
    }

    PythiaTreeReader* pythiaReader = nullptr;
    ThermalLibrary thermalLibrary;

    try {
        pythiaReader = new PythiaTreeReader(pythiaFile);
        thermalLibrary = loadThermalEvents(thermalFiles);
    } 
    catch (const runtime_error& e) {
        cerr << "Initialization Error: " << e.what() << endl;
        return 1;
    }

    // Set number of events to process
    long nEvents = pythiaReader->entries();
    int nThermalPool = (int)thermalLibrary.size();

    if (nThermalPool == 0) {
        cerr << "ERROR: Thermal library is empty. Cannot perform background subtraction." << endl;
        delete pythiaReader;
        return 1;
    }
    
    // --- Term Generation and Histogram Initialization ---
    vector<TermData> terms = generate_all_terms(N);
    vector<TH1D*> termHists;
    // ... (Histogram initialization logic as before, including Sumw2) ...
    for (size_t i = 0; i < terms.size(); i++) {
        string histName = "h_" + varName + "_" + terms[i].termString;
        TH1D* h = new TH1D(histName.c_str(), 
                           (terms[i].termString + " (" + varName + ")").c_str(), 
                           binEdges.size()-1, binEdges.data());
        h->Sumw2();
        termHists.push_back(h);
    }


    // --- Event Loop ---
    for (long int iEvent = iEventStart; iEvent < iEventEnd; iEvent++) {
        if(iEvent % 100 == 0) cout << "Processing Event " << iEvent << "/" << nEvents << endl;
        
        // Load the Pythia (Signal) event 
        const vector<PseudoJet> eventJ = pythiaReader->getEvent(iEvent);

        // Loop Over All Terms
        for (size_t t = 0; t < terms.size(); t++) {
            TermData& term = terms[t];
            
            int nBlocks = term.partition.size();
            vector<vector<int>> thermalCombos;
            vector<int> currentCombo(nBlocks);
            
            if (nBlocks > 0) {
                generateBlockIndexCombos(nBlocks, nThermalPool, thermalCombos, currentCombo);
            } else {
                thermalCombos.push_back({}); // Pure J-term case
            }

            // Iterate over all valid thermal event combinations
            for (const auto& combo : thermalCombos) {
                
                vector<const vector<fastjet::PseudoJet>*> inputVectors;
                vector<int> indices(N);

                // Add J vectors (signal event)
                // We use the const_cast to add the address of the const vector
                for (int k = 0; k < term.jPower; k++) {
                    inputVectors.push_back(&eventJ);
                }
                
                // Add Background vectors
                for(size_t b=0; b < term.partition.size(); b++) {
                    int blockSize = term.partition[b];
                    int eventIndex = combo[b]; 
                    
                    const vector<PseudoJet>& thermalEvent = thermalLibrary.bkgevents[eventIndex];
                    // Push the pointer 'blockSize' times
                    for(int r=0; r<blockSize; r++) {
                        inputVectors.push_back(&thermalEvent);
                    }
                }

                // Run Recursion
                recursiveFillTerm(0, N, inputVectors, indices, termHists[t], var);
            }
        }
    }
    
    // --- Finalization and Background Subtraction ---
    
    // ... (Background Subtraction logic as before) ...
    TH1D* hbkgsub = (TH1D*)termHists[0]->Clone("hbkgsub_Final");
    hbkgsub->Reset();

    for (size_t t = 0; t < terms.size(); t++) {
        double coeff = terms[t].coefficient;
        hbkgsub->Add(termHists[t], coeff); // Corrected: Use the signed coefficient directly
        cout << "Term: " << terms[t].termString << " | Coeff: " << coeff << " | Entries: " << termHists[t]->GetEntries() << endl;
    }

    // Save to file
   string outputFilename = outputBaseName + "_" + varName + "_EEC_N" + to_string(N) 
                            + "_E" + to_string(iEventStart) + "_to_" + to_string(iEventEnd) 
                            + "_Results.root";
    TFile* f = new TFile(outputFilename.c_str(), "RECREATE");
    hbkgsub->Write();
    for(auto h : termHists) h->Write();
    f->Close();

    // Clean up readers
    delete pythiaReader;
    
    cout << "Done. Results saved to " << outputFilename << endl;
    return 0;
}