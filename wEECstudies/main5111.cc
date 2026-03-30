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
  double dphi = abs(p1.phi() - p2.phi());
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
    double dphi = abs(p1.phi() - p2.phi());
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
    
vector<PseudoJet> thermalGen(int nParticles = 1500){
    
    TF1 *mb = new TF1("mb","exp(-[0]*x/[1])",0.3,10); // pt range from .3 to 10 Gev
    mb->SetParameters(1,0.26); // shape and temperature (260 MeV)
    
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

    Double_t topi_div_pi[41] = {2.5000002068509275e-11, 8.27088397770126e-11, 2.736301740746683e-10, 9.052655292052236e-10, 2.9949388524741494e-09, 
                                9.908318221452106e-09, 3.278022497887889e-08, 1.084485892755005e-07, 3.5878630572083736e-07, 1.1869918202234153e-06, 
                                3.926985676583339e-06, 1.2991819933971804e-05, 4.298111116507464e-05, 0.00014219203290743998, 0.0004703701265936777, 
                                0.0015555877026633969, 0.00514027405974854, 0.016938553474117934, 0.05530720928632632, 0.17514843789181533, 
                                0.49999999999999994,  
                                0.8248515621081847, 0.9446927907136737, 0.983061446525882, 0.9948597259402514, 0.9984444122973366, 
                                0.9995296298734063, 0.9998578079670926, 0.9999570188888349, 0.999987008180066, 0.9999960730143234, 
                                0.9999988130081798, 0.9999996412136942, 0.9999998915514108, 0.999999967219775, 0.9999999900916818, 
                                0.9999999970050611, 0.9999999990947345, 0.9999999997263698, 0.9999999999172912, 0.999999999975 };

    Double_t topi1[31] = {1.00000000e-05, 2.22028125e-05, 4.92964884e-05, 1.09452069e-04, 2.43014377e-04,
                        5.39560265e-04, 1.19797554e-03, 2.65984263e-03, 5.90559874e-03, 1.31120902e-02,
                        2.91125280e-02, 6.46380000e-02, 1.43514540e-01, 3.18642642e-01, 7.07476284e-01,
                        1.57079633e+00,
                        2.43411637e+00, 2.82295001e+00, 2.99807811e+00, 3.07695465e+00, 3.11248013e+00,
                        3.12848056e+00, 3.13568705e+00, 3.13893281e+00, 3.14039468e+00, 3.14105309e+00,
                        3.14134964e+00, 3.14148320e+00, 3.14154336e+00, 3.14157045e+00, 3.14158265e+00
        
    };
    
    
    
    
    
    Double_t eecbounds1[31] =   {1.00000000e-05, 2.05714387e-05, 4.23184092e-05, 8.70550563e-05, 1.79084776e-04,
                            3.68403150e-04, 7.57858283e-04, 1.55902353e-03, 3.20713570e-03, 6.59753955e-03,
                            1.35720881e-02, 2.79197379e-02, 5.74349177e-02, 1.18151889e-01, 2.43055435e-01,
                            5.00000000e-01,
                            7.56944565e-01, 8.81848111e-01, 9.42565082e-01, 9.72080262e-01, 9.86427912e-01,
                            9.93402460e-01, 9.96792864e-01, 9.98440976e-01, 9.99242142e-01, 9.99631597e-01,
                            9.99820915e-01, 9.99912945e-01, 9.99957682e-01, 9.99979429e-01, 1.00000000};
    
    Double_t topi_div_pi1[31] = {2.5000002068509275e-11, 1.232411950269352e-10, 6.07535965979622e-10, 2.9949388524741494e-09, 1.4763996747380048e-08, 
                                7.278131802790355e-08, 3.5878630572083736e-07, 1.7686896639501626e-06, 8.718998747880846e-06, 4.298111116507464e-05, 
                                0.00021186985581017614, 0.0010441541372358532, 0.00514027405974854, 0.02516923842387242, 0.11999774832124455, 
                                0.49999999999999994, 
                                0.8800022516787553, 0.9748307615761276, 0.9948597259402514, 0.9989558458627641, 0.9997881301441898, 
                                0.9999570188888349, 0.9999912810012521, 0.999998231310336, 0.9999996412136942, 0.9999999272186819, 
                                0.9999999852360033, 0.9999999970050611, 0.999999999392464, 0.9999999998767588, 0.999999999975};
                                
    Double_t topi2[61] = {1.00000000e-05, 1.49006082e-05, 2.22028125e-05, 3.30835411e-05,
                            4.92964884e-05, 7.34547660e-05, 1.09452069e-04, 1.63090240e-04,
                            2.43014377e-04, 3.62106202e-04, 5.39560265e-04, 8.03977611e-04,
                            1.19797554e-03, 1.78505642e-03, 2.65984263e-03, 3.96332730e-03,
                            5.90559873e-03, 8.79970130e-03, 1.31120901e-02, 1.95378118e-02,
                            2.91125279e-02, 4.33794373e-02, 6.46379999e-02, 9.63145513e-02,
                            1.43514539e-01, 2.13845393e-01, 3.18642641e-01, 4.74796916e-01,
                            7.07476283e-01, 1.05418269e+00, 1.57079633e+00, 
                            2.0874099625242373, 2.4341163709009024, 2.6667957376383757, 2.8229500122560345, 
                            2.927747261033487, 2.9980781141479063, 3.0452781022917392, 3.0769546536504087, 
                            3.098213216320398, 3.1124801256716323, 3.122054841779965, 3.1284805634480732, 
                            3.132792952294784, 3.1356870548603424, 3.137629326292187, 3.1389328109578147, 
                            3.139807597172383, 3.1403946780498013, 3.1407886759785417, 3.14105309332515, 
                            3.1412305473879445, 3.141349639213069, 3.1414295633499876, 3.1414832015208316, 
                            3.1415191988238047, 3.1415433571014066, 3.141559570048725, 3.1415704507772717, 
                            3.141577752981577, 3.141582653589793};
                            
                            
    Double_t eecbounds2[61] = {1e-05, 1.43814e-05, 2.06634e-05, 2.96705e-05, 4.25849e-05, 
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
    
    Double_t topi_div_pi2[61] = {2.5000002068509275e-11, 5.5507043406066714e-11, 1.232411950269352e-10, 2.736301740746683e-10, 6.07535965979622e-10, 
                                1.3489006578382146e-09, 2.9949388524741494e-09, 6.649606587583179e-09, 1.4763996747380048e-08, 3.278022497887889e-08, 
                                7.278131802790355e-08, 1.6159499116596976e-07, 3.5878630572083736e-07, 7.966063917952404e-07, 1.7686896639501626e-06, 
                                3.926985676583339e-06, 8.718998747880846e-06, 1.9358560801097102e-05, 4.298111116507464e-05, 9.54284868925348e-05, 
                                0.00021186985581017614, 0.0004703701265936777, 0.0010441541372358532, 0.002317330974731746, 0.00514027405974854, 
                                0.011388962265783398, 0.02516923842387242, 0.05530720928632632, 0.11999774832124455, 0.2530307286888087, 
                                0.49999999999999994, 
                                0.7469692713111913, 0.8800022516787553, 0.9446927907136737, 0.9748307615761276, 0.9886110377342165, 
                                0.9948597259402514, 0.9976826690252683, 0.9989558458627641, 0.9995296298734063, 0.9997881301441898, 
                                0.9999045715131074, 0.9999570188888349, 0.999980641439199, 0.9999912810012521, 0.9999960730143234, 
                                0.999998231310336, 0.9999992033936083, 0.9999996412136942, 0.9999998384050088, 0.9999999272186819, 
                                0.999999967219775, 0.9999999852360033, 0.9999999933503934, 0.9999999970050611, 0.9999999986510993, 
                                0.999999999392464, 0.9999999997263698, 0.9999999998767588, 0.999999999944493, 0.999999999975 };
    // Initialize histogram'
    
    int bins = 40;
    int bins1 = 30;
    int bins2 = 60;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    

    TH1D EEC_corr("EEC_corr", "Energy Energy Correlator", bins, topi_div_pi); 
    TH1D EEC_theta("EEC_theta", "Energy Energy Correlator", bins, topi); 

    TH1D EEC_ltt("EEC_ltt", "Energy Energy Correlator", 40, 0, M_PI); // mix signal phi
    TH1D EEC_ltz("EEC_ltz", "Energy Energy Correlator", 40, 0,  1); // mix1 mix1 phi 
    
    TH1D EEC_z("EEC_z", "Energy Energy Correlator",bins, eecbounds); // mix1 mix2 phi
    
    TH1D EEC_corr1("EEC_corr1", "Energy Energy Correlator", bins1, topi_div_pi1); 
    TH1D EEC_theta1("EEC_theta1", "Energy Energy Correlator", bins1, topi1); 
    TH1D EEC_z1("EEC_z1", "Energy Energy Correlator",bins1, eecbounds1); // mix1 mix2 phi
    
    TH1D EEC_corr2("EEC_corr2", "Energy Energy Correlator", bins2, topi_div_pi2); 
    TH1D EEC_theta2("EEC_theta2", "Energy Energy Correlator", bins2, topi2);
    TH1D EEC_z2("EEC_z2", "Energy Energy Correlator",bins2, eecbounds2); // mix1 mix2 phi
    TH1D JetSpectrum("JetSpectrum", "Jet p{T} spectrum", 70, 0, 70);
    
    // Initialize Pythia for 200 GeV pp collision
    Pythia pythia;
    
    if (!pythia.readFile(argv[1])) {
    cerr << "Error: Could not read file.cmnd!" << endl;
    return 1;
    }
    
    pythia.init();

    for (int iEvent = 0; iEvent < 1000; ++iEvent) { // 1M events
        if (!pythia.next()) continue;
    
        vector<PseudoJet> thermalParticles = thermalGen(1500);
        
        vector<PseudoJet> m1 = thermalGen(1500);
        vector<PseudoJet> m2 = thermalGen(1500);
        
        for (size_t i = 0; i < m1.size(); ++i) {
            for (size_t k = 0; k < m2.size(); ++k) {
                double mm12_eec = m1[i].pt() * m2[k].pt();
                double mm12_ctheta = costheta(m1[i], m2[k]);
                double mm12_z = (1 - mm12_ctheta) / 2;
                double mm12_delphi = deltaphi(m1[i], m2[k]);
                double theta = acos(mm12_ctheta);
                
                if (theta < 0) theta = -1 * theta;
                
                
                EEC_z.Fill(mm12_z, mm12_eec);
                EEC_corr.Fill(mm12_z, mm12_eec);
                EEC_theta.Fill(theta, mm12_eec);
                EEC_z1.Fill(mm12_z, mm12_eec);
                EEC_corr1.Fill(mm12_z, mm12_eec);
                EEC_theta1.Fill(theta, mm12_eec);
                EEC_z2.Fill(mm12_z, mm12_eec);
                EEC_corr2.Fill(mm12_z, mm12_eec);
                EEC_theta2.Fill(theta, mm12_eec);
                EEC_ltt.Fill(theta, mm12_eec);
                EEC_ltz.Fill(mm12_z, mm12_eec);
    
            }
        }
      
        /*
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
                    EEC_spl.Fill(delphi, eec/pmq2);                }
                    
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
    
        
        vector<PseudoJet> m2 = thermalGen(1500);
        
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
        */
        
    }//! event loop 
    
  
    out->Write();
  
    // Pythia cleanup
    pythia.stat();
    out->Close();
    return 0;
    }
