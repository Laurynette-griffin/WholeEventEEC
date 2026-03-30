#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <iostream>

void drawText(const char *text, double xp, double yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(8);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

using namespace std;
void Linearweec() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_mixedeventb_phi_Dec15_5M.root");

 
    int const number = 1;
    double min = 1e-2; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e1;
    
    double min1 = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max1 = 1e1;
          
    
    vector<const char*> histonames = {"EEC_spl"};
    
    bool rebin = 0;
    bool delphi = 0;
    
    TH1D* histos[number];
    TH1D* prehistos[number];
    TH1D* preEEC_w_p[number];
    TH1D* EEC_w_p[number];
    TH1D* mbssignal;
    TH1D* mbssignalfinal;
    TH1D* mbssignal1;
    
    string outputFilename = "../plots/delphilineardec15.C";
    string outputFilename1 = "../plots/delphilineardec15.png";
    const char* pttrackcut = "charged constituent p_{T} > 0.2 GeV";
    
    int bins = 40;
    for (int i = 0; i < number; ++i) {
    TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
        // to keep it separte and not mess up your root file
        prehistos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
        prehistos[i]->SetDirectory(0);
        
        histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
        histos[i]->SetDirectory(0);
        
        histos[i]->Scale(1.0/ histos[i]->Integral("width"));
        
        EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", bins, 0, M_PI);
        
        for (int bin = 1; bin <= bins; ++bin) {
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
            double newcontent = histos[i]->GetBinContent(bin);
            //cout << newcontent << " " << bin << endl;
            //if (newcontent > 0 && newcontent < min) min = newcontent;
            //if (newcontent > max) max = newcontent;
        }
    }
    /*        
    mbssignal = (TH1D*)prehistos[0]->Clone();
    mbssignal->Add(prehistos[1], -1);
    mbssignal->Add(prehistos[2], -1);
    mbssignal->Add(prehistos[3], +1);
    
    for (int bin = 1; bin <= mbssignal->GetNbinsX(); ++bin) {
            double width = mbssignal->GetBinWidth(bin);
            double content = mbssignal->GetBinContent(bin);
            double error = mbssignal->GetBinError(bin);
            mbssignal->SetBinContent(bin, content / width);
            mbssignal->SetBinError(bin, error);
    }
    // double log to linear for mbs 
    
    mbssignalfinal = new TH1D(Form("mbssignalfinal"), "", bins, 0, 1);
    for (int bin = 1; bin <= bins; ++bin) {
        mbssignalfinal->SetBinContent(bin, mbssignal->GetBinContent(bin));
        mbssignalfinal->SetBinError(bin, mbssignal->GetBinError(bin));
        double newcontent = mbssignal->GetBinContent(bin);
        //cout << newcontent << " " << bin << endl;
        //if (newcontent > 0 && newcontent < min) min = newcontent;
        //if (newcontent > max) max = newcontent;
    }
    //mbssignal->Scale(1.0/mbssignal->Integral("width"));
    */
    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 600);
    canvas1->cd();
    canvas1->SetLogy();
    canvas1->SetLogx();
    canvas1->SetLeftMargin(0.125); 
    
    // empty hist since we are rebinning on the edges 
    TH1D* EEC_empty = new TH1D("EEC_empty", "", bins, 0, M_PI);
    EEC_empty->SetStats(0);
    //EEC_empty->GetXaxis()->SetLabelOffset(999);
    //EEC_empty->GetXaxis()->SetTitleOffset(999);
    //EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1 / (#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2} EEC");
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    for (int i = 0; i < number; ++i) {
        //marker size changes & style assignments 

        EEC_w_p[i]->SetStats(0);
        //20 is filled in circle 21 us
        if (i < 2) {
            EEC_w_p[i]->SetMarkerStyle(20 + i);   
    
            EEC_w_p[i]->SetMarkerSize(.75);
        }
        else if (i == 2){
            EEC_w_p[i]->SetMarkerStyle(23);   
            EEC_w_p[i]->SetMarkerSize(1);
    
        }
        else { 
            EEC_w_p[i]->SetMarkerStyle(33);
            EEC_w_p[i]->SetMarkerSize(1.25);
        }
        
        //color assignments 418 is green+2, 616 is kmagenta 
        if (i == 2) {
            EEC_w_p[i]->SetLineColor(418);
            EEC_w_p[i]->SetMarkerColor(418);
        }
        
        else if(i == 1) {
            EEC_w_p[i]->SetLineColor(616);
            EEC_w_p[i]->SetMarkerColor(616);
        }
        
        else if(i == 4) {
            EEC_w_p[i]->SetLineColor(kViolet+7);
            EEC_w_p[i]->SetMarkerColor(kViolet+7);
        }
    
        else { // 0 and 3
            EEC_w_p[i]->SetLineColor(i + 1);
            EEC_w_p[i]->SetMarkerColor(i + 1);
            }
        // just plot whole & signal 
        if (i < number){
        EEC_w_p[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }
    }
    //mbssignalfinal->Draw("P E1 SAME");
    
    drawText("#Delta#phi", 0.90, .03, 20);
    drawText("PYTHIA8 ", .15, .36, 17);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.32, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.285, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.24, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.20, 15);
    drawText(pttrackcut, 0.15, 0.16, 15);

    
    //legend->AddEntry(EEC_w_p[4], "Signal + Signal", "pl");
   // legend->AddEntry(mbssignal1final, "Mixed background subtraction post bin", "pl");
    //legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl");
    //legend->Draw();

    //canvas1->SaveAs();
    canvas1->SaveAs(outputFilename1.c_str());

}

    
    
    
    
    
    
    
    
    
    
            
            
