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
void MixedBackgroundSub() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
   /* const char* filename = {
       "../outfiles/pythia_mixedevent_03_15pthat50_May27_50m.root"};
       //"../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_May22_50m.root"};
       */
    TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_mixedeventb_50M_Jul31.root");

    
    int const number = 5;
    double min = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e1;
    
    double min1 = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max1 = 1e1;
          
    
    
    //mixed event 
    //vector<const char*> histonames = { "EEC_sp"};
    vector<const char*> histonames = {"EEC_w", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_s"};
    
    bool rebin = 0;
    bool delphi = 0;
    
    TH1F* histos[number];
    TH1F* prehistos[number];
    TH1F* preEEC_w_p[number];
    TH1F* EEC_w_p[number];
    TH1F* mbssignal;
    TH1F* mbssignal1final;
    TH1F* mbssignal1;
    TH1F* mbssignalfinal;
    
    string outputFilename = "../plots/mbscentrality03Aug4.C";
    string outputFilename1 = "../plots/mbscentrality03Aug4.png";
    
    int bins = 60;
    for (int i = 0; i < number; ++i) {
    TH1F* htemp = dynamic_cast<TH1F*>(file->Get(histonames[i]));
            // to keep it separte and not mess up your root file
            histos[i] = (TH1F*)htemp->Clone(Form("EEC_clone_%d", i));
            histos[i]->SetDirectory(0);
            
            // histos that will be subtracted 
            prehistos[i] = (TH1F*)htemp->Clone(Form("EEC_clone_%d", i));
            prehistos[i]->SetDirectory(0);
           
            // gonna do the subtraction before and after bin width and see what the heck goes on. 
            // Subtraction before normalization but idk which one fr
            
            // Normalize to bin width
            for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
                double width = histos[i]->GetBinWidth(bin);
                double content = histos[i]->GetBinContent(bin);
                double error = histos[i]->GetBinError(bin);
                histos[i]->SetBinContent(bin, content / width);
                histos[i]->SetBinError(bin, error / width);
            }
            //mbs stuff only
    }
   /* cout <<"seg check 1" << endl;
    // post bin width subtraction :)
    mbssignal1 =(TH1F*)histos[0]->Clone();
    mbssignal1->Add(histos[1], -1);
    mbssignal1->Add(histos[2], -1);
    mbssignal1->Add(histos[3],  +1);
    cout <<"seg check 2" << endl;
    //scale by integram
    mbssignal1->Scale(1.0/ mbssignal1->Integral("width"));
    
    //log to linear
    mbssignal1final = new TH1F(Form("mbssignal1final"), "", 60, 0, 60);
    cout <<"seg check 3" << endl;    
    for (int bin = 1; bin <= 60; ++bin) {
        mbssignal1final->SetBinContent(bin, mbssignal1->GetBinContent(bin));
        mbssignal1final->SetBinError(bin, mbssignal1->GetBinError(bin));
        double newcontent1 = mbssignal1->GetBinContent(bin);
        //cout << newcontent << " " << bin << endl;
        if (newcontent1 > 0 && newcontent1 < min1) min1 = newcontent1;
        if (newcontent1 > max1) max1 = newcontent1;
    }
    cout <<"seg check 4" << endl;
    */
    
    for (int i = 0; i < number; ++i) {        
        // regular histos[i] loop
        histos[i]->Scale(1.0/ histos[i]->Integral("width"));
        
        EEC_w_p[i] = new TH1F(Form("EEC_w_p_%d", i), "", 60, 0, 60);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
            double newcontent = histos[i]->GetBinContent(bin);
            //cout << newcontent << " " << bin << endl;
            if (newcontent > 0 && newcontent < min) min = newcontent;
            if (newcontent > max) max = newcontent;
        }
    cout <<"seg check 5" << endl;
    } //closes i = loop

//pre bin width SUBTRACTION IS GONNA OCCUR NOW 
//Subhistos are subtracted then normalized by bin width then integral 
// okay so EEC_w - EEC_ms - EEC_mm + EEC_m2 = EEC_s allegedly and thats what this code does :)
    mbssignal = (TH1F*)prehistos[0]->Clone();
    mbssignal->Add(prehistos[1], -1);
    mbssignal->Add(prehistos[2], -1);
    mbssignal->Add(prehistos[3], +1);
    
    //mbs has occured now to normaliza just mbssignal by bin width
    for (int bin = 1; bin <= mbssignal->GetNbinsX(); ++bin) {
            double width = mbssignal->GetBinWidth(bin);
            double content = mbssignal->GetBinContent(bin);
            double error = mbssignal->GetBinError(bin);
                mbssignal->SetBinContent(bin, content / width);
                mbssignal->SetBinError(bin, error / width);
    }
    
    mbssignal->Scale(1.0/mbssignal->Integral("width"));
    
    // double log to linear for mbs 
    //reinitiallizing min and max 
    
    mbssignalfinal = new TH1F(Form("mbssignalfinal"), "", 60, 0, 60);
    for (int bin = 1; bin <= 60; ++bin) {
        mbssignalfinal->SetBinContent(bin, mbssignal->GetBinContent(bin));
        mbssignalfinal->SetBinError(bin, mbssignal->GetBinError(bin));
        double newcontent = mbssignal->GetBinContent(bin);
        //cout << newcontent << " " << bin << endl;
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
 // pre bin width subtraction + bin width & integral normalization + log to linear shift (final histogram )  :)
    
    //plot :)
    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 600);
    canvas1->cd();
    canvas1->SetLogy();
    canvas1->SetLeftMargin(0.125); 
    cout << "seg check 1" << endl;
    // empty hist since we are rebinning on the edges 
    TH1F* EEC_empty = new TH1F("EEC_empty", "", 60, 0, 60);
    EEC_empty->SetStats(0);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1 / (#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2} EEC");
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    for (int i = 0; i < number; ++i) {
        //marker size changes & style assignments 
        cout << i << endl;
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
        if (i == 0 || i == 4){
        EEC_w_p[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }
    }
    /*
	mbssignal1final->SetMarkerStyle(41);
	mbssignal1final->SetMarkerColor(kOrange+1);
	mbssignal1final->SetLineColor(kOrange+1);
    mbssignal1final->Draw("P E1 SAME");
    */
    mbssignalfinal->SetMarkerStyle(20);
	mbssignalfinal->SetMarkerColor(kCyan+3);
	mbssignalfinal->SetLineColor(kCyan+3);
	mbssignalfinal->Draw("P E1 SAME");
    
    TGaxis *axis1 = new TGaxis(0, 0, 30, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(60, 0, 30, 0, .00001 , .5 , 510, "-G");
    axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    if(delphi){
        axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-4}" );
        axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-2}" );
    }
    else{
        axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
        axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
    }
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->SetLabelOffset(0.045);
    axis1->SetLabelSize(.03);
    axis2->SetLabelSize(.03);
    axis1->SetTickSize(0.05);
    axis2->SetTickSize(0.05);
    
 
    TGaxis *axis3 = new TGaxis( 60, min , 60, max , min, max, 510,"+LG");
    axis3->SetLabelFont(42);
    axis3->SetLabelSize(.035);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(0, max, 30, max, .00001, .5, 510, "-G"); 
    TGaxis *axis5 = new TGaxis(60, max, 30, max, .00001 , .5 , 510, "+G");
    axis4->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis4->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis4->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis5->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis5->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    if(delphi){
        axis5->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-4}" );
        axis5->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-2}" );
        }
    else{
        axis5->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
        axis5->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
        }
    axis5->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis4->SetLabelOffset(0.015);
    axis5->SetLabelOffset(-0.025);
    axis5->SetLabelFont(42);
    axis4->SetLabelSize(.03);
    axis5->SetLabelSize(.03);
    axis4->SetTickSize(0.05);
    axis5->SetTickSize(0.05);
    
    axis1->Draw();
    axis2->Draw();
    axis3->Draw();
    axis4->Draw();
    axis5->Draw();
    
    if(delphi){
        drawText("#pi/2", 0.505, .065, 17 );
        drawText("#pi/2", 0.505, .91, 17 );
        drawText("#Delta#phi", 0.8, .025, 16);
    }
    else{
        drawText("0.5", 0.505, .065, 17 );
        drawText("0.5", 0.50, .91, 17 );
        drawText("z = #frac{(1 - cos(#theta))}{2}", 0.8, .025, 16);
    }
    if(delphi){
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.38, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.34, 15);
    drawText("15 GeV < #hat{p}_{T} #leq 50 GeV", 0.15, 0.30, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.26, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.22, 15);
    //drawText("jet p_{T2} #leq 9.4 GeV", 0.15, 0.26, 15);
    drawText("charged constituent p_{T} > 0.2 GeV", 0.15, 0.18, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.14, 15);
    canvas1->Update();
    }
    else{
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.35, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.32, 15);
    drawText("15 GeV < #hat{p}_{T} #leq 50 GeV", 0.15, 0.28, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.24, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.20, 15);
    //drawText("jet p_{T2} #leq 9.4 GeV", 0.15, 0.26, 15);
    drawText("charged constituent p_{T} > 0.2 GeV", 0.15, 0.16, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.12, 15);
    canvas1->Update();
    }
    
    TLegend* legend = new TLegend(0.55, 0.125, 0.875, 0.325);
    legend->AddEntry(EEC_w_p[0], "wEEC", "pl");
    //legend->AddEntry(EEC_w_p[1], "Thermal + Thermal", "pl");
    //legend->AddEntry(EEC_w_p[2], "Signal + Mixed", "pl");
    //legend->AddEntry(EEC_w_p[3], "Mixed + Mixed", "pl");
    legend->AddEntry(EEC_w_p[4], "Signal + Signal", "pl");
   // legend->AddEntry(mbssignal1final, "Mixed background subtraction post bin", "pl");
    legend->AddEntry(mbssignalfinal, "Mixed background subtraction pre bin", "pl");
    
    legend->Draw();
 
    canvas1->SaveAs(outputFilename.c_str());
    canvas1->SaveAs(outputFilename1.c_str());
    
}

    
    
    
        