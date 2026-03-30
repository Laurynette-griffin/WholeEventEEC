#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <iostream>

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(8);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}



// code takes in 4 files with the histograms of the same name and plots them on 1 canvas :)

void plotmaker4to1() {
    TH1::SetDefaultSumw2();
    
    const int number = 4;
    double min = 5e-3; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e1;
    
    const char* filenames[number] = {
       "../outfiles/pythia_pp_rhic_45GeV_fullevent_44pthat45_on_April19.root",
        "../outfiles/pythia_pp_rhic_45GeV_fullevent_44pthat45_isroff_April19.root",
        "../outfiles/pythia_pp_rhic_45GeV_fullevent_44pthat45_mpioff_April19.root",
        "../outfiles/pythia_pp_rhic_45GeV_fullevent_44pthat45_bothoff_April19.root",
    };
    
    // things changed from plot to plot 
    const char* histoname = "EEC_w_b"; 
    string outputFilename = "../plots/MPIISRToggledMay30delphirebin.png";
    bool rebin = 1; // bool for rebinning with 3 bins in small and large angle for fluctuations  
    bool delphi = 1; // bool for delta phi plots for axes and stuff
    
    
    
    TFile* files[number];
    TH1F* histos[number];
    TH1F* EEC_w_p[number];
    TH1F* hrebin[number];
    TH1F* hrebin2[number];
    
    for (int i = 0; i < number; ++i) {
    files[i] = TFile::Open(filenames[i]);
    if (!files[i] || files[i]->IsZombie()) {
        std::cerr << "Error opening file: " << filenames[i] << std::endl;
        return;
    }

    TH1F* htemp = dynamic_cast<TH1F*>(files[i]->Get(histoname));
    if (!htemp) {
        std::cerr << "Histogram not found in: " << filenames[i] << std::endl;
        return;
    }

    // Clone and rename to keep each histogram separate
    histos[i] = (TH1F*)htemp->Clone(Form("EEC_clone_%d", i));
    histos[i]->SetDirectory(0);  // Detach from file
    files[i]->Close();
    
    
    // Normalize to bin width
    for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
    float width = histos[i]->GetBinWidth(bin);
    float content = histos[i]->GetBinContent(bin);
    float error = histos[i]->GetBinError(bin);
        if (width > 0) {
        histos[i]->SetBinContent(bin, content / width);
        histos[i]->SetBinError(bin, error / width);
        }
    }
    // 
    
    histos[i]->Scale(1.0/ histos[i]->Integral("width"));
    
    EEC_w_p[i] = new TH1F(Form("EEC_w_p_%d", i), "", 60, 0, 60);
    
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
        EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
        double newcontent = histos[i]->GetBinContent(bin);
        //incase min and max dont contain all the data points 
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
        }
    // figuring out min and max here because I am lazy and dont want to have to put it in two conditionals. It should work 
    
    
   if(rebin){

    hrebin[i] = (TH1F*)EEC_w_p[i]->Clone(Form("hrebin%d", i));
    hrebin[i]->Rebin(6);
    hrebin[i]->Scale(1./6);
    hrebin2[i] = (TH1F*)hrebin[i]->Clone(Form("hrebin2_%d", i));
    hrebin[i]->GetXaxis()->SetRange(1,3);
    EEC_w_p[i]->GetXaxis()->SetRange(19,42);
    hrebin2[i]->GetXaxis()->SetRange(8,10);
   }
}

TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 600);
canvas1->cd();
canvas1->SetLogy();
canvas1->SetLeftMargin(0.125); 

// empty hist since we are rebinning on the edges 
TH1F* EEC_empty = new TH1F("EEC_w_m", "", 60, 0, 60);
EEC_empty->SetStats(0);
//EEC_empty->SetMinimum(1e-3);
//EEC_empty->SetMaximum(1e3);
EEC_empty->GetXaxis()->SetLabelOffset(999);
EEC_empty->GetXaxis()->SetTitleOffset(999);
EEC_empty->GetXaxis()->SetTickLength(0);
EEC_empty->GetYaxis()->SetTitle("#frac{1}{W_{pairs}} #frac{1}{(#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2}} EEC");
EEC_empty->SetStats(0);
EEC_empty->SetMinimum(min);
EEC_empty->SetMaximum(max);
EEC_empty->Draw("P E1");



// Axis embellishments
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


axis1->Draw();
axis2->Draw(); 

canvas1->Update();

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

axis4->Draw();
axis5->Draw();


//EEC_w_p[0]->SetMarkerSize(2);
//drawText("0.5", 0.4875, .06, 17 );
//drawText("z = (1 - cos(#theta))/2", 0.8, .025, 16);
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
drawText("Work in Progress", 0.15, 0.38, 16);
drawText("PYTHIA8 p+p #sqrt{s} = 200 GeV", 0.15, 0.34, 16);
drawText("43 GeV < jet p_{T1} #leq 45 GeV", 0.15, 0.30, 16);
drawText("43 GeV < jet p_{T2} #leq 45 GeV", 0.15, 0.26, 16);
drawText("44.75 GeV < #hat{p}_{T} #leq 45.25 GeV", 0.15, 0.22, 16);
drawText("charged constituent p_{T} > 0.2 GeV", 0.15, 0.18, 16);
drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 7#pi/8", 0.15, 0.14, 16);
//drawText("Deltaphi", 0.60, .88, 16 );


if(rebin){
    
// horribly done very innefficient but it works :)
    for (int i = 0; i < number; ++i) {
        //marker size changes & style assignments 
        EEC_w_p[i]->SetStats(0);
        hrebin[i]->SetStats(0);
        hrebin2[i]->SetStats(0);
        
        //20 is filled in circle 21 us
        if (i < 2) {
            EEC_w_p[i]->SetMarkerStyle(20 + i);   
            hrebin[i]->SetMarkerStyle(20 + i);  
            hrebin2[i]->SetMarkerStyle(20 + i); 
            
            EEC_w_p[i]->SetMarkerSize(.75);
            hrebin[i]->SetMarkerSize(.75);
            hrebin2[i]->SetMarkerSize(.75);
        }
        else if (i == 2){
            EEC_w_p[i]->SetMarkerStyle(23);   
            hrebin[i]->SetMarkerStyle(23);  
            hrebin2[i]->SetMarkerStyle(23); 
            
            EEC_w_p[i]->SetMarkerSize(1);
            hrebin[i]->SetMarkerSize(1);
            hrebin2[i]->SetMarkerSize(1);
        }
        else { 
            EEC_w_p[i]->SetMarkerStyle(33);
            hrebin[i]->SetMarkerStyle(33);
            hrebin2[i]->SetMarkerStyle(33);
            
            EEC_w_p[i]->SetMarkerSize(1.25);
            hrebin[i]->SetMarkerSize(1.25);
            hrebin2[i]->SetMarkerSize(1.25);
        }
        
        //color assignments 418 is green+2, 616 is kmagenta 
        if (i == 2) {
            EEC_w_p[i]->SetLineColor(418);
            EEC_w_p[i]->SetMarkerColor(418);
            hrebin[i]->SetLineColor(418);
            hrebin[i]->SetMarkerColor(418);
            hrebin2[i]->SetLineColor(418);
            hrebin2[i]->SetMarkerColor(418);
        }
        
        else if(i == 1) {
             EEC_w_p[i]->SetLineColor(616);
            EEC_w_p[i]->SetMarkerColor(616);
            hrebin[i]->SetLineColor(616);
            hrebin[i]->SetMarkerColor(616);
            hrebin2[i]->SetLineColor(616);
            hrebin2[i]->SetMarkerColor(616);
        }

        else { // 0 and 3
            EEC_w_p[i]->SetLineColor(i + 1);
            EEC_w_p[i]->SetMarkerColor(i + 1);
            hrebin[i]->SetLineColor(i + 1);
            hrebin[i]->SetMarkerColor(i + 1);
            hrebin2[i]->SetLineColor(i + 1);
            hrebin2[i]->SetMarkerColor(i + 1);
            }
        
        hrebin[i]->Draw("P E1 SAME");
        EEC_w_p[i]->Draw("P E1 SAME");
        hrebin2[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        
        TLegend* legend = new TLegend(0.65, 0.125, 0.85, 0.325);
        if (i ==0){
        legend->AddEntry(EEC_w_p[0], "MPI on ISR on", "pl");
        }
        if (i ==1){
        legend->AddEntry(EEC_w_p[0], "MPI on ISR on", "pl");
        legend->AddEntry(EEC_w_p[1], "MPI on ISR off", "pl");
        }
        if (i ==2){
        legend->AddEntry(EEC_w_p[0], "MPI on ISR on", "pl");
        legend->AddEntry(EEC_w_p[1], "MPI on ISR off", "pl");
        legend->AddEntry(EEC_w_p[2], "MPI off ISR on", "pl");
        }
        //legend->AddEntry(EEC_w_p[2], "MPI on ISR off", "pl");
        if (i ==3){
        legend->AddEntry(EEC_w_p[0], "MPI on ISR on", "pl");
        legend->AddEntry(EEC_w_p[1], "MPI on ISR off", "pl");
        legend->AddEntry(EEC_w_p[2], "MPI off ISR on", "pl");
        legend->AddEntry(EEC_w_p[3], "MPI off ISR off", "pl");
        }
        legend->Draw();
        TString filename = Form("MPIrebin_%d.png", i);
        canvas1->SaveAs(filename);
    }
}

if(!rebin){
    for (int i = 0; i < number; ++i) {
        EEC_w_p[i]->SetStats(0);
        
        if (i < 2) { // smaller markers for the square and circle
            EEC_w_p[i]->SetMarkerStyle(20 + i);  
            EEC_w_p[i]->SetMarkerSize(.75);
        }
        else if (i == 2) {// bigger markers 
            EEC_w_p[i]->SetMarkerStyle(23);
            EEC_w_p[i]->SetMarkerSize(1);
        }
        else {
            EEC_w_p[i]->SetMarkerStyle(33);
            EEC_w_p[i]->SetMarkerSize(1.25);
        }
        if (i == 2) {
        EEC_w_p[i]->SetLineColor(418);
        EEC_w_p[i]->SetMarkerColor(418);
        }
        else if(i == 1) {
        EEC_w_p[i]->SetLineColor(616);
        EEC_w_p[i]->SetMarkerColor(616);
        }
        else{
        EEC_w_p[i]->SetLineColor(i + 1);
        EEC_w_p[i]->SetMarkerColor(i + 1);
        }
        EEC_w_p[i]->GetYaxis()->SetLabelOffset(999);
        EEC_w_p[i]->Draw("P E1 SAME");

    }
}//end of else loop 

canvas1->Update();
// Legend
TLegend* legend = new TLegend(0.65, 0.125, 0.85, 0.325);
legend->AddEntry(EEC_w_p[0], "MPI on ISR on", "pl");
legend->AddEntry(EEC_w_p[2], "MPI on ISR off", "pl");
legend->AddEntry(EEC_w_p[1], "MPI off ISR on", "pl");
//legend->AddEntry(EEC_w_p[2], "MPI on ISR off", "pl");
legend->AddEntry(EEC_w_p[3], "MPI off ISR off", "pl");
legend->Draw();
canvas1->Update();

canvas1->SaveAs(outputFilename.c_str());
}
    