#include <TFile.h>
#include <TH1D.h>
#include <TH1.h>
#include <TH2.h>
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
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

using namespace std;
void threeto1() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    TFile* fin = TFile::Open("../outfiles/pythia_mixedevent_4045_15pthat50_May27_50m.root");
    
    const int number = 5;
    double min = 1e-2; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e1;
    
    TH1D* histos[number];
    TH1D* EEC_w_p[number];
    TH1D* hrebin[number];
    TH1D* hrebin2[number];
    
    //vector<const char*> histonames = {"EEC_wp", "EEC_tp", "EEC_mp", "EEC_sp"};
    //vector<const char*> histonames = {"EEC_w"};
    // these are for variable pt1 and pt2 as of April23
    //vector<const char*> histonames = {"EEC_w", "EEC_t", "EEC_m", "EEC_s"};
    
    vector<const char*> histonames = {"EEC_w", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_s"};
    
    //vector<const char*> histonames = {"EEC_w", "EEC_w_low", "EEC_w_mid", "EEC_w_high"};
    //vector<const char*> histonames = {"EEC_w_phi", "EEC_w_lowphi", "EEC_w_midphi", "EEC_w_hiphi"};
    //vector<const char*> histonames = {"EEC_d"};
    
    bool rebin = 0;
    bool delphi = 0; 
    bool mbs = 1;
    
    string outputFilename = "../plots/tryjune20.png";
    
    int bins = 60;
    
    for (int i = 0; i < number; ++i) {
        cout << "i = " << i << endl;
        TH1D* htemp = dynamic_cast<TH1D*>(fin->Get(histonames[i]));
        // to keep it separte and not mess up your root file
        histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
        histos[i]->SetDirectory(0);
        
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
        if(mbs){
            
            if(i > 0){
               
            double norm = histos[0]->Integral();
            histos[i]->Scale(norm / histos[i]->Integral());
            }
        }
     
        
        EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", 60, 0, 60);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
            double newcontent = histos[i]->GetBinContent(bin);
            //if (newcontent > 0 && newcontent < min) min = newcontent;
           // if (newcontent > max) max = newcontent;
            }
        // figuring out min and max here because I am lazy and dont want to have to put it in two conditionals. It should work 
        
        cout << min << " " << max << endl;
        if(rebin){
    
            hrebin[i] = (TH1D*)EEC_w_p[i]->Clone(Form("hrebin%d", i));
            hrebin[i]->Rebin(6);
            hrebin[i]->Scale(1./6);
            hrebin2[i] = (TH1D*)hrebin[i]->Clone(Form("hrebin2_%d", i));
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
    TH1D* EEC_empty = new TH1D("EEC_empty", "", 60, 0, 60);
    EEC_empty->SetStats(0);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("#frac{1}{W_{pairs}} #frac{1}{(#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2}} EEC");
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    
    if(rebin){
        
    // horribly done very innefficient but it works :)
        for (int i = 0; i < number; ++i) {
            hrebin[i] = (TH1D*)EEC_w_p[i]->Clone(Form("hrebin%d", i));
            hrebin[i]->Rebin(6);
            hrebin[i]->Scale(1./6);
            hrebin2[i] = (TH1D*)hrebin[i]->Clone(Form("hrebin2_%d", i));
            hrebin[i]->GetXaxis()->SetRange(1,3);
            EEC_w_p[i]->GetXaxis()->SetRange(19,42);
            hrebin2[i]->GetXaxis()->SetRange(8,10);
            
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
        }
    }// end of rebin loop 
    
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
                else if (i == 0){
                    EEC_w_p[i]->SetLineColor(kOrange + 7);
                    EEC_w_p[i]->SetMarkerColor(kOrange + 7);
                }
                else{
                    EEC_w_p[i]->SetLineColor(i + 1);
                    EEC_w_p[i]->SetMarkerColor(i + 1);
                }
                
                EEC_w_p[i]->GetYaxis()->SetLabelOffset(999);
                EEC_w_p[i]->Draw("P E1 SAME");
            
            }
        }//end of not rebin loop 
            
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
        /*
        drawText("Work in Progress", 0.15, 0.38, 16);
        drawText("PYTHIA8 p+p #sqrt{s} = 200 GeV", 0.15, 0.345, 16);
        drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.31, 16);
        drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.26, 16);
        drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.22, 16);
        //drawText("jet p_{T2} #geq 9.4 GeV", 0.15, 0.18, 16);
        //drawText("5 GeV < #hat{p}_{T} #leq 60 GeV", 0.15, 0.14, 16);
        drawText("charged constituent p_{T} > 0.2 GeV", 0.15, 0.18, 16);
        drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.14, 16);
        canvas1->Update();
        
        */
        
        TLegend* legend = new TLegend(0.65, 0.125, 0.875, 0.325);
        legend->AddEntry(EEC_w_p[0], "All contributions", "pl");
        legend->AddEntry(EEC_w_p[1], "Thermal + Thermal", "pl");
        legend->AddEntry(EEC_w_p[2], "Thermal + Signal", "pl");
        legend->AddEntry(EEC_w_p[3], "Signal + Signal", "pl");
        
        legend->Draw();
        
       /*
        //PT2 VARIABLE BINS
        TLegend* legend = new TLegend(0.60, 0.125, 0.88, 0.375);
        legend->AddEntry(EEC_w_p[0], "9.4 GeV #leq jet p_{T2} < 60.8 GeV", "pl");
        legend->AddEntry(EEC_w_p[2], "9.4 GeV #leq jet p_{T2} < 20.9 GeV", "pl");
        legend->AddEntry(EEC_w_p[1], "20.9 GeV #leq jet p_{T2} < 27.3 GeV", "pl");
        legend->AddEntry(EEC_w_p[3], "27.3 GeV #leq jet p_{T2} < 31.2 GeV", "pl");
        legend->Draw();
        
    
        
        // THIS IS FOR PT1 BINS FROM DAN'S QM POSTER 
        // Legend
        TLegend* legend = new TLegend(0.60, 0.125, 0.88, 0.375);
        legend->AddEntry(EEC_w_p[0], "20.9 GeV #leq jet p_{T1} < 60.8 GeV", "pl");
        legend->AddEntry(EEC_w_p[1], "20.9 GeV #leq jet p_{T1} < 31.2 GeV", "pl");
        legend->AddEntry(EEC_w_p[2], "31.2 GeV #leq jet p_{T1} < 40.7 GeV", "pl");
        legend->AddEntry(EEC_w_p[3], "40.7 GeV #leq jet p_{T1} < 60.8 GeV", "pl");
        legend->Draw();

        */
        canvas1->SaveAs(outputFilename.c_str());
    }

        