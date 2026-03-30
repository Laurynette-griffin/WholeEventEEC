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
#include <TStyle.h>
#include <TLine.h>
#include <iostream>
#include <vector>

using namespace std;

void drawText(const char *text, double xp, double yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(8);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //! tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

void mbsplotter() {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetPadTickY(1);

    //! Open the input file
    //TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy010_psi2match_trackptcut1gev_copycopy_Mar17_5k.root");
    TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_mixedeventb_phi_1e-3_Dec16_5M.root");
    if (!file || file->IsZombie()) {
        cout << "Error opening file!" << endl;
        return;
    }

    //! Configuration variables
    int const number = 7;
    int bins = 40;
    
    //! Initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double min = 1e-1; 
    double max = 1e8;
    
    double min1 = 1e-1;
    double max1 = 1e1;
    
    double minr = 0.0;
    double maxr = 2.0;

    bool rebin = 0;
    bool delphi = 1;
    bool integrate = 0;
    bool linear = 0;
    bool scale = 0;
    bool all = 0;
    bool debug = 0;
    bool ratio = 1;
    
    vector<const char*> histonames;
    if (delphi && linear) {
        histonames = {"EEC_wpl", "EEC_spl", "EEC_mpl", "EEC_tpl", "EEC_mspl", "EEC_mmpl", "EEC_m2pl"};
    } else if (delphi) {
        histonames = {"EEC_wp", "EEC_sp", "EEC_mp", "EEC_tp", "EEC_msp", "EEC_mmp", "EEC_m2p"};
    } else {
        histonames = {"EEC_w", "EEC_s", "EEC_m", "EEC_t", "EEC_ms", "EEC_mm", "EEC_m2"};
    }

    TH1D* histos[number];
    TH1D* prehistos[number];
    TH1D* EEC_w_p[number];
    TH1D* mbssignal;
    TH1D* mbssignalfinal;
    TH1D* rh;
    TH1D* sm;
    TH1D* mbsstep1;
    TH1D* mbsstep2;
    TH1D* mbsstep1final;
    TH1D* mbsstep2final;
    string baseName = "../plots/weec_mbtherm_regenMar18_phi";
    
    string outputFilename  = baseName + ".C";
    string outputFilename1 = baseName + ".png";
    string outputFiledebug = baseName + ".root";
    
    const char* pttrackcut = "charged constituent p_{T} > 0.2 GeV";

    //! Histogram extraction and processing
    if (linear) {
        for (int i = 0; i < number; ++i) {
            TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
            if (!htemp) {
                cout << "Histogram not found: " << histonames[i] << endl;
                return;
            }
            prehistos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_p_%d", i));
            prehistos[i]->SetDirectory(0);
            
            histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            histos[i]->SetDirectory(0);
            
            EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", bins, 0, bins);
        
            for (int bin = 1; bin <= bins; ++bin) {
                EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
                EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
                double newcontent = histos[i]->GetBinContent(bin);
                if (newcontent > 0 && newcontent < min) min = newcontent;
                if (newcontent > max) max = newcontent;
            }
        }
    } else {    
        for (int i = 0; i < number; ++i) {
            TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
            
            //! To keep it separate and not mess up your root file
            histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            histos[i]->SetDirectory(0);
            
            //! Histos that will be subtracted 
            prehistos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            prehistos[i]->SetDirectory(0);
            
            //! Normalize to bin width
            for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
                double width = histos[i]->GetBinWidth(bin);
                double content = histos[i]->GetBinContent(bin);
                double error = histos[i]->GetBinError(bin);
                histos[i]->SetBinContent(bin, content / width);
                
                if (i == 0 || i == 4 || i == 5 || i == 6) {
                    cout << "Bin content of prehistos "<< i << " bin number = " << bin << " is = to "<< prehistos[i]->GetBinContent(bin)  << endl;
                }
            }
        }
        for (int i = 0; i < number; ++i) {        
            if (integrate) {
                histos[i]->Scale(1.0/ histos[i]->Integral("width"));
            }
    
            EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", bins, 0, bins);
            
            for (int bin = 1; bin <= bins; ++bin) {
                EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
                EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
                double newcontent = histos[i]->GetBinContent(bin);
                cout << "Bin content of histos" << i << " "<< bin << " is = to "<< newcontent  << endl;
                if (newcontent > 0 && newcontent < min) min = newcontent;
                if (newcontent > max) max = newcontent;
            }
        }
    }

    //! Pre bin width SUBTRACTION IS GONNA OCCUR NOW 
    //! Subhistos are subtracted then normalized by bin width then integral 
    //! EEC_w - EEC_ms - EEC_mm + EEC_m2 = EEC_s allegedly and thats what this code does :)
    
    mbssignal = (TH1D*)prehistos[0]->Clone();
    mbssignal->Add(prehistos[4], -1.0);
    mbsstep1 = (TH1D*)mbssignal->Clone();
    mbssignal->Add(prehistos[5], -1.0);
    mbsstep2 = (TH1D*)mbssignal->Clone();
    mbssignal->Add(prehistos[6], +1.0); //! Opposite sign of incorrect background modeling in sm1 so this IS the right sign :) (and the only way it works)

    //! MBS has occured now to normaliza just mbssignal by bin width
    mbssignalfinal = new TH1D(Form("mbssignalfinal"), "", bins, 0, bins);
    mbsstep1final = new TH1D(Form("mbsstep1final"), "", bins, 0, bins);
    mbsstep2final = new TH1D(Form("mbsstep2final"), "", bins, 0, bins);
    
    for (int bin = 1; bin <= mbssignal->GetNbinsX(); ++bin) {
        double width = mbssignal->GetBinWidth(bin);
        double content = mbssignal->GetBinContent(bin);
        double error = mbssignal->GetBinError(bin);
        mbssignalfinal->SetBinContent(bin, content / width);
        mbssignalfinal->SetBinError(bin, error/width);
        double newcontent = mbssignalfinal->GetBinContent(bin);
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
    if (integrate) {
        mbssignalfinal->Scale(1.0/mbssignalfinal->Integral("width"));
    }
    
    for (int bin = 1; bin <= mbsstep1->GetNbinsX(); ++bin) {
        double width = mbsstep1->GetBinWidth(bin);
        double content = mbsstep1->GetBinContent(bin);
        double error = mbsstep1->GetBinError(bin);
        mbsstep1final->SetBinContent(bin, content / width);
        mbsstep1final->SetBinError(bin, error/width);
        double newcontent = mbsstep1final->GetBinContent(bin);
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
    if (integrate) {
        mbsstep1final->Scale(1.0/mbsstep1final->Integral("width"));
    }
    
    //! MBS step 2
    for (int bin = 1; bin <= mbsstep2->GetNbinsX(); ++bin) {
        double width = mbsstep2->GetBinWidth(bin);
        double content = mbsstep2->GetBinContent(bin);
        double error = mbsstep2->GetBinError(bin);
        mbsstep2final->SetBinContent(bin, content / width);
        mbsstep2final->SetBinError(bin, error/width);
        double newcontent = mbsstep2final->GetBinContent(bin);
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
    if (integrate) {
        mbsstep2final->Scale(1.0/mbsstep2final->Integral("width"));
    }

    TFile *out = new TFile(outputFiledebug.c_str(), "RECREATE");
    if (!out || out->IsZombie()) {
        cerr << "Error: Could not open output file:" << outputFiledebug.c_str() << endl;
        return;
    }
    out->cd();
    mbssignalfinal->Write();
    mbsstep1final->Write();
    mbsstep2final->Write();
    out->Close();
    
    //! Plotting
    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 800);
    canvas1->Divide(1, 2); //! Split into 2 vertical pads
    canvas1->SetLeftMargin(0.125); 
    
    TPad* pad1 = (TPad*)canvas1->cd(1);
    pad1->SetPad(0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0.0);
    pad1->SetLogy();
    
    TPad* pad2 = (TPad*)canvas1->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.30); 
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    
    pad1->cd();
    
    //! Empty hist since we are rebinning on the edges 
    TH1D* EEC_empty = new TH1D("EEC_empty", "", bins, 0, bins);
    EEC_empty->SetStats(0);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1 / (#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2} EEC");
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    for (int i = 0; i < 7; ++i) {
        EEC_w_p[i]->SetStats(0);
        //! Marker size changes & style assignments 
        //! 20 is filled in circle 21 us
        if (i < 2) {
            EEC_w_p[i]->SetMarkerStyle(20 + i);   
            if (i == 1) EEC_w_p[i]->SetMarkerSize(1.25);
        } else if (i == 2 ) {
            EEC_w_p[i]->SetMarkerStyle(23);   
            EEC_w_p[i]->SetMarkerSize(1);
        } else { 
            EEC_w_p[i]->SetMarkerStyle(33);
            EEC_w_p[i]->SetMarkerSize(1.25);
        }
        
        //! Color assignments 418 is green+2, 616 is kmagenta 
        if (i == 2) {
            EEC_w_p[i]->SetLineColor(418);
            EEC_w_p[i]->SetMarkerColor(418);
        } else if(i == 1) {
            EEC_w_p[i]->SetLineColor(kOrange+1);
            EEC_w_p[i]->SetMarkerColor(kOrange+1);
        } else if(i == 4) {
            EEC_w_p[i]->SetLineColor(kYellow+1);
            EEC_w_p[i]->SetMarkerColor(kYellow+1);
        } else if(i == 5) {
            EEC_w_p[i]->SetLineColor(46);
            EEC_w_p[i]->SetMarkerColor(46);
        } else if(i == 6) {
            EEC_w_p[i]->SetLineColor(kViolet-1);
            EEC_w_p[i]->SetMarkerColor(kViolet-1);
        } else { //! 0 and 3
            EEC_w_p[i]->SetLineColor(i + 1);
            EEC_w_p[i]->SetMarkerColor(i + 1);
        }

        if (all || i <= 1) {
            EEC_w_p[i]->Draw("P E1 SAME");
            EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
            EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
            EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }
        
        
    }
    
    mbssignalfinal->SetMarkerStyle(107);
    mbssignalfinal->SetMarkerColor(618);
    mbssignalfinal->SetLineColor(618);
    mbssignalfinal->Draw("P E1 SAME");
    
    if (debug) {
        mbsstep1final->SetMarkerStyle(108);
        mbsstep2final->SetMarkerStyle(109);
        mbsstep1final->SetMarkerColor(860);
        mbsstep2final->SetMarkerColor(880);
        mbsstep1final->SetLineColor(860);
        mbsstep2final->SetLineColor(880);
        mbsstep1final->Draw("P E1 SAME");
        mbsstep2final->Draw("P E1 SAME");
    }
    
    pad2->cd();
    rh = (TH1D*)EEC_w_p[1]->Clone();
    sm = (TH1D*)mbssignalfinal->Clone();
    rh->Divide(sm);
    
    rh->SetTitle("");
    rh->GetYaxis()->SetTitle("Pythia / Bkd Sub");
    rh->GetYaxis()->SetLabelSize(.05);
    rh->SetMarkerStyle(20);
    rh->SetMarkerSize(0.75);
    rh->SetStats(0);
    rh->SetLineColor(kBlack);
    rh->SetMarkerColor(kBlack);
    rh->GetXaxis()->SetLabelOffset(999);
    rh->GetXaxis()->SetTitleOffset(999);
    rh->GetXaxis()->SetTickLength(0);
    
    rh->SetMinimum(minr);
    rh->SetMaximum(maxr);
    
    rh->Draw("P E1");
    TLine *line = new TLine(0, 1, bins, 1);
    
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kOrange +2);
    line->Draw();
    
    TGaxis *axis7 = new TGaxis(0, minr, bins/2, minr, .0001, .5, 510, "G"); 
    TGaxis *axis8 = new TGaxis(bins, minr, bins/2, minr, .0001 , .5 , 510, "-G");
    axis7->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    if (delphi) {
        axis8->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-4}" );
        axis8->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-2}" );
    } else {
        axis8->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
        axis8->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
    }
    axis7->ChangeLabel(2, -1, -1, -1, -1, 62, "10^{-4}" );
    axis7->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis7->ChangeLabel(4, -1, -1, -1, -1, 62, "10^{-2}" );
    axis7->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis8->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis8->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis8->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis8->SetLabelOffset(0.125);
    axis7->SetLabelSize(.075);
    axis8->SetLabelSize(.075);
    axis7->SetTickSize(0.5);
    axis8->SetTickSize(0.5);
    
    axis7->Draw();
    axis8->Draw();
    
    if (delphi) {
        drawText("#pi/2", 0.4850, .21, 17 );
        drawText("#Delta#phi", 0.8, .025, 16);
    } else {
        drawText("0.5", 0.4850, .21, 16);
        drawText("z = #frac{(1 - cos(#theta))}{2}", .8, .02, 16);
    }
    
    pad1->cd();
    /*
    TGaxis *axis1 = new TGaxis(0, 0, bins/2, 0, .0001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(bins, 0, bins/2, 0, .0001 , .5 , 510, "-G");
    axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(2, -1, 0, -1, -1, 62, "10^{-4}" );
    axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(4, -1, 0, -1, -1, 62, "10^{-2}" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    
    if (delphi) {
        axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-5}" );
        axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-3}" );
    } else {
        axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
        axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
    }
    axis2->SetLabelOffset(0.045);
    axis1->SetLabelSize(.03);
    axis2->SetLabelSize(.03);
    axis1->SetTickSize(0.05);
    axis2->SetTickSize(0.05);
    
    TGaxis *axis3 = new TGaxis( bins, min , bins, max , min, max, 510,"+LG");
    axis3->SetLabelFont(42);
    axis3->SetLabelSize(.035);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(0, max, bins/2, max, .0001, .5, 510, "-G"); 
    TGaxis *axis5 = new TGaxis(bins, max, bins/2, max, .0001 , .5 , 510, "+G");
    axis4->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis4->ChangeLabel(2, -1, -1, -1, -1, 62, "10^{-4}" );
    axis4->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis4->ChangeLabel(4, -1, -1, -1, -1, 62, "10^{-2}" );
    axis4->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis5->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis5->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    
    if (delphi) {
        axis5->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-4}" );
        axis5->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-2}" );
    } else {
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
*/
    TGaxis *axis1 = new TGaxis(0, 0, bins/2, 0, .001, 1.57079633, 510, "G"); 
    TGaxis *axis2 = new TGaxis(bins, 0, bins/2, 0, .001 , 1.57079633 , 510, "-G");
    
    axis1->ChangeLabel(2, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(1, -1, -1, -1, -1, 62, "10^{-3}" );
    axis1->ChangeLabel(4, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(3, -1, -1, -1, -1, 62, "10^{-1}" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    //axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(2, -1, 0, -1, -1, 62, "-" );
   // axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(4, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    if(delphi){
        axis2->ChangeLabel(1, -1, -1, -1, -1, 62, "#pi - 10^{-3}" );
        axis2->ChangeLabel(3, -1, -1, -1, -1, 62, "#pi - 10^{-1}" );
    }
    else{
        axis2->ChangeLabel(1, -1, -1, -1, -1, 62, "1 - 10^{-3}" );
        axis2->ChangeLabel(3, -1, -1, -1, -1, 62, "1 - 10^{-1}" );
    }
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "" );
    
    axis2->SetLabelOffset(0.055);
    axis1->SetLabelOffset(0.015);
    axis1->SetLabelSize(.03);
    axis2->SetLabelSize(.03);
    axis1->SetTickSize(0.05);
    axis2->SetTickSize(0.05);
    
 
    TGaxis *axis3 = new TGaxis( bins, min , bins, max , min, max, 510,"+LG");
    axis3->SetLabelFont(42);
    axis3->SetLabelOffset(.005);
    axis3->SetLabelSize(.035);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(0, max, bins/2, max, .001, 1.5707933, 510, "-G"); 
    TGaxis *axis5 = new TGaxis(bins, max, bins/2, max, .001 , 1.5707933 , 510, "+G");
    
    axis4->ChangeLabel(2, -1, 0, -1, -1, 62, "" );
    axis4->ChangeLabel(3, -1, -1, -1, -1, 62, "10^{-1}" );
    axis4->ChangeLabel(4, -1, 0, -1, -1, 62, "" );
    axis4->ChangeLabel(1, -1, -1, -1, -1, 62, "10^{-3}" );
    axis4->ChangeLabel(5, -1, 0, -1, -1, 62, "" );
    axis5->ChangeLabel(4, -1, 0, -1, -1, 62, "" );
    axis5->ChangeLabel(2, -1, 0, -1, -1, 62, "" );
    if(delphi){
        axis5->ChangeLabel(1, -1, -1, -1, -1, 62, "#pi - 10^{-3}" );
        axis5->ChangeLabel(3, -1, -1, -1, -1, 62, "#pi - 10^{-1}" );
        }
    else{
        axis5->ChangeLabel(1, -1, -1, -1, -1, 62, "#pi - 10^{-1}" );
        axis5->ChangeLabel(3, -1, -1, -1, -1, 62, "#pi - 10^{-3}" );
        }
    axis5->ChangeLabel(5, -1, 0, -1, -1, 62, "" );
    
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
    
    
    TGaxis *axis9 = new TGaxis(0, 0, bins, 0, .001, M_PI , 2, "C"); 
    TGaxis *axis0 = new TGaxis(0, max, bins, max, .001, M_PI , 2, "-C");
    
    axis9->ChangeLabel(1, -1, 0, -1, -1, 62, "" );
    axis9->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi/2" );
    axis9->ChangeLabel(3, -1, 0, -1, -1, 62, "" );
    //axis7->ChangeLabel(4, -1, -1, -1, -1, 62, "-" );
    //axis7->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis0->ChangeLabel(1, -1, 0, -1, -1, 62, "" );
    axis0->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi/2" );
    axis0->ChangeLabel(3, -1, 0, -1, -1, 62, "" );
    /*
    if(delphi){
        axis8->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi/2" );
        axis8->ChangeLabel(4, -1, -1, -1, -1, 62, "" );
        }
    else{
        axis8->ChangeLabel(3, -1, -1, -1, -1, 62, "" );
        axis8->ChangeLabel(4, -1, -1, -1, -1, 62, "" );
        }
        */
    axis9->SetLabelOffset(0.020);
    axis0->SetLabelOffset(0.0025);
    axis9->SetLabelSize(.03);
    axis0->SetLabelSize(.03);
    axis9->SetTickSize(0.5);
    axis0->SetTickSize(0.5);
    
    axis9->Draw();
    axis0->Draw();
    

    if (delphi) {
        if(!ratio) drawText("#pi/2", 0.505, .065, 17 );
        drawText("#pi/2", 0.505, .91, 17 );
        drawText("#Delta#phi", 0.8, .025, 16);
    } else {
        drawText("0.5", 0.490, .91, 16);
        drawText("z = #frac{(1 - cos(#theta))}{2}", .8, .02, 16);
    }
    
    drawText("PYTHIA8 preliminary", .15, .28, 15);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.24, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.20, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.16, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.12, 15);
    drawText(pttrackcut, 0.15, 0.08, 15);
    //drawText("Centrality 0-10%", 0.15, 0.04, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.04, 15);
    canvas1->Update();

    TLegend* legend = new TLegend(0.65, 0.1, 0.89, 0.25);
    if (all && !debug) {
        legend->AddEntry(EEC_w_p[0], "wEEC (Comb. Pythia + Hydro)", "pl");
        legend->AddEntry(EEC_w_p[1], "Pythia", "pl");
        legend->AddEntry(EEC_w_p[2], "Pythia + Hydro", "pl");
        legend->AddEntry(EEC_w_p[3], "Hydro + Hydro", "pl");
        legend->AddEntry(EEC_w_p[4], "Jet Event + Non Jet E", "pl");
        legend->AddEntry(EEC_w_p[5], "Non Jet E+ Non Jet E", "pl");
        legend->AddEntry(EEC_w_p[6], "Non Jet E 1 + + Non Jet E 2", "pl");
        legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl"); 
    } else if (debug && !all) {
        legend->AddEntry(EEC_w_p[0], "wEEC (Comb. Pythia + Hydro)", "pl");
        legend->AddEntry(EEC_w_p[1], "Pythia", "pl");
        legend->AddEntry(mbsstep1final, "wEEC - HP ", "pl"); 
        legend->AddEntry(mbsstep2final, "wEEC - HP - H1H2 1", "pl"); 
        legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl"); 
    } else if (debug && all) {
        legend->AddEntry(EEC_w_p[0], "wEEC (Comb. Pythia + Hydro)", "pl");
        legend->AddEntry(EEC_w_p[1], "Pythia", "pl");
        legend->AddEntry(EEC_w_p[2], "Pythia + Hydro", "pl");
        legend->AddEntry(EEC_w_p[3], "Hydro + Hydro", "pl");
        legend->AddEntry(EEC_w_p[4], "Jet Event + Non Jet E", "pl");
        legend->AddEntry(EEC_w_p[5], "Non Jet E+ Non Jet E", "pl");
        legend->AddEntry(EEC_w_p[6], "Non Jet E 1 + + Non Jet E 2", "pl");
        legend->AddEntry(mbsstep1final, "wEEC - HP ", "pl"); 
        legend->AddEntry(mbsstep2final, "wEEC - HP - H1H2 1", "pl"); 
        legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl"); 
    } else {
        legend->AddEntry(EEC_w_p[0], "wEEC ( Combined Pythia & Thermal Background )", "pl");
        legend->AddEntry(EEC_w_p[1], "Signal(Pythia)", "pl");
        legend->AddEntry(mbssignalfinal, "Mixed Background Subtraction", "pl"); 
    }

    legend->Draw();
    canvas1->SaveAs(outputFilename.c_str());
    canvas1->SaveAs(outputFilename1.c_str());
}