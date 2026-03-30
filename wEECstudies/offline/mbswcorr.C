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
void mbswcorr() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetPadTickY(1);
   /* const char* filename = {
       "../outfiles/pythia_mixedevent_03_15pthat50_May27_50m.root"};
       //"../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_May22_50m.root"};
       pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_mbs_thermalcopy_Aug28_1M.root
       */
    TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_mbs_4thermalevents_100k_Nov19_2m.root");
    //TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_mbs_2dhists_Sep22_100k.root");
    //TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy010_psi2match_trackptcut05gev_20centbins_Mar6_5k.root");
    if (!file || file->IsZombie()) {
        cout << "Error opening file!" << endl;
        return;
    }
    int const number = 17;
    
    double min = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e8;
    
    double min1 = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max1 = 1e1;
    
    double minr = 0.0;
    double maxr = 2.0;
    
    bool rebin = 0;
    bool delphi = 0;
    bool integrate = 0;
    bool linear = 0;
    bool scale = 0;
    bool all = 0;
    bool debug = 0;
    int bins = 60;
 
    vector<const char*> histonames;
 
    if(scale){
         histonames = {"EEC_w", "EEC_s", "EEC_m", "EEC_t", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_mt", 
        "EEC_pt1", "EEC_pt2","EEC_pt3", "EEC_t1t1", "EEC_t2t2", "EEC_t3t3", "EEC_t1t2", "EEC_t1t3", "EEC_t2t3"};
    }
    else if(delphi && linear){
        histonames = {"EEC_wpl", "EEC_spl", "EEC_mpl", "EEC_tpl", "EEC_mspl", "EEC_mmpl", "EEC_m2pl"};
    }
    else if(delphi){
        histonames = {"EEC_wp", "EEC_sp", "EEC_mp", "EEC_tp", "EEC_msp", "EEC_mmp", "EEC_m2p"};
    }
    else{
        histonames = {"EEC_w", "EEC_s", "EEC_m", "EEC_t", "EEC_ms", "EEC_mm", "EEC_m2"};
    }
    
    
    TH1D* histos[number];
    TH1D* prehistos[number];
    TH1D* preEEC_w_p[number];
    TH1D* EEC_w_p[number];
    TH1D* mbssignal;
    TH1D* mbssignal1final;
    TH1D* mbssignal1;
    TH1D* mbssignalfinal;
    TH1D* rh;
    TH1D* rh1;
    TH1D* rh2;
    TH1D* msms;
    TH1D* t12t23;
    TH1D* sm;
    TH1D* t11t22;
    TH1D* t22t33;
    TH1D* mbsstep1;
    TH1D* mbsstep2;
    TH1D* mbsstep1final;
    TH1D* mbsstep2final;
    string baseName = "../plots/weec_mbthermal_regen_scale_Mar18";
    
    string outputFilename  = baseName + ".C";
    string outputFilename1 = baseName + ".png";
    string outputFiledebug = baseName + ".root";
    
    const char* pttrackcut = "charged constituent p_{T} > 0.5 GeV";
    
    if (linear){
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
            
            //histos[i]->Scale(1.0/ histos[i]->Integral("width"));
            
            EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", bins, 0, bins);
        
            for (int bin = 1; bin <= bins; ++bin) {
                EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
                EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
                double newcontent = histos[i]->GetBinContent(bin);
                //cout << newcontent << " " << bin << endl;
                if (newcontent > 0 && newcontent < min) min = newcontent;
                if (newcontent > max) max = newcontent;
            }
        }
    }
    else{     
        for (int i = 0; i < number; ++i) {
            TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
            // to keep it separte and not mess up your root file
            histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            histos[i]->SetDirectory(0);
            
            // histos that will be subtracted 
            prehistos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            prehistos[i]->SetDirectory(0);
            // Normalize to bin width
            for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
                double width = histos[i]->GetBinWidth(bin);
                double content = histos[i]->GetBinContent(bin);
                double error = histos[i]->GetBinError(bin);
                histos[i]->SetBinContent(bin, content / width);
                
                if (i == 0 || i == 4 || i == 5 || i == 6){
                    
                    cout << "Bin content of prehistos "<< i << " bin number = " << bin << " is = to "<< prehistos[i]->GetBinContent(bin)  << endl;
                    }
                }
            
        //mbs stuff only
        }
        for (int i = 0; i < number; ++i) {        
            if (integrate){
            histos[i]->Scale(1.0/ histos[i]->Integral("width"));
            }
    
            EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", bins, 0, bins);
            
            for (int bin = 1; bin <= bins; ++bin) {
                EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
                EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
                double newcontent = histos[i]->GetBinContent(bin);
                cout << "Bin content of histos" << i << " "<< bin << " is = to "<< newcontent  << endl;
                //cout << newcontent << " " << bin << endl;
                if (newcontent > 0 && newcontent < min) min = newcontent;
                if (newcontent > max) max = newcontent;
            }
        } //closes i = loop
    }
    //pre bin width SUBTRACTION IS GONNA OCCUR NOW 
    //Subhistos are subtracted then normalized by bin width then integral 
    // okay so EEC_w - EEC_ms - EEC_mm + EEC_m2 = EEC_s allegedly and thats what this code does :)

    t11t22 = (TH1D*)EEC_w_p[13]->Clone();
    t11t22->Divide(EEC_w_p[12]);
    
    t12t23 =  (TH1D*)EEC_w_p[14]->Clone(); // ask rithya why is this "better than the opposite which feels more correct "
    t12t23->Divide(EEC_w_p[16]);
    
    msms = (TH1D*)prehistos[4]->Clone();
    msms->Add(prehistos[10], +1.0);
    msms->Add(prehistos[15], +1.0);
    msms->Scale(1.0/2.0);
    
    t22t33= (TH1D*)prehistos[12]->Clone();
    t22t33->Add(prehistos[13], +1.0);
    t22t33->Scale(1.0/2.0);
    
    
    mbssignal = (TH1D*)prehistos[0]->Clone();

    mbssignal->Add(prehistos[4], -1.0);
    mbsstep1 = (TH1D*)mbssignal->Clone();
    mbssignal->Add(prehistos[5], -1.0);
    mbsstep2 = (TH1D*)mbssignal->Clone();
    mbssignal->Add(prehistos[6], +1.0); //opposite sign of incorrect background modeling in sm1 so this IS the right sign :) (and the only way it works)

    
    if (scale){
        for (int i = 1; i <= bins; i++) {
            double content = prehistos[0]->GetBinContent(i);
            
            double scalect11t22 = t11t22->GetBinContent(i);
            double scalect12t23 = t12t23->GetBinContent(i);
            //if (i == 30 )cout << "scalect11t33: " << scalect11t33 << "    scalect12t23:   " << scalect12t23 << endl;
            content -= prehistos[4]->GetBinContent(i); // this is pythia + thermal so I cannot scale this
            content -= scalect11t22 * (prehistos[5]->GetBinContent(i)); // scale these though 
            //content += scalect12t23 * (prehistos[6]->GetBinContent(i));
            content += (prehistos[6]->GetBinContent(i));
            
            double scalert11t22 = t11t22->GetBinError(i);
            double scalert12t23 = t12t23->GetBinError(i);
            
            double err0 = prehistos[0]->GetBinError(i);
            double err4 = prehistos[4]->GetBinError(i); 
            double err5 = scalert11t22 * (prehistos[5]->GetBinError(i)); //and these 
            double err6 = scalert12t23 * (prehistos[6]->GetBinError(i));
    
            double error = std::sqrt(err0 * err0 + err4 * err4 + err5 * err5 + err6 * err6);
    
            mbssignal->SetBinContent(i, content);
            mbssignal->SetBinError(i, error);
        }
    }
    
    //mbs has occured now to normaliza just mbssignal by bin width
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
        //cout << newcontent << " " << bin << endl;
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
    if (integrate){
        mbssignalfinal->Scale(1.0/mbssignalfinal->Integral("width"));
        }
    
    if(debug){
        //! MBS step 1
        for (int bin = 1; bin <= mbsstep1->GetNbinsX(); ++bin) {
            double width = mbsstep1->GetBinWidth(bin);
            double content = mbsstep1->GetBinContent(bin);
            double error = mbsstep1->GetBinError(bin);
            mbsstep1final->SetBinContent(bin, content / width);
            mbsstep1final->SetBinError(bin, error/width);
            double newcontent = mbsstep1final->GetBinContent(bin);
            //cout << newcontent << " " << bin << endl;
            if (newcontent > 0 && newcontent < min) min = newcontent;
            if (newcontent > max) max = newcontent;
        }
        if (integrate){
            mbsstep1final->Scale(1.0/mbsstep1final->Integral("width"));
            }
        
         for (int bin = 1; bin <= mbsstep2->GetNbinsX(); ++bin) {
            double width = mbsstep2->GetBinWidth(bin);
            double content = mbsstep2->GetBinContent(bin);
            double error = mbsstep2->GetBinError(bin);
            mbsstep2final->SetBinContent(bin, content / width);
            mbsstep2final->SetBinError(bin, error/width);
            double newcontent = mbsstep2final->GetBinContent(bin);
            //cout << newcontent << " " << bin << endl;
            if (newcontent > 0 && newcontent < min) min = newcontent;
            if (newcontent > max) max = newcontent;
        }
        if (integrate){
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
    }
    
    //plot :)
    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 800);
    canvas1->Divide(1, 2); // Split into 2 vertical pads
    canvas1->SetLeftMargin(0.125); 
    
    TPad* pad1 = (TPad*)canvas1->cd(1);
    pad1->SetPad(0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0.0);
    pad1->SetLogy();
    
    TPad* pad2 = (TPad*)canvas1->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.30); //
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    
    pad1->cd();
    
    // empty hist since we are rebinning on the edges 
    TH1D* EEC_empty = new TH1D("EEC_empty", "", bins, 0, bins);
    EEC_empty->SetStats(0);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1 / (#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2} EEC");
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    for (int i = 0; i < 7; ++i) {
        if (i == 2 || i == 3) continue;
        //marker size changes & style assignments 
        EEC_w_p[i]->SetStats(0);
        //20 is filled in circle 21 us
        if (i < 2) {
            EEC_w_p[i]->SetMarkerStyle(20 + i);   
            if (i == 1){
            EEC_w_p[i]->SetMarkerSize(1.25);
            }
        }
        else if (i == 2 ){
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
            EEC_w_p[i]->SetLineColor(kOrange+1);
            EEC_w_p[i]->SetMarkerColor(kOrange+1);
        }
        else if(i == 4) {
            EEC_w_p[i]->SetLineColor(kYellow+1);
            EEC_w_p[i]->SetMarkerColor(kYellow+1);
        }
         else if(i == 5) {
            EEC_w_p[i]->SetLineColor(46);
            EEC_w_p[i]->SetMarkerColor(46);
        }
        else if(i == 6) {
            EEC_w_p[i]->SetLineColor(kViolet-1);
            EEC_w_p[i]->SetMarkerColor(kViolet-1);
        }
        else { // 0 and 3
            EEC_w_p[i]->SetLineColor(i + 1);
            EEC_w_p[i]->SetMarkerColor(i + 1);
            }
        if(all){
            if (i >= 0 ){
                EEC_w_p[i]->Draw("P E1 SAME");
                EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
                EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
                EEC_w_p[i]->GetXaxis()->SetTickLength(0);
            }
        }
        else{
            if (i <= 1 ){
                EEC_w_p[i]->Draw("P E1 SAME");
                EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
                EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
                EEC_w_p[i]->GetXaxis()->SetTickLength(0);
            }
        }
    }
    
    //mbssignal->Draw("P E1 SAME");
    mbssignalfinal->SetMarkerStyle(107);
	mbssignalfinal->SetMarkerColor(618);
	mbssignalfinal->SetLineColor(618);
	mbssignalfinal->Draw("P E1 SAME");
	if(debug){
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
    //sm = (TH1D*)EEC_w_p[1]->Clone();
    sm = (TH1D*)mbssignalfinal->Clone();
    rh->Divide(sm);
    
    rh->SetTitle("");
    rh->GetYaxis()->SetTitle("Pythia / Bkd Sub");
    //rh->GetYaxis()->SetTitle("EEC_corr / EEC_theta");
    //rh->GetYaxis()->SetTitleSize(10.);
    //rh->GetYaxis()->SetTitleOffset(999);
    rh->GetYaxis()->SetLabelSize(.05);
    //rh->GetYaxis()->SetTickLength(0);
    //rh->GetXaxis()->SetLabelSize(0.1);
    rh->SetMarkerStyle(20);
    rh->SetMarkerSize(0.75);
    rh->SetStats(0);
    rh->SetLineColor(kBlack);
    rh->SetMarkerColor(kBlack);
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
    axis7->ChangeLabel(2, -1, -1, -1, -1, 62, "10^{-4}" );
    axis7->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis7->ChangeLabel(4, -1, -1, -1, -1, 62, "10^{-2}" );
    axis7->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis8->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis8->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis8->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
    axis8->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
    axis8->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis8->SetLabelOffset(0.125);
    axis7->SetLabelSize(.075);
    axis8->SetLabelSize(.075);
    axis7->SetTickSize(0.5);
    axis8->SetTickSize(0.5);
    
    //TGaxis *axis9 = new TGaxis(0, .95, 0, 1.05, .95, 1.05, 1005, ""); 
    //TGaxis *axis10 = new TGaxis(bins-1, .95, bins-1, 1.05, .95, 1.05, 1005, ""); 
 
    axis7->Draw();
    axis8->Draw();
    //axis9->Draw();
    //axis10->Draw();
    // Axis limits
    drawText("0.5", 0.4850, .21, 16 );
    /*
    TLegend* legend1 = new TLegend(0.75, 0.125, 0.975, 0.375);
    legend1->AddEntry(rh, " Pythia / mbs ", "pl");
    legend1->Draw();
    */
    pad1->cd();
    
    TGaxis *axis1 = new TGaxis(0, 0, bins/2, 0, .0001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(bins, 0, bins/2, 0, .0001 , .5 , 510, "-G");
    axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(2, -1, 0, -1, -1, 62, "10^{-4}" );
    axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(4, -1, 0, -1, -1, 62, "10^{-2}" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(2, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(4, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    if(delphi){
        axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-5}" );
        axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-3}" );
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
        //drawText("0.5", 0.490, .025, 17 );
        drawText("0.5", 0.490, .91, 16);
        //drawText("0.5", 0.490, .01, 17 );
        drawText("z = #frac{(1 - cos(#theta))}{2}", .8, .02, 16);
    }
    
    /*
    TGaxis *axis1 = new TGaxis(0, 0, bins/2, 0, .0001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(bins, 0, bins/2, 0, .0001 , .5 , 510, "-G");
    axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(2, -1, 0, -1, -1, 62, "10^{-4}" );
    axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(4, -1, 0, -1, -1, 62, "10^{-2}" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(2, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(4, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    if(delphi){
        axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-5}" );
        axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-3}" );
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
        //drawText("0.5", 0.490, .025, 17 );
        drawText("0.5", 0.490, .91, 16);
        //drawText("0.5", 0.490, .01, 17 );
        drawText("z = #frac{(1 - cos(#theta))}{2}", .8, .02, 16);
    }
    */
    
    /*if(delphi){
    //rawText("NO MPI", 0.15, 0.36, 15);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.32, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.28, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.24, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.20, 15);
    drawText(pttrackcut, 0.15, 0.16, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.12, 15);
    canvas1->Update();
    }*/
    
    drawText("PYTHIA8 + JETSCAPE soft sector preliminary", .15, .28, 15);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.24, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.20, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.16, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.12, 15);
    drawText(pttrackcut, 0.15, 0.08, 15);
    drawText("Centrality 0-10%", 0.15, 0.04, 15);
    canvas1->Update();

    TLegend* legend = new TLegend(0.65, 0.1, 0.89, 0.25);
    if(all && !debug){
        legend->AddEntry(EEC_w_p[0], "wEEC (Comb. Pythia + Hydro)", "pl");
        legend->AddEntry(EEC_w_p[1], "Pythia", "pl");
        legend->AddEntry(EEC_w_p[2], "Pythia + Hydro", "pl");
        legend->AddEntry(EEC_w_p[3], "Hydro + Hydro", "pl");
        legend->AddEntry(EEC_w_p[4], "Jet Event + Non Jet E", "pl");
        legend->AddEntry(EEC_w_p[5], "Non Jet E+ Non Jet E", "pl");
        legend->AddEntry(EEC_w_p[6], "Non Jet E 1 + + Non Jet E 2", "pl");
        legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl"); 
    }
    else if(debug && !all){
        legend->AddEntry(EEC_w_p[0], "wEEC (Comb. Pythia + Hydro)", "pl");
        legend->AddEntry(EEC_w_p[1], "Pythia", "pl");
        legend->AddEntry(mbsstep1final, "wEEC - HP ", "pl"); 
        legend->AddEntry(mbsstep2final, "wEEC - HP - H1H2 1", "pl"); 
        legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl"); 
        
    }
    else if(debug && all){
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
    }
    else{
        legend->AddEntry(EEC_w_p[0], "wEEC (Comb. Pythia + soft sector hydro)", "pl");
        legend->AddEntry(EEC_w_p[1], "Pythia", "pl");
        legend->AddEntry(mbssignalfinal, "MBS [wEEC - HP - H1H2 + H1H1]", "pl"); 
    }

    legend->Draw();
    canvas1->SaveAs(outputFilename.c_str());
    canvas1->SaveAs(outputFilename1.c_str());
}

    
    
    
        
