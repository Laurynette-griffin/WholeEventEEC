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
void MixedBackgrounfSubdoubles() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
   /* const char* filename = {
       "../outfiles/pythia_mixedevent_03_15pthat50_May27_50m.root"};
       //"../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_May22_50m.root"};
       */
    TFile* file = TFile::Open("../outfiles/pythia_pp_fullevent_5pthat60_mbs_40bins_5M.root");

    vector<const char*> histonames = {"EEC_wpl", "EEC_mspl", "EEC_mmpl", "EEC_m2pl", "EEC_spl", "EEC_tpl", "EEC_mpl"};
    cout <<"is this thing onnnnn" << endl;   
       
    int const number = 7;
    double min = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e6; 
    
    bool rebin = 0;
    bool delphi = 0;
    
    TH1D* histos[number];
    TH1D* prehistos[number];
    TH1D* preEEC_w_p[number];
    TH1D* EEC_w_p[number];
    TH1D* mbssignal;
    TH1D* mbssignal1final;
    TH1D* mbssignal1;
    TH1D* mbssignalfinal;
    
    string outputFilename = "../plots/weeccentrality03sept8-40binlp.C";
    string outputFilename1 = "../plots/weecentrality03sept8-40binlp.png";
    const char* pttrackcut = "charged constituent p_{T} > 0.2 GeV";
    
    int bins = 40;

    for (int i = 0; i < number; ++i) {
        TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
        // to keep it separte and not mess up your root file
        prehistos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
        prehistos[i]->SetDirectory(0);
        
        histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
        histos[i]->SetDirectory(0);
        for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
            double width = histos[i]->GetBinWidth(bin);
            double content = histos[i]->GetBinContent(bin);
            double error = histos[i]->GetBinError(bin);
            histos[i]->SetBinContent(bin, content / width);
            histos[i]->SetBinError(bin, error / width);
        }
        
        EEC_w_p[i] = (TH1D*)histos[i]->Clone(Form("EEC_clone_%d", i));
        EEC_w_p[i]->SetDirectory(0);
        
        //EEC_w_p[i]->Scale(1.0/ EEC_w_p[i]->Integral("width"));
 
    }
        
    mbssignal = (TH1D*)prehistos[0]->Clone();
    mbssignal->Add(prehistos[1], -1);
    mbssignal->Add(prehistos[2], -0.1);
    mbssignal->Add(prehistos[3], +0.1);
    
    
        for (int bin = 1; bin <= mbssignal->GetNbinsX(); ++bin) {
            double width = mbssignal->GetBinWidth(bin);
            double content = mbssignal->GetBinContent(bin);
            double error = mbssignal->GetBinError(bin);
                //mbssignal->SetBinContent(bin, content / width);
                //mbssignal->SetBinError(bin, error / width);
    }
    
    //mbssignal->Scale(1.0/mbssignal->Integral("width"));
    mbssignalfinal = (TH1D*)mbssignal->Clone("mbssignalfinal");
    
    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 600);
    canvas1->cd();
    canvas1->SetLogy();
    canvas1->SetLeftMargin(0.125); 
 
    // empty hist since we are rebinning on the edges 
    TH1D* EEC_empty = new TH1D("EEC_empty", "", 40, 0, M_PI);
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
            EEC_w_p[i]->SetLineColor(kYellow+1);
            EEC_w_p[i]->SetMarkerColor(kYellow+1);
        }
         else if(i == 5) {
            EEC_w_p[i]->SetLineColor(46);
            EEC_w_p[i]->SetMarkerColor(46);
        }
        
        else if(i == 6) {
            EEC_w_p[i]->SetLineColor(kYellow+1);
            EEC_w_p[i]->SetMarkerColor(kYellow+1);
        }

        else { // 0 and 3
            EEC_w_p[i]->SetLineColor(i + 1);
            EEC_w_p[i]->SetMarkerColor(i + 1);
            }
        // just plot whole & signal 
        if (i < 6 ){
        //if(i == 0 || i ==4){
        EEC_w_p[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }
        
    }
    mbssignalfinal->SetMarkerStyle(20);
	mbssignalfinal->SetMarkerColor(kCyan+3);
	mbssignalfinal->SetLineColor(kCyan+3);
	mbssignalfinal->Draw("P E1 SAME");
	
    TGaxis *axis1 = new TGaxis(0, 0, 20, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(40, 0, 20, 0, .00001 , .5 , 510, "-G");
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
    
 
    TGaxis *axis3 = new TGaxis(M_PI, min , M_PI, max , min, max, 510,"+LG");
    axis3->SetLabelFont(42);
    axis3->SetLabelSize(.035);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(0, max, M_PI, max, 0, M_PI, 510, "-"); 
    TGaxis *axis5 = new TGaxis(0, min, M_PI, min, 0, M_PI, 510, "-"); 
    
     //axis4->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
   // axis4->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    //axis4->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
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
    
    //axis1->Draw();
    //axis2->Draw();
    axis4->Draw();
    //axis5->Draw();
	
    drawText("#Delta#phi", 0.8, .025, 16);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.32, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.28, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.24, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.20, 15);
    drawText(pttrackcut, 0.15, 0.16, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.12, 15);
    canvas1->Update();
    
    TLegend* legend = new TLegend(0.55, 0.125, 0.875, 0.375);
    legend->AddEntry(EEC_w_p[0], "wEEC", "pl");
    legend->AddEntry(EEC_w_p[1], "Mixed1 + Signal", "pl");
    legend->AddEntry(EEC_w_p[3], "Mixed1 + Mixed2", "pl");
    legend->AddEntry(EEC_w_p[2], "Mixed1 + Mixed1 ", "pl");
    
    legend->AddEntry(EEC_w_p[4], "Signal + Signal", "pl");
    legend->AddEntry(EEC_w_p[5], "Thermal + Thermal", "pl");
    //legend->AddEntry(EEC_w_p[2], "Signal + Mixed", "pl");
    //1legend->AddEntry(EEC_w_p[3], "Mixed + Mixed", "pl");
    //legend->AddEntry(EEC_w_p[6], "Signal + Signal", "pl");
   // legend->AddEntry(mbssignal1final, "Mixed background subtraction post bin", "pl");
    legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl");
    
    legend->Draw();
 
    canvas1->SaveAs(outputFilename.c_str());
    canvas1->SaveAs(outputFilename1.c_str());

}





    