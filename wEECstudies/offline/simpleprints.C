#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}
using namespace std;
void simpleprints() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
  
    TFile* fin = TFile::Open("../../outfiles/ztotheta_2025-10-08_1.root");
    
    //TH2F* histetaphi = (TH2F*)fin->Get("hEtaPhi");
    //TH2F* histq2pmq2 = (TH2F*)fin->Get("comparisonpmq2");
    //TH2F* histqpmq = (TH2F*)fin->Get("comparisonpthat");
    //TH1F* histphi = (TH1F*)fin->Get("Deltaphi");
    TH1D* EEC_corr = (TH1D*)fin->Get("EEC_z");
    TH1D* EEC_theta = (TH1D*)fin->Get("EEC_theta");
    //TH1D* histsubleadingjetspec = (TH1F*)fin->Get("SubleadingJetSpectrum");
    //TH2D* hzvsdelphi = (TH2D*)fin->Get("zvsdelphi");
    //TH2D* hzvstheta = (TH2D*)fin->Get("zvstheta");
 
    //TH1D* theta = (TH1D*)fin->Get("htheta");
    TH1D* sm;
    TH1D* EEC_w_p;
    TH1D* EEC_w_p1;
    TCanvas* canvas1 = new TCanvas("canvas1", "z vs delphi", 800, 600);
    TCanvas* canvas2 = new TCanvas("canvas2", "z vs theta", 800, 600);
    TCanvas* canvas3 = new TCanvas("canvas3", "theta", 800, 600);
    //TCanvas* canvas3 = new TCanvas("canvas3", "Jet Spectrum", 800, 600);
    TCanvas* canvas4 = new TCanvas("canvas4", "Jet Spectrum leading", 800, 600);
    //TCanvas* canvas5 = new TCanvas("canvas5", "Jet Spectrum subleading", 800, 600);

    
    /*
    canvas1->cd();
    
    hzvsdelphi->SetStats(0);
    hzvsdelphi->GetXaxis()->SetRangeUser(0.0, M_PI);
    //hzvsdelphi>GetYaxis()->SetRangeUser(0.0, 6.28);
    //hzvsdelphi->GetZaxis()->SetRangeUser(0.0, 6.28);
    hzvsdelphi->GetYaxis()->SetTitle("#Delta#phi");
    hzvsdelphi->GetXaxis()->SetTitle("z = (1-cos(#theta)/2)");
    canvas1->SetLogz();
    hzvsdelphi->Draw("COLZ");
    //drawText("1M events", 0.1, .9, 18 );
    canvas1->SaveAs("hzvsdelphi.png");
    
    
    canvas2->cd();
    
    hzvstheta->SetStats(0);
    hzvstheta->GetYaxis()->SetTitle("#theta");
    hzvstheta->GetXaxis()->SetTitle("z = (1-cos(#theta)/2)");
    //hzvstheta->GetXaxis()->SetTitleOffset(1.2);
    canvas2->SetLogz();
    hzvstheta->Draw("COLZ");
    //drawText("1M events", 0.1, .9, 18 );
    
    canvas2->SaveAs("zvstheta.png");
    */
    canvas3->cd();
    /*
    TPad* pad1 = (TPad*)canvas1->cd(1);
    pad1->SetPad(0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0.0);
    
    TPad* pad2 = (TPad*)canvas1->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.30); //
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    
    pad1->cd();
    */
    double min = 0;
    double max = 1e10;
    int bins = 100;
    
    canvas3->Divide(1, 2); // Split into 2 vertical pads
    canvas3->SetLeftMargin(0.125); 
 
    TPad* pad1 = (TPad*)canvas3->cd(1);
    pad1->SetPad(0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0.0);
    pad1->SetLogy();
    
    TPad* pad2 = (TPad*)canvas3->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.30); //
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    
    pad1->cd();
    EEC_theta->SetStats(0);
    EEC_theta->GetYaxis()->SetTitle("EEC");
    EEC_theta->GetXaxis()->SetTitle("#theta ");
    //EEC_theta->GetXaxis()->SetRangeUser(0, 2* M_PI);
    //canvas3->SetLogz();
    EEC_theta->SetLineColor(kPink);
    EEC_theta->SetMarkerColor(kPink);
    EEC_theta->SetMarkerStyle(20);
    EEC_theta->SetMarkerStyle(22);
    
    EEC_w_p = new TH1D("EEC_w_p", "", 100, 0, 100);
    EEC_w_p1 = new TH1D("EEC_w_p1", "", 100, 0, 100);
    EEC_w_p->SetLineColor(kPink);
    EEC_w_p->SetMarkerColor(kPink);
    //EEC_w_p->SetMarkerStyle(20);
    for (int bin = 1; bin <= bins; ++bin) {
        EEC_w_p->SetBinContent(bin-30, EEC_theta->GetBinContent(bin));
        EEC_w_p->SetBinError(bin-30, EEC_theta->GetBinError(bin));
        double newcontent = EEC_theta->GetBinContent(bin);

    }
    for (int bin = 1; bin <= bins; ++bin) {
        float content = EEC_corr->GetBinContent(bin);
        //float bc = acos(1-2*content);
        EEC_w_p1->SetBinContent(bin-35, content);
        EEC_w_p1->SetBinError(bin-35, EEC_corr->GetBinError(bin));
       
    }
    EEC_w_p->Draw("P E1");
    EEC_w_p1->SetMarkerStyle(22);
    EEC_w_p1->Draw("P E1 SAME");
    
    TLegend* legend = new TLegend(0.65, 0.1, 0.89, 0.25);
    legend->AddEntry(EEC_w_p1, "EEC(acos(1-2*z))", "pl");
    legend->AddEntry(EEC_w_p, "EEC(#theta)", "pl");
    legend->Draw();

    pad2->cd();
    sm = (TH1D*)EEC_w_p->Clone();
    sm->SetMinimum(0.8);
    sm->SetMaximum(1.2);
    sm->Divide(EEC_w_p1);
    sm->Draw("P E1");

    /*
    pad2->cd();
    
    sm = (TH1D*)EEC_theta->Clone();
    sm->Divide(EEC_corr);
    //sm->Draw("P E1");
    
    sm->SetMinimum(0.8);
    sm->SetMaximum(1.2);
    sm->GetYaxis()->SetTitle(" EEC(#theta) / EEC(acos(1-2z))");
    */
    /*TLine *line = new TLine(0, 1, 3, 1);
    
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kOrange+2);
    line->Draw();
    */
    //drawText("1M events", 0.1, .9, 18 );
    
    canvas3->SaveAs("notlineareeccompoct13.png");

    canvas4->cd();
    sm = (TH1D*)EEC_theta->Clone();
    sm->SetMinimum(0.8);
    sm->SetMaximum(1.2);
    sm->Divide(EEC_corr);
    sm->Draw("P E1");
    canvas4->SaveAs("ratiothetatocorrzoct13.png");
    
    }
