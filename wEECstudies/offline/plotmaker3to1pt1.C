#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(42);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

//takes 3/4 hists from a file and plots them on one canvas 
void plotmaker3to1() {
    
    bool delphi = false;
    bool rebin = false;
    
    TH1::SetDefaultSumw2();
    
    TFile* fin = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_pt1_on_April23_50mil.root");
    
    const char* histNames[] = {"EEC_w", "EEC_w_high", "EEC_w_mid", "EEC_w_low"};
    for (int i = 0; i < nHists; ++i) {
    TH1F* htemp = dynamic_cast<TH1F*>(fin->Get(histNames[i]));
    TH1F* histEECwhole = (TH1F*)fin->Get("EEC_w");
    TH1F* histEEChigh = (TH1F*)fin->Get("EEC_w_high"); 
    TH1F* histEECmid = (TH1F*)fin->Get("EEC_w_mid");
    TH1F* histEEClow = (TH1F*)fin->Get("EEC_w_low"); 

    TCanvas* canvas5 = new TCanvas("canvas5", "Comparison", 800, 600);
    
    //axis stuff 
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


TH1F* histo1 = (TH1F*)histEECwhole->Clone("histo1");

// Normalize to bin width
for (int bin = 1; bin <= histo1->GetNbinsX(); ++bin) {
        float width = histo1->GetBinWidth(bin);
        float content = histo1->GetBinContent(bin);
        float error = histo1->GetBinError(bin);
        histo1->SetBinContent(bin, content / width);
        histo1->SetBinError(bin, error / width);
    }
    
    histo1->Scale(1.0 / histo1->Integral("width"));
    
    TH1F* EEC_w_p = new TH1F("EEC_w_p", "WEEC", 60, 0, 60);
    
    for (int bin = 1; bin <= 60; ++bin) {
        
        EEC_w_p->SetBinContent(bin, histo1->GetBinContent(bin));
        EEC_w_p->SetBinError(bin, histo1->GetBinError(bin));
    }
    
    
    EEC_w_p->SetStats(0);
    EEC_w_p->Draw("P E1");
    EEC_w_p->GetXaxis()->SetLabelOffset(999);
    EEC_w_p->GetXaxis()->SetTitleOffset(999);
    EEC_w_p->GetXaxis()->SetTickLength(0);
    EEC_w_p->SetStats(0);
    EEC_w_p->GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");

    
    TH1F* histo2 = (TH1F*)histEEChigh->Clone("histo2");
    for (int bin = 1; bin <= histo2->GetNbinsX(); ++bin) {
        float width = histo2->GetBinWidth(bin);
        float content = histo2->GetBinContent(bin);
        float error = histo2->GetBinError(bin);
        histo2->SetBinContent(bin, content / width);
        histo2->SetBinError(bin, error);
    }
    
    histo2->Scale(1.0 / histo2->Integral("width"));
    TH1F* EEC_w_h = new TH1F("EEC_w_h", "WEEC", 60, 0, 60);
    
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_h->SetBinContent(bin, histo2->GetBinContent(bin));
        EEC_w_h->SetBinError(bin, histo2->GetBinError(bin));
    }
    
    
    TH1F* histo3 = (TH1F*)histEECmid->Clone("histo3");
    
    for (int bin = 1; bin <= histo3->GetNbinsX(); ++bin) {
        float width = histo3->GetBinWidth(bin);
        float content = histo3->GetBinContent(bin);
        float error = histo3->GetBinError(bin);
        histo3->SetBinContent(bin, content / width);
        histo3->SetBinError(bin, error);
    }
    histo3->Scale(1.0 / histo3->Integral("width"));
    TH1F* EEC_w_m = new TH1F("EEC_w_m", "WEEC", 60, 0, 60);
  
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_m->SetBinContent(bin, histo3->GetBinContent(bin));
        EEC_w_m->SetBinError(bin, histo3->GetBinError(bin));
    }


    TH1F* histo4 = (TH1F*)histEEClow->Clone("histo4");
    
    for (int bin = 1; bin <= histo4->GetNbinsX(); ++bin) {
        float width = histo4->GetBinWidth(bin);
        float content = histo4->GetBinContent(bin);
        float error = histo4->GetBinError(bin);
        histo4->SetBinContent(bin, content / width);
        histo4->SetBinError(bin, error / width);
    }
    
    histo4->Scale(1.0 / histo4->Integral("width"));
    TH1F* EEC_w_l = new TH1F("EEC_w_l", "WEEC", 60, 0, 60);
 
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_l->SetBinContent(bin, histo4->GetBinContent(bin));
        EEC_w_l->SetBinError(bin, histo4->GetBinError(bin));
    }
    

    axis1->Draw();
    axis2->Draw();
    
    if(delphi){
        drawText("#pi/2", 0.4875, .06, 17 );
        drawText("#Delta#phi", 0.8, .025, 16);
    }
    else{
        drawText("0.5", 0.4875, .06, 17 );
        drawText("z = #frac{(1 - cos(#theta))}{2}", 0.8, .025, 16);
    }
    drawText("p+p #sqrt{s} = 91.2 GeV", 0.15, 0.28, 16);
    drawText("43 GeV < jet p_{T1} <= 45 GeV", 0.15, 0.25, 16);
    drawText("43 GeV < jet p_{T2} <= 35 GeV", 0.15, 0.22, 16);
    drawText("44.60 GeV < #hat{p}_{T} <= 45.25 GeV", 0.15, 0.19, 16);
    drawText("charged constituent p_{T} > 0.2 GeV", 0.15, 0.16, 16);
    drawText("#Delta#phi > 7#pi/8", 0.15, 0.13, 16);
    
    canvas5->cd();
    canvas5->SetLogy();
    
    EEC_w_l ->SetLineColor(kRed);
    EEC_w_m->SetLineColor(kBlue);
    EEC_w_h->SetLineColor(kMagenta);
    EEC_w_p->SetLineColor(kGreen+2);

    EEC_w_l->SetMarkerStyle(20);
    EEC_w_l->SetMarkerColor(kRed);
    
    EEC_w_m->SetMarkerStyle(21);
    EEC_w_m->SetMarkerColor(kBlue);
    
    EEC_w_h->SetMarkerStyle(23);
    EEC_w_h->SetMarkerColor(kMagenta);

    EEC_w_p->SetMarkerStyle(24);
    EEC_w_p->SetMarkerColor(kGreen+2);

    EEC_w_l->SetMarkerSize(0.5);
    EEC_w_m->SetMarkerSize(0.5);
    EEC_w_h->SetMarkerSize(0.5);
    EEC_w_p->SetMarkerSize(0.5);
    
    EEC_w_l->SetMinimum(1e-1);
    EEC_w_l->SetMaximum(100);
    
    EEC_w_l->Draw("P E1");
    EEC_w_m->Draw("P E1 SAME");
    EEC_w_h->Draw("P E1 SAME");
    //EEC_w_p->Draw("P E1 SAME");
    
    axis1->Draw();
    axis2->Draw();
    
    drawText("0.5", 0.4875, .06, 17 );
    //drawText("z = (1 - cos(#theta))/2", 0.8, .025, 16);
    drawText("#Delta#Phi", 0.8, .025, 18 );
    drawText("p+p #sqrt{s} = 200 GeV", 0.20, 0.88, 16);
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.84, 16);
    drawText("jet p_{T2} >= 9.4 GeV", 0.20, 0.80, 16);
    drawText("charged constituent p_{T} > .2 GeV", 0.20, .76, 16 );
    drawText("#Delta#Phi > 3/4 ", 0.20, .72, 16 );
    
    TLegend* legend = new TLegend(0.6, 0.13, 0.9, 0.35);
    legend->AddEntry(EEC_w_p, "20.9 Gev < jet p_{T1} < 60.8 GeV", "pl");
    legend->AddEntry(EEC_w_l, "20.9 Gev < jet p_{T1} < 31.2 GeV", "pl");
    legend->AddEntry(EEC_w_m, "31.2 GeV <= jet p_{T1} < 40.7 GeV", "pl");
    legend->AddEntry(EEC_w_h, "40.7 GeV <= jet p_{T1} < 60.8 GeV", "pl");

    legend->Draw();

    canvas5->SaveAs("Variablept1phiMay16.png");
}
