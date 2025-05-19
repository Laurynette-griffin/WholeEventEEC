#include <TFile.h>
#include <TH1F.h>
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

void plotmaker_eec15() {
    
    TH1::SetDefaultSumw2();
    //TFile* fin = TFile::Open("../offline/pythia_pp_rhic_15GeV_fullevent_10pthat60_March3_mpioff.root"); // mpi off
    TFile* fin = TFile::Open("../outfiles/pythia_pp_rhic_8GeV_fullevent_5pthat60_April18.root");  
    //TFile* fin = TFile::Open("../offline/pythia_pp_rhic_15GeV_fullevent_10pthat60_March3_hadoff.root");
    
    TH1F* histEECwhole = (TH1F*)fin->Get("EEC_w");
    
    TH1F* histEEChigh = (TH1F*)fin->Get("EEC_w_high"); // aj < .1
    TH1F* histEECmid = (TH1F*)fin->Get("EEC_w_mid"); // aj < .667
    TH1F* histEEClow = (TH1F*)fin->Get("EEC_w_low"); // aj < .15003
    
    TH1F* histjetspec = (TH1F*)fin->Get("JetSpectrum");
    TH1F* histleadingjetspec = (TH1F*)fin->Get("LeadingJetSpectrum");
    TH1F* histsubleadingjetspec = (TH1F*)fin->Get("SubleadingJetSpectrum");
    
    
    TCanvas* canvas1 = new TCanvas("canvas1", "Whole event EEC", 800, 600);
    TCanvas* canvas2 = new TCanvas("canvas2", "Whole event EEC high", 800, 600);
    TCanvas* canvas3 = new TCanvas("canvas3", "Whole event EEC mid", 800, 600);
    TCanvas* canvas4 = new TCanvas("canvas4", "Whole event EEC low", 800, 600);
    //TCanvas* canvas5 = new TCanvas("canvas5", "Jet Spectrum", 800, 600);
    //TCanvas* canvas6 = new TCanvas("canvas6", "Jet Spectrum leading", 800, 600);
    //TCanvas* canvas7 = new TCanvas("canvas7", "Jet Spectrum subleading", 800, 600);
    TCanvas* canvas5 = new TCanvas("canvas5", "Comparison", 800, 600);
    
    //axis stuff 
    TGaxis *axis1 = new TGaxis(0, 0, 30, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(60, 0, 30, 0, .00001 , .5 , 510, "-G");
    axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
    axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->SetLabelOffset(0.045);
    axis1->SetLabelSize(.03);
    axis2->SetLabelSize(.03);
    axis1->SetTickSize(0.05);
    axis2->SetTickSize(0.05);

canvas1->cd();
canvas1->SetLogy();

// Clone and rename to keep each histogram separate
TH1F* histo1 = (TH1F*)histEECwhole->Clone("histo1");

// Normalize to bin width
for (int bin = 1; bin <= histo1->GetNbinsX(); ++bin) {
        float width = histo1->GetBinWidth(bin);
        float content = histo1->GetBinContent(bin);
        float error = histo1->GetBinContent(bin);
        histo1->SetBinContent(bin, content / width);
        histo1->SetBinError(bin, error / width);
    }
    
    TH1F* EEC_w_p = new TH1F("EEC_w_p", "WEEC", 60, 0, 60);
    
    for (int bin = 1; bin <= 60; ++bin) {
        
        EEC_w_p->SetBinContent(bin, histo1->GetBinContent(bin));
        EEC_w_p->SetBinError(bin, histo1->GetBinError(bin));
    }
    
    EEC_w_p->Scale(1.0 / EEC_w_p->Integral("width"));
    
    EEC_w_p->SetStats(0);
    EEC_w_p->Draw("P E1");
    EEC_w_p->GetXaxis()->SetLabelOffset(999);
    EEC_w_p->GetXaxis()->SetTitleOffset(999);
    EEC_w_p->GetXaxis()->SetTickLength(0);
    EEC_w_p->SetStats(0);
    
    axis1->Draw();
    axis2->Draw();
    
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
    drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
    //drawText("jet p_{T1} >= 20.9 GeV", 0.20, 0.79, 18);
    drawText("31.2 Gev < jet p_{T1} >= 40.7 GeV", 0.20, 0.79, 18);
    drawText("9.4 GeV < jet p_{T2} < 20.9 GeV", 0.20, 0.75, 18);
    //drawText("jet p_{T1} >= 9.4 GeV", 0.20, 0.75, 18);
    drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
    canvas1->SaveAs("WholeEventEECallApril17.png");
    
    canvas2->cd();
    canvas2->SetLogy();
    
    TH1F* histo2 = (TH1F*)histEEChigh->Clone("histo2");
    for (int bin = 1; bin <= histo2->GetNbinsX(); ++bin) {
        float width = histo2->GetBinWidth(bin);
        float content = histo2->GetBinContent(bin);
        float error = histo2->GetBinContent(bin);
        histo2->SetBinContent(bin, content / width);
        histo2->SetBinError(bin, error);
    }
    
    TH1F* EEC_w_h = new TH1F("EEC_w_h", "WEEC", 60, 0, 60);
    
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_h->SetBinContent(bin, histo2->GetBinContent(bin));
        EEC_w_h->SetBinError(bin, histo2->GetBinError(bin));
    }
    
    EEC_w_h->Scale(1.0 / EEC_w_h->Integral("width"));
    
    EEC_w_h->SetStats(0);
    EEC_w_h->Draw("P E1");
    EEC_w_h->GetXaxis()->SetLabelOffset(999);
    EEC_w_h->GetXaxis()->SetTitleOffset(999);
    EEC_w_h->GetXaxis()->SetTickLength(0);
    
    axis1->Draw();
    axis2->Draw();
    
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
    drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
   // drawText("40.7 Gev < jet p_{T1} >= 60.8 GeV", 0.20, 0.79, 18);
    //drawText("jet p_{T2} >= 9.4 GeV", 0.20, 0.75, 18);
    drawText("31.2 Gev < jet p_{T1} <= 40.7 GeV", 0.20, 0.79, 18);
    drawText("27.3 GeV <= jet p_{T2} < 31.2 GeV", 0.20, 0.75, 18);
    drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
    canvas2->SaveAs("WholeEventEEClowpthighApril17.png");
    
    canvas3->cd();
    canvas3->SetLogy();
    
    TH1F* histo3 = (TH1F*)histEECmid->Clone("histo3");
    
    for (int bin = 1; bin <= histo3->GetNbinsX(); ++bin) {
        float width = histo3->GetBinWidth(bin);
        float content = histo3->GetBinContent(bin);
        float error = histo3->GetBinContent(bin);
        histo3->SetBinContent(bin, content / width);
        histo3->SetBinError(bin, error);
    }
    
    TH1F* EEC_w_m = new TH1F("EEC_w_m", "WEEC", 60, 0, 60);
  
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_m->SetBinContent(bin, histo3->GetBinContent(bin));
        EEC_w_m->SetBinError(bin, histo3->GetBinError(bin));
    }
    
        
    EEC_w_m->Scale(1.0 / EEC_w_m->Integral("width"));
    
    EEC_w_m->SetStats(0);
    EEC_w_m->Draw("P E1");
    EEC_w_m->GetXaxis()->SetLabelOffset(999);
    EEC_w_m->GetXaxis()->SetTitleOffset(999);
    EEC_w_m->GetXaxis()->SetTickLength(0);
    
    axis1->Draw();
    axis2->Draw();
    
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
    drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
    //drawText("31.2 Gev < jet p_{T1} >= 40.7 GeV", 0.20, 0.79, 18);
    //drawText("jet p_{T2} >= 9.4 GeV", 0.20, 0.75, 18);
    drawText("31.2 Gev < jet p_{T1} <= 40.7 GeV", 0.20, 0.79, 18);
    drawText("20.9 GeV <= jet p_{T2} < 27.3 GeV", 0.20, 0.75, 18);
    drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
    canvas3->SaveAs("WholeEventEEClowptmidApril17.png");
    
    canvas4->cd();
    canvas4->SetLogy();
    
    TH1F* histo4 = (TH1F*)histEEClow->Clone("histo4");
    
    for (int bin = 1; bin <= histo4->GetNbinsX(); ++bin) {
        float width = histo4->GetBinWidth(bin);
        float content = histo4->GetBinContent(bin);
        float error = histo4->GetBinContent(bin);
        histo4->SetBinContent(bin, content / width);
        histo4->SetBinError(bin, error);
    }
    
    TH1F* EEC_w_l = new TH1F("EEC_w_l", "WEEC", 60, 0, 60);
 
    for (int bin = 1; bin <= 60; ++bin) {
        EEC_w_l->SetBinContent(bin, histo4->GetBinContent(bin));
        EEC_w_l->SetBinError(bin, histo4->GetBinError(bin));
    }
    
    EEC_w_l->Scale(1.0 / EEC_w_l->Integral("width"));
    
    EEC_w_l->SetStats(0);
    EEC_w_l->Draw("P E1");
    EEC_w_l->GetXaxis()->SetLabelOffset(999);
    EEC_w_l->GetXaxis()->SetTitleOffset(999);
    EEC_w_l->GetXaxis()->SetTickLength(0);
    
    axis1->Draw();
    axis2->Draw();
    
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
    drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
    //drawText("20.9 Gev < jet p_{T1} >= 31.2 GeV", 0.20, 0.79, 18);
    //drawText("jet p_{T2} >= 9.4 GeV", 0.20, 0.75, 18);
    drawText("31.2 Gev < jet p_{T1} <= 40.7 GeV", 0.20, 0.79, 18);
    drawText("9.4 Gev < jet p_{T2} < 20.9 GeV", 0.20, 0.75, 18);
    drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
    canvas4->SaveAs("WholeEventEEClowptlowApril17.png");
    
    canvas5->cd();
    canvas5->SetLogy();
    
    EEC_w_l ->SetLineColor(kRed);
    EEC_w_m->SetLineColor(kBlue);
    EEC_w_h->SetLineColor(kMagenta);

    EEC_w_l->SetMarkerStyle(20);
    EEC_w_l->SetMarkerColor(kRed);
    
    EEC_w_m->SetMarkerStyle(21);
    EEC_w_m->SetMarkerColor(kBlue);
    
    EEC_w_h->SetMarkerStyle(23);
    EEC_w_h->SetMarkerColor(kMagenta);

    EEC_w_l->SetMarkerSize(0.5);
    EEC_w_m->SetMarkerSize(0.5);
    EEC_w_h->SetMarkerSize(0.5);
    
    EEC_w_l->SetMinimum(1e-5);
    EEC_w_l->SetMaximum(1);
    
    EEC_w_l->Draw("P E1");
    EEC_w_m->Draw("P E1 SAME");
    EEC_w_h->Draw("P E1 SAME");
    
    axis1->Draw();
    axis2->Draw();
    
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
    drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
    drawText("31.2 Gev < jet p_{T1} >= 40.7 GeV", 0.20, 0.79, 18);
    drawText("charged constituent p_{T} > .2 GeV", 0.20, .75, 18 );
    
    TLegend* legend = new TLegend(0.6, 0.12, 0.9, 0.35);
    legend->AddEntry(EEC_w_l, "9.4 Gev < jet p_{T2} < 20.9 GeV", "pl");
    legend->AddEntry(EEC_w_m, "20.9 GeV <= jet p_{T2} < 27.3 GeV", "pl");
    legend->AddEntry(EEC_w_h, "27.3 GeV <= jet p_{T2} < 31.2 GeV", "pl");

    legend->Draw();

    canvas5->SaveAs("LowptdijetwholeeventcomparisonApril19.png");
}
