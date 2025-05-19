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

void simpleprints() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
  
    TFile* fin = TFile::Open("../outfiles/Thermaltest10kl.root");
    
    //TH2F* histetaphi = (TH2F*)fin->Get("hEtaPhi");
    //TH2F* histq2pmq2 = (TH2F*)fin->Get("comparisonpmq2");
    //TH2F* histqpmq = (TH2F*)fin->Get("comparisonpthat");
    //TH1F* histphi = (TH1F*)fin->Get("Deltaphi");
    //TH1F* histjetspec = (TH1F*)fin->Get("JetSpectrum");
    //TH1F* histleadingjetspec = (TH1F*)fin->Get("LeadingJetspectrum");
    //TH1F* histsubleadingjetspec = (TH1F*)fin->Get("SubleadingJetSpectrum");
    TH3F* histptetaphi = (TH3F*)fin->Get("ptvetaphi");
    
    TCanvas* canvas1 = new TCanvas("canvas1", "pt eta phi", 800, 600);
    TCanvas* canvas2 = new TCanvas("canvas2", "Comparison q2 vs q2 approximation", 800, 600);
    TCanvas* canvas3 = new TCanvas("canvas3", "Comparison q vs q approximation", 800, 600);
    TCanvas* canvas4 = new TCanvas("canvas4", "ThermalPhi", 800, 600);
    //TCanvas* canvas3 = new TCanvas("canvas3", "Jet Spectrum", 800, 600);
    //TCanvas* canvas4 = new TCanvas("canvas4", "Jet Spectrum leading", 800, 600);
    //TCanvas* canvas5 = new TCanvas("canvas5", "Jet Spectrum subleading", 800, 600);

    
    
    canvas1->cd();
    
    histptetaphi->SetStats(0);
    histptetaphi->GetXaxis()->SetRangeUser(0.0, 6.28);
    //histptetaphi->GetYaxis()->SetRangeUser(0.0, 6.28);
    histptetaphi->GetZaxis()->SetRangeUser(0.0, 6.28);
    //canvas1->SetLogz();
    histptetaphi->Draw("LEGO2");
    //drawText("1M events", 0.1, .9, 18 );
    canvas1->SaveAs("thermalptetaphi.png");
    
    /*
    canvas2->cd();
    
    histq2pmq2->SetStats(0);
    histq2pmq2->GetYaxis()->SetTitle("Q^{2} Gev^{2}");
    histq2pmq2->GetXaxis()->SetTitle("#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}^{2} GeV^{2}");
    histq2pmq2->GetXaxis()->SetTitleOffset(1.2);
    histq2pmq2->Draw();
    //drawText("1M events", 0.1, .9, 18 );
    
    canvas2->SaveAs("q2study.png");
    
    canvas3->cd();
    
    histqpmq->SetStats(0);
    histqpmq->GetYaxis()->SetTitle("#hat{p}_{T} GeV");
    histqpmq->GetXaxis()->SetTitle("#frac{jet p_{T, leading} + jet p_{T, subleading}}{2} GeV ");
    histqpmq->GetXaxis()->SetTitleOffset(1.2);
    canvas3->SetLogz();
    histqpmq->Draw("COLZ");
    //drawText("1M events", 0.1, .9, 18 );
    
    canvas3->SaveAs("qstudyMay13.png");
    
    canvas4->cd();
    
    histphi->SetStats(0);
    histphi-> Scale(1./ histphi->Integral("width"));
    histphi->GetYaxis()->SetRangeUser(0.0, 1);
    drawText("1M events", 0.1, .9, 18 );
    //histphi->GetYaxis()->SetRangeUser(0.0, .1);
    histphi->Draw();
    canvas4->SaveAs("thermalbackgroundphi.png");
    
    */
    }
