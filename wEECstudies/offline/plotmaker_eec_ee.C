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

void plotmaker_eec_ee() {

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  //TFile* fin = TFile::Open("../offline/pythia_wholeeventEEC_ee_March5.root"); // e+e-
  TFile* fin = TFile::Open("../outfiles/pythia_ee_aleph_fullevent_June03_50m.root"); 

  TH1F* histEECwhole = (TH1F*)fin->Get("EEC_w");

  TCanvas* canvas1 = new TCanvas("canvas1", "Whole event EEC", 800, 600);

  
 // make this a function 
  canvas1->cd();
  canvas1->SetLogy();
  for (int i = 0; i <= histEECwhole->GetNbinsX(); ++i) {
      double binWidth = histEECwhole->GetBinWidth(i);
      double binContent = histEECwhole->GetBinContent(i);
      histEECwhole->SetBinContent(i, binContent / binWidth);
  }
  histEECwhole->Scale(1.0/ histEECwhole->Integral("width"));
  
  TH1F EEC_w_p("EEC_w_p", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 0; i <= histEECwhole->GetNbinsX(); ++i){
      
    EEC_w_p.SetBinContent(i, histEECwhole->GetBinContent(i));
    EEC_w_p.SetBinError(i, histEECwhole->GetBinError(i));
  }
        
  EEC_w_p.SetStats(0);
  EEC_w_p.Draw("P E1");
  
  EEC_w_p.GetXaxis()->SetLabelOffset(999);
  EEC_w_p.GetXaxis()->SetTitleOffset(999);
  EEC_w_p.GetXaxis()->SetTickLength(0);
  EEC_w_p.GetYaxis()->SetTitle("");
  
  //thankfully this only needs to be done once since all of the WEEEC have the same x axis 
  TGaxis *axis1 = new TGaxis(0, 0, 30, 0, .00001, .5, 510, "G"); 
  TGaxis *axis2 = new TGaxis(60, 0, 30, 0, .00001 , .5 , 510, "-G");
  axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" ); //turn off
  axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" ); //turn off
  axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" ); //turn off
  axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" ); //turn off
  axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
  axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" ); //turn off
  axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
  axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" ); //turn off
  axis2->SetLabelOffset(0.045);
  axis1->SetLabelSize(.03);
  axis2->SetLabelSize(.03);
  axis1->SetTickSize(0.05);
  axis2->SetTickSize(0.05);
  axis1->Draw();
  axis2->Draw(); 
  
  drawText("0.5", 0.4875, .06, 17 );
  //drawText("#pi/2", 0.505, .065, 17 );
  //drawText("#Delta#phi", 0.8, .025, 18);
  drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
  drawText("e+e- #srqt{s} = 91.2 GeV", 0.20, 0.85, 18); 
  drawText("event p_{T} > 15 GeV", 0.20, 0.81, 18);
  drawText("charged particle p_{T} > .2 GeV", 0.20, .77, 18 );
  
  canvas1->SaveAs("WholeEventeeJune03.png");


}
