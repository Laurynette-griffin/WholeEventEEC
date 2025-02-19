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

void plotmaker_eec() {

 TH1::SetDefaultSumw2();
  TFile* fin = TFile::Open("pythia_pp_rhic_10pthat60_eec_star+fullevent_feb5.root");

  //TFile* fin = TFile::Open("pythia_wholeeventEEC_starEEC_Feb17.out");  

  TH1F* histEECwhole = (TH1F*)fin->Get("EEC_w");
  
  TH1F* histEECb = (TH1F*)fin->Get("EEC_w_b"); // aj <.167
  TH1F* histEECb2 = (TH1F*)fin->Get("EEC_w_b2"); // aj < .333
  TH1F* histEECub2 = (TH1F*)fin->Get("EEC_w_ub2"); // aj > .167
  TH1F* histEECub = (TH1F*)fin->Get("EEC_w_ub"); // aj > .333
  
  TH1F* histEEC = (TH1F*)fin->Get("EEC_star");
  TH1F* histAjspec = (TH1F*)fin->Get("Aj_spectrum");


  TCanvas* canvas1 = new TCanvas("canvas1", "Whole event EEC", 800, 600);
  TCanvas* canvas2 = new TCanvas("canvas2", "EEC", 800, 600);
  TCanvas* canvas3 = new TCanvas("canvas3", "Aj spectrum", 800, 600);
  TCanvas* canvas4 = new TCanvas("canvas4", "Whole event EEC balanced", 800, 600);
  TCanvas* canvas5 = new TCanvas("canvas5", "Whole event EEC less balanced", 800, 600);
  TCanvas* canvas6 = new TCanvas("canvas6", "Whole event EEC unbalanced", 800, 600);
  TCanvas* canvas7 = new TCanvas("canvas7", "Whole event EEC less unbalanced", 800, 600);
  
 // make this a function 
  canvas1->cd();
  canvas1->SetLogy();
  for (int i = 0; i <= histEECwhole->GetNbinsX(); ++i) {
      double binWidth = histEECwhole->GetBinWidth(i);
      double binContent = histEECwhole->GetBinContent(i);
      histEECwhole->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_p("EEC_w_p", "Energy Energy Correlator", 50, 0, 50);

  for(int i = 0; i <= histEECwhole->GetNbinsX(); ++i){
      
    EEC_w_p.SetBinContent(i, histEECwhole->GetBinContent(i));
    EEC_w_p.SetBinError(i, histEECwhole->GetBinError(i));
  }
        
  EEC_w_p.SetStats(0);
  EEC_w_p.Draw("P E1");
  
  EEC_w_p.GetXaxis()->SetLabelOffset(999);
  EEC_w_p.GetXaxis()->SetTitleOffset(999);
  EEC_w_p.GetXaxis()->SetTickLength(0);

  TGaxis *axis1 = new TGaxis(0, 0, 25, 0, .0001, .5, 510, "G"); 
  TGaxis *axis2 = new TGaxis(50, 0, 25, 0, .0001 , .5 , 510, "-G");
  //axis2->SetTitle("z = (1 - cos(#Deltar))/2");
  axis1->ChangeLabel(2, -1, 0, -1, -1, 62, "1 - 10^{-4}" );
  axis1->ChangeLabel(4, -1, 0, -1, -1, 62, "1 - 10^{-4}" );
  //cout << "axis font" << axis2->GetLabelFont() << endl; 
  axis2->ChangeLabel(1, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
  axis2->ChangeLabel(2, -1, 0, -1, -1, 62, "1 - 10^{-3}" );
  axis2->ChangeLabel(3, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
  axis2->ChangeLabel(4, -1, 0, -1, -1, 62, "1 - 10^{-1}" );
  axis2->SetLabelOffset(0.045);
  axis1->SetLabelSize(.03);
  axis2->SetLabelSize(.03);
  axis1->SetTickSize(0.05);
  axis2->SetTickSize(0.05);
  axis1->Draw();
  axis2->Draw(); 
  
  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#Deltar))/2", 0.8, .025, 18 );
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 15 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("constituent p_{T} > .5 GeV", 0.20, .69, 18 );
  drawText("z = (1 - cos(#Deltar))/2", 0.8, .025, 18 );
  
  canvas1->SaveAs("WholeEventEEC.png");
    
  // Draw second histogram
  canvas2->cd();
  canvas2->SetLogx();
  canvas2->SetLogy();
  histEEC->Scale(1.0 / histEEC->Integral());
  histEEC->GetXaxis()->SetRangeUser(0.015, 1);
  histEEC->SetStats(0);
  histEEC->Draw("P E0");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.37, 18);
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.6|", 0.20, 0.33, 18);
  drawText("constituent p_{T} > 2 GeV", 0.20, 0.29, 18);;
  drawText("15 < jet p_{T} < 20", 0.20, 0.25, 18);

  canvas2->SaveAs("EEC.png");
 
  // Draw third histogram     
  canvas3->cd();
  histAjspec->Scale(1.0 / histAjspec->Integral());
  histAjspec->SetStats(0);
  histAjspec->Draw("P E0");
 
  canvas3->SaveAs("AjSpectrum.png");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.37, 18);
  drawText("jet p_{T, leading} > 15 GeV", 0.20, 0.25, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.21, 18);
  
  canvas4->cd();
  canvas4->SetLogy();
  for (int i = 1; i <= histEECb->GetNbinsX(); ++i) {
      double binWidth = histEECb->GetBinWidth(i);
      double binContent = histEECb->GetBinContent(i);
      histEECb->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pb("EEC_w_pb", "Energy Energy Correlator", 50, 0, 50);

  for(int i = 1; i <= histEECb->GetNbinsX(); ++i){
    EEC_w_pb.SetBinContent(i, histEECb->GetBinContent(i));
    EEC_w_pb.SetBinError(i, histEECb->GetBinError(i));
  }
        
  EEC_w_pb.SetStats(0);
  EEC_w_pb.Draw("P E0");
  
 // EEC_w_pb.GetXaxis()->SetLabelOffset(999);
 // EEC_w_pb.GetXaxis()->SetTitleOffset(999);
 
  TGaxis *axis3 = new TGaxis(0, 10e-4, 25, 10e-4, .0001, .5, 510, "G");
  TGaxis *axis4 = new TGaxis(50, 10e-4, 25, 10e-4, .0001, 1, 510, "-G");
  axis4->SetTitle("z = 1/(cos(#Deltar))");
  axis3->Draw();
  axis4->Draw(); 
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 15 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} < 0.667", 0.20, 0.69, 18);
  drawText("constituent p_{T} > .5 GeV", 0.20, .65, 18 );
  
  canvas4->SaveAs("WholeEventEECbalanced.png");
    
  canvas5->cd();
  canvas5->SetLogy();
  for (int i = 1; i <= histEECb2->GetNbinsX(); ++i) {
      double binWidth = histEECb2->GetBinWidth(i);
      double binContent = histEECb2->GetBinContent(i);
      histEECb2->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pb2("EEC_w_pb2", "Energy Energy Correlator", 50, 0, 50);

  for(int i = 1; i <= histEECb2->GetNbinsX(); ++i){
    EEC_w_pb2.SetBinContent(i, histEECb2->GetBinContent(i));
    EEC_w_pb2.SetBinError(i, histEECb2->GetBinError(i));
  }
        
  EEC_w_pb2.SetStats(0);
  EEC_w_pb2.Draw("P E0");
  
 // EEC_w_pb2.GetXaxis()->SetLabelOffset(999);
 // EEC_w_pb2.GetXaxis()->SetTitleOffset(999);
  
  TGaxis *axis5 = new TGaxis(0, 10e-4, 25, 10e-4, .0001, .5, 510, "G");
  TGaxis *axis6 = new TGaxis(50, 10e-4, 25, 10e-4, .0001, 1, 510, "-G");
  axis6->SetTitle("z = 1/(cos(#Deltar))");
  axis5->Draw();
  axis6->Draw(); 
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 15 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} < 0.1", 0.20, 0.69, 18);
  drawText("constituent p_{T} > .5 GeV", 0.20, .65, 18 );
  
  canvas5->SaveAs("WholeEventEEClessbalanced.png");
  
  canvas6->cd();
  canvas6->SetLogy();
  for (int i = 1; i <= histEECub2->GetNbinsX(); ++i) {
      double binWidth = histEECub2->GetBinWidth(i);
      double binContent = histEECub2->GetBinContent(i);
      histEECub2->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pub2("EEC_w_pb2", "Energy Energy Correlator", 50, 0, 50);

  for(int i = 1; i <= histEECub2->GetNbinsX(); ++i){
    EEC_w_pub2.SetBinContent(i, histEECub2->GetBinContent(i));
    EEC_w_pub2.SetBinError(i, histEECub2->GetBinError(i));
  }
        
  EEC_w_pub2.SetStats(0);
  EEC_w_pub2.Draw("P E0");
  
 // EEC_w_pub2.GetXaxis()->SetLabelOffset(999);
 // EEC_w_pub2.GetXaxis()->SetTitleOffset(999);
  
  TGaxis *axis7 = new TGaxis(0, 10e-4, 25, 10e-4, .0001, .5, 510, "G");
  TGaxis *axis8 = new TGaxis(50, 10e-4, 25, 10e-4, .0001, 1, 510, "-G");
  axis8->SetTitle("z = 1/(cos(#Deltar))");
  axis7->Draw();
  axis8->Draw(); 
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 15 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} > 0.1", 0.20, 0.69, 18);
  drawText("constituent p_{T} > .5 GeV", 0.20, .65, 18 );
  
  canvas6->SaveAs("WholeEventEECunbalanced.png");
  
  canvas7->cd();
  canvas7->SetLogy();
  for (int i = 1; i <= histEECub->GetNbinsX(); ++i) {
      double binWidth = histEECub->GetBinWidth(i);
      double binContent = histEECub->GetBinContent(i);
      histEECub->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pub("EEC_w_pb2", "Energy Energy Correlator", 50, 0, 50);

  for(int i = 1; i <= histEECub->GetNbinsX(); ++i){
    EEC_w_pub.SetBinContent(i, histEECub->GetBinContent(i));
    EEC_w_pub.SetBinError(i, histEECub->GetBinError(i));
  }
        
  EEC_w_pub.SetStats(0);
  EEC_w_pub.Draw("P E0");
  
 // EEC_w_pub.GetXaxis()->SetLabelOffset(999);
 // EEC_w_pub.GetXaxis()->SetTitleOffset(999);
  
  TGaxis *axis9 = new TGaxis(0, 10e-4, 25, 10e-4, .0001, .5, 510, "G");
  TGaxis *axis0 = new TGaxis(50, 10e-4, 25, 10e-4, .0001, 1, 510, "-G");
  axis0->SetTitle("z = 1/(cos(#Deltar))");
  axis9->Draw();
  axis0->Draw(); 
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 15 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} > 0.667", 0.20, 0.69, 18);
  drawText("constituent p_{T} > .5 GeV", 0.20, .65, 18 );
  
  canvas7->SaveAs("WholeEventEEClessunbalanced.png");
}
