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

void plotmaker_eec8() {

 TH1::SetDefaultSumw2();
  TFile* fin = TFile::Open("pythia_pp_rhic_8GeV_fullevent_5pthat60_Apr9_dr_thirdjet.root"); // 8 GeV q2
  
  
  TH1F* histEECwhole = (TH1F*)fin->Get("EEC_w");
  
  TH1F* histEECb = (TH1F*)fin->Get("EEC_w_b"); // aj < .1
  TH1F* histEECb2 = (TH1F*)fin->Get("EEC_w_b2"); // aj < .667
  TH1F* histEECb3 = (TH1F*)fin->Get("EEC_w_b3"); // aj < .15003
  TH1F* histEECub3 = (TH1F*)fin->Get("EEC_w_ub3"); // aj > .15003
  TH1F* histEECub2 = (TH1F*)fin->Get("EEC_w_ub2"); // aj > .1
  TH1F* histEECub = (TH1F*)fin->Get("EEC_w_ub"); // aj > .667
  
  
  TH1F* histEEC = (TH1F*)fin->Get("EEC_star");
  TH1F* histAjspec = (TH1F*)fin->Get("Aj_spectrum");
  TH1F* histcosspec = (TH1F*)fin->Get("Costheta_spectrum");
  TH1F* histjetspec = (TH1F*)fin->Get("JetSpectrum");
  TH2F* histetaphi = (TH2F*)fin->Get("etaphi_spectrum");

  TCanvas* canvas1 = new TCanvas("canvas1", "Whole event EEC", 800, 600);
  TCanvas* canvas2 = new TCanvas("canvas2", "Cosine Theta Spectrum", 800, 600);
  TCanvas* canvas3 = new TCanvas("canvas3", "Aj spectrum", 800, 600);
  TCanvas* canvas4 = new TCanvas("canvas4", "Whole event EEC balanced", 800, 600);
  TCanvas* canvas5 = new TCanvas("canvas5", "Whole event EEC less balanced", 800, 600);
  TCanvas* canvas6 = new TCanvas("canvas6", "Whole event EEC unbalanced", 800, 600);
  TCanvas* canvas7 = new TCanvas("canvas7", "Whole event EEC less unbalanced", 800, 600);
  TCanvas* canvas8 = new TCanvas("canvas8", "Whole event EEC least balanced", 800, 600);
  TCanvas* canvas9 = new TCanvas("canvas9", "Whole event EEC most unbalanced", 800, 600);
  TCanvas* canvas0 = new TCanvas("canvas0", "Jet Spectrum", 800, 600);
  TCanvas* canvas10 = new TCanvas("canvas10", "EtaPhi spectrum", 800, 600);
  TCanvas* canvas11 = new TCanvas("canvas11", "#Deltar spectrum", 800, 600);
  
  
  
  
 // make this a function 
  canvas1->cd();
  canvas1->SetLogy();
  for (int i = 0; i <= histEECwhole->GetNbinsX(); ++i) {
      double binWidth = histEECwhole->GetBinWidth(i);
      double binContent = histEECwhole->GetBinContent(i);
      histEECwhole->SetBinContent(i, binContent / binWidth);
  }
  
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
  EEC_w_p.GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");
  
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
  drawText("z = (1 - cos(#Deltar))/2", 0.8, .025, 18 );
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("charged constituent p_{T} > .2 GeV", 0.20, .69, 18 );
  
  canvas1->SaveAs("WholeEventEEC8.png");
    
  // Draw second histogram
  canvas2->cd();
  //histcosspec->Scale(1.0 / histAjspec->Integral());
  histcosspec->SetStats(0);
  histcosspec->GetXaxis()->SetTitle("cos(#Deltar)");
  histcosspec->Draw("P E1");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.77, 18);
  
  canvas0->SaveAs("CosSpectrum8.png");
 
  // Draw third histogram     
  canvas3->cd();
  histAjspec->Scale(1.0 / histAjspec->Integral());
  histAjspec->SetStats(0);
  histAjspec->GetXaxis()->SetTitle("(jet p_{T{0}} - jet p_{T{1}})/(jet p_{T{0}} + jet p_{T{1}})");
  histAjspec->Draw("P E0");
  
 
  canvas3->SaveAs("AjSpectrum8.png");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.37, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.25, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.21, 18);
  
  canvas4->cd();
  canvas4->SetLogy();
  for (int i = 1; i <= histEECb->GetNbinsX(); ++i) {
      double binWidth = histEECb->GetBinWidth(i);
      double binContent = histEECb->GetBinContent(i);
      histEECb->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pb("EEC_w_pb", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 1; i <= histEECb->GetNbinsX(); ++i){
    EEC_w_pb.SetBinContent(i, histEECb->GetBinContent(i));
    EEC_w_pb.SetBinError(i, histEECb->GetBinError(i));
  }
        
  EEC_w_pb.SetStats(0);
  EEC_w_pb.GetXaxis()->SetTickLength(0);
  EEC_w_pb.GetXaxis()->SetLabelOffset(999);
  EEC_w_pb.GetXaxis()->SetTitleOffset(999);
  EEC_w_pb.Draw("P E1");

  axis1->Draw();
  axis2->Draw();
  
  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#Deltar))/2", 0.8, .025, 18 );
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} < 0.0667", 0.20, 0.69, 18);
  drawText("charged constituent p_{T} > .5 GeV", 0.20, .65, 18 );
  
  canvas4->SaveAs("WholeEventEEC8dtp.png");
    
  canvas5->cd();
  canvas5->SetLogy();
  for (int i = 1; i <= histEECb2->GetNbinsX(); ++i) {
      double binWidth = histEECb2->GetBinWidth(i);
      double binContent = histEECb2->GetBinContent(i);
      histEECb2->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pb2("EEC_w_pb2", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 1; i <= histEECb2->GetNbinsX(); ++i){
    EEC_w_pb2.SetBinContent(i, histEECb2->GetBinContent(i));
    EEC_w_pb2.SetBinError(i, histEECb2->GetBinError(i));
  }
        
        
  EEC_w_pb2.SetStats(0);
  EEC_w_pb2.GetXaxis()->SetTickLength(0);
  EEC_w_pb2.GetXaxis()->SetLabelOffset(999);
  EEC_w_pb2.GetXaxis()->SetTitleOffset(999);
  EEC_w_pb2.Draw("P E1");
  
  axis1->Draw();
  axis2->Draw();

  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} < 0.1", 0.20, 0.69, 18);
  drawText("charged constituent p_{T} > .2 GeV", 0.20, .65, 18 );
  
  canvas5->SaveAs("WholeEventEEClessbalanced8.png");
  /*
  canvas6->cd();
  canvas6->SetLogy();
  for (int i = 1; i <= histEECub2->GetNbinsX(); ++i) {
      double binWidth = histEECub2->GetBinWidth(i);
      double binContent = histEECub2->GetBinContent(i);
      histEECub2->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pub2("EEC_w_pub2", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 1; i <= histEECub2->GetNbinsX(); ++i){
    EEC_w_pub2.SetBinContent(i, histEECub2->GetBinContent(i));
    EEC_w_pub2.SetBinError(i, histEECub2->GetBinError(i));
  }
        
  EEC_w_pub2.SetStats(0);
  EEC_w_pub2.GetXaxis()->SetLabelOffset(999);
  EEC_w_pub2.GetXaxis()->SetTitleOffset(999);
  EEC_w_pub2.GetXaxis()->SetTickLength(0);
  EEC_w_pub2.Draw("P E1");
  
  axis1->Draw();
  axis2->Draw(); 

  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );

  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} > 0.1", 0.20, 0.69, 18);
  drawText("charged constituent p_{T} > .2 GeV", 0.20, .65, 18 );
  
  canvas6->SaveAs("WholeEventEECunbalanced8.png");
  
  canvas7->cd();
  canvas7->SetLogy();
  for (int i = 1; i <= histEECub->GetNbinsX(); ++i) {
      double binWidth = histEECub->GetBinWidth(i);
      double binContent = histEECub->GetBinContent(i);
      histEECub->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pub("EEC_w_pub", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 1; i <= histEECub->GetNbinsX(); ++i){
    EEC_w_pub.SetBinContent(i, histEECub->GetBinContent(i));
    EEC_w_pub.SetBinError(i, histEECub->GetBinError(i));
  }
        
  EEC_w_pub.SetStats(0);
  EEC_w_pub.GetXaxis()->SetLabelOffset(999);
  EEC_w_pub.GetXaxis()->SetTitleOffset(999);
  EEC_w_pub.GetXaxis()->SetTickLength(0);
  EEC_w_pub.Draw("P E1");
  
  axis1->Draw();
  axis2->Draw(); 

  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} > 0.0667", 0.20, 0.69, 18);
  drawText("charged constituent p_{T} > .2 GeV", 0.20, .65, 18 );
  
  canvas7->SaveAs("WholeEventEEClessunbalanced8.png");
  
  canvas8->cd();
  canvas8->SetLogy();
  for (int i = 1; i <= histEECub3->GetNbinsX(); ++i) {
      double binWidth = histEECub3->GetBinWidth(i);
      double binContent = histEECub3->GetBinContent(i);
      histEECub3->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pub3("EEC_w_pub3", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 1; i <= histEECub3->GetNbinsX(); ++i){
    EEC_w_pub3.SetBinContent(i, histEECub3->GetBinContent(i));
    EEC_w_pub3.SetBinError(i, histEECub3->GetBinError(i));
  }
        
  EEC_w_pub3.SetStats(0);
  EEC_w_pub3.GetXaxis()->SetLabelOffset(999);
  EEC_w_pub3.GetXaxis()->SetTitleOffset(999);
  EEC_w_pub3.GetXaxis()->SetTickLength(0);
  EEC_w_pub3.Draw("P E1");
  
  axis1->Draw();
  axis2->Draw(); 

  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} > 0.15003", 0.20, 0.69, 18); 
  drawText("charged constituent p_{T} > .2 GeV", 0.20, .65, 18 );
  
  canvas8->SaveAs("WholeEventEECmostunbalanced8.png");
  
  canvas9->cd();
  canvas9->SetLogy();
  for (int i = 1; i <= histEECb3->GetNbinsX(); ++i) {
      double binWidth = histEECb3->GetBinWidth(i);
      double binContent = histEECb3->GetBinContent(i);
      histEECb3->SetBinContent(i, binContent / binWidth);
  }
  
  TH1F EEC_w_pb3("EEC_w_pb3", "Energy Energy Correlator", 60, 0, 60);

  for(int i = 1; i <= histEECb3->GetNbinsX(); ++i){
    EEC_w_pb3.SetBinContent(i, histEECb3->GetBinContent(i));
    EEC_w_pb3.SetBinError(i, histEECb3->GetBinError(i));
  }
        
  EEC_w_pb3.SetStats(0);
  EEC_w_pb3.GetXaxis()->SetLabelOffset(999);
  EEC_w_pb3.GetXaxis()->SetTitleOffset(999);
  EEC_w_pb3.GetXaxis()->SetTickLength(0);
  EEC_w_pb3.Draw("P E1");
  
  axis1->Draw();
  axis2->Draw(); 

  drawText("0.5", 0.4875, .06, 17 );
  drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.85, 18); 
  drawText("anti-k_{T} R_{jet} = 0.4 |#eta_{jet} < 0.7|", 0.20, 0.81, 18);
  drawText("jet p_{T, leading} > 8 GeV", 0.20, 0.77, 18);
  drawText("jet p_{T, subleading} > 8 GeV", 0.20, 0.73, 18);
  drawText("jet A_{j} < 0.15003", 0.20, 0.69, 18);
  drawText("charged constituent p_{T} > .2 GeV", 0.20, .65, 18 );
  
  canvas9->SaveAs("WholeEventEECleastbalanced8.png");
  */
  canvas0->cd();
  histjetspec->Scale(1.0 / histjetspec->Integral());
  histjetspec->SetStats(0);
  histjetspec->GetXaxis()->SetTitle("p{T}");
  histjetspec->Draw("P E1");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.77, 18);
  
  canvas0->SaveAs("JetSpectrum8.png");
  
  canvas10->cd();
  histetaphi->SetStats(0);
  histetaphi->GetXaxis()->SetTitle("#phi");
  histetaphi->GetYaxis()->SetTitle("#eta");
  histetaphi->Draw("COLZ");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.77, 18);
  
  canvas10->SaveAs("EtaPhi8.png");
  
  canvas11->cd();
  histcosspec->SetStats(0);
  histcosspec->GetXaxis()->SetTitle("cos(#theta)");
  histcosspec->Draw("P E1");
  
  drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.77, 18);
  
  canvas2->SaveAs("thetadistribution8.png");
}
