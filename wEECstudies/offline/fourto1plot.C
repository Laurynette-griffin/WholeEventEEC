#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
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

void fourto1plot(){
    
    TH1::SetDefaultSumw2();
    TFile* fin = TFile::Open("files/pythia_pp_rhic_fullevent_thermal_1mill_mixedeventstudy_centrality0-3_2025-05-14m.root");

    
    TH1F* histEECwhole = (TH1F*)fin->Get("EEC_w_phi");
    TH1F* histEEChigh = (TH1F*)fin->Get("EEC_w_hiphi"); // aj < .1
    TH1F* histEECmid = (TH1F*)fin->Get("EEC_w_midphi"); // aj < .667
    TH1F* histEEClow = (TH1F*)fin->Get("EEC_w_lowphi"); // aj < .15003
    
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
    //TCanvas* canvas5 = new TCanvas("canvas5", "Comparison", 800, 600);
    
    //axis stuff 
    TGaxis *axis1 = new TGaxis(0, 0, 30, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(60, 0, 30, 0, .00001 , .5, 510, "-G");
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
    
        TH1F* EEC_empty = new TH1F("EEC_empty", "WEEC", 60, 0, 60);
        EEC_empty->SetStats(0);
        EEC_empty->SetMinimum(1e-4);
        EEC_empty->SetMaximum(.5);
        EEC_empty->GetXaxis()->SetLabelOffset(999);
        EEC_empty->GetXaxis()->SetTitleOffset(999);
        EEC_empty->GetXaxis()->SetTickLength(0);
        EEC_empty->Draw("P E1");
        
        TH1F* EEC_w_p = new TH1F("EEC_w_p", "WEEC", 60, 0, 60);
        
        for (int bin = 1; bin <= 60; ++bin) {
            
            EEC_w_p->SetBinContent(bin, histo1->GetBinContent(bin));
            EEC_w_p->SetBinError(bin, histo1->GetBinError(bin));
        }
        
        
        TH1F* hrebintest = (TH1F*)EEC_w_p->Clone("hrebintest");
        hrebintest->Rebin(6);
        hrebintest->Scale(1./6);
        
        //hrebintest->SetMarkerStyle(20);
        //hrebintest->SetMarkerColor(kGreen+10);
        //hrebintest->SetLineColor(kGreen+10);
        //hrebintest->SetMarkerSize(2);
        hrebintest->GetXaxis()->SetRange(1,3);
        
        TH1F* hrebintest2 = (TH1F*)hrebintest->Clone("hrebintest2");
        hrebintest2->GetXaxis()->SetRange(8,10);
        EEC_w_p->SetMarkerSize(0.75);
        //EEC_w_p->GetXaxis()->SetRange(19,42);
        //hrebintest->Draw("p E1 same");
        EEC_w_p->Draw("P E1 same");
        //hrebintest2->Draw("p E1 same");
        
        
        EEC_w_p->GetXaxis()->SetLabelOffset(999);
        EEC_w_p->GetXaxis()->SetTitleOffset(999);
        EEC_w_p->GetXaxis()->SetTickLength(0);
        EEC_w_p->SetStats(0);
        EEC_w_p->GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");
        
        axis1->Draw();
        axis2->Draw();
        
        drawText("#pi/2", 0.4875, .06, 17 );
        //drawText("0.5", 0.4875, .06, 17 );
        //drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
        drawText("#Delta#pi", 0.8, .025, 18 );
        //drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
        drawText("Thermal Background", 0.20, 0.87, 18);
        drawText(".3 GeV< charged constituent p_{T} < 10 GeV", 0.20, .83, 18 );
        drawText("1M events", 0.20, 0.79, 18);
        //drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
        //drawText("20.9 Gev < jet p_{T1} >= 60.8 GeV", 0.20, 0.79, 18);
        //drawText("jet p_{T2} > 9.4 GeV", 0.20, 0.75, 18);
        //drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
        canvas1->SaveAs("Centrality03.png");
        

        canvas2->cd();
        canvas2->SetLogy();

        TH1F* histo2 = (TH1F*)histEEChigh->Clone("histo2");
        for (int bin = 1; bin <= histo2->GetNbinsX(); ++bin) {
            float width = histo2->GetBinWidth(bin);
            float content = histo2->GetBinContent(bin);
            float error = histo2->GetBinError(bin);
            histo2->SetBinContent(bin, content / width);
            histo2->SetBinError(bin, error);
        }
        
        histo2->Scale(1.0 / histo2->Integral("width"));
        
        TH1F* EEC_emptyh = new TH1F("EEC_emptyh", "WEEC", 60, 0, 60);
        EEC_emptyh->SetStats(0);
        EEC_emptyh->SetMinimum(1e-4);
        EEC_emptyh->SetMaximum(.5);
        EEC_emptyh->GetXaxis()->SetLabelOffset(999);
        EEC_emptyh->GetXaxis()->SetTitleOffset(999);
        EEC_emptyh->GetXaxis()->SetTickLength(0);
        EEC_emptyh->SetStats(0);
        EEC_emptyh->Draw("P E1");
        
        TH1F* EEC_w_h = new TH1F("EEC_w_h", "WEEC", 60, 0, 60);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_h->SetBinContent(bin, histo2->GetBinContent(bin));
            EEC_w_h->SetBinError(bin, histo2->GetBinError(bin));
        }
        

        
        TH1F * hrebintesth = (TH1F*)EEC_w_h->Clone("hrebintesth");
        hrebintesth->Rebin(6);
        hrebintesth->Scale(1./6);
    
        hrebintesth->GetXaxis()->SetRange(1,3);
        
        TH1F* hrebintest2h = (TH1F*)hrebintesth->Clone("hrebintest2h");
        hrebintest2h->GetXaxis()->SetRange(8,10);
        //EEC_w_h->GetXaxis()->SetRange(19,42);
        //hrebintesth->Draw("p E1 same");
        EEC_w_h->Draw("P E1 same");
        //hrebintest2h->Draw("p E1 same");
        
        EEC_w_h->SetStats(0);
        //EEC_w_h->Draw("P E1");
        EEC_w_h->GetXaxis()->SetLabelOffset(999);
        EEC_w_h->GetXaxis()->SetTitleOffset(999);
        EEC_w_h->GetXaxis()->SetTickLength(0);
        EEC_w_h->GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");
        
        axis1->Draw();
        axis2->Draw();
        
        drawText("0.5", 0.4875, .06, 17 );
        drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
        drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
        drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
        drawText("40.7 Gev < jet p_{T1} <= 60.8 GeV", 0.20, 0.79, 18);
        drawText("jet p_{T2} > 9.4 GeV", 0.20, 0.75, 18);
        drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
        canvas2->SaveAs("Centrality1015.png");
        
  
        canvas3->cd();
        canvas3->SetLogy();
        
        TH1F* histo3 = (TH1F*)histEECmid->Clone("histo3");
        
        for (int bin = 1; bin <= histo3->GetNbinsX(); ++bin) {
            float width = histo3->GetBinWidth(bin);
            float content = histo3->GetBinContent(bin);
            float error = histo3->GetBinError(bin);
            histo3->SetBinContent(bin, content / width);
            histo3->SetBinError(bin, error);
        }
        
        histo3->Scale(1.0 / histo3->Integral("width"));
        
        TH1F* EEC_emptym = new TH1F("EEC_w_m", "WEEC", 60, 0, 60);
        EEC_emptym->SetStats(0);
        EEC_emptym->SetMinimum(1e-4);
        EEC_emptym->SetMaximum(.5);
        EEC_emptym->GetXaxis()->SetLabelOffset(999);
        EEC_emptym->GetXaxis()->SetTitleOffset(999);
        EEC_emptym->GetXaxis()->SetTickLength(0);
        EEC_emptym->SetStats(0);
        EEC_emptym->Draw("P E1");
        
        TH1F* EEC_w_m = new TH1F("EEC_w_m", "WEEC", 60, 0, 60);
      
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_m->SetBinContent(bin, histo3->GetBinContent(bin));
            EEC_w_m->SetBinError(bin, histo3->GetBinError(bin));
        }

        TH1F * hrebintestm = (TH1F*)EEC_w_m->Clone("hrebintesth");
        hrebintestm->Rebin(6);
        hrebintestm->Scale(1./6);
    
        hrebintestm->GetXaxis()->SetRange(1,3);
        
        TH1F* hrebintest2m = (TH1F*)hrebintestm->Clone("hrebintest2m");
        hrebintest2m->GetXaxis()->SetRange(8,10);
        //EEC_w_m->GetXaxis()->SetRange(19,42);
        //hrebintestm->Draw("p E1 same");
        EEC_w_m->Draw("P E1 SAME");
        //hrebintest2m->Draw("p E1 same");
        
        EEC_w_m->SetStats(0);
        EEC_w_m->Draw("P E1");
        EEC_w_m->GetXaxis()->SetLabelOffset(999);
        EEC_w_m->GetXaxis()->SetTitleOffset(999);
        EEC_w_m->GetXaxis()->SetTickLength(0);
        EEC_w_m->GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");
        
        axis1->Draw();
        axis2->Draw();
        
        drawText("0.5", 0.4875, .06, 17 );
        drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
        drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
        drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
        drawText("31.2 Gev < jet p_{T1} <= 40.7 GeV", 0.20, 0.79, 18);
        drawText("jet p_{T2} > 9.4 GeV", 0.20, 0.75, 18);
        drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
        canvas3->SaveAs("Centrality2530.png");

        canvas4->cd();
        canvas4->SetLogy();
        
        TH1F* histo4 = (TH1F*)histEEClow->Clone("histo4");
        
        for (int bin = 1; bin <= histo4->GetNbinsX(); ++bin) {
            float width = histo4->GetBinWidth(bin);
            float content = histo4->GetBinContent(bin);
            float error = histo4->GetBinError(bin);
            histo4->SetBinContent(bin, content / width);
            histo4->SetBinError(bin, error / width);
        }
        
        histo4->Scale(1.0 / histo4->Integral("width"));
        
        TH1F* EEC_emptyl = new TH1F("EEC_w_m", "WEEC", 60, 0, 60);
        EEC_emptyl->SetStats(0);
        EEC_emptyl->SetMinimum(1e-4);
        EEC_emptyl->SetMaximum(.5);
        EEC_emptyl->GetXaxis()->SetLabelOffset(999);
        EEC_emptyl->GetXaxis()->SetTitleOffset(999);
        EEC_emptyl->GetXaxis()->SetTickLength(0);
        EEC_emptyl->SetStats(0);
        EEC_emptyl->Draw("P E1");
        
        TH1F* EEC_w_l = new TH1F("EEC_w_l", "WEEC", 60, 0, 60);
     
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_l->SetBinContent(bin, histo4->GetBinContent(bin));
            EEC_w_l->SetBinError(bin, histo4->GetBinError(bin));
        }
        
     
        
        TH1F * hrebintestl = (TH1F*)EEC_w_l->Clone("hrebintestl");
        hrebintestl->Rebin(6);
        hrebintestl->Scale(1./6);
    
        hrebintestl->GetXaxis()->SetRange(1,3);
        
        TH1F* hrebintest2l = (TH1F*)hrebintestl->Clone("hrebintest2l");
        hrebintest2l->GetXaxis()->SetRange(8,10);
        //EEC_w_l->GetXaxis()->SetRange(19,42);
        //hrebintestl->Draw("p E1 same");
        EEC_w_l->Draw("P E1 SAME");
        //hrebintest2l->Draw("p E1 same");
        
        
        EEC_w_l->SetStats(0);
        EEC_w_l->Draw("P E1");
        EEC_w_l->GetXaxis()->SetLabelOffset(999);
        EEC_w_l->GetXaxis()->SetTitleOffset(999);
        EEC_w_l->GetXaxis()->SetTickLength(0);
        EEC_w_l->GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");
        
        axis1->Draw();
        axis2->Draw();
        
        drawText("0.5", 0.4875, .06, 17 );
        drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18 );
        drawText("p+p #srqt{s} = 200 GeV", 0.20, 0.87, 18);
        drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.83, 18);
        drawText("20.9 Gev < jet p_{T1} <= 31.2 GeV", 0.20, 0.79, 18);
        drawText("jet p_{T2} > 9.4 GeV", 0.20, 0.75, 18);
        drawText("charged constituent p_{T} > .2 GeV", 0.20, .71, 18 );
        canvas4->SaveAs("Centrality4045.png");

        canvas5->cd();
        canvas5->SetLogy();

        TH1F* EEC_empty5 = new TH1F("EEC_empty5", "WEEC", 60, 0, 60);
        EEC_empty5->GetXaxis()->SetTickLength(0);
        EEC_empty5->SetMinimum(1e-1);
        EEC_empty5->SetMaximum(1);
        EEC_empty5->SetStats(0);
        EEC_empty5->Draw("P E1");
  
        EEC_w_h->SetLineColor(kGreen+2);
        EEC_w_h->SetMarkerStyle(24);
        EEC_w_h->SetMarkerColor(kGreen+2);
        EEC_w_h->SetMarkerSize(0.5);
  

        EEC_w_l->SetLineColor(kRed);

        EEC_w_m->SetLineColor(kBlue);
        
        EEC_w_h->SetLineColor(kMagenta);

        EEC_w_l->SetMarkerStyle(20);
        EEC_w_l->SetMarkerColor(kRed);
   
        EEC_w_m->SetMarkerStyle(21);
        EEC_w_m->SetMarkerColor(kBlue);
        
        EEC_w_h->SetMarkerStyle(25);
        EEC_w_h->SetMarkerColor(kMagenta);
       
        
        EEC_w_l->SetMarkerSize(0.75);
        EEC_w_m->SetMarkerSize(0.75);
        EEC_w_h->SetMarkerSize(0.75);
        
        EEC_w_l->Draw("P E1 SAME");
        EEC_w_m->Draw("P E1 SAME");
        EEC_w_h->Draw("P E1 SAME");
      
        cout <<"seg" << endl;
        axis1->Draw();
        axis2->Draw();
        
        drawText("0.5", 0.4875, .06, 17 );
        drawText("z = (1 - cos(#theta))/2", 0.8, .025, 16);
        //drawText("#Delta#Phi", 0.8, .025, 16);
        drawText("p+p #sqrt{s} = 200 GeV", 0.20, 0.88, 16);
        drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.20, 0.84, 16);
        drawText("31.2 Gev < jet p_{T1} <= 40.7 GeV", 0.20, 0.80, 16);
        drawText("27.3 Gev < jet p_{T2} <= 31.2 GeV", 0.20, 0.76, 16);
        drawText("charged constituent p_{T} > .2 GeV", 0.20, .72, 16 );
        drawText("#Delta#Phi > 3/4 ", 0.20, .68, 16 );
        
        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.88);
        legend->AddEntry(EEC_w_p, "", "pl");
        legend->AddEntry(EEC_w_l, "Signal + Signal", "pl");
        legend->AddEntry(EEC_w_m, "Thermal + Signal", "pl");
        legend->AddEntry(EEC_w_h, "Thermal + Thermal", "pl");
        
        legend->Draw();
        cout <<"seg" << endl;
        canvas5->SaveAs("centralitystudyWEECphiMay16.png");
        
    //delete EEC_w_p;
    delete EEC_w_h;
    delete EEC_w_m;
    delete EEC_w_l;
    
    }
