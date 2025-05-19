#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <iostream>

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}

void plotmaker4to1() {
    TH1::SetDefaultSumw2();

    const int number = 4;
    
    const char* filenames[number] = {
       "files/pythia_pp_rhic_45GeV_fullevent_44pthat45_on_April19.root",
        "files/pythia_pp_rhic_45GeV_fullevent_44pthat45_mpioff_April19.root",
        "files/pythia_pp_rhic_45GeV_fullevent_44pthat45_isroff_April19.root",
        "files/pythia_pp_rhic_45GeV_fullevent_44pthat45_bothoff_April19.root",
    };
    
        const char* histoname = "EEC_w"; 
    
        TFile* files[number];
        TH1F* histos[number];
        TH1F* EEC_w_p[number];
        TH1F* hrebin[number];
        TH1F* hrebin2[number];
        
       for (int i = 0; i < number; ++i) {
        files[i] = TFile::Open(filenames[i]);
        if (!files[i] || files[i]->IsZombie()) {
            std::cerr << "Error opening file: " << filenames[i] << std::endl;
            return;
        }

        TH1F* htemp = dynamic_cast<TH1F*>(files[i]->Get(histoname));
        if (!htemp) {
            std::cerr << "Histogram not found in: " << filenames[i] << std::endl;
            return;
        }

        // Clone and rename to keep each histogram separate
        histos[i] = (TH1F*)htemp->Clone(Form("EEC_clone_%d", i));
        histos[i]->SetDirectory(0);  // Detach from file
        files[i]->Close();
        
        
        // Normalize to bin width
        for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
        float width = histos[i]->GetBinWidth(bin);
        float content = histos[i]->GetBinContent(bin);
        float error = histos[i]->GetBinError(bin);
            if (width > 0) {
            histos[i]->SetBinContent(bin, content / width);
            histos[i]->SetBinError(bin, error / width);
            }
        }
        // 
        histos[i]->Scale(1.0/ histos[i]->Integral("width"));
        
        EEC_w_p[i] = new TH1F(Form("EEC_w_p_%d", i), "", 60, 0, 60);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
            }
    
        hrebin[i] = (TH1F*)EEC_w_p[i]->Clone(Form("hrebin%d", i));
        hrebin[i]->Rebin(6);
        hrebin[i]->Scale(1./6);
        hrebin2[i] = (TH1F*)hrebin[i]->Clone(Form("hrebin2_%d", i));
        hrebin[i]->GetXaxis()->SetRange(1,3);
        EEC_w_p[i]->GetXaxis()->SetRange(19,42);
        hrebin2[i]->GetXaxis()->SetRange(8,10);
    }

    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 600);
    canvas1->cd();
    canvas1->SetLogy();

    // empty hist since we are rebinning on the edges 
    TH1F* EEC_empty = new TH1F("EEC_w_m", "WEEC", 60, 0, 60);
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(1e-3);
    EEC_empty->SetMaximum(1e3);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1/#frac{jet p_{T, leading} + jet p_{T, subleading}}{2}}^2 EEC ");
    EEC_empty->SetStats(0);
    EEC_empty->Draw("P E1");
    
    // horribly done very innefficient but it works :)
    for (int i = 0; i < number; ++i) {
        EEC_w_p[i]->SetStats(0);
        EEC_w_p[i]->SetMarkerSize(0.60);
        if (i == 1) EEC_w_p[i]->SetMarkerStyle(21 + i);   
        else EEC_w_p[i]->SetMarkerStyle(22 + i);
        hrebin[i]->SetStats(0);
        hrebin[i]->SetMarkerSize(0.60);
        if (i == 1 || i == 2) hrebin[i]->SetMarkerStyle(20 + i);   
        else hrebin[i]->SetMarkerStyle(21 + i);
        hrebin2[i]->SetStats(0);
        hrebin2[i]->SetMarkerSize(0.60);
        if (i == 1 || i == 2) hrebin2[i]->SetMarkerStyle(20 + i);   
        else hrebin2[i]->SetMarkerStyle(21 + i);
        if (i == 2) {
        EEC_w_p[i]->SetLineColor(418);
        EEC_w_p[i]->SetMarkerColor(418);
        hrebin[i]->SetLineColor(418);
        hrebin[i]->SetMarkerColor(418);
        hrebin2[i]->SetLineColor(418);
        hrebin2[i]->SetMarkerColor(418);
        }
        else {
        EEC_w_p[i]->SetLineColor(i + 1);
        EEC_w_p[i]->SetMarkerColor(i + 1);
        hrebin[i]->SetLineColor(i + 1);
        hrebin[i]->SetMarkerColor(i + 1);
        hrebin2[i]->SetLineColor(i + 1);
        hrebin2[i]->SetMarkerColor(i + 1);
        }
        hrebin[i]->Draw("P E1 SAME");
        EEC_w_p[i]->Draw("P E1 SAME");
        hrebin2[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
    }


    // Axis embellishments
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
    axis1->Draw();
    axis2->Draw(); 
    
    //EEC_w_p[0]->SetMarkerSize(2);
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 16);
    //drawText("#pi/2", 0.4875, .06, 17 );
    //drawText("#Delta#phi", 0.8, .025, 16);
    drawText("p+p #sqrt{s} = 91.2 GeV", 0.60, 0.88, 16);
    drawText("43 Gev < jet p_{T1} <= 45 GeV", 0.60, 0.85, 16);
    drawText("43 Gev < jet p_{T2} <= 35 GeV", 0.60, 0.82, 16);
    drawText("44.60 Gev < #hat{p}_{T} <= 45.25 GeV", 0.60, 0.79, 16);
    drawText("charged constituent p_{T} > .2 GeV", 0.60, .76, 16 );
    drawText("#Delta#phi > 7#pi/8 ", 0.60, .73, 16 );
    //drawText("Deltaphi", 0.60, .88, 16 );
    canvas1->Update();
    // Legend
    TLegend* legend = new TLegend(0.6, 0.125, 0.9, 0.225);
    legend->AddEntry(EEC_w_p[0], "MPI on ISR on", "pl");
    legend->AddEntry(EEC_w_p[2], "MPI on ISR off", "pl");
    legend->AddEntry(EEC_w_p[1], "MPI off ISR on", "pl");
    //legend->AddEntry(EEC_w_p[2], "MPI on ISR off", "pl");
    legend->AddEntry(EEC_w_p[3], "MPI off ISR off", "pl");
    legend->Draw();
    canvas1->Update();
    
    canvas1->SaveAs("plots/MPIISRToggledMay16.png");
}
