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

    const char* filenames[4] = {
        "pythia_pp_rhic_45GeV_fullevent_10pthat60_Apr9_isroff.root",
        "pythia_pp_rhic_45GeV_fullevent_10pthat60_Apr9_mpioff.root",
        "pythia_pp_rhic_45GeV_fullevent_10pthat60_Apr9_bothoff.root",
        "pythia_pp_rhic_15GeV_fullevent_10pthat60_april15_deltheta.root"
    };

    const char* histoname = "EEC_w";

    TFile* files[4];
    TH1F* histos[4];
    TH1F* EEC_w_p[4];
    
    for (int i = 0; i < 4; ++i) {
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
        double width = histos[i]->GetBinWidth(bin);
        double content = histos[i]->GetBinContent(bin);
        double error = histos[i]->GetBinError(bin);
            if (width > 0) {
            histos[i]->SetBinContent(bin, content / width);
            histos[i]->SetBinError(bin, error / width);
            }
        }
        // 
        
        EEC_w_p[i] = new TH1F(Form("EEC_w_p_%d", i), "", 60, 0, 60);
    
        EEC_w_p[i]->Fill(3.5, histos[i]->GetBinContent(1));
        EEC_w_p[i]->SetBinError(3.5, histos[i]->GetBinError(1));
        EEC_w_p[i]->Fill(10.5, histos[i]->GetBinContent(2)); 
        EEC_w_p[i]->Fill(16.5, histos[i]->GetBinContent(3));
        EEC_w_p[i]->Fill(44, histos[i]->GetBinContent(26));
        EEC_w_p[i]->Fill(50, histos[i]->GetBinContent(27));
        EEC_w_p[i]->Fill(56.5, histos[i]->GetBinContent(28));
                
        for (int bin = 20; bin <= 41; ++bin) {
            int origbin = bin - 16 ;
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(origbin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(origbin));
            }
    
        
        // Normalize to unit integral
        double integral = EEC_w_p[i]->Integral("width");  // use width-aware integral
        if (integral > 0) EEC_w_p[i]->Scale(1.0 / integral);
    }

    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 600);
    canvas1->SetLogy();

    EEC_w_p[0]->SetLineColor(kRed);
    EEC_w_p[1]->SetLineColor(kBlue);
    EEC_w_p[2]->SetLineColor(kGreen + 2);
    EEC_w_p[3]->SetLineColor(kMagenta);

    EEC_w_p[0]->SetMarkerStyle(20);
    EEC_w_p[1]->SetMarkerStyle(21);
    EEC_w_p[2]->SetMarkerStyle(22);
    EEC_w_p[3]->SetMarkerStyle(23);
    
    EEC_w_p[0]->SetMarkerSize(0.5);
    EEC_w_p[0]->SetMinimum(1e-5);
    EEC_w_p[0]->SetMaximum(1);
    EEC_w_p[0]->Draw("P");
    EEC_w_p[0]->GetXaxis()->SetLabelOffset(999);
    EEC_w_p[0]->GetXaxis()->SetTitleOffset(999);
    EEC_w_p[0]->GetXaxis()->SetTickLength(0);
    EEC_w_p[0]->SetStats(0);
    
    for (int i = 1; i < 4; ++i) {
        EEC_w_p[i]->SetMarkerSize(0.5);
        EEC_w_p[i]->SetMinimum(1e-5);
        EEC_w_p[i]->SetMaximum(1);
        EEC_w_p[i]->Draw("P SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        EEC_w_p[i]->SetStats(0);
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
    
    // Text
    drawText("0.5", 0.4875, .06, 17 );
    drawText("z = (1 - cos(#theta))/2", 0.8, .025, 18);
    drawText("p+p #sqrt{s} = 200 GeV", 0.40, 0.85, 12); 
    drawText("anti-k_{T} R_{jet} = 0.4 |#eta < 0.7|", 0.40, 0.81, 12);
    drawText("44.75 GeV < #hat{p}_{T} < 45.25 GeV", 0.40, 0.77, 12);
    drawText("43 GeV < jet p_{T, leading} < 45 GeV", 0.40, 0.73, 12);
    drawText("43 GeV < jet p_{T, subleading} < 45 GeV", 0.40, 0.69, 12);
    drawText("charged constituent p_{T} > 0.2 GeV", 0.40, .65, 12);
    
    // Legend
    TLegend* legend = new TLegend(0.7, 0.15, 0.9, 0.35);
    legend->AddEntry(EEC_w_p[0], "ISR off ", "l");
    legend->AddEntry(EEC_w_p[1], "MPI off", "l");
    legend->AddEntry(EEC_w_p[2], "ISR & MPI off", "l");
    legend->AddEntry(EEC_w_p[3], "ISR & MPI on", "l");
    legend->Draw();

    canvas1->SaveAs("MPIISRdifference.png");
}
