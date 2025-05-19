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

using namespace std;

TH1::SetDefaultSumw2();
TH2::SetDefaultSumw2();
    
void threeto1() {
    TFile* fin = TFile::Open("../outfiles/thermalbackgroundphistudy_10k_2025-05-16_1.root");
    
    const char* histname1 = "EEC_bl1";
    const char* histname2 = "EEC_bl2";
    const char* histname3 = "EEC_bl3";
    
    bool linear = false;
    bool logplot = true;
    bool reverse = false;
    
    int bins = 60;
    
    TGaxis *axis1 = new TGaxis(0, 0, 30, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(60, 0, 30, 0, .00001 , .5, 510, "-G");
    axis1->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis1->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-4}" );
    axis2->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis2->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-2}" );
    axis2->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis2->SetLabelOffset(0.045);
    axis1->SetLabelSize(.03);
    axis2->SetLabelSize(.03);
    axis1->SetTickSize(0.05);
    axis2->SetTickSize(0.05);
    
    TH1F* histEEC1 = (TH1F*)fin->Get(histname1);
    TH1F* histEEC2 = (TH1F*)fin->Get(histname2);
    TH1F* histEEC3 = (TH1F*)fin->Get(histname3);
    
    TCanvas* canvas1 = new TCanvas("canvas1", "Whole event EEC1", 800, 600);
    TCanvas* canvas2 = new TCanvas("canvas2", "Whole event EEC2", 800, 600);
    TCanvas* canvas3 = new TCanvas("canvas3", "Whole event EEC3", 800, 600);
    
    if (linear){
        canvas1->cd();
        TH1F* histo1 = (TH1F*)histEEC1->Clone("histo1");
        TH1F* histo2 = (TH1F*)histEEC2->Clone("histo2");
        TH1F* histo3 = (TH1F*)histEEC3->Clone("histo3");
        
        TH1F* EEC_empty = new TH1F("EEC_empty", "Linear Delta Phi 240 bins", bins, 0, TMath::Pi());
        EEC_empty->SetStats(0);
        //EEC_empty->GetXaxis()->SetLabelOffset(999);
        //EEC_empty->GetXaxis()->SetTitleOffset(999);
        //EEC_empty->GetXaxis()->SetTickLength(0);
        EEC_empty->Draw("P E1");
        
        histo1->SetStats(0);
        histo1->GetXaxis()->SetLabelOffset(999);
        histo1->GetXaxis()->SetTitleOffset(999);
        histo1->GetXaxis()->SetTickLength(0);
        histo1->Scale(1.0 / histo1->Integral("width"));
        histo1->SetLineColor(kMagenta);
        histo1->SetMarkerColor(kMagenta);
        histo1->SetMarkerStyle(20);
        histo1->SetMarkerSize(0.75);
        
        histo2->SetStats(0);
        histo2->GetXaxis()->SetLabelOffset(999);
        histo2->GetXaxis()->SetTitleOffset(999);
        histo2->GetXaxis()->SetTickLength(0);
        histo2->Scale(1.0 / histo2->Integral("width"));
        histo2->SetLineColor(1);
        histo2->SetMarkerColor(1);
        histo2->SetMarkerStyle(2);
        histo2->SetMarkerSize(0.75);
        
        histo3->SetStats(0);
        histo3->GetXaxis()->SetLabelOffset(999);
        histo3->GetXaxis()->SetTitleOffset(999);
        histo3->GetXaxis()->SetTickLength(0);
        histo3->Scale(1.0 / histo3->Integral("width"));
        histo3->SetLineColor(kBlue);
        histo3->SetMarkerColor(kBlue);
        histo3->SetMarkerStyle(36);
        histo3->SetMarkerSize(0.75);  
        
        histo1->Draw("P E1 SAME");
        histo2->Draw("P E1 SAME");
        histo3->Draw("P E1 SAME");
        
        TLegend* legend = new TLegend(0.7, 0.725, 0.9, 0.875);
        legend->AddEntry(histo1, "100 events", "pl");
        legend->AddEntry(histo2, "1000 events", "pl");
        legend->AddEntry(histo3, "10000 events", "pl");
        legend->Draw();
        
        
    canvas1->SaveAs("../plots/BigLinearDeltaPhi3to1.png");
    }// end of linear == true
    
    if (logplot){
        TH1F* histo1 = (TH1F*)histEEC1->Clone("histo1");
        TH1F* histo2 = (TH1F*)histEEC2->Clone("histo2");
        TH1F* histo3 = (TH1F*)histEEC3->Clone("histo3");
        
         //first histogram
        for (int bin = 1; bin <= histo1->GetNbinsX(); ++bin) {
        float width = histo1->GetBinWidth(bin);
        float content = histo1->GetBinContent(bin);
        float error = histo1->GetBinError(bin);
        histo1->SetBinContent(bin, content / width);
        histo1->SetBinError(bin, error / width);
        }
        
        histo1->Scale(1.0/histo1->Integral("width"));
        TH1F*EEC_1 = new TH1F("EEC_1","Double log 100", bins , 0, bins);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_1->SetBinContent(bin, histo1->GetBinContent(bin));
            EEC_1->SetBinError(bin, histo1->GetBinError(bin));
            }
        

         //second histogram
        for (int bin = 2; bin <= histo2->GetNbinsX(); ++bin) {
        float width = histo2->GetBinWidth(bin);
        float content = histo2->GetBinContent(bin);
        float error = histo2->GetBinError(bin);
        histo2->SetBinContent(bin, content / width);
        histo2->SetBinError(bin, error / width);
        }
        
        histo2->Scale(1.0/histo2->Integral("width"));
        TH1F* EEC_2 = new TH1F("EEC_2","Double log 1000", bins, 0, bins);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_2->SetBinContent(bin, histo2->GetBinContent(bin));
            EEC_2->SetBinError(bin, histo2->GetBinError(bin));
            }
        
        

         //third histogram
        for (int bin = 1; bin <= histo3->GetNbinsX(); ++bin) {
        float width = histo3->GetBinWidth(bin);
        float content = histo3->GetBinContent(bin);
        float error = histo3->GetBinError(bin);
        histo3->SetBinContent(bin, content / width);
        histo3->SetBinError(bin, error / width);
        }
        
        histo3->Scale(1.0/histo3->Integral("width"));
        TH1F* EEC_3 = new TH1F("EEC_3","Double Log 10k", bins, 0,bins);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_3->SetBinContent(bin, histo3->GetBinContent(bin));
            EEC_3->SetBinError(bin, histo3->GetBinError(bin));
            }
        
       
        
        canvas1->cd();
        canvas1->SetLogy();
        
        // plotting 
        TH1F* EEC_empty = new TH1F("EEC_empty", "WEEC", bins, 0, TMath::Pi());
        EEC_empty->SetStats(0);
        EEC_empty->GetXaxis()->SetLabelOffset(999);
        EEC_empty->GetXaxis()->SetTitleOffset(999);
        EEC_empty->GetXaxis()->SetTickLength(0);
        //EEC_empty->Draw("P E1");
        
        EEC_1->SetStats(0);
        EEC_1->GetXaxis()->SetLabelOffset(999);
        EEC_1->GetXaxis()->SetTitleOffset(999);
        EEC_1->GetXaxis()->SetTickLength(0);
        EEC_1->SetLineColor(kMagenta);
        EEC_1->SetMarkerColor(kMagenta);
        EEC_1->SetMarkerStyle(25);
        EEC_1->SetMarkerSize(0.75);
        EEC_1->SetMinimum(1e-4);
        EEC_1->SetMaximum(1);
        
        EEC_2->SetStats(0);
        EEC_2->GetXaxis()->SetLabelOffset(999);
        EEC_2->GetXaxis()->SetTitleOffset(999);
        EEC_2->GetXaxis()->SetTickLength(0);
        EEC_2->SetLineColor(1);
        EEC_2->SetMarkerColor(1);
        EEC_2->SetMarkerStyle(24);
        EEC_2->SetMarkerSize(0.75);
        
        EEC_3->SetStats(0);
        EEC_3->GetXaxis()->SetLabelOffset(999);
        EEC_3->GetXaxis()->SetTitleOffset(999);
        EEC_3->GetXaxis()->SetTickLength(0);
        EEC_3->SetLineColor(kBlue);
        EEC_3->SetMarkerColor(kBlue);
        EEC_3->SetMarkerStyle(21);
        EEC_3->SetMarkerSize(0.75);
        
        EEC_1->Draw("P E1");
        EEC_2->Draw("P E1 SAME");
        EEC_3->Draw("P E1 SAME");
        
        axis1->Draw();
        axis2->Draw();
        
        drawText("#pi/2", 0.4875, .06, 17 );
        
        TLegend* legend = new TLegend(0.7, 0.725, 0.9, 0.875);
        legend->AddEntry(EEC_1, "100 events", "pl");
        legend->AddEntry(EEC_2, "1000 events", "pl");
        legend->AddEntry(EEC_3, "10000 events", "pl");
        legend->Draw();
        
        canvas1->SaveAs("../plots/BiglinearDeltaPhi3to1.png");
    }// end of logplot
    
    if (reverse){
    
        Double_t topi[61] =   {1.00000000e-05, 1.49006082e-05, 2.22028125e-05, 3.30835411e-05,
                                4.92964884e-05, 7.34547660e-05, 1.09452069e-04, 1.63090240e-04,
                                2.43014377e-04, 3.62106202e-04, 5.39560265e-04, 8.03977611e-04,
                                1.19797554e-03, 1.78505642e-03, 2.65984263e-03, 3.96332730e-03,
                                5.90559873e-03, 8.79970130e-03, 1.31120901e-02, 1.95378118e-02,
                                2.91125279e-02, 4.33794373e-02, 6.46379999e-02, 9.63145513e-02,
                                1.43514539e-01, 2.13845393e-01, 3.18642641e-01, 4.74796916e-01,
                                7.07476283e-01, 1.05418269e+00, 1.57079633e+00, 2.08761365,
                                2.43432005, 2.66799942, 2.82315369, 2.92795094, 
                                2.99828179, 3.04548178, 3.07715833, 3.09841689, 
                                3.11268380, 3.12225852, 3.12868424, 3.13299663, 
                                3.13589073, 3.13783300, 3.13813649, 3.13981128, 
                                3.14059876, 3.14099235, 3.14115677, 3.14123423, 
                                3.14135331, 3.14143323, 3.14148688, 3.14155555, 
                                3.14156300, 3.14157005, 3.14157413, 3.14157740, 
                                3.14158265};
        // histogram 1
        TH1F* tolog_1 = new TH1F("tolog_1", "240 bin linear delta phi to 60 bin double log", 60, topi);
        
        //histEEC1->Scale(1.0/histEEC1->Integral("width"));
        // basically trying to figure out where the contents of each old bin goes 
         for (int old_bin = 1; old_bin <= bins; ++old_bin) {
            double width  = histEEC1->GetBinWidth(old_bin);
            double content = width * histEEC1->GetBinContent(old_bin);
            double error = width * histEEC1->GetBinError(old_bin);
            double center = histEEC1->GetBinCenter(old_bin);
        
            int target_bin = tolog_1->FindBin(center);
        
            // Get current content and error squared in the new bin
            double existing_content = tolog_1->GetBinContent(target_bin);
            double existing_error2 = pow(tolog_1->GetBinError(target_bin), 2);
        
            // Add new content and error in quadrature
            double updated_content = existing_content + content;
            double updated_error2 = existing_error2 + pow(error,2);
        
            tolog_1->SetBinContent(target_bin, updated_content);
            tolog_1->SetBinError(target_bin, sqrt(updated_error2));
        }
       
       
        // Normalize by bin width (after tolog is filled)
        for (int bin = 1; bin <= tolog_1->GetNbinsX(); ++bin) {
            float width = tolog_1->GetBinWidth(bin);
            double content = tolog_1->GetBinContent(bin);
            double error = tolog_1->GetBinError(bin);
        
            tolog_1->SetBinContent(bin, content / width);
            tolog_1->SetBinError(bin, error / width);
        }
        
        tolog_1->Scale(1.0/tolog_1->Integral("width"));
        TH1F*EEC_ul1 = new TH1F("EEC_ul1","Double log 100", 60 , 0, TMath::Pi());
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_ul1->SetBinContent(bin, tolog_1->GetBinContent(bin));
            EEC_ul1->SetBinError(bin, tolog_1->GetBinError(bin));
            }
        
        //EEC_ul1->Scale(1.0/EEC_ul1->Integral("width"));
         
         
         // histogram 2
        TH1F* tolog_2 = new TH1F("tolog_2", "240 bin linear delta phi to 60 bin double log", 60, topi);
        
        for (int old_bin = 1; old_bin <= bins; ++old_bin) {
            double width  = histEEC2->GetBinWidth(old_bin);
            double content = width * histEEC2->GetBinContent(old_bin);
            double error = width * histEEC2->GetBinError(old_bin);
            double center = histEEC2->GetBinCenter(old_bin);
        
            int target_bin = tolog_2->FindBin(center);
        
            // Get current content and error squared in the new bin
            double existing_content = tolog_2->GetBinContent(target_bin);
            double existing_error2 = pow(tolog_2->GetBinError(target_bin), 2);
        
            // Add new content and error in quadrature
            double updated_content = existing_content + content;
            double updated_error2 = existing_error2 + pow(error,2);
        
            tolog_2->SetBinContent(target_bin, updated_content);
            tolog_2->SetBinError(target_bin, sqrt(updated_error2));
        }
       
        // Normalize by bin width (after tolog is filled)
        for (int bin = 1; bin <= tolog_2->GetNbinsX(); ++bin) {
            float width = tolog_2->GetBinWidth(bin);
            double content = tolog_2->GetBinContent(bin);
            double error = tolog_2->GetBinError(bin);
        
                if (content != content || width == 0){
                cout << bin <<" "<< width<< " "<< content << " "<< error <<endl;

            }
            tolog_2->SetBinContent(bin, content / width);
            tolog_2->SetBinError(bin, error / width);
        }
        
        tolog_2->Scale(1.0/tolog_2->Integral("width"));
        TH1F*EEC_ul2 = new TH1F("EEC_ul2","Double log 100", 60 , 0, TMath::Pi());
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_ul2->SetBinContent(bin, tolog_2->GetBinContent(bin));
            EEC_ul2->SetBinError(bin, tolog_2->GetBinError(bin));
            }
   
          // histogram 3
        TH1F* tolog_3 = new TH1F("tolog_3", "240 bin linear delta phi to 60 bin double log", 60, topi);
        
        // basically trying to figure out where the contents of each old bin goes 
        for (int old_bin = 1; old_bin <= bins; ++old_bin) {
            double width  = histEEC3->GetBinWidth(old_bin);
            double content = width * histEEC3->GetBinContent(old_bin);
            double error = width * histEEC3->GetBinError(old_bin);
            double center = histEEC3->GetBinCenter(old_bin);
        
            int target_bin = tolog_3->FindBin(center);
        
            // Get current content and error squared in the new bin
            double existing_content = tolog_3->GetBinContent(target_bin);
            double existing_error2 = pow(tolog_3->GetBinError(target_bin), 2);
        
            // Add new content and error in quadrature
            double updated_content = existing_content + content;
            double updated_error2 = existing_error2 + pow(error,2);
        
            tolog_3->SetBinContent(target_bin, updated_content);
            tolog_3->SetBinError(target_bin, sqrt(updated_error2));
        }
       
        // Normalize by bin width (after tolog is filled)
        for (int bin = 1; bin <= tolog_3->GetNbinsX(); ++bin) {
            float width = tolog_3->GetBinWidth(bin);
            double content = tolog_3->GetBinContent(bin);
            double error = tolog_3->GetBinError(bin);
            if (content != content || width == 0){
                cout << bin <<" "<< width<< " "<< content << " "<< error <<endl;

            }
            
            tolog_3->SetBinContent(bin, content / width);
            tolog_3->SetBinError(bin, error / width);
        }
        
        tolog_3->Scale(1.0/tolog_3->Integral("width"));
        TH1F*EEC_ul3 = new TH1F("EEC_ul3","Double log 100", 60 , 0, TMath::Pi());
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_ul3->SetBinContent(bin, tolog_3->GetBinContent(bin));
            EEC_ul3->SetBinError(bin, tolog_3->GetBinError(bin));
            }
    
        
        canvas1->cd();
        canvas1->SetLogy();
       
        // plotting 
        TH1F* EEC_emptyd = new TH1F("EEC_emptyd", "WEEC", 60, 0, TMath::Pi());
        EEC_emptyd->SetStats(0);
        EEC_emptyd->GetXaxis()->SetLabelOffset(999);
        EEC_emptyd->GetXaxis()->SetTitleOffset(999);
        EEC_emptyd->GetXaxis()->SetTickLength(0);
        //EEC_emptyd->Draw("P E1");
        
        EEC_ul1->SetStats(0);
        EEC_ul1->GetXaxis()->SetLabelOffset(999);
        EEC_ul1->GetXaxis()->SetTitleOffset(999);
        EEC_ul1->GetXaxis()->SetTickLength(0);
        EEC_ul1->SetLineColor(kMagenta);
        EEC_ul1->SetMarkerColor(kMagenta);
        EEC_ul1->SetMarkerStyle(25);
        EEC_ul1->SetMarkerSize(0.75);
        //EEC_ul1->SetMinimum(1e-4);
        //tolog_1->SetMaximum(1);
        
        EEC_ul2->SetStats(0);
        EEC_ul2->GetXaxis()->SetLabelOffset(999);
        EEC_ul2->GetXaxis()->SetTitleOffset(999);
        EEC_ul2->GetXaxis()->SetTickLength(0);
        EEC_ul2->SetLineColor(1);
        EEC_ul2->SetMarkerColor(1);
        EEC_ul2->SetMarkerStyle(24);
        EEC_ul2->SetMarkerSize(0.75);
        
        EEC_ul3->SetStats(0);
        EEC_ul3->GetXaxis()->SetLabelOffset(999);
        EEC_ul3->GetXaxis()->SetTitleOffset(999);
        EEC_ul3->GetXaxis()->SetTickLength(0);
        EEC_ul3->SetLineColor(kBlue);
        EEC_ul3->SetMarkerColor(kBlue);
        EEC_ul3->SetMarkerStyle(29);
        EEC_ul3->SetMarkerSize(2.5);
       
        EEC_ul1->Draw("P E1");
        EEC_ul2->Draw("P E1 SAME");
        EEC_ul3->Draw("P E1 SAME");
        
        axis1->Draw();
        axis2->Draw();

        drawText("#pi/2", 0.4875, .06, 17 );
 
        TLegend* legend = new TLegend(0.7, 0.725, 0.9, 0.875);
        legend->AddEntry(EEC_ul1, "100 events", "pl");
        legend->AddEntry(EEC_ul2, "1000 events", "pl");
        legend->AddEntry(EEC_ul3, "10000 events", "pl");
        legend->Draw();

        canvas1->Draw();
        canvas1->SaveAs("../plots/edoublelogDeltaPhi3to1.png");
    
    }// end of linear == true 
    
    
    
}
    
    
    