#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <iostream>

void drawText(const char *text, double xp, double yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(8);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

TH1D* makeRatio(TH1D* h1, TH1D* h2, Color_t color, const char* name) {
    TH1D* ratio = (TH1D*)h1->Clone(name);
    ratio->Divide(h2);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(0.75);
    ratio->SetLineColor(color);
    ratio->SetMarkerColor(color);
    ratio->SetStats(0);
    return ratio;
}

void plotComparison(TH1D* h1, TH1D* h2, TH1D* h3,
                    const char* label1, const char* label2, const char* label3,
                    const char* saveNameC, const char* saveNameP) {

    TCanvas* c = new TCanvas(Form("c_%s", saveNameC), "Comparison", 800, 800);
    c->Divide(1, 2);
    
    TPad* pad1 = (TPad*)c->cd(1);
    pad1->SetPad(0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.0);
    pad1->SetLogy();

    TPad* pad2 = (TPad*)c->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);

    // Top pad histograms
    pad1->cd();
    h1->SetMarkerStyle(20); h1->SetMarkerColor(kBlack); h1->Draw("P E1");
    h2->SetMarkerStyle(21); h2->SetMarkerColor(kRed);   h2->Draw("P E1 SAME");
    h3->SetMarkerStyle(22); h3->SetMarkerColor(kBlue);  h3->Draw("P E1 SAME");

    TLegend* leg = new TLegend(0.65, 0.05, 0.89, 0.25);
    leg->AddEntry(h1, label1, "pl");
    leg->AddEntry(h2, label2, "pl");
    leg->AddEntry(h3, label3, "pl");
    leg->Draw();

    // Bottom pad ratios
    pad2->cd();
    TH1D* r12 = makeRatio(h1, h2, kBlack, "ratio12");
    TH1D* r13 = makeRatio(h1, h3, kRed,   "ratio13");
    TH1D* r23 = makeRatio(h2, h3, kBlue,  "ratio23");

    r12->GetYaxis()->SetTitle("Ratio");
    r12->SetMinimum(0.9);
    r12->SetMaximum(1.1);
    r12->Draw("P E1");
    r13->Draw("P E1 SAME");
    r23->Draw("P E1 SAME");

    TLine* line = new TLine(0, 1, r12->GetNbinsX(), 1);
    line->SetLineStyle(2); line->SetLineWidth(2); line->Draw();

    // Save
    c->SaveAs(saveNameC);
    c->SaveAs(saveNameP);
}
using namespace std;
void checkermbs() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
   /* const char* filename = {
       "../outfiles/pythia_mixedevent_03_15pthat50_May27_50m.root"};
       //"../outfiles/pythia_pp_pt12ic_fullevent_5pthat60_CentralityStudy03_May22_50m.root"};
       pythia_pp_pt12ic_fullevent_5pthat60_CentralityStudy03_mbs_thermalcopy_Aug28_1M.root
       */
    TFile* file = TFile::Open("../outfiles/triplechecker_mbs_xtrahists_100k_centrality03_Sept12_5M.root");
    //TFile* file = TFile::Open("../../outfiles/triplechecker_mbs_xtrahists_100k_morestatistics_2025-09-12_1.root");
    
    int const number = 9;
    double min = 1e7; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e8;
    
    double min1 = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max1 = 1e1;
          
    
    
    //mixed event 
    //vector<const char*> histonames = { "EEC_sp"};
    //vector<const char*> histonames = {"EEC_w", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_s"};
    vector<const char*> histonames = {"EEC_pt1", "EEC_pt2","EEC_pt3", "EEC_t1t1", "EEC_t2t2", "EEC_t3t3", "EEC_t1t2", "EEC_t1t3", "EEC_t2t3"};
    //vector<const char*> histonames = { "EEC_t1t1", "EEC_tt"};
    //vector<const char*> histonames = {"EEC_w", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_s", "EEC_t", "EEC_m"};
    //vector<const char*> histonames = {"EEC_wp", "EEC_msp", "EEC_mmp", "EEC_m2p", "EEC_sp", "EEC_tp", "EEC_mp"};
    
    bool rebin = 0;
    bool delphi = 0;
    bool integrate = 0;
    
    TH1D* histos[number];
    TH1D* EEC_w_p[number];
 
    TH1D* pt12;
    TH1D* pt13;
    TH1D* pt23;
    TH1D* p2;
    TH1D* p3;
    
    TH1D* t12t12;
    TH1D* t13t13;
    TH1D* t23t23;
    TH1D* t2;
    TH1D* t3;
    
    TH1D* t12t23;
    TH1D* t11t23;
    TH1D* t21t33;
    TH1D* t23;
    TH1D* t13;
    
    string outputFilename = "../plots/weecmathcheckSep12mixedthermal.C";
    string outputFilename1 = "../plots/weecmathcheckSep12mixedthermal.png";
    const char* pttrackcut = "charged constituent p_{T} > 0.2 GeV";

    
    int bins = 40;
    
    for (int i = 0; i < number; ++i) {
    TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
            // to keep it separte and not mess up your root file
            histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            histos[i]->SetDirectory(0);

            // Normalize to bin width
            for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
                double width = histos[i]->GetBinWidth(bin);
                double content = histos[i]->GetBinContent(bin);
                double error = histos[i]->GetBinError(bin);
                histos[i]->SetBinContent(bin, content / width);
            }
    }
   
    for (int i = 0; i < number; ++i) {        
        if (integrate){
        histos[i]->Scale(1.0/ histos[i]->Integral("width"));
        }

        EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", bins, 0, bins);
        
        for (int bin = 1; bin <= bins; ++bin) {
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
            double newcontent = histos[i]->GetBinContent(bin);
            //cout << newcontent << " " << bin << endl;
            if (newcontent > 0 && newcontent < min) min = newcontent;
            if (newcontent > max) max = newcontent;
        }
    } //closes i = loop


 // pre bin width subtraction + bin width & integral normalization + log to linear shift (final histogram )  :)
    
    //plot :)
    TCanvas* canvas1 = new TCanvas("c", "Comparison", 800, 800);
    canvas1->Divide(1, 2); // Split into 2 vertical pads
    canvas1->SetLeftMargin(0.125); 
    
    TPad* pad1 = (TPad*)canvas1->cd(1);
    pad1->SetPad(0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0.0);
    pad1->SetLogy();
    
    TPad* pad2 = (TPad*)canvas1->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.3); //
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    
    pad1->cd();

    // empty hist since we are rebinning on the edges 
    TH1D* EEC_empty = new TH1D("EEC_empty", "", bins, 0, bins);
    EEC_empty->SetStats(0);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1 / (#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2} EEC");
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    for (int i = 6; i < 9; ++i) {
        //marker size changes & style assignments 
        cout << i << endl;
        EEC_w_p[i]->SetStats(0);
        //20 is filled in circle 21 us
        if (i < 2) {
            EEC_w_p[i]->SetMarkerStyle(20 + i);   
    
            EEC_w_p[i]->SetMarkerSize(.75);
        }
        else if (i == 2){
            EEC_w_p[i]->SetMarkerStyle(23);   
            EEC_w_p[i]->SetMarkerSize(1);
    
        }
        else { 
            EEC_w_p[i]->SetMarkerStyle(33);
            EEC_w_p[i]->SetMarkerSize(1.25);
        }
        
        //color assignments 418 is green+2, 616 is kmagenta 
        if (i == 2) {
            EEC_w_p[i]->SetLineColor(418);
            EEC_w_p[i]->SetMarkerColor(418);
        }
        
        else if(i == 1) {
            EEC_w_p[i]->SetLineColor(616);
            EEC_w_p[i]->SetMarkerColor(616);
        }
        
        else if(i == 4) {
            EEC_w_p[i]->SetLineColor(kYellow+1);
            EEC_w_p[i]->SetMarkerColor(kYellow+1);
        }
         else if(i == 5) {
            EEC_w_p[i]->SetLineColor(46);
            EEC_w_p[i]->SetMarkerColor(46);
        }
        
        else if(i == 6) {
            EEC_w_p[i]->SetLineColor(kViolet-1);
            EEC_w_p[i]->SetMarkerColor(kViolet-1);
        }
        else { // 0 and 3
            EEC_w_p[i]->SetLineColor(i + 1);
            EEC_w_p[i]->SetMarkerColor(i + 1);
            }
        // just plot whole & signal 
        /*
        if(i == 0 || i == 4){
        EEC_w_p[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }

        */
        if (i >=0 ){
        EEC_w_p[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }
    }
  
    
    TGaxis *axis1 = new TGaxis(0, 0, 20, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(40, 0, 20, 0, .00001 , .5 , 510, "-G");
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
    
 
    TGaxis *axis3 = new TGaxis( 40, min , 40, max , min, max, 510,"+LG");
    axis3->SetLabelFont(42);
    axis3->SetLabelSize(.035);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(0, max, 20, max, .00001, .5, 510, "-G"); 
    TGaxis *axis5 = new TGaxis(40, max, 20, max, .00001 , .5 , 510, "+G");
    axis4->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis4->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    axis4->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis5->ChangeLabel(1, -1, 0, -1, -1, 62, "-" );
    axis5->ChangeLabel(3, -1, 0, -1, -1, 62, "-" );
    if(delphi){
        axis5->ChangeLabel(2, -1, -1, -1, -1, 62, "#pi - 10^{-4}" );
        axis5->ChangeLabel(4, -1, -1, -1, -1, 62, "#pi - 10^{-2}" );
        }
    else{
        axis5->ChangeLabel(2, -1, -1, -1, -1, 62, "1 - 10^{-4}" );
        axis5->ChangeLabel(4, -1, -1, -1, -1, 62, "1 - 10^{-2}" );
        }
    axis5->ChangeLabel(5, -1, 0, -1, -1, 62, "-" );
    axis4->SetLabelOffset(0.015);
    axis5->SetLabelOffset(-0.025);
    axis5->SetLabelFont(42);
    axis4->SetLabelSize(.03);
    axis5->SetLabelSize(.03);
    axis4->SetTickSize(0.05);
    axis5->SetTickSize(0.05);
    
    axis1->Draw();
    axis2->Draw();
    axis3->Draw();
    axis4->Draw();
    axis5->Draw();
    
    
    if(delphi){
        drawText("#pi/2", 0.505, .065, 17 );
        drawText("#pi/2", 0.505, .91, 17 );
        drawText("#Delta#phi", 0.8, .025, 16);
    }
    else{
        drawText("0.5", 0.505, .065, 17 );
        drawText("0.5", 0.50, .91, 17 );
        drawText("z = #frac{(1 - cos(#theta))}{2}", 0.8, .025, 16);
    }
    if(delphi){
    //rawText("NO MPI", 0.15, 0.36, 15);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.32, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.28, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.24, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.20, 15);
    drawText(pttrackcut, 0.15, 0.16, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.12, 15);
    canvas1->Update();
    }
    else{
    //drawText("NO MPI", 0.15, 0.36, 15);
    drawText("p+p #sqrt{s} = 200 GeV", 0.15, 0.32, 15);
    drawText("|#phi_{jet_{1}} - #phi_{jet_{2}}| > 3#pi/4", 0.15, 0.28, 15);
    drawText("31.2 #leq jet p_{T1} < 40.7 GeV", 0.15, 0.24, 15);
    drawText("27.3 GeV #leq jet p_{T2} < 31.2 GeV", 0.15, 0.20, 15);
    drawText(pttrackcut, 0.15, 0.16, 15);
    drawText("1500 thermal particles || Centrality 0-3%", 0.15, 0.12, 15);
    canvas1->Update();
    }
    /*
    TLegend* legend = new TLegend(0.55, 0.125, 0.875, 0.375);
    legend->AddEntry(EEC_w_p[0], "wEEC", "pl");
    legend->AddEntry(EEC_w_p[4], "Signal + Signal", "pl");
    legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl");
    legend->Draw();


    */
    
    TLegend* legend = new TLegend(0.65, 0.05, 0.89, 0.25);
    
    /*
    legend->AddEntry(EEC_w_p[3], "Thermal_{1} + Thermal_{1}", "pl");
    legend->AddEntry(EEC_w_p[4], "Thermal_{2} + Thermal_{2}", "pl");
    legend->AddEntry(EEC_w_p[5], "Thermal_{3} + Thermal_{3}", "pl");
     
    */
    legend->AddEntry(EEC_w_p[6], "Thermal_{1} + Thermal_{2}", "pl");
    legend->AddEntry(EEC_w_p[7], "Thermal_{1} + Thermal_{3}", "pl");
    legend->AddEntry(EEC_w_p[8], "Thermal_{2} + Thermal_{3}", "pl");
    
    /*
    legend->AddEntry(EEC_w_p[0], "Pythia + Thermal_{1}", "pl");
    legend->AddEntry(EEC_w_p[1], "Pythia + Thermal_{2}", "pl");
    legend->AddEntry(EEC_w_p[2], "Pythia + Thermal_{3}", "pl");
    */
   // legend->Draw();
    
    pad2->cd();
    pt12 = (TH1D*)EEC_w_p[6]->Clone();
    p2 = (TH1D*)EEC_w_p[7]->Clone();
  
    pt12->Divide(p2);
    
    pt12->SetTitle("");
    pt12->GetYaxis()->SetTitle("Ratio");
    pt12->GetYaxis()->SetTitleSize(0.1);
    pt12->GetYaxis()->SetTitleOffset(0.4);
    pt12->GetYaxis()->SetLabelSize(0.08);
    pt12->GetXaxis()->SetLabelSize(0.08);
    pt12->GetXaxis()->SetTitleSize(0.1);
    pt12->GetXaxis()->SetTitleOffset(1.0);
    pt12->SetMarkerStyle(20);
    pt12->SetMarkerSize(0.75);
    pt12->SetStats(0);
    pt12->SetLineColor(kBlack);
    pt12->SetMarkerColor(kBlack);
    pt12->SetMarkerColor(kBlack);
    pt12->GetXaxis()->SetLabelOffset(999);
    pt12->GetXaxis()->SetTitleOffset(999);
    pt12->GetXaxis()->SetTickLength(0);
    
    pt13 = (TH1D*)EEC_w_p[6]->Clone();
    p3 = (TH1D*)EEC_w_p[8]->Clone();
    pt13->Divide(p3);
    
    pt13->SetTitle("");
    pt13->GetYaxis()->SetTitle("Ratio1");
    pt13->GetYaxis()->SetTitleSize(0.1);
    pt13->GetYaxis()->SetTitleOffset(0.4);
    pt13->GetYaxis()->SetLabelSize(0.08);
    pt13->GetXaxis()->SetLabelSize(0.08);
    pt13->GetXaxis()->SetTitleSize(0.1);
    pt13->GetXaxis()->SetTitleOffset(1.0);
    pt13->SetMarkerStyle(20);
    pt13->SetMarkerSize(0.75);
    pt13->SetStats(0);
    pt13->SetLineColor(kCyan+3);
    pt13->SetMarkerColor(kCyan+3);
    pt13->SetMarkerColor(kCyan+3);
    pt13->GetXaxis()->SetLabelOffset(999);
    pt13->GetXaxis()->SetTitleOffset(999);
    pt13->GetXaxis()->SetTickLength(0);
    
    pt23 = (TH1D*)EEC_w_p[7]->Clone();
    p3 = (TH1D*)EEC_w_p[8]->Clone();
    pt23->Divide(p3);
    
    pt23->SetTitle("");
    pt23->GetYaxis()->SetTitle("Ratio1");
    pt23->GetYaxis()->SetTitleSize(0.1);
    pt23->GetYaxis()->SetTitleOffset(0.4);
    pt23->GetYaxis()->SetLabelSize(0.08);
    pt23->GetXaxis()->SetLabelSize(0.08);
    pt23->GetXaxis()->SetTitleSize(0.1);
    pt23->GetXaxis()->SetTitleOffset(1.0);
    pt23->SetMarkerStyle(20);
    pt23->SetMarkerSize(0.75);
    pt23->SetStats(0);
    pt23->SetLineColor(kViolet);
    pt23->SetMarkerColor(kViolet);
    pt23->SetMarkerColor(kViolet);
    pt23->GetXaxis()->SetLabelOffset(999);
    pt23->GetXaxis()->SetTitleOffset(999);
    pt23->GetXaxis()->SetTickLength(0);    
    
    // Axis limits
    pt12->SetMinimum(0.995);
    pt12->SetMaximum(1.005);
    
    
    pt12->Draw("P E1");
    pt13->Draw("P E1 SAME");
    pt23->Draw("P E1 SAME");
    TLegend* legend1 = new TLegend(0.75, 0.125, 0.975, 0.375);
    /*
    legend1->AddEntry(pt12, " P+T_{1}/P+T_{2}", "pl");
    legend1->AddEntry(pt13, " P+T_{1}/P+T_{3}", "pl");
    legend1->AddEntry(pt23, " P+T_{2}/P+T_{3}", "pl"); 
    
    
    */
    legend1->AddEntry(pt12, " T_{1}+T_{2}/T_{1}+T_{3}", "pl");
    legend1->AddEntry(pt13, " T_{1}+T_{2}/T_{2}+T_{3}", "pl");
    legend1->AddEntry(pt23, " T_{2}+T_{3}/T_{1}+T_{3}", "pl"); 
    
    /*
    legend1->AddEntry(pt12, " T_{1}+T_{1}/T_{2}+T_{2}", "pl");
    legend1->AddEntry(pt13, " T_{1}+T_{1}/T_{3}+T_{3}", "pl");
    legend1->AddEntry(pt23, " T_{2}+T_{2}/T_{3}+T_{3}", "pl"); 
*/
    
    legend1->Draw();
    
    TLine *line = new TLine(0, 1, bins, 1);
    
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw();
    
    axis1->Draw();
    axis2->Draw();
    axis3->Draw();
    axis4->Draw();
    axis5->Draw();
    
    
    canvas1->SaveAs(outputFilename.c_str());
    canvas1->SaveAs(outputFilename1.c_str());
    

}

    
    
    
        
