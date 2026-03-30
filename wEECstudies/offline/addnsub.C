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


using namespace std;
void addnsub() {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
   /* const char* filename = {
       "../outfiles/pythia_mixedevent_03_15pthat50_May27_50m.root"};
       //"../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_May22_50m.root"};
       */
    TFile* file = TFile::Open("../outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy03_mixedeventb_50m_Aug4.root");
    cout <<"1" << endl;   
    
    int const number = 7;
    double min = 1e-10; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max = 1e10;
    
    double min1 = 1e-1; // initializing min and max for y axis because its easier to do right axis later. Def a better way to do this
    double max1 = 1e1;
          
    
    
    //mixed event 
    //vector<const char*> histonames = { "EEC_sp"};
    //vector<const char*> histonames = {"EEC_w", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_s"};
    vector<const char*> histonames = {"EEC_w", "EEC_s", "EEC_m", "EEC_t", "EEC_ms", "EEC_mm", "EEC_m2"};
    //vector<const char*> histonames = {"EEC_w", "EEC_ms", "EEC_mm", "EEC_m2", "EEC_s", "EEC_t", "EEC_m"};
    //vector<const char*> histonames = {"EEC_wp", "EEC_msp", "EEC_mmp", "EEC_m2p", "EEC_sp", "EEC_tp", "EEC_mp"};
    
    cout <<"1" << endl;   
    
    bool rebin = 0;
    bool delphi = 1;
    bool integrate = 0;
    
    TH1D* histos[number];
    TH1D* prehistos[number];
    TH1D* preEEC_w_p[number];
    TH1D* EEC_w_p[number];
    TH1D* mbssignal;
    TH1D* mbssignalfinal;
    TH1D* mbssignal1;
    TH1D* mbssignalfinal1;
    
    string outputFilename = "../plots/TESTweeccrecreations.C";
    string outputFilename1 = "../plots/TESTweeccreations.png";
    const char* pttrackcut = "charged constituent p_{T} > 0.2 GeV";
    
    
    int bins = 60;
        cout <<"1" << endl;   
    
    for (int i = 0; i < number; ++i) {
    TH1D* htemp = dynamic_cast<TH1D*>(file->Get(histonames[i]));
            // to keep it separte and not mess up your root file
            histos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            histos[i]->SetDirectory(0);
            
            // histos that will be subtracted 
            prehistos[i] = (TH1D*)htemp->Clone(Form("EEC_clone_%d", i));
            prehistos[i]->SetDirectory(0);
           
            // gonna do the subtraction before and after bin width and see what the heck goes on. 
            // Subtraction before normalization but idk which one fr
                cout <<"1" << endl;   
    
            // Normalize to bin width
            for (int bin = 1; bin <= histos[i]->GetNbinsX(); ++bin) {
                double width = histos[i]->GetBinWidth(bin);
                double content = histos[i]->GetBinContent(bin);
                double error = histos[i]->GetBinError(bin);
                histos[i]->SetBinContent(bin, content / width);
                histos[i]->SetBinError(bin, error / width);
            }
            //mbs stuff only
    }
    
    for (int i = 0; i < number; ++i) {        
        // regular histos[i] loop
        if (integrate){
        histos[i]->Scale(1.0/ histos[i]->Integral("width"));
        }
            cout <<"2" << endl;   
    
        EEC_w_p[i] = new TH1D(Form("EEC_w_p_%d", i), "", 60, 0, 60);
        
        for (int bin = 1; bin <= 60; ++bin) {
            EEC_w_p[i]->SetBinContent(bin, histos[i]->GetBinContent(bin));
            EEC_w_p[i]->SetBinError(bin, histos[i]->GetBinError(bin));
            double newcontent = histos[i]->GetBinContent(bin);
            //cout << newcontent << " " << bin << endl;
            if (newcontent > 0 && newcontent < min) min = newcontent;
            if (newcontent > max) max = newcontent;
        }
    } //closes i = loop

//pre bin width SUBTRACTION IS GONNA OCCUR NOW 
//Subhistos are subtracted then normalized by bin width then integral 
// okay so EEC_w - EEC_ms - EEC_mm + EEC_m2 = EEC_s allegedly and thats what this code does :)
    mbssignal = (TH1D*)prehistos[1]->Clone();
    mbssignal->Add(prehistos[2], +1);
    mbssignal->Add(prehistos[3], +1);
    //mbssignal->Add(prehistos[3], +1); //opposite sign of incorrect background modeling in sm1 so this IS the right sign :) (and the only way it works)
    
    //mbs has occured now to normaliza just mbssignal by bin width
    for (int bin = 1; bin <= mbssignal->GetNbinsX(); ++bin) {
            double width = mbssignal->GetBinWidth(bin);
            double content = mbssignal->GetBinContent(bin);
            double error = mbssignal->GetBinError(bin);
                mbssignal->SetBinContent(bin, content / width);
                mbssignal->SetBinError(bin, error / width);
    }
    
    if (integrate){
    mbssignal->Scale(1.0/mbssignal->Integral("width"));
    }
    
    // double log to linear for mbs 
    //reinitiallizing min and max 
    
    mbssignalfinal = new TH1D(Form("mbssignalfinal"), "", 60, 0, 60);
    for (int bin = 1; bin <= 60; ++bin) {
        mbssignalfinal->SetBinContent(bin, mbssignal->GetBinContent(bin));
        mbssignalfinal->SetBinError(bin, mbssignal->GetBinError(bin));
        double newcontent = mbssignal->GetBinContent(bin);
        //cout << newcontent << " " << bin << endl;
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
 // pre bin width subtraction + bin width & integral normalization + log to linear shift (final histogram )  :)
    
    mbssignal1 = (TH1D*)prehistos[1]->Clone();
    mbssignal1->Add(prehistos[4], +1);
    mbssignal1->Add(prehistos[5], -1);
    mbssignal1->Add(prehistos[6], +1);
    //mbssignal->Add(prehistos[3], +1); //opposite sign of incorrect background modeling in sm1 so this IS the right sign :) (and the only way it works)
    
    //mbs has occured now to normaliza just mbssignal by bin width
    for (int bin = 1; bin <= mbssignal1->GetNbinsX(); ++bin) {
            double width = mbssignal1->GetBinWidth(bin);
            double content = mbssignal1->GetBinContent(bin);
            double error = mbssignal1->GetBinError(bin);
                mbssignal1->SetBinContent(bin, content / width);
                mbssignal1->SetBinError(bin, error / width);
    }
    
    if (integrate){
    mbssignal1->Scale(1.0/mbssignal1->Integral("width"));
    }
    
    // double log to linear for mbs 
    //reinitiallizing min and max 
    
    mbssignalfinal1 = new TH1D(Form("mbssignalfinal1"), "", 60, 0, 60);
    for (int bin = 1; bin <= 60; ++bin) {
        mbssignalfinal1->SetBinContent(bin, mbssignal1->GetBinContent(bin));
        mbssignalfinal1->SetBinError(bin, mbssignal1->GetBinError(bin));
        double newcontent = mbssignal1->GetBinContent(bin);
        //cout << newcontent << " " << bin << endl;
        if (newcontent > 0 && newcontent < min) min = newcontent;
        if (newcontent > max) max = newcontent;
    }
    
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
    TH1D* EEC_empty = new TH1D("EEC_empty", "", 60, 0, 60);
    EEC_empty->SetStats(0);
    EEC_empty->GetXaxis()->SetLabelOffset(999);
    EEC_empty->GetXaxis()->SetTitleOffset(999);
    EEC_empty->GetXaxis()->SetTickLength(0);
    EEC_empty->GetYaxis()->SetTitle("1 / (#frac{jet_{1} p_{T} + jet_{2} p_{T}}{2})^{2} EEC");
    EEC_empty->SetStats(0);
    EEC_empty->SetMinimum(min);
    EEC_empty->SetMaximum(max);
    EEC_empty->Draw("P E1");
    
    for (int i = 0; i < number; ++i) {
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
            EEC_w_p[i]->SetLineColor(kYellow+1);
            EEC_w_p[i]->SetMarkerColor(kYellow+1);
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
        if (i > 0 && i < 6 ){
        EEC_w_p[i]->Draw("P E1 SAME");
        EEC_w_p[i]->GetXaxis()->SetLabelOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTitleOffset(999);
        EEC_w_p[i]->GetXaxis()->SetTickLength(0);
        }
        
    }
    
    
    
    mbssignalfinal->SetMarkerStyle(20);
	mbssignalfinal->SetMarkerColor(kCyan+3);
	mbssignalfinal->SetLineColor(kCyan+3);
	mbssignalfinal->Draw("P E1 SAME");
	
	mbssignalfinal1->SetMarkerStyle(21);
	mbssignalfinal1->SetMarkerColor(kOrange+3);
	mbssignalfinal1->SetLineColor(kOrange+3);
	mbssignalfinal1->Draw("P E1 SAME");
    
    TGaxis *axis1 = new TGaxis(0, 0, 30, 0, .00001, .5, 510, "G"); 
    TGaxis *axis2 = new TGaxis(60, 0, 30, 0, .00001 , .5 , 510, "-G");
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
    
 
    TGaxis *axis3 = new TGaxis( 60, min , 60, max , min, max, 510,"+LG");
    axis3->SetLabelFont(42);
    axis3->SetLabelSize(.035);
    axis3->Draw();
    
    TGaxis *axis4 = new TGaxis(0, max, 30, max, .00001, .5, 510, "-G"); 
    TGaxis *axis5 = new TGaxis(60, max, 30, max, .00001 , .5 , 510, "+G");
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

    /*
    */
    //
    TLegend* legend = new TLegend(0.65, 0.05, 0.89, 0.25);
    legend->AddEntry(EEC_w_p[0], "wEEC", "pl");
    legend->AddEntry(mbssignalfinal, "ss + sb + bb", "pl");
    legend->AddEntry(mbssignalfinal1, "ss + sm + mm", "pl");
    //legend->AddEntry(EEC_w_p[2], "Mixed1 + Mixed1 ", "pl");
    //legend->AddEntry(EEC_w_p[4], "Signal + Signal", "pl");
    //legend->AddEntry(EEC_w_p[5], "Thermal + Thermal", "pl");
    //legend->AddEntry(EEC_w_p[2], "Signal + Mixed", "pl");
    //legend->AddEntry(EEC_w_p[3], "Mixed + Mixed", "pl");
    //legend->AddEntry(EEC_w_p[4], "Signal + Signal", "pl");
    //legend->AddEntry(mbssignal1final, "Mixed background subtraction post bin", "pl");
    //legend->AddEntry(mbssignalfinal, "Mixed background subtraction", "pl"); 
    legend->Draw();
    
    pad2->cd();

    TH1D* ratioHist = (TH1D*)mbssignalfinal->Clone("ratioHist");
    ratioHist->Divide(mbssignalfinal1);
    
    // Styling
    ratioHist->SetTitle("");
    ratioHist->GetYaxis()->SetTitle("Ratio");
    ratioHist->GetYaxis()->SetTitleSize(0.1);
    ratioHist->GetYaxis()->SetTitleOffset(0.4);
    ratioHist->GetYaxis()->SetLabelSize(0.08);
    ratioHist->GetXaxis()->SetLabelSize(0.08);
    ratioHist->GetXaxis()->SetTitleSize(0.1);
    ratioHist->GetXaxis()->SetTitleOffset(1.0);
    ratioHist->SetMarkerStyle(20);
    ratioHist->SetMarkerSize(0.75);
    ratioHist->SetStats(0);
    ratioHist->SetLineColor(kBlack);
    ratioHist->SetMarkerColor(kBlack);
    ratioHist->GetXaxis()->SetLabelOffset(999);
    ratioHist->GetXaxis()->SetTitleOffset(999);
    ratioHist->GetXaxis()->SetTickLength(0);
    
    // Axis limits
    ratioHist->SetMinimum(-3.0);
    ratioHist->SetMaximum(3.0); 
    
    ratioHist->Draw("P E1");
    
    axis1->Draw();
    axis2->Draw();
    axis3->Draw();
    axis4->Draw();
    axis5->Draw();
    //axis6->Draw()
    
    canvas1->SaveAs(outputFilename.c_str());
    canvas1->SaveAs(outputFilename1.c_str());
 
}

    
    
    
        