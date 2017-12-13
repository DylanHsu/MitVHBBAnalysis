#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TList.h>

#include <cassert>

#include "vhbbPlot.h"

using namespace vhbbPlot;

//TCanvas* finalPlot2018(
TList *finalPlot2018(
  vhbbPlot::selectionType selType, 
  TString inputFileName,
  TString histoName="", 
  TString plotTitle="",
  TString extraText="",
  bool isBlinded=false
) {
  const bool plotQCD=false;
  system("mkdir -p MitVHBBAnalysis/plots");
  TFile *inputFile = TFile::Open(inputFileName, "READ"); assert(inputFile);
  string rawName; {
    size_t lastDot = string(inputFileName).find_last_of(".");
    size_t lastSlash = string(inputFileName).find_last_of("/");
    rawName = string(inputFileName).substr(lastSlash+1, lastDot-lastSlash-1);
  }
  
  TList *listOfHistoNames=inputFile->GetListOfKeys();
  TList *listOfCanvases=new TList();
  for(unsigned iHisto=0; iHisto<=listOfHistoNames->LastIndex(); iHisto++) {
    TString theHistoName = listOfHistoNames->At(iHisto)->GetName();
    TString outPdf = Form("MitVHBBAnalysis/plots/%s_%s.pdf", rawName.c_str(), theHistoName.Data()); 
    TString outPng = Form("MitVHBBAnalysis/plots/%s_%s.png", rawName.c_str(), theHistoName.Data()); 
    TString xlabel=""; TString plotName;
    TH1D *histos[nPlotCategories], *hTotalBkg=0;
    for(int iCat=kPlotData; iCat!=nPlotCategories; iCat++) {
      plotCategory i = static_cast<plotCategory>(iCat);
      if(!plotQCD && i==kPlotQCD) continue;
      // Construct histograms
      if (i==kPlotVH) {
        if(selType==kWHSR) plotName="WH(125)";
        else if(selType==kWHLightFlavorCR || selType==kWHHeavyFlavorCR || selType==kWH2TopCR) plotName="WH(125)x100";
        else if(selType==kZnnHSR || selType==kZllHSR) plotName="ZH(125)";
        else plotName="ZH(125)x100";
      } else plotName=plotNames[i]; 

      // Fill histograms
      histos[i]=(TH1D*)gDirectory->Get(Form("%s/histo%d",theHistoName.Data(),iCat));
      histos[i]->SetDirectory(0);
      if(xlabel=="") xlabel=histos[i]->GetTitle();
      histos[i]->SetName(plotName);
      histos[i]->SetTitle(plotTitle);
      
      // Scaling
      //if(i!=kPlotData) histos[i]->Scale(theLumi);
      if(i==kPlotQCD) histos[i]->Scale(1.0);
      if(i==kPlotVH && selType!=kWHSR && selType!=kZnnHSR && selType!=kZllHSR) {
        histos[i]->Scale(100.);
      }
      
      // Colors
      if(i==kPlotData) {
        histos[i]->SetMarkerColor(plotColors[i]); 
        histos[i]->SetLineColor(plotColors[i]); 
        histos[i]->SetMarkerStyle(20);
      } else if(i==kPlotVH) {
        histos[i]->SetLineColor(plotColors[i]); 
        histos[i]->SetLineWidth(3);
        histos[i]->SetFillStyle(0);
      } else {
        histos[i]->SetLineColor(kBlack);
        histos[i]->SetFillColor(plotColors[i]);
        histos[i]->SetMarkerColor(plotColors[i]);
        histos[i]->SetFillStyle(1001);
      }

      // Summing Up
      if(!hTotalBkg) {
        hTotalBkg=(TH1D*)histos[i]->Clone("hTotal");
        hTotalBkg->Reset(); hTotalBkg->Clear();
        hTotalBkg->SetDirectory(0);
        hTotalBkg->SetTitle("hTotal");
      }
      if(hTotalBkg && i!=kPlotData && i!=kPlotVH) {
        hTotalBkg->Add(histos[i]);
      }
    }
    string units=""; {
      size_t lastOpenBrace = string(xlabel).find_last_of("[");
      size_t lastCloseBrace = string(xlabel).find_last_of("]");
      if(lastOpenBrace!=string::npos && lastCloseBrace!=string::npos) 
        units = string(xlabel).substr(lastOpenBrace+1, lastCloseBrace-lastOpenBrace-1);
    }

    // Start plotting stuff
    gStyle->SetOptStat(0);
    TH1D *hRatio = (TH1D*)histos[kPlotData]->Clone("hRatio");
    hRatio->SetDirectory(0);
    for(int nb=1; nb<=hRatio->GetNbinsX(); nb++) {
      if(hTotalBkg->GetBinContent(nb)>0 && histos[kPlotData]->GetBinContent(nb)>0) {
        hRatio->SetBinContent( nb, hRatio->GetBinContent(nb) / hTotalBkg->GetBinContent(nb));
        hRatio->SetBinError( nb, histos[kPlotData]->GetBinError(nb) / histos[kPlotData]->GetBinContent(nb) / hRatio->GetBinContent(nb)); 
      } else {
        hRatio->SetBinContent( nb, 1);
        hRatio->SetBinError(nb, 1);
      }
    } 
 
    TH1D *hErrorBand = (TH1D*)hTotalBkg->Clone("hErrorBand");
    hErrorBand->SetFillColor(kBlack); hErrorBand->SetFillStyle(3254);
    hErrorBand->SetMarkerColor(kBlack); hErrorBand->SetMarkerSize(0);

    TH1D *hRatioBand = (TH1D*)hErrorBand->Clone("hRatioBand");
    for(int nb=1; nb<=hRatio->GetNbinsX(); nb++) {
      hRatioBand->SetBinContent(nb,1);
      hRatioBand->SetBinError(nb, hTotalBkg->GetBinError(nb) / hTotalBkg->GetBinContent(nb) / hRatio->GetBinContent(nb));
    }
    THStack *hs = new THStack("hs",plotTitle!=""?plotTitle:selectionNames[selType]); assert(hs);
    hs->Add(histos[ kPlotVZbb ] );   
    hs->Add(histos[ kPlotVVLF ] );   
    hs->Add(histos[ kPlotZLF  ] );    
    hs->Add(histos[ kPlotZb   ] );   
    hs->Add(histos[ kPlotZbb  ] );   
    hs->Add(histos[ kPlotWLF  ] );   
    hs->Add(histos[ kPlotWb   ] );   
    hs->Add(histos[ kPlotWbb  ] );   
    hs->Add(histos[ kPlotTop  ] );   
    hs->Add(histos[ kPlotTT   ] );   
    if(plotQCD) hs->Add(histos[ kPlotQCD  ] );   
 
    TCanvas *canvas = new TCanvas(Form("canvas_%s",theHistoName.Data()),theHistoName,600,480);
    canvas->SetTopMargin(0.0); 
    canvas->SetBottomMargin(0); 
    canvas->SetRightMargin(0.02);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.,1);
    pad1->SetTopMargin(0.1);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.04);
    pad1->SetBottomMargin(0.03); 
    //pad1->SetGridx();         
    pad1->Draw();             
    pad1->cd();              
    hs->Draw("HIST");
    int nbins=hs->GetXaxis()->GetNbins();
    float xmin=hs->GetXaxis()->GetBinLowEdge(1);
    float xmax=hs->GetXaxis()->GetBinLowEdge(nbins+1);
    hs->GetXaxis()->SetTitle("");
    hs->GetXaxis()->SetLabelSize(0);
    hs->GetYaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetTitle(Form("Events / %.1f %s",(xmax-xmin)/(float)nbins, units.c_str()));
    hs->GetYaxis()->SetLabelSize(0.05);
    double plotMax=1.4;
    if(hTotalBkg->GetMean() > xmin + (xmax-xmin)/4.) plotMax=2.;
    hs->SetMaximum( plotMax*TMath::Max( hs->GetMaximum(), histos[kPlotData]->GetMaximum() ));
    histos[kPlotVH]->Draw("HIST SAME");
    hErrorBand->Draw("E2 same");
    if(!isBlinded) histos[kPlotData]->Draw("P E0 SAME");
    TLegend *legend1,*legend2;
    if(plotQCD) legend1=new TLegend(0.47,0.50,0.71,0.88);
    else        legend1=new TLegend(0.47,0.555,0.71,0.88); 
    legend2=new TLegend(0.71,0.555,0.95,.88);
    legend1->AddEntry(histos[ kPlotData ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotData)] ,"lp");
    if(plotQCD) legend1->AddEntry(histos[ kPlotQCD  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotQCD )] ,"f");
    legend1->AddEntry(histos[ kPlotTT   ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotTT  )] ,"f");
    legend1->AddEntry(histos[ kPlotTop  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotTop )] ,"f");
    legend1->AddEntry(histos[ kPlotWbb  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotWbb )] ,"f");
    legend1->AddEntry(histos[ kPlotWb   ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotWb  )] ,"f");
    legend1->AddEntry(histos[ kPlotWLF  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotWLF )] ,"f");
    legend2->AddEntry(histos[ kPlotZbb  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotZbb )] ,"f");
    legend2->AddEntry(histos[ kPlotZb   ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotZb  )] ,"f");
    legend2->AddEntry(histos[ kPlotZLF  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotZLF )] ,"f");
    legend2->AddEntry(histos[ kPlotVVLF ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotVVLF)] ,"f");
    legend2->AddEntry(histos[ kPlotVZbb ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotVZbb)] ,"f");
    legend2->AddEntry(histos[ kPlotVH   ], histos[kPlotVH]->GetName(), "lp");
    legend1->SetFillColorAlpha(kWhite, .5); legend2->SetFillColorAlpha(kWhite, 0.5);
    legend1->SetBorderSize(0); legend2->SetBorderSize(0);
    legend1->Draw("same"); legend2->Draw("SAME");
    canvas->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetLeftMargin(0.15);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.04);
    //pad2->SetGridx(); 
    pad2->Draw();
    pad2->cd();
    hRatioBand->SetTitle("");
    hRatioBand->GetXaxis()->SetTitleSize(0.16);
    hRatioBand->GetXaxis()->SetTitleOffset(1.05);
    hRatioBand->GetYaxis()->SetNdivisions(503);
    hRatioBand->GetXaxis()->SetLabelSize(0.2);
    hRatioBand->GetXaxis()->SetTitle(xlabel);
    hRatioBand->GetYaxis()->SetTitle("Data/MC");
    hRatioBand->GetYaxis()->CenterTitle();
    hRatioBand->GetYaxis()->SetTitleOffset(0.37);
    hRatioBand->GetYaxis()->SetTitleSize(0.12);
    hRatioBand->GetYaxis()->SetLabelSize(0.15);
    hRatioBand->SetFillStyle(3254);
    hRatioBand->SetFillColor(kBlack);
    hRatioBand->SetMinimum(.3);
    hRatioBand->SetMaximum(1.7);
    hRatioBand->Draw("E2");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.8);
    hRatio->Draw("P E0 x0 SAME");
    TLine *baseline = new TLine(xmin,1,xmax,1);
    baseline->SetLineStyle(kSolid); baseline->Draw("SAME");
    canvas->Print(outPdf);
    system(Form("gs -sDEVICE=png16m -dTextAlphaBits=4 -g1800x1440 -dUseCropBox -dFIXEDMEDIA -dPDFFitPage -o %s %s >/dev/null 2>&1",outPng.Data(), outPdf.Data()));
    listOfCanvases->Add(canvas);
  }
  return listOfCanvases;
  //return canvas;
}

