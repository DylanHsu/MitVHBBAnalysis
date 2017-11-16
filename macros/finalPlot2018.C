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

#include <cassert>

#include "vhbbPlot.h"

using namespace vhbbPlot;

TCanvas* finalPlot2018(
  vhbbPlot::selectionType selType,
  TString plotTreeFileName,
  bool isBlinded,
  TString varexp,
  TString selection,
  int nbins,
  double xmin,
  double xmax,
  TString xlabel="",
  TString units="",
  TString extraString=""
) {
  system("mkdir -p MitVHBBAnalysis/plots");
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ"); assert(plotTreeFile && plotTreeFile->IsOpen());
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);

  TH1D *histos[nPlotCategories], *hTotalBkg=0;
  for(int iCat=kPlotData; iCat!=nPlotCategories; iCat++) {
    plotCategory i = static_cast<plotCategory>(iCat);
    // Construct histograms
    TString plotTitle=plotNames[i];
    TString histoName=Form("h%s",plotBaseNames[i].Data());
    if     (i==kPlotVH && selType<=kWHSR  ) plotTitle="WH(125)";
    else if(i==kPlotVH && selType<=kZnnHSR) plotTitle="ZH(125)";
    else if(i==kPlotVH && selType<=kZllHSR) plotTitle="ZH(125)";
    //histos[i] = new TH1D(histoName, plotTitle, nbins, xmin, xmax); histos[i]->Sumw2(); histos[i]->SetDirectory(0);

    // Fill histograms
    if(isBlinded && i==kPlotData) continue;
    TString cutString = Form("(theCategory==%d)*weight*(%s)",i,selection!=""?selection.Data():"1");
    Long64_t nentries = plotTree->Draw(
      Form("%s >> %s(%d,%f,%f)", varexp.Data(), histoName.Data(), nbins, xmin, xmax), 
      cutString.Data(),
      "e goff"
    );
    printf("%s %s: %lld entries\n", varexp.Data(), cutString.Data(), nentries);
    histos[i]=(TH1D*)gDirectory->Get(histoName); histos[i]->SetDirectory(0); histos[i]->SetTitle(plotTitle);
    
    // Scaling
    if(i!=kPlotData) histos[i]->Scale(theLumi);
    if(i==kPlotQCD) histos[i]->Scale(1.0);
    if(i==kPlotVH && selType!=kWHSR && selType!=kZnnHSR && selType!=kZllHSR) {
      histos[i]->Scale(100.);
      histos[i]->SetTitle(plotTitle+"x100");
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

  // Start plotting stuff
  gStyle->SetOptStat(0);
  TH1D *hRatio = (TH1D*)histos[kPlotData]->Clone("hRatio");
  hRatio->SetDirectory(0);
  for(int nb=1; nb<=hRatio->GetNbinsX(); nb++) {
    hRatio->SetBinContent( nb, hRatio->GetBinContent(nb) / hTotalBkg->GetBinContent(nb));
    hRatio->SetBinError( nb, histos[kPlotData]->GetBinError(nb) / histos[kPlotData]->GetBinContent(nb) / hRatio->GetBinContent(nb)); 
  } 
 
  TH1D *hErrorBand = (TH1D*)hTotalBkg->Clone("hErrorBand");
  hErrorBand->SetFillColor(kBlack); hErrorBand->SetFillStyle(3254);
  hErrorBand->SetMarkerColor(kBlack); hErrorBand->SetMarkerSize(0);

  TH1D *hRatioBand = (TH1D*)hErrorBand->Clone("hRatioBand");
  for(int nb=1; nb<=hRatio->GetNbinsX(); nb++) {
    hRatioBand->SetBinContent(nb,1);
    hRatioBand->SetBinError(nb, hTotalBkg->GetBinError(nb) / hTotalBkg->GetBinContent(nb) / hRatio->GetBinContent(nb));
  }
  THStack *hs = new THStack("hs",selectionNames[selType]);
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
  hs->Add(histos[ kPlotQCD  ] );   
 
  TCanvas *canvas = new TCanvas("canvas","canvas", 600, 480);
  canvas->SetTopMargin(0.0); 
  canvas->SetBottomMargin(0); 
  canvas->SetRightMargin(0.02);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.,1);
  pad1->SetTopMargin(0.1);
  pad1->SetLeftMargin(0.15);
  pad1->SetRightMargin(0.04);
  pad1->SetBottomMargin(0.03); 
  pad1->SetGridx();         
  pad1->Draw();             
  pad1->cd();              
  hs->Draw("HIST");
  hs->GetXaxis()->SetLabelSize(0);
  hs->GetXaxis()->SetTitle("");
  hs->GetYaxis()->SetTitleOffset(1.2);
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitle(Form("Events / %.1f %s",(xmax-xmin)/(float)nbins, units.Data()));
  hs->GetYaxis()->SetLabelSize(0.05);
  hs->SetMaximum( 1.4*TMath::Max( hs->GetMaximum(), histos[kPlotData]->GetMaximum() ));
  histos[kPlotVH]->Draw("HIST SAME");
  hErrorBand->Draw("E2 same");
  if(!isBlinded) histos[kPlotData]->Draw("P E0 SAME");
  TLegend *legend1=new TLegend(0.47,0.50,0.71,0.88);
  TLegend *legend2=new TLegend(0.71,0.55,0.95,.88);
  legend1->AddEntry(histos[ kPlotData ], histos[ kPlotData ]->GetTitle(),"lp");
  legend1->AddEntry(histos[ kPlotVH   ], histos[ kPlotVH   ]->GetTitle(),"lp");
  legend1->AddEntry(histos[ kPlotQCD  ], histos[ kPlotQCD  ]->GetTitle(),"f");
  legend1->AddEntry(histos[ kPlotTT   ], histos[ kPlotTT   ]->GetTitle(),"f");
  legend1->AddEntry(histos[ kPlotTop  ], histos[ kPlotTop  ]->GetTitle(),"f");
  legend1->AddEntry(histos[ kPlotVZbb ], histos[ kPlotVZbb ]->GetTitle(),"f");
  legend1->AddEntry(histos[ kPlotVVLF ], histos[ kPlotVVLF ]->GetTitle(),"f");
  legend2->AddEntry(histos[ kPlotZLF  ], histos[ kPlotZLF  ]->GetTitle(),"f");
  legend2->AddEntry(histos[ kPlotZb   ], histos[ kPlotZb   ]->GetTitle(),"f");
  legend2->AddEntry(histos[ kPlotZbb  ], histos[ kPlotZbb  ]->GetTitle(),"f");
  legend2->AddEntry(histos[ kPlotWLF  ], histos[ kPlotWLF  ]->GetTitle(),"f");
  legend2->AddEntry(histos[ kPlotWb   ], histos[ kPlotWb   ]->GetTitle(),"f");
  legend2->AddEntry(histos[ kPlotWbb  ], histos[ kPlotWbb  ]->GetTitle(),"f");
  legend1->SetFillColorAlpha(kWhite, .5); legend2->SetFillColorAlpha(kWhite, 0.5);
  legend1->SetBorderSize(0); legend2->SetBorderSize(0);
  legend1->Draw("same"); legend2->Draw("SAME");
  canvas->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.15);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin(0.04);
  pad2->SetGridx(); 
  pad2->Draw();
  pad2->cd();
  hRatioBand->SetTitle("");
  hRatioBand->GetXaxis()->SetTitleSize(0.16);
  hRatioBand->GetXaxis()->SetTitleOffset(1.05);
  hRatioBand->GetYaxis()->SetNdivisions(503);
  hRatioBand->GetXaxis()->SetLabelSize(0.2);
  hRatioBand->GetXaxis()->SetTitle( units!=""? Form("%s [%s]", xlabel.Data(), units.Data()) : xlabel);
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
  canvas->Print(Form("MitVHBBAnalysis/plots/%s_%s%s.pdf",selectionNames[selType].Data(),varexp.Data(), extraString.Data()));
  system(Form("gs -sDEVICE=png16m -dTextAlphaBits=4 -g1800x1440 -dUseCropBox -dFIXEDMEDIA -dPDFFitPage -o MitVHBBAnalysis/plots/%s_%s%s.png MitVHBBAnalysis/plots/%s_%s%s.pdf >/dev/null 2>&1",selectionNames[selType].Data(),varexp.Data(),extraString.Data(),selectionNames[selType].Data(),varexp.Data(),extraString.Data()));
  return canvas;
}

