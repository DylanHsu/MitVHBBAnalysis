#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
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
  const char* varexp,
  const char* selection,
  double nbins,
  double xmin,
  double xmax
) {
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ"); assert(plotTreeFile && plotTreeFile->IsOpen());
  TTree *plotTree = (TTree*)plotTreeFile->Get("events"); assert(plotTree);

  TH1D *histos[nPlotCategories];
  for(int i=0; i<nPlotCategories; i++) {
    TString plotTitle=plotNames[i];
    TString histoName=Form("h%s",plotBaseNames[i].Data());
    if     (i==kPlotVH && selType<=kWHSR  ) plotTitle="WH(125)";
    else if(i==kPlotVH && selType<=kZnnHSR) plotTitle="ZH(125)";
    else if(i==kPlotVH && selType<=kZllHSR) plotTitle="ZH(125)";
    histos[i] = new TH1D(histoName, plotTitle, nbins, xmin, xmax); histos[i]->Sumw2();
    if(isBlinded && i==kPlotData) continue;
    plotTree->Draw(
      Form("%s >> +%s", varexp, histoName.Data()),
      selection,
      "e goff"
    );
    //histos[i]=gDirectory->Get(histoName); histos[i]->SetDirectory(0);
    if(i!=kPlotData) histos[i]->Scale(theLumi);
  }
  TCanvas *canvas;
  return canvas;
}

