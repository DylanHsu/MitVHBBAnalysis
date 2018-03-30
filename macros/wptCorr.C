#include "vhbbPlot.h"
#include <TBranch.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TFile.h> 
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h> 
#include <TMath.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TString.h> 
#include <TStyle.h> 
#include <TTree.h> 

//RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
//#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
//#include "RooExtendPdf.h"
//#include "RooMCStudy.h"
#include "RooGenericPdf.h"
using namespace vhbbPlot;
using namespace RooFit;

//const float sf_TT=0.9,sf_WLF=1.2,sf_Wb=1.5,sf_Wbb=1.5;
const float sf_TT=0.91,sf_WLF=1.14,sf_Wb=1.66,sf_Wbb=1.49;
const int nCat=5;// 0=data, 1=ttbar, 2=W+LF, 3=W+HF and singletop, 4=other
const int nRegion=3; // 0=W+LF,1=ttbar,2=W+HF
vector<TString> regionName = {selectionNames[kWHLightFlavorCR], selectionNames[kWH2TopCR], selectionNames[kWHHeavyFlavorCR]};
const int nBinsPt=30;
const float ptMin=100, ptMax=400;

void wptCorr(
  TString inputFileName, // inclusive resolved file
  TString outputFileName,
  int theLepSel=-1,
  int prescale=1
) {
  
  system("mkdir -p MitVHBBAnalysis/plots/wptCorr");
  TFile *plotTreeFile = TFile::Open(inputFileName, "READ");
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);
  unsigned char typeLepSel, theCategory;
  unsigned selectionBits, selectionBits_jesUp, selectionBits_jesDown, nMinusOneBits;
  int nJet, nSoft2, nSoft5, nSoft10, nIsojet;
  float sumEtSoft1;
  float pfmet, pfmetphi, pfmetsig;
  float lepton1Pt, lepton1Eta, lepton1Phi, lepton1RelIso;
  float topWBosonPt,topWBosonPt_jesUp,topWBosonPt_jesDown;
  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_QCDr1f2, weight_QCDr1f5, weight_QCDr2f1, weight_QCDr2f2, weight_QCDr5f1, weight_QCDr5f5, weight_lepSFUp;

  plotTree->SetBranchStatus("*",0);
  plotTree->SetBranchStatus("selectionBits"          , 1); plotTree->SetBranchAddress("selectionBits"          , &selectionBits         ); 
  plotTree->SetBranchStatus("selectionBits_jesUp"    , 1); plotTree->SetBranchAddress("selectionBits_jesUp"    , &selectionBits_jesUp   ); 
  plotTree->SetBranchStatus("selectionBits_jesDown"  , 1); plotTree->SetBranchAddress("selectionBits_jesDown"  , &selectionBits_jesDown ); 
  plotTree->SetBranchStatus("theCategory"            , 1); plotTree->SetBranchAddress("theCategory"            , &theCategory           ); 
  plotTree->SetBranchStatus("typeLepSel"             , 1); plotTree->SetBranchAddress("typeLepSel"             , &typeLepSel            ); 
  plotTree->SetBranchStatus("weight"                 , 1); plotTree->SetBranchAddress("weight"                 , &weight                ); 
  //plotTree->SetBranchStatus("weight_pdfUp"           , 1); plotTree->SetBranchAddress("weight_pdfUp"           , &weight_pdfUp          ); 
  //plotTree->SetBranchStatus("weight_pdfDown"         , 1); plotTree->SetBranchAddress("weight_pdfDown"         , &weight_pdfDown        ); 
  //plotTree->SetBranchStatus("weight_QCDr1f2"         , 1); plotTree->SetBranchAddress("weight_QCDr1f2"         , &weight_QCDr1f2        ); 
  //plotTree->SetBranchStatus("weight_QCDr1f5"         , 1); plotTree->SetBranchAddress("weight_QCDr1f5"         , &weight_QCDr1f5        ); 
  //plotTree->SetBranchStatus("weight_QCDr2f1"         , 1); plotTree->SetBranchAddress("weight_QCDr2f1"         , &weight_QCDr2f1        ); 
  //plotTree->SetBranchStatus("weight_QCDr2f2"         , 1); plotTree->SetBranchAddress("weight_QCDr2f2"         , &weight_QCDr2f2        ); 
  //plotTree->SetBranchStatus("weight_QCDr5f1"         , 1); plotTree->SetBranchAddress("weight_QCDr5f1"         , &weight_QCDr5f1        ); 
  //plotTree->SetBranchStatus("weight_QCDr5f5"         , 1); plotTree->SetBranchAddress("weight_QCDr5f5"         , &weight_QCDr5f5        ); 
  plotTree->SetBranchStatus("topWBosonPt"            , 1); plotTree->SetBranchAddress("topWBosonPt"            , &topWBosonPt           ); 
  //plotTree->SetBranchStatus("topWBosonPt_jesUp"      , 1); plotTree->SetBranchAddress("topWBosonPt_jesUp"      , &topWBosonPt_jesUp     ); 
  //plotTree->SetBranchStatus("topWBosonPt_jesDown"    , 1); plotTree->SetBranchAddress("topWBosonPt_jesDown"    , &topWBosonPt_jesDown   ); 
  
  // Create histograms
  TH1::SetDefaultSumw2();
  TH1F *histo_nominal  [nRegion][nCat];
  TH1F *histo_pdfUp    [nRegion][nCat];
  TH1F *histo_pdfDown  [nRegion][nCat];
  TH1F *histo_QCDr1f2  [nRegion][nCat];
  TH1F *histo_QCDr1f5  [nRegion][nCat];
  TH1F *histo_QCDr2f1  [nRegion][nCat];
  TH1F *histo_QCDr2f2  [nRegion][nCat];
  TH1F *histo_QCDr5f1  [nRegion][nCat];
  TH1F *histo_QCDr5f5  [nRegion][nCat];
  TH1F *histo_jesUp    [nRegion][nCat];
  TH1F *histo_jesDown  [nRegion][nCat];
  for(int i=0;i<nRegion;i++) for(int j=0; j<nCat; j++) {
    histo_nominal [i][j] = new TH1F(Form("histo%d_%s_nominal", j, regionName[i].Data()),Form("histo%d_%s_nominal", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_nominal  [i][j]->SetDirectory(0);
    histo_pdfUp   [i][j] = new TH1F(Form("histo%d_%s_pdfUp"  , j, regionName[i].Data()),Form("histo%d_%s_pdfUp"  , j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_pdfUp    [i][j]->SetDirectory(0);
    histo_pdfDown [i][j] = new TH1F(Form("histo%d_%s_pdfDown", j, regionName[i].Data()),Form("histo%d_%s_pdfDown", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_pdfDown  [i][j]->SetDirectory(0);
    histo_QCDr1f2 [i][j] = new TH1F(Form("histo%d_%s_QCDr1f2", j, regionName[i].Data()),Form("histo%d_%s_QCDr1f2", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_QCDr1f2  [i][j]->SetDirectory(0);
    histo_QCDr1f5 [i][j] = new TH1F(Form("histo%d_%s_QCDr1f5", j, regionName[i].Data()),Form("histo%d_%s_QCDr1f5", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_QCDr1f5  [i][j]->SetDirectory(0);
    histo_QCDr2f1 [i][j] = new TH1F(Form("histo%d_%s_QCDr2f1", j, regionName[i].Data()),Form("histo%d_%s_QCDr2f1", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_QCDr2f1  [i][j]->SetDirectory(0);
    histo_QCDr2f2 [i][j] = new TH1F(Form("histo%d_%s_QCDr2f2", j, regionName[i].Data()),Form("histo%d_%s_QCDr2f2", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_QCDr2f2  [i][j]->SetDirectory(0);
    histo_QCDr5f1 [i][j] = new TH1F(Form("histo%d_%s_QCDr5f1", j, regionName[i].Data()),Form("histo%d_%s_QCDr5f1", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_QCDr5f1  [i][j]->SetDirectory(0);
    histo_QCDr5f5 [i][j] = new TH1F(Form("histo%d_%s_QCDr5f5", j, regionName[i].Data()),Form("histo%d_%s_QCDr5f5", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_QCDr5f5  [i][j]->SetDirectory(0);
    histo_jesUp   [i][j] = new TH1F(Form("histo%d_%s_jesUp"  , j, regionName[i].Data()),Form("histo%d_%s_jesUp"  , j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_jesUp    [i][j]->SetDirectory(0);
    histo_jesDown [i][j] = new TH1F(Form("histo%d_%s_jesDown", j, regionName[i].Data()),Form("histo%d_%s_jesDown", j, regionName[i].Data()),nBinsPt,ptMin,ptMax); histo_jesDown  [i][j]->SetDirectory(0);
  }
 
  // begin plot tree loop
  Long64_t nentries = plotTree->GetEntries();
  Long64_t oneTenth = nentries/10;
  unsigned selection = vhbbPlot::kWHLightFlavorCR | vhbbPlot::kWH2TopCR | vhbbPlot::kWHHeavyFlavorCR;
  for(Long64_t ientry=0; ientry<nentries; ientry++) {
    if(ientry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries); 
    if(ientry%prescale!=0) continue;
    plotTree->GetBranch("typeLepSel")->GetEntry(ientry); if(theLepSel!=-1 && typeLepSel!=theLepSel) continue;
    plotTree->GetBranch("selectionBits")->GetEntry(ientry);
    plotTree->GetBranch("selectionBits_jesUp")->GetEntry(ientry);
    plotTree->GetBranch("selectionBits_jesDown")->GetEntry(ientry);
    bool passFullSel         = (selectionBits & selection)!=0;
    bool passFullSel_jesUp   = (selectionBits_jesUp & selection)!=0; 
    bool passFullSel_jesDown = (selectionBits_jesDown & selection)!=0; 
    //printf("selectionBits=%d, passFullSel=%d\n", selectionBits,passFullSel);
    if(!passFullSel && !passFullSel_jesUp && !passFullSel_jesDown) continue;
    
    // pick category
    plotTree->GetBranch("theCategory")->GetEntry(ientry);
    int iCat,iReg;
    switch(theCategory) {
      case kPlotData:
        iCat=0; break;
      case kPlotTT   :
        iCat=1; break;
      case kPlotWLF  :
        iCat=2; break;
      case kPlotWbb  :
      case kPlotWb   :
      case kPlotTop  :
        iCat=3; break;
      case kPlotVZbb :
      case kPlotVVLF :
      case kPlotZbb  :
      case kPlotZb   :
      case kPlotZLF  :
        iCat=4; break; 
      case kPlotQCD  :
      default:
        iCat=-1; break;
    }
    if(iCat<0 || iCat>=nCat) continue;
    if     ((selectionBits & kWHLightFlavorCR)!=0) iReg=0;
    else if((selectionBits & kWH2TopCR       )!=0) iReg=1;
    else if((selectionBits & kWHHeavyFlavorCR)!=0) iReg=2;
    else continue;
    
    plotTree->GetEntry(ientry);

    // Apply overall normalization SFs from the fit
    if     (theCategory==kPlotTT ) weight*=sf_TT;
    else if(theCategory==kPlotWLF) weight*=sf_WLF;
    else if(theCategory==kPlotWb ) weight*=sf_Wb;
    else if(theCategory==kPlotWbb) weight*=sf_Wbb;

    if(passFullSel) {
      histo_nominal[iReg][iCat]->Fill(topWBosonPt, weight);
      //printf("filling region %d cat %d\n",iReg,iCat);
    }
      /*if(!nominalOnly) {
        if(theCategory!=kPlotData) {
          weight_cmvaJESUpAllBins     =weight; weight_cmvaJESDownAllBins     =weight;
          weight_cmvaLFUpAllBins      =weight; weight_cmvaLFDownAllBins      =weight;
          weight_cmvaHFUpAllBins      =weight; weight_cmvaHFDownAllBins      =weight;
          weight_cmvaHFStats1UpAllBins=weight; weight_cmvaHFStats1DownAllBins=weight;
          weight_cmvaHFStats2UpAllBins=weight; weight_cmvaHFStats2DownAllBins=weight;
          weight_cmvaLFStats1UpAllBins=weight; weight_cmvaLFStats1DownAllBins=weight;
          weight_cmvaLFStats2UpAllBins=weight; weight_cmvaLFStats2DownAllBins=weight;
          weight_cmvaCErr1UpAllBins   =weight; weight_cmvaCErr1DownAllBins   =weight;
          weight_cmvaCErr2UpAllBins   =weight; weight_cmvaCErr2DownAllBins   =weight;
          if(weight!=0) for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
            weight_cmvaJESUpAllBins        *= weight_cmvaJESUp[iPt][iEta]       /weight;
            weight_cmvaLFUpAllBins         *= weight_cmvaLFUp[iPt][iEta]        /weight;
            weight_cmvaHFUpAllBins         *= weight_cmvaHFUp[iPt][iEta]        /weight;
            weight_cmvaHFStats1UpAllBins   *= weight_cmvaHFStats1Up[iPt][iEta]  /weight;
            weight_cmvaHFStats2UpAllBins   *= weight_cmvaHFStats2Up[iPt][iEta]  /weight;
            weight_cmvaLFStats1UpAllBins   *= weight_cmvaLFStats1Up[iPt][iEta]  /weight;
            weight_cmvaLFStats2UpAllBins   *= weight_cmvaLFStats2Up[iPt][iEta]  /weight;
            weight_cmvaCErr1UpAllBins      *= weight_cmvaCErr1Up[iPt][iEta]     /weight;
            weight_cmvaCErr2UpAllBins      *= weight_cmvaCErr2Up[iPt][iEta]     /weight;
            weight_cmvaJESDownAllBins      *= weight_cmvaJESDown[iPt][iEta]     /weight;
            weight_cmvaLFDownAllBins       *= weight_cmvaLFDown[iPt][iEta]      /weight;
            weight_cmvaHFDownAllBins       *= weight_cmvaHFDown[iPt][iEta]      /weight;
            weight_cmvaHFStats1DownAllBins *= weight_cmvaHFStats1Down[iPt][iEta]/weight;
            weight_cmvaHFStats2DownAllBins *= weight_cmvaHFStats2Down[iPt][iEta]/weight;
            weight_cmvaLFStats1DownAllBins *= weight_cmvaLFStats1Down[iPt][iEta]/weight;
            weight_cmvaLFStats2DownAllBins *= weight_cmvaLFStats2Down[iPt][iEta]/weight;
            weight_cmvaCErr1DownAllBins    *= weight_cmvaCErr1Down[iPt][iEta]   /weight;
            weight_cmvaCErr2DownAllBins    *= weight_cmvaCErr2Down[iPt][iEta]   /weight;
          }
          histo_pdfUp           [iCat]->Fill(topWBosonPt, weight_pdfUp                   ); 
          histo_pdfDown         [iCat]->Fill(topWBosonPt, weight_pdfDown                 ); 
          histo_QCDr1f2         [iCat]->Fill(topWBosonPt, weight_QCDr1f2                 ); 
          histo_QCDr1f5         [iCat]->Fill(topWBosonPt, weight_QCDr1f5                 ); 
          histo_QCDr2f1         [iCat]->Fill(topWBosonPt, weight_QCDr2f1                 ); 
          histo_QCDr2f2         [iCat]->Fill(topWBosonPt, weight_QCDr2f2                 ); 
          histo_QCDr5f1         [iCat]->Fill(topWBosonPt, weight_QCDr5f1                 ); 
          histo_QCDr5f5         [iCat]->Fill(topWBosonPt, weight_QCDr5f5                 ); 
          histo_cmvaJESUp       [iCat]->Fill(topWBosonPt, weight_cmvaJESUpAllBins        ); 
          histo_cmvaLFUp        [iCat]->Fill(topWBosonPt, weight_cmvaLFUpAllBins         ); 
          histo_cmvaHFUp        [iCat]->Fill(topWBosonPt, weight_cmvaHFUpAllBins         ); 
          histo_cmvaHFStats1Up  [iCat]->Fill(topWBosonPt, weight_cmvaHFStats1UpAllBins   ); 
          histo_cmvaHFStats2Up  [iCat]->Fill(topWBosonPt, weight_cmvaHFStats2UpAllBins   ); 
          histo_cmvaLFStats1Up  [iCat]->Fill(topWBosonPt, weight_cmvaLFStats1UpAllBins   ); 
          histo_cmvaLFStats2Up  [iCat]->Fill(topWBosonPt, weight_cmvaLFStats2UpAllBins   ); 
          histo_cmvaCErr1Up     [iCat]->Fill(topWBosonPt, weight_cmvaCErr1UpAllBins      ); 
          histo_cmvaCErr2Up     [iCat]->Fill(topWBosonPt, weight_cmvaCErr2UpAllBins      ); 
          histo_cmvaJESDown     [iCat]->Fill(topWBosonPt, weight_cmvaJESDownAllBins      ); 
          histo_cmvaLFDown      [iCat]->Fill(topWBosonPt, weight_cmvaLFDownAllBins       ); 
          histo_cmvaHFDown      [iCat]->Fill(topWBosonPt, weight_cmvaHFDownAllBins       ); 
          histo_cmvaHFStats1Down[iCat]->Fill(topWBosonPt, weight_cmvaHFStats1DownAllBins ); 
          histo_cmvaHFStats2Down[iCat]->Fill(topWBosonPt, weight_cmvaHFStats2DownAllBins ); 
          histo_cmvaLFStats1Down[iCat]->Fill(topWBosonPt, weight_cmvaLFStats1DownAllBins ); 
          histo_cmvaLFStats2Down[iCat]->Fill(topWBosonPt, weight_cmvaLFStats2DownAllBins ); 
          histo_cmvaCErr1Down   [iCat]->Fill(topWBosonPt, weight_cmvaCErr1DownAllBins    ); 
          histo_cmvaCErr2Down   [iCat]->Fill(topWBosonPt, weight_cmvaCErr2DownAllBins    ); 
        } else {
          histo_pdfUp           [iCat]->Fill(topWBosonPt, 1); 
          histo_pdfDown         [iCat]->Fill(topWBosonPt, 1); 
          histo_QCDr1f2         [iCat]->Fill(topWBosonPt, 1); 
          histo_QCDr1f5         [iCat]->Fill(topWBosonPt, 1); 
          histo_QCDr2f1         [iCat]->Fill(topWBosonPt, 1); 
          histo_QCDr2f2         [iCat]->Fill(topWBosonPt, 1); 
          histo_QCDr5f1         [iCat]->Fill(topWBosonPt, 1); 
          histo_QCDr5f5         [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaJESUp       [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaLFUp        [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaHFUp        [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaHFStats1Up  [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaHFStats2Up  [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaLFStats1Up  [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaLFStats2Up  [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaCErr1Up     [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaCErr2Up     [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaJESDown     [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaLFDown      [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaHFDown      [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaHFStats1Down[iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaHFStats2Down[iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaLFStats1Down[iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaLFStats2Down[iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaCErr1Down   [iCat]->Fill(topWBosonPt, 1); 
          histo_cmvaCErr2Down   [iCat]->Fill(topWBosonPt, 1); 
        }
      }*/
    //if(passFullSel_jesUp && !nominalOnly) 
    //  histo_jesUp  [iCat]->Fill(topWBosonPt, weight);
    //if(passFullSel_jesDown && !nominalOnly) 
    //  histo_jesDown[iCat]->Fill(topWBosonPt, weight);
  }
  plotTreeFile->Close();
  TFile *outputFile = TFile::Open(outputFileName,"recreate");
  for(int i=0;i<nRegion;i++) for(int j=0; j<nCat; j++) {
    histo_nominal [i][j]->Write();
    histo_pdfUp   [i][j]->Write();
    histo_pdfDown [i][j]->Write();
    histo_QCDr1f2 [i][j]->Write();
    histo_QCDr1f5 [i][j]->Write();
    histo_QCDr2f1 [i][j]->Write();
    histo_QCDr2f2 [i][j]->Write();
    histo_QCDr5f1 [i][j]->Write();
    histo_QCDr5f5 [i][j]->Write();
    histo_jesUp   [i][j]->Write();
    histo_jesDown [i][j]->Write();
  }
}
void fitWptCorr(TString outputFileName) {
  TFile *outputFile=TFile::Open(outputFileName,"update"); assert(outputFile);
  TH1F *histo_nominal  [nRegion][nCat];
  TH1F *histo_pdfUp    [nRegion][nCat];
  TH1F *histo_pdfDown  [nRegion][nCat];
  TH1F *histo_QCDr1f2  [nRegion][nCat];
  TH1F *histo_QCDr1f5  [nRegion][nCat];
  TH1F *histo_QCDr2f1  [nRegion][nCat];
  TH1F *histo_QCDr2f2  [nRegion][nCat];
  TH1F *histo_QCDr5f1  [nRegion][nCat];
  TH1F *histo_QCDr5f5  [nRegion][nCat];
  TH1F *histo_jesUp    [nRegion][nCat];
  TH1F *histo_jesDown  [nRegion][nCat];
  for(int i=0;i<nRegion;i++) for(int j=0; j<nCat; j++) {
    histo_nominal [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_nominal", j, regionName[i].Data()));
    histo_pdfUp   [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_pdfUp"  , j, regionName[i].Data()));
    histo_pdfDown [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_pdfDown", j, regionName[i].Data()));
    histo_QCDr1f2 [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_QCDr1f2", j, regionName[i].Data()));
    histo_QCDr1f5 [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_QCDr1f5", j, regionName[i].Data()));
    histo_QCDr2f1 [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_QCDr2f1", j, regionName[i].Data()));
    histo_QCDr2f2 [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_QCDr2f2", j, regionName[i].Data()));
    histo_QCDr5f1 [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_QCDr5f1", j, regionName[i].Data()));
    histo_QCDr5f5 [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_QCDr5f5", j, regionName[i].Data()));
    histo_jesUp   [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_jesUp"  , j, regionName[i].Data()));
    histo_jesDown [i][j] = (TH1F*)outputFile->Get(Form("histo%d_%s_jesDown", j, regionName[i].Data()));
    histo_nominal [i][j]->SetDirectory(0);
    histo_pdfUp   [i][j]->SetDirectory(0);
    histo_pdfDown [i][j]->SetDirectory(0);
    histo_QCDr1f2 [i][j]->SetDirectory(0);
    histo_QCDr1f5 [i][j]->SetDirectory(0);
    histo_QCDr2f1 [i][j]->SetDirectory(0);
    histo_QCDr2f2 [i][j]->SetDirectory(0);
    histo_QCDr5f1 [i][j]->SetDirectory(0);
    histo_QCDr5f5 [i][j]->SetDirectory(0);
    histo_jesUp   [i][j]->SetDirectory(0);
    histo_jesDown [i][j]->SetDirectory(0);
  }
  outputFile->Close();
 
  RooRealVar WpT("WpT","W boson pT",ptMin,ptMax);
  WpT.setBins(400);
  RooRealVar *frac[nRegion][nCat-1]; // relative contribution of each process type
  RooRealVar *p0[nCat-2],*p1[nCat-2]; //fit parameters
  p0[0] = new RooRealVar("p0_ttbar" ,"p0_ttbar" ,1);//,-10,10,"");
  p0[1] = new RooRealVar("p0_WLF"   ,"p0_WLF"   ,1);//,-10,10,"");
  p0[2] = new RooRealVar("p0_WHFtop","p0_WHFtop",1);//,-10,10,"");
  p1[0] = new RooRealVar("p1_ttbar" ,"p1_ttbar" ,0,-.1,.1,"GeV^{-1}");
  p1[1] = new RooRealVar("p1_WLF"   ,"p1_WLF"   ,0,-.1,.1,"GeV^{-1}");
  p1[2] = new RooRealVar("p1_WHFtop","p1_WHFtop",0,-.1,.1,"GeV^{-1}");
  // Linear correction PDFs
  RooGenericPdf *pol1_ttbar  = new RooGenericPdf("pol1_ttbar" ,"pol1_ttbar" ,"p0_ttbar  + p1_ttbar  * (WpT - 100)",RooArgList(WpT,*p0[0],*p1[0]));
  RooGenericPdf *pol1_WLF    = new RooGenericPdf("pol1_WLF"   ,"pol1_WLF"   ,"p0_WLF    + p1_WLF    * (WpT - 100)",RooArgList(WpT,*p0[1],*p1[1]));
  RooGenericPdf *pol1_WHFtop = new RooGenericPdf("pol1_WHFtop","pol1_WHFtop","p0_WHFtop + p1_WHFtop * (WpT - 100)",RooArgList(WpT,*p0[2],*p1[2]));
      
  RooDataHist * data_RDHs[nRegion], *ttbar_RDHs[nRegion], *WLF_RDHs[nRegion], *WHFtop_RDHs[nRegion], *other_RDHs[nRegion];
  RooHistPdf *ttbar_RHPs[nRegion], *WLF_RHPs[nRegion], *WHFtop_RHPs[nRegion], *other_RHPs[nRegion];
  RooCategory sample("sample",""); // Define categories
  RooProdPdf *correctedBkg[nRegion][nCat-2];
  RooAddPdf *modelBkg[nRegion];
  float sumBkgs[nRegion];
  for(int i=0;i<nRegion;i++) {
    sample.defineType(regionName[i],i+1); // indexing starts at 1 for RooCategory types
    
    // Construct RooDataHists
    data_RDHs[i]   = new RooDataHist(Form("hist_data_%s"  ,regionName[i].Data()),Form("hist_data_%s"  ,regionName[i].Data()),RooArgSet(WpT), histo_nominal[i][0]);
    ttbar_RDHs[i]  = new RooDataHist(Form("hist_ttbar_%s" ,regionName[i].Data()),Form("hist_ttbar_%s" ,regionName[i].Data()),RooArgSet(WpT), histo_nominal[i][1]);
    WLF_RDHs[i]    = new RooDataHist(Form("hist_WLF_%s"   ,regionName[i].Data()),Form("hist_WLF_%s"   ,regionName[i].Data()),RooArgSet(WpT), histo_nominal[i][2]);
    WHFtop_RDHs[i] = new RooDataHist(Form("hist_WHFtop_%s",regionName[i].Data()),Form("hist_WHFtop_%s",regionName[i].Data()),RooArgSet(WpT), histo_nominal[i][3]);
    other_RDHs[i]  = new RooDataHist(Form("hist_other_%s" ,regionName[i].Data()),Form("hist_other_%s" ,regionName[i].Data()),RooArgSet(WpT), histo_nominal[i][4]);
    
    // Construct RooHistPdfs
    ttbar_RHPs[i]  = new RooHistPdf(Form("ttbar_%s" ,regionName[i].Data()),Form("ttbar_%s" ,regionName[i].Data()),WpT,*ttbar_RDHs[i] ,1);
    WLF_RHPs[i]    = new RooHistPdf(Form("WLF_%s"   ,regionName[i].Data()),Form("WLF_%s"   ,regionName[i].Data()),WpT,*WLF_RDHs[i]   ,1);
    WHFtop_RHPs[i] = new RooHistPdf(Form("WHFtop_%s",regionName[i].Data()),Form("WHFtop_%s",regionName[i].Data()),WpT,*WHFtop_RDHs[i],1);
    other_RHPs[i]  = new RooHistPdf(Form("other_%s" ,regionName[i].Data()),Form("other_%s" ,regionName[i].Data()),WpT,*other_RDHs[i] ,1);
    sumBkgs[i]=0;
    for(int j=1;j<nCat;j++)
      sumBkgs[i]+=histo_nominal[i][j]->Integral();
    frac[i][0] = new RooRealVar(Form("frac_ttbar_%s" ,regionName[i].Data()),Form("frac_ttbar_%s" ,regionName[i].Data()),histo_nominal[i][1]->Integral(1,nBinsPt)/sumBkgs[i]);//,.8*histo_nominal[i][1]->Integral()/sumBkgs[i],1.2*histo_nominal[i][1]->Integral()/sumBkgs[i]);
    frac[i][1] = new RooRealVar(Form("frac_WLF_%s"   ,regionName[i].Data()),Form("frac_WLF_%s"   ,regionName[i].Data()),histo_nominal[i][2]->Integral(1,nBinsPt)/sumBkgs[i]);//,.8*histo_nominal[i][2]->Integral()/sumBkgs[i],1.2*histo_nominal[i][2]->Integral()/sumBkgs[i]);
    frac[i][2] = new RooRealVar(Form("frac_WHFtop_%s",regionName[i].Data()),Form("frac_WHFtop_%s",regionName[i].Data()),histo_nominal[i][3]->Integral(1,nBinsPt)/sumBkgs[i]);//,.8*histo_nominal[i][3]->Integral()/sumBkgs[i],1.2*histo_nominal[i][3]->Integral()/sumBkgs[i]);
    frac[i][3] = new RooRealVar(Form("frac_other_%s" ,regionName[i].Data()),Form("frac_other_%s" ,regionName[i].Data()),histo_nominal[i][4]->Integral(1,nBinsPt)/sumBkgs[i]);
    correctedBkg[i][0] = new RooProdPdf(Form("corr_ttbar_%s" ,regionName[i].Data()), Form("ttbar_%s * pol1_ttbar"  ,regionName[i].Data()), RooArgList(*ttbar_RHPs[i] , *pol1_ttbar ));
    correctedBkg[i][1] = new RooProdPdf(Form("corr_WLF_%s"   ,regionName[i].Data()), Form("WLF_%s * pol1_WLF"      ,regionName[i].Data()), RooArgList(*WLF_RHPs[i]   , *pol1_WLF   ));
    correctedBkg[i][2] = new RooProdPdf(Form("corr_WHFtop_%s",regionName[i].Data()), Form("WHFtop_%s * pol1_WHFtop",regionName[i].Data()), RooArgList(*WHFtop_RHPs[i], *pol1_WHFtop));
    modelBkg[i] = new RooAddPdf(Form("model_%s",regionName[i].Data()),Form("model_%s",regionName[i].Data()),
      RooArgList(*correctedBkg[i][0],*correctedBkg[i][1],*correctedBkg[i][2],*other_RHPs[i]),
      RooArgList(*frac[i][0],*frac[i][1],*frac[i][2],*frac[i][3]));
  }
 
  RooDataHist *dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(WpT),RooFit::Index(sample),
    RooFit::Import(regionName[0],*((RooDataHist*)data_RDHs[0])),
    RooFit::Import(regionName[1],*((RooDataHist*)data_RDHs[1]))//,
    //RooFit::Import(regionName[2],*((RooDataHist*)data_RDHs[2]))
  );
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelBkg[0],regionName[0]);
  totalPdf.addPdf(*modelBkg[1],regionName[1]);
  //totalPdf.addPdf(*modelBkg[2],regionName[2]);

  //WpT.setBins(400);
  RooFitResult *fitResult=0;
  fitResult = totalPdf.fitTo(*dataCombined,
    RooFit::Extended(),
    RooFit::Strategy(2),
    //RooFit::Minos(RooArgSet(eff)),
    RooFit::NumCPU(10),
    RooFit::Save());
  
  RooPlot *WpTframe[nRegion];
  TCanvas *canvas[nRegion];
  for(int i=0;i<nRegion;i++) {
    WpTframe[i] = WpT.frame(Bins(nBinsPt));
    WpTframe[i]->SetTitle(Form("WpT correction in %s",regionName[i].Data()));
    data_RDHs[i]->plotOn(WpTframe[i],MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
    modelBkg[i]->plotOn(WpTframe[i], LineStyle(kDashed),DrawOption("l"),LineColor(kRed));
    //modelBkg[i]->plotOn(WpTframe[i],Components(Form("corr_ttbar_%s" ,regionName[i].Data())),VisualizeError(*fitResult,1,kFALSE), FillStyle(3005), FillColor(kRed));
    modelBkg[i]->plotOn(WpTframe[i],Components(Form("corr_ttbar_%s" ,regionName[i].Data())),LineStyle(kSolid),LineColor(plotColors[kPlotTT]));
    modelBkg[i]->plotOn(WpTframe[i],Components(Form("corr_WLF_%s" ,regionName[i].Data())),LineStyle(kSolid),LineColor(plotColors[kPlotWLF]));
    modelBkg[i]->plotOn(WpTframe[i],Components(Form("corr_WHFtop_%s" ,regionName[i].Data())),LineStyle(kSolid),LineColor(plotColors[kPlotWbb]));
    canvas[i]=new TCanvas;
    WpTframe[i]->Draw();
  }
  /*
  TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  
  std::vector<TH1F**> wptSpectra    ; wptSpectra    .reserve(64);
  std::vector<TF1*  > wptCorrections; wptCorrections.reserve(64);
  wptSpectra.push_back(histo_nominal         );
  if(!nominalOnly) {
    wptSpectra.push_back(histo_pdfUp           );
    wptSpectra.push_back(histo_pdfDown         );
    wptSpectra.push_back(histo_QCDr1f2         );
    wptSpectra.push_back(histo_QCDr1f5         );
    wptSpectra.push_back(histo_QCDr2f1         );
    wptSpectra.push_back(histo_QCDr2f2         );
    wptSpectra.push_back(histo_QCDr5f1         );
    wptSpectra.push_back(histo_QCDr5f5         );
    wptSpectra.push_back(histo_cmvaJESUp       );
    wptSpectra.push_back(histo_cmvaLFUp        );
    wptSpectra.push_back(histo_cmvaHFUp        );
    wptSpectra.push_back(histo_cmvaHFStats1Up  );
    wptSpectra.push_back(histo_cmvaHFStats2Up  );
    wptSpectra.push_back(histo_cmvaLFStats1Up  );
    wptSpectra.push_back(histo_cmvaLFStats2Up  );
    wptSpectra.push_back(histo_cmvaCErr1Up     );
    wptSpectra.push_back(histo_cmvaCErr2Up     );
    wptSpectra.push_back(histo_cmvaJESDown     );
    wptSpectra.push_back(histo_cmvaLFDown      );
    wptSpectra.push_back(histo_cmvaHFDown      );
    wptSpectra.push_back(histo_cmvaHFStats1Down);
    wptSpectra.push_back(histo_cmvaHFStats2Down);
    wptSpectra.push_back(histo_cmvaLFStats1Down);
    wptSpectra.push_back(histo_cmvaLFStats2Down);
    wptSpectra.push_back(histo_cmvaCErr1Down   );
    wptSpectra.push_back(histo_cmvaCErr2Down   );
    wptSpectra.push_back(histo_jesUp           );
    wptSpectra.push_back(histo_jesDown         );
  }
  
  // Derive the WpT corrections for each set of spectra
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  for(unsigned iS=0; iS<wptSpectra.size(); iS++) {
    TString name= wptSpectra[iS][0]->GetName();
    TString theVariation(name(7,name.Length()-7));
    wptSpectra[iS][3]->Add(wptSpectra[iS][0]);
    wptSpectra[iS][3]->Add(wptSpectra[iS][2],-1);
    wptSpectra[iS][3]->Divide(wptSpectra[iS][1]);
    for(int nb=1; nb<=nBinsPt; nb++) {
      float nData=wptSpectra[iS][0]->GetBinContent(nb);
      float nW   =wptSpectra[iS][1]->GetBinContent(nb);
      float nBkg =wptSpectra[iS][2]->GetBinContent(nb);
      //printf("nb=%d, nData=%f, nW=%f, nBkg=%f\n", nb, nData, nW, nBkg);
    }
    wptSpectra[iS][3]->SetTitle(Form("Correction to W p_{T} spectrum (%s)", theVariation.Data()));
    
    TCanvas *cCorr=new TCanvas("cCorr", "cCorr",1024,640);
    cCorr->SetRightMargin(0.3);
    wptSpectra[iS][3]->GetXaxis()->SetTitle("Reconstructed W p_{T} [GeV]");
    wptSpectra[iS][3]->GetXaxis()->SetMoreLogLabels();
    wptSpectra[iS][3]->GetYaxis()->SetTitle("(Data - Bkg.) / W+jets");
    wptSpectra[iS][3]->SetMarkerStyle(20);
    wptSpectra[iS][3]->SetMarkerSize(0.8);
    wptSpectra[iS][3]->Draw("P E1");
    wptSpectra[iS][3]->GetYaxis()->SetRangeUser(0.2,1.8);
    //wptSpectra[iS][3]->GetXaxis()->SetRangeUser(hTotal_met->GetXaxis()->GetBinLowEdge(1), hTotal_met->GetXaxis()->GetBinUpEdge(hTotal_met->GetNbinsX())-.001);
    //TLine *oneline = new TLine(hTotal_met->GetXaxis()->GetBinLowEdge(1), 1, hTotal_met->GetXaxis()->GetBinUpEdge(hTotal_met->GetNbinsX()), 1);
    TF1 *theCorr; TString corrName=Form("wptCorr_%s",theVariation.Data());
    if     (fitFunction==1) { theCorr = new TF1(corrName, "pol5", ptMin, ptMax); theCorr->SetLineColor(kGreen); } 
    else if(fitFunction==2) { theCorr = new TF1(corrName, "pol1", ptMin, ptMax); theCorr->SetLineColor(kMagenta);}
    if(fitFunction==1) {
    } else if(fitFunction==2) {
      //theCorr->SetParNames("k", "E_{c}","#varepsilon_{obs}");
      //theCorr->SetParameter(0,1);  theCorr->SetParLimits(0,0.5,1.5);
      //theCorr->SetParameter(1,0);  theCorr->SetParLimits(1,-0.1,0.1);
      //theCorr->SetParameter(2,0);  theCorr->SetParLimits(2,-0.1,0.1); 
      //theCorr->SetParameter(2,1);  theCorr->SetParLimits(2,1,1); 
    } else {
      printf("bad fit function\n");
      return;
    }
    TFitResultPtr theFitResultPtr = wptSpectra[iS][3]->Fit(theCorr, "EMS","", ptMin, ptMax);
    TFitResult *theFitResult = theFitResultPtr.Get();
    TPaveText *formula = new TPaveText(0.72,0.8,0.96,0.9, "NDC");
    formula->SetFillColorAlpha(0,0.4);
    formula->SetFillStyle(1001);
    formula->SetLineColor(kBlack);
    formula->SetBorderSize(1);
    if     (fitFunction==1) formula->AddText("f(x) = #Sigma_{n=0}^{5} p_{n} x^{n}");
    else if(fitFunction==2) formula->AddText("f(x) = p_{0} + p_{1} x");
    formula->Draw("SAME");
    TPaveText *selectionPave = new TPaveText(0.12,0.77,.4,0.88, "NDC");
    selectionPave->SetFillColorAlpha(0,0.4);
    selectionPave->SetFillStyle(1001);
    selectionPave->SetTextColor(14);
    selectionPave->SetBorderSize(0);
    //if     (mode=="monojet") selectionPave->AddText("#bf{#it{#splitline{\"Monojet\" topology}{Central jet, p_{T}>100, r_{CH}>10%, r_{NH}<80%}}}");
    //else if(mode=="vbf"    ) selectionPave->AddText("#bf{#it{#splitline{VBF topology}{Two jets p_{T}>80(40), leading jet r_{CH}>10%, r_{NH}<80% if central}}}");
    //selectionPave->Draw("SAME");
    cCorr->Modified(); cCorr->Update(); 
    TPaveStats *ps;
    ps  = (TPaveStats *)cCorr->GetPrimitive("stats");
    if(!ps) { printf("bad pointer\n"); return; }
    ps->SetX1NDC(0.72);
    ps->SetX2NDC(0.96);
    ps->SetY1NDC(0.4);
    ps->SetY2NDC(0.78);
    cCorr->Modified(); cCorr->Update(); 
    cCorr->Print(Form("MitVHBBAnalysis/plots/wptCorr/wptCorr_%s.pdf",theVariation.Data()));
    cCorr->Print(Form("MitVHBBAnalysis/plots/wptCorr/wptCorr_%s.png",theVariation.Data()));
    //wptCorrections.push_back(theCorr);
    theCorr->Write(theCorr->GetName());
    TH1F *theCorrError = (TH1F*)theCorr->GetHistogram();
    for(int i=1; i<=theCorrError->GetNbinsX(); i++) {
      double x[1]={ theCorrError->GetBinCenter(i) };
      double err[1];
      theFitResult->GetConfidenceIntervals(1,1,1,x,err,0.683,false);
      theCorrError->SetBinError(i,err[0]);
    }
    theCorrError->Write(Form("histWithStatUnc_%s",theCorr->GetName()));
  }
  //TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  //for(unsigned corr=0; corr<wptCorrections.size(); corr++) {
  //  wptCorrections[corr]->Write(wptCorrections[corr]->GetName());
  //}
  */

}


