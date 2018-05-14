#include "PandaAnalysis/Flat/interface/GeneralTree.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "PandaAnalysis/Utilities/src/CSVHelper.cc"
#include <Compression.h>
#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <cassert>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unistd.h>
#include "MitAnalysisRunII/panda/macros/80x/auxiliar.h"
#include "MitAnalysisRunII/panda/macros/80x/common.h"
#include "vhbbPlot.h"
#include "PandaAnalysis/Flat/interface/Common.h"
#include "formulas.h"


const bool useHtBinnedVJetsKFactor=true;
const int NJES = 1; //(int)shiftjes::N; // Number of JES variations, currently only look at the nominal for performance reason
// Location of the PA ntuples or skims
const TString ntupleDir = "/data/t3home000/dhsu/dylansVHSkims/2016/dilepton/";
//const TString ntupleDir = "/mnt/hadoop/scratch/dhsu/dylansVHSkims/zllTest/";
const int nLepSel=3; // Number of lepton selections

using namespace vhbbPlot;

// useBoostedCategory:
// False means you will use all possible events to do the resolved Z(ll)H(bb)
// True means the events with Z back to back with 250 GeV fatjet are reserved for boosted ZH

void zllhAnalysis(
  TString dataCardDir,
  bool useBoostedCategory=false,
  int MVAVarType=3,
  unsigned year=2016,
  unsigned debug=0 
) {
  
  /////////////////////////////
  // List of Samples
  vector<pair<TString,vhbbPlot::sampleType>> samples;
  samples.emplace_back("DoubleEG"           , vhbbPlot::kData   );
  samples.emplace_back("DoubleMuon"         , vhbbPlot::kData   );
//samples.emplace_back("SingleTop_tT"       , vhbbPlot::kTop    );       
//samples.emplace_back("SingleTop_tTbar"    , vhbbPlot::kTop    );       
  samples.emplace_back("SingleTop_tW"       , vhbbPlot::kTop    );       
  samples.emplace_back("SingleTop_tbarW"    , vhbbPlot::kTop    );       
  samples.emplace_back("TTTo2L2Nu"          , vhbbPlot::kTT     );       
  samples.emplace_back("WWTo2L2Nu"          , vhbbPlot::kWW     );       
//samples.emplace_back("WWTo4Q"             , vhbbPlot::kWW     );       
//samples.emplace_back("WWToLNuQQ"          , vhbbPlot::kWW     );       
//samples.emplace_back("WZTo1L1Nu2Q"        , vhbbPlot::kVZ     );       
  samples.emplace_back("WZTo2L2Q"           , vhbbPlot::kVZ     );       
  samples.emplace_back("ZZTo2L2Q"           , vhbbPlot::kVZ     );       
//samples.emplace_back("ZZTo4Q"             , vhbbPlot::kVZ     );       
  samples.emplace_back("ZJets_ht100to200"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_ht200to400"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_ht400to600"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_ht600to800"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_ht800to1200"  , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_ht1200to2500" , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_ht2500toinf"  , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_pt50to100"    , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_pt100to250"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_pt250to400"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_pt400to650"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZJets_pt650toinf"   , vhbbPlot::kZjets  );       
  samples.emplace_back("ZllHbb_mH125"       , vhbbPlot::kZH     );       
  samples.emplace_back("ggZllHbb_mH125"     , vhbbPlot::kZH     );       
  // End List of Samples
  /////////////////////////////
  
  // Load Shared Objects for ACLIC
  gSystem->Load("libPandaAnalysisFlat.so");
  
  assert(dataCardDir!="");
  if(dataCardDir!="") system(Form("mkdir -p MitVHBBAnalysis/datacards/%s",dataCardDir.Data()));
  double theLumi= (year==2016)? 35900.:41300;
  unsigned long whichTriggers = (1<<pa::kSingleEleTrig) | (1<<pa::kDoubleEleTrig) | (1<<pa::kSingleMuTrig) | (1<<pa::kDoubleMuTrig) | (1<<pa::kEMuTrig);
 // Load Pileup Weights
  TString puPath =
    (year==2016)? "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root" :
    "MitAnalysisRunII/data/90x/puWeights_90x.root"; 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data())); assert(fPUFile);
  TH1D *puWeights     = (TH1D*)(fPUFile->Get("puWeights"    )); assert(puWeights    );
  TH1D *puWeightsUp   = (TH1D*)(fPUFile->Get("puWeightsUp"  )); assert(puWeightsUp  );
  TH1D *puWeightsDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(puWeightsDown);
  puWeights    ->SetDirectory(0); 
  puWeightsUp  ->SetDirectory(0); 
  puWeightsDown->SetDirectory(0); 
  delete fPUFile;
  // Done Loading Pileup Weights

  // Choice of the MVA variable type, binning, and the name
  // This can be different for each control region
  std::map<selectionType, vector<float>> MVAbins;
  std::map<selectionType, TString> MVAVarName, shapeType; // title of histos, and simple name
  
  vector<selectionType> selections = { // List of enums for the selections we care about, defined in vhbbPlot.h
    kZllHLightFlavorCR   ,
    kZllHHeavyFlavorCR   ,
    kZllH2TopCR          ,
    kZllHSR              ,
    kZllHPresel          ,
    kZllHLightFlavorFJCR ,
    kZllHHeavyFlavorFJCR ,
    kZllHTT2bFJCR        ,
    kZllHTT1bFJCR        ,
    kZllHFJSR            ,
    kZllHFJPresel        
  };
  vector<TString> leptonStrings={"mm","ee","em"};
 
  for(unsigned iSel=0; iSel<selections.size(); iSel++) { // Define the shape variable
    selectionType sel = selections[iSel];
    if(MVAVarType==1) {
      // 1 - simple pT variable
      MVAVarName[sel]="Higgs p_{T} classifier [GeV]";
      if(sel>=kZllHLightFlavorCR && sel<=kZllHPresel) {
        MVAbins[sel]={100,120,140,160,180,200,250,300,350};
        MVAVarName[sel]="H(bb) pT";
        shapeType[sel]="ptShape";
      } else if(sel>=kZllHLightFlavorFJCR && sel<=kZllHFJPresel) {
        MVAbins[sel]={250,300,350,400,450,500,550,600};
        MVAVarName[sel]="H(bb) pT";
        shapeType[sel]="ptShape";
      } 
    //} else if(MVAVarType==2) {
      // 2 - multiclass BDT in SR, subleading CMVA in CR (not implemented)
    } else if(MVAVarType==3) {
      // 3 - normal BDT in SR, subleading CMVA in CR
      if(sel==kZllHLightFlavorCR) {
        MVAbins[sel]={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000};
        MVAVarName[sel]="Subleading H(bb) CMVA";
        shapeType[sel]="lesserCMVAShape";
      } else if(sel==kZllHHeavyFlavorCR || sel==kZllH2TopCR || sel==kZllHPresel) {
        MVAbins[sel]={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
        MVAVarName[sel]="Subleading H(bb) CMVA";
        shapeType[sel]="lesserCMVAShape";
      } else if(sel==kZllHSR) {
        MVAbins[sel]={-1,-0.5,0, 0.20,0.40,0.60,0.70,0.80,0.9,1.00};
        MVAVarName[sel]="BDT Output";
        shapeType[sel]="singleClassBDTShape"; 
      } else if((sel>=kZllHLightFlavorFJCR && sel<kZllHFJSR) || sel==kZllHFJPresel) {
        if(sel==kZllHHeavyFlavorFJCR)
          MVAbins[sel]={40,45,50,55,60,65,70,75,80};
        else
          MVAbins[sel]={40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,160,180,200};
        MVAVarName[sel]="Fatjet soft drop mass [GeV]";
        shapeType[sel]="softDropMassShape";
      } else if(sel==kZllHFJSR) {
        MVAbins[sel]={-0.6,-0.2,0,0.2,0.3,0.6};
        MVAVarName[sel]="BDT Output";
        shapeType[sel]="singleClassBDTShape"; 
      }
    } else throw std::runtime_error("bad MVAVarType");
  }
  
  // Declare histograms for plotting
  printf("Building plotting histograms, please wait...\n");
  const int nPlots=50;
  vector<float> xmin(nPlots), xmax(nPlots);  vector<int> nbins(nPlots);
  vector<TString> histoNames(nPlots), histoTitles(nPlots);
  { int p=0;
    histoNames[p]="MVAVar"                  ; histoTitles[p]=""                      ;                                                p++;
    histoNames[p]="lepton1Pt"               ; histoTitles[p]="Lepton 1 p_{T} [GeV]"  ; nbins[p]=  23; xmin[p]=    20; xmax[p]=   250; p++; 
    histoNames[p]="lepton2Pt"               ; histoTitles[p]="Lepton 2 p_{T} [GeV]"  ; nbins[p]=  23; xmin[p]=    20; xmax[p]=   250; p++; 
    histoNames[p]="lepton1Eta"              ; histoTitles[p]="Lepton 1 #eta"         ; nbins[p]=  25; xmin[p]=  -2.5; xmax[p]=   2.5; p++; 
    histoNames[p]="lepton2Eta"              ; histoTitles[p]="Lepton 2 #eta"         ; nbins[p]=  25; xmin[p]=  -2.5; xmax[p]=   2.5; p++; 
    histoNames[p]="ZBosonPt"                ; histoTitles[p]="Z boson pT [GeV]"      ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="ZBosonEta"               ; histoTitles[p]="Z boson #eta"          ; nbins[p]=  25; xmin[p]=  -2.5; xmax[p]=   2.5; p++; 
    histoNames[p]="ZBosonPhi"               ; histoTitles[p]="Z boson phi"           ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="ZBosonM"                 ; histoTitles[p]="Z boson mass"          ; nbins[p]=  60; xmin[p]=     0; xmax[p]=   120; p++; 
    histoNames[p]="ZBosonLep1CosThetaCS"    ; histoTitles[p]="Z(ll) cos#theta^{CS}"  ; nbins[p]=  20; xmin[p]=    -1; xmax[p]=    1.; p++; 
    histoNames[p]="ZBosonLep1CosThetaStar"  ; histoTitles[p]="cos#theta* Z(ll)+jj"   ; nbins[p]=  20; xmin[p]=    -1; xmax[p]=    1.; p++; 
    histoNames[p]="ZBosonLep1CosThetaStarFJ"; histoTitles[p]="cos#theta* Z(ll)+FJ"   ; nbins[p]=  20; xmin[p]=    -1; xmax[p]=    1.; p++; 
    histoNames[p]="Mjj"                     ; histoTitles[p]="Dijet mass [GeV]"      ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="pTjj"                    ; histoTitles[p]="Dijet pT [GeV]"        ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="bjet1Pt"                 ; histoTitles[p]="B-jet 1 pT [GeV]"      ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="bjet2Pt"                 ; histoTitles[p]="B-jet 2 pT [GeV]"      ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="bjet1btag"               ; histoTitles[p]="B-jet 1 btag"          ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="bjet2btag"               ; histoTitles[p]="B-jet 2 btag"          ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="nJet"                    ; histoTitles[p]="N central AK4CHS jets" ; nbins[p]=   8; xmin[p]=    0.; xmax[p]=    8.; p++; 
    histoNames[p]="deltaPhiZH"              ; histoTitles[p]="#Delta#phi(Z,H) [Rad]" ; nbins[p]=  20; xmin[p]= 1.571; xmax[p]= 3.142; p++; 
    histoNames[p]="ptBalanceZH"             ; histoTitles[p]="|H pT / Z pT|"         ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=    3.; p++; 
    histoNames[p]="sumEtSoft1"              ; histoTitles[p]="#sum E_{T}(soft 1)"    ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  300.; p++; 
    histoNames[p]="nSoft2"                  ; histoTitles[p]="N^{soft}_{2}"          ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=   25.; p++; 
    histoNames[p]="nSoft5"                  ; histoTitles[p]="N^{soft}_{5}"          ; nbins[p]=  12; xmin[p]=    0.; xmax[p]=   12.; p++; 
    histoNames[p]="nSoft10"                 ; histoTitles[p]="N^{soft}_{10}"         ; nbins[p]=   8; xmin[p]=    0.; xmax[p]=    8.; p++; 
    histoNames[p]="bdtValue"                ; histoTitles[p]="BDT Output"            ; nbins[p]=  40; xmin[p]=    -1; xmax[p]=    1.; p++; 
  }
  
  TH1F *histos[nLepSel][selections.size()][nPlots][nPlotCategories];
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(unsigned iSel=0; iSel<selections.size(); iSel++)
  for(int p=0; p<nPlots; p++) {
    selectionType sel = selections[iSel];
    if(histoNames[p]=="") continue;
    for(unsigned char ic=0; ic<nPlotCategories; ic++) {
      if(histoNames[p]=="MVAVar") {
        histos[lep][iSel][p][ic] = new TH1F(
          Form("histo%d_Z%sH%s_%s",
            ic,
            leptonStrings[lep].Data(),
            selectionNames[(int)sel].Data(),
            histoNames[p].Data()
          ), MVAVarName[sel], MVAbins[sel].size()-1, MVAbins[sel].data());
      } else {
        histos[lep][iSel][p][ic] = new TH1F(
          Form("histo%d_Z%sH%s_%s",
            ic,
            leptonStrings[lep].Data(),
            selectionNames[(int)sel].Data(),
            histoNames[p].Data()
          ), histoTitles[p], nbins[p], xmin[p], xmax[p]);
      }
      if(debug) printf("\tMade plot %s\n", histos[lep][iSel][p][ic]->GetName());
      histos[lep][iSel][p][ic]->Sumw2();
      histos[lep][iSel][p][ic]->SetDirectory(0);
    }
  }

  printf("Done building plotting histograms\n");

  ////////////////////////////////////////////////////////////////////////
  // Instantiate GeneralTree object for reading the ntuples
  GeneralTree gt;
  gt.is_hbb         = true;
  gt.is_fatjet      = true;
  gt.is_leptonic    = true;
  gt.btagWeights    = true;
  gt.useCMVA        = true;
  // Branches not in GeneralTree;
  std::map<TString, void*> extraAddresses;
  float normalizedWeight; extraAddresses["normalizedWeight"] = &normalizedWeight;
  // Map of the branches
  std::map<TString, TBranch*> b;
  // Dummy tree for a list of branch names
  TTree *dummyTree = new TTree("dummyTree","dummyTree");
  gt.WriteTree(dummyTree);
  // Done making GeneralTree object
  ////////////////////////////////////////////////////////////////////////

  printf("Building uncertainty histograms, please wait...\n");
  // Declare shape uncertainty histograms
  TH1F *histo_Baseline                           [nLepSel][selections.size()][nPlotCategories];  
  TH1F *histo_PDFUp                              [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_PDFDown                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_VHCorrUp                           [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_VHCorrDown                         [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDr1f2                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDr1f5                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDr2f1                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDr2f2                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDr5f1                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDr5f5                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDrScaleUp                        [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDrScaleDown                      [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDfScaleUp                        [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_QCDfScaleDown                      [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_eleSFUp                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_eleSFDown                          [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_muSFUp                             [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_muSFDown                           [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_btag[GeneralTree::nCsvShifts][5][3][nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_VGluUp                             [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_VGluDown                           [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_jes[(int)shiftjes::N]              [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_jesDown                            [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_doubleBUp                          [nLepSel][selections.size()][nPlotCategories];
  TH1F *histo_doubleBDown                        [nLepSel][selections.size()][nPlotCategories];
  
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(unsigned iSel=0; iSel<selections.size(); iSel++)
  for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
    histo_Baseline     [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s"              , plotBaseNames[ic].Data()));
    histo_VHCorrUp     [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_VHCorrUp"     , plotBaseNames[ic].Data()));
    histo_VHCorrDown   [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_VHCorrDown"   , plotBaseNames[ic].Data()));
    histo_PDFUp        [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_PDFUp"        , plotBaseNames[ic].Data()));
    histo_PDFDown      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_PDFDown"      , plotBaseNames[ic].Data()));
    histo_QCDr1f2      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDr1f2"      , plotBaseNames[ic].Data()));
    histo_QCDr1f5      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDr1f5"      , plotBaseNames[ic].Data()));
    histo_QCDr2f1      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDr2f1"      , plotBaseNames[ic].Data()));
    histo_QCDr2f2      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDr2f2"      , plotBaseNames[ic].Data()));
    histo_QCDr5f1      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDr5f1"      , plotBaseNames[ic].Data()));
    histo_QCDr5f5      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDr5f5"      , plotBaseNames[ic].Data()));
    histo_QCDrScaleUp  [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDrScaleUp"  , plotBaseNames[ic].Data()));
    histo_QCDrScaleDown[lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDrScaleDown", plotBaseNames[ic].Data()));
    histo_QCDfScaleUp  [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDfScaleUp"  , plotBaseNames[ic].Data()));
    histo_QCDfScaleDown[lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_QCDfScaleDown", plotBaseNames[ic].Data()));
    histo_eleSFUp      [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_eleSFUp"      , plotBaseNames[ic].Data()));
    histo_eleSFDown    [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_eleSFDown"    , plotBaseNames[ic].Data()));
    histo_muSFUp       [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_muSFUp"       , plotBaseNames[ic].Data()));
    histo_muSFDown     [lep][iSel][ic] = (TH1F*)histos[lep][iSel][0][ic]->Clone(Form("histo_%s_muSFDown"     , plotBaseNames[ic].Data()));
    for(unsigned iPt=0; iPt<5; iPt++)
    for(unsigned iEta=0; iEta<3; iEta++)
    for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
      GeneralTree::csvShift shift = gt.csvShifts[iShift];
      if (shift==GeneralTree::csvCent) continue;
      histo_btag[GeneralTree::csvJESup][iPt][iEta][lep][iSel][ic] =
        (TH1F*)histos[lep][iSel][0][ic]->Clone(
          Form("histo_%s_btag%s_pt%d_eta%d",
          plotBaseNames[ic].Data(),btagShiftName(shift),iPt,iEta)
        );
    }
  }
  printf("Done building uncertainty histograms\n");


  ////////////////////////////////////////////////////////////////////////
  // Load corrections to apply offline
  // CMVA reweighting for 2016, DeepCSV reweighting for 2017
  BTagCalibrationReader *deepcsvSFs=0; BTagCalibration *deepcsvCalib=0; CSVHelper *cmvaReweighter=0;
  std::vector<std::string> btagSystNames;
  if(year==2016) {
    cmvaReweighter = new CSVHelper(
      "PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root", 
      "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root", 
       5);
  } else if(year==2017) {
    // get the official syst name strings
    btagSystNames.reserve(GeneralTree::nCsvShifts);
    for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
      GeneralTree::csvShift shift = gt.csvShifts[iShift];
      if (shift==GeneralTree::csvCent) continue;
      btagSystNames.push_back(GeneralTree::csvShiftName(shift).Data());
    }    
    deepcsvSFs = new BTagCalibrationReader(
      BTagEntry::OP_RESHAPING,
      GeneralTree::csvShiftName(GeneralTree::csvCent).Data(), 
      btagSystNames);
    BTagCalibration *deepcsvCalib = new BTagCalibration(
      "DeepCSV", 
      "PandaAnalysis/data/csv/DeepCSV_94XSF_V2_B_F.csv");
    deepcsvSFs->load(*(deepcsvCalib), BTagEntry::FLAV_B, "iterativeFit");
    deepcsvSFs->load(*(deepcsvCalib), BTagEntry::FLAV_C, "iterativeFit");
    deepcsvSFs->load(*(deepcsvCalib), BTagEntry::FLAV_UDSG, "iterativeFit");
    
  }
  vector<double> ZjetsEWKCorr = EWKCorrPars(kZjets);

  // CMVA jet kinematic decorrelated weight nuisances
  vector<double> jetPts[5][3], jetPtsUp[5][3], jetPtsDown[5][3], jetEtas[5][3], jetBtags[5][3];
  vector<int> jetFlavors[5][3];
  for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
    jetPts    [iPt][iEta].reserve(20);
    jetPtsUp  [iPt][iEta].reserve(20);
    jetPtsDown[iPt][iEta].reserve(20);
    jetEtas   [iPt][iEta].reserve(20);
    jetBtags  [iPt][iEta].reserve(20);
    jetFlavors[iPt][iEta].reserve(20);
  }
  // Done loading offline corrections
  ////////////////////////////////////////////////////////////////////////
  // Declare local analysis variables
  // Isojets: For boosted categories, the number of 30 GeV AK4 jets not in the fat jet
  vector<int> nIsojet(NJES);
  vector<vector<unsigned char>> isojets(NJES);
  vector<unsigned char> isojetNBtags(NJES);
  // Selection bits
  unsigned selectionBits[NJES], nMinusOneBits;
  unsigned char typeLepSel;
  unsigned char category; int countB;
  // Declare variables for the variations of the weights
  float weight;
  float weight_QCDr1f2, weight_QCDr1f5, weight_QCDr2f1, weight_QCDr2f2, weight_QCDr5f1, weight_QCDr5f5;
  float weight_lepSFUp;
  float weight_btag[GeneralTree::nCsvShifts][5][3];
  float weight_VHCorrUp, weight_VHCorrDown; 
  // Done declaring weight variables
  ////////////////////////////////////////////////////////////////////////

  // Analysis Cuts
  float isojetBtagCut = (year==2017)? deepcsvLoose : cmvaLoose;
  std::map<selectionType,vector<TString>> cuts;
  cuts[kZllHLightFlavorCR  ] = {"ZpT","pTjj","bveto","Zmass"                                };
  cuts[kZllHHeavyFlavorCR  ] = {"ZpT","pTjj","btag" ,"ZmassTight","lowMET","dPhiZH","mjjSB" };
  cuts[kZllH2TopCR         ] = {"ZpT","pTjj","btag" ,"ZmassSB"                              };
  cuts[kZllHSR             ] = {"ZpT","pTjj","btag" ,"Zmass"              ,"dPhiZH","mjj"   };
  cuts[kZllHPresel         ] = {"ZpT","pTjj"                                                };
  cuts[kZllHLightFlavorFJCR] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "0ijb", "Zmass"               };
  cuts[kZllHHeavyFlavorFJCR] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD_SB", "0ijb", "ZmassTight", "lowMET"};
  cuts[kZllHTT1bFJCR       ] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "1ijb", "Zmass"               };
  cuts[kZllHTT2bFJCR       ] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "1ijb", "Zmass"               };
  cuts[kZllHFJSR           ] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD_SR", "0ijb", "Zmass"               };
  cuts[kZllHFJPresel       ] = {"ZpTFJ","pTFJ","dPhiZHFJ"                                         };

  ////////////////////////////////////////////////////////////////////////
  // Begin Chain Loop
  TTree *events=0; 
  for(auto const &sample: samples) {
    TString sampleName = sample.first; sampleType type = sample.second;
    TString inputFileName = Form("%s%s.root",ntupleDir.Data(),sampleName.Data());
    printf("### Opening sample %s ###\n",sampleName.Data());
    TFile *inputFile = TFile::Open(inputFileName,"READ");
    if(!inputFile) { throw std::runtime_error("Problem opening file"); }
    events = (TTree*)inputFile->Get("events");
    if(!events) { throw std::runtime_error("Problem loading tree"); }
    
    // Sample properties
    bool isNLOZjets = sampleName.Contains("ZJets_pt"); // Only use events with HT<100 for NLO pt binned samples in 2016
    // End sample properties
    
    
    // Nasty hack code to get the tree branches and their addresses, have to do it for each file
    TObjArray *listOfBranches = events->GetListOfBranches();
    for(unsigned iB=0; iB<(unsigned)listOfBranches->GetEntries(); iB++) {
      TBranch *branch = (TBranch*)listOfBranches->At(iB);
      TString branchName = branch->GetName();
      if(branchName.Contains("jotBReg_")) continue;
      bool isExtraBranch=false;
      auto x = extraAddresses.find(branchName);
      isExtraBranch= (x!=extraAddresses.end());
      if(isExtraBranch) {
        b[x->first] = events->GetBranch(x->first);
        if(!b[x->first]) { throw std::runtime_error(Form("Extra branch %s could not be found in events tree\n", x->first.Data())); }
        b[x->first]->SetAddress(x->second);
        if(debug>=2) printf("\tBooking extra branch \"%s\"\n", x->first.Data());
        continue;
      }
      TBranch *dummyBranch = dummyTree->GetBranch(branchName);
      if(!dummyBranch) { throw std::runtime_error(Form("WARNING: Couldn't find branch \"%s\" in the dummyTree, skipping", branchName.Data())); continue; }
      void* address=dummyBranch->GetAddress();
      b[branchName] = events->GetBranch(branchName);
      if(!b[branchName]) { throw std::runtime_error(Form("Branch \"%s\" could not be found in events tree", x->first.Data())); }
      b[branchName]->SetAddress(address);
      if(debug>=2) printf("\tBooking GeneralTree branch \"%s\"\n", branchName.Data());
    }
    // End Book Branches

    Long64_t nentries = events->GetEntries();
    //nentries = TMath::Min(Long64_t(1e5),nentries);
    // Begin Event Loop
    for (Long64_t ientry=0; ientry<nentries; ientry++) {
      if(debug && ientry!=0) usleep(2e5);
      if(debug || ientry%100000==0) printf("> Reading entry %lld/%lld ...\n",ientry,nentries);
      
      // Stitching for low pT NLO Z+jets in 2016
      if(isNLOZjets && year==2016) {
        bLoad(b["lheHT"],ientry);
        if(gt.lheHT>=100) continue;
      }

      //////////////////////
      // Clear variables
      typeLepSel=99; // 0: Z(mm), 1: Z(ee), 2: e-mu

      // Clear the jet pts for the b tag decorrelation
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        jetPts    [iPt][iEta].clear();
        jetPtsUp  [iPt][iEta].clear();
        jetPtsDown[iPt][iEta].clear();
        jetEtas   [iPt][iEta].clear();
        jetBtags  [iPt][iEta].clear();
        jetFlavors[iPt][iEta].clear();
      }
      for(unsigned iJES=0; iJES<NJES; iJES++) {
        nIsojet[iJES]=0;
        isojetNBtags[iJES]=0;
        isojets[iJES].clear();
      }
      // Done clearing variables
      //////////////////////
      
      // Basic cuts
      bLoad(b["nLooseLep"],ientry);
      if(gt.nLooseLep!=2) continue; 
      if(debug) printf("  Passed lepton multiplicity\n");

      // Trigger
      if(type==kData) {
        bLoad(b["trigger"],ientry);
        bool passTrigger = (gt.trigger & whichTriggers) !=0;
        if(!passTrigger) continue;
        if(debug) printf("  Passed trigger\n");
      }

      // Lepton ID and isolation
      bLoad(b["nLooseElectron"],ientry);
      bLoad(b["nTightElectron"],ientry);
      bLoad(b["nLooseMuon"],ientry);
      bLoad(b["nTightMuon"],ientry);
      bLoad(b["muonSelBit"],ientry);
      bLoad(b["muonPt"],ientry);
      bLoad(b["muonPdgId"],ientry);
      bLoad(b["electronSelBit"],ientry);
      bLoad(b["electronPt"],ientry);
      bLoad(b["electronPdgId"],ientry);

      if(gt.nLooseMuon==2 && gt.nLooseElectron==0 &&
        gt.muonPt[0]>25 && gt.muonPt[1]>10 && 
        gt.muonPdgId[0]+gt.muonPdgId[1]==0 &&
        (gt.muonSelBit[0] & 1<<3)!=0 &&
        (gt.muonSelBit[1] & 1<<3)!=0
      ) typeLepSel=0;
      else if(gt.nLooseElectron==2 && gt.nLooseMuon==0 &&
        gt.electronPt[0]>25 && gt.electronPt[1]>15 && 
        gt.electronPdgId[0]+gt.electronPdgId[1]==0 &&
        (gt.electronSelBit[0] & 1<<5)!=0 &&
        (gt.electronSelBit[1] & 1<<5)!=0
      ) typeLepSel=1;
      else if(gt.nLooseMuon==1 && gt.nLooseElectron==1 &&
        ((gt.muonPt[0]>25 && gt.electronPt[0]>15) || (gt.muonPt[0]>15 && gt.electronPt[0]>25)) && 
        gt.muonPdgId[0]*gt.electronPdgId[0]<0 &&
        (gt.muonSelBit[0]     & 1<<3)!=0 &&
        (gt.electronSelBit[0] & 1<<5)!=0
      ) typeLepSel=2;
      if(typeLepSel!=0 && typeLepSel!=1 && typeLepSel!=2) continue;
      if(debug) printf("  Passed lepton ID/iso multiplicity\n");

      // Lepton kinematics
      float lepton1Pt,lepton1Eta,lepton1Phi,lepton1RelIso,lepton1D0,lepton1DZ,
            lepton2Pt,lepton2Eta,lepton2Phi,lepton2RelIso,lepton2D0,lepton2DZ;
      if (typeLepSel==0) {
        bLoad(b["muonEta"],ientry);
        bLoad(b["muonPhi"],ientry);
        bLoad(b["muonD0"],ientry);
        bLoad(b["muonDZ"],ientry);
        bLoad(b["muonCombIso"],ientry);
        bLoad(b["muonPdgId"],ientry);
        lepton1Pt     = gt.muonPt[0];
        lepton1Eta    = gt.muonEta[0];
        lepton1Phi    = gt.muonPhi[0]; 
        lepton1D0     = gt.muonD0[0];
        lepton1DZ     = gt.muonDZ[0];
        lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
        lepton2Pt     = gt.muonPt[1];
        lepton2Eta    = gt.muonEta[1];
        lepton2Phi    = gt.muonPhi[1]; 
        lepton2D0     = gt.muonD0[1];
        lepton2DZ     = gt.muonDZ[1];
        lepton2RelIso = gt.muonCombIso[1]/gt.muonPt[1];
      } else if(typeLepSel==1) {
        bLoad(b["electronEta"],ientry);
        bLoad(b["electronPhi"],ientry);
        bLoad(b["electronD0"],ientry);
        bLoad(b["electronDZ"],ientry);
        bLoad(b["electronCombIso"],ientry);
        bLoad(b["electronPdgId"],ientry);
        lepton1Pt     = gt.electronPt[0]; 
        lepton1Eta    = gt.electronEta[0];
        lepton1Phi    = gt.electronPhi[0];
        lepton1RelIso = gt.electronCombIso[0]/gt.electronPt[0];
        lepton1D0     = gt.electronD0[0];
        lepton1DZ     = gt.electronDZ[0];
        lepton2Pt     = gt.electronPt[1]; 
        lepton2Eta    = gt.electronEta[1];
        lepton2Phi    = gt.electronPhi[1];
        lepton2RelIso = gt.electronCombIso[1]/gt.electronPt[1];
        lepton2D0     = gt.electronD0[1];
        lepton2DZ     = gt.electronDZ[1];
      } else if (typeLepSel==2) {
        bLoad(b["muonEta"],ientry);
        bLoad(b["muonPhi"],ientry);
        bLoad(b["muonD0"],ientry);
        bLoad(b["muonDZ"],ientry);
        bLoad(b["muonCombIso"],ientry);
        bLoad(b["muonPdgId"],ientry);
        bLoad(b["electronEta"],ientry);
        bLoad(b["electronPhi"],ientry);
        bLoad(b["electronD0"],ientry);
        bLoad(b["electronDZ"],ientry);
        bLoad(b["electronCombIso"],ientry);
        bLoad(b["electronPdgId"],ientry);
        lepton1Pt     = gt.muonPt[0];
        lepton1Eta    = gt.muonEta[0];
        lepton1Phi    = gt.muonPhi[0]; 
        lepton1D0     = gt.muonD0[0];
        lepton1DZ     = gt.muonDZ[0];
        lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
        lepton2Pt     = gt.electronPt[0]; 
        lepton2Eta    = gt.electronEta[0];
        lepton2Phi    = gt.electronPhi[0];
        lepton2D0     = gt.electronD0[0];
        lepton2DZ     = gt.electronDZ[0];
        lepton2RelIso = gt.electronCombIso[0]/gt.electronPt[0];
      }
      if(debug) printf("  Passed lepton kinematics\n");

      bLoad(b["ZBosonPt"],ientry);
      bLoad(b["ZBosonEta"],ientry);
      bLoad(b["ZBosonPhi"],ientry);
      bLoad(b["ZBosonM"],ientry);
      if(gt.ZBosonPt<30) continue;
      if(gt.ZBosonM<10) continue;
      if(debug) printf("  Passed Z boson reconstruction\n");
      
      // Jet multiplicity
      bool isBoostedCategory=false;
      if(useBoostedCategory) { 
        bLoad(b["nFatjet"],ientry);
        bLoad(b["fjMSD_corr"],ientry);
        bLoad(b["fjPt"],ientry);
        bLoad(b["fjEta"],ientry);
        bLoad(b["fjPhi"],ientry);
        float temp_deltaPhiVH;
        if(gt.nFatjet>0) temp_deltaPhiVH = fabs(TVector2::Phi_mpi_pi(gt.ZBosonPhi-gt.fjPhi));
        if(
          gt.nFatjet != 0 && 
          gt.fjPt[0] >= 250 && 
          gt.fjMSD_corr[0] >= 40 &&
          fabs(gt.fjEta) < 2.4 &&
          gt.ZBosonPt >= 250 &&
          temp_deltaPhiVH >= 2.5
        ) isBoostedCategory=true;
      }
      // Category Assignment for Plotting and Datacards
      if(type!=kData) {
        if(isBoostedCategory) {
          bLoad(b["fjGenNumB"],ientry);
          countB = gt.fjGenNumB;
        } else {
          bLoad(b["nBGenJets"],ientry);
          countB = gt.nBGenJets;
        }
      }
      if     (type==vhbbPlot::kData ) category=kPlotData ;
      else if(type==vhbbPlot::kQCD  ) category=kPlotQCD  ;
      else if(type==vhbbPlot::kWW   ) category=kPlotVVLF ;
      else if(type==vhbbPlot::kTT   ) category=kPlotTT   ; 
      else if(type==vhbbPlot::kTop  ) category=kPlotTop  ;
      else if(type==vhbbPlot::kWH   ) category=kPlotWH   ;
      else if(type==vhbbPlot::kZH   ) category=kPlotZH   ;
      else if(type==vhbbPlot::kWjets) {
        if(countB>1) category=kPlotWbb;
        else if(countB>0) category=kPlotWb;
        else category=kPlotWLF;
      } else if(type==vhbbPlot::kZjets) {
        if(countB>1) category=kPlotZbb;
        else if(countB>0) category=kPlotZb;
        else category=kPlotZLF;
      } else if(type==kVZ) {
        if(countB>1) category=kPlotVZbb;
        else category=kPlotVVLF;
      } else throw std::runtime_error("category problem!");

      // Jet kinematics
      if(isBoostedCategory) {
        // No checks here? 
      } else { 
        bLoad(b["nJet"],ientry);
        bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        bLoad(b["hbbpt"],ientry);
        bLoad(b["hbbm_reg"],ientry);
        if(
          gt.nJet[0]<2 || 
          gt.hbbpt[0]<50 || 
          gt.hbbm_reg[0]<0 ||
          gt.hbbm_reg[0]>250) 
          continue;
      }
      if(debug) printf("  Passed jet kinematics\n");
      if(debug) printf("Passed preselection!\n");
      
      // Met
      bLoad(b["pfmet"],ientry);
      //bLoad(b["pfmetphi"],ientry);
      
      // Jet B-tagging
      bool bjet1IsTight, bjet2IsLoose;
      float bjet1btag, bjet2btag;
      if(year==2016) {
        bLoad(b["jotCMVA"],ientry);
        bjet1IsTight = gt.jotCMVA[gt.hbbjtidx[0][0]] > cmvaTight;
        bjet2IsLoose = gt.jotCMVA[gt.hbbjtidx[0][1]] > cmvaLoose;
        bjet1btag = gt.jotCMVA[gt.hbbjtidx[0][0]];
        bjet2btag = gt.jotCMVA[gt.hbbjtidx[0][1]];
      } else if(year==2017) {
        bLoad(b["jotCSV"],ientry);
        bjet1IsTight = gt.jotCSV[gt.hbbjtidx[0][0]] > deepcsvTight;
        bjet2IsLoose = gt.jotCSV[gt.hbbjtidx[0][1]] > deepcsvLoose;
        bjet1btag = gt.jotCMVA[gt.hbbjtidx[0][0]];
        bjet2btag = gt.jotCMVA[gt.hbbjtidx[0][1]];
      }

      bLoad(b["jotPt"],ientry);
      bLoad(b["jotEta"],ientry);
      bLoad(b["jotPhi"],ientry);
      float bjet1Pt = gt.jotPt[0][gt.hbbjtidx[0][0]];
      float bjet2Pt = gt.jotPt[0][gt.hbbjtidx[0][1]];
      
      if(category!=kPlotData && category!=kPlotQCD) { if(isBoostedCategory) {
        // Isojets for boosted category
        bLoad(b["fjEta"],ientry);
        bLoad(b["fjPhi"],ientry);
        for(unsigned iJES=0; iJES<NJES; iJES++) 
        for(unsigned char iJ=0; iJ<gt.nJot[iJES]; iJ++) {
          if(fabs(gt.jotEta[iJ])>2.4) continue;
          float dR2JetFatjet=pow(gt.jotEta[iJ]-gt.fjEta,2)+pow(TVector2::Phi_mpi_pi(gt.jotPhi[iJ]-gt.fjPhi),2);
          if(dR2JetFatjet<0.64) continue;
          
          float isojetBtag = (year==2017)? gt.jotCMVA[iJ] : gt.jotCSV[iJ];
          if(gt.jotPt[iJES][iJ]>30) {
            isojets[iJES].push_back(iJ);
            if(isojetBtag>isojetBtagCut) isojetNBtags[iJES]++;
          }
          
          // Decorrelated b-tag nuisances for the isojets
          // Nominal JES scenario only
          if(iJES!=0) continue;
          int iPt=-1, iEta=-1;
          double jetAbsEta=fabs(gt.jotEta[iJ]);
          //if      (gt.jetPt[iJ] >= 19.99 && gt.jetPt[iJ] < 30 ) iPt = 0;
          //else if (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
          if      (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
          else if (gt.jetPt[iJ] >= 40    && gt.jetPt[iJ] < 60 ) iPt = 2;
          else if (gt.jetPt[iJ] >= 60    && gt.jetPt[iJ] < 100) iPt = 3;
          else if (gt.jetPt[iJ] >= 100                        ) iPt = 4;
          if      (jetAbsEta >= 0   && jetAbsEta < 0.8  ) iEta = 0;
          else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6  ) iEta = 1;
          else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41 ) iEta = 2;
          
          if(iPt>=0 && iEta>=0) {
            jetPts    [iPt][iEta].push_back(gt.jotPt[iJES][iJ]);
            jetEtas   [iPt][iEta].push_back(gt.jotEta[iJ]);
            jetFlavors[iPt][iEta].push_back(gt.jotFlav[iJ]);
            jetBtags  [iPt][iEta].push_back(isojetBtag);
          }
        }
      } else {
        // In resolved case, need decorrelated b-tag nuisances for all the ak4 jets
        // Nominal JES scenario only
        for(unsigned char iJ=0; iJ<gt.nJot[0]; iJ++) {
          int iPt=-1, iEta=-1;
          double jetAbsEta=fabs(gt.jotEta[iJ]);
          if      (gt.jetPt[iJ] >= 19.99 && gt.jetPt[iJ] < 30 ) iPt = 0;
          else if (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
          if      (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
          else if (gt.jetPt[iJ] >= 40    && gt.jetPt[iJ] < 60 ) iPt = 2;
          else if (gt.jetPt[iJ] >= 60    && gt.jetPt[iJ] < 100) iPt = 3;
          else if (gt.jetPt[iJ] >= 100                        ) iPt = 4;
          if      (jetAbsEta >= 0   && jetAbsEta < 0.8  ) iEta = 0;
          else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6  ) iEta = 1;
          else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41 ) iEta = 2;
          if(iPt>=0 && iEta>=0) {
            float btag = (year==2017)? gt.jotCMVA[iJ] : gt.jotCSV[iJ];
            jetPts    [iPt][iEta].push_back(gt.jotPt[0][iJ]);
            jetEtas   [iPt][iEta].push_back(gt.jotEta[iJ]);
            jetFlavors[iPt][iEta].push_back(gt.jotFlav[iJ]);
            jetBtags  [iPt][iEta].push_back(btag);
          }
        }
      }}
      
      // Load branches and calculate stuff for the cuts
      bLoad(b["ZBosonPt"],ientry);
      bLoad(b["ZBosonPhi"],ientry);
      bLoad(b["hbbpt_reg"],ientry);
      bLoad(b["hbbphi"],ientry);
      bLoad(b["hbbm_reg"],ientry);
      bLoad(b["fjPt"],ientry);
      bLoad(b["fjPhi"],ientry);
      bLoad(b["fjMSD_corr"],ientry);
      float deltaPhiZH    = !isBoostedCategory? fabs(TVector2::Phi_mpi_pi(gt.hbbphi[0] - gt.ZBosonPhi)) : -1;
      float deltaPhiZHFJ  = isBoostedCategory? fabs(TVector2::Phi_mpi_pi(gt.fjPhi - gt.ZBosonPhi)) : -1;
      float ptBalanceZH   = !isBoostedCategory? gt.hbbpt_reg[0] /  gt.ZBosonPt : -1;
      float ptBalanceZHFJ = isBoostedCategory? gt.fjPt[0] / gt.ZBosonPt : -1;

      // Set Selection Bits
      std::map<TString, bool> cut;
      for(unsigned iJES=0; iJES<NJES; iJES++) {
        if(iJES==0) {
          cut["ZpT"       ] = gt.ZBosonPt > 50;
          cut["ZpTFJ"     ] = gt.ZBosonPt > 250;
          cut["dPhiZHFJ"  ] = fabs(gt.fjPhi - gt.ZBosonPhi) > 2.5;
          cut["btag"      ] = bjet1IsTight && bjet2IsLoose;
          cut["bveto"     ] = !bjet1IsTight && !bjet2IsLoose;
          cut["btagFJ"    ] = gt.fjDoubleCSV > doubleBCut;
          cut["bvetoFJ"   ] = gt.fjDoubleCSV < doubleBCut;
          cut["Zmass"     ] = gt.ZBosonM >= 75 && gt.ZBosonM < 105;
          cut["ZmassTight"] = gt.ZBosonM >= 85 && gt.ZBosonM < 97;
          cut["ZmassSB"   ] = (gt.ZBosonM >= 10 && gt.ZBosonM < 75) || (gt.ZBosonM >= 105 && gt.ZBosonM < 120);
          cut["lowMET"    ] = gt.pfmet[0] < 60;
        } 
        if(iJES==(int)shiftjes::kJESTotalUp  ) cut["lowMET"] = gt.pfmet[1] < 60;
        if(iJES==(int)shiftjes::kJESTotalDown) cut["lowMET"] = gt.pfmet[2] < 60;
        if(isBoostedCategory) {
          cut["mSD"     ] = gt.fjMSD_corr[iJES] >= 40;
          cut["mSD_SR"  ] = gt.fjMSD_corr[iJES] >= 80 && gt.fjMSD_corr[iJES]<150;
          cut["mSD_SB"  ] = cut["mSD"] && gt.fjMSD_corr[iJES]<80;
          cut["pTFJ"    ] = gt.fjPt[iJES] > 250;
          cut["0ijb"    ] = isojetNBtags[iJES]==0;
          cut["1ijb"    ] = !cut["0ijb"];
        } else {
          cut["dPhiZH"  ] = fabs(gt.hbbphi[iJES] - gt.ZBosonPhi) > 2.5;
          cut["pTjj"    ] = gt.hbbpt_reg[iJES] > 100;
          cut["mjj"     ] = gt.hbbm_reg[iJES] >= 90 && gt.hbbm_reg[iJES] < 150;
          cut["mjjSB"   ] = !cut["mjj"] && gt.hbbm_reg[iJES]<250;
        } 
        selectionBits[iJES]=0; nMinusOneBits=0;
        for(unsigned iSel=0; iSel<selections.size(); iSel++) { 
          selectionType sel = selections[iSel];
          if(passAllCuts(cut, cuts[sel])) {
            selectionBits[iJES] |= sel;
            if(debug>=1) printf("  Passed cuts for region %s\n",vhbbPlot::selectionNames[sel].Data());
          }
          if(iJES==0 && passNMinusOne(cut, cuts[sel]))
            nMinusOneBits |= sel;
        }
      
      }

      // Begin Event weighting
      if(type==kData) {
        weight=1;
      } else {
        bLoad(b["normalizedWeight"],ientry);
        bLoad(b["pu"],ientry);
        float puWeight = nPUScaleFactor(puWeights, gt.pu);
        weight = normalizedWeight * theLumi * puWeight; 
        
        if(type==kWjets || type==kZjets) {
          if(useHtBinnedVJetsKFactor) {
            bLoad(b["trueGenBosonPt"],ientry);
            if(isNLOZjets) {
              gt.sf_qcdV=1;
            } else {
              bLoad(b["lheHT"],ientry);
              gt.sf_qcdV=qcdKFactor(type, gt.lheHT);
            }
            gt.sf_ewkV=ZjetsEWKCorr[0]+ZjetsEWKCorr[1]*
              (TMath::Power((gt.trueGenBosonPt+ZjetsEWKCorr[2]),ZjetsEWKCorr[3]));
          } else {
            bLoad(b["sf_qcdV"],ientry);
            bLoad(b["sf_ewkV"],ientry);
          }
          weight *= gt.sf_qcdV * gt.sf_ewkV;
        }
        
        if (typeLepSel==0) {
          bLoad(b["muonSfReco"],ientry);
          bLoad(b["muonSfTight"],ientry);
          bLoad(b["muonSfUnc"],ientry);
          bLoad(b["sf_muTrig"],ientry);
          weight *= 1;//gt.sf_muTrig;
          weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
          weight *= gt.muonSfReco[1] * gt.muonSfTight[1];
        } else if(typeLepSel==1) {
          bLoad(b["electronSfReco"],ientry);
          bLoad(b["electronSfMvaWP90"],ientry);
          bLoad(b["electronSfUnc"],ientry);
          bLoad(b["sf_eleTrig"],ientry);
          weight *= 1;//gt.sf_eleTrig;
          weight *= gt.electronSfReco[0] * gt.electronSfMvaWP90[0];
          weight *= gt.electronSfReco[1] * gt.electronSfMvaWP90[1];
        } else if (typeLepSel==2) {
          bLoad(b["muonSfReco"],ientry);
          bLoad(b["muonSfTight"],ientry);
          bLoad(b["muonSfUnc"],ientry);
          bLoad(b["electronSfReco"],ientry);
          bLoad(b["electronSfMvaWP90"],ientry);
          bLoad(b["electronSfUnc"],ientry);
          bLoad(b["sf_eleTrig"],ientry);
          weight *= 1.0; // no trigger data/MC SF
          weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
          weight *= gt.electronSfReco[0] * gt.electronSfMvaWP90[0];
        }

        float recorrect_vhEWK=1, recorrect_vhEWKUp=1, recorrect_vhEWKDown=1;
        if(type==kVZ) {
          bLoad(b["sf_wz"],ientry);
          bLoad(b["sf_zz"],ientry);
          weight *= gt.sf_wz * gt.sf_zz;
        } else if(type==kZH) {
          bLoad(b["sf_vh"],ientry);
          bLoad(b["sf_vhUp"],ientry);
          bLoad(b["sf_vhDown"],ientry);
          //TString vhChannel=lepton1Charge>0?"WplusH":"WminusH";
          recorrect_vhEWK     = vhEWKCorr("ZllH",gt.sf_vh    );
          recorrect_vhEWKUp   = vhEWKCorr("ZllH",gt.sf_vhUp  );
          recorrect_vhEWKDown = vhEWKCorr("ZllH",gt.sf_vhDown);
          weight *= recorrect_vhEWK;
        }
        //bLoad(b["sf_cmvaWeight_Cent"],ientry);
        //weight *= gt.sf_csvWeights[GeneralTree::csvCent];
      }

      // Hack for the central weights 
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        if(year==2016) {
          // in 2016 we can use the CSVHelper to calculate the total shift for each jet kinematic bin
          double cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
          weight *= cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetBtags[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCent, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
        } else if(year==2017) {
          // in 2017, we have to calculate the reshape factor for each jet in each kinematic bin
          unsigned iShift=GeneralTree::csvCent;
          GeneralTree::csvShift theShift = gt.csvShifts[iShift];
          for(unsigned iJ=0; iJ<jetPts[iPt][iEta].size(); iJ++) {
            unsigned absid = abs(jetFlavors[iPt][iEta][iJ]);
            BTagEntry::JetFlavor flav = absid == 5 ? BTagEntry::FLAV_B : 
              (absid == 4 ? BTagEntry::FLAV_C : BTagEntry::FLAV_UDSG);
            float reshapeFactor = deepcsvSFs->eval_auto_bounds(
              GeneralTree::csvShiftName(theShift).Data(),
              flav,
              jetEtas[iPt][iEta][iJ], 
              jetPts[iPt][iEta][iJ],
              jetBtags[iPt][iEta][iJ]
            );
            if(reshapeFactor>0.001) weight *= reshapeFactor;
          }
        }
      }
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        if(year==2016) {
          // in 2016 we can use the CSVHelper to calculate the total shift for each jet kinematic bin
          double cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
          double centralWeight = cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetBtags[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCent, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
          for(unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
            GeneralTree::csvShift theShift = gt.csvShifts[iShift];
            weight_btag[iShift][iPt][iEta] = weight*cmvaReweighter->getCSVWeight(
              jetPts[iPt][iEta], jetEtas[iPt][iEta], jetBtags[iPt][iEta], jetFlavors[iPt][iEta],
              theShift,
              cmvaWgtHF, cmvaWgtLF, cmvaWgtCF
            )/centralWeight; 
          }
        } else if(year==2017) {
          // in 2017, we have to calculate the reshape factor for each jet in each kinematic bin
          for(unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
            GeneralTree::csvShift theShift = gt.csvShifts[iShift];
            weight_btag[iShift][iPt][iEta] = weight;
            for(unsigned iJ=0; iJ<jetPts[iPt][iEta].size(); iJ++) {
              unsigned absid = abs(jetFlavors[iPt][iEta][iJ]);
              BTagEntry::JetFlavor flav = absid == 5 ? BTagEntry::FLAV_B : 
                (absid == 4 ? BTagEntry::FLAV_C : BTagEntry::FLAV_UDSG);
              float reshapeFactor = deepcsvSFs->eval_auto_bounds(
                GeneralTree::csvShiftName(theShift).Data(),
                flav,
                jetEtas[iPt][iEta][iJ], 
                jetPts[iPt][iEta][iJ],
                jetBtags[iPt][iEta][iJ]
              );
              if(reshapeFactor>0.001) weight_btag[iShift][iPt][iEta] *= reshapeFactor;
            }
          }
        }
      }

      // Calculate the shape variable in all JES scenarios for all the regions
      float MVAVar[selections.size()][(int)shiftjes::N], bdtValue[(int)shiftjes::N];
      for(unsigned iSel=0; iSel<selections.size(); iSel++) for(unsigned iJES=0; iJES<NJES; iJES++) {
        selectionType sel = selections[iSel];
        switch(MVAVarType) {
          case 1:
          default:
            if(sel==kZllHLightFlavorCR ||
              sel==kZllHHeavyFlavorCR ||
              sel==kZllH2TopCR        ||
              sel==kZllHSR            ||
              sel==kZllHPresel        )
              MVAVar[iSel][iJES]=gt.hbbpt_reg[iJES];
            else
              MVAVar[iSel][iJES]=gt.fjPt[iJES];
            break;
          case 3:
            if(sel==kZllHSR || sel==kZllHFJSR)
              MVAVar[iSel][iJES]=bdtValue[iJES];
            else if(sel==kZllHLightFlavorCR ||
              sel==kZllHHeavyFlavorCR       || 
              sel==kZllH2TopCR)
              MVAVar[iSel][iJES]=bjet2btag;
            else if(sel==kZllHLightFlavorFJCR ||
              sel==kZllHHeavyFlavorFJCR       ||
              sel==kZllHTT2bFJCR              ||
              sel==kZllHTT1bFJCR)
              MVAVar[iSel][iJES]=gt.fjMSD_corr[iJES];
            break;
        }
      }

      // Fill the plotting histograms
      bLoad(b["ZBosonLep1CosThetaCS"],ientry);
      bLoad(b["ZBosonLep1CosThetaStar"],ientry);
      bLoad(b["ZBosonLep1CosThetaStarFJ"],ientry);
      bLoad(b["sumEtSoft1"],ientry);
      bLoad(b["nSoft2"],ientry);
      bLoad(b["nSoft5"],ientry);
      bLoad(b["nSoft10"],ientry);
      for(unsigned iSel=0; iSel<selections.size(); iSel++) {
        if(debug>=3) printf("\tTrying to fill for iSel=%d\n", iSel);
        selectionType sel = selections[iSel];
        bool passFullSel = (selectionBits[0] & sel) != 0;
        if(passFullSel && debug>=3) printf("\tPassed this sel\n");
        float theVar;
        for(int p=0; p<nPlots; p++) { 
          bool makePlot=false;
          // Variables -- change the makePlot for n-1 later
          if      (histoNames[p]=="MVAVar"                  ) { theVar = MVAVar[iSel][0]            ; makePlot = passFullSel; }
          else if (histoNames[p]=="lepton1Pt"               ) { theVar = lepton1Pt                  ; makePlot = passFullSel; }
          else if (histoNames[p]=="lepton2Pt"               ) { theVar = lepton2Pt                  ; makePlot = passFullSel; }
          else if (histoNames[p]=="lepton1Eta"              ) { theVar = lepton1Eta                 ; makePlot = passFullSel; }
          else if (histoNames[p]=="lepton2Eta"              ) { theVar = lepton2Eta                 ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonPt"                ) { theVar = gt.ZBosonPt                ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonEta"               ) { theVar = gt.ZBosonEta               ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonPhi"               ) { theVar = gt.ZBosonPhi               ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonM"                 ) { theVar = gt.ZBosonM                 ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonLep1CosThetaCS"    ) { theVar = gt.ZBosonLep1CosThetaCS    ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonLep1CosThetaStar"  ) { theVar = gt.ZBosonLep1CosThetaStar  ; makePlot = passFullSel; }
          else if (histoNames[p]=="ZBosonLep1CosThetaStarFJ") { theVar = gt.ZBosonLep1CosThetaStarFJ; makePlot = passFullSel; }
          else if (histoNames[p]=="Mjj"                     ) { theVar = gt.hbbm_reg[0]             ; makePlot = passFullSel; }
          else if (histoNames[p]=="pTjj"                    ) { theVar = gt.hbbpt_reg[0]            ; makePlot = passFullSel; }
          else if (histoNames[p]=="bjet1Pt"                 ) { theVar = bjet1Pt                    ; makePlot = passFullSel; }
          else if (histoNames[p]=="bjet2Pt"                 ) { theVar = bjet2Pt                    ; makePlot = passFullSel; }
          else if (histoNames[p]=="bjet1btag"               ) { theVar = bjet1btag                  ; makePlot = passFullSel; }
          else if (histoNames[p]=="bjet2btag"               ) { theVar = bjet2btag                  ; makePlot = passFullSel; }
          else if (histoNames[p]=="nJet"                    ) { theVar = gt.nJet[0]                 ; makePlot = passFullSel; }
          else if (histoNames[p]=="deltaPhiZH"              ) { theVar = deltaPhiZH                 ; makePlot = passFullSel; }
          else if (histoNames[p]=="deltaPhiZHFJ"            ) { theVar = deltaPhiZHFJ               ; makePlot = passFullSel; }
          else if (histoNames[p]=="ptBalanceZH"             ) { theVar = ptBalanceZH                ; makePlot = passFullSel; }
          else if (histoNames[p]=="ptBalanceZHFJ"           ) { theVar = ptBalanceZHFJ              ; makePlot = passFullSel; }
          else if (histoNames[p]=="sumEtSoft1"              ) { theVar = gt.sumEtSoft1              ; makePlot = passFullSel; }
          else if (histoNames[p]=="nSoft2"                  ) { theVar = gt.nSoft2                  ; makePlot = passFullSel; }
          else if (histoNames[p]=="nSoft5"                  ) { theVar = gt.nSoft5                  ; makePlot = passFullSel; }
          else if (histoNames[p]=="nSoft10"                 ) { theVar = gt.nSoft10                 ; makePlot = passFullSel; }
          else if (histoNames[p]=="bdtValue"                ) { theVar = bdtValue[0]                ; makePlot = passFullSel; }
          if(!makePlot) continue;
          if(debug>=3) printf("\t\tFilling %s\n",histoNames[p].Data());
          histos[typeLepSel][iSel][p][category]->Fill(theVar, weight);
        }
      }    
    
    } // End Event Loop
    inputFile->Close();
  } // End Chain Loop
  delete dummyTree; // no longer needed
  if(deepcsvSFs) delete deepcsvSFs;
  if(deepcsvCalib) delete deepcsvCalib;
  if(cmvaReweighter) delete cmvaReweighter;
  
  // Write plots
  char regionName[128];
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(unsigned iSel=0; iSel<selections.size(); iSel++) {
    selectionType sel = selections[iSel];
    sprintf(regionName, "Z%sH%s",leptonStrings[lep].Data(),selectionNames[sel].Data());
    TString plotFileName = Form("MitVHBBAnalysis/datacards/%s/plots_%s.root",dataCardDir.Data(),regionName);
    TFile *outputPlots = new TFile(plotFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
    for(int p=0; p<nPlots; p++) {
      if(histoNames[p]=="") continue;
      TDirectory *plotDir = outputPlots->mkdir(histoNames[p]); plotDir->cd();
      for(unsigned ic=kPlotData; ic!=nPlotCategories; ic++) 
        histos[lep][iSel][p][ic]->Write(Form("histo%d",ic));
      
    }
    outputPlots->Close();
  }


}
