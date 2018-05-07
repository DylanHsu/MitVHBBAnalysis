#include <Compression.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom3.h>

#include <cassert>
#include <iostream>
#include <fstream>

#include "vhbbPlot.h"
using namespace vhbbPlot;
void measureDoubleBEfficiency(TString inputFile="/data/t3home000/dhsu/whbbPlots/root/plotTreeBoosted.root") {
  TFile *plotTreeFile = TFile::Open(inputFile,"READ");
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);

  unsigned long int time_now = static_cast<unsigned long int>(time(NULL));
  unsigned int randomToySeed=(time_now-731178000); // random seed based on Dylan's age in seconds
  TRandom3 toymaker(randomToySeed);
  
  int runNumber = -1;
  int lumiNumber = -1;
  ULong64_t eventNumber = -1;
  unsigned char typeLepSel, theCategory;
  unsigned char isojetNBtags;
  unsigned selectionBits;
  float topWBosonPt;
  float weight;
  float fj1MSD;
  float fj1Pt;
  float fj1DoubleCSV;
  plotTree->SetBranchStatus("*",0);
  //plotTree->SetBranchStatus("runNumber"    ,1); plotTree->SetBranchAddress("runNumber"    , &runNumber    );
  //plotTree->SetBranchStatus("lumiNumber"   ,1); plotTree->SetBranchAddress("lumiNumber"   , &lumiNumber   );
  //plotTree->SetBranchStatus("eventNumber"  ,1); plotTree->SetBranchAddress("eventNumber"  , &eventNumber  );
  plotTree->SetBranchStatus("selectionBits",1); plotTree->SetBranchAddress("selectionBits", &selectionBits);
  plotTree->SetBranchStatus("theCategory"  ,1); plotTree->SetBranchAddress("theCategory"  , &theCategory  );
  plotTree->SetBranchStatus("weight"       ,1); plotTree->SetBranchAddress("weight"       , &weight       );
  plotTree->SetBranchStatus("topWBosonPt"  ,1); plotTree->SetBranchAddress("topWBosonPt"  , &topWBosonPt  );
  plotTree->SetBranchStatus("typeLepSel"   ,1); plotTree->SetBranchAddress("typeLepSel"   , &typeLepSel   );
  plotTree->SetBranchStatus("isojetNBtags" ,1); plotTree->SetBranchAddress("isojetNBtags" , &isojetNBtags );
  plotTree->SetBranchStatus("fj1MSD_corr"  ,1); plotTree->SetBranchAddress("fj1MSD_corr"  , &fj1MSD_corr  );
  plotTree->SetBranchStatus("fj1Pt"        ,1); plotTree->SetBranchAddress("fj1Pt"        , &fj1Pt        );
  plotTree->SetBranchStatus("fj1DoubleCSV" ,1); plotTree->SetBranchAddress("fj1DoubleCSV" , &fj1DoubleCSV );   
  
  //////////////////////////////////////////////////////////////////////////// 
  //                                //                                      //
  //  B: fail 2B, 1+ iso b tags     //  D: pass 2B, 1+ iso b tags           //
  //                                //                                      //
  //                                //                                      //
  //////////////////////////////////////////////////////////////////////////// 
  //                                //                                      //
  //                                //                                      //
  //  A: fail 2B, 0 iso b tags      //  C: pass 2B, 0 iso b tags            //
  //                                //                                      //
  //////////////////////////////////////////////////////////////////////////// 
  unsigned char processIdx, regionIdx;
  float yield[4][4], sumw2[4][4]; unsigned N[4][4];
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) { yield[i][j]=0; sumw2[i][j]=0; N[i][j]=0;}
  // first index is process type: 0=TT,1=Wbb,2=Wb,3=WLF
  // second index is region: 0=A, 1=B, 2=C, 3=D
  
  Long64_t nentries = plotTree->GetEntries();
  Long64_t oneTenth = nentries/10;
  for(Long64_t ientry=0; ientry<nentries; ientry++) {
    if(ientry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries); 
    plotTree->GetBranch("selectionBits")->GetEntry(ientry);
    unsigned selection = kWHLightFlavorFJCR | kWHHeavyFlavorFJCR | kWHTT2bFJCR | kWHTT1bFJCR | kWHFJSR; 
    if((selectionBits&selection)==0) continue;
    plotTree->GetBranch("theCategory")->GetEntry(ientry); 
    if(theCategory!=kPlotTT && theCategory!=kPlotWbb && theCategory!=kPlotWb && theCategory!=kPlotWLF) continue;
    plotTree->GetBranch("fj1MSD_corr")->GetEntry(ientry);
    if(fj1MSD_corr<40) continue;
    plotTree->GetEntry(ientry);

    if     (theCategory==kPlotTT ) processIdx=0; 
    else if(theCategory==kPlotWbb) processIdx=1; 
    else if(theCategory==kPlotWb ) processIdx=2; 
    else if(theCategory==kPlotWLF) processIdx=3; 
    else continue;

    if     ((selectionBits & kWHLightFlavorFJCR)!=0   ) regionIdx=0; //A
    else if((selectionBits & kWHTT1bFJCR       )!=0   ) regionIdx=1; //B
    else if((selectionBits & kWHHeavyFlavorFJCR)!=0 ||  
            (selectionBits & kWHFJSR           )!=0   ) regionIdx=2; //C
    else if((selectionBits & kWHTT2bFJCR       )!=0   ) regionIdx=3; //D
    else continue;
    //float wptCorrFactor=theWptCorr->Eval(TMath::Min(topWBosonPt,(float)499.99));
    //weight *= wptCorrFactor;
    yield[processIdx][regionIdx] += weight;
    sumw2[processIdx][regionIdx] += weight*weight;
    N    [processIdx][regionIdx]++;
  }
  float rms[4][4];
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) {
    :q
    vi//rms[i][j] = sqrt(sumw2[i][j]/(float)N[i][j]);
    rms[i][j] = sqrt(sumw2[i][j]);
    //rms[i][j] = yield[i][j]/sqrt(float(N[i][j]));
  }
  float eff[4][2], err[4][2];
  float toyPass,toyFail,toyNPass,toyNFail,NPass,NFail,NsigErrPass,NsigErrFail;
  TH1F *toyEffs = new TH1F("toyEffs","toyEffs", 1000, 0,1);
  for(int i=0;i<4;i++) {
    // 0 isojet btags
    eff[i][0] = yield[i][2]/(yield[i][2]+yield[i][0]);
    toyEffs->Reset(); toyEffs->Clear();
    NPass=yield[i][2]; NFail=yield[i][0];
    NsigErrPass = rms[i][2];
    NsigErrFail = rms[i][0];
    for(int iToy=0; iToy<1000; iToy++) {
      toyPass=toymaker.Gaus(-1,1);
      toyFail=toymaker.Gaus(-1,1);
      toyNPass = TMath::Max(float(0.),NPass + toyPass * NsigErrPass);
      toyNFail = TMath::Max(float(0.),NFail + toyFail * NsigErrFail);
      toyEffs->Fill(toyNPass/(toyNPass+toyNFail));
    }
    err[i][0] = toyEffs->GetStdDev();
    // 1+ isojet btags
    eff[i][1] = yield[i][3]/(yield[i][3]+yield[i][1]);
    toyEffs->Reset(); toyEffs->Clear();
    NPass=yield[i][3]; NFail=yield[i][1];
    NsigErrPass = rms[i][3];
    NsigErrFail = rms[i][1];
    for(int iToy=0; iToy<1000; iToy++) {
      toyPass=toymaker.Gaus(-1,1);
      toyFail=toymaker.Gaus(-1,1);
      toyNPass = TMath::Max(float(0.),NPass + toyPass * NsigErrPass);
      toyNFail = TMath::Max(float(0.),NFail + toyFail * NsigErrFail);
      toyEffs->Fill(toyNPass/(toyNPass+toyNFail));
    }
    err[i][1] = toyEffs->GetStdDev();
  }
  printf("%16s%16s%8s%16s%8s\n","Process","Eff.(0ijB)","+/-","Eff.(1+ijB)","+/-");
  printf("-----------------------------------------------------------------------------------\n");
  printf("%16s%16.3f%8.5f%16.3f%8.5f\n","ttbar"  ,eff[0][0],err[0][0],eff[0][1],err[0][1]);
  printf("%16s%16.3f%8.5f%16.3f%8.5f\n","W+bb"   ,eff[1][0],err[1][0],eff[1][1],err[1][1]);
  printf("%16s%16.3f%8.5f%16.3f%8.5f\n","W+b"    ,eff[2][0],err[2][0],eff[2][1],err[2][1]);
  printf("%16s%16.3f%8.5f%16.3f%8.5f\n","W+udcsg",eff[3][0],err[3][0],eff[3][1],err[3][1]);

 


}


