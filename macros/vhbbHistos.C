#include <Compression.h>
#include <TCanvas.h>
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

#include <cassert>

#include "vhbbPlot.h"

const int nPlotsWlnHbb=10;

using namespace vhbbPlot;
void vhbbHistos(
  TString plotTreeFileName,
  TString outputFileName,
  vhbbPlot::selectionType selection,
  int theLepSel=-1, // -1: any; 0, e-mu; 1, mu only; 2, ele only
  bool isBlinded=true
) {
  int nPlots;
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) nPlots=nPlotsWlnHbb;

  system("mkdir -p MitVHBBAnalysis/plots");
  system("mkdir -p MitVHBBAnalysis/datacards");
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ");
  assert(plotTreeFile && plotTreeFile->IsOpen());
  
  // Initialize the tree
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);
  unsigned char typeLepSel, theCategory, nJet;
  unsigned selectionBits, nMinusOneBits;
  float pfmet, pfmetphi, pfmetsig;
  float lepton1Pt, lepton1Eta, lepton1Phi;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi;
  float vectorBosonPt, hbbDijetPt, hbbDijetMass;
  float bDiscrMin, bDiscrMax;
  float deltaPhiMetLep1, deltaPhiVH;
  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_scaleUp, weight_scaleDown;
  float weight_lepSFUp;
  float weight_cmvaUp, weight_cmvaDown; 
  
  plotTree->SetBranchAddress("selectionBits"   , &selectionBits   );
  plotTree->SetBranchAddress("nMinusOneBits"   , &nMinusOneBits   );
  plotTree->SetBranchAddress("theCategory"     , &theCategory     );
  plotTree->SetBranchAddress("nJet"            , &nJet            );
  plotTree->SetBranchAddress("pfmet"           , &pfmet           );
  plotTree->SetBranchAddress("pfmetsig"        , &pfmetsig        );
  plotTree->SetBranchAddress("pfmetphi"        , &pfmetphi        );
  plotTree->SetBranchAddress("hbbDijetPt"      , &hbbDijetPt      );
  plotTree->SetBranchAddress("hbbDijetMass"    , &hbbDijetMass    );
  plotTree->SetBranchAddress("bDiscrMin"       , &bDiscrMin       );
  plotTree->SetBranchAddress("bDiscrMax"       , &bDiscrMax       );
  plotTree->SetBranchAddress("hbbJet1Pt"       , &hbbJet1Pt       );
  plotTree->SetBranchAddress("hbbJet1Eta"      , &hbbJet1Eta      );
  plotTree->SetBranchAddress("hbbJet1Phi"      , &hbbJet1Phi      );
  plotTree->SetBranchAddress("hbbJet2Pt"       , &hbbJet2Pt       );
  plotTree->SetBranchAddress("hbbJet2Eta"      , &hbbJet2Eta      );
  plotTree->SetBranchAddress("hbbJet2Phi"      , &hbbJet2Phi      );
  plotTree->SetBranchAddress("weight"          , &weight          );
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    plotTree->SetBranchAddress("typeLepSel"      , &typeLepSel      );
    plotTree->SetBranchAddress("lepton1Pt"       , &lepton1Pt       );
    plotTree->SetBranchAddress("lepton1Eta"      , &lepton1Eta      );
    plotTree->SetBranchAddress("lepton1Phi"      , &lepton1Phi      );
    plotTree->SetBranchAddress("vectorBosonPt"   , &vectorBosonPt   );
    plotTree->SetBranchAddress("deltaPhiMetLep1" , &deltaPhiMetLep1 );
    plotTree->SetBranchAddress("deltaPhiVH"      , &deltaPhiVH      );
    plotTree->SetBranchAddress("weight_pdfUp"    , &weight_pdfUp    );
    plotTree->SetBranchAddress("weight_pdfDown"  , &weight_pdfDown  );
    plotTree->SetBranchAddress("weight_scaleUp"  , &weight_scaleUp  );
    plotTree->SetBranchAddress("weight_scaleDown", &weight_scaleDown);
    plotTree->SetBranchAddress("weight_lepSFUp"  , &weight_lepSFUp  );
    plotTree->SetBranchAddress("weight_cmvaUp"   , &weight_cmvaUp   );
    plotTree->SetBranchAddress("weight_cmvaDown" , &weight_cmvaDown );
  }
  
  // construct histograms
  TString histoTitle; float xmin, xmax; int nbins;
  TH1D *histos[99][nPlotCategories];
  TString histoNames[99];
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) for(int thePlot=0; thePlot<nPlotsWlnHbb; thePlot++) {
    switch(thePlot) {
      case  0: histoNames[thePlot]="pTH"             ; histoTitle="Higgs daughter dijet p_{T} [GeV]"           ; nbins=  18; xmin=    50; xmax=   350; break; 
      case  1: histoNames[thePlot]="mH"              ; histoTitle="Higgs daughter dijet mass [GeV]"            ; nbins=  25; xmin=     0; xmax=   250; break; 
      case  2: histoNames[thePlot]="WpT"             ; histoTitle="W boson p_{T} [GeV]"                        ; nbins=  18; xmin=    50; xmax=   350; break; 
      case  3: histoNames[thePlot]="deltaPhiVH"      ; histoTitle="#Delta#phi(H,W) [Rad]"                      ; nbins=  32; xmin=     0; xmax= 3.142; break; 
      case  4: histoNames[thePlot]="pTBalanceDijetW" ; histoTitle="Higgs daughter dijet p_{T} / W boson p_{T}" ; nbins=  30; xmin=     0; xmax=     2; break; 
      case  5: histoNames[thePlot]="lep1Pt"          ; histoTitle="Lepton p_{T} [GeV]"                         ; nbins=  23; xmin=    20; xmax=   250; break; 
      case  6: histoNames[thePlot]="met"             ; histoTitle="p_{T}^{miss} [GeV]"                         ; nbins=  25; xmin=     0; xmax=   250; break; 
      case  7: histoNames[thePlot]="mTW"             ; histoTitle="W transverse mass [GeV]"                    ; nbins=  20; xmin=     0; xmax=   200; break; 
      case  8: histoNames[thePlot]="Hbjet1Pt"        ; histoTitle="Leading H daughter jet p_{T} [GeV]"         ; nbins=  38; xmin=    20; xmax=   400; break; 
      case  9: histoNames[thePlot]="Hbjet2Pt"        ; histoTitle="Trailing H daughter jet p_{T} [GeV]"        ; nbins=  38; xmin=    20; xmax=   400; break; 
      default: throw std::runtime_error("bad plot number"); break;
    }
    for(theCategory=0; theCategory<nPlotCategories; theCategory++) {
      histos[thePlot][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitle, nbins, xmin, xmax);
      histos[thePlot][theCategory]->Sumw2();
      histos[thePlot][theCategory]->SetDirectory(0);
    }
  }

  Long64_t nentries = plotTree->GetEntries();
  Long64_t oneTenth = nentries/10;
  for(Long64_t ientry=0; ientry<nentries; ientry++) {
    if(ientry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries); 
    plotTree->GetBranch("typeLepSel")->GetEntry(ientry); if(theLepSel!=-1 && typeLepSel!=theLepSel) continue;
    plotTree->GetBranch("nMinusOneBits")->GetEntry(ientry); if((nMinusOneBits & selection)==0) continue;
    plotTree->GetBranch("theCategory")->GetEntry(ientry); if(theCategory>=nPlotCategories) continue;
    if(isBlinded && theCategory==kPlotData && (selection==kWHSR)) continue;
    plotTree->GetEntry(ientry);
    float mT = TMath::Sqrt(2.*pfmet*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-pfmetphi))));
    float balanceVH = hbbDijetPt/vectorBosonPt;
    float theVar; bool passFullSel = (selectionBits & selection)!=0; bool makePlot;
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) for(int thePlot=0; thePlot<nPlotsWlnHbb; thePlot++) {
      // Variables -- change the makePlot for n-1 later
      switch(thePlot) {
        case  0: makePlot = passFullSel; theVar = TMath::Min(hbbDijetPt   , (float)  350); break; 
        case  1: makePlot = passFullSel; theVar = TMath::Min(hbbDijetMass , (float)  250); break; 
        case  2: makePlot = passFullSel; theVar = TMath::Min(vectorBosonPt, (float)  350); break; 
        case  3: makePlot = passFullSel; theVar = TMath::Min(deltaPhiVH   , (float)3.142); break; 
        case  4: makePlot = passFullSel; theVar = TMath::Min(balanceVH    , (float)    2); break; 
        case  5: makePlot = passFullSel; theVar = TMath::Min(lepton1Pt    , (float)  250); break; 
        case  6: makePlot = passFullSel; theVar = TMath::Min(pfmet        , (float)  250); break; 
        case  7: makePlot = passFullSel; theVar = TMath::Min(mT           , (float)  200); break; 
        case  8: makePlot = passFullSel; theVar = TMath::Min(hbbJet1Pt    , (float)  400); break; 
        case  9: makePlot = passFullSel; theVar = TMath::Min(hbbJet2Pt    , (float)  400); break; 
        default: throw std::runtime_error("bad plot number"); break;
      }
      if(makePlot) histos[thePlot][theCategory]->Fill(theVar, weight);
    }
  }
  TFile *output_plots = new TFile(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  for(int thePlot=0; thePlot<nPlots; thePlot++) {
    TDirectory *plotDir = output_plots->mkdir(histoNames[thePlot]);
    plotDir->cd();
    for(theCategory=kPlotData; theCategory!=nPlotCategories; theCategory++)
      histos[thePlot][theCategory]->Write();

  }
  output_plots->Close();

}

