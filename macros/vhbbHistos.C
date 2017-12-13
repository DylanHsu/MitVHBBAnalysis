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

const int nPlotsWlnHbb=15;
const int MVAVarType=1; //1=mass; 2=BDT

using namespace vhbbPlot;
void vhbbHistos(
  TString plotTreeFileName,
  TString outputFileName,
  vhbbPlot::selectionType selection,
  int theLepSel=-1, // -1: any; 0, e-mu; 1, mu only; 2, ele only
  bool isBlinded=true
) {
  vector<double> MVAbins; TString MVAVarName;
  if(MVAVarType==1) {
    MVAVarName="Higgs mass classifier [GeV]";
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) 
      MVAbins={0,60,90,120,150,200,250};
    if(selection==kWHSR) 
      MVAbins={90,100,110,120,130,140,150};
  } else throw std::runtime_error("bad MVAVarType");
  
  int nPlots;
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) nPlots=nPlotsWlnHbb;
  nPlots=99;

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
  vector<float> xmin(nPlots), xmax(nPlots);  vector<int> nbins(nPlots);
  vector<TString> histoNames(nPlots), histoTitle(nPlots);
  
  TH1D *histos[99][nPlotCategories];
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    int p=0;
    histoNames[p]="pTH"             ; histoTitle[p]="Higgs daughter dijet p_{T} [GeV]"           ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="mH"              ; histoTitle[p]="Higgs daughter dijet mass [GeV]"            ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="WpT"             ; histoTitle[p]="W boson p_{T} [GeV]"                        ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="deltaPhiVH"      ; histoTitle[p]="#Delta#phi(H,W) [Rad]"                      ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="pTBalanceDijetW" ; histoTitle[p]="Higgs daughter dijet p_{T} / W boson p_{T}" ; nbins[p]=  30; xmin[p]=     0; xmax[p]=     2; p++; 
    histoNames[p]="lepton1Pt"       ; histoTitle[p]="Lepton p_{T} [GeV]"                         ; nbins[p]=  23; xmin[p]=    20; xmax[p]=   250; p++; 
    histoNames[p]="met"             ; histoTitle[p]="p_{T}^{miss} [GeV]"                         ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="mTW"             ; histoTitle[p]="W transverse mass [GeV]"                    ; nbins[p]=  20; xmin[p]=     0; xmax[p]=   200; p++; 
    histoNames[p]="Hbjet1Pt"        ; histoTitle[p]="Leading H daughter jet p_{T} [GeV]"         ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="Hbjet2Pt"        ; histoTitle[p]="Trailing H daughter jet p_{T} [GeV]"        ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="dPhiHbbJet1MET"  ; histoTitle[p]="#Delta#phi(H daughter jet 1, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhiHbbJet2MET"  ; histoTitle[p]="#Delta#phi(H daughter jet 2, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="bDiscrMin"       ; histoTitle[p]="Greater CMVA of Higgs daughter jets"        ; nbins[p]=  20; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="bDiscrMax"       ; histoTitle[p]="Lesser CMVA of Higgs daughters jets"        ; nbins[p]=  20; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    //case nPlotsWlnHbb-1: histoNames[p]="MVAVar"; histoTitle=MVAVarName; break;
  }  
  
  for(int p=0; p<nPlots; p++) {
    if(histoNames[p]=="") continue;
    for(theCategory=0; theCategory<nPlotCategories; theCategory++) {
      //if(thePlot==nPlotsWlnHbb-1)
      //  histos[thePlot][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitle, MVAbins.size()-1, MVAbins.data());
      //else
        histos[p][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitle[p], nbins[p], xmin[p], xmax[p]);
      histos[p][theCategory]->Sumw2();
      histos[p][theCategory]->SetDirectory(0);
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
    // physics calculations
    float mT = TMath::Sqrt(2.*pfmet*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-pfmetphi))));
    float balanceVH = hbbDijetPt/vectorBosonPt;
    float dPhiHbbJet1MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet1Phi - pfmetphi  ));
    float dPhiHbbJet2MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet2Phi - pfmetphi  ));

    float theVar; bool passFullSel = (selectionBits & selection)!=0; bool makePlot;
    // hack
    //if(typeLepSel==1 && theCategory==kPlotTop) weight*=-1.;
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) for(int p=0; p<nPlots; p++) {
      makePlot=false;
      // Variables -- change the makePlot for n-1 later
      if      (histoNames[p]=="pTH"             ) { theVar = TMath::Min(hbbDijetPt    , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="mH"              ) { theVar = TMath::Min(hbbDijetMass  , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="WpT"             ) { theVar = TMath::Min(vectorBosonPt , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="deltaPhiVH"      ) { theVar = TMath::Min(deltaPhiVH    , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar>=2.5); }
      else if (histoNames[p]=="pTBalanceDijetW" ) { theVar = TMath::Min(balanceVH     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Pt"       ) { theVar = TMath::Min(lepton1Pt     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmet"           ) { theVar = TMath::Min(pfmet         , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetsig"        ) { theVar = TMath::Min(pfmet         , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=2.); }
      else if (histoNames[p]=="mTW"             ) { theVar = TMath::Min(mT            , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet1Pt"        ) { theVar = TMath::Min(hbbJet1Pt     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet2Pt"        ) { theVar = TMath::Min(hbbJet2Pt     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet1MET"  ) { theVar = TMath::Min(dPhiHbbJet1MET, xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet2MET"  ) { theVar = TMath::Min(dPhiHbbJet2MET, xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="bDiscrMin"       ) { theVar = TMath::Min(bDiscrMin     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="bDiscrMax"       ) { theVar = TMath::Min(bDiscrMax     , xmax[p]); makePlot = passFullSel; }
      else continue;
      if(makePlot) histos[p][theCategory]->Fill(theVar, weight);
    }
  }
  TFile *output_plots = new TFile(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  for(int p=0; p<nPlots; p++) {
    if(histoNames[p]=="") continue;
    TDirectory *plotDir = output_plots->mkdir(histoNames[p]);
    plotDir->cd();
    for(theCategory=kPlotData; theCategory!=nPlotCategories; theCategory++)
      histos[p][theCategory]->Write();

  }
  output_plots->Close();

}

