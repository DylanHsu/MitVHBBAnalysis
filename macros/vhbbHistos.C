#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
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
#include <iostream>
#include <fstream>

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
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) MVAbins={0,30,50,70,90,100,110,120,130,140,150,170,190,210,230,250};
    else if(selection==kWHSR)                                                              MVAbins={90,100,110,120,130,140,150};
    else throw std::runtime_error("bad selection argument");
  } else throw std::runtime_error("bad MVAVarType");
  
  // Load Files
  gSystem->Load("libCondFormatsJetMETObjects.so");
  system("mkdir -p MitVHBBAnalysis/datacards");
  system("mkdir -p MitVHBBAnalysis/plots");
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ");
  assert(plotTreeFile && plotTreeFile->IsOpen());
  string theJecUncertainties = "PandaAnalysis/data/jec/23Sep2016V4/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt";
  JetCorrectionUncertainty *jcuTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(theJecUncertainties, "Total")));
  
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
  int nPlots=99;
  vector<float> xmin(nPlots), xmax(nPlots);  vector<int> nbins(nPlots);
  vector<TString> histoNames(nPlots), histoTitles(nPlots);
  
  int pMVAVar=-1;
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    int p=0;
    histoNames[p]="pTH"             ; histoTitles[p]="Higgs daughter dijet p_{T} [GeV]"           ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="mH"              ; histoTitles[p]="Higgs daughter dijet mass [GeV]"            ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="WpT"             ; histoTitles[p]="W boson p_{T} [GeV]"                        ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="deltaPhiVH"      ; histoTitles[p]="#Delta#phi(H,W) [Rad]"                      ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="pTBalanceDijetW" ; histoTitles[p]="Higgs daughter dijet p_{T} / W boson p_{T}" ; nbins[p]=  30; xmin[p]=     0; xmax[p]=     2; p++; 
    histoNames[p]="lepton1Pt"       ; histoTitles[p]="Lepton p_{T} [GeV]"                         ; nbins[p]=  23; xmin[p]=    20; xmax[p]=   250; p++; 
    histoNames[p]="pfmet"           ; histoTitles[p]="p_{T}^{miss} [GeV]"                         ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="mTW"             ; histoTitles[p]="W transverse mass [GeV]"                    ; nbins[p]=  20; xmin[p]=     0; xmax[p]=   200; p++; 
    histoNames[p]="Hbjet1Pt"        ; histoTitles[p]="Leading H daughter jet p_{T} [GeV]"         ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="Hbjet2Pt"        ; histoTitles[p]="Trailing H daughter jet p_{T} [GeV]"        ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="dPhiHbbJet1MET"  ; histoTitles[p]="#Delta#phi(H daughter jet 1, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhiHbbJet2MET"  ; histoTitles[p]="#Delta#phi(H daughter jet 2, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="bDiscrMin"       ; histoTitles[p]="Lesser CMVA of Higgs daughter jets"         ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="bDiscrMax"       ; histoTitles[p]="Greater CMVA of Higgs daughter jets"        ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="MVAVar"          ; histoTitles[p]=MVAVarName; pMVAVar=p; p++;
    //case nPlotsWlnHbb-1: histoNames[p]="MVAVar"; histoTitles=MVAVarName; break;
  }  
  assert(pMVAVar!=-1);

  TH1D *histos[99][nPlotCategories];
  for(int p=0; p<nPlots; p++) {
    if(histoNames[p]=="") continue;
    for(theCategory=0; theCategory<nPlotCategories; theCategory++) {
      if(histoNames[p]=="MVAVar") {
        histos[p][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitles[p], MVAbins.size()-1, MVAbins.data());
      } else
        histos[p][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitles[p], nbins[p], xmin[p], xmax[p]);
      histos[p][theCategory]->Sumw2();
      histos[p][theCategory]->SetDirectory(0);
    }
  }
  TH1D *histo_pdfUp    [nPlotCategories];
  TH1D *histo_pdfDown  [nPlotCategories];
  TH1D *histo_scaleUp  [nPlotCategories];
  TH1D *histo_scaleDown[nPlotCategories];
  TH1D *histo_eleSFUp  [nPlotCategories];
  TH1D *histo_muSFUp   [nPlotCategories];
  TH1D *histo_cmvaUp   [nPlotCategories];
  TH1D *histo_cmvaDown [nPlotCategories];
  TH1D *histo_jesUp   [nPlotCategories];
  TH1D *histo_jesDown [nPlotCategories];
  for(theCategory=0; theCategory<nPlotCategories; theCategory++) {
    histo_pdfUp    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_pdfUp"    , theCategory)); histo_pdfUp    [theCategory]->SetDirectory(0);
    histo_pdfDown  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_pdfDown"  , theCategory)); histo_pdfDown  [theCategory]->SetDirectory(0);
    histo_scaleUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_scaleUp"  , theCategory)); histo_scaleUp  [theCategory]->SetDirectory(0);
    histo_scaleDown[theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_scaleDown", theCategory)); histo_scaleDown[theCategory]->SetDirectory(0);
    histo_eleSFUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_eleSFUp"  , theCategory)); histo_eleSFUp  [theCategory]->SetDirectory(0);
    histo_muSFUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_muSFUp"   , theCategory)); histo_muSFUp   [theCategory]->SetDirectory(0);
    histo_cmvaUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaUp"   , theCategory)); histo_cmvaUp   [theCategory]->SetDirectory(0);
    histo_cmvaDown [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaDown" , theCategory)); histo_cmvaDown [theCategory]->SetDirectory(0);
    histo_jesUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_jesUp"   , theCategory)); histo_jesUp   [theCategory]->SetDirectory(0);
    histo_jesDown [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_jesDown" , theCategory)); histo_jesDown [theCategory]->SetDirectory(0);
  }
  // done constructing histograms

  // begin plot tree loop
  Long64_t nentries = plotTree->GetEntries();
  Long64_t oneTenth = nentries/10;
  for(Long64_t ientry=0; ientry<nentries; ientry++) {
    if(ientry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries); 
    plotTree->GetBranch("typeLepSel")->GetEntry(ientry); if(theLepSel!=-1 && typeLepSel!=theLepSel) continue;
    plotTree->GetBranch("nMinusOneBits")->GetEntry(ientry); if((nMinusOneBits & selection)==0) continue;
    plotTree->GetBranch("theCategory")->GetEntry(ientry); if(theCategory>=nPlotCategories) continue;
    if(isBlinded && theCategory==kPlotData && (selection==kWHSR)) continue;
    plotTree->GetEntry(ientry);
    if(weight>500.) continue;
    // physics calculations
    float mT = TMath::Sqrt(2.*pfmet*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-pfmetphi))));
    float balanceVH = hbbDijetPt/vectorBosonPt;
    float dPhiHbbJet1MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet1Phi - pfmetphi  ));
    float dPhiHbbJet2MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet2Phi - pfmetphi  ));
    
    float weight_jesUp=weight, weight_jesDown=weight;
    jcuTotal->setJetPt(hbbJet1Pt); jcuTotal->setJetEta(hbbJet1Eta);
    weight_jesUp   *= (1+jcuTotal->getUncertainty(true));
    jcuTotal->setJetPt(hbbJet1Pt); jcuTotal->setJetEta(hbbJet1Eta);
    weight_jesDown *= (1-jcuTotal->getUncertainty(true));
    jcuTotal->setJetPt(hbbJet2Pt); jcuTotal->setJetEta(hbbJet2Eta);
    weight_jesUp   *= (1+jcuTotal->getUncertainty(false));
    jcuTotal->setJetPt(hbbJet2Pt); jcuTotal->setJetEta(hbbJet2Eta);
    weight_jesDown *= (1-jcuTotal->getUncertainty(false));
    
    float MVAVar;
    if(MVAVarType==1) MVAVar=hbbDijetMass;
    float theVar; bool passFullSel = (selectionBits & selection)!=0; bool makePlot;
    // hack
    //if(typeLepSel==1 && theCategory==kPlotTop) weight*=-1.;
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) for(int p=0; p<nPlots; p++) {
      makePlot=false;
      // Variables -- change the makePlot for n-1 later
      if      (histoNames[p]=="pTH"             ) { theVar = TMath::Min(hbbDijetPt    , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="mH"              ) { theVar = TMath::Min(hbbDijetMass  , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="WpT"             ) { theVar = TMath::Min(vectorBosonPt , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="deltaPhiVH"      ) { theVar = TMath::Min(deltaPhiVH    , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<2.5); }
      else if (histoNames[p]=="pTBalanceDijetW" ) { theVar = TMath::Min(balanceVH     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Pt"       ) { theVar = TMath::Min(lepton1Pt     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmet"           ) { theVar = TMath::Min(pfmet         , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetsig"        ) { theVar = TMath::Min(pfmet         , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=2.); }
      else if (histoNames[p]=="mTW"             ) { theVar = TMath::Min(mT            , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet1Pt"        ) { theVar = TMath::Min(hbbJet1Pt     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet2Pt"        ) { theVar = TMath::Min(hbbJet2Pt     , xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet1MET"  ) { theVar = TMath::Min(dPhiHbbJet1MET, xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet2MET"  ) { theVar = TMath::Min(dPhiHbbJet2MET, xmax[p]); makePlot = passFullSel; }
      else if (histoNames[p]=="bDiscrMin"       ) { theVar = TMath::Min(bDiscrMin     , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && (selection==kWHSR && theVar<bDiscrLoose )); }
      else if (histoNames[p]=="bDiscrMax"       ) { theVar = TMath::Min(bDiscrMax     , xmax[p]); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && ((selection==kWHLightFlavorCR && (theVar<bDiscrLoose || theVar>bDiscrMedium) ) || (selection!=kWHLightFlavorCR && theVar<bDiscrTight)) ); }
      else if (histoNames[p]=="MVAVar"          ) { theVar = MVAVar                             ; makePlot = passFullSel; }
      else continue;
      if(!makePlot) continue;
      histos[p][theCategory]->Fill(theVar, weight);
    }
    // Fill systematic shapes
    if(passFullSel && theCategory!=kPlotData) {
      histo_pdfUp    [theCategory]->Fill(MVAVar, weight_pdfUp    ); 
      histo_pdfDown  [theCategory]->Fill(MVAVar, weight_pdfDown  ); 
      histo_scaleUp  [theCategory]->Fill(MVAVar, weight_scaleUp  ); 
      histo_scaleDown[theCategory]->Fill(MVAVar, weight_scaleDown); 
      histo_muSFUp   [theCategory]->Fill(MVAVar, typeLepSel==1?weight_lepSFUp  : weight); 
      histo_eleSFUp  [theCategory]->Fill(MVAVar, typeLepSel==2?weight_lepSFUp : weight); 
      histo_cmvaUp   [theCategory]->Fill(MVAVar, weight_cmvaUp   ); 
      histo_cmvaDown [theCategory]->Fill(MVAVar, weight_cmvaDown );
      histo_jesUp   [theCategory]->Fill(MVAVar, weight_jesUp   ); 
      histo_jesDown [theCategory]->Fill(MVAVar, weight_jesDown );
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

  // Prepare the datacard C-('.' Q)
  char leptonChar='l';
  if(theLepSel==1) leptonChar='m';
  if(theLepSel==2) leptonChar='e';
  
  for(int nb=1; nb<=histos[pMVAVar][kPlotData]->GetNbinsX(); nb++) {
    float bgSum=0;
    for(int ic=0; ic<nPlotCategories; ic++)
      if(ic!=kPlotVH) bgSum+=histos[pMVAVar][ic]->GetBinContent(nb);
    if(bgSum<1e-3 && histos[pMVAVar][kPlotVH]->GetBinContent(nb)<1e-3) continue;
    char outputLimitsShape[512];
    //adapt the "shape" term later
    sprintf(outputLimitsShape,"MitVHBBAnalysis/datacards/histo_limits_W%cnHbb_%s_massShape_bin%d.txt",leptonChar,selectionNames[selection].Data(), nb-1);
    ofstream card; card.open(outputLimitsShape);
    card << Form("imax 1 number of channels\n");
    card << Form("jmax * number of background\n");
    card << Form("kmax * number of nuisance parameters\n");
    card << Form("Observation %d\n", isBlinded ? 0 : (int)round(histos[pMVAVar][kPlotData]->GetBinContent(nb))); 
    card << Form("bin     "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9s ",Form("W%cn%db%d",leptonChar,selection,nb-1)); card<<std::endl;
    card << Form("process "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9s ",plotBaseNames[ic].Data()); card<<std::endl;
    card << Form("process "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9d ", ic==kPlotVH?0:ic); card<<std::endl;
    card << Form("rate    "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9.3f ", histos[pMVAVar][ic]->GetBinContent(nb)); card<<std::endl;
    // Systematics - Shape Uncertainties
    // PDF Acceptance
    card<<Form("pdf_qqbar_ACCEPT lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_pdfDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1,
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_pdfUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;
    // QCD Scale
    for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
      card<<Form("QCDscale_%s lnN ",plotBaseNames[ic].Data());
      for(int jc=kPlotData+1; jc<nPlotCategories; jc++) 
        card << (jc!=ic? "- ":Form("%7.5f/%7.5f ",
          histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_scaleDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1,
          histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_scaleUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
        ));
      card<<std::endl;
    }
    // BTag Shape Uncertainty
    card<<Form("CMS_btagShape lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,histo_cmvaDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1,
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_cmvaUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;
    // Total Jet Energy Scale Shape Uncertainty
    card<<Form("CMS_scale_j lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_jesDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1,
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_jesUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;
    // Muon SF Shape
    card<<Form("CMS_eff2016_m lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_muSFUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;
    // Electron SF Shape
    card<<Form("CMS_eff2016_e lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_eleSFUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;

    // Systematics -- Normalization Uncertainties
    // Muon Energy Scale Norm
    card<<Form("CMS_eff2016_m lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.01 "; card <<std::endl;
    // Electron Energy Scale Norm
    card<<Form("CMS_eff2016_e lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.01 "; card<<std::endl;
    // Trigger
    card<<Form("CMS_trigger2016 lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.01 "; card << std::endl;
    // Lumi
    card<<Form("lumi_13TeV2016 lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.025 "; card << std::endl;
    // VH Cross Section
    card<<Form("pdf_qqbar lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.04 "; card << std::endl;
    // Single Top Cross Section
    card<<Form("CMS_VH_TopNorm lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< (ic==kPlotTop? "1.30 ":"- "); card<<std::endl;
    // Diboson Cross Section
    card<<Form("CMS_VH_VVNorm lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< ((ic==kPlotVZbb||ic==kPlotVVLF)? "1.30 ":"- "); card<<std::endl;

    // Statistical Uncertainties
    for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
      card<<Form("W%cn%s_%sStatBounding_13TeV2016_Bin%d lnN ",leptonChar,selectionNames[selection].Data(), plotBaseNames[ic].Data(), nb-1);
      for(int jc=kPlotData+1; jc<nPlotCategories; jc++) 
        card << (jc!=ic? "- ":Form("%7.5f ",histos[pMVAVar][ic]->GetBinContent(nb)>0? 1.+histos[pMVAVar][ic]->GetBinError(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1));
      card<<std::endl;
    }
 
  } 
}

