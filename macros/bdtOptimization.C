#include <TROOT.h>
#include <TFile.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <TList.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TMath.h>
#include <map>
#include <vector>
#include <TPaletteAxis.h>
const unsigned nVarsMax=100;
//const unsigned nVarsMax=2;
const float ln2=0.3010;
//https://arxiv.org/pdf/1408.3122.pdf
void mutualInfTruth(
  TString inputFileName,
  //float f=0.00542, //WH: 0.00542
  TString trainingName="BDT_multiClass__dec22_test1",
  TString bdtClass="Signal",
  TString varTransform="Id",
  bool debug=false
) {
  TFile *inputFile = TFile::Open(inputFileName,"READ");
  assert(inputFile && inputFile->IsOpen());

  TString histoSuffix = "__" + bdtClass + "_" + varTransform;
  
  // Get input variable histograms
  printf("Loading pre-binned input variable histograms ...\n");
  TDirectory *histoDir = (TDirectory*)inputFile->Get("InputVariables_"+varTransform);
  assert(histoDir);
  TList *histoList = histoDir->GetListOfKeys();
  TH1F *histos[nVarsMax];
  TString varNames[nVarsMax], varTitles[nVarsMax];
  unsigned nVars; { unsigned h=0; for(unsigned i=0; i<=(unsigned)histoList->LastIndex(); i++) {
    TString histoName = histoList->At(i)->GetName();
    Ssiz_t pos = histoName.Index(histoSuffix);
    if(pos==-1) continue;
    varNames[h] = TString( histoName(0,pos));
    histos[h] = (TH1F*) histoDir->Get(histoName);
    assert(histos[h]);
    histos[h]->SetDirectory(0); histos[h]->SetName(varNames[h]);
    printf("Loaded variable #%d: \"%s\"\n", h, varNames[h].Data());
    varTitles[h] = histoList->At(i)->GetTitle();
    h++;
    if(h>=nVarsMax) break;
  } nVars=h; }
  printf("Done loading input variable histograms.\n");

  // Get BDT test histograms
  printf("Getting BDT test histogram ...\n");
  TDirectory *testHistoDir = (TDirectory*)inputFile->Get("Method_BDT/"+trainingName);
  assert(testHistoDir);
  //TList *testHistoList = testHistoDir->GetListOfKeys();
  // "MVA_BDT_multiClass__dec22_test1_Test_Signal_prob_for_Signal"
  TH1F *bdtTestHisto = (TH1F*)testHistoDir->Get("MVA_"+trainingName+"_Test_"+bdtClass+"_prob_for_"+bdtClass);
  assert(bdtTestHisto); bdtTestHisto->SetDirectory(0);
  printf("Done getting BDT test histogram.\n");

  // Get joint distributions of the variables among themselves and with the BDT (truth)
  // Z axis is truth info (1=signal, 0=background)
  TH3F *jointDistsAmongVars[nVarsMax][nVarsMax], *jointDistsVarsBdt[nVarsMax];
  // Need to use the TestTree for relative class normalizations
  // all the classes have the same sum of weights in the TrainTree
  // To do: get the relative norms from TestTree 
  // and use TrainTree with the norms in the cut string
  TTree *testTree=(TTree*)inputFile->Get("TestTree"); assert(testTree);
  printf("Computing joint distributions among variables ...\n");
  for(unsigned h1=0; h1<nVars; h1++) {
    if(h1+1<nVars) for(unsigned h2=h1+1; h2<nVars; h2++) {
      testTree->SetBranchStatus("*",0);
      testTree->SetBranchStatus("className",1);
      testTree->SetBranchStatus("weight",1);
      testTree->SetBranchStatus(varNames[h1],1);
      testTree->SetBranchStatus(varNames[h2],1);
      if(debug) printf("\tJoint distribution \"%s\" vs \"%s\":\n",varNames[h2].Data(), varNames[h1].Data());
      // Z axis is signal(0|1)
      jointDistsAmongVars[h1][h2] = new TH3F(varNames[h2]+"_vs_"+varNames[h1],varNames[h2]+"_vs_"+varNames[h1],
        histos[h1]->GetNbinsX(),
        histos[h1]->GetBinLowEdge(1),
        histos[h1]->GetBinLowEdge(histos[h1]->GetNbinsX()+1),
        histos[h2]->GetNbinsX(),
        histos[h2]->GetBinLowEdge(1),
        histos[h2]->GetBinLowEdge(histos[h2]->GetNbinsX()+1),
        2,0,2 );
      
      // documentation for TTree::Draw is wrong, I think: expression is Z:Y:X not X:Y:Z
      TString varExpr=Form("(className==\"%s\"):%s:%s>>+%s_vs_%s",bdtClass.Data(),varNames[h2].Data(),varNames[h1].Data(),varNames[h2].Data(), varNames[h1].Data());
      TString selection = "weight";
      testTree->Draw(varExpr, selection, "e goff");
      if(debug) {
        printf("\tTTree::Draw expression: \"%s\"\n", varExpr.Data());
        printf("\t\tSum of weights signal: %.2f\n", jointDistsAmongVars[h1][h2]->Integral( 1, jointDistsAmongVars[h1][h2]->GetNbinsX(), 1, jointDistsAmongVars[h1][h2]->GetNbinsY(), 2, 2));
        printf("\t\tSum of weights bkg: %.2f\n", jointDistsAmongVars[h1][h2]->Integral( 1, jointDistsAmongVars[h1][h2]->GetNbinsX(), 1, jointDistsAmongVars[h1][h2]->GetNbinsY(), 1, 1));
      }
    }
  }
  printf("Done computing joint distributions among variables.\n");
  
  printf("Computing signal fraction ...\n");
  testTree->SetBranchStatus("*",0);
  testTree->SetBranchStatus("className",1);
  testTree->SetBranchStatus("weight",1);
  testTree->Draw("0>>signalIntegral(1,0,1)", "weight*(className==\"Signal\")", "e goff");
  testTree->Draw("0>>totalIntegral(1,0,1)", "weight", "e goff");
  TEfficiency signalEfficiency( *((TH1F*)gDirectory->Get("signalIntegral")), *((TH1F*)gDirectory->Get("totalIntegral")));
  signalEfficiency.SetStatisticOption(TEfficiency::kFNormal); // kFNormal
  float f = signalEfficiency.GetEfficiency(1); // Signal fraction
  float df = signalEfficiency.GetEfficiencyErrorUp(1);
  float omf=1.-f; // One Minus f
  printf("Signal fraction = %.3e +/- %.2e \n",f,df);
  
  printf("Computing Mutual Information with the truth...\n");
  TH2F *pairwiseMutualInf = new TH2F("pairwiseMutualInf","Pairwise mutual information with the truth", nVars,0,nVars, nVars,0,nVars);
  std::vector<std::pair<unsigned, float> > varMutualInfs;
  float varMutualInfErrors[nVarsMax];
  varMutualInfs.reserve(nVars);
  for(unsigned h1=0; h1<nVars; h1++) {
    TH1F *pdf1d[2]; // array is truth[0|1]
    if(h1+1<nVars) {
      pdf1d[0] = (TH1F*)jointDistsAmongVars[h1][h1+1]->ProjectionX(varNames[h1]+"_background1D",0,-1,1,1,"e");
      pdf1d[1] = (TH1F*)jointDistsAmongVars[h1][h1+1]->ProjectionX(varNames[h1]+"_signal1D",0,-1,2,2,"e");
    } else {
      pdf1d[0] = (TH1F*)jointDistsAmongVars[h1-1][h1]->ProjectionY(varNames[h1]+"_background1D",0,-1,1,1,"e");
      pdf1d[1] = (TH1F*)jointDistsAmongVars[h1-1][h1]->ProjectionY(varNames[h1]+"_signal1D",0,-1,2,2,"e");
    }
    pdf1d[0]->Scale(1./pdf1d[0]->Integral());
    pdf1d[1]->Scale(1./pdf1d[1]->Integral());
    float mutualInfWithTruth=0, mutualError2=0;
    for(unsigned i=1; i<=(unsigned)pdf1d[0]->GetNbinsX(); i++) {
      float pSig=pdf1d[1]->GetBinContent(i), d_pSig=pdf1d[1]->GetBinError(i);
      float pBkg=pdf1d[0]->GetBinContent(i), d_pBkg=pdf1d[0]->GetBinError(i);
      float pTot=f*pSig+omf*pBkg;
      float d_pTot = sqrt(df*df*pSig*pSig + f*f*d_pSig*d_pSig + df*df*pBkg*pBkg + omf*omf*d_pBkg*d_pBkg);
      if(pTot==0) { if(debug) printf("\tsig+bkg=0 in bin #%d for variable \"%s\"\n", i, varNames[h1].Data()); continue;}
      float dMutual=
        f * pSig * TMath::Log2(pSig/pTot) +
        omf * pBkg * TMath::Log2(pBkg/pTot);
      if(dMutual!=dMutual) { if(debug) printf("\tincrement to mutual information is NaN for bin #%d, skipping\n", i); continue; }
      // To do: Check this error calculation
      mutualError2 += pow(d_pSig*(f*TMath::Log2(pSig/pTot) + f/ln2*(1-f*pSig/pTot) - omf/ln2*(-f*pBkg/pTot)), 2); // incremental error from pSig
      mutualError2 += pow(d_pBkg*(omf*TMath::Log2(pBkg/pTot) + omf/ln2*(1-omf*pBkg/pTot) - f/ln2*(omf*pSig/pTot)), 2); // incremental error from pSig
      if(debug) printf("\t%.1e * %.2f * log2(%.4f/%.4f) + %.1e * %.2f * log2(%.4f/%.4f) = %.3f\n", f, pSig, pSig,pTot, omf, pBkg,pBkg,pTot, dMutual);
      mutualInfWithTruth+=dMutual;

    }
    if(debug) printf("Mutual information with truth for variable  \"%s\": %.2e\n", varNames[h1].Data(), mutualInfWithTruth);
    varMutualInfs.push_back( std::pair<unsigned, float> (h1,mutualInfWithTruth));
    float mutualError=sqrt(mutualError2);
    varMutualInfErrors[h1]=mutualError;
    pairwiseMutualInf->SetBinContent(pairwiseMutualInf->FindBin(h1,h1), mutualInfWithTruth);
    pairwiseMutualInf->SetBinError(pairwiseMutualInf->FindBin(h1,h1), mutualError);
    
  }
  sort(varMutualInfs.begin(), varMutualInfs.end(), [=](std::pair<unsigned, float>& a, std::pair<unsigned, float>& b) {
    return a.second > b.second;
  });
  printf("Rank     I(T;A)       Error  Name(A)\n");
  printf("---------------------------------------\n");
  for(unsigned i=0; i<varMutualInfs.size(); i++)
    printf("%4d %10.2e  %10.2e  %s\n", i+1, varMutualInfs[i].second, varMutualInfErrors[varMutualInfs[i].first], varNames[varMutualInfs[i].first].Data());
  
  // pairwise correlations
  for(unsigned h1=0; h1<nVars; h1++) {
    if(h1+1<nVars) for(unsigned h2=h1+1; h2<nVars; h2++) {
      TH2F *pdf2d[2]; // array is truth[0|1]
      jointDistsAmongVars[h1][h2]->GetZaxis()->SetRange(1,1);
      pdf2d[0] = (TH2F*)jointDistsAmongVars[h1][h2]->Project3D("yxe");
      pdf2d[0]->SetName(varNames[h2]+"_vs_"+varNames[h1]+"_background2D");
      jointDistsAmongVars[h1][h2]->GetZaxis()->SetRange(2,2);
      pdf2d[1] = (TH2F*)jointDistsAmongVars[h1][h2]->Project3D("yxe");
      pdf2d[1]->SetName(varNames[h2]+"_vs_"+varNames[h1]+"_signal2D");
      jointDistsAmongVars[h1][h2]->GetZaxis()->SetRange(1,2);
      pdf2d[0]->Scale(1./pdf2d[0]->Integral());
      pdf2d[1]->Scale(1./pdf2d[1]->Integral());
      float mutualInfWithTruth=0, mutualError2=0;
      for(unsigned i=1; i<=(unsigned)pdf2d[0]->GetNbinsX(); i++) for(unsigned j=1; j<=(unsigned)pdf2d[0]->GetNbinsY(); j++) {
        int nb2d=pdf2d[0]->GetBin(i,j);
        float pSig=pdf2d[1]->GetBinContent(nb2d), d_pSig=pdf2d[1]->GetBinError(nb2d);
        float pBkg=pdf2d[0]->GetBinContent(nb2d), d_pBkg=pdf2d[0]->GetBinError(nb2d);
        float pTot=f*pSig+omf*pBkg;
        float d_pTot = sqrt(df*df*pSig*pSig + f*f*d_pSig*d_pSig + df*df*pBkg*pBkg + omf*omf*d_pBkg*d_pBkg);
        if(pTot==0) { if(debug) printf("\tsig+bkg=0 in bin (%d,%d) for variable pair (\"%s\",\"%s\")\n", i,j,varNames[h1].Data(),varNames[h2].Data()); continue;}
        float dMutual=
          f * pSig * TMath::Log2(pSig/pTot) +
          omf * pBkg * TMath::Log2(pBkg/pTot);
        if(dMutual!=dMutual) { if(debug) printf("\tincrement to mutual information is NaN for bin (%d,%d), skipping\n", i,j); continue; }
        mutualError2 += pow(d_pSig*(f*TMath::Log2(pSig/pTot) + f/ln2*(1-f*pSig/pTot) - omf/ln2*(-f*pBkg/pTot)), 2); // incremental error from pSig
        mutualError2 += pow(d_pBkg*(omf*TMath::Log2(pBkg/pTot) + omf/ln2*(1-omf*pBkg/pTot) - f/ln2*(omf*pSig/pTot)), 2); // incremental error from pSig
        if(debug) printf("\t%.1e * %.2f * log2(%.4f/%.4f) + %.1e * %.2f * log2(%.4f/%.4f) = %.3f\n", f, pSig, pSig,pTot, omf, pBkg,pBkg,pTot, dMutual);
        mutualInfWithTruth+=dMutual;
      }
      float mutualError=sqrt(mutualError2);
      if(debug) printf("Mutual information with truth for variable pair (\"%s\",\"%s\"): %.2e\n", varNames[h1].Data(), varNames[h2].Data(), mutualInfWithTruth);
      pairwiseMutualInf->SetBinContent(pairwiseMutualInf->FindBin(h1,h2), mutualInfWithTruth);
      pairwiseMutualInf->SetBinError(pairwiseMutualInf->FindBin(h1,h2), mutualError);
      pairwiseMutualInf->SetBinContent(pairwiseMutualInf->FindBin(h2,h1), mutualInfWithTruth);
      pairwiseMutualInf->SetBinError(pairwiseMutualInf->FindBin(h2,h1), mutualError);

    }
  }
  // Set axis labels
  pairwiseMutualInf->GetXaxis()->SetLabelSize(0.025);
  pairwiseMutualInf->GetYaxis()->SetLabelSize(0.025);
  for(unsigned h1=0; h1<nVars; h1++) {
    pairwiseMutualInf->GetXaxis()->SetBinLabel(h1+1, varTitles[h1]);
    pairwiseMutualInf->GetYaxis()->SetBinLabel(h1+1, varTitles[h1]);
  }

  // Start Drawing stuff
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlackBody);
  TCanvas *c_pairwiseMutualInf = new TCanvas("c_pairwiseMutualInf","c_pairwiseMutualInf",1000,800);
  c_pairwiseMutualInf->SetRightMargin(0.18);
  c_pairwiseMutualInf->SetBottomMargin(0.15);
  c_pairwiseMutualInf->SetLeftMargin(0.15);
  c_pairwiseMutualInf->SetLogz();
  pairwiseMutualInf->LabelsOption("v","x");
  pairwiseMutualInf->Draw("COLZ");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)pairwiseMutualInf->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.84);
  palette->SetX2NDC(0.9);
  palette->SetY1NDC(0.15);
  palette->SetY2NDC(0.9);
  gPad->Modified();
  gPad->Update();
  
}
/*
void bdt_toys(int nuisance=1) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  //const int MVAVarType = 3; const int nBinMVA = 15; Double_t xbins[nBinMVA+1] =  {-2, -1, 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.4}; TString addChan = "3";
  const int MVAVarType = 3; const int nBinMVA = 12; Double_t xbins[nBinMVA+1] =  {-2, -1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}; TString addChan = "3";
  //const int MVAVarType = 3; const int nBinMVA = 13; Double_t xbins[nBinMVA+1] =  {-2, -1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.575, 0.6}; TString addChan = "3";
  
  TFile *input_file_syst   = TFile::Open("zll3hinvll1j_sm_BDTsyst_13TeV2016.root", "READ");
  TFile *input_file_histos = TFile::Open("zll3hinvll1j_sm_input_13TeV2016.root", "READ");
  //TFile *input_file_syst   = TFile::Open("/home/dhsu/cms/cmssw/045/CMSSW_8_0_12/src/MitZHAnalysis/plots/BDT_multiClass_sm_met130_ntrees400_nodesize5pc_maxdepth4_shrinkage0p5_baggedfrac0p5_ncuts1000/zll3hinvll1j_sm_BDTsyst_13TeV2016.root", "READ");
  //TFile *input_file_histos = TFile::Open("/home/dhsu/cms/cmssw/045/CMSSW_8_0_12/src/MitZHAnalysis/plots/BDT_multiClass_sm_met130_ntrees400_nodesize5pc_maxdepth4_shrinkage0p5_baggedfrac0p5_ncuts1000/zll3hinvll1j_sm_input_13TeV2016.root", "READ");
  TH1F *nominal_shape_VVV        = (TH1F*) input_file_histos->Get("histo_VVV");
  TH1F *nominal_shape_WZ         = (TH1F*) input_file_histos->Get("histo_WZ");
  TH1F *nominal_shape_ZZ         = (TH1F*) input_file_histos->Get("histo_ZZ");
  TH1F *nominal_shape_ZH_hinv_sm = (TH1F*) input_file_histos->Get("histo_ZH_hinv_sm");
  TH1F *nominal_shape_ggZH_hinv  = (TH1F*) input_file_histos->Get("histo_ggZH_hinv");
  TH1F *stat_band_VVV        = (TH1F*) nominal_shape_VVV       ->Clone("stat_band_VVV"       ); stat_band_VVV       ->SetDirectory(0);  
  TH1F *stat_band_WZ         = (TH1F*) nominal_shape_WZ        ->Clone("stat_band_WZ"        ); stat_band_WZ        ->SetDirectory(0); 
  TH1F *stat_band_ZZ         = (TH1F*) nominal_shape_ZZ        ->Clone("stat_band_ZZ"        ); stat_band_ZZ        ->SetDirectory(0); 
  TH1F *stat_band_ZH_hinv_sm = (TH1F*) nominal_shape_ZH_hinv_sm->Clone("stat_band_ZH_hinv_sm"); stat_band_ZH_hinv_sm->SetDirectory(0); 
  TH1F *stat_band_ggZH_hinv  = (TH1F*) nominal_shape_ggZH_hinv ->Clone("stat_band_ggZH_hinv" ); stat_band_ggZH_hinv ->SetDirectory(0); 
  double envelope_max_VVV        = (1.1 * nominal_shape_VVV       ->GetMaximum()); 
  double envelope_max_WZ         = (1.1 * nominal_shape_WZ        ->GetMaximum()); 
  double envelope_max_ZZ         = (1.1 * nominal_shape_ZZ        ->GetMaximum()); 
  double envelope_max_ZH_hinv_sm = (1.1 * nominal_shape_ZH_hinv_sm->GetMaximum()); 
  double envelope_max_ggZH_hinv  = (1.1 * nominal_shape_ggZH_hinv ->GetMaximum()); 
  TH2F *toy_shapes_VVV       ;
  TH2F *toy_shapes_WZ        ;
  TH2F *toy_shapes_ZZ        ;
  TH2F *toy_shapes_ZH_hinv_sm;
  TH2F *toy_shapes_ggZH_hinv ;
  if(nuisance==1) {
    toy_shapes_VVV        = (TH2F*) input_file_syst->Get("histo_bdt_toys_electronScale_VVV");
    toy_shapes_WZ         = (TH2F*) input_file_syst->Get("histo_bdt_toys_electronScale_WZ");
    toy_shapes_ZZ         = (TH2F*) input_file_syst->Get("histo_bdt_toys_electronScale_ZZ");
    toy_shapes_ZH_hinv_sm = (TH2F*) input_file_syst->Get("histo_bdt_toys_electronScale_ZH_hinv_sm");
    toy_shapes_ggZH_hinv  = (TH2F*) input_file_syst->Get("histo_bdt_toys_electronScale_ggZH_hinv");
  } else if(nuisance==2) {
    toy_shapes_VVV        = (TH2F*) input_file_syst->Get("histo_bdt_toys_muonScale_VVV");
    toy_shapes_WZ         = (TH2F*) input_file_syst->Get("histo_bdt_toys_muonScale_WZ");
    toy_shapes_ZZ         = (TH2F*) input_file_syst->Get("histo_bdt_toys_muonScale_ZZ");
    toy_shapes_ZH_hinv_sm = (TH2F*) input_file_syst->Get("histo_bdt_toys_muonScale_ZH_hinv_sm");
    toy_shapes_ggZH_hinv  = (TH2F*) input_file_syst->Get("histo_bdt_toys_muonScale_ggZH_hinv");
  } else if(nuisance==3) {
    toy_shapes_VVV        = (TH2F*) input_file_syst->Get("histo_bdt_toys_METScale_VVV");
    toy_shapes_WZ         = (TH2F*) input_file_syst->Get("histo_bdt_toys_METScale_WZ");
    toy_shapes_ZZ         = (TH2F*) input_file_syst->Get("histo_bdt_toys_METScale_ZZ");
    toy_shapes_ZH_hinv_sm = (TH2F*) input_file_syst->Get("histo_bdt_toys_METScale_ZH_hinv_sm");
    toy_shapes_ggZH_hinv  = (TH2F*) input_file_syst->Get("histo_bdt_toys_METScale_ggZH_hinv");
  }
  TString the_nuisance_in_english, the_nuisance_basename;
  if     (nuisance==1) { the_nuisance_in_english = "electron scale"     ; the_nuisance_basename = "electron"     ;}
  else if(nuisance==2) { the_nuisance_in_english = "muon scale"         ; the_nuisance_basename = "muon"         ;}
  else if(nuisance==3) { the_nuisance_in_english = "E_{T}^{miss} scale" ; the_nuisance_basename = "MET"          ;}
  TH2F *toy_envelope_VVV        = new TH2F("toy_envelope_VVV"       , "toy_envelope_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
  TH2F *toy_envelope_WZ         = new TH2F("toy_envelope_WZ"        , "toy_envelope_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  TH2F *toy_envelope_ZZ         = new TH2F("toy_envelope_ZZ"        , "toy_envelope_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  TH2F *toy_envelope_ZH_hinv_sm = new TH2F("toy_envelope_ZH_hinv_sm", "toy_envelope_ZH_hinv_sm", nBinMVA, xbins, 1000 , 0, 4); 
  TH2F *toy_envelope_ggZH_hinv  = new TH2F("toy_envelope_ggZH_hinv" , "toy_envelope_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
  unsigned int num_toys=toy_shapes_VVV->GetNbinsY();
  
  TCanvas *c_binyields_VVV[nBinMVA], *c_binyields_WZ[nBinMVA], *c_binyields_ZZ[nBinMVA], *c_binyields_ZH_hinv_sm[nBinMVA], *c_binyields_ggZH_hinv[nBinMVA];
  TH1F *binyields_VVV[nBinMVA], *binyields_WZ[nBinMVA], *binyields_ZZ[nBinMVA], *binyields_ZH_hinv_sm[nBinMVA], *binyields_ggZH_hinv[nBinMVA];
  for(int nb=1; nb<=nBinMVA; nb++) {
    stat_band_VVV        -> SetBinContent(nb, 1); 
    stat_band_WZ         -> SetBinContent(nb, 1); 
    stat_band_ZZ         -> SetBinContent(nb, 1); 
    stat_band_ZH_hinv_sm -> SetBinContent(nb, 1); 
    stat_band_ggZH_hinv  -> SetBinContent(nb, 1); 
    stat_band_VVV        -> SetBinError(nb, nominal_shape_VVV       ->GetBinError(nb) / nominal_shape_VVV       ->GetBinContent(nb));
    stat_band_WZ         -> SetBinError(nb, nominal_shape_WZ        ->GetBinError(nb) / nominal_shape_WZ        ->GetBinContent(nb));
    stat_band_ZZ         -> SetBinError(nb, nominal_shape_ZZ        ->GetBinError(nb) / nominal_shape_ZZ        ->GetBinContent(nb));
    stat_band_ZH_hinv_sm -> SetBinError(nb, nominal_shape_ZH_hinv_sm->GetBinError(nb) / nominal_shape_ZH_hinv_sm->GetBinContent(nb));
    stat_band_ggZH_hinv  -> SetBinError(nb, nominal_shape_ggZH_hinv ->GetBinError(nb) / nominal_shape_ggZH_hinv ->GetBinContent(nb));
    for(unsigned int i_toy=1; i_toy<=num_toys; i_toy++) {
      if(nominal_shape_VVV       ->GetBinContent(nb)!=0) toy_envelope_VVV       ->Fill( xbins[nb-1], toy_shapes_VVV       ->GetBinContent( toy_shapes_VVV       ->GetBin(nb, i_toy)) / nominal_shape_VVV       ->GetBinContent(nb));  
      if(nominal_shape_WZ        ->GetBinContent(nb)!=0) toy_envelope_WZ        ->Fill( xbins[nb-1], toy_shapes_WZ        ->GetBinContent( toy_shapes_WZ        ->GetBin(nb, i_toy)) / nominal_shape_WZ        ->GetBinContent(nb));  
      if(nominal_shape_ZZ        ->GetBinContent(nb)!=0) toy_envelope_ZZ        ->Fill( xbins[nb-1], toy_shapes_ZZ        ->GetBinContent( toy_shapes_ZZ        ->GetBin(nb, i_toy)) / nominal_shape_ZZ        ->GetBinContent(nb));  
      if(nominal_shape_ZH_hinv_sm->GetBinContent(nb)!=0) toy_envelope_ZH_hinv_sm->Fill( xbins[nb-1], toy_shapes_ZH_hinv_sm->GetBinContent( toy_shapes_ZH_hinv_sm->GetBin(nb, i_toy)) / nominal_shape_ZH_hinv_sm->GetBinContent(nb));  
      if(nominal_shape_ggZH_hinv ->GetBinContent(nb)!=0) toy_envelope_ggZH_hinv ->Fill( xbins[nb-1], toy_shapes_ggZH_hinv ->GetBinContent( toy_shapes_ggZH_hinv ->GetBin(nb, i_toy)) / nominal_shape_ggZH_hinv ->GetBinContent(nb));  
  
    }
    if(nb>2) {
      // Make toy distribution in individual shape bins for VVV
      binyields_VVV[nb-1] = (TH1F*) toy_envelope_VVV->ProjectionY(Form("binyields_VVV_%d", nb),nb,nb);
      c_binyields_VVV[nb-1] = new TCanvas(Form("c_binyields_VVV_%d", nb), Form("c_binyields_VVV_%d", nb));
      binyields_VVV[nb-1]->SetTitle(Form("Distribution of toy yields for %s in VVV process, BDT bin #%d", the_nuisance_in_english.Data(), nb-2));
      binyields_VVV[nb-1]->Draw("HIST");
      binyields_VVV[nb-1]->GetXaxis()->SetRangeUser(0.2,1.8);
      binyields_VVV[nb-1]->GetXaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
      c_binyields_VVV[nb-1]->Update();
      c_binyields_VVV[nb-1]->Print(Form("syst_BDT_VVV_toyhisto_bin%d_%s.pdf", nb-2, the_nuisance_basename.Data()));
      delete c_binyields_VVV[nb-1];
      // Make toy distribution in individual shape bins for WZ
      binyields_WZ[nb-1] = (TH1F*) toy_envelope_WZ->ProjectionY(Form("binyields_WZ_%d", nb),nb,nb);
      c_binyields_WZ[nb-1] = new TCanvas(Form("c_binyields_WZ_%d", nb), Form("c_binyields_WZ_%d", nb));
      binyields_WZ[nb-1]->SetTitle(Form("Distribution of toy yields for %s in WZ process, BDT bin #%d", the_nuisance_in_english.Data(), nb-2));
      binyields_WZ[nb-1]->Draw("HIST");
      binyields_WZ[nb-1]->GetXaxis()->SetRangeUser(0.8,1.199);
      binyields_WZ[nb-1]->GetXaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
      c_binyields_WZ[nb-1]->Update();
      c_binyields_WZ[nb-1]->Print(Form("syst_BDT_WZ_toyhisto_bin%d_%s.pdf", nb-2, the_nuisance_basename.Data()));
      delete c_binyields_WZ[nb-1];
      // Make toy distribution in individual shape bins for ZZ
      binyields_ZZ[nb-1] = (TH1F*) toy_envelope_ZZ->ProjectionY(Form("binyields_ZZ_%d", nb),nb,nb);
      c_binyields_ZZ[nb-1] = new TCanvas(Form("c_binyields_ZZ_%d", nb), Form("c_binyields_ZZ_%d", nb));
      binyields_ZZ[nb-1]->SetTitle(Form("Distribution of toy yields for %s in ZZ process, BDT bin #%d", the_nuisance_in_english.Data(), nb-2));
      binyields_ZZ[nb-1]->Draw("HIST");
      binyields_ZZ[nb-1]->GetXaxis()->SetRangeUser(0.8,1.199);
      binyields_ZZ[nb-1]->GetXaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
      c_binyields_ZZ[nb-1]->Update();
      c_binyields_ZZ[nb-1]->Print(Form("syst_BDT_ZZ_toyhisto_bin%d_%s.pdf", nb-2, the_nuisance_basename.Data()));
      delete c_binyields_ZZ[nb-1];
      // Make toy distribution in individual shape bins for ZH_hinv_sm
      binyields_ZH_hinv_sm[nb-1] = (TH1F*) toy_envelope_ZH_hinv_sm->ProjectionY(Form("binyields_ZH_hinv_sm_%d", nb),nb,nb);
      c_binyields_ZH_hinv_sm[nb-1] = new TCanvas(Form("c_binyields_ZH_hinv_sm_%d", nb), Form("c_binyields_ZH_hinv_sm_%d", nb));
      binyields_ZH_hinv_sm[nb-1]->SetTitle(Form("Distribution of toy yields for %s in #font[12]{qq}#rightarrowZH process, BDT bin #%d", the_nuisance_in_english.Data(), nb-2));
      binyields_ZH_hinv_sm[nb-1]->Draw("HIST");
      binyields_ZH_hinv_sm[nb-1]->GetXaxis()->SetRangeUser(0.8,1.199);
      binyields_ZH_hinv_sm[nb-1]->GetXaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
      c_binyields_ZH_hinv_sm[nb-1]->Update();
      c_binyields_ZH_hinv_sm[nb-1]->Print(Form("syst_BDT_ZH_hinv_sm_toyhisto_bin%d_%s.pdf", nb-2, the_nuisance_basename.Data()));
      delete c_binyields_ZH_hinv_sm[nb-1];
      // Make toy distribution in individual shape bins for ggZH_hinv
      binyields_ggZH_hinv[nb-1] = (TH1F*) toy_envelope_ggZH_hinv->ProjectionY(Form("binyields_ggZH_hinv_%d", nb),nb,nb);
      c_binyields_ggZH_hinv[nb-1] = new TCanvas(Form("c_binyields_ggZH_hinv_%d", nb), Form("c_binyields_ggZH_hinv_%d", nb));
      binyields_ggZH_hinv[nb-1]->SetTitle(Form("Distribution of toy yields for %s in #font[12]{gg}#rightarrowZH process, BDT bin #%d", the_nuisance_in_english.Data(), nb-2));
      binyields_ggZH_hinv[nb-1]->Draw("HIST");
      binyields_ggZH_hinv[nb-1]->GetXaxis()->SetRangeUser(0.8,1.199);
      binyields_ggZH_hinv[nb-1]->GetXaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
      c_binyields_ggZH_hinv[nb-1]->Update();
      c_binyields_ggZH_hinv[nb-1]->Print(Form("syst_BDT_ggZH_hinv_toyhisto_bin%d_%s.pdf", nb-2, the_nuisance_basename.Data()));
      delete c_binyields_ggZH_hinv[nb-1];
    }
  }
  TCanvas *c_envelope_VVV = new TCanvas("c_envelope_VVV","envelope VVV");
  toy_envelope_VVV->Rebin2D(1,5);
  toy_envelope_VVV->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  toy_envelope_VVV->GetYaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
  toy_envelope_VVV->GetYaxis()->SetRangeUser(0.5,1.5);
  toy_envelope_VVV->GetXaxis()->SetRangeUser(xbins[2],xbins[nBinMVA]-.01);
  toy_envelope_VVV->SetTitle(Form("Toy envelope / nominal bin yields (%s, VVV processes)", the_nuisance_in_english.Data()));
  toy_envelope_VVV->Draw("COLZ");
  stat_band_VVV->SetFillStyle(3004);
  stat_band_VVV->SetFillColor(1);
  stat_band_VVV->Draw("E2 SAME");
  c_envelope_VVV->Print(Form("syst_BDT_VVV_toyenvelope_%s.pdf", the_nuisance_basename.Data()));
  TCanvas *c_envelope_WZ = new TCanvas("c_envelope_WZ","envelope WZ");
  toy_envelope_WZ->Rebin2D(1,2);
  toy_envelope_WZ->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  toy_envelope_WZ->GetYaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
  toy_envelope_WZ->GetXaxis()->SetRangeUser(xbins[2],xbins[nBinMVA]-.01);
  toy_envelope_WZ->GetYaxis()->SetRangeUser(0.8,1.2);
  toy_envelope_WZ->SetTitle(Form("Toy envelope / nominal bin yields (%s, WZ processes)", the_nuisance_in_english.Data()));
  toy_envelope_WZ->Draw("COLZ");
  stat_band_WZ->SetFillStyle(3004);
  stat_band_WZ->SetFillColor(1);
  stat_band_WZ->Draw("E2 SAME");
  c_envelope_WZ->Print(Form("syst_BDT_WZ_toyenvelope_%s.pdf", the_nuisance_basename.Data()));
  TCanvas *c_envelope_ZZ = new TCanvas("c_envelope_ZZ","envelope ZZ");
  toy_envelope_ZZ->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  toy_envelope_ZZ->GetYaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
  toy_envelope_ZZ->GetXaxis()->SetRangeUser(xbins[2],xbins[nBinMVA]-.01);
  toy_envelope_ZZ->GetYaxis()->SetRangeUser(0.9,1.1);
  toy_envelope_ZZ->SetTitle(Form("Toy envelope / nominal bin yields (%s, ZZ processes)", the_nuisance_in_english.Data()));
  toy_envelope_ZZ->Draw("COLZ");
  stat_band_ZZ->SetFillStyle(3004);
  stat_band_ZZ->SetFillColor(1);
  stat_band_ZZ->Draw("E2 SAME");
  c_envelope_ZZ->Print(Form("syst_BDT_ZZ_toyenvelope_%s.pdf", the_nuisance_basename.Data()));
  TCanvas *c_envelope_ZH_hinv_sm = new TCanvas("c_envelope_ZH_hinv_sm","envelope ZH_hinv_sm");
  toy_envelope_ZH_hinv_sm->Rebin2D(1,2);
  toy_envelope_ZH_hinv_sm->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  toy_envelope_ZH_hinv_sm->GetYaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
  toy_envelope_ZH_hinv_sm->GetXaxis()->SetRangeUser(xbins[2],xbins[nBinMVA]-.01);
  toy_envelope_ZH_hinv_sm->GetYaxis()->SetRangeUser(0.8,1.2);
  toy_envelope_ZH_hinv_sm->SetTitle(Form("Toy envelope / nominal bin yields (%s, #font[12]{qq}#rightarrowZH process)",the_nuisance_in_english.Data()));
  toy_envelope_ZH_hinv_sm->Draw("COLZ");
  stat_band_ZH_hinv_sm->SetFillStyle(3004);
  stat_band_ZH_hinv_sm->SetFillColor(1);
  stat_band_ZH_hinv_sm->Draw("E2 SAME");
  c_envelope_ZH_hinv_sm->Print(Form("syst_BDT_ZH_hinv_sm_toyenvelope_%s.pdf", the_nuisance_basename.Data()));
  TCanvas *c_envelope_ggZH_hinv = new TCanvas("c_envelope_ggZH_hinv","envelope ggZH_hinv");
  toy_envelope_ggZH_hinv->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  toy_envelope_ggZH_hinv->GetYaxis()->SetTitle("Toy shape bin yield / nominal bin yield");
  toy_envelope_ggZH_hinv->GetXaxis()->SetRangeUser(xbins[2],xbins[nBinMVA]-.01);
  toy_envelope_ggZH_hinv->GetYaxis()->SetRangeUser(0.9,1.1);
  toy_envelope_ggZH_hinv->SetTitle(Form("Toy envelope / nominal bin yields (%s, #font[12]{gg}#rightarrowZH process)",the_nuisance_in_english.Data()));
  toy_envelope_ggZH_hinv->Draw("COLZ");
  stat_band_ggZH_hinv->SetFillStyle(3004);
  stat_band_ggZH_hinv->SetFillColor(1);
  stat_band_ggZH_hinv->Draw("E2 SAME");
  c_envelope_ggZH_hinv->Print(Form("syst_BDT_ggZH_hinv_toyenvelope_%s.pdf", the_nuisance_basename.Data()));
   

  TH1F *bdt_VVV_nuisanceUp_ratio   = (TH1F*) nominal_shape_VVV->Clone("bdt_VVV_nuisanceUp_ratio");    bdt_VVV_nuisanceUp_ratio  ->SetDirectory(0); bdt_VVV_nuisanceUp_ratio  ->Scale(0);
  TH1F *bdt_VVV_nuisanceDown_ratio = (TH1F*) nominal_shape_VVV->Clone("bdt_VVV_nuisanceDown_ratio");  bdt_VVV_nuisanceDown_ratio->SetDirectory(0); bdt_VVV_nuisanceDown_ratio->Scale(0);
  TH1F *bdt_WZ_nuisanceUp_ratio   = (TH1F*) nominal_shape_WZ->Clone("bdt_WZ_nuisanceUp_ratio");    bdt_WZ_nuisanceUp_ratio  ->SetDirectory(0); bdt_WZ_nuisanceUp_ratio  ->Scale(0);
  TH1F *bdt_WZ_nuisanceDown_ratio = (TH1F*) nominal_shape_WZ->Clone("bdt_WZ_nuisanceDown_ratio");  bdt_WZ_nuisanceDown_ratio->SetDirectory(0); bdt_WZ_nuisanceDown_ratio->Scale(0);
  TH1F *bdt_ZZ_nuisanceUp_ratio   = (TH1F*) nominal_shape_ZZ->Clone("bdt_ZZ_nuisanceUp_ratio");    bdt_ZZ_nuisanceUp_ratio  ->SetDirectory(0); bdt_ZZ_nuisanceUp_ratio  ->Scale(0);
  TH1F *bdt_ZZ_nuisanceDown_ratio = (TH1F*) nominal_shape_ZZ->Clone("bdt_ZZ_nuisanceDown_ratio");  bdt_ZZ_nuisanceDown_ratio->SetDirectory(0); bdt_ZZ_nuisanceDown_ratio->Scale(0);
  TH1F *bdt_ZH_hinv_sm_nuisanceUp_ratio   = (TH1F*) nominal_shape_ZH_hinv_sm->Clone("bdt_ZH_hinv_sm_nuisanceUp_ratio");    bdt_ZH_hinv_sm_nuisanceUp_ratio  ->SetDirectory(0); bdt_ZH_hinv_sm_nuisanceUp_ratio  ->Scale(0);
  TH1F *bdt_ZH_hinv_sm_nuisanceDown_ratio = (TH1F*) nominal_shape_ZH_hinv_sm->Clone("bdt_ZH_hinv_sm_nuisanceDown_ratio");  bdt_ZH_hinv_sm_nuisanceDown_ratio->SetDirectory(0); bdt_ZH_hinv_sm_nuisanceDown_ratio->Scale(0);
  TH1F *bdt_ggZH_hinv_nuisanceUp_ratio   = (TH1F*) nominal_shape_ggZH_hinv->Clone("bdt_ggZH_hinv_nuisanceUp_ratio");    bdt_ggZH_hinv_nuisanceUp_ratio  ->SetDirectory(0); bdt_ggZH_hinv_nuisanceUp_ratio  ->Scale(0);
  TH1F *bdt_ggZH_hinv_nuisanceDown_ratio = (TH1F*) nominal_shape_ggZH_hinv->Clone("bdt_ggZH_hinv_nuisanceDown_ratio");  bdt_ggZH_hinv_nuisanceDown_ratio->SetDirectory(0); bdt_ggZH_hinv_nuisanceDown_ratio->Scale(0);

  for(int nb=3; nb<=nBinMVA; nb++) {
    double toy_mean, toy_rms, yield_error;
    double quantileProbs[3]={0.159,0.5,0.841};
    double theQuantiles[3];
    // Compute the syst. for VVV
    binyields_VVV[nb-1]->GetQuantiles(3,theQuantiles,quantileProbs);
    bdt_VVV_nuisanceUp_ratio  ->SetBinContent(nb, theQuantiles[2]);
    bdt_VVV_nuisanceDown_ratio->SetBinContent(nb, theQuantiles[0]);
    printf("bin %d quantiles for VVV: %f, %f, %f\n", nb-2, theQuantiles[0],theQuantiles[1], theQuantiles[2]);
    // Compute the syst. for WZ
    binyields_WZ[nb-1]->GetQuantiles(3,theQuantiles,quantileProbs);
    bdt_WZ_nuisanceUp_ratio  ->SetBinContent(nb, theQuantiles[2]);
    bdt_WZ_nuisanceDown_ratio->SetBinContent(nb, theQuantiles[0]);
    printf("bin %d quantiles for WZ: %f, %f, %f\n", nb-2, theQuantiles[0],theQuantiles[1], theQuantiles[2]);
    // Compute the syst. for ZZ
    binyields_ZZ[nb-1]->GetQuantiles(3,theQuantiles,quantileProbs);
    bdt_ZZ_nuisanceUp_ratio  ->SetBinContent(nb, theQuantiles[2]);
    bdt_ZZ_nuisanceDown_ratio->SetBinContent(nb, theQuantiles[0]);
    printf("bin %d quantiles for ZZ: %f, %f, %f\n", nb-2, theQuantiles[0],theQuantiles[1], theQuantiles[2]);
    // Compute the syst. for ZH_hinv_sm
    binyields_ZH_hinv_sm[nb-1]->GetQuantiles(3,theQuantiles,quantileProbs);
    bdt_ZH_hinv_sm_nuisanceUp_ratio  ->SetBinContent(nb, theQuantiles[2]);
    bdt_ZH_hinv_sm_nuisanceDown_ratio->SetBinContent(nb, theQuantiles[0]);
    printf("bin %d quantiles for qqZH: %f, %f, %f\n", nb-2, theQuantiles[0],theQuantiles[1], theQuantiles[2]);
    // Compute the syst. for ggZH_hinv
    binyields_ggZH_hinv[nb-1]->GetQuantiles(3,theQuantiles,quantileProbs);
    bdt_ggZH_hinv_nuisanceUp_ratio  ->SetBinContent(nb, theQuantiles[2]);
    bdt_ggZH_hinv_nuisanceDown_ratio->SetBinContent(nb, theQuantiles[0]);
    printf("bin %d quantiles for ggZH: %f, %f, %f\n", nb-2, theQuantiles[0],theQuantiles[1], theQuantiles[2]);
  }


  TLine *oneline = new TLine;
  oneline->SetLineColor(13);
  oneline->SetLineWidth(2);
  
  // Plot the systematic bands
  // VVV
  TCanvas *c_VVV_nuisance = new TCanvas("c_VVV_nuisance","c_VVV_nuisance");
  bdt_VVV_nuisanceUp_ratio  ->SetLineWidth( 3 );
  bdt_VVV_nuisanceDown_ratio->SetLineWidth( 3 );
  bdt_VVV_nuisanceUp_ratio  ->SetLineColor( kViolet-1 );
  bdt_VVV_nuisanceDown_ratio->SetLineColor( kBlue+1   );
  bdt_VVV_nuisanceUp_ratio->SetMinimum(0.8);
  bdt_VVV_nuisanceUp_ratio->SetMaximum(1.2);
  bdt_VVV_nuisanceUp_ratio->SetTitle(Form("BDT syst. unc. from %s toy study in VVV", the_nuisance_in_english.Data()));
  bdt_VVV_nuisanceUp_ratio->Draw("HIST");
  bdt_VVV_nuisanceUp_ratio->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  bdt_VVV_nuisanceUp_ratio->GetYaxis()->SetTitle("Ratio to nominal");
  bdt_VVV_nuisanceUp_ratio->GetXaxis()->SetRangeUser(0.2,xbins[nBinMVA]-.0001);
  
  bdt_VVV_nuisanceDown_ratio->Draw("HIST SAME");
  stat_band_VVV->Draw("E2 SAME");
  oneline->DrawLine(0.2,1,xbins[nBinMVA]-.0001,1);
  TLegend *l_bdt_VVV_nuisance = new TLegend(.35,.15,.55,.3);
  l_bdt_VVV_nuisance->SetFillColor(0);
  l_bdt_VVV_nuisance->AddEntry(bdt_VVV_nuisanceUp_ratio, Form("%c%s up",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_VVV_nuisance->AddEntry(bdt_VVV_nuisanceDown_ratio, Form("%c%s down",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_VVV_nuisance->AddEntry(stat_band_VVV, "Stat. Unc. on Bins", "f");
  l_bdt_VVV_nuisance->Draw("SAME");
  c_VVV_nuisance->Print(Form("syst_BDT_VVV_toys_%s.pdf",the_nuisance_basename.Data()));

  // WZ
  TCanvas *c_WZ_nuisance = new TCanvas("c_WZ_nuisance","c_WZ_nuisance");
  bdt_WZ_nuisanceUp_ratio  ->SetLineWidth( 3 );
  bdt_WZ_nuisanceDown_ratio->SetLineWidth( 3 );
  bdt_WZ_nuisanceUp_ratio  ->SetLineColor( kViolet-1 );
  bdt_WZ_nuisanceDown_ratio->SetLineColor( kBlue+1   );
  bdt_WZ_nuisanceUp_ratio->SetMinimum(0.9);
  bdt_WZ_nuisanceUp_ratio->SetMaximum(1.1);
  bdt_WZ_nuisanceUp_ratio->SetTitle(Form("BDT syst. unc. from %s toy study in WZ", the_nuisance_in_english.Data()));
  bdt_WZ_nuisanceUp_ratio->Draw("HIST");
  bdt_WZ_nuisanceUp_ratio->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  bdt_WZ_nuisanceUp_ratio->GetYaxis()->SetTitle("Ratio to nominal");
  bdt_WZ_nuisanceUp_ratio->GetXaxis()->SetRangeUser(0.2,xbins[nBinMVA]-.0001);
  
  bdt_WZ_nuisanceDown_ratio->Draw("HIST SAME");
  stat_band_WZ->Draw("E2 SAME");
  oneline->DrawLine(0.2,1,xbins[nBinMVA]-.0001,1);
  TLegend *l_bdt_WZ_nuisance = new TLegend(.35,.15,.55,.3);
  l_bdt_WZ_nuisance->SetFillColor(0);
  l_bdt_WZ_nuisance->AddEntry(bdt_WZ_nuisanceUp_ratio, Form("%c%s up",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_WZ_nuisance->AddEntry(bdt_WZ_nuisanceDown_ratio, Form("%c%s down",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_WZ_nuisance->AddEntry(stat_band_WZ, "Stat. Unc. on Bins", "f");
  l_bdt_WZ_nuisance->Draw("SAME");
  c_WZ_nuisance->Print(Form("syst_BDT_WZ_toys_%s.pdf",the_nuisance_basename.Data()));

  // ZZ
  TCanvas *c_ZZ_nuisance = new TCanvas("c_ZZ_nuisance","c_ZZ_nuisance");
  bdt_ZZ_nuisanceUp_ratio  ->SetLineWidth( 3 );
  bdt_ZZ_nuisanceDown_ratio->SetLineWidth( 3 );
  bdt_ZZ_nuisanceUp_ratio  ->SetLineColor( kViolet-1 );
  bdt_ZZ_nuisanceDown_ratio->SetLineColor( kBlue+1   );
  bdt_ZZ_nuisanceUp_ratio->SetMinimum(0.95);
  bdt_ZZ_nuisanceUp_ratio->SetMaximum(1.05);
  bdt_ZZ_nuisanceUp_ratio->SetTitle(Form("BDT syst. unc. from %s toy study in ZZ", the_nuisance_in_english.Data()));
  bdt_ZZ_nuisanceUp_ratio->Draw("HIST");
  bdt_ZZ_nuisanceUp_ratio->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  bdt_ZZ_nuisanceUp_ratio->GetYaxis()->SetTitle("Ratio to nominal");
  bdt_ZZ_nuisanceUp_ratio->GetXaxis()->SetRangeUser(0.2,xbins[nBinMVA]-.0001);
  
  bdt_ZZ_nuisanceDown_ratio->Draw("HIST SAME");
  stat_band_ZZ->Draw("E2 SAME");
  oneline->DrawLine(0.2,1,xbins[nBinMVA]-.0001,1);
  TLegend *l_bdt_ZZ_nuisance = new TLegend(.35,.15,.55,.3);
  l_bdt_ZZ_nuisance->SetFillColor(0);
  l_bdt_ZZ_nuisance->AddEntry(bdt_ZZ_nuisanceUp_ratio, Form("%c%s up",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_ZZ_nuisance->AddEntry(bdt_ZZ_nuisanceDown_ratio, Form("%c%s down",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_ZZ_nuisance->AddEntry(stat_band_ZZ, "Stat. Unc. on Bins", "f");
  l_bdt_ZZ_nuisance->Draw("SAME");
  c_ZZ_nuisance->Print(Form("syst_BDT_ZZ_toys_%s.pdf",the_nuisance_basename.Data()));

  // ZH_hinv_sm
  TCanvas *c_ZH_hinv_sm_nuisance = new TCanvas("c_ZH_hinv_sm_nuisance","c_ZH_hinv_sm_nuisance");
  bdt_ZH_hinv_sm_nuisanceUp_ratio  ->SetLineWidth( 3 );
  bdt_ZH_hinv_sm_nuisanceDown_ratio->SetLineWidth( 3 );
  bdt_ZH_hinv_sm_nuisanceUp_ratio  ->SetLineColor( kViolet-1 );
  bdt_ZH_hinv_sm_nuisanceDown_ratio->SetLineColor( kBlue+1   );
  bdt_ZH_hinv_sm_nuisanceUp_ratio->SetMinimum(0.9);
  bdt_ZH_hinv_sm_nuisanceUp_ratio->SetMaximum(1.1);
  bdt_ZH_hinv_sm_nuisanceUp_ratio->SetTitle(Form("BDT syst. unc. from %s toy study in #font[12]{qq}#rightarrowZH(125)", the_nuisance_in_english.Data()));
  bdt_ZH_hinv_sm_nuisanceUp_ratio->Draw("HIST");
  bdt_ZH_hinv_sm_nuisanceUp_ratio->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  bdt_ZH_hinv_sm_nuisanceUp_ratio->GetYaxis()->SetTitle("Ratio to nominal");
  bdt_ZH_hinv_sm_nuisanceUp_ratio->GetXaxis()->SetRangeUser(0.2,xbins[nBinMVA]-.0001);
  
  bdt_ZH_hinv_sm_nuisanceDown_ratio->Draw("HIST SAME");
  stat_band_ZH_hinv_sm->Draw("E2 SAME");
  oneline->DrawLine(0.2,1,xbins[nBinMVA]-.0001,1);
  TLegend *l_bdt_ZH_hinv_sm_nuisance = new TLegend(.35,.15,.55,.3);
  l_bdt_ZH_hinv_sm_nuisance->SetFillColor(0);
  l_bdt_ZH_hinv_sm_nuisance->AddEntry(bdt_ZH_hinv_sm_nuisanceUp_ratio, Form("%c%s up",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_ZH_hinv_sm_nuisance->AddEntry(bdt_ZH_hinv_sm_nuisanceDown_ratio, Form("%c%s down",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_ZH_hinv_sm_nuisance->AddEntry(stat_band_ZH_hinv_sm, "Stat. Unc. on Bins", "f");
  l_bdt_ZH_hinv_sm_nuisance->Draw("SAME");
  c_ZH_hinv_sm_nuisance->Print(Form("syst_BDT_ZH_hinv_sm_toys_%s.pdf",the_nuisance_basename.Data()));

  // ggZH_hinv
  TCanvas *c_ggZH_hinv_nuisance = new TCanvas("c_ggZH_hinv_nuisance","c_ggZH_hinv_nuisance");
  bdt_ggZH_hinv_nuisanceUp_ratio  ->SetLineWidth( 3 );
  bdt_ggZH_hinv_nuisanceDown_ratio->SetLineWidth( 3 );
  bdt_ggZH_hinv_nuisanceUp_ratio  ->SetLineColor( kViolet-1 );
  bdt_ggZH_hinv_nuisanceDown_ratio->SetLineColor( kBlue+1   );
  bdt_ggZH_hinv_nuisanceUp_ratio->SetMinimum(0.95);
  bdt_ggZH_hinv_nuisanceUp_ratio->SetMaximum(1.05);
  bdt_ggZH_hinv_nuisanceUp_ratio->SetTitle(Form("BDT syst. unc. from %s toy study in #font[12]{gg}#rightarrowZH(125)", the_nuisance_in_english.Data()));
  bdt_ggZH_hinv_nuisanceUp_ratio->Draw("HIST");
  bdt_ggZH_hinv_nuisanceUp_ratio->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  bdt_ggZH_hinv_nuisanceUp_ratio->GetYaxis()->SetTitle("Ratio to nominal");
  bdt_ggZH_hinv_nuisanceUp_ratio->GetXaxis()->SetRangeUser(0.2,xbins[nBinMVA]-.0001);
  
  bdt_ggZH_hinv_nuisanceDown_ratio->Draw("HIST SAME");
  stat_band_ggZH_hinv->Draw("E2 SAME");
  oneline->DrawLine(0.2,1,xbins[nBinMVA]-.0001,1);
  TLegend *l_bdt_ggZH_hinv_nuisance = new TLegend(.35,.15,.55,.3);
  l_bdt_ggZH_hinv_nuisance->SetFillColor(0);
  l_bdt_ggZH_hinv_nuisance->AddEntry(bdt_ggZH_hinv_nuisanceUp_ratio, Form("%c%s up",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_ggZH_hinv_nuisance->AddEntry(bdt_ggZH_hinv_nuisanceDown_ratio, Form("%c%s down",toupper(the_nuisance_in_english[0]),(string(the_nuisance_in_english.Data()).substr(1)).c_str()), "l");
  l_bdt_ggZH_hinv_nuisance->AddEntry(stat_band_ggZH_hinv, "Stat. Unc. on Bins", "f");
  l_bdt_ggZH_hinv_nuisance->Draw("SAME");
  c_ggZH_hinv_nuisance->Print(Form("syst_BDT_ggZH_hinv_toys_%s.pdf",the_nuisance_basename.Data()));

}

void bdt_ewk() {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  //const int MVAVarType = 3; const int nBinMVA = 15; Double_t xbins[nBinMVA+1] =  {-2, -1, 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.4}; TString addChan = "3";
  const int MVAVarType = 3; const int nBinMVA = 12; Double_t xbins[nBinMVA+1] =  {-2, -1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}; TString addChan = "3";
  
  TFile *input_file_histos = TFile::Open("/home/dhsu/cms/cmssw/045/CMSSW_8_0_12/src/MitZHAnalysis/plots/BDT_multiClass_sm_met130_ntrees400_nodesize5pc_maxdepth4_shrinkage0p5_baggedfrac0p5_ncuts1000/zll3hinvll1j_sm_input_13TeV2016.root", "READ");
  TH1F *nominal_shape_WZ         = (TH1F*) input_file_histos->Get("histo_WZ");
  TH1F *nominal_shape_ZZ         = (TH1F*) input_file_histos->Get("histo_ZZ");
  TH1F *stat_band_WZ         = (TH1F*) nominal_shape_WZ        ->Clone("stat_band_WZ"        ); stat_band_WZ        ->SetDirectory(0); 
  TH1F *stat_band_ZZ         = (TH1F*) nominal_shape_ZZ        ->Clone("stat_band_ZZ"        ); stat_band_ZZ        ->SetDirectory(0); 
  TH1F *shape_WZ_EWKCorrUp   = (TH1F*) input_file_histos->Get("histo_WZ_EWKCorrUp");    shape_WZ_EWKCorrUp   ->SetDirectory(0);
  TH1F *shape_WZ_EWKCorrDown = (TH1F*) input_file_histos->Get("histo_WZ_EWKCorrDown");  shape_WZ_EWKCorrDown ->SetDirectory(0);
  TH1F *shape_ZZ_EWKCorrUp   = (TH1F*) input_file_histos->Get("histo_ZZ_EWKCorrUp");    shape_ZZ_EWKCorrUp   ->SetDirectory(0);
  TH1F *shape_ZZ_EWKCorrDown = (TH1F*) input_file_histos->Get("histo_ZZ_EWKCorrDown");  shape_ZZ_EWKCorrDown ->SetDirectory(0);
  TH1F *shape_signalRegion_ZZWZ_EWKCorrUp   = (TH1F*) shape_WZ_EWKCorrUp->Clone("shape_signalRegion_ZZWZ_EWKCorrUp");    shape_signalRegion_ZZWZ_EWKCorrUp   ->SetDirectory(0);
  TH1F *shape_signalRegion_ZZWZ_EWKCorrDown = (TH1F*) shape_WZ_EWKCorrUp->Clone("shape_signalRegion_ZZWZ_EWKCorrDown");  shape_signalRegion_ZZWZ_EWKCorrDown ->SetDirectory(0);
  TH1F *shape_VVratio_ZZWZ_EWKCorrUp        = (TH1F*) shape_WZ_EWKCorrUp->Clone("shape_VVratio_ZZWZ_EWKCorrUp");         shape_VVratio_ZZWZ_EWKCorrUp        ->SetDirectory(0);
  TH1F *shape_VVratio_ZZWZ_EWKCorrDown      = (TH1F*) shape_WZ_EWKCorrUp->Clone("shape_VVratio_ZZWZ_EWKCorrDown");       shape_VVratio_ZZWZ_EWKCorrDown      ->SetDirectory(0);

  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(2,   1.06531);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(3,   1.06178);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(4,   1.06921);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(5,   1.07653);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(6,   1.08430);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(7,   1.08688);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(8,   1.09050);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(9,   1.09535);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(10,  1.10661);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(11,  1.13385);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetBinContent(12,  1.16505);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(2,   2-1.06531);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(3,   2-1.06178);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(4,   2-1.06921);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(5,   2-1.07653);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(6,   2-1.08430);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(7,   2-1.08688);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(8,   2-1.09050);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(9,   2-1.09535);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(10,  2-1.10661);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(11,  2-1.13385);
  shape_signalRegion_ZZWZ_EWKCorrDown->SetBinContent(12,  2-1.16505);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(2,   1.05371);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(3,   1.06081);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(4,   1.05449);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(5,   1.06399);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(6,   1.07482);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(7,   1.07423);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(8,   1.07493);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(9,   1.08014);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(10,  1.09298);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(11,  1.10782);
  shape_VVratio_ZZWZ_EWKCorrUp->SetBinContent(12,  1.13152);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(2,   2-1.05371);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(3,   2-1.06081);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(4,   2-1.05449);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(5,   2-1.06399);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(6,   2-1.07482);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(7,   2-1.07423);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(8,   2-1.07493);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(9,   2-1.08014);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(10,  2-1.09298);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(11,  2-1.10782);
  shape_VVratio_ZZWZ_EWKCorrDown->SetBinContent(12,  2-1.13152);
  for(int nb=1; nb<=nBinMVA; nb++) {
    stat_band_WZ         -> SetBinContent(nb, 1); 
    stat_band_ZZ         -> SetBinContent(nb, 1); 
    stat_band_WZ         -> SetBinError(nb, nominal_shape_WZ        ->GetBinError(nb) / nominal_shape_WZ        ->GetBinContent(nb));
    stat_band_ZZ         -> SetBinError(nb, nominal_shape_ZZ        ->GetBinError(nb) / nominal_shape_ZZ        ->GetBinContent(nb));
  }


  TLine *oneline = new TLine;
  oneline->SetLineColor(13);
  oneline->SetLineWidth(2);
  TH1F* stat_band_ZZWZ=(TH1F*)stat_band_ZZ->Clone("stat_band_ZZWZ");
  stat_band_ZZWZ->Divide(stat_band_WZ);
  stat_band_ZZWZ->SetFillStyle(3004);
  stat_band_ZZWZ->SetFillColor(1);
  
  // Plot the systematic bands
  // WZ
  TCanvas *c_ZZWZ_EWK_signalRegion = new TCanvas("c_ZZWZ_EWK_signalRegion","c_ZZWZ_EWK_signalRegion");
  shape_signalRegion_ZZWZ_EWKCorrUp  ->SetLineWidth( 3 );
  shape_signalRegion_ZZWZ_EWKCorrDown->SetLineWidth( 3 );
  shape_signalRegion_ZZWZ_EWKCorrUp  ->SetLineColor( kViolet-1 );
  shape_signalRegion_ZZWZ_EWKCorrDown->SetLineColor( kBlue+1   );
  shape_signalRegion_ZZWZ_EWKCorrUp->SetMinimum(0.8);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetMaximum(1.2);
  shape_signalRegion_ZZWZ_EWKCorrUp->SetTitle("BDT syst. unc. from EWK VV ratio in signal region");
  shape_signalRegion_ZZWZ_EWKCorrUp->Draw("HIST");
  shape_signalRegion_ZZWZ_EWKCorrUp->GetXaxis()->SetTitle("BDT[ #font[12]{qq/gg}#rightarrowZH(125) ]");
  shape_signalRegion_ZZWZ_EWKCorrUp->GetYaxis()->SetTitle("Ratio to nominal");
  shape_signalRegion_ZZWZ_EWKCorrUp->GetXaxis()->SetRangeUser(0.2,0.899);
  
  shape_signalRegion_ZZWZ_EWKCorrDown->Draw("HIST SAME");
  stat_band_ZZWZ->Draw("E2 SAME");
  oneline->DrawLine(0.2,1,0.899,1);
  TLegend *l_bdt_ZZWZ_EWK_signalRegion = new TLegend(.35,.15,.55,.3);
  l_bdt_ZZWZ_EWK_signalRegion->SetFillColor(0);
  l_bdt_ZZWZ_EWK_signalRegion->AddEntry(shape_signalRegion_ZZWZ_EWKCorrUp, "Electroweak up", "l");
  l_bdt_ZZWZ_EWK_signalRegion->AddEntry(shape_signalRegion_ZZWZ_EWKCorrDown, "Electroweak down", "l");
  l_bdt_ZZWZ_EWK_signalRegion->AddEntry(stat_band_ZZWZ, "Stat. Unc. on Bins", "f");
  l_bdt_ZZWZ_EWK_signalRegion->Draw("SAME");
  c_ZZWZ_EWK_signalRegion->Print("syst_BDT_ZZWZ_EWK_signalRegion.pdf");


}
void METvsBDT(TString input_filename="zll3hinvll1j_sm_input_13TeV2016.root") {
  gStyle->SetPalette(kViridis);
  double xmin=0.1, xmax=0.9;
  TFile *input_file=TFile::Open(input_filename,"READ");
  TLine *line_training=new TLine(xmin,130,xmax,130);
  TLine *line_sr1=new TLine(0.2,100,0.2,300);
  TLine *line_sr2=new TLine(0.2,100,xmax,100);
  line_training->SetLineColor(kViolet-4);
  line_training->SetLineWidth(2);
  line_sr1->SetLineColor(kOrange+8);
  line_sr1->SetLineWidth(2);
  line_sr2->SetLineColor(kOrange+8);
  line_sr2->SetLineWidth(2);
  TPaveText *text_training=new TPaveText(xmin+.7*xmax,130,xmax,142,"NB");
  text_training->AddText("Training Region");
  text_training->SetFillStyle(0);
  text_training->SetTextColor(kViolet-4);
  text_training->SetBorderSize(0);
  TPaveText *text_sr=new TPaveText(xmin+0.7*xmax,100,xmax,112,"NB");
  text_sr->AddText("Signal Region");
  text_sr->SetFillStyle(0);
  text_sr->SetTextColor(kOrange+8);
  text_sr->SetBorderSize(0);

  TH2D *METvsBDT_[7]; TCanvas *canvas_[7];
  METvsBDT_[0] = (TH2D*) input_file->Get("METvsBDT_Zjets");        canvas_[0] = new TCanvas("canvas0","E_{T}^{miss} vs. BDT for Z+jets events", 800,600); 
  METvsBDT_[1] = (TH2D*) input_file->Get("METvsBDT_VVV");          canvas_[1] = new TCanvas("canvas1","E_{T}^{miss} vs. BDT for VVV events"   , 800,600); 
  METvsBDT_[2] = (TH2D*) input_file->Get("METvsBDT_WZ");           canvas_[2] = new TCanvas("canvas2","E_{T}^{miss} vs. BDT for WZ events"    , 800,600); 
  METvsBDT_[3] = (TH2D*) input_file->Get("METvsBDT_ZZ");           canvas_[3] = new TCanvas("canvas3","E_{T}^{miss} vs. BDT for ZZ events"    , 800,600); 
  METvsBDT_[4] = (TH2D*) input_file->Get("METvsBDT_EM");           canvas_[4] = new TCanvas("canvas4","E_{T}^{miss} vs. BDT for EM events"    , 800,600); 
  METvsBDT_[5] = (TH2D*) input_file->Get("METvsBDT_ggZH_hinv");    canvas_[5] = new TCanvas("canvas5","E_{T}^{miss} vs. BDT for ggZH events"  , 800,600); 
  METvsBDT_[6] = (TH2D*) input_file->Get("METvsBDT_ZH_hinv_sm");   canvas_[6] = new TCanvas("canvas6","E_{T}^{miss} vs. BDT for qqZH events"  , 800,600); 
  METvsBDT_[0]->SetDirectory(0);
  METvsBDT_[1]->SetDirectory(0);
  METvsBDT_[2]->SetDirectory(0);
  METvsBDT_[3]->SetDirectory(0);
  METvsBDT_[4]->SetDirectory(0);
  METvsBDT_[5]->SetDirectory(0);
  METvsBDT_[6]->SetDirectory(0);

  for(int i=0; i<7; i++) {
    canvas_[i]->cd();
    canvas_[i]->SetLogz();
    canvas_[i]->MoveOpaque();
    canvas_[i]->SetRightMargin(0.35);
    canvas_[i]->SetBottomMargin(0.2);
    for(int j=1; j<= METvsBDT_[i]->GetNbinsX(); j++) { for(int k=1; k<= METvsBDT_[i]->GetNbinsY(); k++) {
      if(METvsBDT_[i]->GetBinContent(METvsBDT_[i]->GetBin(j,k)) < 0.1) METvsBDT_[i]->SetBinContent(METvsBDT_[i]->GetBin(j,k),0.1);
    }}
    METvsBDT_[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    //METvsBDT_[i]->SetMinimum(TMath::Min(METvsBDT_[i]->GetMinimum(),-0.01));
    METvsBDT_[i]->SetMinimum(TMath::Min(METvsBDT_[i]->GetMinimum(),0.1));
    METvsBDT_[i]->SetMaximum(110);
    METvsBDT_[i]->GetXaxis()->SetTitle("#splitline{Multiclass Gradient BDT (400 Trees, Depth 4, }{Shrinkage 0.5, Bagging 0.5, Training E_{T}^{miss}>130}");
    METvsBDT_[i]->GetXaxis()->SetTitleOffset(1.7);
    METvsBDT_[i]->GetYaxis()->SetTitle("E_{T}^{miss} [GeV]");
    METvsBDT_[i]->GetYaxis()->SetTitleOffset(1.3);
    METvsBDT_[i]->Draw("COLZ");
    line_training->Draw("SAME");
    line_sr1->Draw("SAME");
    line_sr2->Draw("SAME");
    text_sr->Draw("SAME");
    text_training->Draw("SAME");
    canvas_[i]->Print(Form("%s.png",METvsBDT_[i]->GetName()));
  }
}

void mutual_information(TString input_filename="~/cms/cmssw/045/CMSSW_8_0_12/src/MitZHAnalysis/mva/training_result_BDT_multiClass_sm__met130_ntrees400_nodesize5pc_maxdepth4_shrinkage0p5_baggedfrac0p5_ncuts1000_14vars.root",TString theClass="Signal") {
  TFile *training_file=TFile::Open(input_filename,"READ");
  vector<TString> cleaned_var_names_, nice_var_names_;
  //cleaned_var_names_.push_back("mva_balance");                     nice_var_names_.push_back("Balance"                        );  
  //cleaned_var_names_.push_back("mva_cos_theta_star_l1");           nice_var_names_.push_back("cos #theta^{*}_{l1}"            );  
  cleaned_var_names_.push_back("TMath_Abs_mva_cos_theta_CS_l1_");  nice_var_names_.push_back("|cos #theta^{CS}_{l1}|"       );  
  cleaned_var_names_.push_back("mva_delphi_ptll_MET");             nice_var_names_.push_back("#Delta#Phi(p^{ll}, p^{miss})"   );  
  cleaned_var_names_.push_back("mva_deltaR_ll");                   nice_var_names_.push_back("#DeltaR(l1, l2)"                );  
  cleaned_var_names_.push_back("TMath_Abs_mva_etal1_");            nice_var_names_.push_back("|#eta_{l1}|"                    );  
  cleaned_var_names_.push_back("TMath_Abs_mva_etal2_");            nice_var_names_.push_back("|#eta_{l2}|"                    );  
  cleaned_var_names_.push_back("mva_MET");                         nice_var_names_.push_back("E_{T}^{miss}"                   );  
  cleaned_var_names_.push_back("mva_mll_minus_mZ");                nice_var_names_.push_back("|m_{ll} - m_{Z}|"               );  
  cleaned_var_names_.push_back("mva_mTl1MET");                     nice_var_names_.push_back("m_{T}(p_{T}^{l1}, p_{T}^{miss})");  
  cleaned_var_names_.push_back("mva_mTl2MET");                     nice_var_names_.push_back("m_{T}(p_{T}^{l2}, p_{T}^{miss})");  
  cleaned_var_names_.push_back("mva_ptl1");                        nice_var_names_.push_back("p_{T}^{l1}"                     );  
  cleaned_var_names_.push_back("mva_ptl2");                        nice_var_names_.push_back("p_{T}^{l2}"                     );  
  cleaned_var_names_.push_back("mva_ptll");                        nice_var_names_.push_back("p_{T}^{ll}"                     );  
  //cleaned_var_names_.push_back("mva_ptl1mptl2_over_ptll");         nice_var_names_.push_back("Lepton Balance"                 );  
  unsigned int num_vars = (unsigned int)cleaned_var_names_.size();
  TH2D *histo_mutual = new TH2D("histo_mutual", Form("Mutual information of variables (%s) #int#int p(x,y) log_{2} #frac{p(x,y)}{p(x)p(y)} dx dy",theClass.Data()), num_vars,0,num_vars, num_vars,0,num_vars);
  TString theTransformation="Id"; //"Gauss_Deco"
  //TString theTransformation="Gauss_Deco"; 
  //printf("%40s %40s Mutual Information\n","Variable 1", "Variable 2");
  for(unsigned int i=0; i<num_vars; i++) { for(unsigned int j=i; j<num_vars; j++) {
    if(i==j) {
      double entropy=0;
      TH1F *the_pdf = (TH1F*) training_file->Get(Form("InputVariables_%s/%s__%s_%s", theTransformation.Data(), cleaned_var_names_[i].Data(), theClass.Data(), theTransformation.Data()));
      the_pdf->Scale(1./the_pdf->Integral());
      for(unsigned int xbin=1; xbin<=(unsigned int)the_pdf->GetNbinsX(); xbin++) { 
        double px  = the_pdf->GetBinContent(xbin);
        if(px != 0) entropy -= px * TMath::Log2( px );
      }
      histo_mutual->SetBinContent( histo_mutual->GetBin(i+1,i+1), entropy);
      continue;
    }
    TH2F *joint_distribution = (TH2F*) training_file->Get(Form("InputVariables_%s/CorrelationPlots/scat_%s_vs_%s_%s_%s", theTransformation.Data(), cleaned_var_names_[j].Data(), cleaned_var_names_[i].Data(), theClass.Data(), theTransformation.Data()));
    if(!joint_distribution) { 
      printf("Null pointer, trying InputVariables_%s/CorrelationPlots/scat_%s_vs_%s_Signal_%s\n", theTransformation.Data(), cleaned_var_names_[i].Data(), cleaned_var_names_[j].Data(), theTransformation.Data());
      joint_distribution = (TH2F*) training_file->Get(Form("InputVariables_%s/CorrelationPlots/scat_%s_vs_%s_Signal_%s", theTransformation.Data(), cleaned_var_names_[i].Data(), cleaned_var_names_[j].Data(), theTransformation.Data()));
    }
    unsigned int xbins = (unsigned int) joint_distribution->GetNbinsX();
    unsigned int ybins = (unsigned int) joint_distribution->GetNbinsY();
    joint_distribution->Scale(1./joint_distribution->Integral());
    TH1D *projectionX = joint_distribution->ProjectionX(Form("marginal_%s_in_%s", cleaned_var_names_[i].Data(), cleaned_var_names_[j].Data()),0,-1,"e");
    TH1D *projectionY = joint_distribution->ProjectionX(Form("marginal_%s_in_%s", cleaned_var_names_[j].Data(), cleaned_var_names_[i].Data()),0,-1,"e");
    double mutual = 0;
    for(unsigned int xbin=1; xbin<=xbins; xbin++) { for(unsigned int ybin=1; ybin<=ybins; ybin++) {
      double pxy = joint_distribution->GetBinContent(joint_distribution->GetBin(xbin,ybin));
      double px  = projectionX->GetBinContent(xbin);
      double py  = projectionY->GetBinContent(ybin);
      if(pxy == 0 || px==0 || py==0) continue;
      double d_mutual = pxy * TMath::Log2( pxy/(px*py));
      if(d_mutual!=d_mutual) continue;
      mutual += d_mutual;
    }}
    //printf("%40s %40s %f\n",cleaned_var_names_[j].Data(),cleaned_var_names_[i].Data(),mutual);
    histo_mutual->SetBinContent( histo_mutual->GetBin(i+1, j+1), mutual);
    histo_mutual->SetBinContent( histo_mutual->GetBin(j+1, i+1), mutual);
  }}
  for(unsigned int i=0; i<num_vars; i++) {
    histo_mutual->GetXaxis()->SetBinLabel(i+1, nice_var_names_[i].Data());
    histo_mutual->GetYaxis()->SetBinLabel(i+1, nice_var_names_[i].Data());
  }
  TCanvas *canvas = new TCanvas("canvas","canvas", 800,800);
  canvas->SetLeftMargin(0.2);
  canvas->SetBottomMargin(0.18);
  canvas->SetTopMargin(0.13);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("3.2f");
  histo_mutual->LabelsOption("v","x");
  histo_mutual->SetMinimum(0.);
  histo_mutual->SetMaximum(12.);
  histo_mutual->Draw("TEXT COLZ");
  canvas->Print(Form("mutual_information_MVAvars_%s.pdf",theClass.Data()));
  canvas->Print(Form("mutual_information_MVAvars_%s.png",theClass.Data()));
  
  unsigned int max_vars=12;
  vector<unsigned int> the_vars;
  the_vars.push_back(5);
  //the_vars.push_back(7);
  while(the_vars.size() < max_vars) {
    printf("Good var set: { ");
    for(unsigned int j=0; j<the_vars.size(); j++) printf("%s ", cleaned_var_names_[the_vars[j]].Data());
    printf("}\n");
    int best_next_var=-1;
    double lowest_mutual_information=10000;
    for(unsigned int i=0; i<num_vars; i++) {
      // Skip this var if it is already in the list of good vars
      if(std::find(the_vars.begin(), the_vars.end(), i) != the_vars.end()) continue;
      double sum_mutual_information=0;
      // Sum the mutual information of this var with the vars already in the "good var" set
      for(unsigned int j=0; j<the_vars.size(); j++) {
        sum_mutual_information += histo_mutual->GetBinContent(histo_mutual->GetBin(i+1,the_vars[j]+1));
      }
      printf("Trying var %s, sum mutual information is %f with existing set\n", cleaned_var_names_[i].Data(), sum_mutual_information);
      if(sum_mutual_information < lowest_mutual_information) {
        best_next_var=i;
        lowest_mutual_information = sum_mutual_information;
      }
    }
    assert(best_next_var!=-1);
    printf("adding variable %s, sum of mutual information with existing set is %f\n", cleaned_var_names_[best_next_var].Data(), lowest_mutual_information);
    the_vars.push_back(best_next_var);
  }
  printf("Final good var set:\n");
  for(unsigned int j=0; j<the_vars.size(); j++) printf("%d. %s\n", j+1, cleaned_var_names_[the_vars[j]].Data());
}













*/
