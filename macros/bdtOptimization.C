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
const unsigned nVarsMax=200;
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
    
  }
  sort(varMutualInfs.begin(), varMutualInfs.end(), [=](std::pair<unsigned, float>& a, std::pair<unsigned, float>& b) {
    return a.second > b.second;
  });
  printf("Rank     I(T;A)       Error  Name(A)\n");
  printf("---------------------------------------\n");
  for(unsigned i=0; i<varMutualInfs.size(); i++) {
    printf("%4d %10.2e  %10.2e  %s\n", i+1, varMutualInfs[i].second, varMutualInfErrors[varMutualInfs[i].first], varNames[varMutualInfs[i].first].Data());
    pairwiseMutualInf->SetBinContent(pairwiseMutualInf->FindBin(i,i), varMutualInfs[i].second);
    pairwiseMutualInf->SetBinError(pairwiseMutualInf->FindBin(i,i), varMutualInfErrors[varMutualInfs[i].first]);
  }
  
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
      unsigned k1,k2;
      for(unsigned i=0; i<varMutualInfs.size(); i++) {
        if(varMutualInfs[i].first==h1) k1=i;
        if(varMutualInfs[i].first==h2) k2=i;
      }
      pairwiseMutualInf->SetBinContent(pairwiseMutualInf->FindBin(k1,k2), mutualInfWithTruth);
      pairwiseMutualInf->SetBinError(pairwiseMutualInf->FindBin(k1,k2), mutualError);
      pairwiseMutualInf->SetBinContent(pairwiseMutualInf->FindBin(k2,k1), mutualInfWithTruth);
      pairwiseMutualInf->SetBinError(pairwiseMutualInf->FindBin(k2,k1), mutualError);

    }
  }
  // Set axis labels
  pairwiseMutualInf->GetXaxis()->SetLabelSize(0.025);
  pairwiseMutualInf->GetYaxis()->SetLabelSize(0.025);
  for(unsigned i=0; i<varMutualInfs.size(); i++) {
    pairwiseMutualInf->GetXaxis()->SetBinLabel(i+1, varTitles[varMutualInfs[i].first]);
    pairwiseMutualInf->GetYaxis()->SetBinLabel(i+1, varTitles[varMutualInfs[i].first]);
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
