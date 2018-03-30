#include "TCanvas.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
void cutsByHand(TString inputFileName, TCut cutString, int targetCategory, TString plotVar="fj1MSD", int nbins=30, float xmin=0, float xmax=300) {
  TFile *inputFile = TFile::Open(inputFileName,"read");
  TTree *tree = (TTree*)inputFile->Get("plotTree");
  tree->Draw(
    Form("%s>>hBkg(%d,%f,%f)", plotVar.Data(), nbins, xmin, xmax),
    TCut("weight") * (cutString && TCut(Form("theCategory!=%d && theCategory!=12 && theCategory>0",targetCategory))),
    "e0 hist"
  );
  TH1F *hBkg = (TH1F*)gDirectory->Get("hBkg"); hBkg->SetDirectory(0);
  int entriesBkg = hBkg->GetEntries();
  float Nbkg = hBkg->Integral();
  //hBkg->Scale(1./Nbkg);
  tree->Draw(
    Form("%s>>hSignal(%d,%f,%f)", plotVar.Data(), nbins, xmin, xmax),
    TCut("weight") * (cutString && TCut(Form("weight*(theCategory==%d)",targetCategory))),
    "e0 hist same"
  );
  TH1F *hSignal = (TH1F*)gDirectory->Get("hSignal"); hSignal->SetDirectory(0);
  int entriesSignal = hSignal->GetEntries();
  float Nsig = hSignal->Integral();
  //hSignal->Scale(1./Nsig);
  hSignal->SetLineColor(kRed);

  printf("Signal entries = %d, bkg entries = %d\n", entriesSignal, entriesBkg);
  printf("Nsig %.1f, Nbkg %.1f\n", Nsig, Nbkg);
  printf("significance = %.3e\n", Nsig/sqrt(Nsig+Nbkg));
  printf("purity       = %.3e\n", Nsig/(Nsig+Nbkg));


}
