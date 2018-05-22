#include <Compression.h>
#include "TBranch.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
void npnloAddBranch(TString inputFileName, TString lookupFileName) {
  unsigned char npnlo;
  ULong64_t eventNumber;
  int lumiNumber;
  TFile *lookupFile = TFile::Open(lookupFileName,"read"); assert(lookupFile);
  TTree *npnloLookup = (TTree*) lookupFile->Get("npnloLookup");
  npnloLookup->BuildIndex("lumiNumber","eventNumber");
  npnloLookup->SetBranchAddress("npnlo",&npnlo);

  TFile *inputFile = TFile::Open(inputFileName,"update");
  TTree *events = (TTree*)inputFile->Get("events"); assert(events);
  events->SetBranchAddress("eventNumber",&eventNumber);
  events->SetBranchAddress("lumiNumber",&lumiNumber);
  TBranch *branch_npnlo = events->Branch("npnlo",&npnlo);
  Long64_t nEntries = events->GetEntries();
  for(Long64_t i=0; i<nEntries; i++) {
    int readLookup = npnloLookup->GetEntryWithIndex(lumiNumber,eventNumber);
    if(readLookup<=0) npnlo=255;
    branch_npnlo->Fill();
  }
  events->Write();
  inputFile->Close();
  lookupFile->Close();
}


