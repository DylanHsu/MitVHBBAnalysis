#include <Compression.h>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

void npnloLookup(TString lfn, TString outputName) {
  TString outDir="/data/t3home000/dhsu/npnloLookup/miniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/";
  system(Form("mkdir -p %s",outDir.Data()));
  TString inputFileName="root://cmsxrootd.fnal.gov/"+lfn;
  TFile *inputFile=0;
  int retries=0;
  while(true) {
    inputFile = TFile::Open(inputFileName,"read");
    if(inputFile && inputFile->IsOpen()) break;
    retries++;
    if(retries>100) { throw std::runtime_error("Error opening input file"); return; }
  }
  
  
  edm::Wrapper<LHEEventProduct> *lheWrapper=0;
  edm::EventAuxiliary *aux=0;

  TTree *inputTree = (TTree*)inputFile->Get("Events");
  inputTree->SetBranchStatus("*",0);
  inputTree->SetBranchStatus("LHEEventProduct*",true);
  inputTree->SetBranchStatus("EventAuxiliary",true);

  inputTree->SetBranchAddress("LHEEventProduct_externalLHEProducer__SIM.",&lheWrapper);
  inputTree->SetBranchAddress("EventAuxiliary",&aux);

  unsigned char npnlo;
  ULong64_t eventNumber;
  unsigned int lumiNumber;
  system("mkdir -p /tmp/dhsu");
  TFile *outputFile = TFile::Open("/tmp/dhsu/"+outputName,"recreate","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  TTree *outputTree = new TTree("npnloLookup","npnloLookup");
  outputTree->Branch("eventNumber", &eventNumber);
  outputTree->Branch("lumiNumber" , &lumiNumber );
  outputTree->Branch("npnlo"      , &npnlo      );

  for(Long64_t iEntry=0; iEntry<inputTree->GetEntries(); iEntry++) {
    inputTree->GetEntry(iEntry);
    const LHEEventProduct* theLHE = lheWrapper->product();
    npnlo=theLHE->npNLO();
    eventNumber=aux->event();
    lumiNumber=aux->luminosityBlock();
    //printf("%d %llu %u\n", npnlo, eventNumber,lumiNumber);
    outputTree->Fill();
  }
  outputTree->Write();
  outputFile->Close();
  inputFile->Close();
  system(Form("mv %s %s",("/tmp/dhsu/"+outputName).Data(),outDir.Data()));

}
