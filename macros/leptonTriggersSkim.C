#ifndef leptonTriggersSkim_cc
#define leptonTriggersSkim_cc
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include "PandaTree/Objects/interface/Event.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <sstream>
#include "PandaAnalysis/Utilities/src/RoccoR.cc"
//#include <TRandom3.h>
typedef std::map<UInt_t,std::vector<std::pair <UInt_t, UInt_t> > > MapType;
bool leptonTriggersSkim(
  TString inputFileName, 
  TString outputFileName, 
  TString flavor="muons",
  bool debug=false
) {
  gSystem->Load("libPandaTreeObjects.so"); 
  assert(inputFileName!=outputFileName);
  assert(flavor=="muons" || flavor=="electrons");
  int lepSelType=-1;
  if(flavor=="muons") lepSelType=1;
  else if(flavor=="electrons") lepSelType=2;

  // Open input file
  printf("Opening file \"%s\"\n", inputFileName.Data());
  TFile *inputFile = TFile::Open(inputFileName,"READ");
  assert(inputFile && inputFile->IsOpen());
  TTree *inputTree = (TTree*)inputFile->Get("events"); assert(inputTree);
  
  //Read json file into boost property tree
  string jsonFile = "PandaAnalysis/data/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
//string jsonFile = "PandaAnalysis/data/certs/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt";
  MapType fMap;
  boost::property_tree::ptree jsonTree;
  boost::property_tree::read_json(jsonFile.c_str(),jsonTree);
  
  //Loop through boost property tree and fill the MapType structure with the list of good lumi
  //ranges for each run
  for (boost::property_tree::ptree::const_iterator it = jsonTree.begin(); it!=jsonTree.end(); ++it) {
    UInt_t runNumber = boost::lexical_cast<UInt_t>(it->first);
    MapType::mapped_type &lumiPairList = fMap[runNumber];
    boost::property_tree::ptree lumiPairListTree = it->second;
    for (boost::property_tree::ptree::const_iterator jt = lumiPairListTree.begin(); jt!=lumiPairListTree.end(); ++jt) {
      boost::property_tree::ptree lumiPairTree = jt->second;
      if (lumiPairTree.size()==2) {
        UInt_t firstLumi = boost::lexical_cast<UInt_t>(lumiPairTree.begin()->second.data());
        UInt_t lastLumi = boost::lexical_cast<UInt_t>((++lumiPairTree.begin())->second.data());
        lumiPairList.push_back(std::pair<UInt_t,UInt_t>(firstLumi,lastLumi));
      }
    }
  }

  //Load rochester corrections
  RoccoR rochesterCorrection("PandaAnalysis/data/rcdata.2016.v3");
  
  // Set branch addresses
  panda::Event event; event.setStatus(*inputTree, {"!*"});
  switch(lepSelType) {
    case 1: event.setAddress(*inputTree, {"muons", "triggers", "triggerObjects", "runNumber", "lumiNumber", "eventNumber"}); break;
    case 2: event.setAddress(*inputTree, {"electrons", "triggers", "triggerObjects", "runNumber", "lumiNumber", "eventNumber"}); break;
    default: assert(0); break;
  }

  // Register triggers
  unsigned refTriggerToken;
  vector<int> testTriggerEnums, refTriggerEnums;
  switch(lepSelType) {
    case 1:
      refTriggerToken=event.registerTrigger("HLT_IsoMu24");
      refTriggerEnums.push_back( panda::Muon::fIsoMu24 );
      testTriggerEnums.push_back( panda::Muon::fIsoMu24          );
      //testTriggerEnums.push_back( panda::Muon::fMu17Mu8FirstLeg  );
      //testTriggerEnums.push_back( panda::Muon::fMu17Mu8SecondLeg );
      //testTriggerEnums.push_back( panda::Muon::fIsoMu22er        );
      //testTriggerEnums.push_back( panda::Muon::fIsoTkMu22er      );
      //testTriggerEnums.push_back( panda::Muon::fIsoTkMu24        );  
      break;
    case 2:
      refTriggerToken=event.registerTrigger("HLT_Ele27_WPTight_Gsf");
      refTriggerEnums.push_back( panda::Electron::fEl27Tight );
      testTriggerEnums.push_back( panda::Electron::fEl27Tight );
      break;
    default: assert(0); break;
  }


  //declare output variables
  unsigned int out_runNum, // event ID
    out_lumiSec,
    out_evtNum,
    out_npv, // number of primary vertices
    pass; // whether probe passes requirements
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  float        met;                             // missing ET
  int          njets;                           // number of jets
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 
  TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  TTree *outputTree = new TTree("events", "Trigger skim");
  outputTree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
  outputTree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
  outputTree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
  outputTree->Branch("npv",      &out_npv,      "npv/i"      );  
  outputTree->Branch("pass",     &pass,     "pass/i"     );  
  outputTree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
  outputTree->Branch("mass",     &mass,     "mass/F"     );  
  outputTree->Branch("qtag",     &qtag,     "qtag/I"     );  
  outputTree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
  //outputTree->Branch("njets",    &njets,    "njets/I"   );  
  //outputTree->Branch("met",      &met,      "met/F"   );  
  outputTree->Branch("tag",   "TLorentzVector", &p4_tag   );  
  outputTree->Branch("probe", "TLorentzVector", &p4_probe );          
  
  long iEntry = 0, nPass=0, nFail=0;
  long nentries=inputTree->GetEntries();
  if(nentries>1000 && debug) nentries=1000;
  // Begin event loop
  while (event.getEntry(*inputTree, iEntry++) > 0 && iEntry<=nentries) {
    if(debug) printf("######## Reading entry %ld/%ld ########################################################\n",iEntry,nentries);
    else if(iEntry%100000==0) printf("######## Reading entry %ld/%ld ########################################################\n",iEntry,nentries);
    if(debug) printf("Processing runNum %d, LS %d, eventNum %llu\n", event.runNumber, event.lumiNumber, event.eventNumber);
    
    // Check data certification
    bool certifiedEvent=false;
    std::pair<unsigned int, unsigned int> runLumi(event.runNumber, event.lumiNumber);      
    MapType::const_iterator it = fMap.find(runLumi.first);
    if (it!=fMap.end()) {
      //check lumis
      const MapType::mapped_type &lumiPairList = it->second;
      for (MapType::mapped_type::const_iterator jt = lumiPairList.begin(); jt<lumiPairList.end(); ++jt) {
        if (runLumi.second >= jt->first && runLumi.second <= jt->second) {
          //found lumi in accepted range
          certifiedEvent=true;
        }
      }
    }
    if(!certifiedEvent) { if(debug) printf("Event failed data certification\n"); continue; }
    
    // Check the reference trigger
    if(!event.triggerFired(refTriggerToken)) { if(debug) printf("Event did not pass the reference trigger\n");  continue; }
    
    vector<panda::Electron*> tagElectrons, probeElectrons;
    vector<panda::Muon*> tagMuons, probeMuons;
    vector<bool> passTestTriggers_;
    // Loop over muons
    if(lepSelType==1 && event.muons.size()>0) for(unsigned iM = 0; iM != event.muons.size(); ++iM) {
      panda::Muon *muon = &event.muons[iM];
      if(!muon->tight) continue;
      double ptCorrection=rochesterCorrection.kScaleDT(muon->charge, muon->pt(), muon->eta(), muon->phi(), 0, 0);
      double pt=muon->pt()*ptCorrection;
      if(pt<25.) continue;
      muon->setPtEtaPhiM(pt,muon->eta(),muon->phi(),0.106);
      bool passRefTrigger=false, passTestTriggers=false;
      for(unsigned i=0; i<refTriggerEnums.size(); i++) passRefTrigger|=muon->triggerMatch[refTriggerEnums[i]];
      for(unsigned i=0; i<testTriggerEnums.size(); i++) passTestTriggers|=muon->triggerMatch[testTriggerEnums[i]];
      if(passRefTrigger) tagMuons.push_back(muon);
      probeMuons.push_back(muon);
      passTestTriggers_.push_back(passTestTriggers);
    }
    // Loop over electrons
    if(lepSelType==2 && event.electrons.size()>0) for(unsigned iE = 0; iE != event.electrons.size(); ++iE) {
      panda::Electron *electron = &event.electrons[iE];
      if(!electron->tight) continue;
      if(electron->smearedPt < 25.) continue;
      bool passRefTrigger=false, passTestTriggers=false;
      for(unsigned i=0; i<refTriggerEnums.size(); i++) passRefTrigger|=electron->triggerMatch[refTriggerEnums[i]];
      for(unsigned i=0; i<testTriggerEnums.size(); i++) passTestTriggers|=electron->triggerMatch[testTriggerEnums[i]];
      if(passRefTrigger) tagElectrons.push_back(electron);
      probeElectrons.push_back(electron);
      passTestTriggers_.push_back(passTestTriggers);
    }
    if(debug) printf("Found %lu tag and %lu probe leptons\n", tagElectrons.size()+tagMuons.size(), probeElectrons.size()+probeMuons.size());
    
    // Event level info
    out_runNum=event.runNumber;
    out_lumiSec=event.lumiNumber;
    out_evtNum=event.eventNumber;
    out_npv=event.npv;
    // Pair association
    if(lepSelType==1) for(unsigned iTag=0; iTag<tagMuons.size(); iTag++) for(unsigned iProbe=0; iProbe<probeMuons.size(); iProbe++) {
      panda::Muon *tagMuon = tagMuons[iTag]; panda::Muon *probeMuon = probeMuons[iProbe];
      if(tagMuon==probeMuon) continue;
      // Pair-level info
      *p4_tag = tagMuon->p4(); *p4_probe = probeMuon->p4();
      qtag = tagMuon->charge; qprobe = probeMuon->charge;
      pass=passTestTriggers_[iProbe];
      TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
      mass = systemP4.M();
      if(debug) printf("Filling a muon pair (pass=%d, mass=%.1f)\n", pass, mass);
      outputTree->Fill();
      if(pass) nPass++; else nFail++;
    }
  } // End event loop
  inputFile->Close();
  outputFile->cd(); outputTree->Write(); outputFile->Close();
  printf("Passing pairs %ld/%ld\n", nPass, nPass+nFail);
  return true;
}
#endif
