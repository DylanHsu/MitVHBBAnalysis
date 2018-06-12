#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TString.h>
#include <iostream>
void splitPandaExpress(string inputFile, Long64_t nEvtsPerFile=500000) {
  printf("Opening file \"%s\" and splitting it into files with %lld events...\n", inputFile.c_str(), nEvtsPerFile);
  TFile *f = TFile::Open(inputFile.c_str(),"update"); assert(f);
  TTree *events=0, *treebuffer;
  
  // Drop extra cycle numbers so rooteventselector can do its job
  printf("\tLooking for extra cycle numbers to drop...\n");
  TIter keyList(f->GetListOfKeys());
  vector<short> extraCycleNumbers; short highestCycleNumber=-1;
  TKey *key;
  while ((key = (TKey*)keyList())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree")) continue;
    treebuffer=(TTree*)key->ReadObj();
    if(strcmp(treebuffer->GetName(),"events")!=0) continue;
    short cycleNumber=key->GetCycle();
    printf("\t%s;%d\n",treebuffer->GetName(),cycleNumber);
    if(cycleNumber > highestCycleNumber) {
      // replace the highest cycle number, putting the previous highest in the list of extras if real
      if(highestCycleNumber>0) extraCycleNumbers.push_back(highestCycleNumber);
      highestCycleNumber=cycleNumber;
    } else {
      // this is an extra cycle number
      extraCycleNumbers.push_back(cycleNumber);
    }
  }
  std::cout << "Extra cycle numbers: ";
  for (std::vector<short>::const_iterator i = extraCycleNumbers.begin(); i != extraCycleNumbers.end(); ++i)
    std::cout << *i << ", ";
  std::cout << std::endl;
  for (std::vector<short>::const_iterator i = extraCycleNumbers.begin(); i != extraCycleNumbers.end(); ++i)
    f->Delete(Form("events;%d",*i));
  // Done dropping extra cycle numbers
 
  
  events = (TTree*)f->Get(Form("events;%d",highestCycleNumber)); assert(events);
  Long64_t nEntries = events->GetEntries();
  unsigned nSplit = ceil( float(nEntries) / nEvtsPerFile);
  f->Close();
  printf("The file had %lld entries, splitting into %d subfiles\n", nEntries, nSplit);
  size_t lastDot = inputFile.find_last_of(".");
  size_t lastSlash = inputFile.find_last_of("/");
  string splitDir = inputFile.substr(0, lastSlash) + string("/split/");
  string rawName = splitDir + inputFile.substr(lastSlash+1, lastDot-lastSlash-1);
  system(Form("mkdir -p %s",splitDir.c_str()));
  for(unsigned i=0; i<nSplit; i++) {
    string splitName;
    if(lastDot!=string::npos) splitName = rawName + string(Form("_%d.root", i));
    else                      splitName = rawName + string(Form("_%d", i));
    Long64_t firstEvt = i*nEvtsPerFile;
    Long64_t lastEvt  = (i+1)*nEvtsPerFile - 1;
    const char * command = Form("rooteventselector --recreate -f %lld -l %lld %s:events %s", firstEvt, lastEvt, inputFile.c_str(), splitName.c_str());
    printf("\tSubfile #%d: writing event # %lld-%lld to \"%s\"\n", i+1, firstEvt, lastEvt, splitName.c_str());
    system(command);
    //printf("\t%s\n",command);
  }
  printf("Done splitting file \"%s\"\n", inputFile.c_str());
}
