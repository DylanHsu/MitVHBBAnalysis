// does not compile in Aclic, use in CINT only
void splitPandaExpress(string inputFile, Long64_t nEvtsPerFile=100000) {
  printf("Opening file \"%s\" and splitting it into files with %lld events...\n", inputFile.c_str(), nEvtsPerFile);
  TFile *f = TFile::Open(inputFile.c_str(),"read"); assert(f);
  TTree *events = (TTree*)f->Get("events"); assert(events);
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
