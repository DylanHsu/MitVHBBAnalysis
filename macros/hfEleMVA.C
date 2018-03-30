#include <TROOT.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>

void hfEleMVA(
  TString inputFileName, 
  TString extraString="", 
  bool useGaussDeco=false
) {
  gROOT->ProcessLine("TMVA::gConfig().GetVariablePlotting().fMaxNumOfAllowedVariablesForScatterPlots = 50");
  TFile *output_file;
  TMVA::Factory *factory;
  
  // Determine the input trees
  TFile *inputFile = TFile::Open(inputFileName, "READ"); assert(inputFile);
  TTree *inputTree = (TTree*)inputFile->Get("tree"); assert(inputTree);
  
  // Initialize the factory
  output_file=TFile::Open(TString("/data/t3home000/dhsu/mva/whbb/trainingResult_hfEleMva_"+extraString+".root"), "RECREATE");
  TString factoryOptions="!V:!Silent:DrawProgressBar";
  factoryOptions+=":AnalysisType=Classification";
  if(useGaussDeco)  factoryOptions += ":Transformations=G,D";
  factory = new TMVA::Factory("hfEleMva", output_file, factoryOptions);
  factory->AddTree(inputTree,"Signal"    ,1.0,"ele_matchedGHFE!=0","test train");
  factory->AddTree(inputTree,"Background",1.0,"ele_matchedGHFE==0","test train");
  factory->SetWeightExpression("1", "Signal");
  factory->SetWeightExpression("1", "Background");
  //factory->SetWeightExpression("weight", "Signal");
  //factory->SetWeightExpression("weight", "Background");
  
  TCut preselectionCut ="ele_conversionVeto && ele_ooEmooP<10 && ele_ooEmooP>=0 && ele_dxy<10 && ele_dz<10 && ele_r9>=0 && ele_r9<=1 ";
  //factory->AddVariable("dEta_eleFatjet"     , "#Delta#eta(e,fj)"           , ""     , 'F');
  //factory->AddVariable("dPhi_eleFatjet"     , "#Delta#phi(e,fj)"           , "Rad"  , 'F');
  factory->AddVariable("ele_eta"            , "eta"                        , ""     , 'F');
  factory->AddVariable("ele_ecalIso/ele_pt" , "Rel ECAL iso"               , ""     , 'F');
  factory->AddVariable("ele_hcalIso/ele_pt" , "Rel HCAL iso"               , ""     , 'F');
  factory->AddVariable("ele_trackIso/ele_pt", "Rel trk iso"                , ""     , 'F');
  factory->AddVariable("ele_chIso/ele_pt"   , "Rel CH iso"                 , ""     , 'F');
  factory->AddVariable("ele_nhIso/ele_pt"   , "Rel NH iso"                 , ""     , 'F');
  factory->AddVariable("ele_phIso/ele_pt"   , "Rel PH iso"                 , ""     , 'F');
  factory->AddVariable("ele_dxy"            , "dxy"                        , "cm"   , 'F');
  factory->AddVariable("ele_dz"             , "dz"                         , "cm"   , 'F');
  factory->AddVariable("ele_sieie"          , "#sigmai#etai#eta"           , ""     , 'F');
  factory->AddVariable("ele_sipip"          , "#sigmai#phii#phi"           , ""     , 'F');
  factory->AddVariable("ele_r9"             , "R9"                         , ""     , 'F');
  factory->AddVariable("ele_dEtaInSeed"     , "#Delta#eta(SC,trk)"         , ""     , 'F');
  factory->AddVariable("ele_dPhiIn"         , "#Delta#phi(SC,trk)"         , ""     , 'F');
  factory->AddVariable("ele_hOverE"         , "H/E"                        , ""     , 'F');
  factory->AddVariable("ele_ooEmooP"        , "|1-E_{SC}/p_{trk}|/E_{ECAL}", "1/GeV", 'F');
  factory->AddVariable("ele_nMissingHits"   , "N missing hits"             , ""     , 'I');
  factory->AddVariable("ele_tripleCharge"   , "triple charge"              , ""     , 'I');
  
  TString prepareOptions="NormMode=NumEvents";
    prepareOptions+=":SplitMode=Random";
    prepareOptions+=":MixMode=Random";
  factory->PrepareTrainingAndTestTree(preselectionCut, prepareOptions);
  vector<TString> hyperparameters_ = {
    //"NTrees=500:MinNodeSize=2.5%:MaxDepth=4:SeparationType=CrossEntropy:BoostType=Grad:Shrinkage=0.1:nCuts=1000",
    //"NTrees=500:MinNodeSize=2.5%:MaxDepth=4:SeparationType=CrossEntropy:BoostType=Grad:Shrinkage=0.2:nCuts=1000",
    //"NTrees=500:MinNodeSize=2.5%:MaxDepth=4:SeparationType=CrossEntropy:BoostType=Grad:Shrinkage=0.3:nCuts=1000",
    "NTrees=500:MinNodeSize=2.5%:MaxDepth=5:SeparationType=CrossEntropy:BoostType=Grad:Shrinkage=0.1:nCuts=1000"
    //"NTrees=500:MinNodeSize=2.5%:MaxDepth=5:SeparationType=CrossEntropy:BoostType=Grad:Shrinkage=0.2:nCuts=1000",
    //"NTrees=500:MinNodeSize=2.5%:MaxDepth=5:SeparationType=CrossEntropy:BoostType=Grad:Shrinkage=0.3:nCuts=1000"
  };
  for(unsigned hp=0; hp<hyperparameters_.size();hp++) {
    TString hyperparameters = hyperparameters_[hp];
    hyperparameters+="!H:!V:NegWeightTreatment=Pray";
    if(useGaussDeco) hyperparameters += ":VarTransform=G,D";
    factory->BookMethod(TMVA::Types::kBDT, Form("hfEleMva_hp%d_%s",hp,extraString.Data()), hyperparameters);
  }
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  output_file->Close();
}
