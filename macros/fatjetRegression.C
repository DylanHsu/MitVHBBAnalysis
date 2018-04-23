#include <TROOT.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include "vhbbPlot.h"


void fatjetRegression(TString inputFileName, TString extraString="") {
  TString outfileName(Form("/data/t3home000/dhsu/trainingResult_fatjetRegression_%s.root",extraString.Data()) );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  //TMVA::Factory *factory = new TMVA::Factory( TString("fatjetRegression"+extraString), outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression:Transformations=G,D" );
  TMVA::Factory *factory = new TMVA::Factory( TString("fatjetRegression"+extraString), outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
  //TString fname1 = "/data/t3home000/dhsu/fatjetRegTree/merged/Hbb.root";
  //TString fname1 = "/data/t3home000/dhsu/fatjetRegTree/merged/QCD.root";
  //TString fname = "/data/t3home000/dhsu/fatjetRegTree/merged/TTbar_Powheg.root";
  //TString fname1= "/data/t3home000/dhsu/fatjetRegTree/flattened/QCD_merged_skimmed.root";
  //TString fname2= "/data/t3home000/dhsu/fatjetRegTree/flattened/WHbb_skimmed.root";
  TString fname1 = inputFileName;
  TFile *input1 = TFile::Open( fname1); // check if file in local directory exists
  //TFile *input2 = TFile::Open( fname2); // check if file in local directory exists
  if (!input1) {// || !input2) {
     std::cout << "ERROR: could not open data file" << std::endl;
     exit(1);
  }
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("MitVHBBAnalysis");
  dataloader->AddRegressionTree( (TTree*)input1->Get("regTree"), 1.0 );
  //dataloader->AddRegressionTree( (TTree*)input2->Get("regTree"), 1.0 );
  dataloader->AddVariable("fj1_pt"                , "fj1_pt"                  , "", 'F');
  dataloader->AddVariable("fj1_MSD"               , "fj1_MSD"                 , "", 'F');
  dataloader->AddVariable("fj1_nhf"               , "fj1_nhf"                 , "", 'F');
  dataloader->AddVariable("fj1_chf"               , "fj1_chf"                 , "", 'F');
  dataloader->AddVariable("fj1_nef"               , "fj1_nef"                 , "", 'F');
  dataloader->AddVariable("fj1_cef"               , "fj1_cef"                 , "", 'F');
  dataloader->AddVariable("fj1_emFraction"        , "fj1_emFraction"          , "", 'F');
  dataloader->AddVariable("fj1_nLep"              , "fj1_nLep"                , "", 'I');
  dataloader->AddVariable("ele1_pt"               , "ele1_pt"                 , "", 'F');
  dataloader->AddVariable("ele1_ptRelFj1"         , "ele1_ptRelFj1"           , "", 'F');
  dataloader->AddVariable("ele1_dRFj1"            , "ele1_dRFj1"              , "", 'F');
  dataloader->AddVariable("ele1_ptRelSv1"         , "ele1_ptRelSv1"           , "", 'F');
  dataloader->AddVariable("ele1_dRSv1"            , "ele1_dRSv1"              , "", 'F');
  dataloader->AddVariable("ele1_ptRelSv2"         , "ele1_ptRelSv2"           , "", 'F');
  dataloader->AddVariable("ele1_dRSv2"            , "ele1_dRSv2"              , "", 'F');
  dataloader->AddVariable("ele1_dxy"              , "ele1_dxy"                , "", 'F');
  dataloader->AddVariable("ele1_dz"               , "ele1_dz"                 , "", 'F');
  dataloader->AddVariable("ele1_hfMva"            , "ele1_hfMva"              , "", 'F');
  dataloader->AddVariable("mu1_pt"                , "mu1_pt"                  , "", 'F');
  dataloader->AddVariable("mu1_ptRelFj1"          , "mu1_ptRelFj1"            , "", 'F');
  dataloader->AddVariable("mu1_dRFj1"             , "mu1_dRFj1"               , "", 'F');
  dataloader->AddVariable("mu1_ptRelSv1"          , "mu1_ptRelSv1"            , "", 'F');
  dataloader->AddVariable("mu1_dRSv1"             , "mu1_dRSv1"               , "", 'F');
  dataloader->AddVariable("mu1_ptRelSv2"          , "mu1_ptRelSv2"            , "", 'F');
  dataloader->AddVariable("mu1_dRSv2"             , "mu1_dRSv2"               , "", 'F');
  dataloader->AddVariable("mu1_dxy"               , "mu1_dxy"                 , "", 'F');
  dataloader->AddVariable("mu1_dz"                , "mu1_dz"                  , "", 'F');
  dataloader->AddVariable("fj1_sv1_pt"            , "fj1_sv1_pt"              , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_ptRelFj1"      , "fj1_sv1_ptRelFj1"        , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_dRFj1"         , "fj1_sv1_dRFj1"           , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_m"             , "fj1_sv1_m"               , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_logMRel"       , "fj1_sv1_logMRel"         , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_vtx3DVal"      , "fj1_sv1_vtx3DVal"        , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_vtx3DeVal"     , "fj1_sv1_vtx3DeVal"       , "", 'F'); 
  dataloader->AddVariable("fj1_sv1_ntrk"          , "fj1_sv1_ntrk"            , "", 'I'); 
  dataloader->AddVariable("fj1_sv1_significance"  , "fj1_sv1_significance"    , "", 'F');  
  dataloader->AddVariable("fj1_sv2_pt"            , "fj1_sv2_pt"              , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_ptRelFj1"      , "fj1_sv2_ptRelFj1"        , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_dRFj1"         , "fj1_sv2_dRFj1"           , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_m"             , "fj1_sv2_m"               , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_logMRel"       , "fj1_sv2_logMRel"         , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_vtx3DVal"      , "fj1_sv2_vtx3DVal"        , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_vtx3DeVal"     , "fj1_sv2_vtx3DeVal"       , "", 'F'); 
  dataloader->AddVariable("fj1_sv2_ntrk"          , "fj1_sv2_ntrk"            , "", 'I'); 
  dataloader->AddVariable("fj1_sv2_significance"  , "fj1_sv2_significance"    , "", 'F');  
  
  // Add the variable carrying the regression target
  //dataloader->AddTarget("fj1_genMnu","fj1_genMnu", "", 'F');
  //dataloader->AddTarget("TMath::Log10(fj1_sv1_ptWithNu/fj1_sv1_pt)","log(sv1 pt corr)", "", 'F');
  dataloader->AddTarget("log(fj1_genMnu/fj1_MSD)","log(kM)", "", 'F');
  // It is also possible to declare additional targets for multi-dimensional regression, ie:
  //     dataloader->AddTarget( "fvalue2" );
  // BUT: this is currently ONLY implemented for MLP
  // Read training and test data (see TMVAClassification for reading ASCII files)
  // load the signal and background event samples from ROOT trees
  //dataloader->SetWeightExpression( "normalizedWeight*steamroller", "Regression" );
  dataloader->SetWeightExpression( "normalizedWeight", "Regression" );
  TCut preselectionCut = "fj1_genMnu/fj1_MSD>0 && fj1_genMnu/fj1_MSD<10 && fj1_genMnu<250";
  //TCut preselectionCut = "fj1_genMnu<250 && fj1_sv1_ptWithNu/fj1_sv1_pt<10";
  TString prepareOptions="NormMode=NumEvents";
    prepareOptions+=":SplitMode=Random"; 
    prepareOptions+=":MixMode=Random";
    //prepareOptions+=":nTrain_Regression=1000000:nTest_Regression=200000";
    //prepareOptions+=":nTrain_Regression=300000:nTest_Regression=120000";

  dataloader->PrepareTrainingAndTestTree(preselectionCut, prepareOptions);
  //Neural network (MLP)
  //factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );
  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NegWeightTreatment=Pray:BoostType=Grad:Shrinkage=0.1:SeparationType=RegressionVariance:nCuts=10000:MaxDepth=3:MinNodeSize=0.2%:NTrees=2000" );
  //factory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:BoostType=Bagging:SeparationType=RegressionVariance:nCuts=10000:PruneMethod=CostComplexity:PruneStrength=20:MaxDepth=15:NTrees=100" );
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  outputFile->Close();
}
