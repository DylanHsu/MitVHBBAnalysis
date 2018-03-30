#include <TROOT.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include "vhbbPlot.h"
void fatjetRegression() {
  TString outfileName( "/data/t3home000/dhsu/fatjetRegression.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
  TString fname = "/data/t3home000/dhsu/fatjetRegTree/WHbb.root";
  TFile *input = TFile::Open( fname ); // check if file in local directory exists
  if (!input) {
     std::cout << "ERROR: could not open data file" << std::endl;
     exit(1);
  }
  float sum_weights = ((TH1F*)input->Get("sum_weights"))->Integral();
  TTree *regTree = (TTree*)input->Get("regTree");
  Double_t regWeight  = 1.0;
  factory->AddRegressionTree( regTree, regWeight );
  factory->AddVariable("fj1_pt"                , "fj1_pt"                  , "", 'F');
  factory->AddVariable("fj1_eta"               , "fj1_eta"                 , "", 'F');
  factory->AddVariable("fj1_phi"               , "fj1_phi"                 , "", 'F');
  //factory->AddVariable("fj1_MSD"               , "fj1_MSD"                 , "", 'F');
  //factory->AddVariable("fj1_Mnu"               , "fj1_Mnu"                 , "", 'F');
  //factory->AddVariable("fj1_genM"              , "fj1_genM"                , "", 'F');
  //factory->AddVariable("fj1_genMnu"            , "fj1_genMnu"              , "", 'F');
  factory->AddVariable("fj1_nhf"               , "fj1_nhf"                 , "", 'F');
  factory->AddVariable("fj1_chf"               , "fj1_chf"                 , "", 'F');
  factory->AddVariable("fj1_nef"               , "fj1_nef"                 , "", 'F');
  factory->AddVariable("fj1_cef"               , "fj1_cef"                 , "", 'F');
  factory->AddVariable("ele1_pt"               , "ele1_pt"                 , "", 'F');
  factory->AddVariable("ele1_dEtaFj1"          , "ele1_dEtaFj1"            , "", 'F');
  factory->AddVariable("ele1_dPhiFj1"          , "ele1_dPhiFj1"            , "", 'F');
  factory->AddVariable("ele1_dxy"              , "ele1_dxy"                , "", 'F');
  factory->AddVariable("ele1_dz"               , "ele1_dz"                 , "", 'F');
  factory->AddVariable("ele1_hfMva"            , "ele1_hfMva"              , "", 'F');
  factory->AddVariable("ele1_sv1_pt"           , "ele1_sv1_pt"             , "", 'F'); 
  factory->AddVariable("ele1_sv1_dEta"         , "ele1_sv1_dEta"           , "", 'F');  
  factory->AddVariable("ele1_sv1_dPhi"         , "ele1_sv1_dPhi"           , "", 'F');  
  factory->AddVariable("ele1_sv1_m"            , "ele1_sv1_m"              , "", 'F'); 
  factory->AddVariable("ele1_sv1_vtx3DVal"     , "ele1_sv1_vtx3DVal"       , "", 'F'); 
  factory->AddVariable("ele1_sv1_vtx3DeVal"    , "ele1_sv1_vtx3DeVal"      , "", 'F'); 
  factory->AddVariable("ele1_sv1_ntrk"         , "ele1_sv1_ntrk"           , "", 'I'); 
  factory->AddVariable("ele1_sv1_ndof"         , "ele1_sv1_ndof"           , "", 'F');  
  factory->AddVariable("ele1_sv1_chi2"         , "ele1_sv1_chi2"           , "", 'F');  
  factory->AddVariable("ele1_sv1_significance" , "ele1_sv1_significance"   , "", 'F');  
  factory->AddVariable("mu1_pt"                , "mu1_pt"                  , "", 'F');
  factory->AddVariable("mu1_dEtaFj1"           , "mu1_dEtaFj1"             , "", 'F');
  factory->AddVariable("mu1_dPhiFj1"           , "mu1_dPhiFj1"             , "", 'F');
  factory->AddVariable("mu1_dxy"               , "mu1_dxy"                 , "", 'F');
  factory->AddVariable("mu1_dz"                , "mu1_dz"                  , "", 'F');
  factory->AddVariable("mu1_sv1_pt"            , "mu1_sv1_pt"              , "", 'F'); 
  factory->AddVariable("mu1_sv1_dEta"          , "mu1_sv1_dEta"            , "", 'F');  
  factory->AddVariable("mu1_sv1_dPhi"          , "mu1_sv1_dPhi"            , "", 'F');  
  factory->AddVariable("mu1_sv1_m"             , "mu1_sv1_m"               , "", 'F'); 
  factory->AddVariable("mu1_sv1_vtx3DVal"      , "mu1_sv1_vtx3DVal"        , "", 'F'); 
  factory->AddVariable("mu1_sv1_vtx3DeVal"     , "mu1_sv1_vtx3DeVal"       , "", 'F'); 
  factory->AddVariable("mu1_sv1_ntrk"          , "mu1_sv1_ntrk"            , "", 'I'); 
  factory->AddVariable("mu1_sv1_ndof"          , "mu1_sv1_ndof"            , "", 'F');  
  factory->AddVariable("mu1_sv1_chi2"          , "mu1_sv1_chi2"            , "", 'F');  
  factory->AddVariable("mu1_sv1_significance"  , "mu1_sv1_significance"    , "", 'F');  
  factory->AddVariable("fj1_sv1_pt"            , "fj1_sv1_pt"              , "", 'F'); 
  factory->AddVariable("fj1_sv1_dEta"          , "fj1_sv1_dEta"            , "", 'F');  
  factory->AddVariable("fj1_sv1_dPhi"          , "fj1_sv1_dPhi"            , "", 'F');  
  factory->AddVariable("fj1_sv1_m"             , "fj1_sv1_m"               , "", 'F'); 
  factory->AddVariable("fj1_sv1_vtx3DVal"      , "fj1_sv1_vtx3DVal"        , "", 'F'); 
  factory->AddVariable("fj1_sv1_vtx3DeVal"     , "fj1_sv1_vtx3DeVal"       , "", 'F'); 
  factory->AddVariable("fj1_sv1_ntrk"          , "fj1_sv1_ntrk"            , "", 'I'); 
  factory->AddVariable("fj1_sv1_ndof"          , "fj1_sv1_ndof"            , "", 'F');  
  factory->AddVariable("fj1_sv1_chi2"          , "fj1_sv1_chi2"            , "", 'F');  
  factory->AddVariable("fj1_sv1_significance"  , "fj1_sv1_significance"    , "", 'F');  
  factory->AddVariable("fj1_sv2_pt"            , "fj1_sv2_pt"              , "", 'F'); 
  factory->AddVariable("fj1_sv2_dEta"          , "fj1_sv2_dEta"            , "", 'F');  
  factory->AddVariable("fj1_sv2_dPhi"          , "fj1_sv2_dPhi"            , "", 'F');  
  factory->AddVariable("fj1_sv2_m"             , "fj1_sv2_m"               , "", 'F'); 
  factory->AddVariable("fj1_sv2_vtx3DVal"      , "fj1_sv2_vtx3DVal"        , "", 'F'); 
  factory->AddVariable("fj1_sv2_vtx3DeVal"     , "fj1_sv2_vtx3DeVal"       , "", 'F'); 
  factory->AddVariable("fj1_sv2_ntrk"          , "fj1_sv2_ntrk"            , "", 'I'); 
  factory->AddVariable("fj1_sv2_ndof"          , "fj1_sv2_ndof"            , "", 'F');  
  factory->AddVariable("fj1_sv2_chi2"          , "fj1_sv2_chi2"            , "", 'F');  
  factory->AddVariable("fj1_sv2_significance"  , "fj1_sv2_significance"    , "", 'F');  
  // Add the variable carrying the regression target
  //factory->AddTarget("sumNu_pt"              , "sumNu_pt"                , "", 'F');
  //factory->AddTarget("sumNu_dEtaFj1"         , "sumNu_dEtaFj1"           , "", 'F');
  //factory->AddTarget("sumNu_dPhiFj1"         , "sumNu_dPhiFj1"           , "", 'F');
  factory->AddTarget("sumNu_E"               , "sumNu_E"                 , "", 'F');
  // It is also possible to declare additional targets for multi-dimensional regression, ie:
  //     factory->AddTarget( "fvalue2" );
  // BUT: this is currently ONLY implemented for MLP
  // Read training and test data (see TMVAClassification for reading ASCII files)
  // load the signal and background event samples from ROOT trees
  factory->SetWeightExpression( Form("weight/%.2f",sum_weights), "Regression" );
  TCut preselectionCut = "";
  TString prepareOptions="NormMode=NumEvents";
    prepareOptions+=":SplitMode=Random"; 
    prepareOptions+=":MixMode=Random";
  //  prepareOptions+=":nTrain_Regression=5000:nTest_Regression=5000";

  factory->PrepareTrainingAndTestTree(preselectionCut, prepareOptions);
  //Neural network (MLP)
  //factory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );
  //factory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:BoostType=Bagging:SeparationType=RegressionVariance:nCuts=10000:PruneMethod=CostComplexity:PruneStrength=20:MaxDepth=15:NTrees=100" );
  factory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:BoostType=Grad:Shrinkage=0.2:SeparationType=RegressionVariance:nCuts=1000:MaxDepth=15:MinNodeSize=3%:NTrees=400" );
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  outputFile->Close();
}
