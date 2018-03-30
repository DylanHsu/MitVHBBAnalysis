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
  factory->AddVariable("fatjet_pt"                 ,"fatjet_pt"                 , "", 'F');
  factory->AddVariable("fatjet_eta"                ,"fatjet_eta"                , "", 'F');
  factory->AddVariable("fatjet_phi"                ,"fatjet_phi"                , "", 'F');
  factory->AddVariable("fatjet_nhf"                ,"fatjet_nhf"                , "", 'F');
  factory->AddVariable("fatjet_chf"                ,"fatjet_chf"                , "", 'F');
  factory->AddVariable("fatjet_nef"                ,"fatjet_nef"                , "", 'F');
  //factory->AddVariable("fatjet_cef"                ,"fatjet_cef"                , "", 'F');
  //factory->AddVariable("ele1_pt"                   ,"ele1_pt"                   , "", 'F');
  //factory->AddVariable("ele1_eta"                  ,"ele1_eta"                  , "", 'F');
  //factory->AddVariable("ele1_phi"                  ,"ele1_phi"                  , "", 'F');
  //factory->AddVariable("ele1_chIso"                ,"ele1_chIso"                , "", 'F');
  //factory->AddVariable("ele1_nhIso"                ,"ele1_nhIso"                , "", 'F');
  //factory->AddVariable("ele1_phIso"                ,"ele1_phIso"                , "", 'F');
  //factory->AddVariable("ele1_dxy"                  ,"ele1_dxy"                  , "", 'F');
  //factory->AddVariable("ele1_dz"                   ,"ele1_dz"                   , "", 'F');
  //factory->AddVariable("ele1_sieie"                ,"ele1_sieie"                , "", 'F');
  //factory->AddVariable("ele1_sipip"                ,"ele1_sipip"                , "", 'F');
  //factory->AddVariable("ele1_r9"                   ,"ele1_r9"                   , "", 'F');
  //factory->AddVariable("ele1_dEtaInSeed"           ,"ele1_dEtaInSeed"           , "", 'F');
  //factory->AddVariable("ele1_dPhiIn"               ,"ele1_dPhiIn"               , "", 'F');
  //factory->AddVariable("ele1_eseed"                ,"ele1_eseed"                , "", 'F');
  //factory->AddVariable("ele1_hOverE"               ,"ele1_hOverE"               , "", 'F');
  //factory->AddVariable("ele1_ecalE"                ,"ele1_ecalE"                , "", 'F');
  //factory->AddVariable("ele1_trackP"               ,"ele1_trackP"               , "", 'F');
  //factory->AddVariable("ele2_pt"                   ,"ele2_pt"                   , "", 'F');
  //factory->AddVariable("ele2_eta"                  ,"ele2_eta"                  , "", 'F');
  //factory->AddVariable("ele2_phi"                  ,"ele2_phi"                  , "", 'F');
  //factory->AddVariable("ele2_chIso"                ,"ele2_chIso"                , "", 'F');
  //factory->AddVariable("ele2_nhIso"                ,"ele2_nhIso"                , "", 'F');
  //factory->AddVariable("ele2_phIso"                ,"ele2_phIso"                , "", 'F');
  //factory->AddVariable("ele2_dxy"                  ,"ele2_dxy"                  , "", 'F');
  //factory->AddVariable("ele2_dz"                   ,"ele2_dz"                   , "", 'F');
  //factory->AddVariable("ele2_sieie"                ,"ele2_sieie"                , "", 'F');
  //factory->AddVariable("ele2_sipip"                ,"ele2_sipip"                , "", 'F');
  //factory->AddVariable("ele2_r9"                   ,"ele2_r9"                   , "", 'F');
  //factory->AddVariable("ele2_dEtaInSeed"           ,"ele2_dEtaInSeed"           , "", 'F');
  //factory->AddVariable("ele2_dPhiIn"               ,"ele2_dPhiIn"               , "", 'F');
  //factory->AddVariable("ele2_eseed"                ,"ele2_eseed"                , "", 'F');
  //factory->AddVariable("ele2_hOverE"               ,"ele2_hOverE"               , "", 'F');
  //factory->AddVariable("ele2_ecalE"                ,"ele2_ecalE"                , "", 'F');
  //factory->AddVariable("ele2_trackP"               ,"ele2_trackP"               , "", 'F');
  factory->AddVariable("mu1_pt"                    ,"mu1_pt"                    , "", 'F');
  factory->AddVariable("mu1_eta"                   ,"mu1_eta"                   , "", 'F');
  factory->AddVariable("mu1_phi"                   ,"mu1_phi"                   , "", 'F');
  factory->AddVariable("mu1_chIso"                 ,"mu1_chIso"                 , "", 'F');
  factory->AddVariable("mu1_nhIso"                 ,"mu1_nhIso"                 , "", 'F');
  factory->AddVariable("mu1_phIso"                 ,"mu1_phIso"                 , "", 'F');
  factory->AddVariable("mu1_dxy"                   ,"mu1_dxy"                   , "", 'F');
  factory->AddVariable("mu1_dz"                    ,"mu1_dz"                    , "", 'F');
  factory->AddVariable("mu1_validFraction"         ,"mu1_validFraction"         , "", 'F');
  factory->AddVariable("mu1_nValidMuon"            ,"mu1_nValidMuon"            , "", 'I');
  factory->AddVariable("mu1_nValidPixel"           ,"mu1_nValidPixel"           , "", 'I');
  factory->AddVariable("mu1_trkLayersWithMmt"      ,"mu1_trkLayersWithMmt"      , "", 'I');
  factory->AddVariable("mu1_pixLayersWithMmt"      ,"mu1_pixLayersWithMmt"      , "", 'I');
  factory->AddVariable("mu1_nMatched"              ,"mu1_nMatched"              , "", 'I');
  //factory->AddVariable("mu1_normChi2"              ,"mu1_normChi2"              , "", 'F');
  factory->AddVariable("mu1_chi2LocalPosition"     ,"mu1_chi2LocalPosition"     , "", 'I');
  factory->AddVariable("mu1_trkKink"               ,"mu1_trkKink"               , "", 'I');
  factory->AddVariable("mu1_segmentCompatibility"  ,"mu1_segmentCompatibility"  , "", 'F');
  //factory->AddVariable("mu2_pt"                    ,"mu2_pt"                    , "", 'F');
  //factory->AddVariable("mu2_eta"                   ,"mu2_eta"                   , "", 'F');
  //factory->AddVariable("mu2_phi"                   ,"mu2_phi"                   , "", 'F');
  //factory->AddVariable("mu2_chIso"                 ,"mu2_chIso"                 , "", 'F');
  //factory->AddVariable("mu2_nhIso"                 ,"mu2_nhIso"                 , "", 'F');
  //factory->AddVariable("mu2_phIso"                 ,"mu2_phIso"                 , "", 'F');
  //factory->AddVariable("mu2_dxy"                   ,"mu2_dxy"                   , "", 'F');
  //factory->AddVariable("mu2_dz"                    ,"mu2_dz"                    , "", 'F');
  //factory->AddVariable("mu2_validFraction"         ,"mu2_validFraction"         , "", 'F');
  //factory->AddVariable("mu2_nValidMuon"            ,"mu2_nValidMuon"            , "", 'I');
  //factory->AddVariable("mu2_nValidPixel"           ,"mu2_nValidPixel"           , "", 'I');
  //factory->AddVariable("mu2_trkLayersWithMmt"      ,"mu2_trkLayersWithMmt"      , "", 'I');
  //factory->AddVariable("mu2_pixLayersWithMmt"      ,"mu2_pixLayersWithMmt"      , "", 'I');
  //factory->AddVariable("mu2_nMatched"              ,"mu2_nMatched"              , "", 'I');
  //factory->AddVariable("mu2_normChi2"              ,"mu2_normChi2"              , "", 'F');
  //factory->AddVariable("mu2_chi2LocalPosition"     ,"mu2_chi2LocalPosition"     , "", 'I');
  //factory->AddVariable("mu2_trkKink"               ,"mu2_trkKink"               , "", 'I');
  //factory->AddVariable("mu2_segmentCompatibility"  ,"mu2_segmentCompatibility"  , "", 'F');
  factory->AddVariable("sv1_pt"                    ,"sv1_pt"                    , "", 'F');
  factory->AddVariable("sv1_m"                     ,"sv1_m"                     , "", 'F');
  factory->AddVariable("sv1_3Dval"                 ,"sv1_3Dval"                 , "", 'F');
  factory->AddVariable("sv1_3Derr"                 ,"sv1_3Derr"                 , "", 'F');
  factory->AddVariable("sv1_ntrk"                  ,"sv1_ntrk"                  , "", 'I');
  //factory->AddVariable("sv2_pt"                    ,"sv2_pt"                    , "", 'F');
  //factory->AddVariable("sv2_m"                     ,"sv2_m"                     , "", 'F');
  //factory->AddVariable("sv2_3Dval"                 ,"sv2_3Dval"                 , "", 'F');
  //factory->AddVariable("sv2_3Derr"                 ,"sv2_3Derr"                 , "", 'F');
  //factory->AddVariable("sv2_ntrk"                  ,"sv2_ntrk"                  , "", 'I');
  // Add the variable carrying the regression target
  //factory->AddTarget("sumNu_pt" );
  //factory->AddTarget("sumNu_eta");
  //factory->AddTarget("sumNu_phi");
  //factory->AddTarget("sumNu_E"  );
  factory->AddTarget("sumNu_E/fatjet_pt"  );
  // It is also possible to declare additional targets for multi-dimensional regression, ie:
  //     factory->AddTarget( "fvalue2" );
  // BUT: this is currently ONLY implemented for MLP
  // Read training and test data (see TMVAClassification for reading ASCII files)
  // load the signal and background event samples from ROOT trees
  factory->SetWeightExpression( Form("weight/%.2f",sum_weights), "Regression" );
  TCut preselectionCut = "fatjet_MSD>70&&sumNu_E>0&&sumNu_E<400";
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
