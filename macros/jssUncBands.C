#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <vector>
#include <map>
#include "vhbbPlot.h"

vector<TString> histoVariations = {"Nominal", "FSRFactUp", "FSRFactDown", "FSRRenormUp", "FSRRenormDown", "ISRFactUp", "ISRFactDown", "ISRRenormUp", "ISRRenormDown", "MPIUp", "MPIDown"};
std::map<TString,int> varColors = {
  {"Nominal"       , kBlack     },     
  {"FSRFactUp"     , kRed       },    
  {"FSRFactDown"   , kRed       },    
  {"FSRRenormUp"   , kOrange-2  },    
  {"FSRRenormDown" , kOrange-2  },    
  {"ISRFactUp"     , kBlue      },    
  {"ISRFactDown"   , kBlue      },    
  {"ISRRenormUp"   , kMagenta   },    
  {"ISRRenormDown" , kMagenta   },    
  {"MPIUp"         , kGreen-3   },    
  {"MPIDown"       , kGreen-3   }
};

void jssUncBands(TString mode="WH", TString var="fj1Tau21", bool nbGt0=false) {
  int nbins; float xmin,xmax; TString axisLabel;
  if     (var=="fj1Tau21"       ) { nbins=20; xmin=  0; xmax=  1; axisLabel="#tau_{2}/#tau_{1}"     ;}
  else if(var=="fj1Tau32"       ) { nbins=20; xmin=  0; xmax=  1; axisLabel="#tau_{3}/#tau_{2}"     ;}
  else if(var=="psi020503010502") { nbins=20; xmin=  0; xmax=0.5; axisLabel="#psi(2,0.5,3,1,0.5,2)" ;}
  else if(var=="psi022003012002") { nbins=20; xmin=  0; xmax=0.4; axisLabel="#psi(2,2.0,3,1,2.0,2)" ;}
  else if(var=="psi032003012002") { nbins=20; xmin=  0; xmax=  4; axisLabel="#psi(3,2.0,3,1,2.0,2)" ;}
  else if(var=="psi031003012002") { nbins=20; xmin=  0; xmax=0.4; axisLabel="#psi(3,1.0,3,1,2.0,2)" ;}
  else if(var=="psi030503010502") { nbins=20; xmin=  0; xmax=  3; axisLabel="#psi(3,0.5,3,1,0.5,2)" ;}
  else { printf("Variable not supported yet...\n"); return; }
  
  TString dir="/mnt/hadoop/scratch/dhsu/pythiaVariations/";
  vector<TString> files;  
  files.push_back("WJetsToLNu_HT-100To200_FSRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-100To200_FSRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-100To200_FSRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-100To200_FSRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-100To200_ISRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-100To200_ISRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-100To200_ISRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-100To200_ISRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-100To200_MPIDown.root"         );  
  files.push_back("WJetsToLNu_HT-100To200_MPIUp.root"           );  
  files.push_back("WJetsToLNu_HT-100To200_Nominal.root"         );  
  files.push_back("WJetsToLNu_HT-1200To2500_FSRFactDown.root"   );  
  files.push_back("WJetsToLNu_HT-1200To2500_FSRFactUp.root"     );  
  files.push_back("WJetsToLNu_HT-1200To2500_FSRRenormDown.root" );  
  files.push_back("WJetsToLNu_HT-1200To2500_FSRRenormUp.root"   );  
  files.push_back("WJetsToLNu_HT-1200To2500_ISRFactDown.root"   );  
  files.push_back("WJetsToLNu_HT-1200To2500_ISRFactUp.root"     );  
  files.push_back("WJetsToLNu_HT-1200To2500_ISRRenormDown.root" );  
  files.push_back("WJetsToLNu_HT-1200To2500_ISRRenormUp.root"   );  
  files.push_back("WJetsToLNu_HT-1200To2500_MPIDown.root"       );  
  files.push_back("WJetsToLNu_HT-1200To2500_MPIUp.root"         );  
  files.push_back("WJetsToLNu_HT-1200To2500_Nominal.root"       );  
  files.push_back("WJetsToLNu_HT-200To400_FSRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-200To400_FSRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-200To400_FSRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-200To400_FSRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-200To400_ISRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-200To400_ISRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-200To400_ISRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-200To400_ISRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-200To400_MPIDown.root"         );  
  files.push_back("WJetsToLNu_HT-200To400_MPIUp.root"           );  
  files.push_back("WJetsToLNu_HT-200To400_Nominal.root"         );  
  files.push_back("WJetsToLNu_HT-2500ToInf_FSRFactDown.root"    );  
  files.push_back("WJetsToLNu_HT-2500ToInf_FSRFactUp.root"      );  
  files.push_back("WJetsToLNu_HT-2500ToInf_FSRRenormDown.root"  );  
  files.push_back("WJetsToLNu_HT-2500ToInf_FSRRenormUp.root"    );  
  files.push_back("WJetsToLNu_HT-2500ToInf_ISRFactDown.root"    );  
  files.push_back("WJetsToLNu_HT-2500ToInf_ISRFactUp.root"      );  
  files.push_back("WJetsToLNu_HT-2500ToInf_ISRRenormDown.root"  );  
  files.push_back("WJetsToLNu_HT-2500ToInf_ISRRenormUp.root"    );  
  files.push_back("WJetsToLNu_HT-2500ToInf_MPIDown.root"        );  
  files.push_back("WJetsToLNu_HT-2500ToInf_MPIUp.root"          );  
  files.push_back("WJetsToLNu_HT-2500ToInf_Nominal.root"        );  
  files.push_back("WJetsToLNu_HT-400To600_FSRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-400To600_FSRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-400To600_FSRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-400To600_FSRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-400To600_ISRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-400To600_ISRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-400To600_ISRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-400To600_ISRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-400To600_MPIDown.root"         );  
  files.push_back("WJetsToLNu_HT-400To600_MPIUp.root"           );  
  files.push_back("WJetsToLNu_HT-400To600_Nominal.root"         );  
  files.push_back("WJetsToLNu_HT-600To800_FSRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-600To800_FSRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-600To800_FSRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-600To800_FSRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-600To800_ISRFactDown.root"     );  
  files.push_back("WJetsToLNu_HT-600To800_ISRFactUp.root"       );  
  files.push_back("WJetsToLNu_HT-600To800_ISRRenormDown.root"   );  
  files.push_back("WJetsToLNu_HT-600To800_ISRRenormUp.root"     );  
  files.push_back("WJetsToLNu_HT-600To800_MPIDown.root"         );  
  files.push_back("WJetsToLNu_HT-600To800_MPIUp.root"           );  
  files.push_back("WJetsToLNu_HT-600To800_Nominal.root"         );  
  files.push_back("WJetsToLNu_HT-800To1200_FSRFactDown.root"    );  
  files.push_back("WJetsToLNu_HT-800To1200_FSRFactUp.root"      );  
  files.push_back("WJetsToLNu_HT-800To1200_FSRRenormDown.root"  );  
  files.push_back("WJetsToLNu_HT-800To1200_FSRRenormUp.root"    );  
  files.push_back("WJetsToLNu_HT-800To1200_ISRFactDown.root"    );  
  files.push_back("WJetsToLNu_HT-800To1200_ISRFactUp.root"      );  
  files.push_back("WJetsToLNu_HT-800To1200_ISRRenormDown.root"  );  
  files.push_back("WJetsToLNu_HT-800To1200_ISRRenormUp.root"    );  
  files.push_back("WJetsToLNu_HT-800To1200_MPIDown.root"        );  
  files.push_back("WJetsToLNu_HT-800To1200_MPIUp.root"          );  
  files.push_back("WJetsToLNu_HT-800To1200_Nominal.root"        );  
  files.push_back("WminusH_HToBB_WToLNu_M125_FSRFactDown.root"  );  
  files.push_back("WminusH_HToBB_WToLNu_M125_FSRFactUp.root"    );  
  files.push_back("WminusH_HToBB_WToLNu_M125_FSRRenormDown.root");  
  files.push_back("WminusH_HToBB_WToLNu_M125_FSRRenormUp.root"  );  
  files.push_back("WminusH_HToBB_WToLNu_M125_ISRFactDown.root"  );  
  files.push_back("WminusH_HToBB_WToLNu_M125_ISRFactUp.root"    );  
  files.push_back("WminusH_HToBB_WToLNu_M125_ISRRenormDown.root");  
  files.push_back("WminusH_HToBB_WToLNu_M125_ISRRenormUp.root"  );  
  files.push_back("WminusH_HToBB_WToLNu_M125_MPIDown.root"      );  
  files.push_back("WminusH_HToBB_WToLNu_M125_MPIUp.root"        );  
  files.push_back("WminusH_HToBB_WToLNu_M125_Nominal.root"      );  
  files.push_back("WplusH_HToBB_WToLNu_M125_FSRFactDown.root"   );  
  files.push_back("WplusH_HToBB_WToLNu_M125_FSRFactUp.root"     );  
  files.push_back("WplusH_HToBB_WToLNu_M125_FSRRenormDown.root" );  
  files.push_back("WplusH_HToBB_WToLNu_M125_FSRRenormUp.root"   );  
  files.push_back("WplusH_HToBB_WToLNu_M125_ISRFactDown.root"   );  
  files.push_back("WplusH_HToBB_WToLNu_M125_ISRFactUp.root"     );  
  files.push_back("WplusH_HToBB_WToLNu_M125_ISRRenormDown.root" );  
  files.push_back("WplusH_HToBB_WToLNu_M125_ISRRenormUp.root"   );  
  files.push_back("WplusH_HToBB_WToLNu_M125_MPIDown.root"       );  
  files.push_back("WplusH_HToBB_WToLNu_M125_MPIUp.root"         );  
  files.push_back("WplusH_HToBB_WToLNu_M125_Nominal.root"       );  
  
  // Declare variables
  int nB;
  float normalizedWeight,genFatJetPt;
  float fj1Tau32,fj1Tau21,fj1Tau32SD,fj1Tau21SD;
  float fj1ECFN_2_3_05,fj1ECFN_1_2_05,fj1ECFN_2_3_20,fj1ECFN_1_2_20,fj1ECFN_3_3_20,fj1ECFN_3_3_10,fj1ECFN_3_3_05;
  std::map<TString, float> psi;
  float plotVar;

  TH1F *histos[11], *ratios[11];
  //TH1F *hNominal, *hFSRFactUp, *hFSRFactDown, *hFSRRenormUp, *hFSRRenormDown, *hISRFactUp, *hISRFactDown, *hISRRenormUp, *hISRRenormDown, *hMPIUp, *hMPIDown;
  for(unsigned iV=0; iV<11; iV++) {
    histos[iV] = new TH1F(TString("h"+histoVariations[iV]), histoVariations[iV], nbins, xmin, xmax); histos[iV]->SetDirectory(0); histos[iV]->Sumw2(); 
    ratios[iV] = new TH1F(TString("r"+histoVariations[iV]), histoVariations[iV], nbins, xmin, xmax); ratios[iV]->SetDirectory(0); ratios[iV]->Sumw2(); 
  }
  for(unsigned iFile=0; iFile<files.size(); iFile++) {
    if(mode=="WH" && !(files[iFile].Contains("WplusH")||files[iFile].Contains("WminusH"))) continue;
    if(mode=="Wjets" && !files[iFile].Contains("WJets")) continue;
    
    string inputFileNameStr(files[iFile]);
    /*float theQCDKFactor=1;
    if(mode=="Wjets") {
      // Determine the HT bin
      // Assume the sample names are of the form "WplusH_HToBB_WToLNu_M125_Nominal.root"
      size_t htWord    = inputFileNameStr.find("HT-");
      size_t toWord    = inputFileNameStr.find("To");
      size_t lastUnderscore   = inputFileNameStr.find_last_of("_");
      string htLowStr   = inputFileNameStr.substr(htWord+3, toWord-htWord-3);
      string htHighStr  = inputFileNameStr.substr(toWord+2, lastUnderscore-toWord-2);
      htLow   = atoi(htLowStr.c_str());
      if(htHighStr=="inf"||htHighStr=="Inf")
        htHigh=99999;
      else
        htHigh  = atoi(htHighStr.c_str());
      if(htLow==0 || htHigh==0) {
        throw std::runtime_error(Form("Warning: Error parsing the filename \"%s\", probably it is not of the form \"WplusH_HToBB_WToLNu_M125_Nominal.root\" (go fix that)", files[iFile].Data()));
        return false;
      }
      theQCDKFactor = qcdKFactor(vhbbPlot::kWjets, htLow, htHigh);
    }*/

    // Find the variation in the list
    bool foundVariationInList=false; unsigned varIdx;
    {
      size_t lastUnderscore=inputFileNameStr.find_last_of("_");
      size_t lastDot   = inputFileNameStr.find_last_of(".");
      TString theVariation(inputFileNameStr.substr(lastUnderscore+1, lastDot-lastUnderscore-1).c_str());
      for(int iV=0; iV<11 && !foundVariationInList; iV++) 
        if(histoVariations[iV]==theVariation) { foundVariationInList=true; varIdx = iV; }
    } assert(foundVariationInList);

    TFile *file = TFile::Open(dir+files[iFile], "read"); assert(file);
    TTree *events = (TTree*)file->Get("events"); assert(events);
    events->SetBranchStatus("*",0);
    events->SetBranchStatus("nB"              , 1); events->SetBranchAddress("nB", &nB);
    events->SetBranchStatus("normalizedWeight", 1); events->SetBranchAddress("normalizedWeight", &normalizedWeight);
    events->SetBranchStatus("genFatJetPt"     , 1); events->SetBranchAddress("genFatJetPt", &genFatJetPt);
    events->SetBranchStatus("fj1Tau21"        , 1); events->SetBranchAddress("fj1Tau21", &fj1Tau21);
    events->SetBranchStatus("fj1Tau32"        , 1); events->SetBranchAddress("fj1Tau32", &fj1Tau32);
    events->SetBranchStatus("fj1ECFN_2_3_05"  , 1); events->SetBranchAddress("fj1ECFN_2_3_05", &fj1ECFN_2_3_05);  
    events->SetBranchStatus("fj1ECFN_1_2_05"  , 1); events->SetBranchAddress("fj1ECFN_1_2_05", &fj1ECFN_1_2_05);  
    events->SetBranchStatus("fj1ECFN_2_3_20"  , 1); events->SetBranchAddress("fj1ECFN_2_3_20", &fj1ECFN_2_3_20);  
    events->SetBranchStatus("fj1ECFN_1_2_20"  , 1); events->SetBranchAddress("fj1ECFN_1_2_20", &fj1ECFN_1_2_20);  
    events->SetBranchStatus("fj1ECFN_3_3_20"  , 1); events->SetBranchAddress("fj1ECFN_3_3_20", &fj1ECFN_3_3_20);  
    events->SetBranchStatus("fj1ECFN_3_3_10"  , 1); events->SetBranchAddress("fj1ECFN_3_3_10", &fj1ECFN_3_3_10);  
    events->SetBranchStatus("fj1ECFN_3_3_05"  , 1); events->SetBranchAddress("fj1ECFN_3_3_05", &fj1ECFN_3_3_05);  
    Long64_t nEntries = events->GetEntries();
    float varMax=histos[varIdx]->GetBinLowEdge(histos[varIdx]->GetNbinsX()+1)-0.0001;
    for(Long64_t iEntry=0; iEntry<nEntries; iEntry++) {
      events->GetEntry(iEntry);
      if(genFatJetPt<250) continue;
      if(nbGt0 && nB==0) continue;
      //psi["021004010502"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),4.00); 
      //psi["012004010502"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),4.00); 
      //psi["021003011002"]=fj1ECFN_2_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),2.00); 
      //psi["022004011002"]=fj1ECFN_2_4_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),4.00); 
      psi["020503010502"]=fj1ECFN_2_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),2.00); 
      psi["022003012002"]=fj1ECFN_2_3_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),2.00); 
      //psi["012004020503"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_2_3_05),2.00); 
      //psi["021004020503"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_2_3_05),2.00); 
      psi["032003012002"]=fj1ECFN_3_3_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),3.00); 
      //psi["012003011002"]=fj1ECFN_1_3_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),2.00); 
      //psi["011003010502"]=fj1ECFN_1_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),2.00); 
      //psi["031003011002"]=fj1ECFN_3_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),3.00); 
      psi["031003012002"]=fj1ECFN_3_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),1.50); 
      //psi["030503011002"]=fj1ECFN_3_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),1.50); 
      //psi["022004012002"]=fj1ECFN_2_4_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),2.00); 
      //psi["021004011002"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),2.00); 
      //psi["012004011002"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),2.00); 
      //psi["022004021003"]=fj1ECFN_2_4_20/pow(TMath::Max((float)0.,fj1ECFN_2_3_10),2.00); 
      //psi["012003010502"]=fj1ECFN_1_3_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),4.00); 
      //psi["022003011002"]=fj1ECFN_2_3_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),4.00); 
      //psi["022004031003"]=fj1ECFN_2_4_20/pow(TMath::Max((float)0.,fj1ECFN_3_3_10),1.33); 
      //psi["012004030503"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_3_3_05),1.33); 
      //psi["021004030503"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_3_3_05),1.33); 
      //psi["020504010502"]=fj1ECFN_2_4_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),2.00); 
      //psi["022004030503"]=fj1ECFN_2_4_20/pow(TMath::Max((float)0.,fj1ECFN_3_3_05),2.67); 
      //psi["011004010502"]=fj1ECFN_1_4_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),2.00); 
      psi["030503010502"]=fj1ECFN_3_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),3.00); 
      //psi["021003010502"]=fj1ECFN_2_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),4.00); 
      //psi["012004011003"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_1_3_10),2.00); 
      //psi["012004021004"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_2_4_10),1.00); 
      //psi["021004012004"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_1_4_20),1.00); 
      //psi["011003020503"]=fj1ECFN_1_3_10/pow(TMath::Max((float)0.,fj1ECFN_2_3_05),1.00); 
      //psi["020503011003"]=fj1ECFN_2_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_3_10),1.00); 
      //psi["012003030503"]=fj1ECFN_1_3_20/pow(TMath::Max((float)0.,fj1ECFN_3_3_05),1.33); 
      //psi["030503012003"]=fj1ECFN_3_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_3_20),(float)0.75); 
      //psi["010503010502"]=fj1ECFN_1_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),1.00); 
      //psi["020503011002"]=fj1ECFN_2_3_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),1.00); 
      //psi["011003011002"]=fj1ECFN_1_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),1.00); 
      //psi["020504011002"]=fj1ECFN_2_4_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),1.00); 
      //psi["012004021003"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_2_3_10),1.00); 
      //psi["012004012002"]=fj1ECFN_1_4_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),1.00); 
      //psi["020504020503"]=fj1ECFN_2_4_05/pow(TMath::Max((float)0.,fj1ECFN_2_3_05),1.00); 
      //psi["011004011002"]=fj1ECFN_1_4_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_10),1.00); 
      //psi["021004012002"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),1.00); 
      //psi["021003012002"]=fj1ECFN_2_3_10/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),1.00); 
      //psi["011004020503"]=fj1ECFN_1_4_10/pow(TMath::Max((float)0.,fj1ECFN_2_3_05),1.00); 
      //psi["010504010502"]=fj1ECFN_1_4_05/pow(TMath::Max((float)0.,fj1ECFN_1_2_05),1.00); 
      //psi["021004021003"]=fj1ECFN_2_4_10/pow(TMath::Max((float)0.,fj1ECFN_2_3_10),1.00); 
      //psi["012003012002"]=fj1ECFN_1_3_20/pow(TMath::Max((float)0.,fj1ECFN_1_2_20),1.00); 
      //psi["022004022003"]=fj1ECFN_2_4_20/pow(TMath::Max((float)0.,fj1ECFN_2_3_20),1.00); 
      if     (var=="fj1Tau21"       ) plotVar = fj1Tau21;                
      else if(var=="fj1Tau32"       ) plotVar = fj1Tau32;                
      else if(var=="psi020503010502") plotVar = psi["020503010502"];
      else if(var=="psi022003012002") plotVar = psi["022003012002"];
      else if(var=="psi032003012002") plotVar = psi["032003012002"];
      else if(var=="psi031003012002") plotVar = psi["031003012002"];
      else if(var=="psi030503010502") plotVar = psi["030503010502"];
      histos[varIdx]->Fill(
        TMath::Min(varMax, plotVar),
        normalizedWeight
      );
    }
    file->Close();
  }
  histos[0]->Scale(1./histos[0]->Integral());
  for(unsigned iV=0; iV<11; iV++) {
    if(iV!=0) histos[iV]->Scale(histos[0]->Integral()/histos[iV]->Integral());
    histos[iV]->SetLineColor(varColors[histoVariations[iV]]);
    histos[iV]->SetLineWidth(2);
    ratios[iV]->SetLineColor(varColors[histoVariations[iV]]);
    ratios[iV]->SetLineWidth(2);
    ratios[iV]->SetFillStyle(0);
    for(int nb=1;nb<=histos[0]->GetNbinsX(); nb++) {
      float yield=histos[0]->GetBinContent(nb);
      if(yield>0) ratios[iV]->SetBinContent(nb, histos[iV]->GetBinContent(nb)/yield);
      else ratios[iV]->SetBinContent(nb, 1);
    }
  }
  // Calculate rms band
  TH1F *rmsBand = (TH1F*)histos[0]->Clone("rmsBand");
  TH1F *rmsRatioBand = (TH1F*)ratios[0]->Clone("rmsRatioBand");
  for(int nb=1;nb<=histos[0]->GetNbinsX(); nb++) {
    float N=0; float sum=0,sum2=0,sumR=0,sumR2=0;
    for(unsigned iV=0; iV<11; iV++) {
      sum+=histos[iV]->GetBinContent(nb);
      sumR+=ratios[iV]->GetBinContent(nb);
      N++;
    }
    float mean=sum/N;
    float meanR=sumR/N;
    for(unsigned iV=0; iV<11; iV++) {
      sum2+=pow(histos[iV]->GetBinContent(nb)-mean,2);
      sumR2+=pow(ratios[iV]->GetBinContent(nb)-meanR,2);
    }
    float rms=sqrt(sum2/N);
    float rmsR=sqrt(sumR2/N);
    rmsBand->SetBinContent(nb,mean);
    rmsBand->SetBinError(nb,rms);
    ratios[0]->SetBinContent(nb,meanR);
    rmsRatioBand->SetBinError(nb,rmsR);
  }
  // Start drawing stuff
  gStyle->SetOptStat(0);
  TCanvas *canvas=new TCanvas("canvas","canvas",600,480);
  canvas->SetTopMargin(0.0); 
  canvas->SetBottomMargin(0); 
  canvas->SetRightMargin(0.02);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.,1);
  pad1->SetTopMargin(0.1);
  pad1->SetLeftMargin(0.15);
  pad1->SetRightMargin(0.04);
  pad1->SetBottomMargin(0.03); 
  //pad1->SetGridx();         
  pad1->Draw();             
  pad1->cd();      
  rmsBand->SetFillStyle(3254);
  rmsBand->SetFillColor(kGray+3);
  rmsBand->SetTitle("");
  rmsBand->GetXaxis()->SetLabelSize(0);
  rmsBand->SetMinimum(0);
  rmsBand->SetMaximum(rmsBand->GetMaximum()*1.5);
  rmsBand->Draw("e2");
  histos[0]->Draw("hist e0 same");
  TLegend *legend = new TLegend(.2,.5,.5,.88);
  legend->AddEntry(histos[0], "Nominal", "lp");
  legend->AddEntry(histos[1], "FSR fact. scale", "l");
  legend->AddEntry(histos[3], "FSR renorm. scale", "l");
  legend->AddEntry(histos[5], "ISR fact. scale", "l");
  legend->AddEntry(histos[7], "ISR renorm. scale", "l");
  legend->AddEntry(histos[9], "UE MPIs", "l");
  legend->AddEntry(rmsBand, "RMS band", "f");
  legend->Draw("same");
  canvas->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.15);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin(0.04);
  //pad2->SetGridx(); 
  pad2->Draw();
  pad2->cd();

  rmsRatioBand->SetTitle("");
  rmsRatioBand->GetXaxis()->SetTitleSize(0.18);
  rmsRatioBand->GetXaxis()->SetTitleOffset(1.1);
  rmsRatioBand->GetYaxis()->SetNdivisions(503);
  rmsRatioBand->GetXaxis()->SetLabelSize(0.2);
  rmsRatioBand->GetXaxis()->SetTitle(axisLabel);
  rmsRatioBand->GetYaxis()->SetTitle("Var./Nominal");
  rmsRatioBand->GetYaxis()->CenterTitle();
  rmsRatioBand->GetYaxis()->SetTitleOffset(0.37);
  rmsRatioBand->GetYaxis()->SetTitleSize(0.12);
  rmsRatioBand->GetYaxis()->SetLabelSize(0.15);
  rmsRatioBand->SetFillStyle(3254);
  rmsRatioBand->SetFillColor(kGray+3);
  rmsRatioBand->SetMinimum(.32);
  rmsRatioBand->SetMaximum(1.68);
  rmsRatioBand->Draw("e2");

  for(unsigned iV=1; iV<11; iV++)
    ratios[iV]->Draw("hist same");
  rmsRatioBand->Draw("e2 same");
  ratios[0]->Draw("hist same");
  
}
