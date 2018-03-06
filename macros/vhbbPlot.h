#ifndef vhbbPlot_h
#define vhbbPlot_h

#include <TBranch.h>
#include <map>

namespace vhbbPlot {
  
  float theLumi=35900.;
  //const float bDiscrLoose = 0.5426, bDiscrMedium = 0.8484, bDiscrTight  = 0.9535; //csv
  const float bDiscrLoose = -0.5884, bDiscrMedium = 0.4432, bDiscrTight  = 0.9432; //cmva
  const float doubleBCut = 0.8; // double b tagger
  
  TString trainTreeEventSplitStr="(eventNumber % 10)<3";
  TString testTreeEventSplitStr="(eventNumber % 10)>=3";
  
  //std::map<TString, float> xs_DY, xs_gg, xs_tloop, xs_gamma, delta_EW;
  // NNLO QCD DY W+H cross sections [fb]
  std::map<TString, float> xs_DY = {
    {"WplusH" , 97.18  },
    {"WminusH", 61.51  },
    {"ZllH"   , 26.66  },
    {"ZnnH"   , 158.10 }
  };
  // Non-DY contributions to the XS which aren't susceptible to the EW corrections [fb]
  std::map<TString, float> xs_gg = {
    {"WplusH" , 0    },
    {"WminusH", 0    },
    {"ZllH"   , 4.14 },
    {"ZnnH"   , 24.57}
  };
  std::map<TString, float> xs_tloop = {
    {"WplusH" , 1.20},
    {"WminusH", 0.78},
    {"ZllH"   , 0.31},
    {"ZnnH"   , 1.85}
  };
  std::map<TString, float> xs_gamma = {
    {"WplusH" , 3.09},
    {"WminusH", 2.00},
    {"ZllH"   , 0.11},
    {"ZnnH"   , 0.00}
  };
  // Fractional inclusive EW NLO corrections
  std::map<TString, float> delta_EW = {
    {"WplusH" , -0.074},
    {"WminusH", -0.073},
    {"ZllH"   , -0.053},
    {"ZnnH"   , -0.044}
  }; 
  enum selectionType { 
    kWHLightFlavorCR   = 0x1 <<  0,
    kWHHeavyFlavorCR   = 0x1 <<  1,
    kWH2TopCR          = 0x1 <<  2,
    kWHSR              = 0x1 <<  3,
    kWHPresel          = 0x1 <<  4,
    kWHLightFlavorFJCR = 0x1 <<  5,
    kWHHeavyFlavorFJCR = 0x1 <<  6,
    kWH2TopFJCR        = 0x1 <<  7,
    kWHFJSR            = 0x1 <<  8, 
    kWHFJPresel        = 0x1 <<  9,
    kZnnHLightFlavorCR = 0x1 << 10,
    kZnnHHeavyFlavorCR = 0x1 << 11,
    kZnnH2TopCR        = 0x1 << 12,
    kZnnHMultijetCR    = 0x1 << 13,
    kZnnHSR            = 0x1 << 14,
    kZllHLightFlavorCR = 0x1 << 15,
    kZllHHeavyFlavorCR = 0x1 << 16,
    kZllH2TopCR        = 0x1 << 17,
    kZllHSR            = 0x1 << 18 
  };
  enum sampleType {
    kData       , // 0 
    kQCD        , // 1
    kVZ         , // 2
    kWW         , // 3
    kTT         , // 4
    kTop        , // 5
    kWjets      , // 6
    kZjets      , // 7
    kVH         , // 7
    nSampleTypes 
  };
  enum plotCategory {
    kPlotData , // 0  
    kPlotQCD  , // 1  
    kPlotVZbb , // 2  
    kPlotVVLF , // 3  
    kPlotTT   , // 4  
    kPlotTop  , // 5  
    kPlotWbb  , // 6  
    kPlotWb   , // 7  
    kPlotWLF  , // 8  
    kPlotZbb  , // 9  
    kPlotZb   , //10  
    kPlotZLF  , //11   
    kPlotVH   , //12   
    nPlotCategories
  };
  
  std::map<int, int> plotColors={
    { kPlotData , kBlack      },
    { kPlotQCD  , kGray       },
    { kPlotVZbb , 842         },
    { kPlotVVLF , kAzure-9    },
    { kPlotTT   , kOrange-2   },
    { kPlotTop  , kYellow-6   },
//    { kPlotTop  , kPink+7     },
    { kPlotWbb  , kViolet+2   },
    { kPlotWb   , kViolet+8   },
    { kPlotWLF  , kViolet+6   },
    { kPlotZbb  , kRed-8      },
    { kPlotZb   , kMagenta-10 },
    { kPlotZLF  , kPink+1     },
//    { kPlotVH   , kOrange+9   }
    { kPlotVH   , kRed+1      }
  };
  std::map<plotCategory, TString> plotNames={
    { kPlotData , "Data"     },
    { kPlotQCD  , "QCD"      },
    { kPlotVZbb , "VZ(bb)"   },
    { kPlotVVLF , "VV+LF"    },
    { kPlotTT   , "t#bar{t}" },
    { kPlotTop  , "Top"      },
    { kPlotWbb  , "W+bb"     },
    { kPlotWb   , "W+b"      },
    { kPlotWLF  , "W+udcsg"  },
    { kPlotZbb  , "Z+bb"     },
    { kPlotZb   , "Z+b"      },
    { kPlotZLF  , "Z+udcsg"  },
    { kPlotVH   , "WH(125)"  }
  };
  std::map<int, TString> plotBaseNames={
    { kPlotData , "Data" },
    { kPlotQCD  , "QCD"  },
    { kPlotVZbb , "VZbb" },
    { kPlotVVLF , "VVLF" },
    { kPlotTT   , "TT"   },
    { kPlotTop  , "Top"  },
    { kPlotWbb  , "Wbb"  },
    { kPlotWb   , "Wb"   },
    { kPlotWLF  , "WLF"  },
    { kPlotZbb  , "Zbb"  },
    { kPlotZb   , "Zb"   },
    { kPlotZLF  , "ZLF"  },
    { kPlotVH   , "VH"   }
  }; 
  std::map<int, TString> selectionNames={ 
    { kWHLightFlavorCR    , "WHLightFlavorCR"    },
    { kWHHeavyFlavorCR    , "WHHeavyFlavorCR"    },
    { kWH2TopCR           , "WH2TopCR"           },
    { kWHSR               , "WHSR"               },
    { kWHLightFlavorFJCR  , "WHLightFlavorFJCR"  },
    { kWHHeavyFlavorFJCR  , "WHHeavyFlavorFJCR"  },
    { kWH2TopFJCR         , "WH2TopFJCR"         },
    { kWHFJSR             , "WHFJSR"             },
    { kZnnHLightFlavorCR  , "ZnnHLightFlavorCR"  },
    { kZnnHHeavyFlavorCR  , "ZnnHHeavyFlavorCR"  },
    { kZnnH2TopCR         , "ZnnH2TopCR"         },
    { kZnnHMultijetCR     , "ZnnHMultijetCR"     },
    { kZnnHSR             , "ZnnHSR"             },
    { kZllHLightFlavorCR  , "ZllHLightFlavorCR"  },
    { kZllHHeavyFlavorCR  , "ZllHHeavyFlavorCR"  },
    { kZllH2TopCR         , "ZllH2TopCR"         },
    { kZllHSR             , "ZllHSR"             }
  };
  // This function loads the ith entry of the branch, only if it has not already been loaded
  int bLoad(TBranch *branch, Long64_t ientry) {
    if(!branch) return 0;
    int bytesRead=0;
    Long64_t readEntry = branch->GetReadEntry();
    if(readEntry != ientry) bytesRead = branch->GetEntry(ientry);
    return bytesRead;
  }
  bool passAllCuts( std::map<TString, bool> cutMap, vector<TString> theCuts, bool debug=false ) {
    unsigned nCutsToPass=theCuts.size();
    for(unsigned i=0; i<nCutsToPass; i++)
      if(cutMap.find(theCuts[i])!=cutMap.end()) if(!cutMap[theCuts[i]]) {
        if(debug) printf("\tvhbbPlot::passAllCuts: Failed cut \"%s\"\n", theCuts[i].Data());
        return false;
      }
    return true;
  }
  bool passNMinusOne( std::map<TString, bool> cutMap, vector<TString> theCuts ) {
    unsigned nCutsToPass=theCuts.size();
    unsigned nCutsPassed=0;
    for(unsigned i=0; i<nCutsToPass; i++) {
      if(cutMap.find(theCuts[i])!=cutMap.end())
        nCutsPassed += (unsigned)cutMap[theCuts[i]];
    }
    return(nCutsPassed >= nCutsToPass-1);
  }

  double qcdKFactor(vhbbPlot::sampleType sample, unsigned htLow, unsigned htHigh) {
    double theQCDKFactor=1;
    if(sample==kWjets) {
      if     (htLow>=100 && htHigh<= 200) theQCDKFactor = 1.459;
      else if(htLow>=200 && htHigh<= 400) theQCDKFactor = 1.434;
      else if(htLow>=400 && htHigh<= 600) theQCDKFactor = 1.532;
      else if(htLow>=600                ) theQCDKFactor = 1.004;
    } else if(sample==kZjets) {
      if     (htLow>=100 && htHigh<= 200) theQCDKFactor = 1.588;
      else if(htLow>=200 && htHigh<= 400) theQCDKFactor = 1.438;
      else if(htLow>=400 && htHigh<= 600) theQCDKFactor = 1.494;
      else if(htLow>=600                ) theQCDKFactor = 1.139;
    }
    return theQCDKFactor;
  }

  vector<double> EWKCorrPars(vhbbPlot::sampleType sample) {
    vector<double> theEWKCorrPars(4);
    if(sample==kWjets)
      theEWKCorrPars = {-0.830041, 7.93714, 877.978, -0.213831};
    else if(sample==kZjets)
      theEWKCorrPars = {-0.1808051, 6.04146, 759.098, -0.242556};
    return theEWKCorrPars;
  }
  
  float vhEWKCorr( TString x, float diffPtCorr) {
    if(x!="WplusH" && x!="WminusH" && x!="ZllH" && x!="ZnnH") return 1;
    return 
      (xs_DY[x]*diffPtCorr+xs_gg[x]+xs_tloop[x]+xs_gamma[x])/
      (xs_DY[x]*(1.+delta_EW[x])+xs_gg[x]+xs_tloop[x]+xs_gamma[x]);
  }

}
#endif
