#ifndef vhbbPlot_h
#define vhbbPlot_h

#include <TBranch.h>
#include <map>

namespace vhbbPlot {
  
  float theLumi=35900.;
  //const float bDiscrLoose = 0.5426, bDiscrMedium = 0.8484, bDiscrTight  = 0.9535; //csv
  const float cmvaLoose = -0.5884, cmvaMedium = 0.4432, cmvaTight  = 0.9432; //cmva
  const float deepcsvLoose = 0.1522, deepcsvMedium = 0.4941, deepcsvTight = 0.8001; //cmva
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
    kWHLightFlavorCR     = 0x1 <<  0,
    kWHHeavyFlavorCR     = 0x1 <<  1,
    kWH2TopCR            = 0x1 <<  2,
    kWHSR                = 0x1 <<  3,
    kWHPresel            = 0x1 <<  4,
    kWHLightFlavorFJCR   = 0x1 <<  5,
    kWHHeavyFlavorFJCR   = 0x1 <<  6,
    kWHTT2bFJCR          = 0x1 <<  7,
    kWHTT1bFJCR          = 0x1 <<  8,
    kWHFJSR              = 0x1 <<  9, 
    kWHFJPresel          = 0x1 << 10,
    kZllHLightFlavorCR   = 0x1 << 11,
    kZllHHeavyFlavorCR   = 0x1 << 12,
    kZllH2TopCR          = 0x1 << 13,
    kZllHSR              = 0x1 << 14,
    kZllHPresel          = 0x1 << 15,
    kZllHLightFlavorFJCR = 0x1 << 16,
    kZllHHeavyFlavorFJCR = 0x1 << 17,
    kZllHTT2bFJCR        = 0x1 << 18,
    kZllHTT1bFJCR        = 0x1 << 19,
    kZllHFJSR            = 0x1 << 20,
    kZllHFJPresel        = 0x1 << 21,
    
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
    kWH         , // 8
    kZH         , // 8
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
    kPlotWH   , //12   
    kPlotZH   , //13   
    nPlotCategories
  };
  
  std::map<int, int> plotColors={
    { kPlotData , kBlack      },
    { kPlotQCD  , kGray       },
    { kPlotVZbb , 842         },
    { kPlotVVLF , kAzure-9    },
    { kPlotTT   , kOrange-2   },
    { kPlotTop  , kOrange-1   },
    { kPlotWbb  , kViolet+3   },
    { kPlotWb   , kViolet+8   },
    { kPlotWLF  , kViolet+6   },
    { kPlotZbb  , kRed-8      },
    { kPlotZb   , kMagenta-10 },
    { kPlotZLF  , kPink+1     },
    { kPlotWH   , kRed+1      },
    { kPlotZH   , kRed+3      }
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
    { kPlotWH   , "WH(125)"  },
    { kPlotZH   , "ZH(125)"  }
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
    { kPlotWH   , "WH"   },
    { kPlotZH   , "WH"   }
  }; 
  std::map<int, TString> selectionNames={ 
    { kWHLightFlavorCR    , "WHLightFlavorCR"    },
    { kWHHeavyFlavorCR    , "WHHeavyFlavorCR"    },
    { kWH2TopCR           , "WH2TopCR"           },
    { kWHSR               , "WHSR"               },
    { kWHPresel           , "WHPresel"           },
    { kWHLightFlavorFJCR  , "WHLightFlavorFJCR"  },
    { kWHHeavyFlavorFJCR  , "WHHeavyFlavorFJCR"  },
    { kWHTT2bFJCR         , "WHTT2bFJCR"         },
    { kWHTT1bFJCR         , "WHTT1bFJCR"         },
    { kWHFJSR             , "WHFJSR"             },
    { kWHFJPresel         , "WHFJPresel"         },
    { kZllHLightFlavorCR    , "ZllHLightFlavorCR"    },
    { kZllHHeavyFlavorCR    , "ZllHHeavyFlavorCR"    },
    { kZllH2TopCR           , "ZllH2TopCR"           },
    { kZllHSR               , "ZllHSR"               },
    { kZllHPresel           , "ZllHPresel"           },
    { kZllHLightFlavorFJCR  , "ZllHLightFlavorFJCR"  },
    { kZllHHeavyFlavorFJCR  , "ZllHHeavyFlavorFJCR"  },
    { kZllHTT2bFJCR         , "ZllHTT2bFJCR"         },
    { kZllHTT1bFJCR         , "ZllHTT1bFJCR"         },
    { kZllHFJSR             , "ZllHFJSR"             },
    { kZllHFJPresel         , "ZllHFJPresel"         }
  };
  std::map<int, float> pdfAcceptUncs={
    { kPlotVZbb , 1.02  },
    { kPlotVVLF , 1.03  },
    { kPlotTT   , 1.005 },
    { kPlotWbb  , 1.02  },
    { kPlotWb   , 1.03  },
    { kPlotWLF  , 1.05  },
    { kPlotZbb  , 1.02  },
    { kPlotZb   , 1.03  },
    { kPlotZLF  , 1.05  },
    { kPlotWH   , 1.01  },
    { kPlotZH   , 1.01  }
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

  double qcdKFactor(vhbbPlot::sampleType sample, float lhe_HT) {
    double theQCDKFactor=1;
    if(sample==kWjets) {
      if     (lhe_HT>=100 && lhe_HT<200) theQCDKFactor = 1.459;
      else if(lhe_HT>=200 && lhe_HT<400) theQCDKFactor = 1.434;
      else if(lhe_HT>=400 && lhe_HT<600) theQCDKFactor = 1.532;
      else if(lhe_HT>=600              ) theQCDKFactor = 1.004;
    } else if(sample==kZjets) {
      if     (lhe_HT>=100 && lhe_HT<200) theQCDKFactor = 1.588;
      else if(lhe_HT>=200 && lhe_HT<400) theQCDKFactor = 1.438;
      else if(lhe_HT>=400 && lhe_HT<600) theQCDKFactor = 1.494;
      else if(lhe_HT>=600              ) theQCDKFactor = 1.139;
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
  float BenrichedVPT100   = 0.108880411891;
  float BenrichedVPT200   = 0.0679764068637;
  float BfilterVPT100     = 0.206608510333;
  float BfilterVPT200     = 0.16396633746;

}
#endif
