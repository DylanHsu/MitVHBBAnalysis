#ifndef vhbbPlot_h
#define vhbbPlot_h

#include <TBranch.h>
#include <map>

namespace vhbbPlot {
  
  float theLumi=35900.;
  //const float bDiscrLoose = 0.5426, bDiscrMedium = 0.8484, bDiscrTight  = 0.9535; //csv
  const float bDiscrLoose = -0.5884, bDiscrMedium = 0.4432, bDiscrTight  = 0.9432; //cmva
  const float doubleBCut = 0.8; // double b tagger

  enum selectionType { 
    kWHLightFlavorCR   = 0x1 <<  0,
    kWHHeavyFlavorCR   = 0x1 <<  1,
    kWH2TopCR          = 0x1 <<  2,
    kWHSR              = 0x1 <<  3,
    kWHLightFlavorFJCR = 0x1 <<  4,
    kWHHeavyFlavorFJCR = 0x1 <<  5,
    kWH2TopFJCR        = 0x1 <<  6,
    kWHFJSR            = 0x1 <<  7, 
    kZnnHLightFlavorCR = 0x1 <<  8,
    kZnnHHeavyFlavorCR = 0x1 <<  9,
    kZnnH2TopCR        = 0x1 << 10,
    kZnnHMultijetCR    = 0x1 << 11,
    kZnnHSR            = 0x1 << 12,
    kZllHLightFlavorCR = 0x1 << 13,
    kZllHHeavyFlavorCR = 0x1 << 14,
    kZllH2TopCR        = 0x1 << 15,
    kZllHSR            = 0x1 << 16 
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
  bool passAllCuts( std::map<TString, bool> cutMap, vector<TString> theCuts ) {
    unsigned nCutsToPass=theCuts.size();
    for(unsigned i=0; i<nCutsToPass; i++)
      if(cutMap.find(theCuts[i])!=cutMap.end()) if(!cutMap[theCuts[i]])
        return false;
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


}
#endif
