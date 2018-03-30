#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>

void plotSystShape(TString process="WLF", TString syst="QCDscaleWLF",TString histDir="", TString region="")
{
  vector<TString>filenames,regions;
  if(histDir=="" || region=="") {
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WenWH2TopCR.root"               ); regions.push_back("WenWH2TopCR"               );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WenWHHeavyFlavorCRHighMass.root"); regions.push_back("WenWHHeavyFlavorCRHighMass");
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WenWHHeavyFlavorCRLowMass.root" ); regions.push_back("WenWHHeavyFlavorCRLowMass" );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WenWHLightFlavorCR.root"        ); regions.push_back("WenWHLightFlavorCR"        );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WenWHSR.root"                   ); regions.push_back("WenWHSR"                   );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WmnWH2TopCR.root"               ); regions.push_back("WmnWH2TopCR"               );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WmnWHHeavyFlavorCRHighMass.root"); regions.push_back("WmnWHHeavyFlavorCRHighMass");
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WmnWHHeavyFlavorCRLowMass.root" ); regions.push_back("WmnWHHeavyFlavorCRLowMass" );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WmnWHLightFlavorCR.root"        ); regions.push_back("WmnWHLightFlavorCR"        );
    filenames.push_back("MitVHBBAnalysis/datacards/all_orthodox11vars_singleClass/oldHistos/hists_WmnWHSR.root"                   ); regions.push_back("WmnWHSR"                   );
  } else {
    filenames.push_back(Form("%s/hists_%s.root",histDir.Data(), region.Data())); regions.push_back(region);
  }
  TCanvas *c[99];
  for(unsigned i=0; i<filenames.size(); i++) {
    TFile *file=TFile::Open(filenames[i],"READ");
    TString shapeType;
    if(regions[i].Contains("WHSR")) shapeType="singleClassBDTShape";
    else if(regions[i].Contains("FJCR")) shapeType="softDropMassShape";
    else shapeType = "lesserCMVAShape";
    TString histoNomName  = Form("%s_%s_%s"       , shapeType.Data(), regions[i].Data(), process.Data()             );
    TString histoUpName   = Form("%s_%s_%s_%sUp"  , shapeType.Data(), regions[i].Data(), process.Data(), syst.Data());
    TString histoDownName = Form("%s_%s_%s_%sDown", shapeType.Data(), regions[i].Data(), process.Data(), syst.Data());
    TH1F *histoNom=0, *histoUp=0,*histoDown=0;
    histoNom  = (TH1F*)file->Get(histoNomName ); assert(histoNom ); histoNom ->SetDirectory(0); 
    histoUp   = (TH1F*)file->Get(histoUpName  ); assert(histoUp  ); histoUp  ->SetDirectory(0); 
    histoDown = (TH1F*)file->Get(histoDownName); assert(histoDown); histoDown->SetDirectory(0); 
    c[i]=new TCanvas(Form("c_%s",regions[i].Data()),regions[i]);
    gStyle->SetOptStat(0);
    histoUp->Divide(histoNom);
    histoDown->Divide(histoNom);
    histoUp->SetMinimum(0.8); histoUp->SetMaximum(1.2);
    histoUp->SetTitle(Form("Syst. %s for process %s in region %s", syst.Data(), process.Data(), regions[i].Data()));
    histoUp->SetLineColor(kViolet); histoDown->SetLineColor(kBlue);
    histoUp->SetLineWidth(2); histoDown->SetLineWidth(2);
    histoUp->Draw("HIST"); 
    histoUp->GetXaxis()->SetTitle(shapeType);
    histoUp->GetYaxis()->SetTitle("Relative shape uncertainty");
    histoDown->Draw("HIST SAME");
    c[i]->Print(Form("MitVHBBAnalysis/plots/syst_%s_%s_%s.pdf",regions[i].Data(),process.Data(),syst.Data()));
  }
}
