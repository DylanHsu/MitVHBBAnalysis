#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>

void plotSystShape(TString process="WLF", TString syst="QCDrScaleWLF",TString histDir="", TString region="", float ymin=.8,float ymax=1.2)
{
  TString filename = Form("%s/datacard_%s.root",histDir.Data(), region.Data());
  system(Form("mkdir -p %s/systs", histDir.Data()));
  TFile *file=TFile::Open(filename,"READ");
  vector<TString> systNames;
  TString shapeType;
  if(region.Contains("WHSR")||region.Contains("WHFJSR")) shapeType="singleClassBDTShape";
  else if(region.Contains("FJCR")) shapeType="softDropMassShape";
  else shapeType = "lesserCMVAShape";

  if(syst=="") {
    TList *listOfHistoNames=file->GetListOfKeys();
    for(unsigned iHisto=0; iHisto<=(unsigned)listOfHistoNames->LastIndex(); iHisto++) {
      TString theHistoName = listOfHistoNames->At(iHisto)->GetName();
      //TString token = Form("%s_%s_%s_", shapeType.Data(), region.Data(), process.Data());
      TString token = Form("histo_%s_", process.Data());
      if(!theHistoName.BeginsWith(token) ||
         !theHistoName.EndsWith("Up")) continue;
      
      systNames.push_back(theHistoName(token.Length(), theHistoName.Length()-2-token.Length()));
      printf("Found syst named \"%s\" from histo \"%s\"\n", systNames[systNames.size()-1].Data(), theHistoName.Data());
    }
  } else systNames.push_back(syst);
  for(unsigned iSyst=0; iSyst<systNames.size(); iSyst++) {
    TString theSyst = systNames[iSyst];
    //TString histoNomName  = Form("%s_%s_%s"       , shapeType.Data(), region.Data(), process.Data()             );
    //TString histoUpName   = Form("%s_%s_%s_%sUp"  , shapeType.Data(), region.Data(), process.Data(), theSyst.Data());
    //TString histoDownName = Form("%s_%s_%s_%sDown", shapeType.Data(), region.Data(), process.Data(), theSyst.Data());
    TString histoNomName  = Form("histo_%s"       , process.Data()             );
    TString histoUpName   = Form("histo_%s_%sUp"  , process.Data(), theSyst.Data());
    TString histoDownName = Form("histo_%s_%sDown", process.Data(), theSyst.Data());
    TH1F *histoNom=0, *histoUp=0,*histoDown=0;
    histoNom  = (TH1F*)file->Get(histoNomName ); assert(histoNom ); histoNom ->SetDirectory(0); 
    histoUp   = (TH1F*)file->Get(histoUpName  ); assert(histoUp  ); histoUp  ->SetDirectory(0); 
    histoDown = (TH1F*)file->Get(histoDownName); assert(histoDown); histoDown->SetDirectory(0); 
    TCanvas *c=new TCanvas(Form("c_%s",region.Data()),region);
    gStyle->SetOptStat(0);
    histoUp->SetLineColor(kViolet); histoDown->SetLineColor(kBlue); histoNom->SetLineColor(kBlack);
    histoUp->SetMaximum(1.2*TMath::Max(histoNom->GetMaximum(),TMath::Max(histoUp->GetMaximum(),histoDown->GetMaximum())));
    histoUp->Draw("hist");
    histoNom->Draw("hist same");
    histoDown->Draw("hist same");
    histoUp->SetTitle(Form("Abs. syst. %s for process %s in region %s", theSyst.Data(), process.Data(), region.Data()));
    c->Print(Form("%s/systs/systAbs_%s_%s_%s.pdf",histDir.Data(),region.Data(),process.Data(),theSyst.Data()));
    c->Print(Form("%s/systs/systAbs_%s_%s_%s.png",histDir.Data(),region.Data(),process.Data(),theSyst.Data()));

    histoUp->Divide(histoNom);
    histoDown->Divide(histoNom);
    histoUp->SetMinimum(ymin); histoUp->SetMaximum(ymax);
    histoUp->SetTitle(Form("Rel syst. %s for process %s in region %s", theSyst.Data(), process.Data(), region.Data()));
    histoUp->SetLineWidth(2); histoDown->SetLineWidth(2);
    histoUp->Draw("HIST"); 
    TLine *oneline = new TLine(histoDown->GetBinLowEdge(1),1,histoDown->GetBinLowEdge(histoDown->GetNbinsX()+1),1);
    oneline->SetLineWidth(2);
    oneline->SetLineColor(kBlack);
    oneline->SetLineStyle(kDashed);
    oneline->Draw("same");
    histoUp->Draw("HIST same"); 
    histoUp->GetXaxis()->SetTitle(shapeType);
    histoUp->GetYaxis()->SetTitle("Relative shape uncertainty");
    histoDown->Draw("HIST SAME");
  
    c->Print(Form("%s/systs/syst_%s_%s_%s.pdf",histDir.Data(),region.Data(),process.Data(),theSyst.Data()));
    c->Print(Form("%s/systs/syst_%s_%s_%s.png",histDir.Data(),region.Data(),process.Data(),theSyst.Data()));
  }
}
