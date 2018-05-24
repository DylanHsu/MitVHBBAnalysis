#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TList.h>
#include <cassert>
#include "RooArgSet.h"

#include "vhbbPlot.h"

using namespace vhbbPlot;

//TList *finalPlot2018(
void finalPlot2018(
  vhbbPlot::selectionType selType, 
  TString inputFileName,
  TString plotTitle="",
  TString extraText="",
  bool isBlinded=false,
  bool normSignalToBkg=false,
  TString plotDir="MitVHBBAnalysis/plots",
  TString mlfitResult=""
) {
  // Bkg only SF
  float SF_TT_Wln         = 1;
  float SF_WLF_Wln        = 1;
  float SF_Wb_Wln         = 1;
  float SF_Wbb_Wln        = 1;
  float SF_Top_Wln        = 1;

  int primarySignal = (selType>=kWHLightFlavorCR && selType<=kWHFJPresel)? kPlotWH:kPlotZH;

  system(Form("mkdir -p %s",plotDir.Data()));

  TFile *inputFile = TFile::Open(inputFileName, "READ");
  if(!inputFile) return; // assert(inputFile);
  TFile *mlfit=0;
  string regionName; {
    size_t lastDot = string(inputFileName).find_last_of(".");
    size_t lastSlash = string(inputFileName).find_last_of("/");
    regionName = string(inputFileName).substr(lastSlash+7, lastDot-lastSlash-7);
  }
  
  if(mlfitResult!="") {
    mlfit=TFile::Open(mlfitResult); assert(mlfit);
    RooArgSet *norm_prefit   = (RooArgSet*)mlfit->Get("norm_prefit");
    RooArgSet *norm_postfit  = (RooArgSet*)mlfit->Get("norm_fit_s" );
    inputFile->cd();
    //RooArgSet *norm_postfit  = (RooArgSet*)mlfit->Get("norm_fit_s" );
    SF_TT_Wln  = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotTT ].Data()))/norm_prefit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotTT ].Data()));
    SF_WLF_Wln = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWLF].Data()))/norm_prefit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWLF].Data()));
    SF_Wb_Wln  = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWb ].Data()))/norm_prefit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWb ].Data()));
    SF_Wbb_Wln = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWbb].Data()))/norm_prefit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWbb].Data()));
    SF_Top_Wln = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotTop].Data()))/norm_prefit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotTop].Data()));
    //SF_TT_Wln  = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotTT ].Data()))/((TH1F*)inputFile->Get(Form("MVAVar/histo%d",kPlotTT )))->Integral();
    //SF_WLF_Wln = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWLF].Data()))/((TH1F*)inputFile->Get(Form("MVAVar/histo%d",kPlotWLF)))->Integral();
    //SF_Wb_Wln  = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWb ].Data()))/((TH1F*)inputFile->Get(Form("MVAVar/histo%d",kPlotWb )))->Integral();
    //SF_Wbb_Wln = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotWbb].Data()))/((TH1F*)inputFile->Get(Form("MVAVar/histo%d",kPlotWbb)))->Integral();
    //SF_Top_Wln = norm_postfit->getRealValue(Form("%s/%s",regionName.c_str(),plotBaseNames[kPlotTop].Data()))/((TH1F*)inputFile->Get(Form("MVAVar/histo%d",kPlotTop)))->Integral();
    printf("Renorm. from mlfit:\n");
    printf("k_TT_Wln  = %.3f\n", SF_TT_Wln  );
    printf("k_WLF_Wln = %.3f\n", SF_WLF_Wln );
    printf("k_Wb_Wln  = %.3f\n", SF_Wb_Wln  );
    printf("k_Wbb_Wln = %.3f\n", SF_Wbb_Wln );
    printf("k_Top_Wln = %.3f\n", SF_Top_Wln );
  }

  TList *listOfHistoNames=inputFile->GetListOfKeys();
  //TList *listOfCanvases=new TList();
  for(unsigned iHisto=0; iHisto<=(unsigned)listOfHistoNames->LastIndex(); iHisto++) {
    TString theHistoName = listOfHistoNames->At(iHisto)->GetName();
    //if(theHistoName.Contains("psi")) continue;
    bool isFitShape = theHistoName.Contains("MVAVar");
    bool stackSignal = isFitShape;
    bool doLogPlot = stackSignal;// || theHistoName.Contains("HTT") || theHistoName.Contains("pfmetsig");
    bool doRatioPad = !isBlinded || stackSignal;
    bool plotQCD=true; if(selType==kWHSR || selType==kWHFJSR /* || (selType>=kWHLightFlavorFJCR && selType<=kWHFJSR)*/ || isFitShape) plotQCD=false;
    TString outPdf = Form("%s/%s_%s.pdf", plotDir.Data(), regionName.c_str(), theHistoName.Data()); 
    TString outPng = Form("%s/%s_%s.png", plotDir.Data(), regionName.c_str(), theHistoName.Data()); 
    TString xlabel=""; TString plotName;
    TH1F *histos[nPlotCategories], *hTotalBkg=0, *hTotalBkgPrefit=0;
    for(int iCat=kPlotData; iCat!=nPlotCategories; iCat++) {
      plotCategory i = static_cast<plotCategory>(iCat);
      if(!plotQCD && i==kPlotQCD) continue;
      if (i==kPlotWH) {
        if(selType==kWHSR || selType==kWHFJSR) {
          if(normSignalToBkg) plotName="WH(125)x?";
          else if(!stackSignal) plotName="WH(125)x10"; 
          else plotName="WH(125)";
        } else if(selType>=kWHLightFlavorCR && selType<kWHFJPresel) {
          if(normSignalToBkg) plotName="WH(125)x?";
          else if(!stackSignal) plotName="WH(125)x100";
          else plotName="WH(125)";
        }
      } else if (i==kPlotZH) {
        if(selType==kZllHSR || selType==kZllHFJSR) {
          if(normSignalToBkg) plotName="ZH(125)x?";
          else if(!stackSignal) plotName="ZH(125)x10"; 
          else plotName="ZH(125)";
        } else if(selType>=kZllHLightFlavorCR && selType<kZllHFJPresel) {
          if(normSignalToBkg) plotName="ZH(125)x?";
          else if(!stackSignal) plotName="ZH(125)x100";
          else plotName="ZH(125)";
        }
      } else plotName=plotNames[i]; 

      // Fill histograms
      histos[i]=(TH1F*)inputFile->Get(Form("%s/histo%d",theHistoName.Data(),iCat));
      if(!histos[i]) continue;
      histos[i]->SetDirectory(0);
      // Colors
      if(i==kPlotData) {
        histos[i]->SetMarkerColor(plotColors[i]); 
        histos[i]->SetLineColor(plotColors[i]); 
        histos[i]->SetMarkerStyle(20);
      } else if(i!=kPlotWH && i!=kPlotZH) {
        histos[i]->SetLineColor(plotColors[i]);
        histos[i]->SetFillColor(plotColors[i]);
        histos[i]->SetMarkerColor(plotColors[i]);
        histos[i]->SetFillStyle(1001);
      }
      if(iCat==kPlotGGZH) 
        histos[(int)kPlotZH]->Add(histos[i]);
      if(iCat==kPlotData ||iCat==kPlotQCD || iCat==kPlotWH || iCat==kPlotZH) {
        // do nothing currently
      } else if(!isFitShape || !mlfit) {
        if     (i==kPlotTT  ) histos[i]->Scale( SF_TT_Wln);
        else if(i==kPlotWLF ) histos[i]->Scale(SF_WLF_Wln);
        else if(i==kPlotWb  ) histos[i]->Scale( SF_Wb_Wln);
        else if(i==kPlotWbb ) histos[i]->Scale(SF_Wbb_Wln);
        else if(i==kPlotTop ) histos[i]->Scale(SF_Top_Wln);
      } else {
        // If we have a MLF result, and this is the shape variable,
        // use the shape from background only postfit
        printf("shapes_fit_s/%s/%s\n",regionName.c_str(), plotBaseNames[iCat].Data());
        TH1F *postfitHisto = (TH1F*)mlfit->Get(Form("shapes_fit_s/%s/%s",regionName.c_str(), plotBaseNames[iCat].Data()));
        TH1F *postfitHisto_b = (TH1F*)mlfit->Get(Form("shapes_fit_b/%s/%s",regionName.c_str(), plotBaseNames[iCat].Data()));
        if(!postfitHisto || !postfitHisto_b) continue;
        for(int nb=1; nb<=postfitHisto->GetNbinsX(); nb++) {
          histos[i]->SetBinContent(nb, postfitHisto->GetBinContent(nb));
          histos[i]->SetBinError(nb, postfitHisto->GetBinError(nb));
          //printf("%s bin%d %.2f\n", plotBaseNames[iCat].Data(), nb, postfitHisto->GetBinError(nb)/postfitHisto->GetBinContent(nb));
        }
      }
      if(xlabel=="") xlabel=histos[i]->GetTitle();
      histos[i]->SetName(plotName);
      histos[i]->SetTitle(plotTitle);

    }
    if(mlfit && isFitShape) {
      TH1F *postfitTotalBkg = (TH1F*)mlfit->Get(Form("shapes_fit_s/%s/total_background",regionName.c_str()));
      TH1F *prefitTotalBkg = (TH1F*)mlfit->Get(Form("shapes_prefit/%s/total_background",regionName.c_str()));
      hTotalBkg=(TH1F*)histos[kPlotData]->Clone("hTotal");
      hTotalBkgPrefit=(TH1F*)histos[kPlotData]->Clone("hTotalPrefit");
      for(int nb=1; nb<=histos[kPlotData]->GetNbinsX(); nb++) {
        hTotalBkg->SetBinContent(nb, postfitTotalBkg->GetBinContent(nb));
        hTotalBkg->SetBinError(nb, postfitTotalBkg->GetBinError(nb));
        hTotalBkgPrefit->SetBinContent(nb, prefitTotalBkg->GetBinContent(nb));
        hTotalBkgPrefit->SetBinError(nb, prefitTotalBkg->GetBinError(nb));
      }
    }
    string units=""; {
      size_t lastOpenBrace = string(xlabel).find_last_of("[");
      size_t lastCloseBrace = string(xlabel).find_last_of("]");
      if(lastOpenBrace!=string::npos && lastCloseBrace!=string::npos) 
        units = string(xlabel).substr(lastOpenBrace+1, lastCloseBrace-lastOpenBrace-1);
    }
    // Bin width calculation
    int nbins=histos[kPlotData]->GetXaxis()->GetNbins();
    bool variableWidth=false; float binWidth=-1;
    for(int nb=1; nb<=nbins; nb++) { 
      if(nb==1) { binWidth=histos[kPlotData]->GetXaxis()->GetBinWidth(nb); continue; }
      if(!variableWidth && TMath::Abs(binWidth/((float)histos[kPlotData]->GetXaxis()->GetBinWidth(nb))-1.) > 0.0001) variableWidth=true;
      if(histos[kPlotData]->GetXaxis()->GetBinWidth(nb) < binWidth) binWidth=histos[kPlotData]->GetXaxis()->GetBinWidth(nb);
    }
    if(stackSignal) {
      histos[kPlotWH]->SetFillColor(plotColors[kPlotWH]);
      histos[kPlotWH]->SetLineColor(plotColors[kPlotWH]);
      histos[kPlotZH]->SetFillColor(plotColors[kPlotZH]);
      histos[kPlotZH]->SetLineColor(plotColors[kPlotZH]);
    } else {
      // Scaling
      if(!normSignalToBkg && !stackSignal && selType!=kWHSR && selType!=kWHFJSR /*&& selType!=kZnnHSR && selType!=kZllHSR*/) histos[kPlotWH]->Scale(100.);
      if(!normSignalToBkg && !stackSignal && (selType==kWHSR || selType==kWHFJSR /*|| selType==kZnnHSR || selType==kZllHSR*/)) histos[kPlotWH]->Scale(10.);
      histos[kPlotWH]->SetMarkerColor(plotColors[kPlotWH]); 
      histos[kPlotWH]->SetLineColor(plotColors[kPlotWH]); 
      histos[kPlotWH]->SetLineWidth(3);
      histos[kPlotWH]->SetFillStyle(0);
      if(!normSignalToBkg && !stackSignal && selType!=kZllHSR && selType!=kZllHFJSR /*&& selType!=kZnnHSR && selType!=kZllHSR*/) histos[kPlotZH]->Scale(100.);
      if(!normSignalToBkg && !stackSignal && (selType==kZllHSR || selType==kZllHFJSR /*|| selType==kZnnHSR || selType==kZllHSR*/)) histos[kPlotZH]->Scale(10.);
      histos[kPlotZH]->SetMarkerColor(plotColors[kPlotZH]); 
      histos[kPlotZH]->SetLineColor(plotColors[kPlotZH]); 
      histos[kPlotZH]->SetLineWidth(3);
      histos[kPlotZH]->SetFillStyle(0);
    }
    float xmin=histos[kPlotData]->GetXaxis()->GetBinLowEdge(1);
    float xmax=histos[kPlotData]->GetXaxis()->GetBinLowEdge(nbins+1);
    
    if(mlfit && isFitShape && variableWidth) for(int nb=1; nb<=nbins; nb++) {
      hTotalBkg->SetBinContent(nb, hTotalBkg->GetBinContent(nb) / hTotalBkg->GetXaxis()->GetBinWidth(nb)*binWidth);
      hTotalBkg->SetBinError(nb, hTotalBkg->GetBinError(nb) / hTotalBkg->GetXaxis()->GetBinWidth(nb)*binWidth);
      hTotalBkgPrefit->SetBinContent(nb, hTotalBkgPrefit->GetBinContent(nb) / hTotalBkgPrefit->GetXaxis()->GetBinWidth(nb)*binWidth);
      hTotalBkgPrefit->SetBinError(nb, hTotalBkgPrefit->GetBinError(nb) / hTotalBkgPrefit->GetXaxis()->GetBinWidth(nb)*binWidth);
    }

    for(int iCat=kPlotData; iCat!=nPlotCategories; iCat++) {
      plotCategory i = static_cast<plotCategory>(iCat);
      if(!plotQCD && i==kPlotQCD) continue;
      if(variableWidth) for(int nb=1; nb<=nbins; nb++) { 
        if(binWidth!=histos[i]->GetXaxis()->GetBinWidth(nb)) {
          histos[i]->SetBinContent(nb, histos[i]->GetBinContent(nb) / histos[i]->GetXaxis()->GetBinWidth(nb)*binWidth);
          histos[i]->SetBinError(nb, histos[i]->GetBinError(nb) / histos[i]->GetXaxis()->GetBinWidth(nb)*binWidth);
        }
      }
      // Summing Up
      if(!hTotalBkg && !(mlfit && isFitShape)) {
        hTotalBkg=(TH1F*)histos[i]->Clone("hTotal");
        hTotalBkg->Reset(); hTotalBkg->Clear();
        hTotalBkg->SetDirectory(0);
        hTotalBkg->SetTitle("hTotal");
      }
      if(!(mlfit && isFitShape) && hTotalBkg && i!=kPlotData && i!=kPlotWH && i!=kPlotZH) {
        hTotalBkg->Add(histos[i]);
      }
    }
    if(normSignalToBkg) {
      float signalInflationFactor = hTotalBkg->Integral(1., hTotalBkg->GetNbinsX()) / histos[primarySignal]->Integral(1., histos[primarySignal]->GetNbinsX());
      histos[primarySignal]->Scale(signalInflationFactor);
      // needs to be worked on
      if(primarySignal==kPlotWH) {
        if     (signalInflationFactor>1e5) plotName = Form("WH(125)x%.1e", signalInflationFactor);
        else if(signalInflationFactor>  1) plotName = Form("WH(125)x%d", (int)round(signalInflationFactor));
        else                          plotName = "WH(125)";
      } else {
        if     (signalInflationFactor>1e5) plotName = Form("ZH(125)x%.1e", signalInflationFactor);
        else if(signalInflationFactor>  1) plotName = Form("ZH(125)x%d", (int)round(signalInflationFactor));
        else                          plotName = "ZH(125)";
      }
      histos[primarySignal]->SetName(plotName);
      histos[primarySignal]->SetTitle(plotName);
    }
   

    // Start plotting stuff
    gStyle->SetOptStat(0);
    if(isBlinded) for(int nb=1; nb<=histos[kPlotData]->GetNbinsX(); nb++) {
      if(!stackSignal || histos[kPlotData]->GetBinLowEdge(nb)>=0) { histos[kPlotData]->SetBinContent(nb,-9999); histos[kPlotData]->SetBinError(nb,0); }
    }
    //TH1F *hRatio = (TH1F*) (isBlinded? hTotalBkg->Clone("hRatio") : histos[kPlotData]->Clone("hRatio"));
    TH1F *hRatio = (TH1F*) histos[kPlotData]->Clone("hRatio");
    hRatio->SetDirectory(0);
    for(int nb=1; nb<=hRatio->GetNbinsX(); nb++) {
      if(hTotalBkg->GetBinContent(nb)>0 && hRatio->GetBinContent(nb)>0) {
        hRatio->SetBinContent( nb, hRatio->GetBinContent(nb) / hTotalBkg->GetBinContent(nb));
        hRatio->SetBinError( nb, hRatio->GetBinError(nb) / histos[kPlotData]->GetBinContent(nb) / hRatio->GetBinContent(nb)); 
        if(isBlinded && (!stackSignal||histos[kPlotData]->GetBinLowEdge(nb)>=0)) hRatio->SetBinError( nb, hRatio->GetBinError(nb) / hTotalBkg->GetBinContent(nb) / hRatio->GetBinContent(nb)); 
        //else            hRatio->SetBinError( nb, hRatio->GetBinError(nb) / histos[kPlotData]->GetBinContent(nb) / hRatio->GetBinContent(nb)); 
        else            hRatio->SetBinError( nb, histos[kPlotData]->GetBinError(nb) / hTotalBkg->GetBinContent(nb));
      } else {
        hRatio->SetBinContent( nb, 1);
        hRatio->SetBinError(nb, 1);
      }
    } 
 
    TH1F *hErrorBand = (TH1F*)hTotalBkg->Clone("hErrorBand");
    hErrorBand->SetFillColor(kBlack); hErrorBand->SetFillStyle(3254);
    hErrorBand->SetMarkerColor(kBlack); hErrorBand->SetMarkerSize(0);

    TH1F *hRatioBand = (TH1F*)hErrorBand->Clone("hRatioBand");
    TH1F *hRatioBandPrefit=0;
    if(mlfit && isFitShape) {
      hRatioBandPrefit = (TH1F*)hRatio->Clone("hRatioBandPrefit");
    }
    for(int nb=1; nb<=hRatio->GetNbinsX(); nb++) {
      hRatioBand->SetBinContent(nb,1);
      //hRatioBand->SetBinError(nb, hTotalBkg->GetBinError(nb) / hTotalBkg->GetBinContent(nb) / hRatio->GetBinContent(nb));
      hRatioBand->SetBinError(nb, hTotalBkg->GetBinError(nb) / hTotalBkg->GetBinContent(nb));
      if(hRatioBandPrefit) {
        hRatioBandPrefit->SetBinContent(nb, 1);
        hRatioBandPrefit->SetBinError(nb, hTotalBkgPrefit->GetBinError(nb) / hTotalBkgPrefit->GetBinContent(nb));
      }
    }
    THStack *hs = new THStack("hs",plotTitle!=""?plotTitle:selectionNames[selType]); assert(hs);
    hs->Add(histos[ kPlotVZbb ] );   
    hs->Add(histos[ kPlotVVLF ] );   
    hs->Add(histos[ kPlotZLF  ] );    
    hs->Add(histos[ kPlotZb   ] );   
    hs->Add(histos[ kPlotZbb  ] );   
    hs->Add(histos[ kPlotWLF  ] );   
    hs->Add(histos[ kPlotWb   ] );   
    hs->Add(histos[ kPlotWbb  ] );   
    hs->Add(histos[ kPlotTop  ] );   
    hs->Add(histos[ kPlotTT   ] );   
    if(plotQCD) hs->Add(histos[ kPlotQCD  ] );   
    if(stackSignal && primarySignal==kPlotWH) hs->Add(histos[kPlotWH]);
    if(stackSignal && primarySignal==kPlotZH) hs->Add(histos[kPlotZH]);
 
    TCanvas *canvas = new TCanvas(Form("canvas_%s",theHistoName.Data()),theHistoName,600,480);
    TPad *pad1=0,*pad2=0;
    if(doRatioPad) {
      canvas->SetTopMargin(0.0); 
      canvas->SetRightMargin(0.02);
      canvas->SetBottomMargin(0); 
      pad1= new TPad("pad1", "pad1", 0, 0.3, 1.,1);
      pad1->SetTopMargin(0.1);
      pad1->SetLeftMargin(0.15);
      pad1->SetRightMargin(0.04);
      pad1->SetBottomMargin(0.03); 
      //pad1->SetGridx();         
      pad1->Draw();             
      pad1->cd();              
      if(doLogPlot) pad1->SetLogy();
    }
    hs->Draw("HIST");
    if(doRatioPad) {
      hs->GetXaxis()->SetTitle("");
      hs->GetXaxis()->SetLabelSize(0);
      hs->GetYaxis()->SetTitleOffset(1.2);
      hs->GetYaxis()->SetTitleSize(0.05);
      hs->GetYaxis()->SetLabelSize(0.05);
    } else {
      hs->GetXaxis()->SetTitle(xlabel);
    }
    
    //float binWidth=(xmax-xmin)/(float)nbins;
    
    hs->GetYaxis()->SetTitle(Form( (binWidth>=10.? "Events / %.0f %s" : binWidth>=1.? "Events / %.1f %s": "Events / %.2f %s"), binWidth, units.c_str()));
    float plotMax=1.4;
    if(hTotalBkg->GetMean() > xmin + (xmax-xmin)/4.) plotMax=2.;
    float theMax = TMath::Max(histos[primarySignal]->GetMaximum(),TMath::Max(hs->GetMaximum(), histos[kPlotData]->GetMaximum()));
    if(doLogPlot) { 
      //if(hTotalBkg->GetBinContent(hTotalBkg->GetMinimumBin())>1 && histos[kPlotData]->GetBinContent(histos[kPlotData]->GetMinimumBin())>1)
      if(hTotalBkg->GetBinContent(hTotalBkg->GetMinimumBin())>1 || histos[kPlotData]->GetBinContent(histos[kPlotData]->GetMinimumBin())>1)
        hs->SetMinimum(hTotalBkg->GetBinContent(hTotalBkg->GetMinimumBin())/2.);
      else hs->SetMinimum(0.1);
      float span=hTotalBkg->GetBinContent(hTotalBkg->GetMaximumBin())-hTotalBkg->GetBinContent(hTotalBkg->GetMinimumBin());
      //float span=hTotalBkg->GetBinContent(hTotalBkg->GetMaximumBin());
      hs->SetMaximum(TMath::Max((float)50.0,float(hTotalBkg->GetBinContent(hTotalBkg->GetMaximumBin())+pow(span,1.5))));
    } else {
      hs->SetMaximum( plotMax*theMax);
    }
    if(!stackSignal) histos[primarySignal]->Draw("HIST SAME");
    hErrorBand->Draw("E2 same");
    //if(!isBlinded) histos[kPlotData]->Draw("P E0 SAME");
    histos[kPlotData]->Draw("P E0 SAME");
    TLegend *legend1,*legend2;
    float x1=.59,x2=.77,x3=.95;
    if(!doRatioPad) { x1=.4; x2=.64; x3=.88;}
    if(plotQCD) legend1=new TLegend(x1,0.50 ,x2,0.88);
    else        legend1=new TLegend(x1,0.555,x2,0.88); 
    legend2=new TLegend(x2,0.555,x3,.88);
    legend1->AddEntry(histos[ kPlotData ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotData)] ,"lp");
    if(plotQCD) legend1->AddEntry(histos[ kPlotQCD  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotQCD )] ,"f");
    legend1->AddEntry(histos[ kPlotTT   ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotTT  )] ,"f");
    legend1->AddEntry(histos[ kPlotTop  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotTop )] ,"f");
    legend1->AddEntry(histos[ kPlotWbb  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotWbb )] ,"f");
    legend1->AddEntry(histos[ kPlotWb   ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotWb  )] ,"f");
    legend1->AddEntry(histos[ kPlotWLF  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotWLF )] ,"f");
    legend2->AddEntry(histos[ kPlotZbb  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotZbb )] ,"f");
    legend2->AddEntry(histos[ kPlotZb   ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotZb  )] ,"f");
    legend2->AddEntry(histos[ kPlotZLF  ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotZLF )] ,"f");
    legend2->AddEntry(histos[ kPlotVVLF ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotVVLF)] ,"f");
    legend2->AddEntry(histos[ kPlotVZbb ], vhbbPlot::plotNames[static_cast<plotCategory>(kPlotVZbb)] ,"f");
    legend2->AddEntry(histos[primarySignal], histos[primarySignal]->GetName(), stackSignal?"f":"lp");
    legend1->SetFillColorAlpha(kWhite, .5); legend2->SetFillColorAlpha(kWhite, 0.5);
    legend1->SetBorderSize(0); legend2->SetBorderSize(0);
    legend1->Draw("same"); legend2->Draw("SAME");
    if(doRatioPad) {
      canvas->cd();
      pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
      pad2->SetLeftMargin(0.15);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(0.4);
      pad2->SetRightMargin(0.04);
      //pad2->SetGridx(); 
      pad2->Draw(); pad2->cd();
      hRatioBand->SetTitle("");
      hRatioBand->GetXaxis()->SetTitleSize(0.16);
      hRatioBand->GetXaxis()->SetTitleOffset(1.05);
      hRatioBand->GetYaxis()->SetNdivisions(503);
      hRatioBand->GetXaxis()->SetLabelSize(0.2);
      hRatioBand->GetXaxis()->SetTitle(xlabel);
      hRatioBand->GetYaxis()->SetTitle("Data/MC");
      hRatioBand->GetYaxis()->CenterTitle();
      hRatioBand->GetYaxis()->SetTitleOffset(0.37);
      hRatioBand->GetYaxis()->SetTitleSize(0.12);
      hRatioBand->GetYaxis()->SetLabelSize(0.15);
      hRatioBand->SetFillStyle(3254);
      hRatioBand->SetFillColor(kBlack);
      hRatioBand->SetLineColor(kBlack);
      hRatioBand->SetMinimum(.3);
      hRatioBand->SetMaximum(1.7);
      hRatioBand->Draw("E2");
      if(hRatioBandPrefit) {
        hRatioBandPrefit->SetMarkerColor(kBlack);
        hRatioBandPrefit->SetMarkerSize(0);
        hRatioBandPrefit->SetLineColor(kBlack);
        hRatioBandPrefit->SetFillColorAlpha(kBlue-9,0.5);
        hRatioBandPrefit->SetFillStyle(1001);
        hRatioBandPrefit->Draw("e2 same");
        hRatioBand->Draw("E2 same");
      }
      
      hRatio->SetLineColor(kBlack);
      hRatio->SetMarkerStyle(20);
      hRatio->SetMarkerSize(0.8);
      hRatio->Draw("P E1 x0 SAME");
      TLine *baseline = new TLine(xmin,1,xmax,1);
      baseline->SetLineStyle(kSolid); baseline->Draw("SAME");
      TLegend *ratioLegend=0;
      if(hRatioBandPrefit) {
        ratioLegend = new TLegend(.17,.86,.46,.97);
        ratioLegend->SetNColumns(2);
        ratioLegend->SetFillColorAlpha(kWhite, .5);
        ratioLegend->SetBorderSize(0);
        ratioLegend->AddEntry(hRatioBandPrefit, "Prefit stat.+syst.","f");
        ratioLegend->AddEntry(hRatioBand, "Postfit stat.+syst.","f");
        ratioLegend->Draw("same");
      }
    }
    canvas->Print(outPdf);
    system(Form("gs -sDEVICE=png16m -dTextAlphaBits=4 -g1800x1440 -dUseCropBox -dFIXEDMEDIA -dPDFFitPage -o %s %s >/dev/null 2>&1 &",outPng.Data(), outPdf.Data()));
    delete canvas;
  }


  return;
}

