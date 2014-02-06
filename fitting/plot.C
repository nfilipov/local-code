#if !defined(__CINT__) || defined(__MAKECINT__)
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <TROOT.h>

#include <TH1.h>
#include <TH2D.h>

#include <TBranch.h>
#include <TCanvas.h>
#include "TClonesArray.h"
#include <TDirectory.h>
#include <TFile.h>
#include "TH1.h"
#include <TLatex.h>
#include <TLegend.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "TRandom.h"
#include <TStyle.h>
#include "TPaveLabel.h"
#include "TVirtualFitter.h"
#include <TSystem.h>
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TCut.h"
// miscellaneous  
#include <fstream>
#include <map>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooKeysPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#endif

using namespace RooFit;
using namespace RooStats;
using namespace std;

void plot(
	    ){
  gROOT->Macro("../code/cm/logon.C+");
  const int nDataTrees = 5;
  const int nMCTrees = 4;
  const int nBinsPtRap = 2*6-1;
  const char* outFigsDir      = "test_pdfOutput/nonFit/";
  int nTotalTrees = nDataTrees+nMCTrees;
  int _sampleMC=1;
  int _sampleData=1;
  TString inFileData[nDataTrees] = {"",
				    "../dimuonTree_upsiMiniTree_aa276tev_regitreco_glbglb_Runa_trigBit1_allTriggers0_pt4.root",
				    "../dimuonTree_upsiMiniTree_AA2p76tev_ptmuSpecial_nov25_2013_trigBit1_allTriggers1_testNoCut.root",//dummy
				    "../dimuonTree_upsiMiniTree_pp276tev_5p41_Run211739-211831_trigBit1_allTriggers0_pt4.root", // pp trk trk
				    "../dimuonTree_upsiMiniTree_pp276tev_5p41_ptmu2_Run211739-211831_trigBit1_allTriggers0.root"// pp GlbGlb
  };
  TString inFileMC[nMCTrees] = {"",
				"../dimuonTree_upsiMiniTree_1Spythia2p76_ptmu4_Nov06_2013_RunMC_trigBit1_allTriggers1.root",
				"../dimuonTree_upsiMiniTree_2Spythia2p76_ptmu4_Nov06_2013_RunMC_trigBit1_allTriggers1.root",
				"../dimuonTree_upsiMiniTree_3Spythia2p76_ptmu4_Nov13_2013_RunMC_trigBit1_allTriggers1.root"
  };
  bool BoldCut=1;
  bool EasyCut=0;
  bool Overlay=1;
  bool doAll=1;
  bool doMC =0;
  bool do3p5;
  bool do4;
  if(BoldCut){
    do3p5=true;
    do4=false;
  }
	gROOT->Macro("MCParameters.h");

  char* _sample;
 
  // write out the fitting params MC, data
  //PICK A FILE : MC, DATA
  if(doMC){
    std::string outParametersMC = "testOutMC.txt";
    cout<<"Output file: " << outParametersMC<<endl;
    ofstream outfileFitResultsMC;
    outfileFitResultsMC.open(outParametersMC.c_str(), ios_base::out | ios_base::app);
    TString finput_mc = inFileMC[_sampleMC];
    TFile f(finput_mc,"read"); 
    switch (_sampleMC)
      {
      case 0:
	cout<< "no MC sample! Pick from 1 to 3."<<endl;
	break;
      case 1:
	_sample = "#varUpsilon(1S) [Pythia]";
	break;
      case 2:
	_sample = "#varUpsilon(2S) [Pythia]";
	break;
      case 3:
	_sample = "#varUpsilon(3S) [Pythia]";
	break;
      default:
	cout << "pick a Data sample!"<<endl;
	break;
      }
  }
  if(!doMC) 
    { std::string outParameters = "testOut.txt";
      cout<<"Output file: " << outParameters<<endl;
      ofstream outfileFitResults;
      outfileFitResults.open(outParameters.c_str(), ios_base::out | ios_base::app);
      TString finput = inFileData[_sampleData];
      TFile f(finput, "read");
      switch (_sampleData)
	{
	case 0:
	  cout<< "no Data sample! Pick from 1 to 4."<<endl;
	  break;
	case 1:
	  _sample = "PbPb, GG RegIT";
	  break;
	case 2:
	  _sample = "PbPb, noCut RegIT";
	  break;
	case 3:
	  _sample = "pp Trk-Trk";
	  break;
	case 4:
	  _sample ="pp Glb-Glb";
	  break;
	default:
	  cout << "pick a Data sample!"<<endl;
	  break;
	}
    }
 

  double mass_l = 8.5 ;
  if(doAll)
    {
  double mass_h = 11.5;			     
    }
  else if(!doAll)
    {
      double mass_h = 9.8;
    }
if (doMC)
  {
    double mass_l = 8.5 ;
    double mass_h = 10.1;
  }

//fixed parameters from MC;
 if(do3p5){
       float massRes1_MC_pt[5]={0.0702278,
     			       0.0704769,
     			       0.0723508,
     			       0.072539,
     			       0.0741795};
       float massRes2_MC_rap[5]={0.0987689,
     			    	0.10944,
     			    	0.12396,
     			    	0.149925,
     			    	0.181642};
       float massRes2_MC_pt[5]={0.135703,
     			    	0.139709,
     			    	0.140769,
     			    	0.141848,
     			    	0.14674};
       float massRes1_MC_rap[5]={0.055119,
     			    	 0.0654371,
     			    	 0.0789018,
     			    	 0.0947856,
     			    	 0.113839};
        float npowMC_pt[5]={1.63403,
     			   1.52982,
     			   1.53018,
     			   1.3078,
     			   1.35339};
       float npowMC_rap[5]={1.13925,
     			    1.22514,
     			    1.20287,
     			    1.29962,
     			    1.23529};
       float alphaMC_pt[5]={ 1.75943,
     			     1.87361,
     			     1.82563,
     			     1.93345,
     			     1.93388};
       float alphaMC_rap[5]={1.98463,
     			     1.97219,
     			     2.01558,
     			     1.9974,
     			     2.06177};
       float sigmaFraction_rap[5]={ 0.10950, 0.12203, 0.13997,  0.15214, 0.28186};
       float sigmaFraction_pt[5]={0.68306,0.6819,0.70786,0.70281,0.73227};
 }
 
 
  if (!doMC)  {
    TTree *UpsilonTree = (TTree*)gROOT->FindObject("UpsilonTree");
  } // OS --- all mass
  else if (doMC) {
    TTree *UpsilonTree = (TTree*)gROOT->FindObject("RecoUpsilonTree");
  }// OS --- all mass
  RooRealVar* mass       = new RooRealVar("invariantMass","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");
  RooRealVar* upsPt      = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,60,"GeV");
  RooRealVar* upsEta     = new RooRealVar("upsEta",  "upsEta"  ,-10,10);
  RooRealVar* upsRapidity= new RooRealVar("upsRapidity",  "upsRapidity",-2.4, 2.4);
  RooRealVar* vProb      = new RooRealVar("vProb",  "vProb"  ,0.01,1.00);  
  RooRealVar* muPlusPt   = new RooRealVar("muPlusPt","muPlusPt",0,100);
  RooRealVar* muMinusPt  = new RooRealVar("muMinusPt","muMinusPt",0,100);
  RooRealVar* muPlusEta  = new RooRealVar("muPlusEta","muPlusEta",  -2.4,2.4);
  RooRealVar* muMinusEta = new RooRealVar("muMinusEta","muMinusEta",-2.4,2.4);

  RooDataSet *data0 = new RooDataSet("data0","data0",UpsilonTree,
				     RooArgSet(*mass,*upsPt,*upsRapidity,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta));
  data0->Print();
  //VARIOUS CUTS
  TCut muPlus_10_006_barrel(Form("(muPlusPt>3.4 && (muPlusEta<1 && muPlusEta>-1))"));
  TCut muPlus_10_006_int(Form("(((muPlusEta>1 && muPlusEta <1.6)||(muPlusEta>-1.6 && muPlusEta<-1.)) && muPlusPt>3)"));
  TCut muPlus_10_006_endcap(Form("(((muPlusEta>1.6 && muPlusEta <2.4)||(muPlusEta>-2.4 && muPlusEta<-1.6)) && muPlusPt>2.5 )"));
 
  TCut muMinus_10_006_barrel(Form("(muMinusPt>3.4 && (muMinusEta<1 && muMinusEta>-1))"));
  TCut muMinus_10_006_int(Form("(((muMinusEta>1 && muMinusEta <1.6)||(muMinusEta>-1.6 && muMinusEta<-1.)) && muMinusPt>3)"));
  TCut muMinus_10_006_endcap(Form("(((muMinusEta>1.6 && muMinusEta <2.4)||(muMinusEta>-2.4 && muMinusEta<-1.6)) && muMinusPt>2)"));
 
  TCut muPlus_10_006(muPlus_10_006_barrel+muPlus_10_006_int+muPlus_10_006_endcap); 
  TCut muMinus_10_006(muMinus_10_006_barrel+muMinus_10_006_int+muMinus_10_006_endcap);  // already applied on tree!
 
  TCut muPlus_11_011(Form("(muPlusEta<2.4 && muPlusEta>-2.4) && muPlusPt>4"));
  TCut muMinus_11_011(Form("(muMinusEta<2.4 && muMinusEta>-2.4) && muMinusPt>4"));
  double binw = 0.05;
  int nbins = ceil((mass_h-mass_l)/binw); 
  float vProbCut = 0.01;

  // TString mass_vProb(Form("invariantMass<%.3f&&invariantMass>%.3f&&vProb>%.3f",mass_h,mass_l,vProbCut)); 
  if (doMC)
    {  TString mass_vProb(Form("invariantMass<%.3f&&invariantMass>%.3f",mass_h,mass_l));
    } else {TString mass_vProb(Form("invariantMass<%.3f&&invariantMass>%.3f&&vProb>%.3f",mass_h,mass_l,vProbCut)); }
 
  TCut acc_aa_10_006(muMinus_10_006*muPlus_10_006*mass_vProb);
  TCut acc_aa_11_011((muMinus_11_011*muPlus_11_011)*mass_vProb);

  std::string rapCuts[6]= {"",
			   "(abs(upsRapidity)<0.4) && ",
			   "(abs(upsRapidity)>0.4&&abs(upsRapidity)<0.7) && ",
			   "(abs(upsRapidity)>0.7&&abs(upsRapidity)<1.) && ",
			   "(abs(upsRapidity)>1.&&abs(upsRapidity)<1.5) && ",
			   "(abs(upsRapidity)>1.5&&abs(upsRapidity)<2.4) &&"};
  std::string _binSuffix[nBinsPtRap] = {"",
				   "_maxRap0p4",
				   "_maxRap0p7",
				   "_maxRap1",
				   "_maxRap1p5",
				   "_maxRap2p4",
				   "_maxPt2p5",
				   "_maxPt5",
				   "_maxPt8",
				   "_maxPt12",
				   "_maxPt20",
  };
 char _plotSuffix[nBinsPtRap][1000] = {"","#varUpsilon |y| < 0.4","#varUpsilon 0.4 < |y| < 0.7","#varUpsilon 0.7 < |y| < 1.0","#varUpsilon 1.0 < |y| < 1.5","#varUpsilon 1.5 < |y| < 2.4","p_{T}^{#mu#mu} < 2.5","2.5 < p_{T}^{#mu#mu} < 5","5 < p_{T}^{#mu#mu} < 8","8 < p_{T}^{#mu#mu} < 12","12 < p_{T}^{#mu#mu} < 20"}; 
  std::string ptCuts[6]= {"(upsPt<2.5)&& ",
			  "(upsPt>=2.5&&upsPt<5.)&& ",
			  "(upsPt>=5.&&upsPt<8.)&& ",
			  "(upsPt>=8.&&upsPt<12.)&& ",
			  "(upsPt>=12.&&upsPt<20.)&&",
			  ""};
 TString KinCut="";
  //ROOFIT VARS
  //FIT VARIABLES
  const double M1S = 9.46;   //upsilon 1S pgd mass value
  const double M2S = 10.023;  //upsilon 2S pgd mass value
  const double M3S = 10.355;  //upsilon 3S pgd mass value
   int nEntries = data0->sumEntries();
  RooRealVar  *mean = new RooRealVar("mass1S","#Upsilon mean",M1S-0.1,M1S+0.1);
  RooConstVar *rat2 = new RooConstVar("rat2", "rat2", M2S/M1S);
  RooConstVar *rat3 = new RooConstVar("rat3", "rat3", M3S/M1S);

  // scale mean and resolution by mass ratio
  RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",RooArgList(*mean));
  RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*@1", RooArgList(*mean,*rat2));
  RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0*@1", RooArgList(*mean,*rat3));

  //detector resolution ?? where is this coming from?
  RooRealVar    *sigma1  = new RooRealVar("sigma1","#sigma_{1S}",0.08,0.01,0.2); // 
  RooRealVar    *sigmaGaus = new RooRealVar("sigmaGaus","#sigmaGaus_{1S}",0.08,0.01,0.2); 
  RooFormulaVar *sigma1S = new RooFormulaVar("sigma1S","@0"   ,RooArgList(*sigma1));
  RooFormulaVar *sigma2S = new RooFormulaVar("sigma2S","@0*@1",RooArgList(*sigma1,*rat2));
  RooFormulaVar *sigma3S = new RooFormulaVar("sigma3S","@0*@1",RooArgList(*sigma1,*rat3));
  
  /// to describe final state radiation tail on the left of the peaks
  RooRealVar *alpha  = new RooRealVar("alpha","tail shift",1.2,0.01,4);    // MC 5tev 1S pol2 
  RooRealVar *npow   = new RooRealVar("npow","power order",2.3,1.,20);    // MC 5tev 1S pol2 

  // relative fraction of the two Gaussians components for each CB
  RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.1,0,1.);
  if(!doMC)
    {
  sigmaFraction->setVal(0);
  sigmaFraction->setConstant(kTRUE);
    }
      //ROOFIT PDFDEF
  RooGaussian *g1S_1 = new RooGaussian("g1S_1","gaussian 1S",
				       *mass,*mean1S,*sigmaGaus);
  RooCBShape  *cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
					   *mass,*mean1S,*sigma1,*alpha,*npow);
  RooCBShape  *cb1S_2    = new RooCBShape ("cb1S_2", "FSR cb 1s",
					   *mass,*mean1S,*sigmaGaus,*alpha,*npow);
  if(!doMC)
    {
  RooCBShape  *cb1S_2    = new RooCBShape ("cb1S_2", "FSR cb 1s",
					   *mass,*mean1S,*sigma1,*alpha,*npow);
    }
  RooAddPdf      *sig1S  = new RooAddPdf  ("sig1S","1S mass pdf",
					   RooArgList(*cb1S_1,*cb1S_2),*sigmaFraction);

  // bkg Chebychev
  RooRealVar *nbkgd   = new RooRealVar("n_{Bkgd}","nbkgd",0,nEntries*4);
  RooRealVar *bkg_a1  = new RooRealVar("a1_bkg", "bkg_{a1}", 0, -1.2, 1.2);
  RooRealVar *bkg_a2  = new RooRealVar("a2_Bkg", "bkg_{a2}", 0, -0.8, 0.8);
  RooRealVar *bkg_a3  = new RooRealVar("a3_Bkg", "bkg_{a3}", 0, -0.5, 0.5);

 

  // sig1S is then jsut cb1S2: c*pdf_1+(1-c)*pdf_2

  /// Upsilon 2S
  RooCBShape  *cb2S_1    = new RooCBShape ("cb2S_1", "FSR cb 2s", 
					   *mass,*mean2S,*sigma1,*alpha,*npow); 
  RooCBShape  *cb2S_2    = new RooCBShape ("cb2S_2", "FSR cb 2s", 
					   *mass,*mean2S,*sigma2S,*alpha,*npow); 
  RooAddPdf      *sig2S  = new RooAddPdf  ("sig2S","2S mass pdf",
					   RooArgList(*cb2S_1,*cb2S_2),*sigmaFraction);
  
  /// Upsilon 3S
  RooCBShape  *cb3S_1    = new RooCBShape ("cb3S_1", "FSR cb 3s", 
					   *mass,*mean3S,*sigma1,*alpha,*npow); 
  RooCBShape  *cb3S_2    = new RooCBShape ("cb3S_2", "FSR cb 3s", 
					   *mass,*mean3S,*sigma3S,*alpha,*npow); 
  RooAddPdf      *sig3S  = new RooAddPdf  ("sig3S","3S mass pdf",
					   RooArgList(*cb3S_1,*cb3S_2),*sigmaFraction); // = cb3S1*sigmaFrac + cb3S2*(1-sigmaFrac)
 
  RooRealVar *nsig1f   = new RooRealVar("N_{#Upsilon(1S)}","nsig1S",10,nEntries);
  // if(doMC) 
  //   {nsig1f->setVal(nEntries);
  //     nsig1f->setConstant(kTRUE);}
  RooRealVar *nsig2f   = new RooRealVar("N_{#Upsilon(2S)}","nsig2S",1,nEntries);
  RooRealVar *nsig3f   = new RooRealVar("N_{#Upsilon(3S)}","nsig3S",1,nEntries);

  RooFitResult* fit;// fit results
  //containers for signal and bkg estimates
  double gNB=5;
  double sig_mc[2*6]; sig_mc[0]=0;
  double sig_mcErr[2*6]; sig_mcErr[0]=0;
  double mean_mc[2*6]; 
  double mean_mcErr[2*6];
  double sigma_mc[2*6];
  double sigma_mcErr[2*6];
  double npow_mc[2*6];
  double npow_mcErr[2*6];
  double alpha_mc[2*6];
  double alpha_mcErr[2*6];

  double sig_data[2*6];
  double sig_dataErr[2*6];
  double bkg_data[2*6];
  double bkg_dataErr[2*6];
  double mean_data[2*6]; 
  double mean_dataErr[2*6];
  double sigma_data[2*6];
  double sigma_dataErr[2*6];
  // double npow_data[2*6];
  // double npow_dataErr[2*6];
  // double alpha_data[2*6];
  // double alpha_dataErr[2*6];
  //variables for plotting
  float pt[5]={1.25,3.75,6.75,10,16};
  float pte[5]={1.25,1.25,1.5,2,4};
  float rap[5]={0.2,0.55,0.85,1.25,1.95};
  float rape[5]={0.2,0.15,0.15,0.25,0.45};
  float sig_mcRap[5]={};
  float sig_mcRapErr[5]={};
  float sig_mcPt[5]={};
  float sig_mcPtErr[5]={};
  float sig_dataRap[5]={};
  float sig_dataRapErr[5]={};
  float sig_dataPt[5]={};
  float sig_dataPtErr[5]={};
  float bkg_dataPt[2*6]={};
  float bkg_dataPtErr[2*6]={};
  float bkg_dataRap[2*6]={};
  float bkg_dataRapErr[2*6]={};
  float yieldOverErr[nBinsPtRap]={};
  float yieldOverErr_e[nBinsPtRap]={};
  float ptYieldOverErr[5]={};
  float ptYieldOverErr_e[5]={};
  float rapYieldOverErr[5]={};
  float rapYieldOverErr_e[5]={};
  
  //LOOP

  for(int iLoop = 1 ; iLoop < nBinsPtRap ; iLoop++)
    {
      if(BoldCut)
	{
	  do3p5=true;
	}
         gROOT->Macro("MCParameters.h");
      if(iLoop<6)
	{ KinCut = rapCuts[iLoop];
	  if(!doMC){
	    npow->setVal(npowMC_rap[iLoop-1]);npow->setConstant(kTRUE);  
	    alpha->setVal(alphaMC_rap[iLoop-1]);alpha->setConstant(kTRUE);
	    sigma1->setVal(sigmaFraction_rap[iLoop-1]*massRes2_MC_rap[iLoop-1]+(1-sigmaFraction_rap[iLoop-1])*massRes1_MC_rap[iLoop-1]);
	    sigma1->setConstant(kTRUE);
	    cout<< npowMC_rap[iLoop-1]  << " " << alphaMC_rap[iLoop-1]  << " " << sigmaFraction_rap[iLoop-1]*massRes2_MC_rap[iLoop-1]+(1-sigmaFraction_rap[iLoop-1])*massRes1_MC_rap[iLoop-1] << endl;
	    }
	}
      else if(iLoop>=6)
		{
		  KinCut = ptCuts[iLoop-6];
		  if(!doMC){
		    npow->setVal(npowMC_pt[iLoop-6]);npow->setConstant(kTRUE);
		    alpha->setVal(alphaMC_pt[iLoop-6]);alpha->setConstant(kTRUE);
		    sigma1->setVal(sigmaFraction_pt[iLoop-6]*massRes1_MC_pt[iLoop-6]+(1-sigmaFraction_pt[iLoop-6])*massRes2_MC_pt[iLoop-6]);
		    sigma1->setConstant(kTRUE);
		  }
		}
   
      float muonPtLeg1=3.5;
      float muonPtLeg2=4.0;
      TString bSuffix = _binSuffix[iLoop];
      double muonEtaMax =2.4;
      
      if (BoldCut)
	{KinCut+=" (abs(muPlusEta)<2.4 && abs(muMinusEta)<2.4) && ((muMinusPt>3.5 && muPlusPt>4.0) || (muMinusPt>4. && muPlusPt>3.5))";}
      else if(EasyCut) 
	{ KinCut+=" (abs(muPlusEta)<2.4 && abs(muMinusEta)<2.4) && ((muMinusPt>4 && muPlusPt>4))";} 
      TString cut_app= KinCut;
      cout << "entering fitting loop : " << KinCut << endl;
      RooDataSet *data0_ap = ( RooDataSet*)data0->reduce(Cut(cut_app));
      cout << "yeah?" << endl;
      data0_ap->SetName("data0_ap");
      cout<< "FUCK!!!"<<endl;
      data0_ap->Print();
      RooFitResult* fit;// fit results
      RooFitResult* efit;
      RooFitResult* efit2;
      std::string cName = "testCanvas";
      std::string cTitle ="test canvas";
      cName+=_binSuffix[iLoop];
      cTitle+=_binSuffix[iLoop];
      cout<< "Drawing..." << endl;
      TCanvas *c_uno = new TCanvas(cName.c_str(),cTitle.c_str(),10,10,1010,510);
      c_uno->Divide(2,1);
      c_uno->cd(1);
      RooPlot* frame = mass->frame(Bins(nbins));
      frame->SetTitle("");
      frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
      frame->GetXaxis()->CenterTitle(kTRUE);
      frame->GetYaxis()->SetTitleOffset(1.3);
      data0_ap->plotOn(frame,Name("theData"),MarkerSize(0.8));
      if(!doMC)
      	{
      	  sigmaFraction->setVal(0);
      	  sigmaFraction->setConstant(kTRUE);
      	  // mass->setRange("sb_left",8.81,9.31);  // change "natural" a few lines below"
      	  // mass->setRange("sb_right",9.61,10.11);
	  
      	  RooRealVar nbkgw_sb("nbkgw_sb","nBkgd in side-band window",0,nEntries);
      	  mass->setRange("natural",mass_l,mass_h);
      	  RooAbsPdf  *pdf_combinedbkgd  = new RooChebychev("bkgPdf","bkgPdf",
      	  						   *mass, RooArgList(*bkg_a1,*bkg_a2,*bkg_a3));	
      	  if(doAll)
      	    {
      	      RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",
      							   RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
      							   RooArgList(*nsig1f,*nsig2f,*nsig3f,nbkgw_sb));
      	      efit = pdf->fitTo(*data0_ap,Save(kTRUE),NumCPU(4),Extended(kTRUE));
      	    }
      	  else if(!doAll)
      	    {
      	      RooExtendPdf ebkgd("ebkgd","ebkgd",*pdf_combinedbkgd,nbkgw_sb,"natural");
	      RooExtendPdf esig("esig","esig",*g1S_1,*nsig1f,"natural");
      	      RooAddPdf  pdf_data("pdf_data","total p.d.f.",
      				  RooArgList(esig,ebkgd)); 
      	      efit = pdf_data.fitTo(*data0_ap,Save(kTRUE),NumCPU(4));
      	    }
	}
	  if (doMC){
	    mass->setRange("massWindow",mass_l,mass_h);
	    
	    //	RooRealVar nsigw("nsigw","nsignal in window",0,nEntries*5.);
	    RooAbsPdf  *pdf_mc             = new RooAddPdf ("pdf","total p.d.f.",
							    RooArgList(*sig1S),
							    RooArgList(*nsig1f)); 
            
	    fit     = pdf_mc->fitTo(*data0_ap,Save(kTRUE),Minos(),NumCPU(4));
	    
	    // RooArgSet * pars = pdf->getParameters(data0_ap);
	    
	    sig_mc[iLoop]     = nsig1f->getVal();
	    sig_mcErr[iLoop]  = nsig1f->getError(); 
	    mean_mc[iLoop]    = mean->getVal();
	    mean_mcErr[iLoop] = mean->getError();
	    sigma_mc[iLoop]   = sigma1->getVal();
	    sigma_mcErr[iLoop]= sigma1->getError();
	    npow_mc[iLoop]    = npow->getVal();
	    npow_mcErr[iLoop] = npow->getError();
	    alpha_mc[iLoop]   = alpha->getVal();
	    alpha_mcErr[iLoop]= alpha->getError();
	    if (iLoop<6)
	      { sig_mcRap[iLoop-1]=sig_mc[iLoop];
		sig_mcRap[iLoop-1]=sig_mcErr[iLoop];
		
	      }
	    else if(iLoop>=6)
	      { sig_mcPt[iLoop-6]=sig_mc[iLoop];
		sig_mcPtErr[iLoop-6]=sig_mc[iLoop];
		
	      }  
	    double signalStrength = nsig1f->getVal()/sqrt(nsig1f->getError());
	
	    cout <<" there?" <<endl;
	    pdf_mc->plotOn(frame,Components("cb1S_2"),Name("theFitPDF_mca"),LineColor(kOrange));
	    pdf_mc->plotOn(frame,Components("cb1S_1"),Name("theFitPDF_mcb"),LineColor(kGray));
	    pdf_mc->plotOn(frame,Name("theFitPDF_mc"),LineColor(kBlue)); 
	    double chi2FromRoo = frame->chiSquare(fit->floatParsFinal().getSize());
	    data0_ap->plotOn(frame,Name("theData"),MarkerSize(0.8));
	    cout<<"Writing to file..."<<endl;
	    outfileFitResultsMC <<" "<< finput_mc <<" "<< _binSuffix[iLoop] <<", pt1:"<< muonPtLeg1 << ", pt2:" << muonPtLeg2 << " - "<< sig_mc[iLoop] <<" "<< sig_mcErr[iLoop] <<" "<< mean_mc[iLoop] <<" "<< mean_mcErr[iLoop] <<" "<< sigma_mc[iLoop] <<" "<< sigma_mcErr[iLoop] <<" "<<sigmaGaus->getVal()<<" "<<sigmaGaus->getError()<<" "<< npow_mc[iLoop] <<" "<< npow_mcErr[iLoop] <<" "<< alpha_mc[iLoop] <<" "<< alpha_mcErr[iLoop] <<" "<< fit->minNll() << endl;
	    cout<<"Writing to file...!"<< endl;
	    fit->Print();
	    frame->Draw(); 
	    TLatex latex1;
	    latex1.SetNDC();
	    latex1.SetTextSize(0.03);
	    latex1.DrawLatex(0.15,1.-0.1*1,Form("p_{T}^{#mu(1)} = %.1f,p_{T}^{#mu(2)} = %.1f",muonPtLeg1,muonPtLeg2));
	    latex1.DrawLatex(0.7,1.-0.1*1,Form("%s",_plotSuffix[iLoop]));
	    latex1.DrawLatex(0.15,1.-0.1*1.5,Form("N_{#varUpsilon(1S)} = %.1f #pm %.1f ",nsig1f->getVal(),nsig1f->getError()));
	    latex1.DrawLatex(0.15,1.-0.1*3,Form("n_{CB} = %.2f #pm %.2f",npow->getVal(),npow->getError()));
	    latex1.DrawLatex(0.15,1.-0.1*3.5,Form("#alpha = %.2f #pm %.2f",alpha->getVal(),alpha->getError()));
	    latex1.DrawLatex(0.15,1.-0.1*2,Form("m_{#varUpsilon(1S)} = %.3f #pm %.3f",mean->getVal(),mean->getError()));
	    latex1.DrawLatex(0.15,1.-0.1*2.5,Form("#sigma_{#varUpsilon(1S)}(MeV) = %.1f #pm %.1f",1000*(sigma1->getVal()),1000*(sigma1->getError())));
	    latex1.DrawLatex(0.15,1.-0.1*4,Form("#chi^{2} = %.2f",chi2FromRoo));
	    latex1.DrawLatex(0.7,1.-0.1*1.5,Form("%s",_sample));
	    cout<<"Plot Drawn!"<<endl;
	  }
	  if(!doMC)
	    {
	      	  sig_data[iLoop] = nsig1f->getVal();
	      	  sig_dataErr[iLoop] = nsig1f->getError(); 
	      	  mean_data[iLoop]    = mean->getVal();
	      	  mean_dataErr[iLoop] = mean->getError();
	      	  sigma_data[iLoop]   = sigma1->getVal();
	      	  sigma_dataErr[iLoop]= sigma1->getError();
	      	  bkg_data[iLoop] = nbkgw_sb.getVal();
	      	  bkg_dataErr[iLoop] = nbkgw_sb.getError();
	      	  yieldOverErr[iLoop]= sig_data[iLoop]/sig_dataErr[iLoop];
	      	  yieldOverErr_e[iLoop]=0;
	      	  if (iLoop<6)
	      	    { sig_dataRap[iLoop-1]=sig_data[iLoop];
	      	      sig_dataRapErr[iLoop-1]=sig_dataErr[iLoop];
	      	      bkg_dataRap[iLoop-1]=bkg_data[iLoop];
	      	      bkg_dataRapErr[iLoop-1]=bkg_data[iLoop]; 
	      	      rapYieldOverErr[iLoop-1]=yieldOverErr[iLoop];
	      	      rapYieldOverErr_e[iLoop-1]=yieldOverErr_e[iLoop];
	      	    } else if(iLoop>=6)
	      	    { sig_dataPt[iLoop-6]=sig_data[iLoop];
	      	      sig_dataPtErr[iLoop-6]=sig_data[iLoop]; 
	      	      bkg_dataPt[iLoop-6]=bkg_data[iLoop];
	      	      bkg_dataPtErr[iLoop-6]=bkg_data[iLoop];
	      	      ptYieldOverErr[iLoop-6]=yieldOverErr[iLoop];
	      	      ptYieldOverErr_e[iLoop-6]=yieldOverErr_e[iLoop];
	      	    }
	      	  if(!doAll)
	      	    {
		      pdf_data.plotOn(frame,Components("ebkgd"),Name("theFitPDFbkg"),LineColor(kBlue),DrawOption("LF"),FillStyle(3004), FillColor(kBlue));  
	      	      pdf_data.plotOn(frame,Name("theFitPDF"));
		      double chi2FromRoo = frame->chiSquare(efit->floatParsFinal().getSize());
	      	    }
	      	  if(doAll)
	      	    {
		      pdf->plotOn(frame,Components("bkgPdf"),Name("theFitPDFbkg"),LineColor(kBlue),DrawOption("LF"),FillStyle(3004), FillColor(kBlue));
		      pdf->plotOn(frame,Name("theFitPDF"));
		      double chi2FromRoo = frame->chiSquare(efit->floatParsFinal().getSize());
	      	      //  pdf_data.paramOn(frame,Layout(0.15,0.7,0.9));
	      	      cout <<" there?" <<endl;
	      	    }
	      	  cout<<"Writing to file... "<<endl;
	      	  outfileFitResults<<" "<< finput<<" "<< bSuffix<<", pt1:"<< muonPtLeg1 << ", pt2:" << muonPtLeg2 << " - " << yieldOverErr[iLoop]<<" "<<sig_data[iLoop]<<" "<<sig_dataErr[iLoop]<<" "<<mean_data[iLoop]<<" "<<mean_dataErr[iLoop]<<" "<<sigma_data[iLoop]<<" "<<sigma_dataErr[iLoop]<<" "<<efit->minNll()<<endl;
	      	  cout<<"Writing to file...!"<<endl; 
	      	  efit->Print();
	      	
	      	  data0_ap->plotOn(frame,Name("theData"),MarkerSize(0.8));
	      	  frame->Draw();
	      	  TLatex latex1;
	      	  latex1.SetNDC();
	      	  latex1.SetTextSize(0.032);
	      	  latex1.DrawLatex(0.15,1.-0.1*1,Form("p_{T}^{#mu(1)} = %.1f,p_{T}^{#mu(2)} = %.1f",muonPtLeg1,muonPtLeg2));
	      	  latex1.DrawLatex(0.7,1.-0.1*1,Form("%s",_plotSuffix[iLoop]));
		  if(doMC){
	      	  latex1.DrawLatex(0.15,1.-0.1*3,Form("n_{CB} = %.2f #pm %.2f",npow->getVal(),npow->getError()));
	      	  latex1.DrawLatex(0.15,1.-0.1*3.5,Form("#alpha = %.2f #pm %.2f",alpha->getVal(),alpha->getError()));
		  }
	      	  latex1.DrawLatex(0.15,1.-0.1*2,Form("m_{#varUpsilon(1S)} = %.3f #pm %.3f",mean->getVal(),mean->getError()));
	      	  latex1.DrawLatex(0.15,1.-0.1*2.5,Form("#sigma_{#varUpsilon(1S)}(MeV) = %.1f #pm %.1f",1000*(sigma1->getVal()),1000*(sigma1->getError())));
	      	  latex1.DrawLatex(0.15,1.-0.1*4,Form("#chi^{2} = %.2f",chi2FromRoo));
		  latex1.DrawLatex(0.7,1.-0.1*1.5,Form("%s",_sample));
	      	  latex1.DrawLatex(0.15,1.-0.1*1.5,Form("N_{#varUpsilon(1S)} = %.1f #pm %.1f, #frac{N}{err} = %.1f ",nsig1f->getVal(),nsig1f->getError(),yieldOverErr[iLoop]));
	      	  latex1.DrawLatex(0.15,1.-0.1*3,Form("N_{Bkg} = %.1f #pm %.1f",nbkgw_sb->getVal(),nbkgw_sb->getError()));
	    } 	
    
	  
	  cout<<"Overlaying..."<<endl;
	  if(Overlay){
	    do4 = true;
	    do3p5 = false;
	    gROOT->Macro("MCParameters.h");
	    if(do4){
	      float sigmaFraction_rap[5]={0.07687, 0.097598, 0.169 , 0.12379, 0.27256};
	      float sigmaFraction_pt[5]={0.66387,0.67881,0.723,0.69545,0.71832};
	      float massRes1_MC_pt[5]={0.0683712,
				       0.0694958,
				       0.0717929,
				       0.0713796,
				       0.0734393};
	      float massRes2_MC_rap[5]={0.107213,
					0.112244,
					0.12054,
					0.155516,
					0.18437};
	      float massRes2_MC_pt[5]={0.137299,
				       0.140186,
				       0.142773,
				       0.140888,
				       0.144339};
	      float massRes1_MC_rap[5]={0.0552849,
					0.0651953,
					0.077617,
					0.0962206,
					0.113659};
	      float npowMC_pt[5]={2.07882,
				  1.65775,
				  1.73509,
				  1.44071,
				  1.37867};
	      float npowMC_rap[5]={1.25553,
				   1.34141,
				   1.22328,
				   1.5064,
				   1.48334};
	      float alphaMC_pt[5]={1.77283,
				   1.8701,
				   1.75824,
				   1.87573,
				   1.92024};
	      float alphaMC_rap[5]={1.96474,
				    1.96184,
				    2.06546,
				    1.97356,
				    2.05739};
	    }
	    if(iLoop<6)
	      { KinCut = rapCuts[iLoop];
		if(!doMC){
	      npow->setVal(npowMC_rap[iLoop-1]);npow->setConstant(kTRUE);
	      alpha->setVal(alphaMC_rap[iLoop-1]);alpha->setConstant(kTRUE);
	      sigma1->setVal(sigmaFraction_rap[iLoop-1]*massRes2_MC_rap[iLoop-1]+(1-sigmaFraction_rap[iLoop-1])*massRes1_MC_rap[iLoop-1]);
	      sigma1->setConstant(kTRUE);
	    }
	  }
	else if(iLoop>=6)
	  { KinCut = ptCuts[iLoop-6];
	    if(!doMC){
	      npow->setVal(npowMC_pt[iLoop-6]);npow->setConstant(kTRUE);
	      alpha->setVal(alphaMC_pt[iLoop-6]);alpha->setConstant(kTRUE);
	      sigma1->setVal(sigmaFraction_pt[iLoop-6]*massRes1_MC_pt[iLoop-6]+(1-sigmaFraction_pt[iLoop-6])*massRes2_MC_pt[iLoop-6]);
	      sigma1->setConstant(kTRUE);
	    }
	  }
	
	TString bSuffix = _binSuffix[iLoop];
	double muonEtaMax =2.4;

	float muonPtLeg1_mem = muonPtLeg1;
	float muonPtLeg2_mem = muonPtLeg2;
	muonPtLeg1=4.0;
	muonPtLeg2=4.0;
	
	if (BoldCut)
	  {KinCut+=" (abs(muPlusEta)<2.4 && abs(muMinusEta)<2.4) && ((muMinusPt>4 && muPlusPt>4))";}
	else if(EasyCut) 
	  { KinCut+=" (abs(muPlusEta)<2.4 && abs(muMinusEta)<2.4) && ((muMinusPt>3 && muPlusPt>4) || (muMinusPt>4 && muPlusPt>3))";}  
	TString cut_app= KinCut;
	cout << "entering fitting loop : " << KinCut << endl;
	cout << muonPtLeg1 << "," << muonPtLeg2 <<endl; 
	RooDataSet *data0_ap2 = ( RooDataSet*)data0->reduce(Cut(cut_app));
	cout << "yeah?" << endl;
	data0_ap2->SetName("data0_ap2");
	cout<< "FUCK!!!"<<endl;
	data0_ap2->Print();
	RooFitResult* efit2;
	RooFitResult* fit2;
	c_uno->cd(2);
	//   gPad->SetLogy();
	RooPlot* frame2= mass->frame(Bins(nbins));
	frame2->SetTitle("");
	//frame2->GetYaxis()->SetTitle("Yield, current cut: pt > 4");
	frame2->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	frame2->GetXaxis()->CenterTitle(kTRUE);
	frame2->GetYaxis()->SetTitleOffset(1.3);
	data0_ap2->plotOn(frame2,Name("theData2"),MarkerSize(0.8));
		if(!doMC)
		  {	    
		    RooRealVar nbkgw_sb2("nbkgw_sb2","nBkgd2",0,nEntries/2);
		    mass->setRange("natural",mass_l,mass_h);
	    
		    RooAbsPdf  *pdf_combinedbkgd2  = new RooChebychev("bkgPdf2","bkgPdf2",
								      *mass, RooArgList(*bkg_a1,*bkg_a2,*bkg_a3));	
		    RooExtendPdf ebkgd2("ebkgd2","ebkgd2",*pdf_combinedbkgd2,nbkgw_sb2,"natural");
	    
		    RooExtendPdf esig2("esig2","esig2",*g1S_1,*nsig1f,"natural");
		    if(doAll)
		      {
			RooAbsPdf  *pdf2             = new RooAddPdf ("pdf","total p.d.f.",
								     RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd2),
								     RooArgList(*nsig1f,*nsig2f,*nsig3f,nbkgw_sb2));
		        efit2 = pdf2->fitTo(*data0_ap2,Save(kTRUE),NumCPU(4),Extended(kTRUE));
		      }
		    else if(!doAll)
		      {
			RooExtendPdf ebkgd2("ebkgd2","ebkgd2",*pdf_combinedbkgd2,nbkgw_sb2,"natural");
		
			RooExtendPdf esig2("esig2","esig2",*g1S_1,*nsig1f,"natural");
		    	RooAddPdf  pdf_data2("pdf_data2","total p.d.f.",
		    			     RooArgList(esig2,ebkgd2)); 
			efit2 = pdf_data2.fitTo(*data0_ap2,Save(kTRUE),NumCPU(4));
		      }
		  sig_data[iLoop] = nsig1f->getVal();
		  sig_dataErr[iLoop] = nsig1f->getError(); 
		  mean_data[iLoop]    = mean->getVal();
		  mean_dataErr[iLoop] = mean->getError();
		  sigma_data[iLoop]   = sigma1->getVal();
		  sigma_dataErr[iLoop]= sigma1->getError();
		  bkg_data[iLoop] = nbkgw_sb2.getVal();
		  bkg_dataErr[iLoop] = nbkgw_sb2.getError();
		  yieldOverErr[iLoop]= sig_data[iLoop]/sig_dataErr[iLoop];
		  yieldOverErr_e[iLoop]=0;
		  if (iLoop<6)
		    { sig_dataRap[iLoop-1]=sig_data[iLoop];
		      sig_dataRapErr[iLoop-1]=sig_dataErr[iLoop];
		      bkg_dataRap[iLoop-1]=bkg_data[iLoop];
		      bkg_dataRapErr[iLoop-1]=bkg_data[iLoop]; 
		      rapYieldOverErr[iLoop-1]=yieldOverErr[iLoop];
		      rapYieldOverErr_e[iLoop-1]=yieldOverErr_e[iLoop];
		    } else if(iLoop>=6)
		    { sig_dataPt[iLoop-6]=sig_data[iLoop];
		      sig_dataPtErr[iLoop-6]=sig_data[iLoop]; 
		      bkg_dataPt[iLoop-6]=bkg_data[iLoop];
		      bkg_dataPtErr[iLoop-6]=bkg_dataErr[iLoop];
		      ptYieldOverErr[iLoop-6]=yieldOverErr[iLoop];
		      ptYieldOverErr_e[iLoop-6]=yieldOverErr_e[iLoop];
		    }
		  if(!doAll)
		    {
		      pdf_data2.plotOn(frame2,Components("ebkgd2"),Name("theFit2PDFbkg"),LineColor(kGray),DrawOption("LF"),FillStyle(3005), FillColor(kGray));
		      pdf_data2.plotOn(frame2,Name("theFit2PDF"),LineColor(kGray));
		    }
		  else if(doAll)
		    {	 
		      pdf2->plotOn(frame2,Components("bkgPdf2"),Name("theFit2PDFbkg"),LineColor(kGray),DrawOption("LF"),FillStyle(3005), FillColor(kGray));
		      pdf2->plotOn(frame2,Name("theFit2PDF"),LineColor(kGray));
		    }
		  double chi2FromRoo2 = frame2->chiSquare(efit2->floatParsFinal().getSize());
		  outfileFitResults<<" " << finput <<" " << bSuffix<<", pt1:"<< muonPtLeg1 << " , pt2:" << muonPtLeg2 << " - " << yieldOverErr[iLoop]<<" "<<sig_data[iLoop]<<" "<<sig_dataErr[iLoop]<<" "<<mean_data[iLoop]<<" "<<mean_dataErr[iLoop]<<" "<<sigma_data[iLoop]<<" "<<sigma_dataErr[iLoop]<<" "<<bkg_data[iLoop] << " " << bkg_dataErr[iLoop] << ", "<< efit->minNll()<<endl;
		  cout<<"Writing to file...!"<<endl;
		  }
	if (doMC){
	  mass->setRange("massWindow",9.31,9.61);
	  RooAbsPdf  *pdf_mc2             = new RooAddPdf ("pdf","total p.d.f.",
							   RooArgList(*sig1S),
							   RooArgList(*nsig1f)); 
	  fit2     = pdf_mc2->fitTo(*data0_ap2,Save(kTRUE),Minos(),NumCPU(4));
	  // RooArgSet * pars = pdf->getParameters(data0_ap);
	  sig_mc[iLoop]     = nsig1f->getVal();
	  sig_mcErr[iLoop]  = nsig1f->getError(); 
	  mean_mc[iLoop]    = mean->getVal();
	  mean_mcErr[iLoop] = mean->getError();
	  sigma_mc[iLoop]   = sigma1->getVal();
	  sigma_mcErr[iLoop]= sigma1->getError();
	  npow_mc[iLoop]    = npow->getVal();
	  npow_mcErr[iLoop] = npow->getError();
	  alpha_mc[iLoop]   = alpha->getVal();
	  alpha_mcErr[iLoop]= alpha->getError();
	  if (iLoop<6)
	    { sig_mcRap[iLoop-1]=sig_mc[iLoop];
	      sig_mcRap[iLoop-1]=sig_mcErr[iLoop];
	      
	    }
	  else if(iLoop>=6)
	    { sig_mcPt[iLoop-6]=sig_mc[iLoop];
	      sig_mcPtErr[iLoop-6]=sig_mc[iLoop];
	      
	    }  
	  pdf_mc2->plotOn(frame2,Components("cb1S_2"),Name("theFit2PDF_mca"),LineColor(kOrange));
	  pdf_mc2->plotOn(frame2,Components("cb1S_1"),Name("theFit2PDF_mcb"),LineColor(kGray));
	  pdf_mc2->plotOn(frame2,Name("theFit2PDF_mc"),LineColor(kBlue));
	  double chi2FromRoo2 = frame2->chiSquare(fit2->floatParsFinal().getSize());
	  //	double signalStrength = nsig1f->getVal()/sqrt(nsig1f->getError());
	cout<<"Writing to file..."<<endl;
	outfileFitResultsMC<<" "<< _binSuffix[iLoop]<<", pt1:" << muonPtLeg1 << ", pt2:" << muonPtLeg2 << " - " << sig_mc[iLoop]<<" "<<sig_mcErr[iLoop]<<" "<<mean_mc[iLoop]<<" "<<mean_mcErr[iLoop]<<" "<<sigma_mc[iLoop]<<" "<<sigma_mcErr[iLoop]<<" "<<sigmaGaus->getVal()<<" "<<sigmaGaus->getError()<<" "<<npow_mc[iLoop]<<" "<<npow_mcErr[iLoop]<<" "<<alpha_mc[iLoop]<<" "<<alpha_mcErr[iLoop]<<" "<<fit2->minNll()<<endl;}
	  cout<<"Writing to file...!"<<endl;
	  }

      data0_ap2->plotOn(frame2,Name("theData2"),MarkerSize(0.8));
    
      frame2->Draw();
      TLatex latex2;
      latex2.SetNDC();
      latex2.SetTextSize(0.032);
      latex2.DrawLatex(0.15,1.-0.1*1,Form("p_{T}^{#mu(1)} = %.1f,p_{T}^{#mu(2)} = %.1f",muonPtLeg1,muonPtLeg2));
      latex2.DrawLatex(0.7,1.-0.1*1,Form("%s",_plotSuffix[iLoop]));
      latex2.DrawLatex(0.15,1.-0.1*2,Form("m_{#varUpsilon(1S)} = %.3f #pm %.3f",mean->getVal(),mean->getError()));
      latex2.DrawLatex(0.15,1.-0.1*2.5,Form("#sigma_{#varUpsilon(1S)}(MeV) = %.1f #pm %.1f",1000*(sigma1->getVal()),1000*(sigma1->getError())));
     
      latex2.DrawLatex(0.15,1.-0.1*4,Form("#Chi^{2} = %.2f",chi2FromRoo2));
      latex2.DrawLatex(0.15,1.-0.1*35,Form("#alpha = %.2f #pm %.2f",alpha->getVal(),alpha->getError()));
      if(!doMC) 
	{
	  latex2.DrawLatex(0.15,1.-0.1*1.5,Form("N_{#varUpsilon(1S)} = %.1f #pm %.1f, #frac{N}{err} = %.1f ",nsig1f->getVal(),nsig1f->getError(),yieldOverErr[iLoop]));
	  latex2.DrawLatex(0.7,1.-0.1*1.5,Form("%s",_sample));
	  latex2.DrawLatex(0.15,1.-0.1*3,Form("N_{Bkg} = %.1f #pm %.1f",nbkgw_sb2.getVal(),nbkgw_sb2.getError()));
	  efit2->Print();
	  if(BoldCut)
	    {
	      c_uno->SaveAs(Form("%s_%s_boldCut%d_data_Pt%.1f_wOverlay%d_%d.png",outFigsDir,cName.c_str(),BoldCut,muonPtLeg1_mem,Overlay,_sampleData));
	    }
	  else if (!BoldCut)
	    {
	      c_uno->SaveAs(Form("%s_%s_boldCut%d_data_Pt%.1f_wOverlay%d_%d.png",outFigsDir,cName.c_str(),BoldCut,muonPtLeg2,Overlay,_sampleData));
	    }
	}
      else if(doMC) 
	{
	  latex2.DrawLatex(0.15,1.-0.1*3,Form("n_{CB} = %.2f #pm %.2f",npow->getVal(),npow->getError()));
	  latex2.DrawLatex(0.7,1.-0.1*1.5,Form("%s",_sample));
	  latex2.DrawLatex(0.15,1.-0.1*1.5,Form("N_{#varUpsilon(1S)} = %.1f #pm %.1f ",nsig1f->getVal(),nsig1f->getError()));
	  fit2->Print();
			   c_uno->SaveAs(Form("%s_%s_boldCut%d_MC_Pt%.1f_wOverlay%d_%d_mc2CBs_pulls.png",outFigsDir,cName.c_str(),BoldCut,muonPtLeg1_mem,Overlay,_sampleMC));
	}
      cout<<"Plot saved "<<endl;


      c_uno->Close();
    }  
  
}


      

      ///////////PULLS DRAWN HERE/////////////////////////
      // RooHist *phPullm = frame->pullHist(0,0,true); // this calcualtes the pulls taking the integral of the fit in each bin, instead of the value in the middle of the bid
      // phPullm->SetName("phPullm");
      // double *ypull     = phPullm->GetY();

      // TPad *pPad2 = new TPad("pPad2","pPad2",0.0,0.05,0.97,0.35);
      // pPad2->SetTopMargin(0.0);
      // pPad2->cd();
      // pPad2->Draw();
      // TH1 *phData      = data0_ap->createHistogram("invariantMass",nbins);
      // double Chi2       = 0;
      // int nFullBinsPull = 0;
      // for (int i=0; i < nbins; i++) 
      // 	{
      // 	  if (phData->GetBinContent(i) == 0) continue;
      // 	  nFullBinsPull++;
      // 	  Chi2 = Chi2 + pow(ypull[i],2);
      // 	}
      // if(doMC)
      // 	{
      // 	  int nFitParam     = fit->floatParsFinal().getSize();
      // 	}
      // if(!doMC)
      // 	{
      // 	  int nFitParam     = efit->floatParsFinal().getSize();
      // 	}
      // int Dof           = nFullBinsPull - nFitParam;
      // double UnNormChi2 = Chi2;
      // Chi2             /= (nFullBinsPull - nFitParam);
      
      // cout<<"!!!!! nFullBinsPull="<<nFullBinsPull<<"\tnFitParam="<<nFitParam<<endl;
      // // draw pulls
      // pPad2->cd();
      // double mFrameMax = 0;
      // RooPlot* prpFramePull = mass->frame(Title("Pull"),Bins(nbins),Range(mass_l,mass_h));
      // prpFramePull->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
      // prpFramePull->GetXaxis()->CenterTitle(kTRUE);
      // prpFramePull->GetYaxis()->SetTitleOffset(1.3);
      // prpFramePull->GetYaxis()->SetTitle("Pull");
      // prpFramePull->addPlotable(phPullm,"PX");
      // if (prpFramePull->GetMinimum()*-1 > prpFramePull->GetMaximum()) mFrameMax = prpFramePull->GetMinimum()*-1;
      // else mFrameMax = prpFramePull->GetMaximum();
      // prpFramePull->SetMaximum(mFrameMax); 
      // prpFramePull->SetMinimum(-1*mFrameMax); 
      // prpFramePull->Draw();
    
