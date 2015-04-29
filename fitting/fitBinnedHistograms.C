#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

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
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooChi2Var.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

// Root stuff
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"

#include "TMath.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"


//bkgtable
#include "bkgTableFSR.h"
#include "dataTable.h"
// miscellaneous  
#include <fstream>
#include <new>
#include <iostream>
double mass_l =  8.8;
double mass_h = 10.;
double binw   = 0.01;    //bin width of the histogram

const int nData   = 2;
const char* choseSampleLegend[nData] = {"",
					"CMS Simulation"}; //CMS simulation, pp #sqrt{s} = 2.76 TeV

const char* choseSampleLumi[nData] = {"",
				      ""};

// -------- make some local picks
bool doMinos      = true;     //kFALSE;
bool chisquare = false; // default false: do NLL instead. works as long as you're not in extended likelihood mode.
using namespace std;
using namespace RooFit;
using namespace RooStats;

void fitBinnedHistograms(int choseSample    = 1, //Input data sample. 1: pyquen 1S, 2: 2Spyquen
			 int choseFitParams = 0, //0: (1s, 2s, 3s) 1: (1s, 2s/1s; 3s/1s); 2: (1S, (2s+3s)/1s) //dropped
			 int bkgdModel      = 0, //0: Nothing; 1:LS erf*exp + pol2; 2:LS RookeyPdf + pol2; 3:erf*exp; 4:pol2; 5:erf*exp+pol2 6:pol3
			 int fixFSR         = 0,
			 int fixSigma1      = 0, // 0 free; 1: fix
			 int centralityMin = 0,
			 int centralityMax = 40,
			 float muonEtaMin  = -1.6,
			 float muonEtaMax  = 1.6, 
			 float dimuPtMin    = 0., 
			 float dimuPtMax    = 100,
			 double muonpTcut1  = 3.5, //single muon pT cut
			 double muonpTcut2 = 4, //1 should always be lower than 2!
			 bool plotBkg      = 0, //0: hide LS or trkRot; 1: plot LS or trkRot data points and fit lines;
			 bool doTrkRot     = 0, //0: use LS;   1: use track rotation
			 bool doConstrainFit   = 0,  //1: use constrain method
			 int useRef            = 0, // # 0 none 1: data, 2: MC, 3: old
			 bool plotpars         = 1, //1: plot parameters;   0: plot CMS label
			 const char* choseSampleCase = "Pyquen", //
			 const char* outFigsDir      = "pdfOutput/Syst_sig/08102014_pyquen/",// figs dir outfile location
			 TString outDatsDir          = "txtOutput",// dats dir outfile location
			 const char* outFilePrefix   = "MCpars", // yields, ratios
			 bool narrowMass             = false,
			 float dimuYMin = 0.0,
			 float dimuYMax = 100,
			 bool isHI =1,
			 int signalModel =4,
			 int doRap=1,
			 int doPt=1){
  int doCent=0;
  if(!doPt && !doRap) doCent=1;

  cout << "hi,"<<endl;
  cout<< "check the variable names!! " << endl;
  cout << "you're tryin to bin, ya sucker!"<<endl; 
  cout<< "-----------------------------------"<<endl; 
  cout << "don't get spastic now, chill!"<<endl;
  gROOT->Macro("cm/logon.C+");
  // input files
  std::string finput;
  double muonPtCut_min1 = muonpTcut1;
  double muonPtCut_min2 = muonpTcut2;
 switch (choseSample) 
    {
    case 0:
      finput ="";
      cout << "you don't pick any tree! Fuck I'm gonna crash fo' sure..." << endl;
      break;
    case 1:
      finput = "~/Project/ups2013/upsiMiniTree_Pyquen1S_QQtrigbit1_Trig_Unfolding_postCut_deltaRmatched_withCentrality.root";
      cout << "picked the 1S embedded sample, filtered for only passing events" << endl;
      break;
    default: break;
    }

    bool scanLL =true;
    TString figsDir(Form("%s",outFigsDir)); //output fig location
  
    // if param are on, this gets filled in the naming
    TString paramOn_("");
    // kinematic cuts:
    double muonEtaCut_min = muonEtaMin;
    double muonEtaCut_max = muonEtaMax; 
    double upsYCut_min    = -2.4;
    double upsYCut_max    = 2.4;
    int whatBin;
    if(doRap){
      whatBin=8;
      if(dimuYMin ==0. && dimuYMax==0.4)    {whatBin=0;}	
      if(dimuYMin ==0.4 && dimuYMax==0.8)    {whatBin=1;}	
      if(dimuYMin ==0.8 && dimuYMax==1.2)    {whatBin=2;}	
      if(dimuYMin ==1.2 && dimuYMax==1.6)    {whatBin=3;}	
      if(dimuYMin ==1.6 && dimuYMax==2.0)     {whatBin=4;}
      if(dimuYMin ==2.0 && dimuYMax==2.4)     {whatBin=5;}
      if(dimuYMin ==0.0 && dimuYMax==1.2)     {whatBin=6;}
      if(dimuYMin==1.2 && dimuYMax==2.4)      {whatBin=7;}
    }
    if(doPt){
      whatBin=9;
      if(dimuPtMin == 0. && dimuPtMax ==2.5)   {whatBin=0;}	
      if(dimuPtMin == 2.5 && dimuPtMax ==5.)    {whatBin=1;}	
      if(dimuPtMin == 5. && dimuPtMax == 8.)    {whatBin=2;}	
      if(dimuPtMin == 8 && dimuPtMax ==12)      {whatBin=3;}	
      if(dimuPtMin == 12 && dimuPtMax ==20)      {whatBin=4;}
      if(dimuPtMin == 20 && dimuPtMax ==50)      {whatBin=5;}
      if(dimuPtMin == 0. && dimuPtMax ==6.5)   {whatBin=6;}
      if(dimuPtMin == 6.5 && dimuPtMax ==10)     {whatBin=7;}
      if(dimuPtMin == 10 && dimuPtMax ==20)  {whatBin=8;}
    }
    if(doCent){
      if(centralityMin==0 && centralityMax==2) {whatBin=0;}
      if(centralityMin==2 && centralityMax==4) {whatBin=1;}
      if(centralityMin==4 && centralityMax==8) {whatBin=2;}
      if(centralityMin==8 && centralityMax==12) {whatBin=3;}
      if(centralityMin==12 && centralityMax==16) {whatBin=4;}
      if(centralityMin==16 && centralityMax==20) {whatBin=5;}
      if(centralityMin==20 && centralityMax==28) {whatBin=6;}
      if(centralityMin==28 && centralityMax==40) {whatBin=7;}
      if(centralityMin==20 && centralityMax==40) {whatBin=13;}
    }
    cout<< whatBin << endl;
    cout << "been there" << endl;
    int centrality_max = centralityMax;  ///  
    int centrality_min = centralityMin; 
    if(choseSample!=1){ TString cut_ap("");}
    else { TString cut_ap(Form("  (%.2f<RecoMuPlusEta && RecoMuPlusEta < %.2f) && (%.2f<RecoMuMinusEta && RecoMuMinusEta < %.2f) && (%.2f<RecoUpsPt && RecoUpsPt<%.2f)  && (abs(RecoUpsRap)<%.2f && abs(RecoUpsRap)>%.2f)  &&((RecoMuPlusPt > %.2f && RecoMuMinusPt > %.2f) || (RecoMuPlusPt > %.2f && RecoMuMinusPt > %.2f))",muonEtaCut_min,muonEtaCut_max,muonEtaCut_min,muonEtaCut_max,dimuPtMin,dimuPtMax,dimuYMax,dimuYMin,muonPtCut_min1,muonPtCut_min2,muonPtCut_min2,muonPtCut_min1));}
    

    int centMin = centrality_min;
    int centMax = centrality_max;
    if (centralityMin==28) binw=0.05;
    
    if(choseSample==1) 
      {
	centMin = (int)(centrality_min*2.5);
	centMax = (int)(centrality_max*2.5);
      }
    TString figName_(Form("%s_%s_cent%d-%d_bkgModel%d_sigModel%d_muonEta%.2f%.2f_muonPt%.2f-%.2f_dimuPt%.2f%.2f_dimuY%.2f%.2f_trkRot%d_constrain%d_fsr%d_sigma%d_ref%d",
			  outFilePrefix,choseSampleCase,centMin,centMax,
			  bkgdModel, signalModel , muonEtaCut_min,muonEtaCut_max,muonPtCut_min1,muonPtCut_min2,dimuPtMin,dimuPtMax,dimuYMin,dimuYMax,doTrkRot,doConstrainFit,fixFSR,fixSigma1,useRef));
    figName_.ReplaceAll("-","M");
    figName_.ReplaceAll(".","");
    
    cout<<"Fitting: Pt["<< dimuPtMin <<","<< dimuPtMax <<"] , muPt > ["<< muonPtCut_min1<<","<<muonPtCut_min2<<"] and centrality ["<<centrality_min<<","<<centrality_max<<"]!!!!"<<endl;
    cout << "oniafitter processing"
	 << "\n\tInput:  \t" << finput
	 << "\n\tOutput: \t" << figName_
	 << endl;
    // -----------    
    TFile *f = new TFile(finput.c_str());
    TTree* theTree       = (TTree*)gROOT->FindObject("UpsilonTree"); // OS --- all mass
    RooRealVar* mass       = new RooRealVar("RecoUpsM","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");
    RooRealVar* upsPt      = new RooRealVar("RecoUpsPt","p_{T}(#Upsilon)",0,60,"GeV");
    //  RooRealVar* upsEta     = new RooRealVar("upsEta",  "upsEta"  ,-10,10);
    RooRealVar* upsRapidity= new RooRealVar("RecoUpsRap",  "upsRapidity",-1000, 1000);
    RooRealVar* vProb      = new RooRealVar("_VtxProb",  "vProb"  ,0.01,1.00);
    //   RooRealVar* QQsign     = new RooRealVar("QQsign",  "QQsign"  ,-1,5);
    RooRealVar* Centrality = new RooRealVar("centrality","centrality",0,100);
    RooRealVar* GenUpsPt= new RooRealVar("GenUpsPt","p_{T}(#Upsilon^{GEN})",0,200,"GeV");
    RooRealVar* GenMuPlusPt   = new RooRealVar("GenMuPlusPt","GenMuPlusPt",muonPtCut_min1,100);
    RooRealVar* GenMuMinusPt  = new RooRealVar("GenMuMinusPt","GenMuMinusPt",muonPtCut_min1,100);
    RooRealVar* GenMuPlusEta  = new RooRealVar("GenMuPlusEta","GenMuPlusEta",  -2.4,2.4);
    RooRealVar* GenMuMinusEta = new RooRealVar("GenMuMinusEta","GenMuMinusEta",-2.4,2.4);
    RooRealVar* muPlusPt   = new RooRealVar("RecoMuPlusPt","RecoMuPlusPt",muonPtCut_min1,100);
    RooRealVar* muMinusPt  = new RooRealVar("RecoMuMinusPt","RecoMuMinusPt",muonPtCut_min1,100);
    RooRealVar* muPlusEta  = new RooRealVar("RecoMuPlusEta","RecoMuPlusEta",  -2.4,2.4);
    RooRealVar* muMinusEta = new RooRealVar("RecoMuMinusEta","RecoMuMinusEta",-2.4,2.4);
    RooDataSet *data0 = new RooDataSet("data0","data0",theTree,
				       RooArgSet(*mass,*upsPt,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta,*upsRapidity,*GenUpsPt,*Centrality));
    RooDataSet *redData;
    // RooDataSet *tmp;
    RooFormulaVar wFunc("w","event weight","GenUpsPt>30?0:(GenUpsPt>15.0?(0.0877039/107560.0):(GenUpsPt>12?(0.0977067/117270.0):(GenUpsPt>9?(0.169418/162777.0):(GenUpsPt>6?(0.430737/172165.0):(GenUpsPt>3?(1.1129191/171048.0):(1.0/172208.0))))))",*GenUpsPt);
    RooRealVar *w = (RooRealVar*) data0->addColumn(wFunc);

    data0->Print();
    redData = new RooDataSet(data0->GetName(),data0->GetTitle(),data0,*data0->get(),cut_ap,w->GetName());
    // for weighting 
    redData->SetName("redData");
    cout << "koko..." << endl;
    redData->Print();
    // RooDataSet *data = ( RooDataSet*)redData;
    // data->SetName("data");
    // data->Print("v");

    // *************************************************** signal PDF
    const double M1S = 9.46;   //upsilon 1S pgd mass value
    const double M2S = 10.023;  //upsilon 2S pgd mass value
    const double M3S = 10.355;  //upsilon 3S pgd mass value

    // *************************************************** free param in the fit
    int nt = 10000;
    int nt = redData->sumEntries();
    if(isHI==1 && choseSample==1){  RooRealVar *nsig1f   = new RooRealVar("N_{#Upsilon(1S)}","nsig1S",0,nt);
    }else{ RooRealVar *nsig1f   = new RooRealVar("N_{#Upsilon(1S)}","nsig1S",0,nt*10);
      RooRealVar *nsig2f  = new RooRealVar("N_{#Upsilon(2S)}","nsig2S",   nt*0.25,-1*nt,10*nt);
      RooRealVar *nsig3f  = new RooRealVar("N_{#Upsilon(3S)}","nsig3S",   nt*0.25,-1*nt,10*nt);
    } 
    RooRealVar  *mean = new RooRealVar("mass1S","#Upsilon mean",M1S,M1S-0.1,M1S+0.1);
    RooConstVar *rat2 = new RooConstVar("rat2", "rat2", M2S/M1S);
    RooConstVar *rat3 = new RooConstVar("rat3", "rat3", M3S/M1S);
    // scale mean and resolution by mass ratio
    RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",RooArgList(*mean));
    RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*@1", RooArgList(*mean,*rat2));
    RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0*@1", RooArgList(*mean,*rat3));
    
    //detector resolution ?? where is this coming from?
    RooRealVar    *sigma1  = new RooRealVar("#sigma_{CB1}","#sigma_{CB1}",sigma_min[whatBin],sigma_max[whatBin]); // 
    RooFormulaVar *sigma1S = new RooFormulaVar("sigma1S","@0"   ,RooArgList(*sigma1));
    RooFormulaVar *sigma2S = new RooFormulaVar("sigma2S","@0*@1",RooArgList(*sigma1,*rat2));
    RooFormulaVar *sigma3S = new RooFormulaVar("sigma3S","@0*@1",RooArgList(*sigma1,*rat3));
    RooRealVar *alpha  = new RooRealVar("#alpha_{CB}","tail shift",alpha_min[whatBin],alpha_max[whatBin]);    // MC 5tev 1S pol2 
    RooRealVar *npow   = new RooRealVar("npow","power order",npow_min[whatBin],npow_max[whatBin]);    // MC 5tev 1S pol2 
    RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.,1.);
    // scale the sigmaGaus with sigma1S*scale=sigmaGaus now.
    RooRealVar    *scaleWidth = new RooRealVar("#sigma_{CB2}/#sigma_{CB1}","scaleWidth",1.,2.7);
    RooFormulaVar *sigmaGaus = new RooFormulaVar("sigmaGaus","@0*@1", RooArgList(*sigma1,*scaleWidth));
    RooFormulaVar *sigmaGaus2 = new RooFormulaVar("sigmaGaus","@0*@1*@2", RooArgList(*sigma1,*scaleWidth,*rat2));
    RooFormulaVar *sigmaGaus3 = new RooFormulaVar("sigmaGaus","@0*@1*@2", RooArgList(*sigma1,*scaleWidth,*rat3));
    RooGaussian* gauss1 = new RooGaussian("gaus1s","gaus1s",
					  *nsig1f,
					  *mass,    //mean
					  *sigmaGaus); //sigma
    RooGaussian* gauss1b = new RooGaussian("gaus1sb","gaus1sb",
					   *nsig1f,
					   *mass,    //mean
					   *sigma1); //sigma
    cout << signalModel << "^th model" << endl;
    switch(signalModel){    
    case 1: //crystal boule
      RooCBShape  *sig1S   = new RooCBShape ("cb1S_1", "FSR cb 1s",
					     *mass,*mean1S,*sigma1,*alpha,*npow);
      break;
    case 2: //Gaussein
      RooAbsPdf      *sig1S  = new RooGaussian ("g1", "gaus 1s",
						*mass,*mean1S,*sigma1);
      break;
    case 3: //Gaussein + crystal boule
      RooCBShape  *cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
					       *mass,*mean1S,*sigma1,*alpha,*npow);
      RooAddPdf      *sig1S  = new RooAddPdf ("cbg", "cbgaus 1s",
					      RooArgList(*gauss1,*cb1S_1),*sigmaFraction);
      break;
    case 4: //crystal boules
      RooCBShape  *cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
					       *mass,*mean1S,*sigma1,*alpha,*npow);
       
      RooCBShape  *cb1S_2    = new RooCBShape ("cb1S_2", "FSR cb 1s",
					       *mass,*mean1S,*sigmaGaus,*alpha,*npow);
      RooAddPdf      *sig1S  = new RooAddPdf  ("cbcb","1S mass pdf",
					       RooArgList(*cb1S_1,*cb1S_2),*sigmaFraction);
      // /// Upsilon 2S
      RooCBShape  *cb2S_1    = new RooCBShape ("cb2S_1", "FSR cb 2s", 
					       *mass,*mean2S,*sigma2S,*alpha,*npow); 
      RooCBShape  *cb2S_2    = new RooCBShape ("cb2S_2", "FSR cb 2s", 
					       *mass,*mean2S,*sigmaGaus2,*alpha,*npow); 
      RooAddPdf      *sig2S  = new RooAddPdf  ("sig2S","2S mass pdf",
					       RooArgList(*cb2S_1,*cb2S_2),*sigmaFraction);
      
      // /// Upsilon 3S
      RooCBShape  *cb3S_1    = new RooCBShape ("cb3S_1", "FSR cb 3s", 
					       *mass,*mean3S,*sigma3S,*alpha,*npow); 
      RooCBShape  *cb3S_2    = new RooCBShape ("cb3S_2", "FSR cb 3s", 
					       *mass,*mean3S,*sigmaGaus3,*alpha,*npow); 
      RooAddPdf      *sig3S  = new RooAddPdf  ("sig3S","3S mass pdf",
					       RooArgList(*cb3S_1,*cb3S_2),*sigmaFraction); // = cb3S1*sigmaFrac + cb3S2*(1-sigmaFrac)
      break;
      
    case 5: //deux Gausseins
      RooAddPdf      *sig1S  = new RooAddPdf ("cb1S_1", "cbgaus 1s",
					      RooArgList(*gauss1,*gauss1b),*sigmaFraction);
         
      break;
    }
    // bkg Chebychev
    RooRealVar *nbkgd   = new RooRealVar("n_{Bkgd}","nbkgd",0,nt);
    RooRealVar *bkg_a1  = new RooRealVar("a1_bkg", "bkg_{a1}", 0, -5, 5);
    RooRealVar *bkg_a2  = new RooRealVar("a2_Bkg", "bkg_{a2}", 0, -5, 5);
    RooRealVar *bkg_a3  = new RooRealVar("a3_Bkg", "bkg_{a3}", 0, -2, 2);
    RooAbsPdf  *pdf_combinedbkgd  = new RooChebychev("bkgPdf","bkgPdf",
						     *mass, RooArgList(*bkg_a1,*bkg_a2));
    bkg_a2->setVal(0);
    bkg_a2->setConstant();
    RooRealVar *bkgdFraction = new RooRealVar("bkgdFrac","Background normalisation",0.001,0,1);
    //fitted histo
    int nbins = ceil((mass_h-mass_l)/binw); 
    TH1D *MReco;
    MReco = new TH1D("MReco","Reco di-muon mass",nbins,mass_l,mass_h);
    MReco = (TH1D*) redData->createHistogram("RecoUpsM",*mass);
    MReco->SetBinErrorOption(TH1::kPoisson);
    cout << "kokomo" << endl;
    RooDataHist binnedData ("binnedData","binnedData",*mass,Import(*MReco));
    binnedData.Print("v");  
    RooFitResult *fit_2nd; 
    RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",
						 RooArgList(*pdf_combinedbkgd,*sig1S),
    						 *bkgdFraction);
    if(chisquare){
    //chi2 of a binned dataset.
    RooChi2Var chi2("chi2","chi2",*sig1S,binnedData,SumW2Error(kTRUE)) ;
    RooMinuit m(chi2) ;
    m.migrad() ;
    m.hesse() ;
    //m.minos();
    fit_2nd= (RooFitResult*) m.save();
    }
    else{
    // NLL fit
    RooAbsReal* nll = pdf->createNLL(binnedData,NumCPU(4)) ;
    RooMinuit(*nll).migrad();
    RooMinuit(*nll).hesse();
    fit_2nd = pdf->fitTo(binnedData,Extended(0),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
    }
    //for the plots!
    RooPlot* frame = mass->frame(Bins(nbins),Range(mass_l,mass_h));
    binnedData.plotOn(frame,DataError(RooAbsData::SumW2),Name("theData"),MarkerSize(0.8));
    pdf->plotOn(frame,Components("cb1S_1"),Name("TheFirstCB"),LineColor(kTeal)); 
    pdf->plotOn(frame,Components("cb1S_2"),Name("TheSecondCB"),LineColor(kOrange)); 
    pdf->plotOn(frame,Components("bkgPdf"),Name("TheBackground"),LineColor(kDashed)); 
    pdf->plotOn(frame,Name("thePdf"),LineColor(kBlue)); 

    RooArgSet * pars = pdf->getParameters(binnedData);
    //plot data
    TCanvas c; c.cd();
    binnedData.plotOn(frame,Name("theData"),MarkerSize(0.8)); 
    pdf->plotOn(frame,Components("cb1S_1"),Name("TheFirstCB"),LineColor(kTeal)); 
    pdf->plotOn(frame,Components("cb1S_2"),Name("TheSecondCB"),LineColor(kOrange)); 
    pdf->plotOn(frame,Components("bkgPdf"),Name("TheBackground"),LineColor(kDashed)); 
    pdf->plotOn(frame,Name("thePdf"),LineColor(kBlue)); 
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    frame->GetXaxis()->CenterTitle(kTRUE);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->Draw();
    c.Draw();
    c.SaveAs(figsDir+figName_+paramOn_+".pdf");
    fit_2nd->Print("v");
    int ibin = 102;
    double err_low = MReco->GetBinErrorLow(ibin);
    double err_up = MReco->GetBinErrorUp(ibin);
    cout<<err_low <<endl;
    cout <<"crazy"<<endl;
    cout<<err_up<<endl;
    TCanvas cm("cm","cm");
    cm.cd();
    TLatex latex1;
    latex1.SetNDC();
    TPad *pPad1 = new TPad("pPad1","pPad1",0.01,0.3-0.03,0.96,0.92);
    pPad1->SetBottomMargin(0.03);
    TPad *pPad2 = new TPad("pPad2","pPad2",0.05,0.05,1,0.3);
    pPad2->SetTopMargin(0.0);
    pPad2->SetBottomMargin(-0.1);
    frame->SetMinimum(0.00001);
    pPad1->SetLogy();
    // pPad2->SetBottomMargin(gStyle->GetPadBottomMargin()/0.3);
    // pPad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
    pPad1->Draw();
    pPad1->cd();
    //  sig1S->paramOn(frame,Layout(0.12,0.5,0.38));
    pPad1->Update();
    frame->Draw();
    latex1.SetTextSize(0.032);
    latex1.DrawLatex(0.15,1.-0.05*1.8,Form("%s",choseSampleLegend[choseSample]));
    latex1.DrawLatex(0.15,1.-0.05*4.5,Form("%.2f < |y| < %.2f",dimuYMin,dimuYMax)); 
    latex1.DrawLatex(0.15,1.-0.05*5.5,Form("%.1f < p_{T}^{#Upsilon} < %.1f",dimuPtMin,dimuPtMax));
    latex1.DrawLatex(0.15,1.-0.05*6.5,Form("p_{T}^{#mu1} > %.1f GeV/c",muonPtCut_min1));
    latex1.DrawLatex(0.15,1.-0.05*7.5,Form("p_{T}^{#mu2} > %.1f GeV/c",muonPtCut_min2));

    latex1.DrawLatex(0.78,1.-0.05*3.5,Form("n_{CB} = %.2f",npow->getVal()));
    latex1.DrawLatex(0.78,1.-0.05*4.5,Form("#alpha_{CB} = %.2f",alpha->getVal()));
    latex1.DrawLatex(0.78,1.-0.05*5.5,Form("#sigma_{CB1} = %.2f",sigma1->getVal()));
    latex1.DrawLatex(0.78,1.-0.05*6.5,Form("#sigma_{CB2}/#sigma_{CB1} = %.2f",scaleWidth->getVal()));
    latex1.DrawLatex(0.72,1.-0.05*7.5,Form("normalisation = %.2f",sigmaFraction->getVal()));
    latex1.DrawLatex(0.78,1.-0.05*8.5,Form("m_{0} = %.3f",mean1S->getVal()));
    latex1.DrawLatex(0.78,1.-0.05*9.5,Form("Bkg_{frac} = %.3f",bkgdFraction->getVal()));
    cm.cd(0);
    TPad *pPad2 = new TPad("pPad2","pPad2",0.01,0.05,0.96,0.29);
    pPad2->SetTopMargin(0.0);
    pPad2->SetBottomMargin(-0.1);
    pPad2->Draw();
    pPad2->cd();
    double chi2FromRoo = frame->chiSquare(fit_2nd->floatParsFinal().getSize());
    cout<<"!!!!!!!! chi2 from simple pull= "<<frame->chiSquare()<<"\t chi2 from RooFit= "<<chi2FromRoo <<endl;
    RooHist *phPullm = frame->pullHist(0,0,true); // this calcualtes the pulls taking the integral of the fit in each bin, instead of the value in the middle of the bid
    phPullm->SetName("phPullm");
    double *ypull     = phPullm->GetY();
    double Chi2       = 0;
    int nFullBinsPull = 0;
    for (int i=0; i < nbins; i++) 
      {
        if (MReco->GetBinContent(i) == 0) continue;
        nFullBinsPull++;
        Chi2 = Chi2 + pow(ypull[i],2);
      }
    // for writing on canvas
    int nFitParam     = fit_2nd->floatParsFinal().getSize();
    int Dof           = nFullBinsPull - nFitParam;
    double UnNormChi2 = Chi2;
    Chi2             /= (nFullBinsPull - nFitParam);
    
    cout<<"!!!!! nFullBinsPull="<<nFullBinsPull<<"\tnFitParam="<<nFitParam<<endl;
    // draw pulls
    pPad2->cd();
    double mFrameMax = 0;
  

    RooPlot* prpFramePull = mass->frame(Title("Pull"),Bins(nbins),Range(mass_l,mass_h));
    
    prpFramePull->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    prpFramePull->GetXaxis()->CenterTitle(kTRUE);
    prpFramePull->GetXaxis()->SetTitleSize(0.06);
    prpFramePull->GetXaxis()->SetLabelSize(0.1);
    prpFramePull->GetYaxis()->CenterTitle(kTRUE);
    prpFramePull->GetYaxis()->SetTitleSize(0.08);
    prpFramePull->GetYaxis()->SetLabelSize(0.1);
    prpFramePull->GetYaxis()->SetTitleOffset(0.4);
    prpFramePull->GetXaxis()->SetTitleOffset(0.6);
   
    prpFramePull->GetYaxis()->SetTitle("Pull");
    prpFramePull->addPlotable(phPullm,"PX");
 
    if (prpFramePull->GetMinimum()*-1 > prpFramePull->GetMaximum()) mFrameMax = prpFramePull->GetMinimum()*-1;
    else mFrameMax = prpFramePull->GetMaximum();
    prpFramePull->SetMaximum(mFrameMax); 
    prpFramePull->SetMinimum(-1*mFrameMax); 
    prpFramePull->Draw();

    latex1.SetTextSize(0.085);
    double myChi2 = chi2FromRoo*Dof;
    latex1.DrawLatex(0.7,1.-0.05*3.5,Form("#chi^{2}/ndf = %2.1f/%d",myChi2,Dof));
   
    //cm.SaveAs(figsDir+figName_+paramOn_+"_pulls.png");
    cm.SaveAs(figsDir+figName_+paramOn_+"_pulls.pdf");
    float baseNll = fit_2nd->minNll();
    float estimatedDistance2Minimum = fit_2nd->edm();
    string outParameters_forNote = outDatsDir+"/"+figName_+"forNote.txt";
    cout<<"Output file: " << outParameters_forNote<<endl;
    ofstream outfileFitResults_forNote;
    outfileFitResults_forNote.open(outParameters_forNote.c_str(), ios_base::out);
    outfileFitResults_forNote<<figName_<<" "<<nsig1f->getVal()<<" "<<nsig1f->getError()<<" "<<npow->getVal()<<" "<<npow->getError()<<" "<<alpha->getVal()<<" "<<alpha->getError()<<" "<<sigma1->getVal()<<" "<<sigma1->getError()<<" "<<scaleWidth->getVal()<<" "<<scaleWidth->getError()<<" "<<sigmaFraction->getVal()<< " "<<sigmaFraction->getError()<<" "<<estimatedDistance2Minimum<<" "<<baseNll<<endl;
    outfileFitResults_forNote.close();
    RooWorkspace *ws = new RooWorkspace("ws","workspace with good stuff");
    ws->import(binnedData);
    ws->import(*fit_2nd);
    ws->import(*pdf);
    
    
    string outRooFitResult = outDatsDir+"/WS_"+figName_+"_fitres.root"; 
    ws->writeToFile(Form("%s",outRooFitResult.c_str())); 
}
