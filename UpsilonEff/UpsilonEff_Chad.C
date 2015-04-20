#ifndef __CINT__
#endif
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include "FitFunctions.h"
#include "TObjArray.h"

bool IsAccept(Double_t pt, Double_t eta);
double FindCenWeight(int Bin);
double RError(double A, double eA, double B, double eB);
double PError(double A, double eA, double B, double eB);

double getEffEtaPbPb(Double_t pt, Double_t eta);
double getEffEtaPP(Double_t pt, Double_t eta);


void UpsilonEff_Chad(int genacc=2, int pp=2, int YS=1, int iSpec = 2, int type =2)
{
  // int genacc =2    // 1 for Acceptance, 2 for (acc*)efficiency
  //  int pp=1;       // 1 for pp, 2 for PbPb 
  //  int YS=1;       // 1 for Y1S, 2 for Y2S, 3 for Y3S
  //  int iSpec=1;    // 1 for pT, 2 for rap, 3 for centrality 
  //  int type=1;     // 1 or 2, For Y2S and Y3S only for iSpec=1 and 2


  //See pT Cut
  double PtCut= 0.0;
  if(YS==1) PtCut=3.5;
  if(YS==2) PtCut=3.5;
  //  if(YS==2) PtCut=4.0;
  if(YS==3) PtCut=4.0;

  int isParent=1;  
  
  int Nhits = 10;
  int NPxlLayers = 0;
  double InnerChi = 4.0;
  double Dxy = 3;
  double Dz = 15;
  double GlobalChi = 10;
  int NMuonhits = 0;
  double Prob = 0.01;


  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetOptTitle(0); // at least most of the time
  gStyle->SetOptStat(1); // most of the time, sometimes "nemriou" might be useful to display name, 
  //number of entries, mean, rms, integral, overflow and underflow
  gStyle->SetOptFit(1); // set to 1 only if you want to display fit results
  
  // Pt bin sizes
  int Nptbin=5;
  double pt_bound[100] = {0};

  //  Pt bins (1S) ; 0 to 2.5, 2.5 to 5, 5 to 8, 8 to 12, 12 to 20. 
  //  Rapidity bins (1S) ; |y| < 0.4 , 0.4 < |y| < 0.7 , 0.7 < |y| < 1.0 , 1.0 < |y| < 1.5, 1.5 < |y| < 2.4


  // centrality for double differential
  double centlow = 0.0/2.5;
  double centup = 100.0/2.5;

 
  if(iSpec == 1) {   //pt
    if(type==1){
      Nptbin = 3;
      pt_bound[0] = 0.0;
      pt_bound[1] = 5.0;
      pt_bound[2] = 12.0;
      pt_bound[3] = 20.0;
    }

    if(type==2){ 
      Nptbin = 5;
      pt_bound[0] = 0.0;
      pt_bound[1] = 2.5;
      pt_bound[2] = 5.0;
      pt_bound[3] = 8.0;
      pt_bound[4] = 12.0;
      pt_bound[5] = 20.0;
    }
  }

  if(iSpec == 2) { // rap
    if(type==1) {
      Nptbin = 2;
      pt_bound[0] = 0.0; 
      pt_bound[1] = 1.2; 
      pt_bound[2] = 2.4; 
    }
    
    if(type==2) {
      Nptbin = 6;
      pt_bound[0] = 0.0; 
      pt_bound[1] = 0.4; 
      pt_bound[2] = 0.8; 
      pt_bound[3] = 1.2; 
      pt_bound[4] = 1.6; 
      pt_bound[5] = 2.0; 
      pt_bound[6] = 2.4; 
    }
  }
  
  if(iSpec == 3) { // cent
    Nptbin = 8;
    pt_bound[0] = 0.0;
    pt_bound[1] = 5.0;
    pt_bound[2] = 10.0;
    pt_bound[3] = 20.0;
    pt_bound[4] = 30.0;
    pt_bound[5] = 40.0;
    pt_bound[6] = 50.0;
    pt_bound[7] = 70.0;
    pt_bound[8] = 100.0;
    for(int i=0; i<Nptbin+1; i++) {  pt_bound[i]=pt_bound[i]/2.5;}
  }



  // Plot weights 
  double pmin = 3.0, pmax = 10.0, pstep = 0.1;
  int np = (pmax-pmin)/pstep;
  double ptp[1000], wep[1000];
  for (int i=0.0; i<np; i++) {
    ptp[i] = pmin + pstep*i;
    wep[i] = getEffEtaPbPb(ptp[i], 1.0);
  }
  new TCanvas;
  TGraph *tnp = new TGraph(np, ptp, wep);
  gPad->SetTicks(1);
  tnp->Draw("APL");



  //===============Input Root File============================================//
  char fileName[10][500];

  int nfile =1;

  bool isWeight=0;

  if(pp==1 && YS==1 && isParent ==1 && genacc==2){sprintf(fileName[0],"pp/Upsi1s_ppMC_All.root");}  
  if(pp==1 && YS==2 && isParent ==1 && genacc==2){sprintf(fileName[0],"pp/Upsi2s_ppMC_All.root");}   
  if(pp==1 && YS==3 && isParent ==1 && genacc==2){sprintf(fileName[0],"pp/Upsi3s_ppMC_All.root");}   

  if(YS==1 && isParent ==1 && genacc==1){sprintf(fileName[0]  ,"ForAcc/DimuonOnia2Dplots_Y1SGen_PP.root");}  
  if((YS==2) && isParent ==1 && genacc==1){sprintf(fileName[0],"ForAcc/DimuonOnia2Dplots_Y2SGen_PP.root");}   
  if((YS==3) && isParent ==1 && genacc==1){sprintf(fileName[0],"ForAcc/DimuonOnia2Dplots_Y3SGen_PP.root");}   

  
  if(pp==2 && YS==1 && isParent ==1 && genacc==2){
    nfile =7;
    sprintf(fileName[0],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt03.root");
    sprintf(fileName[1],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt36.root");
    sprintf(fileName[2],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt69.root");
    sprintf(fileName[3],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt912.root");
    sprintf(fileName[4],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt1215.root");
    sprintf(fileName[5],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt1530.root");
    sprintf(fileName[6],"PbPb/Y1S/Upsi1s_PbPbMC_AllPt30.root");
  isWeight = 1;
  }  

  if(pp==2 && YS==2 && isParent ==1 && genacc==2){
    nfile = 6;
    sprintf(fileName[0],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt03.root");
    sprintf(fileName[1],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt36.root");
    sprintf(fileName[2],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt69.root");
    sprintf(fileName[3],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt912.root");
    sprintf(fileName[4],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt1215.root");
    sprintf(fileName[5],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt1530.root");
    // sprintf(fileName[6],"./PbPb/Y2S/Upsi2s_PbPbMC_AllPt30.root");
    // sprintf(fileName[0], "~/Desktop/AnaUpsilon/MCupsilonPbPb/Upsi2SMC/Upsi2s_DimuTree_Pt03_17Apr14.root");
    // sprintf(fileName[1], "~/Desktop/AnaUpsilon/MCupsilonPbPb/Upsi2SMC/Upsi2s_DimuTree_Pt36_17Apr14.root");
    // sprintf(fileName[2], "~/Desktop/AnaUpsilon/MCupsilonPbPb/Upsi2SMC/Upsi2s_DimuTree_Pt69_17Apr14.root");
    // sprintf(fileName[3], "~/Desktop/AnaUpsilon/MCupsilonPbPb/Upsi2SMC/Upsi2s_DimuTree_Pt912_17Apr14.root");
    // sprintf(fileName[4], "~/Desktop/AnaUpsilon/MCupsilonPbPb/Upsi2SMC/Upsi2s_DimuTree_Pt1215_17Apr14.root");  
    // sprintf(fileName[5], "~/Desktop/AnaUpsilon/MCupsilonPbPb/Upsi2SMC/Upsi2s_DimuTree_Pt1530_17Apr14.root"); 
    isWeight = 1;
  }   
 
  if(pp==2 && YS==4 && isParent ==1 && genacc==2) {
    nfile = 6;
    // filter eff/ #events
    sprintf(fileName[0], "~/Desktop/AnaUpsilon/MCJpsiPbPbTrigger/JpsiMC/Jpsi_DimuTree_Pt03_22Oct13.root");
    sprintf(fileName[1], "~/Desktop/AnaUpsilon/MCJpsiPbPbTrigger/JpsiMC/Jpsi_DimuTree_Pt36_24Oct13.root");   
    sprintf(fileName[2], "~/Desktop/AnaUpsilon/MCJpsiPbPbTrigger/JpsiMC/Jpsi_DimuTree_Pt69_24Oct13.root");
    sprintf(fileName[3], "~/Desktop/AnaUpsilon/MCJpsiPbPbTrigger/JpsiMC/Jpsi_DimuTree_Pt912_24Oct13.root");
    sprintf(fileName[4], "~/Desktop/AnaUpsilon/MCJpsiPbPbTrigger/JpsiMC/Jpsi_DimuTree_Pt1215_24Oct13.root");  
    sprintf(fileName[5], "~/Desktop/AnaUpsilon/MCJpsiPbPbTrigger/JpsiMC/Jpsi_DimuTree_Pt1530_24Oct13.root"); 
    isWeight = 1;
  }

  double masslow = 2.5;
  double massup = 3.5;
  if(YS==1 || YS==2 || YS==3) { masslow = 8.0; massup = 12.0;}

  //  double scale[10]={1.0, 1.0, 410.0/1060.0, 410.0/1060.0*530.0/1350.0, 1.0, 1.0};  
  double scale[10]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};  

  if(pp==2 && YS==1){
    //  {172208,171048,172165,162777,117270,107560,107280};
    scale[0] = 1.0;
    scale[1] = 1.112919104;
    scale[2] = 0.430737304;
    scale[3] = 0.169418472;
    scale[4] = 0.097706729;
    scale[5] = 0.087703928;
    scale[6] = 0.004902708;
    /*
    scale[0] = 514414.0/172214;
    scale[1] = 485009.0/171054;
    scale[2] = 205234.0/172166;
    scale[3] = 86123.0/162783;
    scale[4] = 37875.0/117270;
    scale[5] = 32515.0/107563;
    scale[6] = 1987.0/105280;
    */
  }


  if(pp==2 && YS==2){
    scale[0] =  1.0;
    scale[1] =  1.17201;
    scale[2] =  0.40233;
    scale[3] =  0.13702;
    scale[4] =  0.043731;
    scale[5] =  0.032069;
  }

  if(pp==2 && YS==4){
    scale[0] = 1.0/117512;
    scale[1] = 0.915/124802;
    scale[2] = 0.178/172176;
    scale[3] = 0.0379/167336;
    scale[4] = 0.00979/143944;
    scale[5] = 0.00579/166960;
  }


  //X Axis error on Eff graph 
  double PT[100], DelPT[100], mom_err[100];
  for (Int_t ih = 0; ih < Nptbin; ih++) {
    PT[ih] = (pt_bound[ih] + pt_bound[ih+1])/2.0;
    DelPT[ih] = pt_bound[ih+1] - pt_bound[ih];
    mom_err[ih] = DelPT[ih]/2.0;
  }

  if(iSpec==3) {
    for(int ih=0; ih<Nptbin; ih++) { 
      PT[ih] = PT[ih]*2.5;
      DelPT[ih] = DelPT[ih]*2.5;
      mom_err[ih] = DelPT[ih]/2.0;
    }
  }
  
  double genError, recError;
  double gen_pt[100]={0}, gen_ptError[100]={0}; 
  double gen_ptAll[100]={0}, gen_ptAllError[100]={0}; 
  double rec_pt[100]={0}, rec_ptError[100]={0}; 
  double recTP_pt[100]={0}, recTP_ptError[100]={0}; 

  double recT_pt[100]={0}, recT_ptError[100]={0}; 
  double recI_pt[100]={0}, recI_ptError[100]={0}; 

  double Eff[100]={0}, Err_Eff[100]={0};  
  double EffTP[100]={0}, Err_EffTP[100]={0};  
  double Acc[100]={0}, Err_Acc[100]={0};  
  double Trig[100]={0}, Err_Trig[100]={0};  
  double muEff[100]={0}, Err_muEff[100]={0};  
  double AccEff[100]={0}, Err_AccEff[100]={0};  
  double AccEffTP[100]={0}, Err_AccEffTP[100]={0};  

 //================================= Define Histograms====================================//
  char OutTextFile[100]; 
  char OutTextFileGen[100]; 
  char inTextFileGen[100]; 

  if(pp==1) sprintf(OutTextFile,"EffY_PP.tex");
  if(pp==2) sprintf(OutTextFile,"EffY_PbPb.tex");
  ofstream dataFile(Form(OutTextFile),ios::app);
  
  sprintf(OutTextFileGen,"GENY_PP.tex");
  ofstream outFileGen(Form(OutTextFileGen),ios::app);
  
  sprintf(inTextFileGen,"GENY_PP.tex");
  ifstream inFileGen(Form(inTextFileGen));

  // M, pT, Rap, Centrality (Gen)  
  TH1D *diMuonsInvMass_Gen = new TH1D("diMuonsInvMass_Gen","diMuonsInvMass_Gen", 100,masslow,massup);
  diMuonsInvMass_Gen->Sumw2();
  TH1D *diMuonsInvMass_GenAll = new TH1D("diMuonsInvMass_GenAll","diMuonsInvMass_GenAll", 100,masslow,massup);
  diMuonsInvMass_GenAll->Sumw2();
  TH1D *diMuonsPt_Gen = new TH1D("diMuonsPt_Gen","diMuonsPt_Gen", 400,0,100);
  TH1D *diMuonsRap_Gen = new TH1D("diMuonsRap_Gen","diMuonsRap_Gen", 100,-5,5);
  TH1D *Bin_Gen = new TH1D("Bin_Gen","Bin_Gen", 50,0,50);

  // M, pT, Rap, Centrality (Rec)  
  TH1D *diMuonsInvMass_Rec = new TH1D("diMuonsInvMass_Rec","diMuonsInvMass_Rec", 100,masslow,massup);
  diMuonsInvMass_Rec->Sumw2();
  diMuonsInvMass_Rec->SetMarkerStyle(8);
  diMuonsInvMass_Rec->SetMarkerColor(4);
  diMuonsInvMass_Rec->SetLineColor(4);

  TH1D *diMuonsInvMass_RecTP = new TH1D("diMuonsInvMass_RecTP","diMuonsInvMass_RecTP", 100,masslow,massup);
  diMuonsInvMass_RecTP->Sumw2();
  diMuonsInvMass_RecTP->SetMarkerStyle(8);
  diMuonsInvMass_RecTP->SetMarkerColor(4);
  diMuonsInvMass_RecTP->SetLineColor(4);

  TH1D *diMuonsPt_rec = new TH1D("diMuonsPt_rec","diMuonsPt_rec", 400,0,100);
  TH1D *Bin_Rec = new TH1D("Bin_Rec","Bin_Rec", 50,0,50);  // centrality

  // Single Muon Gen
  TH1D *MuonPt_Gen = new TH1D("MuonPt_Gen","MuonPt_Gen", 100,0,20);
  TH1D *MuonEta_Gen = new TH1D("MuonEta_Gen","MuonEta_Gen", 100,-5,5);
  TH1D *MuonPhi_Gen = new TH1D("MuonPhi_Gen","MuonPhi_Gen", 100,-5,5);

  // Single Muon Rec
  TH1D *MuonPt_Rec = new TH1D("MuonPt_Rec","MuonPt_Rec", 100,0,20);
  TH1D *MuonEta_Rec = new TH1D("MuonEta_Rec","MuonEta_Rec", 100,-5,5);
  TH1D *MuonPhi_Rec = new TH1D("MuonPhi_Rec","MuonPhi_Rec", 100,-5,5);



  // Efficiency 2D maps for Single Muons
  // 1 for |RapY|<0.8 --- 2 for 0.8<|RapY|<1.6 -- 3 for 1.6<|RapY|<2.4
  TH2D *MuonPtEta_Gen = new TH2D("MuonPtEta_Gen","",96,-2.4,2.4,96,0,9.6);
  TH2D *MuonPtEta1_Gen = new TH2D("MuonPtEta1_Gen","",96,-2.4,2.4,96,0,9.6); 
  TH2D *MuonPtEta2_Gen = new TH2D("MuonPtEta2_Gen","",96,-2.4,2.4,96,0,9.6);
  TH2D *MuonPtEta3_Gen = new TH2D("MuonPtEta3_Gen","",96,-2.4,2.4,96,0,9.6);

  TH2D *MuonPtEta_Rec = new TH2D("MuonPtEta_Rec","",96,-2.4,2.4,96,0,9.6);
  TH2D *MuonPtEta1_Rec = new TH2D("MuonPtEta1_Rec","",96,-2.4,2.4,96,0,9.6);
  TH2D *MuonPtEta2_Rec = new TH2D("MuonPtEta2_Rec","",96,-2.4,2.4,96,0,9.6);
  TH2D *MuonPtEta3_Rec = new TH2D("MuonPtEta3_Rec","",96,-2.4,2.4,96,0,9.6);

  TH2D *MuonPtEta1_Eff = new TH2D("MuonPtEta1_Eff","",96,-2.4,2.4,96,0,9.6);


  // Histogram arrays
  TH1D *diMuonsInvMass_GenA[10][1000];
  TH1D *diMuonsInvMass_GenB[10][1000];  // without acceptance

  TH1D *diMuonsInvMass_RecA[10][1000];

  TH1D *diMuonsInvMass_RecATP[10][1000];  // with Tag and Probe

  TH1D *diMuonsInvMass_RecT[10][1000];  // without trigger
  TH1D *diMuonsInvMass_RecI[10][1000];  // without mu ID

  TH1D *diMuonsPt_GenA[10][1000];
  TH1D *diMuonsPt_RecA[10][1000];
  char nameGen[10][500], nameRec[10][500], nameGenPt[10][500], nameRecPt[10][500];
 
  for (int ifile = 0; ifile < nfile; ifile++) {
    for (Int_t ih = 0; ih < Nptbin; ih++) {
      sprintf(nameGen[ifile],"DiMuonMassGen_pt_%d_%d",ih,ifile);
      sprintf(nameRec[ifile],"DiMuonMassRec_pt_%d_%d",ih,ifile);
      
      sprintf(nameGenPt[ifile],"DiMuonPtGen_pt_%d_%d",ih,ifile);
      sprintf(nameRecPt[ifile],"DiMuonPtRec_pt_%d_%d",ih,ifile);
      
      diMuonsInvMass_GenA[ifile][ih]= new TH1D(nameGen[ifile],nameGen[ifile],  100,masslow,massup); //for eff Gen;
      diMuonsInvMass_GenA[ifile][ih]->Sumw2();
      diMuonsInvMass_GenA[ifile][ih]->SetMarkerStyle(7);
      diMuonsInvMass_GenA[ifile][ih]->SetMarkerColor(4);
      diMuonsInvMass_GenA[ifile][ih]->SetLineColor(4);

      sprintf(nameGen[ifile],"DiMuonMassGenB_pt_%d_%d",ih,ifile);
      diMuonsInvMass_GenB[ifile][ih]= new TH1D(nameGen[ifile],nameGen[ifile],  100,masslow,massup); // for acc;
      diMuonsInvMass_GenB[ifile][ih]->Sumw2();
      diMuonsInvMass_GenB[ifile][ih]->SetMarkerStyle(7);
      diMuonsInvMass_GenB[ifile][ih]->SetMarkerColor(4);
      diMuonsInvMass_GenB[ifile][ih]->SetLineColor(4);
      
      diMuonsInvMass_RecA[ifile][ih] = new TH1D(nameRec[ifile],nameRec[ifile], 100,masslow,massup); //for eff Rec;
      diMuonsInvMass_RecA[ifile][ih]->Sumw2();
      diMuonsInvMass_RecA[ifile][ih]->SetMarkerStyle(8);
      diMuonsInvMass_RecA[ifile][ih]->SetMarkerColor(4);
      diMuonsInvMass_RecA[ifile][ih]->SetLineColor(4);

      sprintf(nameRec[ifile],"DiMuonMassRecATP_pt_%d_%d",ih,ifile);
      diMuonsInvMass_RecATP[ifile][ih] = new TH1D(nameRec[ifile],nameRec[ifile], 100,masslow,massup); //for eff Rec;
      diMuonsInvMass_RecATP[ifile][ih]->Sumw2();
      diMuonsInvMass_RecATP[ifile][ih]->SetMarkerStyle(8);
      diMuonsInvMass_RecATP[ifile][ih]->SetMarkerColor(4);
      diMuonsInvMass_RecATP[ifile][ih]->SetLineColor(4);

      sprintf(nameRec[ifile],"DiMuonMassRecT_pt_%d_%d",ih,ifile);
      diMuonsInvMass_RecT[ifile][ih] = new TH1D(nameRec[ifile],nameRec[ifile], 100,masslow,massup); //for eff Rec;
      diMuonsInvMass_RecT[ifile][ih]->Sumw2();
      diMuonsInvMass_RecT[ifile][ih]->SetMarkerStyle(8);
      diMuonsInvMass_RecT[ifile][ih]->SetMarkerColor(4);
      diMuonsInvMass_RecT[ifile][ih]->SetLineColor(4);

      sprintf(nameRec[ifile],"DiMuonMassRecI_pt_%d_%d",ih,ifile);
      diMuonsInvMass_RecI[ifile][ih] = new TH1D(nameRec[ifile],nameRec[ifile], 100,masslow,massup); //for eff Rec;
      diMuonsInvMass_RecI[ifile][ih]->Sumw2();
      diMuonsInvMass_RecI[ifile][ih]->SetMarkerStyle(8);
      diMuonsInvMass_RecI[ifile][ih]->SetMarkerColor(4);
      diMuonsInvMass_RecI[ifile][ih]->SetLineColor(4);
      
      diMuonsPt_GenA[ifile][ih]= new TH1D(nameGenPt[ifile],nameGenPt[ifile],  100,0,100); //
      diMuonsPt_RecA[ifile][ih]= new TH1D(nameRecPt[ifile],nameRecPt[ifile],  100,0,100); //
    }
  }
  
  TFile *infile;
  TTree *tree;
  TTree *gentree;



  //===========File loop ======================
  
  for(int ifile =0; ifile<nfile; ifile++){
    cout<<" file loop "<<endl;
    infile=new TFile(fileName[ifile],"R");

    tree=(TTree*)infile->Get("SingleMuonTree");
    gentree=(TTree*)infile->Get("SingleGenMuonTree");

    //Event variables
    int eventNb, runNb, lumiBlock, gbin, rbin;
    bool evTrigger0;
    //Jpsi Variables
    Double_t JpsiMass, JpsiPt, JpsiRap, JpsiCharge;
    Double_t JpsiVprob;
    //2.) muon variables RECO
    double muPosPx, muPosPy, muPosPz,  muPosEta, muPosPt;
    double muNegPx, muNegPy, muNegPz,  muNegEta, muNegPt;

    //(1).Positive Muon                                     
    double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;
    int muPos_found,muPos_pixeLayers,muPos_nValidMuHits,muPos_arbitrated;
    bool muPos_matches,muPos_tracker,muPos_global;
    int muPos_Trigger1, muPos_Trigger2,muPos_Trigger3,muPos_Trigger4,muPos_Trigger5;
    int muPos_Trigger6,muPos_Trigger7,muPos_Trigger8,muPos_Trigger9,muPos_Trigger10; 

    //(2).Negative Muon                                     
    double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;
    int muNeg_found,muNeg_pixeLayers,muNeg_nValidMuHits,muNeg_arbitrated;
    bool muNeg_matches,muNeg_tracker,muNeg_global;
    int muNeg_Trigger1,muNeg_Trigger2,muNeg_Trigger3,muNeg_Trigger4,muNeg_Trigger5;
    int muNeg_Trigger6,muNeg_Trigger7,muNeg_Trigger8,muNeg_Trigger9,muNeg_Trigger10;

    //Gen Parent Variables
    double GenJpsiMassP, GenJpsiPtP, GenJpsiRapP;
    double GenJpsiPxP, GenJpsiPyP, GenJpsiPzP;
    //Gen JPsi Variables
    double GenJpsiMass, GenJpsiPt, GenJpsiRap;
    double GenJpsiPx, GenJpsiPy, GenJpsiPz;
    // Gen muon variables 
    double GenmuPosPx, GenmuPosPy, GenmuPosPz,  GenmuPosEta, GenmuPosPt, GenmuPosPhi;
    double GenmuNegPx, GenmuNegPy, GenmuNegPz,  GenmuNegEta, GenmuNegPt, GenmuNegPhi;

    //Event variables
    tree->SetBranchAddress("eventNb",&eventNb);
    tree->SetBranchAddress("runNb",&runNb);
    tree->SetBranchAddress("lumiBlock",&lumiBlock);
    tree->SetBranchAddress("evTrigger0",&evTrigger0);

    //Jpsi Variables
    tree->SetBranchAddress("JpsiCharge",&JpsiCharge);
    tree->SetBranchAddress("JpsiMass",&JpsiMass);
    tree->SetBranchAddress("JpsiPt",&JpsiPt);
    tree->SetBranchAddress("JpsiRap",&JpsiRap);
    tree->SetBranchAddress("JpsiVprob",&JpsiVprob);
    tree->SetBranchAddress("rbin",&rbin);
    //muon variables
    tree->SetBranchAddress("muPosPx",&muPosPx);
    tree->SetBranchAddress("muPosPy",&muPosPy);
    tree->SetBranchAddress("muPosPz",&muPosPz);
    tree->SetBranchAddress("muPosEta",&muPosEta);
    tree->SetBranchAddress("muNegPx", &muNegPx);
    tree->SetBranchAddress("muNegPy",    &muNegPy);
    tree->SetBranchAddress("muNegPz",    &muNegPz);
    tree->SetBranchAddress("muNegEta",    &muNegEta);
    //1). Positive Muon
    tree->SetBranchAddress("muPos_nchi2In", &muPos_nchi2In);
    tree->SetBranchAddress("muPos_dxy", &muPos_dxy);
    tree->SetBranchAddress("muPos_dz", &muPos_dz);
    tree->SetBranchAddress("muPos_nchi2Gl", &muPos_nchi2Gl);
    tree->SetBranchAddress("muPos_found", &muPos_found);
    tree->SetBranchAddress("muPos_pixeLayers", &muPos_pixeLayers);
    tree->SetBranchAddress("muPos_nValidMuHits", &muPos_nValidMuHits);
    tree->SetBranchAddress("muPos_matches", &muPos_matches);
    tree->SetBranchAddress("muPos_tracker", &muPos_tracker);
    tree->SetBranchAddress("muPos_global",&muPos_global);
    tree->SetBranchAddress("muPos_arbitrated", &muPos_arbitrated);
    
    tree->SetBranchAddress("muPos_Trigger1",&muPos_Trigger1);
    tree->SetBranchAddress("muPos_Trigger2",&muPos_Trigger2);
    tree->SetBranchAddress("muPos_Trigger3",&muPos_Trigger3);
    tree->SetBranchAddress("muPos_Trigger4",&muPos_Trigger4);
    tree->SetBranchAddress("muPos_Trigger5",&muPos_Trigger5);
    tree->SetBranchAddress("muPos_Trigger6",&muPos_Trigger6);
    tree->SetBranchAddress("muPos_Trigger7",&muPos_Trigger7);
    tree->SetBranchAddress("muPos_Trigger8",&muPos_Trigger8);
    tree->SetBranchAddress("muPos_Trigger9",&muPos_Trigger9);
    tree->SetBranchAddress("muPos_Trigger10",&muPos_Trigger10);

    //2). Negative Muon                                                                            
    tree->SetBranchAddress("muNeg_nchi2In", &muNeg_nchi2In);
    tree->SetBranchAddress("muNeg_dxy", &muNeg_dxy);
    tree->SetBranchAddress("muNeg_dz", &muNeg_dz);
    tree->SetBranchAddress("muNeg_nchi2Gl", &muNeg_nchi2Gl);
    tree->SetBranchAddress("muNeg_found", &muNeg_found);
    tree->SetBranchAddress("muNeg_pixeLayers", &muNeg_pixeLayers);
    tree->SetBranchAddress("muNeg_nValidMuHits", &muNeg_nValidMuHits);
    tree->SetBranchAddress("muNeg_matches", &muNeg_matches);
    tree->SetBranchAddress("muNeg_tracker", &muNeg_tracker);
    tree->SetBranchAddress("muNeg_global",&muNeg_global);
    tree->SetBranchAddress("muNeg_arbitrated", &muNeg_arbitrated);

    tree->SetBranchAddress("muNeg_Trigger1",&muNeg_Trigger1);
    tree->SetBranchAddress("muNeg_Trigger2",&muNeg_Trigger2);
    tree->SetBranchAddress("muNeg_Trigger3",&muNeg_Trigger3);
    tree->SetBranchAddress("muNeg_Trigger4",&muNeg_Trigger4);
    tree->SetBranchAddress("muNeg_Trigger5",&muNeg_Trigger5);
    tree->SetBranchAddress("muNeg_Trigger6",&muNeg_Trigger6);
    tree->SetBranchAddress("muNeg_Trigger7",&muNeg_Trigger7);
    tree->SetBranchAddress("muNeg_Trigger8",&muNeg_Trigger8);
    tree->SetBranchAddress("muNeg_Trigger9",&muNeg_Trigger9);
    tree->SetBranchAddress("muNeg_Trigger10",&muNeg_Trigger10);
    //HLT_HIL1DoubleMuOpen_v1  HLT_HIL1DoubleMu0_HighQ_v1  HLT_HIL2Mu3_v1  HLT_HIL2Mu3_NHitQ_v1  HLT_HIL2Mu7_v1
    //HLT_HIL2Mu15_v1  HLT_HIL2DoubleMu0_v1  HLT_HIL2DoubleMu0_NHitQ_v1  HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1  HLT_HIL2DoubleMu3_v1
    
    // Gen Parent Variables
    gentree->SetBranchAddress("GenJpsiMassP",   &GenJpsiMassP);
    gentree->SetBranchAddress("GenJpsiPtP",     &GenJpsiPtP);
    gentree->SetBranchAddress("GenJpsiRapP",    &GenJpsiRapP);
    gentree->SetBranchAddress("GenJpsiPxP",     &GenJpsiPxP);
    gentree->SetBranchAddress("GenJpsiPyP",     &GenJpsiPyP);
    gentree->SetBranchAddress("GenJpsiPzP",     &GenJpsiPzP);

    //Gen Jpsi Variables
    gentree->SetBranchAddress("GenJpsiMass",   &GenJpsiMass);
    gentree->SetBranchAddress("GenJpsiPt",     &GenJpsiPt);
    gentree->SetBranchAddress("GenJpsiRap",    &GenJpsiRap);
    gentree->SetBranchAddress("GenJpsiPx",     &GenJpsiPx);
    gentree->SetBranchAddress("GenJpsiPy",     &GenJpsiPy);
    gentree->SetBranchAddress("GenJpsiPz",     &GenJpsiPz);
    gentree->SetBranchAddress("gbin",&gbin);

    //muon variables
    gentree->SetBranchAddress("GenmuPosPx",    &GenmuPosPx);
    gentree->SetBranchAddress("GenmuPosPy",    &GenmuPosPy);
    gentree->SetBranchAddress("GenmuPosPz",    &GenmuPosPz);
    gentree->SetBranchAddress("GenmuPosEta",    &GenmuPosEta);
    gentree->SetBranchAddress("GenmuNegPx",    &GenmuNegPx);
    gentree->SetBranchAddress("GenmuNegPy",    &GenmuNegPy);
    gentree->SetBranchAddress("GenmuNegPz",    &GenmuNegPz);
    gentree->SetBranchAddress("GenmuNegEta",    &GenmuNegEta);
    

    //============================== Gen tree loop ================================
    int NAccep=0;
    int nGenEntries=gentree->GetEntries();
    cout<<" Entries in Gen Tree for pT range: "<<fileName[ifile]<<"  "<<   nGenEntries<< " ========="<<endl;
    
    for(int i=0; i< nGenEntries; i++)  {	    
      gentree->GetEntry(i);
      
      if(isParent==1){
	GenJpsiMass=GenJpsiMassP; 
	GenJpsiPt=GenJpsiPtP;
	GenJpsiRap=GenJpsiRapP;
	GenJpsiPx=GenJpsiPxP;
	GenJpsiPy=GenJpsiPyP;
	GenJpsiPz=GenJpsiPzP;
      }

      if(i%200000==0){
	cout<<" processing record Gen Tree "<<i<<endl;
	cout<<" Mass "<< GenJpsiMass<< " pT "<< GenJpsiPt << " Y " <<GenJpsiRap<<endl; 
      }
      
      GenmuPosPt= TMath::Sqrt(GenmuPosPx*GenmuPosPx + GenmuPosPy*GenmuPosPy); 
      GenmuNegPt= TMath::Sqrt(GenmuNegPx*GenmuNegPx + GenmuNegPy*GenmuNegPy); 
      
      GenmuPosPhi=TMath::ACos(GenmuPosPx/GenmuPosPt); 
      GenmuNegPhi=TMath::ACos(GenmuNegPx/GenmuNegPt); 
 
      // Acceptance and pT cut
      bool GenPosIn=0,GenNegIn=0,GenPosPass=0,GenNegPass=0;
      if(IsAccept(GenmuPosPt, GenmuPosEta)) {GenPosIn=1;}
      if(IsAccept(GenmuNegPt, GenmuNegEta)) {GenNegIn=1;}

      if(GenmuPosPt> PtCut) {GenPosPass=1;}
      if(GenmuNegPt> PtCut) {GenNegPass=1;}
       
      bool genGen=0;
      genGen = GenJpsiPt > 0.0 && TMath::Abs(GenJpsiRap)<2.4;

      bool acceptedGen=0;
      bool passpt=0;
	
     if(YS==1) passpt = ( (GenmuPosPt> 4.0 && GenmuNegPt> 3.5) || (GenmuPosPt> 3.5 && GenmuNegPt> 4.0) );

      if(YS==2) passpt = ( (GenmuPosPt> 4.0 && GenmuNegPt> 3.5) || (GenmuPosPt> 3.5 && GenmuNegPt> 4.0) );
      //      if(YS==2) passpt = (GenmuPosPt> 4.0 && GenmuNegPt> 4.0);

      if(YS==3) passpt = (GenmuPosPt> 4.0 && GenmuNegPt> 4.0);


      acceptedGen = (GenPosIn==1 && GenNegIn==1) && passpt && TMath::Abs(GenJpsiRap)<2.4;

      // Weights 
      double GenCenWeight=0,GenWeight=0;
      GenCenWeight=FindCenWeight(gbin);
      GenWeight=GenCenWeight*scale[ifile];
      if(isWeight==0) GenWeight=1;

      if(genGen) diMuonsInvMass_GenAll->Fill(GenJpsiMass, GenWeight);

      if(genGen && acceptedGen) {
	diMuonsInvMass_Gen->Fill(GenJpsiMass, GenWeight);
	NAccep++;
      }	

      Bin_Gen->Fill(gbin);      
      diMuonsPt_Gen->Fill(GenJpsiPt, GenWeight);  
      diMuonsRap_Gen->Fill(GenJpsiRap, GenWeight);

      MuonPt_Gen->Fill(GenmuPosPt, GenWeight);
      MuonPt_Gen->Fill(GenmuNegPt, GenWeight);
      
      MuonEta_Gen->Fill(GenmuPosEta, GenWeight);
      MuonEta_Gen->Fill(GenmuNegEta, GenWeight);
	
      MuonPhi_Gen->Fill(GenmuPosPhi, GenWeight);
      MuonPhi_Gen->Fill(GenmuNegPhi, GenWeight);

      if(genGen){
      MuonPtEta_Gen ->Fill(GenmuPosEta,GenmuPosPt, GenWeight);
      MuonPtEta_Gen ->Fill(GenmuNegEta,GenmuNegPt,GenWeight); 
    }

      if(genGen &&acceptedGen && TMath::Abs(GenJpsiRap)<0.8){
      MuonPtEta1_Gen ->Fill(GenmuPosEta,GenmuPosPt, GenWeight);
      MuonPtEta1_Gen ->Fill(GenmuNegEta,GenmuNegPt,GenWeight); 
      }
      if(genGen &&acceptedGen&& TMath::Abs(GenJpsiRap)>0.8 && TMath::Abs(GenJpsiRap)<1.6){
      MuonPtEta2_Gen ->Fill(GenmuPosEta,GenmuPosPt, GenWeight); 
      MuonPtEta2_Gen ->Fill(GenmuNegEta,GenmuNegPt, GenWeight); 

      }
      if(genGen && acceptedGen&&TMath::Abs(GenJpsiRap)>1.6 && TMath::Abs(GenJpsiRap)<2.4){
      MuonPtEta3_Gen ->Fill(GenmuPosEta,GenmuPosPt, GenWeight);
      MuonPtEta3_Gen ->Fill(GenmuNegEta,GenmuNegPt, GenWeight);
 
      }

      bool CentCut = (gbin>=centlow && gbin<centup);
      if(pp==1 || genacc==1) CentCut=1;

      for (Int_t ih = 0; ih < Nptbin; ih++) {
	//adding pT of all pt bins to see diss is cont
	if(iSpec == 1)  if(genGen && CentCut && (GenJpsiPt>pt_bound[ih] && GenJpsiPt<=pt_bound[ih+1])) {
	    if(acceptedGen) diMuonsPt_GenA[ifile][ih]->Fill(GenJpsiPt, GenWeight);
	    if(acceptedGen) diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass, GenWeight);
	    diMuonsInvMass_GenB[ifile][ih]->Fill(GenJpsiMass, GenWeight);
	  }
	//for non symetric plots
	if(iSpec == 2)  if(genGen && CentCut && (TMath::Abs(GenJpsiRap) > pt_bound[ih] && TMath::Abs(GenJpsiRap) <=pt_bound[ih+1] )) {
	    if(acceptedGen) {diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass, GenWeight);}
	    diMuonsInvMass_GenB[ifile][ih]->Fill(GenJpsiMass, GenWeight );
	  }
	if(iSpec == 3)  if(genGen && (gbin>=pt_bound[ih] && gbin<pt_bound[ih+1]) ) {
	    if(acceptedGen) {diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass,GenWeight );}
	    diMuonsInvMass_GenB[ifile][ih]->Fill(GenJpsiMass, GenWeight);
	  }
      }
    }//gen loop end
    
    cout<<" accepted no "<< NAccep<<endl;
    cout<<" End of Gen Tree " <<endl<< endl;

    //=============== Rec Tree Loop ==============================================//
    
    int nRecEntries=tree->GetEntries();
    cout<<"Entries in reco Tree for pT range "<<fileName[ifile]<<"  "<<nRecEntries<< "====="<<endl;
    
    for(int i=0; i<nRecEntries; i++)  {	    
      tree->GetEntry(i);
      
      if(i%200000==0){
	cout<<" processing record Rec Tree"<<i<<endl;
	cout<<" processing Run  " <<runNb <<" event "<<eventNb<<" lum block "<<lumiBlock<<endl;    
	cout<<" Mass "<< JpsiMass<< " pT "<< JpsiPt << " Y " <<JpsiRap<<"  "<<JpsiVprob<<" charge "<<JpsiCharge<<endl; 
      }

      muPosPt= TMath::Sqrt(muPosPx*muPosPx + muPosPy*muPosPy); 
      //      muPosP = TMath::Sqrt(muPosPx*muPosPx + muPosPy*muPosPy+ muPosPz*muPosPz); 
      muNegPt= TMath::Sqrt(muNegPx*muNegPx + muNegPy*muNegPy); 
      //      muNegP = TMath::Sqrt(muNegPx*muNegPx + muNegPy*muNegPy +muNegPz*muNegPz); 

      // acceptance cuts
      bool PosIn=0, NegIn=0, accepted=0;    
      if(IsAccept(muPosPt, muPosEta)){PosIn=1;}
      if(IsAccept(muNegPt, muNegEta)){NegIn=1;}
      
     
      bool passpt=0;
     
     if(YS==1) passpt = ( (muPosPt> 4.0 && muNegPt> 3.5) || (muPosPt> 3.5 && muNegPt> 4.0) );

      if(YS==2) passpt = ( (muPosPt> 4.0 && muNegPt> 3.5) || (muPosPt> 3.5 && muNegPt> 4.0) );

      //      if(YS==2) passpt = (muPosPt> 4.0 && muNegPt> 4.0);

      if(YS==3) passpt = (muPosPt> 4.0 && muNegPt> 4.0);

      accepted = (PosIn==1 && NegIn==1) && passpt && (TMath::Abs(JpsiRap) < 2.4);

      // Quality Cuts
      bool PosPass=0, NegPass=0, muQpassed=0;

      // Quality cuts
      //      track.numberOfValidHits > 10 && track.normalizedChi2 < 4 && track.hitPattern.pixelLayersWithMeasurement > 0 
      //      && abs(track.dxy) < 3 && abs(track.dz) < 15
      //      && isGlobalMuon && isTrackerMuon
      //      globalTrack.normalizedChi2 < 20 && muonID('TrackerMuonArbitrated')"

      if(muPos_found > Nhits  && TMath::Abs(muPos_nchi2In) < InnerChi && muPos_pixeLayers > NPxlLayers    
         && TMath::Abs(muPos_dxy) < Dxy && TMath::Abs(muPos_dz) < Dz
	 && muPos_global==1 && muPos_tracker==1
	 && TMath::Abs(muPos_nchi2Gl) < GlobalChi && muPos_arbitrated==1 ){
	PosPass=1;
      }
      if(muNeg_found > Nhits  && TMath::Abs(muNeg_nchi2In) < InnerChi && muNeg_pixeLayers > NPxlLayers 
         && TMath::Abs(muNeg_dxy) < Dxy && TMath::Abs(muNeg_dz) < Dz 
	 && muNeg_global==1 && muNeg_tracker==1 
	 && TMath::Abs(muNeg_nchi2Gl) < GlobalChi && muNeg_arbitrated==1 ) {
	NegPass=1;
      }
      muQpassed = (PosPass==1 && NegPass==1);

      // Triggered Matched
      bool TrigMatch = 0;

      if(pp==1 && (YS==1 || YS==2 || YS==3)) TrigMatch = (muPos_Trigger2==1 && muNeg_Trigger2==1); // (Upsilon pp) HLT_HIL1DoubleMu0_HighQ_v1
      if(pp==2 && (YS==1 || YS==2 || YS==3)) TrigMatch = (muPos_Trigger1==1 && muNeg_Trigger1==1); // (Upsilon PbPb) HLT_HIL1DoubleMu0_HighQ_v1
      if(YS==4) TrigMatch = (muPos_Trigger10==1 && muNeg_Trigger10==1);    // J/psi

      bool isDimuon =0;
      isDimuon = evTrigger0==1 && (JpsiCharge == 0) && (JpsiVprob > Prob) && (JpsiPt> 0.0  && JpsiPt<1000.0);

      // all cut
      bool muCuts=0;
      if(accepted && muQpassed && TrigMatch){muCuts=1;}

      double RecCenWeight=0,RecWeight=0;
      RecCenWeight=FindCenWeight(rbin);
      RecWeight=RecCenWeight*scale[ifile];
      if(isWeight==0)  RecWeight=1;
      
      Bin_Rec->Fill(rbin);

      double tnpWeight_1=1.0, tnpWeight_2=1.0;

      // Efficiency correction by Tag And Probe
      if(pp==1) {
	tnpWeight_1 = getEffEtaPP(muPosPt, muPosEta);
	tnpWeight_2 = getEffEtaPP(muNegPt, muNegEta);
      }
      if(pp==2) {
	tnpWeight_1 = getEffEtaPbPb(muPosPt, muPosEta);
	tnpWeight_2 = getEffEtaPbPb(muNegPt, muNegEta);
      }
      double tnpWeight = tnpWeight_1*tnpWeight_2;


      if(isDimuon && muCuts && (JpsiMass > masslow && JpsiMass < massup)) {
	diMuonsInvMass_Rec->Fill(JpsiMass, RecWeight );
	diMuonsInvMass_RecTP->Fill(JpsiMass, RecWeight*tnpWeight);

	diMuonsPt_rec->Fill(JpsiPt, RecWeight);
      }
      
      if(TMath::Abs(JpsiRap) < 2.4&&JpsiPt>0. &&muQpassed){
      MuonPtEta_Rec ->Fill(muPosEta,muPosPt, RecWeight);
      MuonPtEta_Rec ->Fill(muNegEta,muNegPt, RecWeight); 
   }

      if(isDimuon && muCuts && TMath::Abs(JpsiRap)<.8 && (JpsiMass > masslow && JpsiMass < massup)){
      MuonPtEta1_Rec ->Fill(muPosEta,muPosPt, RecWeight);
      MuonPtEta1_Rec ->Fill(muNegEta,muNegPt, RecWeight); 
      }
      if(isDimuon && muCuts && TMath::Abs(JpsiRap)>.8 && TMath::Abs(JpsiRap)<1.6 && (JpsiMass > masslow && JpsiMass < massup)){
      MuonPtEta2_Rec ->Fill(muPosEta,muPosPt, RecWeight); 
      MuonPtEta2_Rec ->Fill(muNegEta,muNegPt, RecWeight); 
      }
      if(isDimuon && muCuts &&  TMath::Abs(JpsiRap)>1.6 && TMath::Abs(JpsiRap)<2.4 && (JpsiMass > masslow && JpsiMass < massup)){
      MuonPtEta3_Rec ->Fill(muPosEta,muPosPt, RecWeight); 
      MuonPtEta3_Rec ->Fill(muNegEta,muNegPt, RecWeight);
      }


      bool CentCut = (rbin>=centlow && rbin<centup);
      if(pp==1 || genacc==1) CentCut=1;

      //Rec pt Loop for reco
      if(isDimuon) {
	for (Int_t ih = 0; ih < Nptbin; ih++) {
	  if(iSpec == 1) if( accepted &&  CentCut && (JpsiPt>pt_bound[ih] && JpsiPt<=pt_bound[ih+1])) {
	      if(muQpassed && TrigMatch) diMuonsInvMass_RecATP[ifile][ih]->Fill(JpsiMass, RecWeight*tnpWeight);
	      if(muQpassed && TrigMatch) diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass, RecWeight );
	      if(muQpassed && TrigMatch) diMuonsPt_RecA[ifile][ih]->Fill(JpsiPt,  RecWeight);
	      if(muQpassed)  diMuonsInvMass_RecT[ifile][ih]->Fill(JpsiMass,  RecWeight);
	      if(TrigMatch)  diMuonsInvMass_RecI[ifile][ih]->Fill(JpsiMass,  RecWeight);
	    }
	  //for non symetric plots
	  if(iSpec == 2) if( accepted  && CentCut && (TMath::Abs(JpsiRap) > pt_bound[ih] && TMath::Abs(JpsiRap)<=pt_bound[ih+1]) ){
	      if(muQpassed && TrigMatch) diMuonsInvMass_RecATP[ifile][ih]->Fill(JpsiMass, RecWeight*tnpWeight);
	      if(muQpassed && TrigMatch) diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass,  RecWeight);
	      if(muQpassed)  diMuonsInvMass_RecT[ifile][ih]->Fill(JpsiMass,  RecWeight);
	      if(TrigMatch)  diMuonsInvMass_RecI[ifile][ih]->Fill(JpsiMass,  RecWeight);
	    }
	  if(iSpec == 3) if( accepted  && (rbin>=pt_bound[ih] &&  rbin < pt_bound[ih+1])) {
	      if(muQpassed && TrigMatch) diMuonsInvMass_RecATP[ifile][ih]->Fill(JpsiMass, RecWeight*tnpWeight);
	      if(muQpassed && TrigMatch) diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass,  RecWeight);
	      if(muQpassed)  diMuonsInvMass_RecT[ifile][ih]->Fill(JpsiMass,  RecWeight);
	      if(TrigMatch)  diMuonsInvMass_RecI[ifile][ih]->Fill(JpsiMass,  RecWeight);
	    }

	} // for over pT 
      } // if prob
    } //rec entry loop    
  }  // file loop ends 
  
  ///////////////////////////////////////////////////////////////////

    
  new TCanvas;
  diMuonsInvMass_Gen->Draw();
  gPad->Print("plots/diMuonsInvMass_GenMC.pdf");
  gPad->Print("plots/diMuonsInvMass_GenMC.png");
  new TCanvas;
  gPad->SetLogy(1);
  diMuonsPt_Gen->Draw();
  gPad->Print("plots/diMuonsPt_GenMC.pdf");
  gPad->Print("plots/diMuonsPt_GenMC.png");
  new TCanvas;
  gPad->SetLogy(0);
  diMuonsRap_Gen->Draw();
  gPad->Print("plots/diMuonsRap_GenMC.pdf");
  gPad->Print("plots/diMuonsRap_GenMC.png");
  

  cout<< " adding "<<endl;
  TH1D *diMuonsInvMass_RecA1[100];
  TH1D *diMuonsInvMass_RecATP1[100];
  TH1D *diMuonsInvMass_RecT1[100];
  TH1D *diMuonsInvMass_RecI1[100];
  TH1D *diMuonsInvMass_GenA1[100];
  TH1D *diMuonsInvMass_GenB1[100];
  TH1D *diMuonsPt_RecA1[100];
  TH1D *diMuonsPt_GenA1[100];
  for(Int_t ih = 0; ih < Nptbin; ih++){
    diMuonsInvMass_RecA1[ih] = diMuonsInvMass_RecA[0][ih];
    diMuonsInvMass_RecATP1[ih] = diMuonsInvMass_RecATP[0][ih];
    diMuonsInvMass_RecT1[ih] = diMuonsInvMass_RecT[0][ih];
    diMuonsInvMass_RecI1[ih] = diMuonsInvMass_RecI[0][ih];
    diMuonsInvMass_GenA1[ih] = diMuonsInvMass_GenA[0][ih];
    diMuonsInvMass_GenB1[ih] = diMuonsInvMass_GenB[0][ih];
    diMuonsPt_RecA1[ih] = diMuonsPt_RecA[0][ih];
    diMuonsPt_GenA1[ih] = diMuonsPt_GenA[0][ih];
  }

  for(int ifile=1; ifile < nfile; ifile++) {  
    for(Int_t ih = 0; ih < Nptbin; ih++){
      diMuonsInvMass_RecA1[ih]->Add(diMuonsInvMass_RecA[ifile][ih]);
      diMuonsInvMass_RecATP1[ih]->Add(diMuonsInvMass_RecATP[ifile][ih]);
      diMuonsInvMass_RecT1[ih]->Add(diMuonsInvMass_RecT[ifile][ih]);
      diMuonsInvMass_RecI1[ih]->Add(diMuonsInvMass_RecI[ifile][ih]);
      diMuonsInvMass_GenA1[ih]->Add(diMuonsInvMass_GenA[ifile][ih]);
      diMuonsInvMass_GenB1[ih]->Add(diMuonsInvMass_GenB[ifile][ih]);
      diMuonsPt_RecA1[ih]->Add(diMuonsPt_RecA[ifile][ih]);
      diMuonsPt_GenA1[ih]->Add(diMuonsPt_GenA[ifile][ih]);
    }
  }
    
  //===========================Fitting======================================================//
  // Fit ranges
  double massFitlow, massFithigh;
  double MassUpsilon, WidthUpsilon;
  if(YS==1){MassUpsilon = 9.46; WidthUpsilon = 0.09; massFitlow = 8.0; massFithigh = 10.5;}
  if(YS==2){MassUpsilon = 10.05; WidthUpsilon = 0.09; massFitlow = 9.0; massFithigh = 11.0;}
  if(YS==3){MassUpsilon = 10.35; WidthUpsilon = 0.09; massFitlow = 9.35; massFithigh = 11.35;}

  if(YS==4){MassUpsilon = 3.1; WidthUpsilon = 0.04; massFitlow = 2.7; massFithigh = 3.4;}
  
  double MassLow,MassHigh;
  if(YS==1){MassLow=9.0; MassHigh=10.0;}
  if(YS==2){MassLow=9.5; MassHigh=10.5;}
  if(YS==3){MassLow=9.8; MassHigh=10.8;}
  if(YS==4){MassLow=2.9; MassHigh=3.3;}

  // Fit ranges
  TF1 *backfun_1;
  char namePt_1B[500]; //for bkg func

  // Fit Function crystall ball
  TF1 *GAUSPOL = new TF1("GAUSPOL",CBPol1,masslow,massup,8);
  GAUSPOL->SetParNames("Yield (#varUpsilon)","BinWidth","Mean","Sigma","#alpha","n");
 
  GAUSPOL->SetParameter(2, MassUpsilon);
  GAUSPOL->SetParLimits(2, massFitlow, massFithigh);  

  GAUSPOL->SetParameter(3, WidthUpsilon);
  GAUSPOL->SetParLimits(3, 0.1*WidthUpsilon,5.0*WidthUpsilon);
  
  GAUSPOL->FixParameter(4, 1.6); //1.0
  GAUSPOL->SetParameter(5, 2.5); //2.0
 
  if(YS==1) GAUSPOL->SetLineColor(2);
  if(YS==2) GAUSPOL->SetLineColor(1);
  if(YS==3) GAUSPOL->SetLineColor(1);
 
  if(YS==4) GAUSPOL->SetLineColor(1);

  new TCanvas;
  //  diMuonsInvMass_Rec;
  //  diMuonsInvMass_Rec->Rebin(2);
  //  GAUSPOL->SetParameter(0, diMuonsInvMass_Rec->Integral(0,50));
  //  GAUSPOL->FixParameter(1, diMuonsInvMass_Rec->GetBinWidth(1));

  //  diMuonsInvMass_RecA1[0]->Rebin(2);
  GAUSPOL->SetParameter(0, diMuonsInvMass_RecA1[0]->Integral(0,50));
  GAUSPOL->FixParameter(1, diMuonsInvMass_RecA1[0]->GetBinWidth(1));

  diMuonsInvMass_RecA1[0]->Fit("GAUSPOL","LLMERQ", "", massFitlow, massFithigh);
  diMuonsInvMass_RecA1[0]->Draw();
  gPad->SaveAs("plots/UpsilonInMass_PPMC.pdf");
  gPad->SaveAs("plots/UpsilonInMass_PPMC.png");

  new TCanvas;
  gPad->SetLogy(1);
  diMuonsPt_rec->Draw();
  gPad->SaveAs("plots/UpsilonPt_PPMC.pdf");
  gPad->SaveAs("plots/UpsilonPt_PPMC.png");

  //Cal eff 
  double gen_Total = diMuonsInvMass_Gen->IntegralAndError(1,100,genError);
  double gen_Error_Total= genError;

  double gen_TotalAll = diMuonsInvMass_GenAll->IntegralAndError(1,100,genError);
  double gen_Error_TotalAll= genError;


  //yield by histogram integral

  int binlow1 =diMuonsInvMass_Rec->GetXaxis()->FindBin(MassLow);
  int binhi1 =diMuonsInvMass_Rec->GetXaxis()->FindBin(MassHigh);

  double rec_Total = diMuonsInvMass_Rec->IntegralAndError(binlow1, binhi1,recError);
  double rec_Error_Total = recError;

  double recTP_Total = diMuonsInvMass_RecTP->IntegralAndError(binlow1, binhi1,recError);
  double recTP_Error_Total = recError;

  double Eff_Total = rec_Total/gen_Total; 
  double Err_Eff_Total = RError(rec_Total,rec_Error_Total,gen_Total,gen_Error_Total);  

  double EffTP_Total = recTP_Total/rec_Total; 
  double Err_EffTP_Total = RError(recTP_Total,recTP_Error_Total,rec_Total,rec_Error_Total);  

  double Acc_Total = gen_Total/gen_TotalAll; 
  double Err_Acc_Total = RError(gen_Total,gen_Error_Total,gen_TotalAll,gen_Error_TotalAll);  
  
  if(genacc==1) {
    outFileGen<< Acc_Total<<"  "<<Err_Acc_Total<< endl;
  }
  if(genacc==2) {
    inFileGen >> Acc_Total >>  Err_Acc_Total;
  }
  
  double AccEff_Total = Acc_Total*Eff_Total; 
  double Err_AccEff_Total = PError(Acc_Total,Err_Acc_Total,Eff_Total,Err_Eff_Total);  

  double AccEffTP_Total = AccEff_Total*EffTP_Total; 
  double Err_AccEffTP_Total = PError(AccEff_Total,Err_AccEff_Total,EffTP_Total,Err_EffTP_Total);  


  cout<< " ////////////////////////////////////////////////////////" << endl;
  cout<<"Upsilon eff Total "<< Eff_Total<<" error "<<Err_Eff_Total<< endl;
  cout<<"Upsilon Acc Total "<< Acc_Total<<" error "<<Err_Acc_Total<< endl;
  cout<<"Upsilon Acc*eff Total "<< AccEff_Total<<" error "<<Err_AccEff_Total<< endl;

  //  dataFile<<"Upsilon eff Total = " << Eff_Total<<" $\\pm$ "<<Err_Eff_Total<< endl;
  //  dataFile<<"Upsilon Acc Total = "<< Acc_Total<<" $\\pm$ "<<Err_Acc_Total<< endl;
  //  dataFile<<"Upsilon Acc*eff Total = "<< AccEff_Total<<" $\\pm$ "<<Err_AccEff_Total<< endl;

  dataFile<< "  " <<setprecision(3)<<fixed<< "0.0-100.0 " <<"  &  "<<Eff_Total<<"  $\\pm$  "<<Err_Eff_Total   <<"  &  " <<Acc_Total<<"  $\\pm$  "<<Err_Acc_Total <<"  &  "<<AccEff_Total<<"  $\\pm$  "<<Err_AccEff_Total  <<"  &  "<< EffTP_Total<<"  $\\pm$  "<<Err_EffTP_Total <<"  &  "<< AccEffTP_Total<<"  $\\pm$  "<<Err_AccEffTP_Total << " \\\ "    << endl;


  cout<< " ////////////////////////////////////////////////////////" << endl;


  //=====================Loop for eff==================================================

  TCanvas *CanCentR = new TCanvas("CanCentR","CanCentR",20,20,1100,700);
  CanCentR->Divide(3,2);

  char Spectra[100];
  if(iSpec==1) sprintf(Spectra,"Pt");
  if(iSpec==2) sprintf(Spectra,"Rap");
  if(iSpec==3) sprintf(Spectra,"Cent");

  if(YS==1) dataFile<<"Y(1S) Efficiency "<<Spectra<<" ========= "<<endl;
  if(YS==2) dataFile<<"Y(2S) Efficiency "<<Spectra<<" ========= "<<endl;
  if(YS==3) dataFile<<"Y(3S) Efficiency "<<Spectra<<" ========= "<<endl;
  if(YS==4) dataFile<<"J/psi Efficiency "<<Spectra<<" ========= "<<endl;

  for (Int_t ih = 0; ih < Nptbin; ih++) {
    CanCentR->cd(ih+1);
   
    gen_pt[ih] = diMuonsInvMass_GenA1[ih]->IntegralAndError(1,100,genError);
    gen_ptError[ih]= genError;
    gen_ptAll[ih] = diMuonsInvMass_GenB1[ih]->IntegralAndError(1,100,genError);
    gen_ptAllError[ih]= genError;

    diMuonsInvMass_RecA1[ih]->Rebin(2);
    diMuonsInvMass_RecATP1[ih]->Rebin(2);
    diMuonsInvMass_RecT1[ih]->Rebin(2);
    diMuonsInvMass_RecI1[ih]->Rebin(2);

    GAUSPOL->SetParameter(0, diMuonsInvMass_RecA1[ih]->Integral(0,50));
    GAUSPOL->FixParameter(1, diMuonsInvMass_RecA1[ih]->GetBinWidth(1));
    
    //new TCanvas; 
    gPad->SetLogy(0);
    diMuonsInvMass_RecA1[ih]->Fit("GAUSPOL","LLMERQ", "", massFitlow, massFithigh);
    diMuonsInvMass_RecA1[ih]->Draw();

    double UpsilonMass = GAUSPOL->GetParameter(2);
    double UpsilonWidth = GAUSPOL->GetParameter(3);
    double UpsilonYield = GAUSPOL->GetParameter(0);
    //    double UpsilonYieldError = GAUSPOL->GetParError(0);
    double par[20];
    GAUSPOL->GetParameters(par);
    sprintf(namePt_1B,"pt_1B_%d",ih);
    backfun_1 = new TF1(namePt_1B, Pol1, massFitlow, massFithigh, 2);
    backfun_1->SetParameters(&par[6]);
    backfun_1->SetLineColor(4);
    backfun_1->SetLineWidth(1);
    backfun_1->Draw("same");

    //yield by histogram integral
    binlow1 =diMuonsInvMass_RecA1[ih]->GetXaxis()->FindBin(MassLow);
    binhi1 =diMuonsInvMass_RecA1[ih]->GetXaxis()->FindBin(MassHigh);
    double binwidth=diMuonsInvMass_RecA1[ih]->GetBinWidth(1);

    rec_pt[ih] = diMuonsInvMass_RecA1[ih]->IntegralAndError(binlow1, binhi1,recError);
    rec_ptError[ih]= recError;

    // Eff including T&P
    recTP_pt[ih] = diMuonsInvMass_RecATP1[ih]->IntegralAndError(binlow1, binhi1,recError);
    recTP_ptError[ih]= recError;

    // For Trigger Eff.
    recT_pt[ih] = diMuonsInvMass_RecT1[ih]->IntegralAndError(binlow1, binhi1,recError);
    recT_ptError[ih]= recError;

    // For muon ID Eff.
    recI_pt[ih] = diMuonsInvMass_RecI1[ih]->IntegralAndError(binlow1, binhi1,recError);
    recI_ptError[ih]= recError;

    //yield by full integral
    //    cout<<"Rec dimuons from mass "<<diMuonsInvMass_RecA1[ih]->Integral(1,100)<<endl; 
    //    cout<<"Rec diMuons from Pt histo "<<diMuonsPt_RecA1[ih]->Integral(1,100)<<endl;  
   
    //Cal eff 
    Eff[ih] = rec_pt[ih]/gen_pt[ih]; 
    Err_Eff[ih]= RError(rec_pt[ih],rec_ptError[ih],gen_pt[ih],gen_ptError[ih]);  

    //Cal eff 
    EffTP[ih] = recTP_pt[ih]/rec_pt[ih]; 
    Err_EffTP[ih]= RError(recTP_pt[ih],recTP_ptError[ih],rec_pt[ih],rec_ptError[ih]);  

    //Cal acc 
    Acc[ih] = gen_pt[ih]/gen_ptAll[ih]; 
    Err_Acc[ih]= RError(gen_pt[ih],gen_ptError[ih],gen_ptAll[ih],gen_ptAllError[ih]);  

    double te1, te2;
    if(genacc==1) {
      if(iSpec==1 || iSpec==2) {
	outFileGen<< "  " <<setprecision(4)<<pt_bound[ih]<<"  "<<setprecision(4)<<pt_bound[ih+1] << "  " << Acc[ih] << "  "<< Err_Acc[ih] << endl;
      }
      if(iSpec==3) {
	outFileGen<< "  " <<setprecision(4)<<pt_bound[ih]<<"  "<<setprecision(4)<<pt_bound[ih+1] << "  " << Acc[0] << "  "<< Err_Acc[0] << endl;
      }
    }
    if(genacc==2) {
      inFileGen >> te1 >> te2 >> Acc[ih] >>  Err_Acc[ih];
    }

    //Cal Trigger Eff
    Trig[ih] = rec_pt[ih]/recT_pt[ih]; 
    Err_Trig[ih]= RError(rec_pt[ih],rec_ptError[ih],recT_pt[ih],recT_ptError[ih]);  

    //Cal Trigger Eff
    muEff[ih] = rec_pt[ih]/recI_pt[ih]; 
    Err_muEff[ih]= RError(rec_pt[ih],rec_ptError[ih],recI_pt[ih],recI_ptError[ih]);  

    //Cal acc*eff 
    AccEff[ih] =  Acc[ih]*Eff[ih];
    Err_AccEff[ih]= PError(Acc[ih],Err_Acc[ih],Eff[ih],Err_Eff[ih]);  

    //Cal acc*eff 
    AccEffTP[ih] =  AccEff[ih]*EffTP[ih];
    Err_AccEffTP[ih]= PError(AccEff[ih],Err_AccEff[ih],EffTP[ih],Err_EffTP[ih]);  


    //    cout<<"UpsilonYield by CB yield det: "<< UpsilonYield << endl;
    //    cout<<"Upsilon Yield by Function integral: "<< GAUSPOL->Integral(MassLow,MassHigh)/binwidth <<endl;
    cout<< " //////////////////////////////////// " << endl;
    cout<< pt_bound[ih] << " - " << pt_bound[ih+1] << endl; 
    cout<<" gen_pt_All  "<<  gen_ptAll[ih] <<" +- "<< gen_ptAllError[ih]<<endl;
    cout<<" gen_pt  "<<  gen_pt[ih] <<" +- "<< gen_ptError[ih]<<endl;
    cout<<" rec_pt  "<<  rec_pt[ih] <<" +- "<< rec_ptError[ih]<<endl;
    cout<<" recT_pt "<< recT_pt[ih] <<" +- "<< recT_ptError[ih]<<endl;
    cout<<" recI_pt "<< recI_pt[ih] <<" +- "<< recI_ptError[ih]<<endl;

    cout<<"Upsilon eff "<< Eff[ih]<<" +- "<<Err_Eff[ih]<< endl;

    cout<<"Upsilon eff with T&P"<< EffTP[ih]<<" +- "<<Err_EffTP[ih]<< endl;

    if(iSpec==1 || iSpec ==2){
      dataFile<< "  " <<setprecision(3)<<fixed<<pt_bound[ih]<<"-"<<pt_bound[ih+1]<<" &  "<<Eff[ih]<<"  $\\pm$  "<<Err_Eff[ih]    <<"  &   "<<Acc[ih]<<"  $\\pm$  "<<Err_Acc[ih] <<"  &   "<<AccEff[ih]<<"  $\\pm$  "<<Err_AccEff[ih] <<" &  "<<EffTP[ih]<<"  $\\pm$  "<<Err_EffTP[ih] <<" &  "<<AccEffTP[ih]<<"  $\\pm$  "<<Err_AccEffTP[ih]  << " \\\ "    << endl;
    }
    if(iSpec==3){
      dataFile<< "  " <<setprecision(3)<<fixed<<pt_bound[ih]*2.5<<"-"<<pt_bound[ih+1]*2.5<<" &  "<<Eff[ih]<<"  $\\pm$  "<<Err_Eff[ih]    <<"  &   "<<Acc[ih]<<"  $\\pm$  "<<Err_Acc[ih] <<"  &   "<<AccEff[ih]<<"  $\\pm$  "<<Err_AccEff[ih] <<" &  "<<EffTP[ih]<<"  $\\pm$  "<<Err_EffTP[ih] <<" &  "<<AccEffTP[ih]<<"  $\\pm$  "<<Err_AccEffTP[ih]  << " \\\ "    << endl;
    }
  }

  dataFile<<endl<<endl; 
  dataFile.close();
  
  TFile *outfile;
  if(pp==1) {
    if(YS==1){
      if(iSpec==1)outfile =new TFile("EffUpsilon1spp_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffUpsilon1spp_Rap.root","Recreate");
    }
    if(YS==2){
      if(iSpec==1)outfile =new TFile("EffUpsilon2spp_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffUpsilon2spp_Rap.root","Recreate");
    }
    if(YS==3){
      if(iSpec==1)outfile =new TFile("EffUpsilon3spp_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffUpsilon3spp_Rap.root","Recreate");
    }
    if(YS==4){
      if(iSpec==1)outfile =new TFile("EffJpsipp_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffJpsipp_Rap.root","Recreate");
    }
  }

  if(pp==2) {
    if(YS==1){
      if(iSpec==1)outfile =new TFile("EffUpsilon1sPbPb_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffUpsilon1sPbPb_Rap.root","Recreate");
      if(iSpec==3)outfile =new TFile("EffUpsilon1sPbPb_Cen.root","Recreate");
    }
    if(YS==2){
      if(iSpec==1)outfile =new TFile("EffUpsilon2sPbPb_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffUpsilon2sPbPb_Rap.root","Recreate");
      if(iSpec==3)outfile =new TFile("EffUpsilon2sPbPb_Cen.root","Recreate");
    }
    if(YS==3){
      if(iSpec==1)outfile =new TFile("EffUpsilon3sPbPb_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffUpsilon3sPbPb_Rap.root","Recreate");
      if(iSpec==3)outfile =new TFile("EffUpsilon3sPbPb_Cen.root","Recreate");
    }
    if(YS==4){
      if(iSpec==1)outfile =new TFile("EffJpsiPbPb_Pt.root","Recreate");
      if(iSpec==2)outfile =new TFile("EffJpsiPbPb_Rap.root","Recreate");
      if(iSpec==3)outfile =new TFile("EffJpsiPbPb_Cen.root","Recreate");
    }
  }

  TGraphErrors *Eff_Upsilon = new TGraphErrors(Nptbin, PT, Eff, mom_err, Err_Eff);
  Eff_Upsilon->SetName("Eff_Upsilon");
  Eff_Upsilon->SetMarkerStyle(21);
  Eff_Upsilon->SetMarkerColor(2);
  Eff_Upsilon->GetYaxis()->SetTitle("Reconstruction Efficiency");
 
  if(iSpec==1) Eff_Upsilon->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  if(iSpec==2) Eff_Upsilon->GetXaxis()->SetTitle("rapidity");
  if(iSpec==3) Eff_Upsilon->GetXaxis()->SetTitle("bin");
  Eff_Upsilon->GetYaxis()->SetRangeUser(0,1.0);

  TLegend *legend_GP = new TLegend( 0.70,0.79,0.80,0.89);
  legend_GP->SetBorderSize(0);
  legend_GP->SetFillStyle(0);
  legend_GP->SetFillColor(0);
  legend_GP->SetTextSize(0.032);
  legend_GP->AddEntry(Eff_Upsilon,"PythiaEvtGen", "P");
  
  new TCanvas;
  Eff_Upsilon->Draw("AP");
  legend_GP->Draw("Same");
  Eff_Upsilon->Write();

  if(iSpec==1){ gPad->Print("plots/Eff_Upsilon_PtMC.pdf");}
  if(iSpec==1){ gPad->Print("plots/Eff_Upsilon_PtMC.png");}
  if(iSpec==2){ gPad->Print("plots/Eff_Upsilon_RapMC.pdf");}
  if(iSpec==2){ gPad->Print("plots/Eff_Upsilon_RapMC.png");}
  if(iSpec==3){ gPad->Print("plots/Eff_Upsilon_CentMC.pdf");}
  if(iSpec==3){ gPad->Print("plots/Eff_Upsilon_CentMC.png");}

  // Trigger efficiency
  TGraphErrors *Trig_Upsilon = new TGraphErrors(Nptbin, PT, Trig, mom_err, Err_Trig);
  Trig_Upsilon->SetName("Trig_Upsilon");
  Trig_Upsilon->SetMarkerStyle(21);
  Trig_Upsilon->SetMarkerColor(2);
  Trig_Upsilon->GetYaxis()->SetTitle("Trigger match Efficiency");
 
  if(iSpec==1) Trig_Upsilon->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  if(iSpec==2) Trig_Upsilon->GetXaxis()->SetTitle("rapidity");
  if(iSpec==3) Trig_Upsilon->GetXaxis()->SetTitle("bin");
  Trig_Upsilon->GetYaxis()->SetRangeUser(0,1.0);

  new TCanvas;
  Trig_Upsilon->Draw("AP");
  legend_GP->Draw("Same");
  Trig_Upsilon->Write();
  if(iSpec==1){ gPad->Print("plots/Trig_Upsilon_PtMC.pdf");}
  if(iSpec==1){ gPad->Print("plots/Trig_Upsilon_PtMC.png");}
  if(iSpec==2){ gPad->Print("plots/Trig_Upsilon_RapMC.pdf");}
  if(iSpec==2){ gPad->Print("plots/Trig_Upsilon_RapMC.png");}
  if(iSpec==3){ gPad->Print("plots/Trig_Upsilon_CentMC.pdf");}
  if(iSpec==3){ gPad->Print("plots/Trig_Upsilon_CentMC.png");}

  // muon ID eff
  TGraphErrors *muEff_Upsilon = new TGraphErrors(Nptbin, PT, muEff, mom_err, Err_muEff);
  muEff_Upsilon->SetName("muEff_Upsilon");
  muEff_Upsilon->SetMarkerStyle(21);
  muEff_Upsilon->SetMarkerColor(2);
  muEff_Upsilon->GetYaxis()->SetTitle("muon ID Efficiency");
 
  if(iSpec==1) muEff_Upsilon->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  if(iSpec==2) muEff_Upsilon->GetXaxis()->SetTitle("rapidity");
  if(iSpec==3) muEff_Upsilon->GetXaxis()->SetTitle("bin");
  muEff_Upsilon->GetYaxis()->SetRangeUser(0,1.0);

  new TCanvas;
  muEff_Upsilon->Draw("AP");
  legend_GP->Draw("Same");
  muEff_Upsilon->Write();

  if(iSpec==1){ gPad->Print("plots/muEff_Upsilon_PtMC.pdf");}
  if(iSpec==1){ gPad->Print("plots/muEff_Upsilon_PtMC.png");}
  if(iSpec==2){ gPad->Print("plots/muEff_Upsilon_RapMC.pdf");}
  if(iSpec==2){ gPad->Print("plots/muEff_Upsilon_RapMC.png");}
  if(iSpec==3){ gPad->Print("plots/muEff_Upsilon_CentMC.pdf");}
  if(iSpec==3){ gPad->Print("plots/muEff_Upsilon_CentMC.png");}

  // Acceptance
  TGraphErrors *Acc_Upsilon = new TGraphErrors(Nptbin, PT, Acc, mom_err, Err_Acc);
  Acc_Upsilon->SetName("Acc_Upsilon");
  Acc_Upsilon->SetMarkerStyle(21);
  Acc_Upsilon->SetMarkerColor(2);
  Acc_Upsilon->GetYaxis()->SetTitle("Acceptance");
 
  if(iSpec==1) Acc_Upsilon->GetXaxis()->SetTitle("#Upsilon pT (GeV/c)");
  if(iSpec==2) Acc_Upsilon->GetXaxis()->SetTitle("#Upsilon rapidity");
  if(iSpec==3) Acc_Upsilon->GetXaxis()->SetTitle("bin");
  Acc_Upsilon->GetYaxis()->SetRangeUser(0,1.0);

  new TCanvas;
  Acc_Upsilon->Draw("AP");
  legend_GP->Draw("Same");
  Acc_Upsilon->Write();

  if(iSpec==1){ gPad->Print("plots/Acc_Upsilon_PtMC.pdf");}
  if(iSpec==1){ gPad->Print("plots/Acc_Upsilon_PtMC.png");}
  if(iSpec==2){ gPad->Print("plots/Acc_Upsilon_RapMC.pdf");}
  if(iSpec==2){ gPad->Print("plots/Acc_Upsilon_RapMC.png");}
  if(iSpec==3){ gPad->Print("plots/Acc_Upsilon_CentMC.pdf");}
  if(iSpec==3){ gPad->Print("plots/Acc_Upsilon_CentMC.png");}
 

  // Acceptance*Efficiency
  TGraphErrors *AccEff_Upsilon = new TGraphErrors(Nptbin, PT, AccEff, mom_err, Err_AccEff);
  AccEff_Upsilon->SetName("AccEff_Upsilon");
  AccEff_Upsilon->SetMarkerStyle(21);
  AccEff_Upsilon->SetMarkerColor(2);
  AccEff_Upsilon->GetYaxis()->SetTitle("Acc*Eff");
 
  if(iSpec==1) AccEff_Upsilon->GetXaxis()->SetTitle("#Upsilon pT (GeV/c)");
  if(iSpec==2) AccEff_Upsilon->GetXaxis()->SetTitle("#Upsilon rapidity");
  if(iSpec==3) AccEff_Upsilon->GetXaxis()->SetTitle("bin");
  AccEff_Upsilon->GetYaxis()->SetRangeUser(0,1.0);

  new TCanvas;
  AccEff_Upsilon->Draw("AP");
  legend_GP->Draw("Same");
  AccEff_Upsilon->Write();

  if(iSpec==1){ gPad->Print("plots/AccEff_Upsilon_PtMC.pdf");}
  if(iSpec==1){ gPad->Print("plots/AccEff_Upsilon_PtMC.png");}
  if(iSpec==2){ gPad->Print("plots/AccEff_Upsilon_RapMC.pdf");}
  if(iSpec==2){ gPad->Print("plots/AccEff_Upsilon_RapMC.png");}
  if(iSpec==3){ gPad->Print("plots/AccEff_Upsilon_CentMC.pdf");}
  if(iSpec==3){ gPad->Print("plots/AccEff_Upsilon_CentMC.png");}
 


  MuonPtEta1_Eff -> Divide(MuonPtEta_Rec, MuonPtEta_Gen);
  new TCanvas;
  MuonPtEta_Gen -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta_Gen -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta_Gen->Draw("colz");
   new TCanvas;
  MuonPtEta1_Gen -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta1_Gen -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta1_Gen->Draw("colz");
  new TCanvas;
  MuonPtEta2_Gen -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta2_Gen -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta2_Gen->Draw("colz");
  new TCanvas;
  MuonPtEta3_Gen -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta3_Gen -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta3_Gen->Draw("colz");
   new TCanvas;
  MuonPtEta_Rec -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta_Rec -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta_Rec->Draw("colz");
 new TCanvas;
  MuonPtEta1_Rec -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta1_Rec -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta1_Rec->Draw("colz");
  new TCanvas;
  MuonPtEta2_Rec -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta2_Rec -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta2_Rec->Draw("colz");
  new TCanvas;
  MuonPtEta3_Rec -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta3_Rec -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta3_Rec->Draw("colz");
  new TCanvas;
  MuonPtEta1_Eff -> GetXaxis()->SetTitle("#eta^{#mu}");
  MuonPtEta1_Eff -> GetYaxis()->SetTitle("p_{T}^{#mu}(GeV/c)");
  MuonPtEta1_Eff->Draw("colz");

  MuonPtEta1_Eff ->Write();

  diMuonsPt_Gen->Write();
  diMuonsRap_Gen->Write();
  
  MuonPt_Gen->Write();
  MuonEta_Gen->Write();
  MuonPhi_Gen->Write();

  outfile->Write();
  outfile->Close();

  cout<< " ////////////////////////////////////////////////////////" << endl;
  cout<<"Upsilon eff Total "<< Eff_Total<<" error "<<Err_Eff_Total<< endl;
  cout<< " ////////////////////////////////////////////////////////" << endl;

}


bool IsAccept(Double_t pt, Double_t eta)
{
  return (fabs(eta) < 2.4);
  return (fabs(eta) < 2.4 &&
  	  (    ( fabs(eta) < 1.0 && pt >= 3.4 ) ||
  	       (  1.0 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 5.8-2.4*fabs(eta) ) ||
  	       (  1.5 <= fabs(eta) && pt >= 3.3667 - 7.0/9.0*fabs(eta)) ));
}

double FindCenWeight(int Bin)
{
  float NCollArray[40]=
  {1747.49
   ,1566.92
   ,1393.97
   ,1237.02
   ,1095.03
   ,979.836
   ,863.228
   ,765.968
   ,677.894
   ,594.481
   ,522.453
   ,456.049
   ,399.178
   ,347.174
   ,299.925
   ,258.411
   ,221.374
   ,188.676
   ,158.896
   ,135.117
   ,112.481
   ,93.5697
   ,77.9192
   ,63.2538
   ,52.0938
   ,42.3553
   ,33.7461
   ,27.3213
   ,21.8348
   ,17.1722
   ,13.5661
   ,10.6604
   ,8.31383
   ,6.37662
   ,5.12347
   ,3.73576
   ,3.07268
   ,2.41358
   ,2.10707
   ,1.76851};/*double NCollArray[50]={1563.03, 1370.41, 1204.05, 1063.45, 943.5,
			 834.12, 736.223, 654.913, 576.466, 507.757, 443.05, 386.802, 334.48,
			 290.097, 247.779, 211.762, 179.834, 153.509, 127.75, 106.59, 88.1189,
			 72.3836, 59.1049, 47.3574, 37.7951, 30.1705, 23.6861, 18.6918, 14.2287,
			 10.9705, 8.76148, 6.57459, 5.01557, 3.78525, 2.9123, 2.12377, 1.5,
			 0.922951, 0.581967, 0.503279,};*/
  return(NCollArray[Bin]);
}


//Ratio Error
double RError(double A, double eA, double B, double eB){
  double f=A/B;
  double fA=eA/A;
  double fB=eB/B;
  double eR=  f*sqrt( (fA*fA + fB*fB )) ;
  return eR;
}


//Product Error
double PError(double A, double eA, double B, double eB){
  double f=A*B;
  double fA=eA/A;
  double fB=eB/B;
  double eR=  f*sqrt( (fA*fA + fB*fB )) ;
  return eR;
}

///pbpb
double getEffEtaPbPb(Double_t pt, Double_t eta)
{
  if(fabs(eta)<1.6) return ( 0.9988*TMath::Erf( (pt-1.322)/2.688)/TMath::Erf( (pt-1.796)/2.516) );

  //  if(1.6<fabs(eta)) return ( 1.04*TMath::Erf( (pt-1.41)/3.629)/TMath::Erf( (pt-1.571)/3.752) );
  if(1.6<fabs(eta)) return (0.8299*TMath::Erf((pt-1.2785)/1.8833))/(0.7810*TMath::Erf((pt-1.3609)/2.1231));

  return 1.0;
}

///pp 
double getEffEtaPP(Double_t pt, Double_t eta)
{
  if(fabs(eta)<1.6) return ( 0.998474*TMath::Erf( (pt-1.63555)/(0.797933) ) / TMath::Erf( (pt-0.222866)/(2.95593) ) );

  //  if(1.6<fabs(eta)) return ( 0.6977*TMath::Erf( (pt-1.61309)/(-1.03684) ) + 1.76794);
  if(1.6<fabs(eta)) return (0.7788*TMath::Erf((pt-1.1903)/1.9880))/(0.7364*TMath::Erf((pt-1.2538)/2.2530));

  return 1.0;
}


