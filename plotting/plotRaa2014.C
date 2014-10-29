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

#include "data_raa.h"

// miscellaneous  
#include <fstream>
#include <iostream>

using namespace std;

const bool plotCS =false ;//fully corrected
const bool plotUncorrected = false;
const bool plotFiducial=false;//fiducial plots(not corrected for acc)
const bool plotEffOrAcc=false;//control plots (ratio of efficiencies etc)
const bool plotRAA=true;// centrality, transverse momentum, rapidity,
const bool plotPA=false;
const bool plotTNP=false;
const bool plot2010=true;
float computeRatio(float x, float y) ;
float computeRatioError(float x, float y, float xerr, float yerr);
void plot2010();
void combine_blue(double val1, double err1, double val2, double err2);
void plotRAA_uncorr();
void plotDoubleRatios();
void plotRaa2014()
{
  gROOT->Macro("../code/cm/logon.C+");//it all looks much nicer with this.
  //  gROOT->Macro("data_raa.h");

 

  double RapBinWidth = 4.8;
  double PtBinWidth = 20;
  float pt [nPtBins_2013] = {1.25, 3.75, 6.5, 10., 16.};
  float pt2014 [nPtBins_2014] = {1.25, 3.75, 6.5, 10., 16.,30};
  float pt2014e[nPtBins_2014] = {1.25, 1.25, 1.5, 2., 4.,10.};
  float pte[nPtBins_2013] = {1.25, 1.25, 1.5, 2., 4.};
  float deltaPt[nPtBins_2013]   = {2.5,2.5,3,4,8};
  float pt_2010 [nPtBins_2010] = {3.25,8.25,15.};
  float pte_2010[nPtBins_2010] = {3.25,1.75,5};
  float deltaPt_2010[nPtBins_2010]   = {6.5,3.5,10};
  float deltaRap[nRapBins_2013]  = {0.8,0.6,0.6,1,1.8};
  float deltaRapEven[nRapBins_2014] = {0.8,0.8,0.8,0.8,0.8,0.8};
  float deltaRap2010[nRapBins_2010] = {2.4,2.4};

  float R_A_1S_pt[nPtBins_2013]={};
  float R_A_1S_pte[nPtBins_2013]={};
  float R_A_1S_rap[nRapBins_2014]={};
  float R_A_1S_rape[nRapBins_2014]={}; 
  float R_e_1S_pt[nPtBins_2013]={};
  float R_e_1S_pte[nPtBins_2013]={};
  float R_e_1S_rap[nRapBins_2014]={};
  float R_e_1S_rape[nRapBins_2014]={};

  float CS1S_pp_rap[nRapBins_2013] = {}; // left this one to compare 2013 with 2014 (only difference: larger bin 1.5-2.4). /// not used for now.
  float CS1S_pp_rape[nRapBins_2013] = {};
  float CS1S_pp_pt[nPtBins_2013] = {};
  float CS1S_pp_pte[nPtBins_2013] = {};
  //large things
  float CS1S_pp_ptLarge[nPtBins_2010] = {};
  float CS1S_pp_pteLarge[nPtBins_2010] = {};
  float CS1S_aa_ptLarge[nPtBins_2010] = {};
  float CS1S_aa_pteLarge[nPtBins_2010] = {};
  //RAA pt 2010
  float RAA_1S_ptLarge[nPtBins_2010]={};
  float RAA_1S_pteLarge[nPtBins_2010]={};
  //
  float CS1S_pp_tnp_pt[nPtBins_2013] = {};
  float CS1S_pp_tnp_pte[nPtBins_2013] = {};
  float CS1S_pp_rap2014[nRapBins_2014] = {};
  float CS1S_pp_rap2014e[nRapBins_2014] = {};
  float CS1S_pp_tnp_rap2014[nRapBins_2014] = {};
  float CS1S_pp_tnp_rap2014e[nRapBins_2014] = {};
  float CS1S_aa_pt[nPtBins_2013] = {};
  float CS1S_aa_pte[nPtBins_2013] = {}; //not sure i'll use it for the moment
  float CS1S_aa_rap[nRapBins_2013] = {};
  float CS1S_aa_rape[nRapBins_2013] = {};
  float CS1S_aa_rap2014[nRapBins_2014] = {};
  float CS1S_aa_rap2014e[nRapBins_2014] = {};
  float CS1S_aa_tnp_pt[nPtBins_2013] = {};
  float CS1S_aa_tnp_pte[nPtBins_2013] = {};
  float CS1S_aa_tnp_rap[nRapBins_2013] = {};
  float CS1S_aa_tnp_rape[nRapBins_2013] = {};
  float CS1S_aa_tnp_rap2014[nRapBins_2014] = {};
  float CS1S_aa_tnp_rap2014e[nRapBins_2014] = {};

  float CS1S_pp_ptFiducial[nPtBins_2013]={}; //
  float CS1S_aa_ptFiducial[nPtBins_2013]={}; //
  float CS1S_pp_rap2014Fiducial[nRapBins_2014]={}; //
  float CS1S_aa_rapFiducial[nRapBins_2013]={};//

  float CS2S_pp_pt2013Fiducial[nPtBins_2013]={};//
  // float CS2S_pp_ptFiducial[nPtBins_2010]={};
  // float CS2S_aa_ptFiducial[nPtBins_2010]={};
  // float CS2S_aa_rapFiducial[nRapBins_2010]={};
  // float CS2S_pp_rapFiducial[nRapBins_2010]={};
  float CS2S_pp_rap2014Fiducial[nRapBins_2014]={};

  float CS1S_pp_ptFiduciale[nPtBins_2013]={};
  float CS1S_aa_ptFiduciale[nPtBins_2013]={};
  float CS1S_pp_rap2014Fiduciale[nRapBins_2014]={};
  float CS1S_aa_rapFiduciale[nRapBins_2014]={};

  float CS2S_pp_pt2013Fiduciale[nPtBins_2013]={};
  // float CS2S_pp_ptFiduciale[nPtBins_2010]={};
  // float CS2S_aa_ptFiduciale[nPtBins_2010]={};
  // float CS2S_aa_rapFiduciale[nRapBins_2010]={};
  // float CS2S_pp_rapFiduciale[nRapBins_2010]={};
  float CS2S_pp_rap2014Fiduciale[nRapBins_2014]={};

  float CS2S_pp_pt[nPtBins_2010] = {};
  float CS2S_pp_pt2013[nPtBins_2013] ={};
  float CS2S_pp_pte[nPtBins_2010] = {};
  float CS2S_pp_pt2013e[nPtBins_2013] ={};
  float CS2S_pp_tnp_pt[nPtBins_2010] = {};
  float CS2S_pp_tnp_pt2013[nPtBins_2013] ={};
  float CS2S_pp_tnp_pte[nPtBins_2010] = {};
  float CS2S_pp_tnp_pt2013e[nPtBins_2013] ={};
  float CS2S_pp_rap[nRapBins_2010] = {};
  float CS2S_pp_rap2014[nRapBins_2014] = {};
  float CS2S_pp_rape[nRapBins_2010] = {};
  float CS2S_pp_rap2014e[nRapBins_2014] = {};
  float CS2S_pp_tnp_rap[nRapBins_2010] = {};
  float CS2S_pp_tnp_rape[nRapBins_2010] = {};
  float CS2S_pp_tnp_rap2014[nRapBins_2014] = {};
  float CS2S_pp_tnp_rap2014e[nRapBins_2014] = {}; 
  float CS2S_aa_pt[nPtBins_2010] = {};
  float CS2S_aa_pte[nPtBins_2010] = {};
  float CS2S_aa_rap[nRapBins_2010] = {};
  float CS2S_aa_rape[nRapBins_2010] = {};
  float CS2S_aa_tnp_pt[nPtBins_2010] = {};
  float CS2S_aa_tnp_pte[nPtBins_2010] = {};
  float CS2S_aa_tnp_rap[nRapBins_2010] = {};
  float CS2S_aa_tnp_rape[nRapBins_2010] = {};

  
  float CS3S_pp_pt2013[nPtBins_2013] ={};
  float CS3S_pp_pt2013e[nPtBins_2013] ={};
  float CS3S_pp_rap2014[nRapBins_2014] = {};
  float CS3S_pp_rap2014e[nRapBins_2014] = {};
  float CS3S_pp_tnp_pt2013[nPtBins_2013] ={};
  float CS3S_pp_tnp_pt2013e[nPtBins_2013] ={};
  float CS3S_pp_tnp_rap2014[nRapBins_2014] = {};
  float CS3S_pp_tnp_rap2014e[nRapBins_2014] = {};

  float RAA_1S_pt[nPtBins_2013]={};
  float RAA_1S_rap[nRapBins_2014]={};
  float RAA_1S_pte[nPtBins_2013]={};
  float RAA_1S_rape[nRapBins_2014]={};

  float RAA_2S_pt[nPtBins_2010]={};
  float RAA_2S_rap[nRapBins_2010]={};
  float RAA_2S_pte[nPtBins_2010]={};
  float RAA_2S_rape[nRapBins_2010]={};
  //tnp correction
  float RAA_1S_tnp_pt[nPtBins_2013]={};
  float RAA_1S_tnp_rap[nRapBins_2014]={};
  float RAA_1S_tnp_pte[nPtBins_2013]={};
  float RAA_1S_tnp_rape[nRapBins_2014]={};

  float RAA_2S_tnp_pt[nPtBins_2010]={};
  float RAA_2S_tnp_rap[nRapBins_2010]={};
  float RAA_2S_tnp_pte[nPtBins_2010]={};
  float RAA_2S_tnp_rape[nRapBins_2010]={};

  for(int i = 0 ; i<nPtBins_2013 ; i++)
    {
      //invariant yields.
      //1S
      CS1S_pp_pt[i]= computeRatio( N1S_pp_pt3p5[i] , Ae_1S_pythia_pt[i] );  
      CS1S_aa_pt[i]= computeRatio( N1S_aa_pt3p5[i] , Ae_1S_pyquen_pt[i] );   
      CS1S_pp_pte[i] = computeRatioError( N1S_pp_pt3p5[i] , Ae_1S_pythia_pt[i], N1S_pp_pt3p5e[i] , Ae_1S_pythia_pte[i]);
      CS1S_aa_pte[i] = computeRatioError(  N1S_aa_pt3p5[i] , Ae_1S_pyquen_pt[i] ,  N1S_aa_pt3p5e[i] , Ae_1S_pyquen_pte[i] );
      //2S
      CS2S_pp_pt2013[i]= computeRatio( N2S_pp_pt4_2013[i] , Ae_2S_pythia_pt2013[i] ); 
      CS2S_pp_pt2013e[i] = computeRatioError( N2S_pp_pt4_2013[i] , Ae_2S_pythia_pt2013[i] , N2S_pp_pt4_2013e[i] , Ae_2S_pythia_pt2013e[i] );
      //3S
      CS3S_pp_pt2013[i]= computeRatio( N3S_pp_pt4_2013[i] , Ae_3S_pythia_pt2013[i] ); 
      CS3S_pp_pt2013e[i] = computeRatioError( N3S_pp_pt4_2013[i] , Ae_3S_pythia_pt2013[i] , N3S_pp_pt4_2013e[i] , Ae_3S_pythia_pt2013e[i] );
      //comparison with tnp
      CS1S_pp_tnp_pt[i]= computeRatio( N1S_pp_pt3p5[i] , Aet_1S_pythia_pt[i] ); 
      CS1S_aa_tnp_pt[i]= computeRatio( N1S_aa_pt3p5[i] , Aet_1S_pyquen_pt[i] );
      CS1S_pp_tnp_pte[i] = computeRatioError( N1S_pp_pt3p5[i] , Aet_1S_pythia_pt[i], N1S_pp_pt3p5e[i] , Aet_1S_pythia_pte[i]);
      CS1S_aa_tnp_pte[i] = computeRatioError(  N1S_aa_pt3p5[i] , Aet_1S_pyquen_pt[i] ,  N1S_aa_pt3p5e[i] , Aet_1S_pyquen_pte[i] );
      CS2S_pp_tnp_pt2013[i]= computeRatio( N2S_pp_pt4_2013[i] , Aet_2S_pythia_pt2013[i] ); 
      CS2S_pp_tnp_pt2013e[i] = computeRatioError( N2S_pp_pt4_2013[i] , Aet_2S_pythia_pt2013[i] , N2S_pp_pt4_2013e[i] , Aet_2S_pythia_pt2013e[i] );        
      CS3S_pp_tnp_pt2013[i]= computeRatio( N3S_pp_pt4_2013[i] , Aet_3S_pythia_pt2013[i] ); 
      CS3S_pp_tnp_pt2013e[i] = computeRatioError( N3S_pp_pt4_2013[i] , Aet_3S_pythia_pt2013[i] , N3S_pp_pt4_2013e[i] , Aet_3S_pythia_pt2013e[i] );

    
    
      

      //fiducial shit
      CS2S_pp_pt2013Fiducial[i]= computeRatio( N2S_pp_pt4_2013[i] , e_2S_pythia_pt2013[i] );      
      CS2S_pp_pt2013Fiduciale[i] = computeRatioError( N2S_pp_pt4_2013[i] , e_2S_pythia_pt2013[i] , N2S_pp_pt4_2013e[i] , e_2S_pythia_pt2013e[i] );
      CS1S_aa_ptFiducial[i]= computeRatio( N1S_aa_pt3p5[i] , e_1S_pyquen_pt[i] );
      CS1S_pp_ptFiducial[i]= computeRatio( N1S_pp_pt3p5[i] , e_1S_pythia_pt3p5[i] );
      CS1S_aa_pt[i]= computeRatio( N1S_aa_pt3p5[i] , Ae_1S_pyquen_pt[i] ); 
      CS1S_pp_ptFiducial[i]= computeRatio( N1S_pp_pt3p5[i] , e_1S_pythia_pt3p5[i] ); 
      CS1S_aa_ptFiducial[i]= computeRatio( N1S_aa_pt3p5[i] , e_1S_pyquen_pt[i] );
      CS1S_pp_ptFiduciale[i] = computeRatioError( N1S_pp_pt3p5[i] , e_1S_pythia_pt3p5[i], N1S_pp_pt3p5e[i] , e_1S_pythia_pt3p5e[i]);   
      CS1S_aa_ptFiduciale[i] = computeRatioError(  N1S_aa_pt3p5[i] , e_1S_pyquen_pt[i] ,  N1S_aa_pt3p5e[i] , e_1S_pyquen_pte[i] );
      CS1S_pp_ptFiducial[i]=CS1S_pp_ptFiducial[i]/L_pp;
      CS1S_aa_ptFiducial[i]=CS1S_aa_ptFiducial[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_ptFiduciale[i]=CS1S_pp_ptFiduciale[i]/L_pp;
      CS1S_aa_ptFiduciale[i]=CS1S_aa_ptFiduciale[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_ptFiducial[i]=CS1S_pp_ptFiducial[i]/deltaPt[i];
      CS1S_aa_ptFiducial[i]=CS1S_aa_ptFiducial[i]/deltaPt[i];
      CS1S_pp_ptFiduciale[i]=CS1S_pp_ptFiduciale[i]/deltaPt[i];
      CS1S_aa_ptFiduciale[i]=CS1S_aa_ptFiduciale[i]/deltaPt[i];    
      CS2S_pp_pt2013Fiducial[i]=CS2S_pp_pt2013Fiducial[i]/L_pp;
      CS2S_pp_pt2013Fiduciale[i]=CS2S_pp_pt2013Fiduciale[i]/L_pp;
      CS2S_pp_pt2013Fiducial[i]=CS2S_pp_pt2013Fiducial[i]/deltaPt[i];
      CS2S_pp_pt2013Fiduciale[i]=CS2S_pp_pt2013Fiduciale[i]/deltaPt[i];
      //good ones
      CS1S_pp_pt[i]=CS1S_pp_pt[i]/(L_pp*RapBinWidth*deltaPt[i]);    
      CS1S_aa_pt[i]=CS1S_aa_pt[i]/(N_MB_corr * T_AA_b*(RapBinWidth*deltaPt[i]));
      CS1S_pp_pte[i]=CS1S_pp_pte[i]/(L_pp*RapBinWidth*deltaPt[i]);
      CS1S_aa_pte[i]=CS1S_aa_pte[i]/(N_MB_corr * T_AA_b *(RapBinWidth*deltaPt[i]));
      CS2S_pp_pt2013[i]=CS2S_pp_pt2013[i]/(RapBinWidth*deltaPt[i]*L_pp);
      CS2S_pp_pt2013e[i]=CS2S_pp_pt2013e[i]/(RapBinWidth*deltaPt[i]*L_pp);
      CS3S_pp_pt2013[i]=CS2S_pp_pt2013[i]/(RapBinWidth*deltaPt[i]*L_pp);
      CS3S_pp_pt2013e[i]=CS2S_pp_pt2013e[i]/(RapBinWidth*deltaPt[i]*L_pp);
      //tnp correction
      CS1S_pp_tnp_pt[i]=CS1S_pp_tnp_pt[i]/(L_pp*RapBinWidth*deltaPt[i]);    
      CS1S_aa_tnp_pt[i]=CS1S_aa_tnp_pt[i]/(N_MB_corr * T_AA_b*(RapBinWidth*deltaPt[i]));
      CS1S_pp_tnp_pte[i]=CS1S_pp_tnp_pte[i]/(L_pp*RapBinWidth*deltaPt[i]);
      CS1S_aa_tnp_pte[i]=CS1S_aa_tnp_pte[i]/(N_MB_corr * T_AA_b *(RapBinWidth*deltaPt[i]));
      CS2S_pp_tnp_pt2013[i]=CS2S_pp_tnp_pt2013[i]/(RapBinWidth*deltaPt[i]*L_pp);
      CS2S_pp_tnp_pt2013e[i]=CS2S_pp_tnp_pt2013e[i]/(RapBinWidth*deltaPt[i]*L_pp);
      CS3S_pp_tnp_pt2013[i]=CS2S_pp_tnp_pt2013[i]/(RapBinWidth*deltaPt[i]*L_pp);
      CS3S_pp_tnp_pt2013e[i]=CS2S_pp_tnp_pt2013e[i]/(RapBinWidth*deltaPt[i]*L_pp);
      RAA_1S_pt[i]= computeRatio( CS1S_aa_pt[i] , CS1S_pp_pt[i]);
      RAA_1S_pte[i] =computeRatioError( CS1S_aa_pt[i] , CS1S_pp_pt[i], CS1S_aa_pte[i] , CS1S_pp_pte[i]);
      // cout<<"cs1S (ppPt, ppRap, aaPt, aaRap)"<< CS1S_pp_pt[i] <<", " <<CS1S_pp_rap2014[i] <<", " <<CS1S_aa_pt[i] <<", " <<CS1S_aa_rap[i] <<". " << endl;
      // cout <<"sigma(1S)_pp vs Pt ="<< endl;
      //invariant yields, corrected by taa and lumi corr.//////raaaaaaaa!
      RAA_1S_tnp_pt[i]= computeRatio( CS1S_aa_tnp_pt[i] , CS1S_pp_tnp_pt[i]);
      RAA_1S_tnp_pte[i] =computeRatioError( CS1S_aa_tnp_pt[i] , CS1S_pp_tnp_pt[i], CS1S_aa_tnp_pte[i] , CS1S_pp_tnp_pte[i]);
      R_A_1S_pt[i]=computeRatio(A_1S_pythia_pt3p5[i],A_1S_pyquen_pt[i]);
      R_A_1S_pte[i]=computeRatioError(A_1S_pythia_pt3p5[i],A_1S_pyquen_pt[i],A_1S_pythia_pt3p5e[i],A_1S_pyquen_pte[i]);
      R_e_1S_pt[i]=computeRatio(e_1S_pythia_pt3p5[i],e_1S_pyquen_pt[i]);
      R_e_1S_pte[i]=computeRatioError(e_1S_pythia_pt3p5[i],e_1S_pyquen_pt[i],e_1S_pythia_pt3p5e[i],e_1S_pyquen_pte[i]);

    }
 
  for(int i=0; i < nRapBins_2014; i++)
    {

      //fiducial shit
      CS1S_aa_rapFiducial[i]= computeRatio( N1S_aa_rap3p5_2014[i] , e_1S_pyquen_rap2014[i] );
      CS1S_pp_rap2014Fiduciale[i] = computeRatioError(  N1S_pp_rap3p5_2014[i] , e_1S_pythia_rap3p5[i] ,  N1S_pp_rap3p5_2014e[i] , e_1S_pythia_rap3p5e[i] );
      CS1S_pp_rap2014Fiducial[i]= computeRatio( N1S_pp_rap3p5_2014[i] , e_1S_pythia_rap3p5[i] );  
      CS1S_aa_rapFiduciale[i] = computeRatioError(   N1S_aa_rap3p5_2014[i] , e_1S_pyquen_rap2014[i] ,  N1S_aa_rap3p5_2014e[i] , e_1S_pyquen_rap2014e[i]);
      CS1S_pp_rap2014Fiducial[i]=CS1S_pp_rap2014Fiducial[i]/L_pp;
      CS1S_aa_rapFiducial[i]=CS1S_aa_rapFiducial[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_rap2014Fiduciale[i]=CS1S_pp_rap2014Fiduciale[i]/L_pp;
      CS1S_aa_rapFiduciale[i]=CS1S_aa_rapFiduciale[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_rap2014Fiducial[i]=CS1S_pp_rap2014Fiducial[i]/deltaRapEven[i];
      CS1S_aa_rapFiducial[i]=CS1S_aa_rapFiducial[i]/deltaRapEven[i];
      CS1S_pp_rap2014Fiduciale[i]=CS1S_pp_rap2014Fiduciale[i]/deltaRapEven[i];
      CS1S_aa_rapFiduciale[i]=CS1S_aa_rapFiduciale[i]/deltaRapEven[i];
      ///good ones
      CS1S_pp_rap2014[i]= computeRatio( N1S_pp_rap3p5_2014[i] , Ae_1S_pythia_rap2014[i] ); 
      CS1S_aa_rap[i]= computeRatio( N1S_aa_rap3p5_2014[i] , Ae_1S_pyquen_rap2014[i] );
      CS1S_pp_rap2014e[i] = computeRatioError(  N1S_pp_rap3p5_2014[i] , Ae_1S_pythia_rap2014[i] ,  N1S_pp_rap3p5_2014e[i] , Ae_1S_pythia_rap2014e[i] );
      CS1S_aa_rape[i] = computeRatioError(   N1S_aa_rap3p5_2014[i] , Ae_1S_pyquen_rap2014[i] ,  N1S_aa_rap3p5_2014e[i] , Ae_1S_pyquen_rap2014e[i]); 
      CS1S_pp_rap2014[i]=CS1S_pp_rap2014[i]/L_pp;
      CS1S_aa_rap[i]=CS1S_aa_rap[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_rap2014e[i]=CS1S_pp_rap2014e[i]/L_pp;
      CS1S_aa_rape[i]=CS1S_aa_rape[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_rap2014[i]=CS1S_pp_rap2014[i]/(deltaRapEven[i]);
      CS1S_aa_rap[i]=CS1S_aa_rap[i]/deltaRapEven[i];
      CS1S_pp_rap2014e[i]=CS1S_pp_rap2014e[i]/(deltaRapEven[i]);
      CS1S_aa_rape[i]=CS1S_aa_rape[i]/deltaRapEven[i];
      //2S
      CS2S_pp_rap2014[i]= computeRatio( N2S_pp_rap4_2014[i] , Ae_2S_pythia_rap2014[i] ); 
      CS2S_pp_rap2014e[i]= computeRatioError(  N2S_pp_rap4_2014[i] , Ae_2S_pythia_rap2014[i] ,  N2S_pp_rap4_2014e[i] , Ae_2S_pythia_rap2014e[i] );
      CS2S_pp_rap2014[i]=CS2S_pp_rap2014[i]/(deltaRapEven[i] * L_pp);
      CS2S_pp_rap2014e[i]=CS2S_pp_rap2014e[i]/(deltaRapEven[i] * L_pp);
      //3S
      CS3S_pp_rap2014[i]= computeRatio( N3S_pp_rap4_2014[i] , Ae_3S_pythia_rap2014[i] ); 
      CS3S_pp_rap2014e[i] = computeRatioError( N3S_pp_rap4_2014[i] , Ae_3S_pythia_rap2014[i] , N3S_pp_rap4_2014e[i] , Ae_3S_pythia_rap2014e[i] );
      CS3S_pp_rap2014[i]=      CS3S_pp_rap2014[i]/L_pp;
      CS3S_pp_rap2014e[i]= CS3S_pp_rap2014e[i]/L_pp;
      CS3S_pp_rap2014[i]=      CS3S_pp_rap2014[i]/(deltaRapEven[i]);
      CS3S_pp_rap2014e[i]= CS3S_pp_rap2014e[i]/L_pp;
      CS3S_pp_rap2014e[i]= CS3S_pp_rap2014e[i]/(deltaRapEven[i]);
      //TAG AND PROBE COMPARISON
      //1S
      CS1S_pp_tnp_rap2014[i]= computeRatio( N1S_pp_rap3p5_2014[i] , Aet_1S_pythia_rap2014[i] ); 
      CS1S_aa_tnp_rap[i]= computeRatio( N1S_aa_rap3p5_2014[i] , Aet_1S_pyquen_rap2014[i] );
      CS1S_pp_tnp_rap2014e[i] = computeRatioError(  N1S_pp_rap3p5_2014[i] , Aet_1S_pythia_rap2014[i] ,  N1S_pp_rap3p5_2014e[i] , Aet_1S_pythia_rap2014e[i] );
      CS1S_aa_tnp_rape[i] = computeRatioError(   N1S_aa_rap3p5_2014[i] , Aet_1S_pyquen_rap2014[i] ,  N1S_aa_rap3p5_2014e[i] , Aet_1S_pyquen_rap2014e[i]); 
      CS1S_pp_tnp_rap2014[i]=CS1S_pp_tnp_rap2014[i]/L_pp;
      CS1S_aa_tnp_rap[i]=CS1S_aa_tnp_rap[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_tnp_rap2014e[i]=CS1S_pp_tnp_rap2014e[i]/L_pp;
      CS1S_aa_tnp_rape[i]=CS1S_aa_tnp_rape[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_tnp_rap2014[i]=CS1S_pp_tnp_rap2014[i]/(deltaRapEven[i]);
      CS1S_aa_tnp_rap[i]=CS1S_aa_tnp_rap[i]/deltaRapEven[i];
      CS1S_pp_tnp_rap2014e[i]=CS1S_pp_tnp_rap2014e[i]/(deltaRapEven[i]);
      CS1S_aa_tnp_rape[i]=CS1S_aa_tnp_rape[i]/deltaRapEven[i];
      //2S
      CS2S_pp_tnp_rap2014[i]= computeRatio( N2S_pp_rap4_2014[i] , Aet_2S_pythia_rap2014[i] ); 
      CS2S_pp_tnp_rap2014e[i]= computeRatioError(  N2S_pp_rap4_2014[i] , Aet_2S_pythia_rap2014[i] ,  N2S_pp_rap4_2014e[i] , Aet_2S_pythia_rap2014e[i] );
      CS2S_pp_tnp_rap2014[i]=CS2S_pp_tnp_rap2014[i]/(deltaRapEven[i] * L_pp);
      CS2S_pp_tnp_rap2014e[i]=CS2S_pp_tnp_rap2014e[i]/(deltaRapEven[i] * L_pp);
      //3S
      CS3S_pp_tnp_rap2014[i]= computeRatio( N3S_pp_rap4_2014[i] , Aet_3S_pythia_rap2014[i] ); 
      CS3S_pp_tnp_rap2014e[i] = computeRatioError( N3S_pp_rap4_2014[i] , Aet_3S_pythia_rap2014[i] , N3S_pp_rap4_2014e[i] , Aet_3S_pythia_rap2014e[i] );
      CS3S_pp_tnp_rap2014[i]=      CS3S_pp_tnp_rap2014[i]/L_pp;
      CS3S_pp_tnp_rap2014e[i]= CS3S_pp_tnp_rap2014e[i]/L_pp;
      CS3S_pp_tnp_rap2014[i]=      CS3S_pp_tnp_rap2014[i]/(deltaRapEven[i]);
      CS3S_pp_tnp_rap2014e[i]= CS3S_pp_tnp_rap2014e[i]/L_pp;
      CS3S_pp_tnp_rap2014e[i]= CS3S_pp_tnp_rap2014e[i]/(deltaRapEven[i]);
      //RAAAA and other ratios!
      RAA_1S_rap[i]= computeRatio( CS1S_aa_rap[i] , CS1S_pp_rap2014[i]);
      RAA_1S_rape[i]= computeRatioError( CS1S_aa_rap[i] , CS1S_pp_rap2014[i],  CS1S_aa_rape[i] , CS1S_pp_rap2014e[i]);
      //tag and probe just after
      RAA_1S_tnp_rap[i]= computeRatio( CS1S_aa_tnp_rap[i] , CS1S_pp_tnp_rap2014[i]);
      RAA_1S_tnp_rape[i]= computeRatioError( CS1S_aa_tnp_rap[i] , CS1S_pp_tnp_rap2014[i],  CS1S_aa_tnp_rape[i] , CS1S_pp_tnp_rap2014e[i]);
      R_A_1S_rap[i]=computeRatio(A_1S_pythia_rap3p5[i],A_1S_pyquen_rap2014[i]);
      R_A_1S_rape[i]=computeRatioError(A_1S_pythia_rap3p5[i],A_1S_pyquen_rap2014[i],A_1S_pythia_rap3p5e[i],A_1S_pyquen_rap2014e[i]);
      R_e_1S_rap[i]=computeRatio(e_1S_pythia_rap3p5[i],e_1S_pyquen_rap2014[i]);
      R_e_1S_rape[i]=computeRatioError(e_1S_pythia_rap3p5[i],e_1S_pyquen_rap2014[i],e_1S_pythia_rap3p5e[i],e_1S_pyquen_rap2014e[i]);
      
    }


cout << "  --- 1S Cross section in pp vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_pp = "<< CS1S_pp_rap2014[j] <<" +/- "<<CS1S_pp_rap2014e[j]<<" b" << endl;
   }

cout << "  --- 1S Cross section in pp vs. pt ---" << endl;
 for(int j =0 ; j<nRapBins_2013 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_pp = "<< CS1S_pp_pt[j] <<" +/- "<<CS1S_pp_pte[j]<<" b" << endl;
   }


cout << "  --- 1S Cross section in PbPb vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_PbPb = "<< CS1S_aa_rap[j] <<" +/- "<<CS1S_aa_rape[j]<<" b" << endl;
   }

cout << "  --- 1S Cross section in PbPb vs. pt ---" << endl;
 for(int j =0 ; j<nRapBins_2013 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_PbPb = "<< CS1S_aa_pt[j] <<" +/- "<<CS1S_aa_pte[j]<<" b" << endl;
   }

cout << "  --- 1S RAA vs. p_{T} ---" << endl;
 for(int j =0 ; j<nPtBins_2013 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_1S_pt[j] <<" +/- "<< RAA_1S_pte[j]<< endl;
   }

cout << "  --- 1S RAA vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_1S_rap[j] <<" +/- "<< RAA_1S_rape[j]<< endl;
   }

 for(int i=0 ; i<nRapBins_2010; i++)
   {
     //    CS2S_pp_rap[i]= computeRatio( N2S_pp_rap3p5[i] , Ae_2S_pythia_rap[i] ); 
     CS2S_pp_rap[i]= computeRatio( N2S_pp_rap4_2014Large[i] , Ae_2S_pythia_rap2010[i] ); 
     CS2S_aa_rap[i]= computeRatio( N2S_aa_rap4_2014Large[i] , Ae_2S_pyquen_rap[i] );
     CS2S_pp_rape[i] = computeRatioError(  N2S_pp_rap4_2014Large[i] , Ae_2S_pythia_rap2010[i] ,  N2S_pp_rap4_2014Largee[i] , Ae_2S_pythia_rap2010e[i] );
     CS2S_aa_rape[i] = computeRatioError(   N2S_aa_rap4_2014Large[i] , Ae_2S_pyquen_rap[i] ,  N2S_aa_rap4_2014Largee[i] , Ae_2S_pyquen_rape[i]);
     CS2S_pp_rap[i]=CS2S_pp_rap[i]/(L_pp*deltaRap2010[i]);
     CS2S_aa_rap[i]=CS2S_aa_rap[i]/(N_MB_corr * T_AA_b * deltaRap2010[i]);
     CS2S_pp_rape[i]=CS2S_pp_rape[i]/(L_pp*deltaRap2010[i]);
     CS2S_aa_rape[i]=CS2S_aa_rape[i]/(N_MB_corr * T_AA_b * deltaRap2010[i]);
     RAA_2S_rap[i]= computeRatio( CS2S_aa_rap[i] , CS2S_pp_rap[i]);
     RAA_2S_rape[i]= computeRatioError( CS2S_aa_rap[i] , CS2S_pp_rap[i],  CS2S_aa_rape[i] , CS2S_pp_rape[i]);

     //keeping it just in case i need to copy to compare with tnp in large bins as well.
     // CS2S_pp_rap[i]= computeRatio( N2S_pp_rap4_2014Large[i] , Ae_2S_pythia_rap2010[i] ); 
     // CS2S_aa_rap[i]= computeRatio( N2S_aa_rap4_2014Large[i] , Ae_2S_pyquen_rap[i] );
     // CS2S_pp_rape[i] = computeRatioError(  N2S_pp_rap4_2014Large[i] , Ae_2S_pythia_rap2010[i] ,  N2S_pp_rap4_2014Largee[i] , Ae_2S_pythia_rap2010e[i] );
     // CS2S_aa_rape[i] = computeRatioError(   N2S_aa_rap4_2014Large[i] , Ae_2S_pyquen_rap[i] ,  N2S_aa_rap4_2014Largee[i] , Ae_2S_pyquen_rape[i]);
     // CS2S_pp_rap[i]=CS2S_pp_rap[i]/(L_pp*deltaRap2010[i]);
     // CS2S_aa_rap[i]=CS2S_aa_rap[i]/(N_MB_corr * T_AA_b * deltaRap2010[i]);
     // CS2S_pp_rape[i]=CS2S_pp_rape[i]/(L_pp*deltaRap2010[i]);
     // CS2S_aa_rape[i]=CS2S_aa_rape[i]/(N_MB_corr * T_AA_b * deltaRap2010[i]);
     // RAA_2S_rap[i]= computeRatio( CS2S_aa_rap[i] , CS2S_pp_rap[i]);
     // RAA_2S_rape[i]= computeRatioError( CS2S_aa_rap[i] , CS2S_pp_rap[i],  CS2S_aa_rape[i] , CS2S_pp_rape[i]);
 }
 
 for(int i = 0 ; i<nPtBins_2010 ; i++)
   {
     CS1S_pp_ptLarge[i]=computeRatio(N1S_pp_pt3p5Large[i],0.231);
     CS1S_pp_pteLarge[i]=computeRatioError(N1S_pp_pt3p5Large[i],0.231,N1S_pp_pt3p5eLarge[i],0.0006);
     CS1S_aa_ptLarge[i]=computeRatio(N1S_aa_pt3p5Large[i],0.210);
     CS1S_aa_pteLarge[i]=computeRatioError(N1S_aa_pt3p5Large[i],0.210,N1S_aa_pt3p5eLarge[i],0.0009);
     //   CS1S_pp_ptLarge[i]=computeRatio(N1S_pp_pt3p5Large[i],N1S_aa_pt3p5Large[i]);
     CS1S_pp_ptLarge[i]=CS1S_pp_ptLarge[i]/(RapBinWidth*L_pp);
     CS1S_aa_ptLarge[i]=CS1S_aa_ptLarge[i]/(RapBinWidth*N_MB_corr*T_AA_b);
     CS1S_pp_pteLarge[i]=CS1S_pp_pteLarge[i]/(RapBinWidth*L_pp);
     CS1S_aa_pteLarge[i]=CS1S_aa_pteLarge[i]/(RapBinWidth*N_MB_corr*T_AA_b);
     //invariant yields.
     CS2S_pp_pt[i]= computeRatio( N2S_pp_pt4_2013Large[i] , Ae_2S_pythia_pt2010[i] ); 
     CS2S_aa_pt[i]= computeRatio( N2S_aa_pt4_2013Large[i] , Ae_2S_pyquen_pt[i] );
     CS2S_pp_pte[i]= computeRatioError( N2S_pp_pt4_2013Large[i] , Ae_2S_pythia_pt2010[i], N2S_pp_pt4_2013Largee[i] , Ae_2S_pythia_pt2010e[i]);
     CS2S_aa_pte[i]= computeRatioError(  N2S_aa_pt4_2013Large[i] , Ae_2S_pyquen_pt[i] ,  N2S_aa_pt4_2013Largee[i] , Ae_2S_pyquen_pte[i] );
     CS2S_pp_pt[i]=CS2S_pp_pt[i]/(RapBinWidth*L_pp);
     CS2S_aa_pt[i]=CS2S_aa_pt[i]/(N_MB_corr * T_AA_b*RapBinWidth);
     CS2S_pp_pte[i]=CS2S_pp_pte[i]/(RapBinWidth*L_pp);
     CS2S_aa_pte[i]=CS2S_aa_pte[i]/(N_MB_corr * T_AA_b * RapBinWidth);
     // cout<<"cs1S (ppPt, ppRap, aaPt, aaRap)"<< CS2S_pp_pt[i] <<", " <<CS2S_pp_rap[i] <<", " <<CS2S_aa_pt[i] <<", " <<CS2S_aa_rap[i] <<". " << endl;
     // cout <<"sigma(2S)_pp vs Pt ="<< endl;
     //invariant yields, corrected by taa and lumi corr.
     RAA_1S_ptLarge[i]=computeRatio(CS1S_aa_ptLarge[i],CS1S_pp_ptLarge[i]); 
     RAA_1S_pteLarge[i]=computeRatioError(CS1S_aa_ptLarge[i],CS1S_pp_ptLarge[i],CS1S_aa_pteLarge[i],CS1S_pp_pteLarge[i]); 
     RAA_2S_pt[i]= computeRatio( CS2S_aa_pt[i] , CS2S_pp_pt[i]);
     RAA_2S_pte[i] =computeRatioError( CS2S_aa_pt[i] , CS2S_pp_pt[i], CS2S_aa_pte[i] , CS2S_pp_pte[i]);
   }

 cout << "  --- 2S Cross section in pp vs. y in large bins---" << endl;
 for(int j =0 ; j<nRapBins_2010 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(2S)_pp = "<< CS2S_pp_rap[j] <<" +/- "<<CS2S_pp_rape[j]<<" b" << endl;
   }

cout << "  --- 2S Cross section in pp vs. pt in large bins---" << endl;
 for(int j =0 ; j<nRapBins_2010 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(2S)_pp = "<< CS2S_pp_pt[j] <<" +/- "<<CS2S_pp_pte[j]<<" b" << endl;
   }


cout << "  --- 2S Cross section in PbPb vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2010 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(2S)_PbPb = "<< CS2S_aa_rap[j] <<" +/- "<<CS2S_aa_rape[j]<<" b" << endl;
   }

cout << "  --- 2S Cross section in PbPb vs. pt ---" << endl;
 for(int j =0 ; j<nRapBins_2010 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(2S)_PbPb = "<< CS2S_aa_pt[j] <<" +/- "<<CS2S_aa_pte[j]<<" b" << endl;
   }

cout << "  --- 1S RAA vs. p_{T} in LARGE BINS---" << endl;
 for(int j =0 ; j<nPtBins_2010 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_1S_ptLarge[j] <<" +/- "<< RAA_1S_pteLarge[j]<< endl;
   }

cout << "  --- 2S RAA vs. p_{T} ---" << endl;
 for(int j =0 ; j<nPtBins_2010 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_2S_pt[j] <<" +/- "<< RAA_2S_pte[j]<< endl;
   }

cout << "  --- 2S RAA vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2010 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_2S_rap[j] <<" +/- "<< RAA_2S_rape[j]<< endl;
   }

 //for the cross sections in pp : shorter bins


cout << "  --- 2S Cross section in pp vs. pt short bins---" << endl;
 for(int j =0 ; j<nPtBins_2013 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(2S)_pp = "<< CS2S_pp_pt2013[j] <<" +/- "<<CS2S_pp_pt2013e[j]<<" b" << endl;
   }

cout << "  --- 2S Cross section in pp vs. y short bins---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(2S)_pp = "<< CS2S_pp_rap2014[j] <<" +/- "<<CS2S_pp_rap2014e[j]<<" b" << endl;
   }

 if(plotCS){
 ////////////////////////////////////////////////////////////////
 /// drawing Pt-binned Data
 ////////////////////////////////////////////////////////////////
 
 TCanvas *cpt = new TCanvas("cpt","cpt"); 
 cpt->cd();
 TPad *ppt1 = new TPad("ppt1","ppt1",0.0,0.0,1.0,1.0);
 ppt1->SetBottomMargin(0.12);
 ppt1->SetTopMargin(0.03);
 ppt1->SetRightMargin(0.03);
 ppt1->SetLeftMargin(0.16);
 ppt1->SetLogy();
 ppt1->Draw();
 ppt1->cd();
 TF1 *f4Pt = new TF1("f4Pt","0.000000001",0,21);
 f4Pt->SetLineWidth(0);
 f4Pt->GetYaxis()->SetTitleOffset(2);
 f4Pt->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");		
 // f4Pt->GetYaxis()->SetTitle("#frac{1}{L_{pp,PbPb}}#frac{N}{#alpha #varepsilon #Deltap_{T}} (b/(GeV/c))");
 f4Pt->GetYaxis()->SetTitle("#sigma(#varUpsilon #rightarrow #mu^{+}#mu^{-1}) (b)");
 f4Pt->GetYaxis()->SetTitleSize(0.028);
 f4Pt->GetYaxis()->SetRangeUser(0.000000000001,.000000001);
 f4Pt->GetXaxis()->CenterTitle(kTRUE);
 f4Pt->Draw();
 /// one pad to draw PbPb yields,
 TGraphErrors *gpt1 = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_pt,pte,CS1S_aa_pte);
 gpt1->SetMarkerColor(8);
 gpt1->SetMarkerStyle(33);
 gpt1->SetMarkerSize(2);
 TGraphErrors *gpt1circle = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_pt,pte,CS1S_aa_pte);
 gpt1circle->SetMarkerStyle(27);
 gpt1circle->SetMarkerSize(2);
 gpt1circle->SetLineColor(kBlack);
 gpt1->Draw("pe");
 gpt1circle->Draw("p");
 // f4Pt->Draw("same");
 gPad->RedrawAxis();
 TGraphErrors *gpt1pp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 gpt1pp->SetMarkerColor(kAzure-9);
 gpt1pp->SetMarkerStyle(21);
 gpt1pp->SetMarkerSize(1.2);
 TGraphErrors *gpt1circlepp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 gpt1circlepp->SetMarkerStyle(25);
 gpt1circlepp->SetMarkerSize(1.22);
 gpt1circlepp->SetLineColor(kBlack);
 gpt1pp->Draw("pe");
 gpt1circlepp->Draw("p");
 //f4Pt->Draw("same");
 gPad->RedrawAxis();

 ///tnp comparison
 if(plotTNP){
   TGraphErrors *gpt1TNP = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_tnp_pt,pte,CS1S_aa_tnp_pte);
   gpt1TNP->SetMarkerColor(kOrange-4);
   gpt1TNP->SetMarkerStyle(33);
   gpt1TNP->SetMarkerSize(2);
   TGraphErrors *gpt1TNPcircle = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_tnp_pt,pte,CS1S_aa_tnp_pte);
   gpt1TNPcircle->SetMarkerStyle(27);
   gpt1TNPcircle->SetMarkerSize(2);
   gpt1TNPcircle->SetLineColor(kBlack);
   gpt1TNP->Draw("pe");
   gpt1TNPcircle->Draw("p");
   // f4Pt->Draw("same");
   gPad->RedrawAxis();
   TGraphErrors *gpt1TNPpp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_tnp_pt,pte,CS1S_pp_tnp_pte);
   gpt1TNPpp->SetMarkerColor(kViolet-9);
   gpt1TNPpp->SetMarkerStyle(21);
   gpt1TNPpp->SetMarkerSize(1.2);
   TGraphErrors *gpt1TNPcirclepp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_tnp_pt,pte,CS1S_pp_tnp_pte);
   gpt1TNPcirclepp->SetMarkerStyle(25);
   gpt1TNPcirclepp->SetMarkerSize(1.22);
   gpt1TNPcirclepp->SetLineColor(kBlack);
   gpt1TNPpp->Draw("pe");
   gpt1TNPcirclepp->Draw("p");
   //f4Pt->Draw("same");
   gPad->RedrawAxis(); 
 }
 TLegend *legend = new TLegend(0.682,0.67,0.85,0.8);
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 legend->AddEntry(gpt1pp,"#varUpsilon(1S), pp ","lp");
 legend->AddEntry(gpt1,"#varUpsilon(1S), PbPb ","lp");
 if(plotTNP){
   legend->AddEntry(gpt1TNPpp,"#varUpsilon(1S), pp TNP corrected","lp");
   legend->AddEntry(gpt1TNP,"#varUpsilon(1S), PbPb TNP corrected","lp");
 } 
 legend->Draw();
 // TLatex *l1CMSpt = new TLatex(10,0.0000000008, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
 TLatex *l1CMSpt = new TLatex(12,0.0000000005, "Work in progress");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.04);
 l1CMSpt->Draw();
 TLatex *lyL= new TLatex(2,0.000000000011,"L_{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 lyL->SetTextFont(42);
 lyL->SetTextSize(0.029);
 lyL->DrawLatex(2,0.000000000007,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->DrawLatex(2,0.000000000004,"corrected cross-sections");
 lyL->Draw();

 //Unfolding unfolding!
 // TCanvas *cUnfold = new TCanvas("cUnfold","cUnfold"); 
 // cUnfold->cd();
 // TPad *pUnfold1 = new TPad("pUnfold1","pUnfold1",0.0,0.0,1.0,1.0);
 // pUnfold1->SetBottomMargin(0.12);
 // pUnfold1->SetTopMargin(0.03);
 // pUnfold1->SetRightMargin(0.03);
 // pUnfold1->SetLeftMargin(0.12);
 // pUnfold1->SetLogy();
 // pUnfold1->Draw();
 // pUnfold1->cd();
 // //f4Unfold->GetYaxis()->SetRangeUser(0.01,.09);
 // TGraphErrors *gpt1unfoldpp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 // gpt1unfoldpp->SetMarkerStyle(25);
 // gpt1unfoldpp->SetMarkerSize(1.22);
 // gpt1unfoldpp->SetLineColor(kBlack);
 // gpt1unfoldpp->Draw("p");
 // gPad->RedrawAxis();

 //////////
 //2S+1S pp
 TCanvas *cptpp = new TCanvas("cptpp","cptpp"); 
 cptpp->cd();
 TPad *ppt1pp = new TPad("ppt1pp","ppt1pp",0.0,0.0,1.0,1.0);
 ppt1pp->SetBottomMargin(0.12);
 ppt1pp->SetTopMargin(0.03);
 ppt1pp->SetRightMargin(0.03);
 ppt1pp->SetLeftMargin(0.16);
 ppt1pp->SetLogy();
 ppt1pp->Draw();
 ppt1pp->cd();
 TF1 *f4Pt = new TF1("f4Pt","0.000000001",0,21);
 f4Pt->SetLineWidth(0);
 f4Pt->GetYaxis()->SetTitleOffset(2);
 f4Pt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} (GeV/c)");		
 f4Pt->GetYaxis()->SetTitle("#frac{1}{<L_{pp}>}#frac{N}{#alpha #varepsilon #Deltap_{T}} (b/(GeV/c))");

 f4Pt->GetYaxis()->SetTitleSize(0.028);
 f4Pt->GetYaxis()->SetRangeUser(0.000000000001,0.000000001);
 //f4Pt->GetYaxis()->SetRangeUser(0.01,.09);
 f4Pt->GetXaxis()->CenterTitle(kTRUE);
 f4Pt->Draw();
 TGraphErrors *g2pt = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013,pte,CS2S_pp_pt2013e);
 g2pt->SetMarkerColor(kAzure-8);
 g2pt->SetMarkerStyle(21);
 g2pt->SetMarkerSize(1.2);

 TGraphErrors *g3pt = new TGraphErrors(nPtBins_2013,pt,CS3S_pp_pt2013,pte,CS3S_pp_pt2013e);
 g3pt->SetMarkerColor(kAzure-7);
 g3pt->SetMarkerStyle(21);
 g3pt->SetMarkerSize(1.2);


 // TGraphErrors *g2Ssyst = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013,pte,CS2S_pp_pt2013e);
 // g2Ssyst->SetLineColor(kAzure-9);
 // g2Ssyst->SetFillStyle(0);
 // // g2Ssyst->SetLineWidth(18);
 // g2Ssyst->SetMarkerSize(0);
 // g2Ssyst->Draw("2");


 TGraphErrors *g2circle = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013,pte,CS2S_pp_pt2013e);
 g2circle->SetMarkerStyle(25);
 g2circle->SetMarkerSize(1.3);
 g2circle->SetLineColor(kBlack);
 g2pt->Draw("pe");
 g2circle->Draw("p");

 TGraphErrors *g3circle = new TGraphErrors(nPtBins_2013,pt,CS3S_pp_pt2013,pte,CS3S_pp_pt2013e);
 g3circle->SetMarkerStyle(25);
 g3circle->SetMarkerSize(1.3);
 g3circle->SetLineColor(kBlack);
 g3pt->Draw("pe");
 g3circle->Draw("p");

TGraphErrors *gpt1pp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 gpt1pp->SetMarkerColor(kAzure-9);
 gpt1pp->SetMarkerStyle(21);
 gpt1pp->SetMarkerSize(1.2);
 TGraphErrors *gpt1circlepp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 gpt1circlepp->SetMarkerStyle(25);
 gpt1circlepp->SetMarkerSize(1.22);
 gpt1circlepp->SetLineColor(kBlack);
 gpt1pp->Draw("pe");
 gpt1circlepp->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();
 f4Pt->Draw("same");
 gPad->RedrawAxis();
 TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 legend->AddEntry(gpt1pp,"#varUpsilon(1S), pp ","lp");
 legend->AddEntry(g2pt,"#varUpsilon(2S), pp ","lp");
 legend->AddEntry(g3pt,"#varUpsilon(3S), pp ","lp");
 legend->Draw();
 TLatex *l1CMSpt = new TLatex(2,0.0000000008, "CMS Internal #sqrt{s} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.036);
 l1CMSpt->Draw();
 
 TLatex *lyL= new TLatex(2,0.0000000006,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->SetTextSize(0.029);
 lyL->Draw();
 }


 if(plotFiducial){
 //////////
 //2S+1S pp FIDUCIAL
 TCanvas *cptppFiducial = new TCanvas("cptppFiducial","cptppFiducial"); 
 cptppFiducial->cd();
 TPad *ppt1ppFiducial = new TPad("ppt1ppFiducial","ppt1ppFiducial",0.0,0.0,1.0,1.0);
 ppt1ppFiducial->SetBottomMargin(0.12);
 ppt1ppFiducial->SetTopMargin(0.03);
 ppt1ppFiducial->SetRightMargin(0.03);
 ppt1ppFiducial->SetLeftMargin(0.16);
 ppt1ppFiducial->SetLogy();
 ppt1ppFiducial->Draw();
 ppt1ppFiducial->cd();
 TF1 *f4Pt = new TF1("f4Pt","0.000000001",0,21);
 f4Pt->SetLineWidth(0);
 f4Pt->GetYaxis()->SetTitleOffset(2);
 f4Pt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} (GeV/c)");		
 f4Pt->GetYaxis()->SetTitle("#frac{1}{L_{pp}}#frac{N}{#varepsilon #Deltap_{T}} (b/(GeV/c))");

 f4Pt->GetYaxis()->SetTitleSize(0.028);
 //f4Pt->GetYaxis()->SetRangeUser(0.01,.09);
 f4Pt->GetXaxis()->CenterTitle(kTRUE);
 f4Pt->Draw();
 TGraphErrors *g2ptFiducial = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013Fiducial,pte,CS2S_pp_pt2013Fiduciale);
 g2ptFiducial->SetMarkerColor(kAzure-7);
 g2ptFiducial->SetMarkerStyle(21);
 g2ptFiducial->SetMarkerSize(1.2);


 // TGraphErrors *g2Ssyst = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013,pte,CS2S_pp_pt2013e);
 // g2Ssyst->SetLineColor(kAzure-9);
 // g2Ssyst->SetFillStyle(0);
 // // g2Ssyst->SetLineWidth(18);
 // g2Ssyst->SetMarkerSize(0);
 // g2Ssyst->Draw("2");


 TGraphErrors *g2circleF = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013Fiducial,pte,CS2S_pp_pt2013Fiduciale);
 g2circleF->SetMarkerStyle(25);
 g2circleF->SetMarkerSize(1.2);
 g2circleF->SetLineColor(kBlack);
 g2ptFiducial->Draw("pe");
 g2circleF->Draw("p");

TGraphErrors *gpt1ppF = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_ptFiducial,pte,CS1S_pp_ptFiduciale);
 gpt1ppF->SetMarkerColor(kAzure-9);
 gpt1ppF->SetMarkerStyle(21);
 gpt1ppF->SetMarkerSize(1.2);
 TGraphErrors *gpt1circleFpp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_ptFiducial,pte,CS1S_pp_ptFiduciale);
 gpt1circleFpp->SetMarkerStyle(25);
 gpt1circleFpp->SetMarkerSize(1.22);
 gpt1circleFpp->SetLineColor(kBlack);
 gpt1ppF->Draw("pe");
 gpt1circleFpp->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();
 f4Pt->Draw("same");
 gPad->RedrawAxis();
 TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 legend->AddEntry(gpt1ppF,"#varUpsilon(1S), pp ","lp");
 legend->AddEntry(g2ptFiducial,"#varUpsilon(2S), pp ","lp");
 legend->Draw();
 TLatex *l1CMSpt = new TLatex(2,0.0000000008, "CMS Internal #sqrt{s} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.032);
 l1CMSpt->Draw();
 
 TLatex *lyL= new TLatex(2,0.0000000005,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->SetTextSize(0.04);
lyL->DrawLatex(2,0.000000002,"Fiducial");
 lyL->Draw();


 //comparison pp,pbpb fiducial



TCanvas *cptFiducial = new TCanvas("cptFiducial","cptFiducial"); 
 cptFiducial->cd();
 TPad *ppt1ppFiducial = new TPad("ppt1ppFiducial","ppt1ppFiducial",0.0,0.0,1.0,1.0);
 ppt1ppFiducial->SetBottomMargin(0.12);
 ppt1ppFiducial->SetTopMargin(0.03);
 ppt1ppFiducial->SetRightMargin(0.03);
 ppt1ppFiducial->SetLeftMargin(0.16);
 ppt1ppFiducial->SetLogy();
 ppt1ppFiducial->Draw();
 ppt1ppFiducial->cd();
 TF1 *f4Pt = new TF1("f4Pt","0.000000001",0,21);
 f4Pt->SetLineWidth(0);
 f4Pt->GetYaxis()->SetTitleOffset(2);
 f4Pt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} (GeV/c)");		
 f4Pt->GetYaxis()->SetTitle("#frac{1}{L_{pp}}#frac{N}{#varepsilon #Deltap_{T}} (b/(GeV/c))");

 f4Pt->GetYaxis()->SetTitleSize(0.028);
 //f4Pt->GetYaxis()->SetRangeUser(0.01,.09);
 f4Pt->GetXaxis()->CenterTitle(kTRUE);
 f4Pt->Draw();

TGraphErrors *gpt1F = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_ptFiducial,pte,CS1S_aa_ptFiduciale);
 gpt1F->SetMarkerColor(8);
 gpt1F->SetMarkerStyle(33);
 gpt1F->SetMarkerSize(2);



 TGraphErrors *gpt1circle = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_ptFiducial,pte,CS1S_aa_ptFiduciale);
 gpt1circle->SetMarkerStyle(27);
 gpt1circle->SetMarkerSize(2);
 gpt1circle->SetLineColor(kBlack);
 gpt1F->Draw("pe");
 gpt1circle->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();
    

TGraphErrors *gpt1ppF = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_ptFiducial,pte,CS1S_pp_ptFiduciale);
 gpt1ppF->SetMarkerColor(kAzure-9);
 gpt1ppF->SetMarkerStyle(21);
 gpt1ppF->SetMarkerSize(1.2);
 TGraphErrors *gpt1circleFpp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_ptFiducial,pte,CS1S_pp_ptFiduciale);
 gpt1circleFpp->SetMarkerStyle(25);
 gpt1circleFpp->SetMarkerSize(1.22);
 gpt1circleFpp->SetLineColor(kBlack);
 gpt1ppF->Draw("pe");
 gpt1circleFpp->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();
 f4Pt->Draw("same");
 gPad->RedrawAxis();
 TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 legend->AddEntry(gpt1ppF,"#varUpsilon(1S), pp ","lp");
 legend->AddEntry(gpt1F,"#varUpsilon(1S), PbPb ","lp");
 legend->Draw();
 TLatex *l1CMSpt = new TLatex(2,0.0000000008, "CMS Internal #sqrt{s} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.032);
 l1CMSpt->Draw();
 
 TLatex *lyL= new TLatex(2,0.000000008,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->SetTextSize(0.04);
lyL->DrawLatex(2,0.000000002,"Fiducial");
 lyL->Draw();

}
 ///raa vs pt

 if(plotRAA){
 TCanvas *cRaapt = new TCanvas("cRaapt","cRaapt"); 
 cRaapt->cd();
 TPad *ppt2 = new TPad("ppt2","ppt2",0.0,0.0,1.0,1.0);
 ppt2->SetBottomMargin(0.12);
 ppt2->SetTopMargin(0.03);
 ppt2->SetRightMargin(0.03);
 ppt2->SetLeftMargin(0.16);
 ppt2->Draw();
 ppt2->cd();
 //one pad to draw RaaPt!
 TF1 *f4RaaPt = new TF1("f4RaaPt","1",0,20.5);
 f4RaaPt->SetLineWidth(0);
 f4RaaPt->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
 f4RaaPt->GetYaxis()->SetTitle("R_{AA}");
 f4RaaPt->GetYaxis()->SetTitleOffset(1.8);
 f4RaaPt->GetYaxis()->SetTitleSize(0.028);
 f4RaaPt->GetYaxis()->SetRangeUser(0.,2);
 f4RaaPt->GetXaxis()->CenterTitle(kTRUE);
 f4RaaPt->Draw();
 TGraphErrors *gRaaPt1 = new TGraphErrors(nPtBins_2013,pt,RAA_1S_pt,pte,RAA_1S_pte);
 gRaaPt1->SetMarkerColor(kBlue-4);
 gRaaPt1->SetMarkerStyle(33);
 gRaaPt1->SetMarkerSize(2);
 TGraphErrors *gRaaPt1circle = new TGraphErrors(nPtBins_2013,pt,RAA_1S_pt,pte,RAA_1S_pte);
 gRaaPt1circle->SetMarkerStyle(27);
 gRaaPt1circle->SetMarkerSize(2);
 gRaaPt1circle->SetLineColor(kBlack);
 gRaaPt1->Draw("pe");
 gRaaPt1circle->Draw("p");
 f4RaaPt->Draw("same");
 if(plotTNP){
   TGraphErrors *gRaaPt1TNP = new TGraphErrors(nPtBins_2013,pt,RAA_1S_tnp_pt,pte,RAA_1S_tnp_pte);
   gRaaPt1TNP->SetMarkerColor(kOrange+1);
   gRaaPt1TNP->SetMarkerStyle(33);
   gRaaPt1TNP->SetMarkerSize(2);
   TGraphErrors *gRaaPt1TNPcircle = new TGraphErrors(nPtBins_2013,pt,RAA_1S_tnp_pt,pte,RAA_1S_tnp_pte);
   gRaaPt1TNPcircle->SetMarkerStyle(27);
   gRaaPt1TNPcircle->SetMarkerSize(2);
   gRaaPt1TNPcircle->SetLineColor(kBlack);
   gRaaPt1TNP->Draw("pe");
   gRaaPt1TNPcircle->Draw("p");
   f4RaaPt->Draw("same");
 }
 TLatex *l1CMSpt = new TLatex(12,0.2, "Work in progress");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.04);
 l1CMSpt->Draw();
 TLatex *lyLPT= new TLatex(1,1.8,"2011 L_{int}^{PbPb} = 150 #mub^{-1}; 2013 L_{int}^{pp} = 5.4 pb^{-1};");
 lyLPT->SetTextFont(42);
 lyLPT->SetTextSize(0.027);
 lyLPT->Draw();
 lyLPT->DrawLatex(1,1.6,"2010 L_{int}^{PbPb} = 7.28 #mub^{-1}; L_{int}^{pp} = 225 nb^{-1};");
 ppt2->Update();
 //2010stuff
 TGraphErrors *gpt2010 = new TGraphErrors(nPtBins_2010,pt_2010,raaPt2010,pte_2010,raaPt2010e);
 gpt2010->SetMarkerColor(kTeal+3);
 gpt2010->SetMarkerStyle(33);
 gpt2010->SetMarkerSize(2);
 TGraphErrors *gpt2010s = new TGraphErrors(nPtBins_2010,pt_2010,raaPt2010,centnoErr,raaPt2010s);
 gpt2010s->SetLineColor(8);
 gpt2010s->SetLineWidth(18);
 gpt2010s->SetMarkerSize(0);
 gpt2010s->Draw("e");
 gpt2010->Draw("pe");
 TGraphErrors *gpt2010circle = new TGraphErrors(nPtBins_2010,pt_2010,raaPt2010,pte_2010,raaPt2010e);
 gpt2010circle->SetMarkerStyle(27);
 gpt2010circle->SetMarkerSize(2);
 gpt2010circle->SetLineColor(kBlack);
 gpt2010circle->Draw("p");
 f4RaaPt->Draw("same");
 gPad->RedrawAxis();
 TLegend *legend = new TLegend(0.2,0.65,0.4,0.75);
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 legend->AddEntry(gpt2010,"#varUpsilon(1S) JHEP 05 (2012) 063","lp");
 legend->AddEntry(gRaaPt1,"#varUpsilon(1S) in progress","lp");
  if(plotTNP)
    {
      legend->AddEntry(gRaaPt1TNP,"#varUpsilon(1S) in progress","lp");
     }
  legend->Draw();
  gPad->RedrawAxis();
 }

 if(plot2010){
   TCanvas *cRaapt = new TCanvas("cRaapt","cRaapt"); 
   cRaapt->cd();
   TPad *ppt2 = new TPad("ppt2","ppt2",0.0,0.0,1.0,1.0);
   ppt2->SetBottomMargin(0.12);
   ppt2->SetTopMargin(0.03);
   ppt2->SetRightMargin(0.03);
   ppt2->SetLeftMargin(0.16);
   ppt2->Draw();
   ppt2->cd();
   //one pad to draw RaaPt!
   TF1 *f4RaaPt = new TF1("f4RaaPt","1",0,20.5);
   f4RaaPt->SetLineWidth(0);
   f4RaaPt->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
   f4RaaPt->GetYaxis()->SetTitle("R_{AA}");
   f4RaaPt->GetYaxis()->SetTitleOffset(1.8);
   f4RaaPt->GetYaxis()->SetTitleSize(0.028);
   f4RaaPt->GetYaxis()->SetRangeUser(0.,2);
   f4RaaPt->GetXaxis()->CenterTitle(kTRUE);
   f4RaaPt->Draw();
   TLatex *l1CMSpt = new TLatex(12,0.2, "Work in progress");
   l1CMSpt->SetTextFont(42);
   l1CMSpt->SetTextSize(0.04);
   l1CMSpt->Draw();
   TLatex *lyLPT= new TLatex(1,1.8,"2011 L_{int}^{PbPb} = 150 #mub^{-1}; 2013 L_{int}^{pp} = 5.4 pb^{-1};");
   lyLPT->SetTextFont(42);
   lyLPT->SetTextSize(0.027);
   lyLPT->Draw();
   lyLPT->DrawLatex(1,1.6,"2010 L_{int}^{PbPb} = 7.28 #mub^{-1}; L_{int}^{pp} = 225 nb^{-1};");
   ppt2->Update();
   //2010stuff
   TGraphErrors *gpt2010 = new TGraphErrors(nPtBins_2010,pt_2010,raaPt2010,pte_2010,raaPt2010e);
   gpt2010->SetMarkerColor(kTeal+3);
   gpt2010->SetMarkerStyle(33);
   gpt2010->SetMarkerSize(2);
   TGraphErrors *gpt2010s = new TGraphErrors(nPtBins_2010,pt_2010,raaPt2010,centnoErr,raaPt2010s);
   gpt2010s->SetLineColor(8);
   gpt2010s->SetLineWidth(18);
   gpt2010s->SetMarkerSize(0);
   gpt2010s->Draw("e");
   gpt2010->Draw("pe");
   TGraphErrors *gpt2010circle = new TGraphErrors(nPtBins_2010,pt_2010,raaPt2010,pte_2010,raaPt2010e);
   gpt2010circle->SetMarkerStyle(27);
   gpt2010circle->SetMarkerSize(2);
   gpt2010circle->SetLineColor(kBlack);
   gpt2010circle->Draw("p");
   f4RaaPt->Draw("same");

   TGraphErrors *gRaaPtLarge1 = new TGraphErrors(nPtBins_2010,pt_2010,RAA_1S_ptLarge,pte_2010,RAA_1S_pteLarge);
   gRaaPtLarge1->SetMarkerColor(kBlue-4);
   gRaaPtLarge1->SetMarkerStyle(33);
   gRaaPtLarge1->SetMarkerSize(2);
   TGraphErrors *gRaaPtLarge1circle = new TGraphErrors(nPtBins_2010,pt_2010,RAA_1S_ptLarge,pte_2010,RAA_1S_pteLarge);
   gRaaPtLarge1circle->SetMarkerStyle(27);
   gRaaPtLarge1circle->SetMarkerSize(2);
   gRaaPtLarge1circle->SetLineColor(kBlack);
   gRaaPtLarge1->Draw("pe");
   gRaaPtLarge1circle->Draw("p");
   f4RaaPt->Draw("same");
   gPad->RedrawAxis();
   TLegend *legend = new TLegend(0.2,0.65,0.4,0.75);
   legend->SetTextSize(0.029);
   legend->SetFillStyle(0);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->SetTextFont(42);
   legend->AddEntry(gpt2010,"#varUpsilon(1S) JHEP 05 (2012) 063","lp");
   legend->AddEntry(gRaaPtLarge1,"#varUpsilon(1S) in progress","lp");
   legend->Draw();
   gPad->RedrawAxis();
 }

 ////////////////////////////////////////////////////////////////
 /// drawing Rap-binned Data
 ////////////////////////////////////////////////////////////////
 if(plotCS){
 TCanvas *crap = new TCanvas("crap","crap"); 
 crap->cd();
 TPad *prap1 = new TPad("prap1","prap1",0.0,0.0,1.0,1.0);
 prap1->SetBottomMargin(0.12);
 prap1->SetTopMargin(0.08);
 prap1->SetRightMargin(0.03);
 prap1->SetLeftMargin(0.16);
 prap1->SetLogy();
 prap1->Draw();
 prap1->cd();
 TF1 *f4Rap = new TF1("f4Rap","0.0000001",0,2.45);
 f4Rap->SetLineWidth(0);

 f4Rap->GetYaxis()->SetTitleOffset(2);
 f4Rap->GetYaxis()->SetTitleSize(0.028);
 f4Rap->GetYaxis()->SetRangeUser(0.000000001,0.00000001);
 f4Rap->GetXaxis()->SetTitle("|y^{#Upsilon}| ");
 f4Rap->GetYaxis()->SetTitle("#frac{1}{L_{pp,PbPb}}#frac{N}{#alpha #varepsilon #Deltay} (b)");
 f4Rap->GetXaxis()->CenterTitle(kTRUE);
 f4Rap->Draw();
 //pbpb data
 TGraphErrors *grap1 = new TGraphErrors(nRapBins_2014,rap2014,CS1S_aa_rap,rap2014e,CS1S_aa_rape);
 grap1->SetMarkerColor(8);
 grap1->SetMarkerStyle(33);
 grap1->SetMarkerSize(2);
 TGraphErrors *grap1circle = new TGraphErrors(nRapBins_2014,rap2014,CS1S_aa_rap,rap2014e,CS1S_aa_rape);
 grap1circle->SetMarkerStyle(27);
 grap1circle->SetMarkerSize(2);
 grap1circle->SetLineColor(kBlack);
 grap1->Draw("pe");
 grap1circle->Draw("p");
 f4Rap->Draw("same");
 gPad->RedrawAxis();
 TGraphErrors *grap1pp = new TGraphErrors(nRapBins_2014,rap2014,CS1S_pp_rap2014,rap2014e,CS1S_pp_rap2014e);
 grap1pp->SetMarkerColor(kAzure-9);
 grap1pp->SetMarkerStyle(21);
 grap1pp->SetMarkerSize(1.2);
 TGraphErrors *grap1circlepp = new TGraphErrors(nRapBins_2014,rap2014,CS1S_pp_rap2014,rap2014e,CS1S_pp_rap2014e);
 grap1circlepp->SetMarkerStyle(25);
 grap1circlepp->SetMarkerSize(1.22);
 grap1circlepp->SetLineColor(kBlack);
 grap1pp->Draw("pe");
 grap1circlepp->Draw("p");
 f4Rap->Draw("same");

 if(plotTNP){
   TGraphErrors *grap1TNP = new TGraphErrors(nRapBins_2014,rap2014,CS1S_aa_tnp_rap,rap2014e,CS1S_aa_tnp_rape);
   grap1TNP->SetMarkerColor(kOrange+1);
   grap1TNP->SetMarkerStyle(33);
   grap1TNP->SetMarkerSize(2);
   TGraphErrors *grap1TNPcircle = new TGraphErrors(nRapBins_2014,rap2014,CS1S_aa_tnp_rap,rap2014e,CS1S_aa_tnp_rape);
   grap1TNPcircle->SetMarkerStyle(27);
   grap1TNPcircle->SetMarkerSize(2);
   grap1TNPcircle->SetLineColor(kBlack);
   grap1TNP->Draw("pe");
   grap1TNPcircle->Draw("p");
   f4Rap->Draw("same");
   gPad->RedrawAxis();
   TGraphErrors *grap1TNPpp = new TGraphErrors(nRapBins_2014,rap2014,CS1S_pp_tnp_rap2014,rap2014e,CS1S_pp_tnp_rap2014e);
   grap1TNPpp->SetMarkerColor(kViolet-9);
   grap1TNPpp->SetMarkerStyle(21);
   grap1TNPpp->SetMarkerSize(1.2);
   TGraphErrors *grap1TNPcirclepp = new TGraphErrors(nRapBins_2014,rap2014,CS1S_pp_tnp_rap2014,rap2014e,CS1S_pp_tnp_rap2014e);
   grap1TNPcirclepp->SetMarkerStyle(25);
   grap1TNPcirclepp->SetMarkerSize(1.22);
   grap1TNPcirclepp->SetLineColor(kBlack);
   grap1TNPpp->Draw("pe");
   grap1TNPcirclepp->Draw("p");
   f4Rap->Draw("same");
 }
 TLatex *l1CMSpt = new TLatex(0.15,7e-9, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.032);
 l1CMSpt->Draw();

 
 TLatex *lyL= new TLatex(0.15,3e-9,"L_{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 
 lyL->SetTextSize(0.029);
 lyL->DrawLatex(0.15,2e-9,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->Draw();

TLegend *legendB = new TLegend(0.7,0.84,0.88,0.73);
legendB->SetTextSize(0.029);
legendB->SetFillStyle(0);
legendB->SetFillColor(0);
legendB->SetBorderSize(0);
legendB->SetTextFont(42);
legendB->AddEntry(grap1,"#varUpsilon(1S) PbPb ","lp");
legendB->AddEntry(grap1pp,"#varUpsilon(1S) pp ","lp");
 if(plotTNP){
   legendB->AddEntry(grap1TNP,"#varUpsilon(1S) PbPb TnP corrected","lp");
   legendB->AddEntry(grap1TNPpp,"#varUpsilon(1S) pp TnP corrected","lp");
 }
 legendB->Draw();

 gPad->RedrawAxis();
 //////////////////////////////////////////////
TCanvas *crappp = new TCanvas("crappp","crappp"); 
 crappp->cd();
 TPad *prap1 = new TPad("prap1","prap1",0.0,0.0,1.0,1.0);
 prap1->SetBottomMargin(0.12);
 prap1->SetTopMargin(0.03);
 prap1->SetRightMargin(0.03);
 prap1->SetLeftMargin(0.16);
 prap1->SetLogy();
 prap1->Draw();
 prap1->cd();
 TF1 *f4Rap = new TF1("f4Rap","0.0000001",0,2.45);
 f4Rap->SetLineWidth(0);
 f4Rap->GetYaxis()->SetTitleOffset(2);
 f4Rap->GetYaxis()->SetTitleSize(0.028);
 f4Rap->GetYaxis()->SetRangeUser(0.00000000001,0.00000001);
 f4Rap->GetXaxis()->SetTitle("|y^{#Upsilon}| ");
 f4Rap->GetYaxis()->SetTitle("#frac{1}{L_{pp,PbPb}}#frac{N}{#alpha #varepsilon #Deltay} (b)");
 f4Rap->GetXaxis()->CenterTitle(kTRUE);
 f4Rap->Draw();
 //pbpb data
TGraphErrors *grap2pp = new TGraphErrors(nRapBins_2014,rap2014,CS2S_pp_rap2014,rap2014e,CS2S_pp_rap2014e);
 grap2pp->SetMarkerColor(kAzure-8);
 grap2pp->SetMarkerStyle(21);
 grap2pp->SetMarkerSize(1.2);
TGraphErrors *grap3pp = new TGraphErrors(nRapBins_2014,rap2014,CS3S_pp_rap2014,rap2014e,CS3S_pp_rap2014e);
 grap3pp->SetMarkerColor(kAzure-7);
 grap3pp->SetMarkerStyle(21);
 grap3pp->SetMarkerSize(1.2);


 // TGraphErrors *g2Ssyst = new TGraphErrors(nPtBins_2013,pt,CS2S_pp_pt2013,pte,CS2S_pp_pt2013e);
 // g2Ssyst->SetLineColor(kAzure-9);
 // g2Ssyst->SetFillStyle(0);
 // // g2Ssyst->SetLineWidth(18);
 // g2Ssyst->SetMarkerSize(0);
 // g2Ssyst->Draw("2");


 TGraphErrors *g2rapcircle = new TGraphErrors(nRapBins_2014,rap2014,CS2S_pp_rap2014,rap2014e,CS2S_pp_rap2014e);
 g2rapcircle->SetMarkerStyle(25);
 g2rapcircle->SetMarkerSize(1.3);
 g2rapcircle->SetLineColor(kBlack);
 grap2pp->Draw("pe");
 g2rapcircle->Draw("p");
 TGraphErrors *g3rapcircle = new TGraphErrors(nRapBins_2014,rap2014,CS3S_pp_rap2014,rap2014e,CS3S_pp_rap2014e);
 g3rapcircle->SetMarkerStyle(25);
 g3rapcircle->SetMarkerSize(1.3);
 g3rapcircle->SetLineColor(kBlack);
 grap3pp->Draw("pe");
 g3rapcircle->Draw("p");
 TGraphErrors *grap1pp = new TGraphErrors(nRapBins_2014,rap2014,CS1S_pp_rap2014,rap2014e,CS1S_pp_rap2014e);
 grap1pp->SetMarkerColor(kAzure-9);
 grap1pp->SetMarkerStyle(21);
 grap1pp->SetMarkerSize(1.2);
 TGraphErrors *grap1circlepp = new TGraphErrors(nRapBins_2014,rap2014,CS1S_pp_rap2014,rap2014e,CS1S_pp_rap2014e);
 grap1circlepp->SetMarkerStyle(25);
 grap1circlepp->SetMarkerSize(1.22);
 grap1circlepp->SetLineColor(kBlack);
 grap1pp->Draw("pe");
 grap1circlepp->Draw("p");
 f4Rap->Draw("same");


 TLatex *l1CMSpt = new TLatex(0.15,7e-9, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.032);
 l1CMSpt->Draw();

 
 TLatex *lyL= new TLatex(0.15,3e-9,"L_{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 
 lyL->SetTextSize(0.029);
 lyL->DrawLatex(0.15,2e-9,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->Draw();

TLegend *legendB = new TLegend(0.7,0.84,0.88,0.73);
legendB->SetTextSize(0.029);
legendB->SetFillStyle(0);
legendB->SetFillColor(0);
legendB->SetBorderSize(0);
legendB->SetTextFont(42);
legendB->AddEntry(grap1pp,"#varUpsilon(1S), pp ","lp");
legendB->AddEntry(grap2pp,"#varUpsilon(2S), pp ","lp");
 legendB->Draw();

 gPad->RedrawAxis();

 //////////////////////////////////////////////////////////////////
 }

 if(plotRAA){
 TCanvas *cRaarap = new TCanvas("cRaarap","cRaarap"); 
 cRaarap->cd();
 TPad *prap2 = new TPad("prap2","prap2",0.0,0.0,1.0,1.0);
 prap2->SetBottomMargin(0.12);
 prap2->SetTopMargin(0.03);
 prap2->SetRightMargin(0.03);
 prap2->SetLeftMargin(0.16);
 prap2->Draw();
 prap2->cd();
 //one pad to draw RaaRap!
 TF1 *f4RaaRap = new TF1("f4RaaRap","1",0,2.45);
 f4RaaRap->SetLineWidth(0);
 f4RaaRap->GetXaxis()->SetTitle("|y^{#Upsilon}|");
 f4RaaRap->GetYaxis()->SetTitle("R_{AA}");
 f4RaaRap->GetYaxis()->SetTitleOffset(1.8);
 f4RaaRap->GetYaxis()->SetTitleSize(0.028);
 f4RaaRap->GetYaxis()->SetRangeUser(0.,2);
 f4RaaRap->GetXaxis()->CenterTitle(kTRUE);
 f4RaaRap->Draw();
 TGraphErrors *gRaaRap1 = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_rap,rap2014e,RAA_1S_rape);
 gRaaRap1->SetMarkerColor(8);
 gRaaRap1->SetMarkerStyle(33);
 gRaaRap1->SetMarkerSize(2);
 TGraphErrors *gRaaRap1circle = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_rap,rap2014e,RAA_1S_rape);
 gRaaRap1circle->SetMarkerStyle(27);
 gRaaRap1circle->SetMarkerSize(2);
 gRaaRap1circle->SetLineColor(kBlack);
 gRaaRap1->Draw("pe");
 gRaaRap1circle->Draw("p");
 f4RaaRap->Draw("same");
 gPad->RedrawAxis();
 if(plotTNP)
   {
     TGraphErrors *gRaaRap1TNP = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_tnp_rap,rap2014e,RAA_1S_tnp_rape);
     gRaaRap1TNP->SetMarkerColor(kOrange+1);
     gRaaRap1TNP->SetMarkerStyle(33);
     gRaaRap1TNP->SetMarkerSize(2);
     TGraphErrors *gRaaRap1TNPcircle = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_tnp_rap,rap2014e,RAA_1S_tnp_rape);
     gRaaRap1TNPcircle->SetMarkerStyle(27);
     gRaaRap1TNPcircle->SetMarkerSize(2);
     gRaaRap1TNPcircle->SetLineColor(kBlack);
     gRaaRap1TNP->Draw("pe");
     gRaaRap1TNPcircle->Draw("p");
     f4RaaRap->Draw("same");
     gPad->RedrawAxis();
   }
 TLatex *l1CMSrap = new TLatex(1.2,0.2, "Work in progress");
 l1CMSrap->SetTextFont(42);
 l1CMSrap->SetTextSize(0.04);
 l1CMSrap->Draw();
 TLatex *lyLRAP= new TLatex(0.2,1.8,"2011 L_{int}^{PbPb} = 150 #mub^{-1}; 2013 L_{int}^{pp} = 5.4 pb^{-1};");
 lyLRAP->SetTextFont(42);
 lyLRAP->SetTextSize(0.027);
 lyLRAP->Draw();
lyLRAP->DrawLatex(0.2,1.6,"2010 L_{int}^{PbPb} = 7.28 #mub^{-1}; L_{int}^{pp} = 225 nb^{-1};");
 prap2->Update();
 TGraphErrors *grap2010 = new TGraphErrors(nRapBins_2010,rap2010,raaRap2010,rap2010e,raaRap2010e);
 grap2010->SetMarkerColor(kTeal+3);
 grap2010->SetMarkerStyle(33);
 grap2010->SetMarkerSize(2);
  
  TGraphErrors *grap2010s = new TGraphErrors(nRapBins_2010,rap2010,raaRap2010,centnoErr,raaRap2010s);
  grap2010s->SetLineColor(8);
  grap2010s->SetLineWidth(18);
  grap2010s->SetMarkerSize(0);
  grap2010s->Draw("e");
  
  grap2010->Draw("pe");
  
  TGraphErrors *grap2010circle = new TGraphErrors(nRapBins_2010,rap2010,raaRap2010,rap2010e,raaRap2010e);
  grap2010circle->SetMarkerStyle(27);
  grap2010circle->SetMarkerSize(2);
  grap2010circle->SetLineColor(kBlack);
  grap2010circle->Draw("p");
  f4RaaRap->Draw("same");
  gPad->RedrawAxis(); 
  TLegend *legend = new TLegend(0.2,0.65,0.4,0.75);
  legend->SetTextSize(0.029);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->AddEntry(grap2010,"#varUpsilon(1S) JHEP 05 (2012) 063","lp");
  legend->AddEntry(gRaaRap1,"#varUpsilon(1S) in progress","lp");
  if(plotTNP)  legend->AddEntry(gRaaRap1TNP,"#varUpsilon(1S) in progress","lp");
  legend->Draw();

gPad->RedrawAxis();

}
 plot2010();
 plotRAA_uncorr();
 if(plotPA) plotDoubleRatios();
 if(plotEffOrAcc){
 TCanvas *cAccComp1 = new TCanvas("cAccComp1","cAccComp1"); 
 cAccComp1->cd();
 TPad *pComp1 = new TPad("pComp1","pComp1",0.0,0.0,1.0,1.0);
 pComp1->SetBottomMargin(0.12);
 pComp1->SetTopMargin(0.03);
 pComp1->SetRightMargin(0.03);
 pComp1->SetLeftMargin(0.16);
 pComp1->Draw();
 pComp1->cd();
TF1 *f4CompRap = new TF1("f4CompRap","1",0,2.45);
 f4CompRap->SetLineWidth(0);
 f4CompRap->GetXaxis()->SetTitle("|y^{#Upsilon}|");
 f4CompRap->GetYaxis()->SetTitle("#frac{A_{pp}}{A_{PbPb}}");
 f4CompRap->GetYaxis()->SetTitleOffset(1.8);
 f4CompRap->GetYaxis()->SetTitleSize(0.028);
 f4CompRap->GetYaxis()->SetRangeUser(0.,2);
 f4CompRap->GetXaxis()->CenterTitle(kTRUE);
 f4CompRap->Draw();
 TGraphErrors *gCompRap1 = new TGraphErrors(nRapBins_2014,rap2014,R_A_1S_rap,rap2014e,R_A_1S_rape);
 gCompRap1->SetMarkerColor(kRed);
 gCompRap1->SetMarkerStyle(21);
 gCompRap1->SetMarkerSize(1.2);
 TGraphErrors *gCompRap1circle = new TGraphErrors(nRapBins_2014,rap2014,R_A_1S_rap,rap2014e,R_A_1S_rape);
 gCompRap1circle->SetMarkerStyle(25);
 gCompRap1circle->SetMarkerSize(1.2);
 gCompRap1circle->SetLineColor(kBlack);
 gCompRap1->Draw("pe");
 gCompRap1circle->Draw("p");
 f4CompRap->Draw("same");
 gPad->RedrawAxis();
TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
legend->SetTextSize(0.029);
legend->SetFillStyle(0);
legend->SetFillColor(0);
legend->SetBorderSize(0);
legend->SetTextFont(42);
legend->AddEntry(gCompRap1,"#varUpsilon(1S) A_{pp}/A_{PbPb} ","lp");
 legend->Draw();


 TCanvas *cAccComp2 = new TCanvas("cAccComp2","cAccComp2"); 
 cAccComp2->cd();
 TPad *pComp2 = new TPad("pComp2","pComp2",0.0,0.0,1.0,1.0);
 pComp2->SetBottomMargin(0.12);
 pComp2->SetTopMargin(0.03);
 pComp2->SetRightMargin(0.03);
 pComp2->SetLeftMargin(0.16);
 pComp2->Draw();
 pComp2->cd();
TF1 *f4CompPt = new TF1("f4CompPt","1",0,20);
 f4CompPt->SetLineWidth(0);
 f4CompPt->GetXaxis()->SetTitle("p_{T}^{#Upsilon}");
 f4CompPt->GetYaxis()->SetTitle("#frac{A_{pp}}{A_{PbPb}}");
 f4CompPt->GetYaxis()->SetTitleOffset(1.8);
 f4CompPt->GetYaxis()->SetTitleSize(0.028);
 f4CompPt->GetYaxis()->SetRangeUser(0.,2.);
 f4CompPt->GetXaxis()->CenterTitle(kTRUE);
 f4CompPt->Draw();
 TGraphErrors *gCompPt1 = new TGraphErrors(nPtBins_2013,pt,R_A_1S_pt,pte,R_A_1S_pte);
 gCompPt1->SetMarkerColor(kRed);
 gCompPt1->SetMarkerStyle(21);
 gCompPt1->SetMarkerSize(1.2);
 TGraphErrors *gCompPt1circle = new TGraphErrors(nPtBins_2013,pt,R_A_1S_pt,pte,R_A_1S_pte);
 gCompPt1circle->SetMarkerStyle(25);
 gCompPt1circle->SetMarkerSize(1.2);
 gCompPt1circle->SetLineColor(kBlack);
 gCompPt1->Draw("pe");
 gCompPt1circle->Draw("p");
 f4CompPt->Draw("same");
 gPad->RedrawAxis();
TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
legend->SetTextSize(0.029);
legend->SetFillStyle(0);
legend->SetFillColor(0);
legend->SetBorderSize(0);
legend->SetTextFont(42);
legend->AddEntry(gCompPt1,"#varUpsilon(1S) A_{pp}/A_{PbPb} ","lp");
 legend->Draw();

 TCanvas *cAccEffComp1 = new TCanvas("cAccEffComp1","cAccEffComp1"); 
 cAccEffComp1->cd();
 TPad *pAccEffComp1 = new TPad("pAccEffComp1","pAccEffComp1",0.0,0.0,1.0,1.0);
 pAccEffComp1->SetBottomMargin(0.12);
 pAccEffComp1->SetTopMargin(0.03);
 pAccEffComp1->SetRightMargin(0.03);
 pAccEffComp1->SetLeftMargin(0.16);
 pAccEffComp1->Draw();
 pAccEffComp1->cd();
TF1 *f4CompRap = new TF1("f4CompRap","1",0,2.45);
 f4CompRap->SetLineWidth(0);
 f4CompRap->GetXaxis()->SetTitle("|y^{#Upsilon}|");
 f4CompRap->GetYaxis()->SetTitle("#alpha#cdot#epsilon_{pp,AA}");
 f4CompRap->GetYaxis()->SetTitleOffset(1.8);
 f4CompRap->GetYaxis()->SetTitleSize(0.028);
 f4CompRap->GetYaxis()->SetRangeUser(0.,0.5);
 f4CompRap->GetXaxis()->CenterTitle(kTRUE);
 f4CompRap->Draw();
 TGraphErrors *gCompAccEffRap1 = new TGraphErrors(nRapBins_2014,rap2014,Ae_1S_pythia_rap2014,rap2014e,Ae_1S_pythia_rap2014e);
 gCompAccEffRap1->SetMarkerColor(kRed);
 gCompAccEffRap1->SetMarkerStyle(23);
 gCompAccEffRap1->SetMarkerSize(1.2);
 TGraphErrors *gCompAccEffRap1circle = new TGraphErrors(nRapBins_2014,rap2014,Ae_1S_pythia_rap2014,rap2014e,Ae_1S_pythia_rap2014e);
 gCompAccEffRap1circle->SetMarkerStyle(32);
 gCompAccEffRap1circle->SetMarkerSize(1.2);
 gCompAccEffRap1circle->SetLineColor(kBlack);
 gCompAccEffRap1->Draw("pe");
 gCompAccEffRap1circle->Draw("p");
 f4CompRap->Draw("same");
 gPad->RedrawAxis();

 TGraphErrors *gCompAccEffRap1aa = new TGraphErrors(nRapBins_2014,rap2014,Ae_1S_pyquen_rap2014,rap2014e,Ae_1S_pyquen_rap2014e);
 gCompAccEffRap1aa->SetMarkerColor(kRed-7);
 gCompAccEffRap1aa->SetMarkerStyle(23);
 gCompAccEffRap1aa->SetMarkerSize(1.2);
 TGraphErrors *gCompAccEffRap1aacircle = new TGraphErrors(nRapBins_2014,rap2014,Ae_1S_pyquen_rap2014,rap2014e,Ae_1S_pyquen_rap2014e);
 gCompAccEffRap1aacircle->SetMarkerStyle(32);
 gCompAccEffRap1aacircle->SetMarkerSize(1.2);
 gCompAccEffRap1aacircle->SetLineColor(kBlack);
 gCompAccEffRap1aa->Draw("pe");
 gCompAccEffRap1aacircle->Draw("p");
 f4CompRap->Draw("same");
 gPad->RedrawAxis();
 ///2S

 TCanvas *cAccEffComp2 = new TCanvas("cAccEffComp2","cAccEffComp2"); 
 cAccEffComp2->cd();
 pAccEffComp1->Draw();
 pAccEffComp1->cd();
 f4CompRap->Draw();
 TGraphErrors *gCompAccEffRap2 = new TGraphErrors(nRapBins_2014,rap2014,Ae_2S_pythia_rap2014,rap2014e,Ae_2S_pythia_rap2014e);
 gCompAccEffRap2->SetMarkerColor(kBlue);
 gCompAccEffRap2->SetMarkerStyle(21);
 gCompAccEffRap2->SetMarkerSize(1.2);
 TGraphErrors *gCompAccEffRap2circle = new TGraphErrors(nRapBins_2014,rap2014,Ae_2S_pythia_rap2014,rap2014e,Ae_2S_pythia_rap2014e);
 gCompAccEffRap2circle->SetMarkerStyle(25);
 gCompAccEffRap2circle->SetMarkerSize(1.2);
 gCompAccEffRap2circle->SetLineColor(kBlack);
 gCompAccEffRap2->Draw("pe");
 gCompAccEffRap2circle->Draw("p");
 f4CompRap->Draw("same");
 gPad->RedrawAxis();

 TGraphErrors *gCompAccEffRap2aa = new TGraphErrors(nRapBins_2014,rap2014,Ae_2S_pyquen_rap2014,rap2014e,Ae_2S_pyquen_rap2014e);
 gCompAccEffRap2aa->SetMarkerColor(kBlue-9);
 gCompAccEffRap2aa->SetMarkerStyle(21);
 gCompAccEffRap2aa->SetMarkerSize(1.2);
 TGraphErrors *gCompAccEffRap2aacircle = new TGraphErrors(nRapBins_2014,rap2014,Ae_2S_pyquen_rap2014,rap2014e,Ae_2S_pyquen_rap2014e);
 gCompAccEffRap2aacircle->SetMarkerStyle(25);
 gCompAccEffRap2aacircle->SetMarkerSize(1.2);
 gCompAccEffRap2aacircle->SetLineColor(kBlack);
 gCompAccEffRap2aa->Draw("pe");
 gCompAccEffRap2aacircle->Draw("p");
 f4CompRap->Draw("same");
 gPad->RedrawAxis();




 TCanvas *cEffComp1 = new TCanvas("cEffComp1","cEffComp1"); 
 cEffComp1->cd();
 TPad *pComp1 = new TPad("pComp1","pComp1",0.0,0.0,1.0,1.0);
 pComp1->SetBottomMargin(0.12);
 pComp1->SetTopMargin(0.03);
 pComp1->SetRightMargin(0.03);
 pComp1->SetLeftMargin(0.16);
 pComp1->Draw();
 pComp1->cd();
TF1 *f4CompRap = new TF1("f4CompRap","1",0,2.45);
 f4CompRap->SetLineWidth(0);
 f4CompRap->GetXaxis()->SetTitle("|y^{#Upsilon}|");
 f4CompRap->GetYaxis()->SetTitle("#frac{#epsilon_{pp}}{#epsilon_{PbPb}}");
 f4CompRap->GetYaxis()->SetTitleOffset(1.8);
 f4CompRap->GetYaxis()->SetTitleSize(0.028);
 f4CompRap->GetYaxis()->SetRangeUser(0.,1.5);
 f4CompRap->GetXaxis()->CenterTitle(kTRUE);
 f4CompRap->Draw();
 TGraphErrors *gCompRap1 = new TGraphErrors(nRapBins_2014,rap2014,R_e_1S_rap,rap2014e,R_e_1S_rape);
 gCompRap1->SetMarkerColor(kRed);
 gCompRap1->SetMarkerStyle(21);
 gCompRap1->SetMarkerSize(1.2);
 TGraphErrors *gCompRap1circle = new TGraphErrors(nRapBins_2014,rap2014,R_e_1S_rap,rap2014e,R_e_1S_rape);
 gCompRap1circle->SetMarkerStyle(25);
 gCompRap1circle->SetMarkerSize(1.2);
 gCompRap1circle->SetLineColor(kBlack);
 gCompRap1->Draw("pe");
 gCompRap1circle->Draw("p");
 f4CompRap->Draw("same");
 gPad->RedrawAxis();
TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
legend->SetTextSize(0.029);
legend->SetFillStyle(0);
legend->SetFillColor(0);
legend->SetBorderSize(0);
legend->SetTextFont(42);
legend->AddEntry(gCompRap1,"#varUpsilon(1S) #epsilon_{pp}/#epsilon_{PbPb} ","lp");
 legend->Draw();

 TCanvas *cEffComp2 = new TCanvas("cEffComp2","cEffComp2"); 
 cEffComp2->cd();
 TPad *pComp2 = new TPad("pComp2","pComp2",0.0,0.0,1.0,1.0);
 pComp2->SetBottomMargin(0.12);
 pComp2->SetTopMargin(0.03);
 pComp2->SetRightMargin(0.03);
 pComp2->SetLeftMargin(0.16);
 pComp2->Draw();
 pComp2->cd();
TF1 *f4CompPt = new TF1("f4CompPt","1",0,20);
 f4CompPt->SetLineWidth(0);
 f4CompPt->GetXaxis()->SetTitle("p_{T}^{#Upsilon}");
 f4CompPt->GetYaxis()->SetTitle("#frac{#epsilon_{pp}}{#epsilon_{PbPb}}");
 f4CompPt->GetYaxis()->SetTitleOffset(1.8);
 f4CompPt->GetYaxis()->SetTitleSize(0.028);
 f4CompPt->GetYaxis()->SetRangeUser(0.,1.5);
 f4CompPt->GetXaxis()->CenterTitle(kTRUE);
 f4CompPt->Draw();
 TGraphErrors *gCompPt1 = new TGraphErrors(nPtBins_2013,pt,R_e_1S_pt,pte,R_e_1S_pte);
 gCompPt1->SetMarkerColor(kRed);
 gCompPt1->SetMarkerStyle(21);
 gCompPt1->SetMarkerSize(1.2);
 TGraphErrors *gCompPt1circle = new TGraphErrors(nPtBins_2013,pt,R_e_1S_pt,pte,R_e_1S_pte);
 gCompPt1circle->SetMarkerStyle(25);
 gCompPt1circle->SetMarkerSize(1.2);
 gCompPt1circle->SetLineColor(kBlack);
 gCompPt1->Draw("pe");
 gCompPt1circle->Draw("p");
 f4CompPt->Draw("same");
 gPad->RedrawAxis();
TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
legend->SetTextSize(0.029);
legend->SetFillStyle(0);
legend->SetFillColor(0);
legend->SetBorderSize(0);
legend->SetTextFont(42);
legend->AddEntry(gCompPt1,"#varUpsilon(1S) #epsilon_{pp}/#epsilon_{PbPb} ","lp");
 legend->Draw();
 }

}

 float computeRatio(float x, float y) 
  {
    // pass the yield (x), and Acc*eff (y), and computes the corrected yield. then divides by lumi and delta pt or delta rapidity to get the cross section. in case of pbpb, divides by taa*nMB to get the nColl scaled invariant yield.
    float ratio;
    ratio = x/y;
    
    return ratio;
  }
 
float computeRatioError(float x, float y, float xerr, float yerr) 
 {
   //propagate the error of the ratio
   float err = (xerr*xerr)/(x*x) + (yerr*yerr)/(y*y);
  
 // + 2.*(x.getError()*y.getError())/(x.getVal()*y.getVal())*correlation; // can be needed in case of correlations.
   
   return fabs(computeRatio(x,y))*sqrt(err);
 }


float plot2010()
{

  float CS1S_pp_tot;
  float CS1S_pp_tote;
  float CS1S_aa_tot;
  float CS1S_aa_tote;
  float CS1S_aa_cent[nCentBins_2014] = {};
  float CS1S_aa_cente[nCentBins_2014] = {};
  float RAA_1S_cent[nCentBins_2014]={};
  float RAA_1S_cente[nCentBins_2014]={};

  float RAA_1S_tot;
  float RAA_1S_tote;
  float RAA_2S_tot;
  float RAA_2S_tote;
  float RAA_3S_tot;
  float RAA_3S_tote;
  float CS2S_pp_tot;
  float CS2S_pp_tote;
  float CS2S_aa_tot;
  float CS2S_aa_tote;
  float CS3S_pp_tot;
  float CS3S_pp_tote;
  float CS3S_aa_tot;
  float CS3S_aa_tote;

  float CS2S_aa_cent[bin1] = {};
  float CS2S_aa_cente[bin1] = {};
  float RAA_2S_cent[bin1]={};
  float RAA_2S_cente[bin1]={};
 for(int centi =0 ; centi<bin1 ; centi++)
   {
     taa[centi]=taa[centi]*1000;
     CS2S_aa_cent[centi]= computeRatio( N2S_aa_cent4[centi] , Ae_2S_pyquen_cent2014[centi] );
     CS2S_aa_cente[centi] = computeRatioError( N2S_aa_cent4[centi] , Ae_2S_pyquen_cent2014[centi], N2S_aa_cent4e[centi] , Ae_2S_pyquen_cent2014e[centi]);
     CS2S_aa_cent[centi]=CS2S_aa_cent[centi]/(mb_percentage[centi]*N_MB_corr * taa[centi]);
     CS2S_aa_cente[centi]=CS2S_aa_cente[centi]/(mb_percentage[centi]*N_MB_corr * taa[centi]);
     if(centi==0){
       CS2S_pp_tot = computeRatio(N2S_pp_tot4,Ae_2S_pythia_tot);
       CS2S_pp_tote = computeRatioError(N2S_pp_tot4,Ae_2S_pythia_tot,N2S_pp_tot4e,Ae_2S_pythia_tote);
       CS2S_pp_tot = CS2S_pp_tot/L_pp;
       CS2S_pp_tote=CS2S_pp_tote/L_pp;
       CS3S_pp_tot = computeRatio(N3S_pp_tot4,Ae_3S_pythia_tot);
       CS3S_pp_tote = computeRatioError(N3S_pp_tot4,Ae_3S_pythia_tot,N3S_pp_tot4e,Ae_3S_pythia_tote);
       CS3S_pp_tot = CS3S_pp_tot/L_pp;
       CS3S_pp_tote=CS3S_pp_tote/L_pp;
     }
     RAA_2S_cent[centi]= computeRatio( CS2S_aa_cent[centi] , CS2S_pp_tot);
     RAA_2S_cente[centi]= computeRatioError( CS2S_aa_cent[centi] , CS2S_pp_tot,  CS2S_aa_cente[centi] , CS2S_pp_tote);
   }

 for(int centi =0 ; centi<nCentBins_2014; centi++){
     taa2014[centi]=taa2014[centi]*1000;
     CS1S_aa_cent[centi]= computeRatio( N1S_aa_cent3p5[centi] , Ae_1S_pyquen_cent2014[centi] );
     CS1S_aa_cente[centi] = computeRatioError( N1S_aa_cent3p5[centi] , Ae_1S_pyquen_cent2014[centi], N1S_aa_cent3p5e[centi] , Ae_1S_pyquen_cent2014e[centi]);
     CS1S_aa_cent[centi]=CS1S_aa_cent[centi]/(mb_percentage2014[centi]*N_MB_corr * taa2014[centi]);
     CS1S_aa_cente[centi]=CS1S_aa_cente[centi]/(mb_percentage2014[centi]*N_MB_corr * taa2014[centi]);

if(centi==0){       
  CS1S_pp_tot = computeRatio(N1S_pp_tot3p5,Ae_1S_pythia_tot);
  CS1S_pp_tote = computeRatioError(N1S_pp_tot3p5,Ae_1S_pythia_tot,N1S_pp_tot3p5e,Ae_1S_pythia_tote);
  CS1S_pp_tot = CS1S_pp_tot/L_pp;
  CS1S_pp_tote=CS1S_pp_tote/L_pp;
 }
     RAA_1S_cent[centi]= computeRatio( CS1S_aa_cent[centi] , CS1S_pp_tot);
     RAA_1S_cente[centi]= computeRatioError( CS1S_aa_cent[centi] , CS1S_pp_tot,  CS1S_aa_cente[centi] , CS1S_pp_tote);
 }
 CS1S_aa_tot = computeRatio(N1S_aa_tot3p5,Ae_1S_pyquen_tot);
 CS1S_aa_tote = computeRatioError(N1S_aa_tot3p5,Ae_1S_pyquen_tot,N1S_aa_tot3p5e,Ae_1S_pyquen_tote);
 CS1S_aa_tot = CS1S_aa_tot/L_pp;
 CS1S_aa_tote= CS1S_aa_tote/L_pp;
 CS2S_aa_tot = computeRatio(N2S_aa_tot3p5,Ae_2S_pyquen_tot);
 CS2S_aa_tote = computeRatioError(N2S_aa_tot3p5,Ae_2S_pyquen_tot,N2S_aa_tot3p5e,Ae_2S_pyquen_tote);
 CS2S_aa_tot = CS2S_aa_tot/L_pp;
 CS2S_aa_tote=CS2S_aa_tote/L_pp;
 CS3S_aa_tot = computeRatio(N3S_aa_tot4,Ae_3S_pyquen_tot); // careful here
 CS3S_aa_tote = computeRatioError(N3S_aa_tot4,Ae_3S_pyquen_tot,N3S_aa_tot4e,Ae_3S_pyquen_tote);
 CS3S_aa_tot = CS3S_aa_tot/L_pp;
 CS3S_aa_tote=CS3S_aa_tote/L_pp;

 RAA_1S_tot = computeRatio(CS1S_aa_tot,CS1S_pp_tot);
 RAA_1S_tote = computeRatioError(CS1S_aa_tot,CS1S_pp_tot,CS1S_aa_tote,CS1S_pp_tote);
 RAA_2S_tot = computeRatio(CS2S_aa_tot,CS2S_pp_tot);
 RAA_2S_tote = computeRatioError(CS2S_aa_tot,CS2S_pp_tot,CS2S_aa_tote,CS2S_pp_tote);
 RAA_3S_tot = computeRatio(CS3S_aa_tot,CS3S_pp_tot);
 RAA_3S_tote = computeRatioError(CS3S_aa_tot,CS3S_pp_tot,CS3S_aa_tote,CS3S_pp_tote);
 cout<<" ---####---  total_sigma(1S)_pp = "<<CS1S_pp_tot <<" +/- " <<CS1S_pp_tote<<endl;
cout << "  --- Cross section in PbPb vs. nPart ---" << endl;
 for(int j =0 ; j<nCentBins_2014 ; j++)
   {
     cout <<"bin="<< j << "' ,sigma(1S)_PbPb = "<< CS1S_aa_cent[j] <<" +/- "<<CS1S_aa_cente[j]<<" b" << endl;
   }

cout << "  --- 1S RAA vs. nPart ---" << endl;
 for(int j =0 ; j<nCentBins_2014; j++)
   {
     cout <<"bin="<< j << "' , Raa = "<< RAA_1S_cent[j] <<" +/- "<< RAA_1S_cente[j]<< endl;
   }

 cout<<" ---####---  total_sigma(2S)_pp = "<<CS2S_pp_tot <<" +/- " <<CS2S_pp_tote<<endl;
cout << "  --- Cross section in PbPb vs. nPart ---" << endl;
 for(int j =0 ; j<bin1 ; j++)
   {
     cout <<"bin="<< j << "' ,sigma(2S)_PbPb = "<< CS2S_aa_cent[j] <<" +/- "<<CS2S_aa_cente[j]<<" b" << endl;
   }

cout << "  --- 2S RAA vs. nPart ---" << endl;
 for(int j =0 ; j<bin1 ; j++)
   {
     cout <<"bin="<< j << "' , Raa = "<< RAA_2S_cent[j] <<" +/- "<< RAA_2S_cente[j]<< endl;
   }


//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 28 03:36:31 2014) by ROOT version5.34/02
 if(plotRAA){
   TCanvas *c1 = new TCanvas("c1", "c1",423,55,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   //   c1->Range(-58.8957,-0.2117647,431.9018);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(0);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.12);
   c1->SetRightMargin(0.065);
   c1->SetTopMargin(0.03);
   c1->SetBottomMargin(0.12);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TF1 *f4 = new TF1("f4","1",0,400);
   f4->SetFillColor(19);
   f4->SetFillStyle(0);
   f4->SetMarkerStyle(20);
   f4->SetMarkerSize(0.8);
   f4->SetLineWidth(1);
   f4->GetXaxis()->SetTitle("N_{part}");
   f4->GetXaxis()->CenterTitle(true);
   f4->GetXaxis()->SetLabelFont(42);
   f4->GetXaxis()->SetTitleSize(0.048);
   f4->GetXaxis()->SetTitleOffset(1.15);
   f4->GetXaxis()->SetTitleFont(42);
   f4->GetYaxis()->SetTitle("R_{AA}");
   f4->GetYaxis()->SetLabelFont(42);
   f4->GetYaxis()->SetTitleSize(0.048);
   f4->GetYaxis()->SetTitleFont(42);
   f4->Draw("axis");
   
   TGraphErrors *gre = new TGraphErrors(7);
   // gre->SetName("Graph");
   // gre->SetTitle("Graph");
   gre->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#99ff99");
   gre->SetLineColor(ci);
   gre->SetLineWidth(25);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(0);
   gre->SetPoint(0,22.1,1.005);
   gre->SetPointError(0,0,0.176);
   gre->SetPoint(1,86.3,0.59);
   gre->SetPointError(1,0,0.086);
   gre->SetPoint(2,130,0.681);
   gre->SetPointError(2,0,0.085);
   gre->SetPoint(3,187.1,0.614);
   gre->SetPointError(3,0,0.075);
   gre->SetPoint(4,261.4,0.484);
   gre->SetPointError(4,0,0.049);
   gre->SetPoint(5,329.4,0.432);
   gre->SetPointError(5,0,0.046);
   gre->SetPoint(6,381.3,0.411);
   gre->SetPointError(6,0,0.048);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph_Graph1",100,0,417.22);
   // Graph_Graph1->SetMinimum(0.2812);
   // Graph_Graph1->SetMaximum(1.2628);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->SetMarkerSize(0.8);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.048);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.048);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.048);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   gre->Draw("e");
   
   gre = new TGraphErrors(7);
   gre->SetName("Graph_Graph1");
   gre->SetTitle("Graph_Graph1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(33);
   gre->SetMarkerSize(2);
   gre->SetPoint(0,20.1,1.005);
   gre->SetPointError(0,0,0.121);
   gre->SetPoint(1,84.3,0.59);
   gre->SetPointError(1,0,0.096);
   gre->SetPoint(2,128,0.681);
   gre->SetPointError(2,0,0.069);
   gre->SetPoint(3,185.1,0.614);
   gre->SetPointError(3,0,0.053);
   gre->SetPoint(4,259.4,0.484);
   gre->SetPointError(4,0,0.04);
   gre->SetPoint(5,327.4,0.432);
   gre->SetPointError(5,0,0.048);
   gre->SetPoint(6,379.3,0.411);
   gre->SetPointError(6,0,0.043);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph_Graph2",100,0,417.22);
   // Graph_Graph2->SetMinimum(0.2922);
   // Graph_Graph2->SetMaximum(1.2018);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   Graph_Graph2->SetMarkerStyle(20);
   Graph_Graph2->SetMarkerSize(0.8);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.048);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.048);
   Graph_Graph2->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.048);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph2);
   
   gre->Draw("pe");
   
   gre = new TGraphErrors(7);
   gre->SetName("Graph_Graph2");
   gre->SetTitle("Graph_Graph2");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(27);
   gre->SetMarkerSize(2);
   gre->SetPoint(0,20.1,1.005);
   gre->SetPointError(0,0,0.121);
   gre->SetPoint(1,84.3,0.59);
   gre->SetPointError(1,0,0.096);
   gre->SetPoint(2,128,0.681);
   gre->SetPointError(2,0,0.069);
   gre->SetPoint(3,185.1,0.614);
   gre->SetPointError(3,0,0.053);
   gre->SetPoint(4,259.4,0.484);
   gre->SetPointError(4,0,0.04);
   gre->SetPoint(5,327.4,0.432);
   gre->SetPointError(5,0,0.048);
   gre->SetPoint(6,379.3,0.411);
   gre->SetPointError(6,0,0.043);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph_Graph3",100,0,417.22);
   // Graph_Graph3->SetMinimum(0.2922);
   // Graph_Graph3->SetMaximum(1.2018);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);
   Graph_Graph3->SetMarkerStyle(20);
   Graph_Graph3->SetMarkerSize(0.8);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.048);
   Graph_Graph3->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.048);
   Graph_Graph3->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.048);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph3);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(3);
   gre->SetName("Graph_Graph3");
   gre->SetTitle("Graph_Graph3");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff99ff");
   gre->SetLineColor(ci);
   gre->SetLineWidth(25);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(0);
   gre->SetPoint(0,64.2184,0.677);
   gre->SetPointError(0,0,0.107);
   gre->SetPoint(1,261.4178,0.842);
   gre->SetPointError(1,0,0.13);
   gre->SetPoint(2,355.3528,0.454);
   gre->SetPointError(2,0,0.079);
   

   
 

   TBox *box = new TBox(385,0.864,400,1.136);

   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);
   box->Draw();

   
   TF1 *f4 = new TF1("f4","1",0,400);
   f4->SetFillColor(19);
   f4->SetFillStyle(0);
   f4->SetMarkerStyle(20);
   f4->SetMarkerSize(0.8);
   f4->SetLineWidth(1);
   f4->GetXaxis()->SetTitle("N_{part}");
   f4->GetXaxis()->CenterTitle(true);
   f4->GetXaxis()->SetLabelFont(42);
   f4->GetXaxis()->SetTitleSize(0.048);
   f4->GetXaxis()->SetTitleOffset(1.15);
   f4->GetXaxis()->SetTitleFont(42);
   f4->GetYaxis()->SetTitle("R_{AA}");
   f4->GetYaxis()->SetLabelFont(42);
   f4->GetYaxis()->SetTitleSize(0.048);
   f4->GetYaxis()->SetTitleFont(42);
   f4->Draw("same");
   
   TH1F *Graph = new TH1F("Graph","Graph",100,0,417.22);
   // Graph->SetMinimum(0.2812);
   // Graph->SetMaximum(1.2628);
   Graph->SetDirectory(0);
   Graph->SetStats(0);
   Graph->SetMarkerStyle(20);
   Graph->SetMarkerSize(0.8);
   Graph->GetXaxis()->SetLabelFont(42);
   Graph->GetXaxis()->SetTitleSize(0.048);
   Graph->GetXaxis()->SetTitleOffset(1.15);
   Graph->GetXaxis()->SetTitleFont(42);
   Graph->GetYaxis()->SetLabelFont(42);
   Graph->GetYaxis()->SetTitleSize(0.048);
   Graph->GetYaxis()->SetTitleOffset(1.2);
   Graph->GetYaxis()->SetTitleFont(42);
   Graph->GetZaxis()->SetLabelFont(42);
   Graph->GetZaxis()->SetTitleSize(0.048);
   Graph->GetZaxis()->SetTitleFont(42);
   Graph->Draw("sameaxis");
   TLatex *   tex = new TLatex(220,1.85,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
   tex->DrawLatex(220,1.75,"Work in progress");
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
 



   TGraphErrors *gcent1 = new TGraphErrors(nCentBins_2014,nPart2014,RAA_1S_cent,centnoErr,RAA_1S_cente);
 gcent1->SetMarkerColor(kOrange+7);
 gcent1->SetMarkerStyle(33);
 gcent1->SetMarkerSize(2);
 TGraphErrors *gcent1circle = new TGraphErrors(nCentBins_2014,nPart2014,RAA_1S_cent,centnoErr,RAA_1S_cente);
 gcent1circle->SetMarkerStyle(27);
 gcent1circle->SetMarkerSize(2);
 gcent1circle->SetLineColor(kBlack);
 gcent1->Draw("pe");
 gcent1circle->Draw("p");
 f4->Draw("same");
 // f4->GetYaxis()->SetRangeUser(0.0,1.5);
 gPad->RedrawAxis();

   TGraphErrors *gcent2 = new TGraphErrors(bin1,cent,RAA_2S_cent,centnoErr,RAA_2S_cente);
 gcent2->SetMarkerColor(kAzure-8);
 gcent2->SetMarkerStyle(20);
 gcent2->SetMarkerSize(1.2);
 TGraphErrors *gcent2circle = new TGraphErrors(bin1,cent,RAA_2S_cent,centnoErr,RAA_2S_cente);
 gcent2circle->SetMarkerStyle(24);
 gcent2circle->SetMarkerSize(1.2);
 gcent2circle->SetLineColor(kBlack);
 gcent2->Draw("pe");
 gcent2circle->Draw("p");
 f4->Draw("same");
 gPad->RedrawAxis();

 //x axis for old points is good-2.
   TGraphErrors *gre = new TGraphErrors(7);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   gre = new TGraphErrors(7);
   gre->SetName("Graph_2S");
   gre->SetTitle("Graph_2S");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,20.1,0.3);
   gre->SetPointError(0,0,0.157);
   gre->SetPoint(1,84.3,0.251);
   gre->SetPointError(1,0,0.138);
   gre->SetPoint(2,128,0.237);
   gre->SetPointError(2,0,0.098);
   gre->SetPoint(3,185.1,0.26);
   gre->SetPointError(3,0,0.079);
   gre->SetPoint(4,259.4,0.068);
   gre->SetPointError(4,0,0.053);
   gre->SetPoint(5,327.4,0.044);
   gre->SetPointError(5,0,0.06);
   gre->SetPoint(6,379.3,0.111);
   gre->SetPointError(6,0,0.061);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,417.22);
   Graph_Graph6->SetMinimum(-0.0633);
   Graph_Graph6->SetMaximum(0.5043);
   Graph_Graph6->SetDirectory(0);
   Graph_Graph6->SetStats(0);
   Graph_Graph6->SetMarkerStyle(20);
   Graph_Graph6->SetMarkerSize(0.8);
   Graph_Graph6->GetXaxis()->SetLabelFont(42);
   Graph_Graph6->GetXaxis()->SetTitleSize(0.048);
   Graph_Graph6->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph6->GetXaxis()->SetTitleFont(42);
   Graph_Graph6->GetYaxis()->SetLabelFont(42);
   Graph_Graph6->GetYaxis()->SetTitleSize(0.048);
   Graph_Graph6->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph6->GetYaxis()->SetTitleFont(42);
   Graph_Graph6->GetZaxis()->SetLabelFont(42);
   Graph_Graph6->GetZaxis()->SetTitleSize(0.048);
   Graph_Graph6->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph6);
gre->GetYaxis()->SetRangeUser(0,1.8);
   gre->Draw("p");
   
   TLegend *leg = new TLegend(0.53,0.6,1.,0.8);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.031);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph_Graph1","#varUpsilon(1S) 2011","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);


   ci = TColor::GetColor("#cc00cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);

   entry=leg->AddEntry(gcent1,"#varUpsilon(1S) 2014","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);

   entry=leg->AddEntry("Graph_2S","#varUpsilon(2S) 2011","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);

   entry=leg->AddEntry(gcent2,"#varUpsilon(2S) 2014","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);

   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
 }
 cout << "Raa_1S = "<<RAA_1S_tot<<" +/- "<<RAA_1S_tote<<" stat."<<endl;
 cout << "Raa_2S = "<<RAA_2S_tot<<" +/- "<<RAA_2S_tote<<" stat."<<endl;
 cout << "Raa_3S = "<<RAA_3S_tot<<" +/- "<<RAA_3S_tote<<" stat."<<endl;
}
//only for 1S - > uncorrected RAA
void plotRAA_uncorr(){

  float pt [nPtBins_2013] = {1.25, 3.75, 6.5, 10., 16.};
  float pte[nPtBins_2013] = {1.25, 1.25, 1.5, 2., 4.};
  float deltaPt[nPtBins_2013]   = {2.5,2.5,3,4,8};
  
  float deltaRap[nRapBins_2013]  = {0.8,0.6,0.6,1,1.8};
  float deltaRapEven[nRapBins_2014] = {0.8,0.8,0.8,0.8,0.8,0.8};

  float uncorrRAA_1S_pt[nPtBins_2013] = {};
  float uncorrRAA_1S_pte[nPtBins_2013] = {};
  float uncorrRAA_1S_ptU[nPtBins_2013] = {};
  float uncorrRAA_1S_pteU[nPtBins_2013] = {};
  float uncorrRAA_1S_rap2014[nPtBins_2013] = {};
  float uncorrRAA_1S_rap2014e[nPtBins_2013] = {};
  float CS1S_pp_pt[nPtBins_2013] = {};
  float CS1S_pp_pte[nPtBins_2013] = {};
  float CS1S_pp_ptU[nPtBins_2013] = {};
  float CS1S_pp_pteU[nPtBins_2013] = {};
  float CS1S_pp_rap2014[nRapBins_2014] = {};
  float CS1S_pp_rap2014e[nRapBins_2014] = {};

  float CS1S_aa_pt[nPtBins_2013] = {};
  float CS1S_aa_pte[nPtBins_2013] = {};

  float CS1S_aa_ptU[nPtBins_2013] = {};
  float CS1S_aa_pteU[nPtBins_2013] = {};
  float CS1S_aa_rap2014[nRapBins_2014] = {};
  float CS1S_aa_rap2014e[nRapBins_2014] = {};



  for(int i=0; i<nRapBins_2014;i++)
    {
      CS1S_aa_rap2014[i]=N1S_aa_rap3p5_2014[i]/(N_MB_corr * T_AA_b * deltaRapEven[i]);
      CS1S_aa_rap2014e[i]=N1S_aa_rap3p5_2014e[i]/(N_MB_corr * T_AA_b * deltaRapEven[i]);
      CS1S_pp_rap2014[i]=N1S_pp_rap3p5_2014[i]/(L_pp*deltaRapEven[i]);
      CS1S_pp_rap2014e[i]=N1S_pp_rap3p5_2014e[i]/(L_pp*deltaRapEven[i]);
      uncorrRAA_1S_rap2014[i]=computeRatio(CS1S_aa_rap2014[i] , CS1S_pp_rap2014[i]);
      uncorrRAA_1S_rap2014e[i]=computeRatioError(CS1S_aa_rap2014[i],CS1S_pp_rap2014[i], CS1S_aa_rap2014e[i], CS1S_pp_rap2014e[i]);
    }

  for(int i = 0; i<nPtBins_2013 ; i++)
    {
      CS1S_pp_pt[i]=N1S_pp_pt3p5[i]/(L_pp*deltaPt[i]);
      CS1S_pp_pte[i]=N1S_pp_pt3p5e[i]/(L_pp*deltaPt[i]);
      CS1S_aa_pt[i]=N1S_aa_pt3p5[i]/(N_MB_corr * T_AA_b * deltaPt[i]);
      CS1S_aa_pte[i]=N1S_aa_pt3p5e[i]/(N_MB_corr * T_AA_b * deltaPt[i]);
      
      CS1S_pp_ptU[i]=N1S_pp_pt3p5U[i]/(L_pp*deltaPt[i]);
      CS1S_pp_pteU[i]=N1S_pp_pt3p5eU[i]/(L_pp*deltaPt[i]);
      CS1S_aa_ptU[i]=N1S_aa_pt3p5U[i]/(N_MB_corr * T_AA_b * deltaPt[i]);
      CS1S_aa_pteU[i]=N1S_aa_pt3p5eU[i]/(N_MB_corr * T_AA_b * deltaPt[i]);
      
      uncorrRAA_1S_pt[i]=computeRatio(CS1S_aa_pt[i] , CS1S_pp_pt[i]);
      uncorrRAA_1S_pte[i]=computeRatioError(CS1S_aa_pt[i],CS1S_pp_pt[i], CS1S_aa_pte[i], CS1S_pp_pte[i]);
      uncorrRAA_1S_ptU[i]=computeRatio(CS1S_aa_ptU[i] , CS1S_pp_ptU[i]);
      uncorrRAA_1S_pteU[i]=computeRatioError(CS1S_aa_ptU[i],CS1S_pp_ptU[i], CS1S_aa_pteU[i], CS1S_pp_pteU[i]);
    }
  
//drawing cross sections

  if (plotUncorrected){
 TCanvas *cptu = new TCanvas("cptu","cptu"); 
 cptu->cd();
 TPad *ppt1 = new TPad("ppt1","ppt1",0.0,0.0,1.0,1.0);
 ppt1->SetBottomMargin(0.12);
 ppt1->SetTopMargin(0.03);
 ppt1->SetRightMargin(0.03);
 ppt1->SetLeftMargin(0.16);
 ppt1->SetLogy();
 ppt1->Draw();
 ppt1->cd();
 TF1 *f4Pt = new TF1("f4Pt","0.000000001",0,21);
 f4Pt->SetLineWidth(0);
 f4Pt->GetYaxis()->SetTitleOffset(2);
 f4Pt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} (GeV/c)");		
 f4Pt->GetYaxis()->SetTitle("#frac{1}{L_{pp,PbPb}}#frac{N}{#Deltap_{T}} (b/(GeV/c))");


 f4Pt->GetYaxis()->SetTitleSize(0.028);
 //f4Pt->GetYaxis()->SetRangeUser(0.01,.09);
 f4Pt->GetXaxis()->CenterTitle(kTRUE);
 f4Pt->Draw();
 //one pad to draw PbPb yields,
 TGraphErrors *gpt1 = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_pt,pte,CS1S_aa_pte);
 gpt1->SetMarkerColor(8);
 gpt1->SetMarkerStyle(33);
 gpt1->SetMarkerSize(2);



 TGraphErrors *gpt1circle = new TGraphErrors(nPtBins_2013,pt,CS1S_aa_pt,pte,CS1S_aa_pte);
 gpt1circle->SetMarkerStyle(27);
 gpt1circle->SetMarkerSize(2);
 gpt1circle->SetLineColor(kBlack);
 gpt1->Draw("pe");
 gpt1circle->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();
    
 TGraphErrors *gpt1pp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 gpt1pp->SetMarkerColor(kAzure-9);
 gpt1pp->SetMarkerStyle(21);
 gpt1pp->SetMarkerSize(1.2);
 TGraphErrors *gpt1circlepp = new TGraphErrors(nPtBins_2013,pt,CS1S_pp_pt,pte,CS1S_pp_pte);
 gpt1circlepp->SetMarkerStyle(25);
 gpt1circlepp->SetMarkerSize(1.22);
 gpt1circlepp->SetLineColor(kBlack);

 gpt1pp->Draw("pe");
 gpt1circlepp->Draw("p");
 f4Pt->Draw("same");
 gPad->RedrawAxis();

 TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 legend->AddEntry(gpt1pp,"#varUpsilon(1S), pp ","lp");
 legend->AddEntry(gpt1,"#varUpsilon(1S), PbPb ","lp");

 legend->Draw();
 TLatex *l1CMSpt = new TLatex(8,0.0000000008, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.032);
 l1CMSpt->Draw();

 
 TLatex *lyL= new TLatex(2,0.0000000008,"L_{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 
 lyL->SetTextSize(0.029);
 lyL->DrawLatex(2,0.0000000005,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
 lyL->SetTextSize(0.029);
 lyL->DrawLatex(2,0.000000002,"Raw Yields");
 lyL->Draw();
  }

  if(plotRAA && plotUncorrected){
    TCanvas *cRaaptu = new TCanvas("cRaaptu","cRaaptu"); 
    cRaaptu->cd();
    TPad *ppt2 = new TPad("ppt2","ppt2",0.0,0.0,1.0,1.0);
    ppt2->SetBottomMargin(0.12);
    ppt2->SetTopMargin(0.03);
    ppt2->SetRightMargin(0.03);
    ppt2->SetLeftMargin(0.16);
    ppt2->Draw();
    ppt2->cd();
    //one pad to draw RaaPt!
    TF1 *f4RaaPt = new TF1("f4RaaPt","1",0,21);
    f4RaaPt->SetLineWidth(0);
    f4RaaPt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} ");
    f4RaaPt->GetYaxis()->SetTitle("uncorrected R_{AA}");
    f4RaaPt->GetYaxis()->SetTitleOffset(1.8);
    f4RaaPt->GetYaxis()->SetTitleSize(0.028);
    f4RaaPt->GetYaxis()->SetRangeUser(0.,1.3);
    f4RaaPt->GetXaxis()->CenterTitle(kTRUE);
    f4RaaPt->Draw();
    TGraphErrors *gRaaPt1 = new TGraphErrors(nPtBins_2013,pt,uncorrRAA_1S_pt,pte,uncorrRAA_1S_pte);
    gRaaPt1->SetMarkerColor(8);
    gRaaPt1->SetMarkerStyle(33);
    gRaaPt1->SetMarkerSize(2);
    TGraphErrors *gRaaPt1circle = new TGraphErrors(nPtBins_2013,pt,uncorrRAA_1S_pt,pte,uncorrRAA_1S_pte);
    gRaaPt1circle->SetMarkerStyle(27);
    gRaaPt1circle->SetMarkerSize(2);
    gRaaPt1circle->SetLineColor(kBlack);
    gRaaPt1->Draw("pe");
    gRaaPt1circle->Draw("p");
    f4RaaPt->Draw("same");
    TLatex *l1CMSpt = new TLatex(10,0.8, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
    l1CMSpt->SetTextFont(42);
    l1CMSpt->SetTextSize(0.032);
    l1CMSpt->Draw();
    TLatex *lyLPT= new TLatex(10,0.6,"L_{int}^{PbPb} = 150 #mub^{-1}; L_{int}^{pp} = 5.4 pb^{-1};");
    lyLPT->SetTextSize(0.032);
    lyLPT->Draw();
    lyLPT->SetTextSize(0.04);
    lyLPT->DrawLatex(7,0.15,"Uncorrected");
    ppt2->Update();
    gPad->RedrawAxis();

    TCanvas *cRaaptUnfold = new TCanvas("cRaaptUnfold","cRaaptUnfold"); 
    cRaaptUnfold->cd();
    TPad *ppt2 = new TPad("ppt2","ppt2",0.0,0.0,1.0,1.0);
    ppt2->SetBottomMargin(0.12);
    ppt2->SetTopMargin(0.03);
    ppt2->SetRightMargin(0.03);
    ppt2->SetLeftMargin(0.16);
    ppt2->Draw();
    ppt2->cd();
    TF1 *f4RaaPt = new TF1("f4RaaPt","1",0,21);
    f4RaaPt->SetLineWidth(0);
    f4RaaPt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} ");
    f4RaaPt->GetYaxis()->SetTitle("uncorrected, but unfolded R_{AA}");
    f4RaaPt->GetYaxis()->SetTitleOffset(1.8);
    f4RaaPt->GetYaxis()->SetTitleSize(0.028);
    f4RaaPt->GetYaxis()->SetRangeUser(0.,1.3);
    f4RaaPt->GetXaxis()->CenterTitle(kTRUE);
    f4RaaPt->Draw();
    TGraphErrors *gRaaPt1 = new TGraphErrors(nPtBins_2013,pt,uncorrRAA_1S_pt,pte,uncorrRAA_1S_pte);
    gRaaPt1->SetMarkerColor(8);
    gRaaPt1->SetMarkerStyle(33);
    gRaaPt1->SetMarkerSize(2);
    TGraphErrors *gRaaPt1circle = new TGraphErrors(nPtBins_2013,pt,uncorrRAA_1S_pt,pte,uncorrRAA_1S_pte);
    gRaaPt1circle->SetMarkerStyle(27);
    gRaaPt1circle->SetMarkerSize(2);
    gRaaPt1circle->SetLineColor(kBlack);
    gRaaPt1->Draw("pe");
    gRaaPt1circle->Draw("p");

    TGraphErrors *gRaaPt1U = new TGraphErrors(nPtBins_2013,pt,uncorrRAA_1S_ptU,pte,uncorrRAA_1S_pteU);
    gRaaPt1U->SetMarkerColor(kBlue-3);
    gRaaPt1U->SetMarkerStyle(33);
    gRaaPt1U->SetMarkerSize(2);
    TGraphErrors *gRaaPt1Ucircle = new TGraphErrors(nPtBins_2013,pt,uncorrRAA_1S_ptU,pte,uncorrRAA_1S_pteU);
    gRaaPt1Ucircle->SetMarkerStyle(27);
    gRaaPt1Ucircle->SetMarkerSize(2);
    gRaaPt1Ucircle->SetLineColor(kBlack);
    gRaaPt1U->Draw("pe");
    gRaaPt1Ucircle->Draw("p");
    f4RaaPt->Draw("same");
    TLatex *l1CMSpt = new TLatex(10,0.8, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
    l1CMSpt->SetTextFont(42);
    l1CMSpt->SetTextSize(0.032);
    l1CMSpt->Draw();
    TLatex *lyLPT= new TLatex(10,0.6,"L_{int}^{PbPb} = 150 #mub^{-1}; L_{int}^{pp} = 5.4 pb^{-1};");
    lyLPT->SetTextSize(0.032);
    lyLPT->Draw();
    lyLPT->SetTextSize(0.04);
    lyLPT->DrawLatex(7,0.15,"Uncorrected");
    ppt2->Update();
    TLegend *legend = new TLegend(0.2,0.65,0.4,0.75);
    legend->SetTextSize(0.029);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->AddEntry(gRaaPt1U,"#varUpsilon(1S) unfolded","lp");
    legend->AddEntry(gRaaPt1,"#varUpsilon(1S) raw","lp");
    legend->Draw();
    gPad->RedrawAxis();
  

  }


  if(plotUncorrected){
    TCanvas *crapu = new TCanvas("crapu","crapu"); 
    crapu->cd();
    TPad *prap1 = new TPad("prap1","prap1",0.0,0.0,1.0,1.0);
    prap1->SetBottomMargin(0.12);
    prap1->SetTopMargin(0.03);
    prap1->SetRightMargin(0.03);
    prap1->SetLeftMargin(0.16);
    prap1->SetLogy();
    prap1->Draw();
    prap1->cd();
    TF1 *f4Rap = new TF1("f4Rap","0.000000001",0,2.45);
    f4Rap->SetLineWidth(0);
    f4Rap->GetYaxis()->SetTitleOffset(2);
    f4Rap->GetXaxis()->SetTitle("y^{#Upsilon_{cand.}}");		
    f4Rap->GetYaxis()->SetTitle("#frac{1}{L_{pp,PbPb}}#frac{N}{#Delta_{y}} (b/(GeV/c))");


    f4Rap->GetYaxis()->SetTitleSize(0.028);
    //f4Rap->GetYaxis()->SetRangeUser(0.01,.09);
    f4Rap->GetXaxis()->CenterTitle(kTRUE);
    f4Rap->Draw();
    //one pad to draw PbPb yields,
    TGraphErrors *grap1 = new TGraphErrors(nRapBins_2013,rap,CS1S_aa_rap2014,rape,CS1S_aa_rap2014e);
    grap1->SetMarkerColor(8);
    grap1->SetMarkerStyle(33);
    grap1->SetMarkerSize(2);



    TGraphErrors *grap1circle = new TGraphErrors(nRapBins_2013,rap,CS1S_aa_rap2014,rape,CS1S_aa_rap2014e);
    grap1circle->SetMarkerStyle(27);
    grap1circle->SetMarkerSize(2);
    grap1circle->SetLineColor(kBlack);
    grap1->Draw("pe");
    grap1circle->Draw("p");
    f4Rap->Draw("same");
    gPad->RedrawAxis();
    
    TGraphErrors *grap1pp = new TGraphErrors(nRapBins_2014,rap,CS1S_pp_rap2014,rape,CS1S_pp_rap2014e);
    grap1pp->SetMarkerColor(kAzure-9);
    grap1pp->SetMarkerStyle(21);
    grap1pp->SetMarkerSize(1.2);
    TGraphErrors *grap1circlepp = new TGraphErrors(nRapBins_2014,rap,CS1S_pp_rap2014,rape,CS1S_pp_rap2014e);
    grap1circlepp->SetMarkerStyle(25);
    grap1circlepp->SetMarkerSize(1.22);
    grap1circlepp->SetLineColor(kBlack);
    grap1pp->Draw("pe");
    grap1circlepp->Draw("p");
    f4Rap->Draw("same");
    gPad->RedrawAxis();

    TLegend *legend = new TLegend(0.482,0.84,0.88,0.7);
    legend->SetTextSize(0.029);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->AddEntry(grap1pp,"#varUpsilon(1S), pp ","lp");
    legend->AddEntry(grap1,"#varUpsilon(1S), PbPb ","lp");

    legend->Draw();
    TLatex *l1CMSrap = new TLatex(0.8,0.0000000008, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
    l1CMSrap->SetTextFont(42);
    l1CMSrap->SetTextSize(0.032);
    l1CMSrap->Draw();

 
    TLatex *lyL= new TLatex(1.2,0.0000000008,"L_{PbPb} = 150 #mub^{-1}; |y| < 2.4");
 
    lyL->SetTextSize(0.029);
    lyL->DrawLatex(1.2,0.0000000005,"L_{pp} = 5.4 pb^{-1}; |y| < 2.4");
    lyL->SetTextSize(0.029);
    lyL->DrawLatex(1.2,0.000000002,"Raw Yields");
    lyL->Draw();
  }


 if (plotUncorrected && plotRAA){
 // Uncorrected ratio of Rapidity-binned Cross-Sections
   
   TCanvas *cRaarapu = new TCanvas("cRaarapu","cRaarapu"); 
   cRaarapu->cd();
   TPad *prap2 = new TPad("prap2","prap2",0.0,0.0,1.0,1.0);
   prap2->SetBottomMargin(0.12);
   prap2->SetTopMargin(0.03);
   prap2->SetRightMargin(0.03);
   prap2->SetLeftMargin(0.16);
   prap2->Draw();
   prap2->cd();
   //one pad to draw RaaRap!
   TF1 *f4RaaRap = new TF1("f4RaaRap","1",0,2.45);
   f4RaaRap->SetLineWidth(0);
   f4RaaRap->GetXaxis()->SetTitle("y^{#Upsilon_{cand.}} ");
   f4RaaRap->GetYaxis()->SetTitle("uncorrected R_{AA}");
   f4RaaRap->GetYaxis()->SetTitleOffset(1.8);
   f4RaaRap->GetYaxis()->SetTitleSize(0.028);
   f4RaaRap->GetYaxis()->SetRangeUser(0.,1.3);
   f4RaaRap->GetXaxis()->CenterTitle(kTRUE);
   f4RaaRap->Draw();
   TGraphErrors *gRaaRap1 = new TGraphErrors(nRapBins_2014,rap2014,uncorrRAA_1S_rap2014,rap2014e,uncorrRAA_1S_rap2014e);
   gRaaRap1->SetMarkerColor(8);
   gRaaRap1->SetMarkerStyle(33);
   gRaaRap1->SetMarkerSize(2);
   TGraphErrors *gRaaRap1circle = new TGraphErrors(nRapBins_2014,rap2014,uncorrRAA_1S_rap2014,rap2014e,uncorrRAA_1S_rap2014e);
   gRaaRap1circle->SetMarkerStyle(27);
   gRaaRap1circle->SetMarkerSize(2);
   gRaaRap1circle->SetLineColor(kBlack);
   gRaaRap1->Draw("pe");
   gRaaRap1circle->Draw("p");
   f4RaaRap->Draw("same");
   TLatex *l1CMSrap = new TLatex(1,0.8, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
   l1CMSrap->SetTextFont(42);
   l1CMSrap->SetTextSize(0.032);
   l1CMSrap->Draw();
   TLatex *lyLRAP= new TLatex(1,0.6,"L_{int}^{PbPb} = 150 #mub^{-1}; L_{int}^{pp} = 5.4 pb^{-1};");
   lyLRAP->SetTextSize(0.032);
   lyLRAP->Draw();
   lyLRAP->SetTextSize(0.04);
   lyLRAP->DrawLatex(0.7,0.15,"Uncorrected");
   prap2->Update();
   
 gPad->RedrawAxis();
 }


}





////////
/*2010plots


   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph_Graph4",100,35.10496,384.4663);
   Graph_Graph4->SetMinimum(0.3153);
   Graph_Graph4->SetMaximum(1.0317);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);
   Graph_Graph4->SetMarkerStyle(20);
   Graph_Graph4->SetMarkerSize(0.8);
   Graph_Graph4->GetXaxis()->SetLabelFont(42);
   Graph_Graph4->GetXaxis()->SetTitleSize(0.048);
   Graph_Graph4->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph4->GetXaxis()->SetTitleFont(42);
   Graph_Graph4->GetYaxis()->SetLabelFont(42);
   Graph_Graph4->GetYaxis()->SetTitleSize(0.048);
   Graph_Graph4->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph4->GetYaxis()->SetTitleFont(42);
   Graph_Graph4->GetZaxis()->SetLabelFont(42);
   Graph_Graph4->GetZaxis()->SetTitleSize(0.048);
   Graph_Graph4->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph4);
   
   gre->Draw("e");
   
   gre = new TGraphErrors(3);
   gre->SetName("Graph_Graph4");
   gre->SetTitle("Graph_Graph4");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#cc00cc");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,64.2184,0.677);
   gre->SetPointError(0,0,0.154);
   gre->SetPoint(1,261.4178,0.842);
   gre->SetPointError(1,0,0.212);
   gre->SetPoint(2,355.3528,0.454);
   gre->SetPointError(2,0,0.136);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(20);
   entry=leg->AddEntry("Graph_Graph4","#varUpsilon(1S) 2010","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

  TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph_Graph5",100,35.10496,384.4663);
   // Graph_Graph5->SetMinimum(0.2444);
   // Graph_Graph5->SetMaximum(1.1276);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);
   Graph_Graph5->SetMarkerStyle(20);
   Graph_Graph5->SetMarkerSize(0.8);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.048);
   Graph_Graph5->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetTitleSize(0.048);
   Graph_Graph5->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetTitleSize(0.048);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph5);
   
   gre->Draw("pe");
   
   gre = new TGraphErrors(3);
   gre->SetName("Graph_Graph5");
   gre->SetTitle("Graph_Graph5");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,64.2184,0.677);
   gre->SetPointError(0,0,0.154);
   gre->SetPoint(1,261.4178,0.842);
   gre->SetPointError(1,0,0.212);
   gre->SetPoint(2,355.3528,0.454);
   gre->SetPointError(2,0,0.136);
   

   // TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph_Graph6",100,35.10496,384.4663);
   // // Graph_Graph6->SetMinimum(0.2444);
   // // Graph_Graph6->SetMaximum(1.1276);
   // Graph_Graph6->SetDirectory(0);
   // Graph_Graph6->SetStats(0);
   // Graph_Graph6->SetMarkerStyle(20);
   // Graph_Graph6->SetMarkerSize(0.8);
   // Graph_Graph6->GetXaxis()->SetLabelFont(42);
   // Graph_Graph6->GetXaxis()->SetTitleSize(0.048);
   // Graph_Graph6->GetXaxis()->SetTitleOffset(1.15);
   // Graph_Graph6->GetXaxis()->SetTitleFont(42);
   // Graph_Graph6->GetYaxis()->SetLabelFont(42);
   // Graph_Graph6->GetYaxis()->SetTitleSize(0.048);
   // Graph_Graph6->GetYaxis()->SetTitleOffset(1.2);
   // Graph_Graph6->GetYaxis()->SetTitleFont(42);
   // Graph_Graph6->GetZaxis()->SetLabelFont(42);
   // Graph_Graph6->GetZaxis()->SetTitleSize(0.048);
   // Graph_Graph6->GetZaxis()->SetTitleFont(42);
   // gre->SetHistogram(Graph_Graph6);
   
   // gre->Draw("p");

   // box = new TBox(370,0.94,385,1.06);

   // ci = TColor::GetColor("#ff99ff");
   // box->SetFillColor(ci);
   // box->Draw();
     tex = new TLatex(240,1.1,"|y| < 2.4");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(340,0.62,"0-10%");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(240,0.58,"10-20%");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(40,0.85,"20-100%");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   
*/
void plotDoubleRatios()
{
  float pt2014 [nPtBins_2014] = {1.25, 3.75, 6.5, 10., 16.,30};
  float pt2014e[nPtBins_2014] = {1.25, 1.25, 1.5, 2., 4.,10.};
  float doubleRatio3S1S[nPtBins_2014];
  float doubleRatio3S1Se[nPtBins_2014];
  float doubleRatio2S1S[nPtBins_2014];
  float doubleRatio2S1Se[nPtBins_2014];
  float ppSingleRatio3S1S[nPtBins_2014];
  float ppSingleRatio3S1Se[nPtBins_2014];
  float ppSingleRatio2S1S[nPtBins_2014];
  float ppSingleRatio2S1Se[nPtBins_2014];
  float paSingleRatio3S1S[nPtBins_2014];
  float paSingleRatio3S1Se[nPtBins_2014];
  float paSingleRatio2S1S[nPtBins_2014];
  float paSingleRatio2S1Se[nPtBins_2014];
    
  for(int i=0 ; i<nPtBins_2014 ; i++){
    ppSingleRatio2S1S[i]=computeRatio(N2S_pp_m146p239pt3p5[i],N1S_pp_m146p239pt3p5[i]);
    paSingleRatio2S1S[i]=computeRatio(N2S_pa_pt3p5[i],N1S_pa_pt3p5[i]);
    ppSingleRatio3S1S[i]=computeRatio(N3S_pp_m146p239pt3p5[i],N1S_pp_m146p239pt3p5[i]);
    paSingleRatio3S1S[i]=computeRatio(N3S_pa_pt3p5[i],N1S_pa_pt3p5[i]);
   
    ppSingleRatio2S1Se[i]=    computeRatioError(N2S_pp_m146p239pt3p5[i],N1S_pp_m146p239pt3p5[i],N2S_pp_m146p239pt3p5e[i],N1S_pp_m146p239pt3p5e[i]);
    paSingleRatio2S1Se[i]=    computeRatioError(N2S_pa_pt3p5[i],N1S_pa_pt3p5[i],N2S_pa_pt3p5e[i],N1S_pa_pt3p5e[i]);
    ppSingleRatio3S1Se[i]=    computeRatioError(N3S_pp_m146p239pt3p5[i],N1S_pp_m146p239pt3p5[i],N3S_pp_m146p239pt3p5e[i],N1S_pp_m146p239pt3p5e[i]);
    paSingleRatio3S1Se[i]=    computeRatioError(N3S_pa_pt3p5[i],N1S_pa_pt3p5[i],N3S_pa_pt3p5e[i],N1S_pa_pt3p5e[i]);
 
    doubleRatio2S1S[i]= computeRatio(paSingleRatio2S1S[i],ppSingleRatio2S1S[i]);
    doubleRatio3S1S[i]= computeRatio(paSingleRatio3S1S[i],ppSingleRatio3S1S[i]);
    doubleRatio2S1Se[i]= computeRatioError(paSingleRatio2S1S[i],ppSingleRatio2S1S[i],paSingleRatio2S1Se[i],ppSingleRatio2S1Se[i]);
    doubleRatio3S1Se[i]= computeRatioError(paSingleRatio3S1S[i],ppSingleRatio3S1S[i],paSingleRatio3S1Se[i],ppSingleRatio3S1Se[i]);

    cout <<"pA/pp 2s/1s = "<< doubleRatio2S1S[i] <<" +/- "<< doubleRatio2S1Se[i]<< endl;
    cout <<"pA/pp 3s/1s = "<< doubleRatio3S1S[i] <<" +/- "<< doubleRatio3S1Se[i]<< endl;
  }
  TCanvas *cDRpt = new TCanvas("cDRpt","cDRpt"); 
  cDRpt->cd();
  TPad *ppt2 = new TPad("ppt2","ppt2",0.0,0.0,1.0,1.0);
  ppt2->SetBottomMargin(0.12);
  ppt2->SetTopMargin(0.03);
  ppt2->SetRightMargin(0.03);
  ppt2->SetLeftMargin(0.16);
  ppt2->Draw();
  ppt2->cd();
  //one pad to draw RaaPt!
  TF1 *f4DRPt = new TF1("f4DRPt","1",0,40.5);
  f4DRPt->SetLineWidth(0);
  f4DRPt->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  f4DRPt->GetYaxis()->SetTitle("DOUBLE RATIOS");
  f4DRPt->GetYaxis()->SetTitleOffset(1.8);
  f4DRPt->GetYaxis()->SetTitleSize(0.028);
  f4DRPt->GetYaxis()->SetRangeUser(0.,1.4);
  f4DRPt->GetXaxis()->CenterTitle(kTRUE);
  f4DRPt->Draw();
  TGraphErrors *gDRPt2 = new TGraphErrors(nPtBins_2014,pt2014,doubleRatio2S1S,pt2014e,doubleRatio2S1Se);
  gDRPt2->SetMarkerColor(kBlue+7);
  gDRPt2->SetMarkerStyle(33);
  gDRPt2->SetMarkerSize(2);
  TGraphErrors *gDRPt2circle = new TGraphErrors(nPtBins_2014,pt2014,doubleRatio2S1S,pt2014e,doubleRatio2S1Se);
  gDRPt2circle->SetMarkerStyle(27);
  gDRPt2circle->SetMarkerSize(2);
  gDRPt2circle->SetLineColor(kBlack);
  gDRPt2->Draw("pe");
  gDRPt2circle->Draw("p");
  f4DRPt->Draw("same");
  TGraphErrors *gDRPt3 = new TGraphErrors(nPtBins_2014,pt2014,doubleRatio3S1S,pt2014e,doubleRatio3S1Se);
  gDRPt3->SetMarkerColor(kBlue+2);
  gDRPt3->SetMarkerStyle(33);
  gDRPt3->SetMarkerSize(2);
  TGraphErrors *gDRPt3circle = new TGraphErrors(nPtBins_2014,pt2014,doubleRatio3S1S,pt2014e,doubleRatio3S1Se);
  gDRPt3circle->SetMarkerStyle(27);
  gDRPt3circle->SetMarkerSize(2);
  gDRPt3circle->SetLineColor(kBlack);
  gDRPt3->Draw("pe");
  gDRPt3circle->Draw("p");
  f4DRPt->Draw("same");
  TLatex *l1CMSpt = new TLatex(2,1.32, "CMS - Work in progress");
  l1CMSpt->SetTextFont(42);
  l1CMSpt->SetTextSize(0.032);
  l1CMSpt->Draw();
  l1CMSpt->DrawLatex(2,1.25,"#sqrt{s} = 2.76 TeV, #sqrt{s_{NN}} = 5.02 TeV");
  l1CMSpt->DrawLatex(2,1.18,"-1.47 < |y_{C.M.}| < 2.4");
  TLatex *lyLPT= new TLatex(2,1.1,"2013 L_{int}^{pPb} = 18.4 nb^{-1}, 2013 L_{int}^{pp} = 5.4 pb^{-1}"); // 
  lyLPT->SetTextFont(42);
  lyLPT->SetTextSize(0.027);
  lyLPT->Draw();
  gPad->RedrawAxis();

  TLegend *legend = new TLegend(0.6,0.21,0.8,0.4);
  legend->SetTextSize(0.029);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->AddEntry(gDRPt2,"#varUpsilon(2S)/#varUpsilon(1S) pPb/pp","lp");
  legend->AddEntry(gDRPt3,"#varUpsilon(3S)/#varUpsilon(1S) pPb/pp","lp");
  legend->Draw();

 gPad->RedrawAxis();
 cDRpt->SaveAs("~/Project/ups2013/code/pdfOutput/DoubleRatiosPA.pdf");
 /////////////// pPb vs,RATIOS
 
 TCanvas *cSRpt = new TCanvas("cSRpt","cSRpt"); 
 cSRpt->cd();
 TPad *ppt2 = new TPad("ppt2","ppt2",0.0,0.0,1.0,1.0);
 ppt2->SetBottomMargin(0.12);
 ppt2->SetTopMargin(0.03);
 ppt2->SetRightMargin(0.03);
 ppt2->SetLeftMargin(0.16);
 ppt2->Draw();
 ppt2->cd();

  //one pad to draw RaaPt!
  TF1 *f4SRPt = new TF1("f4SRPt","1",0,40.5);
  f4SRPt->SetLineWidth(0);
  f4SRPt->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  f4SRPt->GetYaxis()->SetTitle("SINGLE ratios");
  f4SRPt->GetYaxis()->SetTitleOffset(1.8);
  f4SRPt->GetYaxis()->SetTitleSize(0.028);
  f4SRPt->GetYaxis()->SetRangeUser(0.,1.);
  f4SRPt->GetXaxis()->CenterTitle(kTRUE);
  f4SRPt->Draw();
  TGraphErrors *gSRPt2pa = new TGraphErrors(nPtBins_2014,pt2014,paSingleRatio2S1S,pt2014e,paSingleRatio2S1Se);
  gSRPt2pa->SetMarkerColor(kBlue);
  gSRPt2pa->SetMarkerStyle(20);
  gSRPt2pa->SetMarkerSize(1.2);
  TGraphErrors *gSRPt2pacircle = new TGraphErrors(nPtBins_2014,pt2014,paSingleRatio2S1S,pt2014e,paSingleRatio2S1Se);
  gSRPt2pacircle->SetMarkerStyle(24);
  gSRPt2pacircle->SetMarkerSize(1.2);
  gSRPt2pacircle->SetLineColor(kBlack);
  gSRPt2pa->Draw("pe");
  gSRPt2pacircle->Draw("p");
  f4SRPt->Draw("same");
  TGraphErrors *gSRPt3pa = new TGraphErrors(nPtBins_2014,pt2014,paSingleRatio3S1S,pt2014e,paSingleRatio3S1Se);
  gSRPt3pa->SetMarkerColor(kGreen);
  gSRPt3pa->SetMarkerStyle(21);
  gSRPt3pa->SetMarkerSize(1.2);
  TGraphErrors *gSRPt3pacircle = new TGraphErrors(nPtBins_2014,pt2014,paSingleRatio3S1S,pt2014e,paSingleRatio3S1Se);
  gSRPt3pacircle->SetMarkerStyle(25);
  gSRPt3pacircle->SetMarkerSize(1.2);
  gSRPt3pacircle->SetLineColor(kBlack);
  gSRPt3pa->Draw("pe");
  gSRPt3pacircle->Draw("p");
  f4SRPt->Draw("same");


  TGraphErrors *gSRPt2pp = new TGraphErrors(nPtBins_2014,pt2014,ppSingleRatio2S1S,pt2014e,ppSingleRatio2S1Se);
  gSRPt2pp->SetMarkerColor(kViolet);
  gSRPt2pp->SetMarkerStyle(20);
  gSRPt2pp->SetMarkerSize(1.2);
  TGraphErrors *gSRPt2ppcircle = new TGraphErrors(nPtBins_2014,pt2014,ppSingleRatio2S1S,pt2014e,ppSingleRatio2S1Se);
  gSRPt2ppcircle->SetMarkerStyle(24);
  gSRPt2ppcircle->SetMarkerSize(1.2);
  gSRPt2ppcircle->SetLineColor(kBlack);
  gSRPt2pp->Draw("pe");
  gSRPt2ppcircle->Draw("p");
  f4SRPt->Draw("same");
  TGraphErrors *gSRPt3pp = new TGraphErrors(nPtBins_2014,pt2014,ppSingleRatio3S1S,pt2014e,ppSingleRatio3S1Se);
  gSRPt3pp->SetMarkerColor(kOrange+1);
  gSRPt3pp->SetMarkerStyle(21);
  gSRPt3pp->SetMarkerSize(1.2);
  TGraphErrors *gSRPt3ppcircle = new TGraphErrors(nPtBins_2014,pt2014,ppSingleRatio3S1S,pt2014e,ppSingleRatio3S1Se);
  gSRPt3ppcircle->SetMarkerStyle(25);
  gSRPt3ppcircle->SetMarkerSize(1.2);
  gSRPt3ppcircle->SetLineColor(kBlack);
  gSRPt3pp->Draw("pe");
  gSRPt3ppcircle->Draw("p");
  f4SRPt->Draw("same");
  TLatex *l1CMSpt = new TLatex(2,0.92, "CMS - Work in progress");
  l1CMSpt->SetTextFont(42);
  l1CMSpt->SetTextSize(0.032);
  l1CMSpt->Draw();
  l1CMSpt->DrawLatex(2,0.87,"#sqrt{s} = 2.76 TeV, #sqrt{s_{NN}} = 5.02 TeV");
  l1CMSpt->DrawLatex(2,0.81,"-1.47 < |y_{C.M.}| < 2.4");
  TLatex *lyLPT= new TLatex(2,0.75,"2013 L_{int}^{pPb} = 18.4 nb^{-1}, 2013 L_{int}^{pp} = 5.4 pb^{-1}"); // 
  lyLPT->SetTextFont(42);
  lyLPT->SetTextSize(0.027);
  lyLPT->Draw();
  gPad->RedrawAxis();

  TLegend *legend = new TLegend(0.2,0.5,0.5,0.7);
  legend->SetTextSize(0.029);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->AddEntry(gSRPt2pp,"#varUpsilon(2S)/#varUpsilon(1S) pp","lp");
  legend->AddEntry(gSRPt2pa,"#varUpsilon(2S)/#varUpsilon(1S) pPb ","lp");
  legend->AddEntry(gSRPt3pp,"#varUpsilon(3S)/#varUpsilon(1S) pp","lp");
  legend->AddEntry(gSRPt3pa,"#varUpsilon(3S)/#varUpsilon(1S) pPb","lp");
  legend->Draw();

 gPad->RedrawAxis();
cSRpt->SaveAs("~/Project/ups2013/code/pdfOutput/SingleRatiosPPandPA.pdf");
}

void combine_blue(double val1, double err1, double val2, double err2)
{
   double w1 = err2*err2/(err1*err1+err2*err2);
   double w2 = err1*err1/(err1*err1+err2*err2);

   cout << w1*val1+w2*val2 << "\\pm" << err1*err2/sqrt(err1*err1+err2*err2) << endl;
}
