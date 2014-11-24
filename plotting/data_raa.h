
  //integers for binning.
  //centrality binning in 11-011

#define nfitvars 5
#define nbkgdvars 2
#define nPtBins_2013  5
#define nPtBins_2014 6
#define nPtBins_2010  3
#define nRapBins_2013 5
#define nRapBins_2010 2
#define nRapBins_2014 6
#define nCentBins_2014 8
#define nRapBins_pPb 9
#define bin1  7
#define bin 8
#define L_pp_invNb 5400
#define L_pp_invNbe 199.8
#define L_pp 5400000000000
#define L_ppe 199800000000 
#define N_MB_uncorr 1.126653312	
#define N_MB_corr   1.138033648
#define T_AA_b 5662
#define T_AA_mb 5.66     

  float cent[bin1]  ={22.05, 86.3, 130.0, 187.1, 261.4, 330.4, 381.3}; // with 40-50 and 50-100 //and the 5-10 bin was 329.5, which is inconsistent with a few lines below.
float nPart2014[nCentBins_2014]  ={8.75, 42.02, 86.3, 130.0, 187.2, 261.4, 330.4, 381.3}; // with 40-50 and 50-100 //and the 5-10 bin was 329.5, which is inconsistent with a few lines below.
float taa[bin1] = {0.486,2.748,5.089,8.782,14.477,20.47,25.901}; // taa of 40-50, 50-100
float taa2014[nCentBins_2014] = {0.125,0.982,2.748,5.089,8.782,14.477,20.47,25.901}; //taa of 40-50, 50-100
float mb_percentage2014[nCentBins_2014] = {0.3,0.2,0.1,0.1,0.1,0.1,0.05,0.05};
float mb_percentage[bin1] = {0.5,0.1,0.1,0.1,0.1,0.05,0.05};
  int binPt = 5;
   int binPt2010 = 3;
   int binRap = 5;
   int binRap2010 = 2;
   int bin2010 = 4;
  
   // float cent1[8] ={14.2, 69.9, 130.0, 187.1, 261.4, 329.4, 381.3}; // with 40-60 and 60-100
   float nPart1[9] ={8.75,42.02, 86.3, 130.1, 187.3, 261.4, 330.3, 381.2}; // with 70-100, 50-70, 40-50, 30-40, 20-30, 10-20, 5-10, 0-5
   float nPart2[7] ={17.8, 69.9, 130.0, 187.1, 261.4, 355.4};//?? 

   float cent2010[5]={308.6,64.24,261.3,355.7};//0-20,20-100,10-20,0-10

   float pt [5] = {1.25, 3.75, 6.5, 10., 16.};
   float pte[5] = {1.25, 1.25, 1.5, 2., 4.};
   float pt2010 [3] = {3.25,8.25,15.};
   float pt2010e[3] = {3.25,1.75,5};

   float rap2010[2]={0.6,1.8};
   float rap2010e[2]={0.6,0.6};
   float rap2010paper[2]={0.64,1.54};
   float rap2010paperel[2]={0.64,0.34};
   float rap2010papereh[2]={0.56,0.86};
   float rap[5] = {0.2 , 0.55, 0.85, 1.25, 1.95};
   float rape[5]= {0.2, 0.15, 0.15, 0.25, 0.45};
float rap2014[6] = {0.2 , 0.6, 1.0, 1.4, 1.8,2.2}; // for the moment with bins of ∆y =0.8 except the last one which is 1.6-2.4
float rap2014e[6] = {0.2,0.2,0.2,0.2,0.2,0.2}; // for the moment with bins of ∆y =0.8 except the last one which is 1.6-2.4

float centErr[bin1]={6,6,6,6,6,6,6};
float centErr2014[nCentBins_2014]={6,6,6,6,6,6,6,6};
float centnoErr[bin1]={0,0,0,0,0,0,0};


   //CONTROL PLOTS
   // mass resolution on the 1S peak in function of pT/rapidity/centrality bins
   float massRes_AA_pT[5]={72.7,100.6,60.8,88.5,75.8};
   float massRes_AA_pTe[5]={6.2,7.4,8.0,9.4,7.8};
   float massRes_pp_pT[5]={89.7,86.7,91.7,95.1,84.1};
   float massRes_pp_pTe[5]={4.2,4.2,3.9,4.4,4.6};
   float massRes_AA_rap[5]={58.6,61.5,86.,100.3,136.0};
   float massRes_AA_rape[5]={5.9,5.9,9.4,8.7,14.8};
   float massRes_pp_rap[5]={62.8,76.1,68.6,104.4,148.2};
   float massRes_pp_rape[5]={2.3,3.4,4.4,3.,5.9};
   float massRes_AA_npart[7]={65.8,88.7,98.6,89.7,77.6,98.8,67.1};
   float massRes_AA_nparte[7]={14.2,11.2,10.4,8.2,7.7,10.6,10.1};
   // mean mass resolution in AA "MB", and in pp "pp"
   float massRes_MB[1]={83.7};
   float massRes_MBe[1]={4.0};
   float massRes_pp[1]={91.6};
   float massRes_ppe[1]={2.1};
 
   //alice shit pT>0 GeV/c 
   #define binA 2

   float centAlice[binA]={72,308}; //20-90, 0-20%
   float centAliceErr[binA]={6,6};
   float centAliceNoErr[binA]={0,0};
   float ratioUpsAlice[binA]    ={0.512,0.351}; // 1S with 'alice' binning
   float ratioUpsAliceStat[binA]={0.036,0.026}; // stat. err.
   float ratioUpsAliceSyst[binA]={0.028,0.029}; // syst. err.

   float ratioUpsAlice2[binA]    ={0.203,0.056}; // 2S with 'alice' binning
   float ratioUpsAliceStat2[binA]={0.055,0.036};
   float ratioUpsAliceSyst2[binA]={0.035,0.050};
 

   //Raa with regit and pp2013

   float ratio1S[bin]       ={1.275,0.699,0.567,0.585,0.418,0.409,0.308,0.387};// with 70-100, 50-70, 40-50, 30-40, 20-30, 10-20, 5-10, 0-5
   float ratio1SstatErr[bin]={0.279,0.081,0.085,0.067,0.045,0.038,0.031,0.049};
   float ratio1SsystErr[bin]={0.186,0.154,0.084,0.037,0.057,0.041,0.050,0.035};
   // float ratio1S[bin1]       ={1.061,0.587,0.585,0.418,0.409,0.308,0.387};// with 60-100, 40-60, 30-40, 20-30, 10-20, 5-10, 0-5 //
   // float ratio1SstatErr[bin1]={0.175,0.079,0.067,0.045,0.038,0.031,0.049}; //
   // float ratio1SsystErr[bin1]={0.   ,0.084,0.037,0.057,0.041,0.050,0.035}; //
   float ratio2S[bin]        ={0.532,0.225,0.353,0.367,0.144,0.130,0.00,0.149};
   float ratio2SstatErr[bin] ={0.333,0.134,0.161,0.109,0.066,0.059,0.00,0.073};
   float ratio2SsystErr[bin] ={  0.0,0.135,0.252,0.043,0.070,0.042,0.038,0.032};
   //  float ratio2S[bin]        ={0.458,0.318,0.367,0.144,0.130,0.00,0.149};// with 60-100, 40-60, 30-40, 20-30, 10-20, 5-10, 0-5 //
   //  float ratio2SstatErr[bin] ={0.280,0.119,0.109,0.066,0.059,0.00,0.073};
   //  float ratio2SsystErr[bin] ={  0.0,0.135,0.043,0.070,0.042,0.038,0.032};

   float ratioJpsi[bin]={0.610, 0.661, 0.506, 0.372, 0.220, 0.202};
   float ratioJpsistatErr[bin]={0.122, 0.131, 0.085, 0.058, 0.037, 0.030};
   float ratioJpsisystErr[bin]={0.097, 0.077, 0.048, 0.029, 0.016, 0.014};

   float tot1SErr[bin];
   float tot2SErr[bin];

   /**/// p_T-BINNED DATA
 
   //------[GeV/c]-----------------0/2.5-2.5/5--5/8--8/12--12/20//	   
   float raaPt1  [nPtBins_2013] = {0.267,0.403,0.385,0.440,0.444};
   float raaPt1e [nPtBins_2013] = {0.022,0.039,0.041,0.053,0.057};
   float raaPt1s [nPtBins_2013] = {0,0,0,0,0};
   //------[GeV/c]-----------------0/6.5-6.5/10-10/20//
   float raaPt2  [nPtBins_2010] = {0.109,0.101,0.143};
   float raaPt2e [nPtBins_2010] = {0.043,0.070,0.068};
   float raaPt2s [nPtBins_2010] = {0,0,0};
   //-------[max.y]-----------------0.4----0.7---1.----1.5---2.4-//
   float raaRap1 [nRapBins_2013] = {0.362,0.384,0.417,0.497,0.414};
   float raaRap1e[nRapBins_2013] = {0.038,0.043,0.055,0.042,0.038};
   float raaRap1s[nRapBins_2013] = {0.039,0.041,0.055,0.026,0.062};
   //-------[max.y]-----------------0.4----0.7---1.----1.5---2.4-//
   float raaRap2 [nRapBins_2013] = {0.071,0.128,0.133,0.09,0.247};
   float raaRap2e[nRapBins_2013] = {0.053,0.060,0.076,0.056,0.069};
   float raaRap2s[nRapBins_2013] = {0.032,0.033,0.038,0.082,0.092};

   float raaPt2010[nPtBins_2010]= {0.43,0.88,1.72}; // 10-006 published this
   float raaPt2010e[nPtBins_2010]={0.1,0.37,0.74};
   float raaPt2010s[nPtBins_2010]={0.07,0.14,0.25};

float raaRap2010 [nRapBins_2010]={0.52,0.83};
float raaRap2010e[nRapBins_2010]={0.12,0.24};
float raaRap2010s[nRapBins_2010]={0.08,0.13};

   //1S data: pp, then PbPb.
   //------[GeV/c]--------------------0/2.5---2.5/5----5/8----8/12----12/20//	   
   float dndyPt1pp [nPtBins_2013] = {0.218  ,0.152  ,0.158  ,0.107  ,0.0583};
   float dndyPt1ppe[nPtBins_2013] = {0.00975,0.00683,0.00667,0.00535,0.00383};
   float dndyPt1pps[nPtBins_2013] = {0,0,0,0,0};
   //------[GeV/c]--------------------0/2.5---2.5/5----5/8----8/12----12/20//	
   float dndyPt1   [nPtBins_2013] = {0.0583 ,0.0615 ,0.0611 ,0.0475 ,0.0259};
   float dndyPt1e  [nPtBins_2013] = {0.00403,0.00541,0.00613,0.00526,0.0029};
   float dndyPt1s  [nPtBins_2013] = {0,0,0,0,0};
   float dndyPt2010[nPtBins_2010] = {0.293,0.093,0.066}; //10-006 published
   float dndyPt2010e[nPtBins_2010]= {0.057,0.028,0.016};
   float dndyPt2010s[nPtBins_2010]= {0.053,0.017,0.012};
   //must make 2S,3S pp yields in both binnings.
   //dndyPt2pp_2013 refers to the *1S* binning of 2013 for *2S* data
   //dndyPt2pp_2010 refers to the *1S* binning of 2010 for *2S* data

/////2S data: pp, in both binnings, then PbPb in the appropriate one.
   //------[GeV/c]------------------------0/2.5---2.5/5----5/8----8/12----12/20//	
   float dndyPt2pp_2013 [nPtBins_2013] = {0.0762 ,0.0454 ,0.0527 ,0.041  ,0.0219};
   float dndyPt2ppe_2013[nPtBins_2013] = {0.00676,0.0044 ,0.00436,0.00366,0.00252};
   float dndyPt2pps_2013[nPtBins_2013] = {0,0,0,0,0};
   //------[GeV/c]------------------------0/6.5---6.5/10--10/20---//	
   float dndyPt2pp_2010 [nPtBins_2010] = {0.15   ,0.054  ,0.0345 };
   float dndyPt2ppe_2010[nPtBins_2010] = {0.00924,0.00425,0.00324};
   float dndyPt2pps_2010[nPtBins_2010] = {0,0,0};
   float dndyPt2        [nPtBins_2010] = {0.0163 ,0.0055 ,0.00497};
   float dndyPt2e       [nPtBins_2010] = {0.00651,0.00378,0.00231};
   float dndyPt2s       [nPtBins_2010] = {0,0,0};

/////3S data: pp, in 2013 *1S* binning
   //------[GeV/c]------------------------0/2.5---2.5/5----5/8----8/12----12/20//	
   float dndyPt3pp_2013 [nPtBins_2013] = {0.033  ,0.0268 ,0.0262 ,0.0204 ,0.0142 };
   float dndyPt3ppe_2013[nPtBins_2013] = {0.00547,0.00386,0.0036 ,0.00298,0.00218};
   float dndyPt3pps_2013[nPtBins_2013] = {0,0,0,0,0};

   /**/// RAPIDITY-BINNED DATA
   float dndyRap1[5]={0.0545,0.04334,0.04645,0.09095,0.0769};
   float dndyRap1e[5]={0.00539,0.00439,0.00574,0.00677,0.00639};
   float dndyRap2[5]={0,0,0,0,0};
   float dndyRap2e[5]={0,0,0,0,0};

  
   //ALice's shit
   float rapAlice[nRapBins_2010]={2.85,3.6};
   float rapAlicee[nRapBins_2010]={0.35,0.4};
   float raaRapAlice [nRapBins_2010]={0.429,0.457};
   float raaRapAlicee[nRapBins_2010]={0.084,0.106};
   float raaRapAlicesup[nRapBins_2010]={0.052,0.100};
   float raaRapAlicesdown[nRapBins_2010]={0.068,0.128};
   float raaRapAliceglob[nRapBins_2010]={0.094,0.1};
   //

   float dndyNP2010[5]={0.468,0.517,0.643,0.347};
   float dndyNP2010e[5]={0.081,0.101,0.144,0.096};
   float dndyNP2010s[5]={0.094,0.101,0.118,0.069};
   
   float dndyNP1[9]={0.989,0.542,0.440,0.454,0.324,0.317,0.238,0.300};
   float dndyNP1e[9]={0.215,0.0625,0.066,0.0515,0.0346,0.0288,0.0238,0.0382};
   float dndyNP1s[9]={0,0,0,0,0,0,0,0};
   float dndyNP2[5]={0.235,0.067,0.0365,0.0183};
   float dndyNP2e[5]={0.0102,0.0147,0.0165,0.0128};
   float dndyNP2s[5]={0,0,0,0};
   

 float dndyRap1pp[5]={0.1533,0.115,0.113,0.186,0.189};
 float dndyRap1ppe[5]={0.0063,0.00569,0.00574,0.00767,0.00776};
 // float dndyRap1pps[5]={};
 float dndyRap2pp[5]={0.0521,0.0408,0.0406,0.0683,0.066};
 float dndyRap2ppe[5]={0.00414,0.0037,0.00394,0.00487,0.00569};
 // float dndyRap2pps[5]={};
 float dndyRap3pp[5]={0.0265,0.0178,0.0247,0.032,0.0383};
 float dndyRap3ppe[5]={0.00333,0.00297,0.00344,0.00389,0.00528};
 // float dndyRap1pps[5]={};
 float dndyRap2010[3]={0.495,0.498};
 float dndyRap2010e[3]={0.091,0.097};
 float dndyRap2010s[3]={0.091,0.092};

 
 // MB stuff !?
  //MB
  
  //pt > 4
  float centMB[1]={0.5};
  float centMB2S[1]={0.3};
  float centMBA[1]={0.4};
  float centMB3S[1]={0.7};
  float cent_limit_err[1]={0.1};
  // float raaMB1S[1]={0.564};
  // float raaMB1SstatErr[1]={0.077};
  // float raaMB1SsystErr[1]={0.071};
  // float raaMB2S[1]={0.120};
  // float raaMB2SstatErr[1]={0.038};
  // float raaMB2SsystErr[1]={0.020};
  // float raaMB3S[1]={0.033};
  // float raaMB3SstatErr[1]={0.035};//{0.035};
  // float raaMB3SsystErr[1]={0.006};
  
  float raaMB1S[1]={0.511};
  float raaMB1SstatErr[1]={0.026};
  float raaMB1SsystErr[1]={0.016};
  float raaMB2S[1]={0.126};
  float raaMB2SstatErr[1]={0.03};
  float raaMB2SsystErr[1]={0.021};
  float raaMB3S[1]={0.033};
//  float raaMB3SstatErr[1]={0.035};//{0.035};
//  float raaMB3SsystErr[1]={0.006};
  
  /// MB result!
  
  float raaMB1 [1]={0.388};
  float raaMB1e[1]={0.018};
  float raaMB1s[1]={0.042};
  
  float raaMB2 [1]={0.106};
  float raaMB2e[1]={0.026};
  float raaMB2s[1]={0.018};
  
  float raaMB3 [1]={0.046};
  float raaMB3e[1]={0.045};
  float raaMB3s[1]={0.028};
  
  float raaMB3Slimit[1]={0.1};
  float raaMBJpsi[1]={0.304};
  float raaMBJpsistatErr[1]={0.025};
  float raaMBJpsisystErr[1]={0.022};
  // float raaMBAlice[1]={0.439};
  // float raaMBAliceStat[1]={0.065};
  // float raaMBAliceSyst[1]={0.028};
  
  float raaMB1StotErr =sqrt(pow(raaMB1e[0],2)+pow(raaMB1s[0],2));
  float raaMB2StotErr =sqrt(pow(raaMB2e[0],2)+pow(raaMB2s[0],2));
  float raaMB3StotErr =sqrt(pow(raaMB3e[0],2)+pow(raaMB3s[0],2));
  
  float ppy[1]={1};
  float ppx[1]={377};
  float ppxEr[1]={0};
  float ppyErSystematic=sqrt(pow(0.06,2)+pow(0.023,2));
  float ppyErSyst[1]={ppyErSystematic};
  float ppyErtol[1]={0}; 

  float Y2Sppx[1]={392};
  float Y2SppyErSystematic=sqrt(pow(0.06,2)+pow(0.033,2));
  float Y2SppyErSyst[1]={Y2SppyErSystematic};
  float Y2SppyErtotal=sqrt(pow(Y2SppyErSystematic,2)+pow(0.202,2));
  float Y2SppyErtol[1]={0};
  float centMBErr[1]={0.15};

// doubledifferential Npart/rapidity
 // !!!! first entry is |y|<1, second entry is 1<|y|<2.4 !!!!
 float raaRN1_c [2] = {0.364,0.407}; // central values
 float raaRN1e_c[2] = {0.036,0.047};
 float raaRN1s_c[2] = {0.033,0.029};

 float raaRN2_c [2] = {0.073,0.127};
 float raaRN2e_c[2] = {0.048,0.061};
 float raaRN2s_c[2] = {0.056,0.058};

 float raaRN3_c [2] = {-0.04,0.108};
 float raaRN3e_c[2] = {0.   ,0.105};
 float raaRN3s_c[2] = {0.   ,0.035};

 float raaRN1_p [2] = {0.423,0.596}; // non-central (20-100) values
 float raaRN1e_p[2] = {0.035,0.047};
 float raaRN1s_p[2] = {0.102,0.033};

 float raaRN2_p [2] = {0.168,0.296};
 float raaRN2e_p[2] = {0.057,0.076};
 float raaRN2s_p[2] = {0.160,0.027};

 float raaRN3_p [2] = {0.145,0.152};
 float raaRN3e_p[2] = {0.107,0.129};
 float raaRN3s_p[2] = {0.064,0.024};

 float rap_short_c[2] = {0.45,1.65};
 float rap_short_p[2] = {0.55,1.75};
 float rap_short_3[2] = {0.58,1.78};
 float rap_shorte[2]= {0.55,0.65};



////EFFICIENCY and acceptance 20th Nov 2014
//A. Pythia sample. 1S //table 6
////JUNE 18th
// total acc*eff 1S = 0.449 +/- 0.0007
float Ae_1S_pythia_pt[nPtBins_2013] = {0.296,0.190,0.185,0.259,0.376};
float Ae_1S_pythia_pte[nPtBins_2013] = {0.001,0.001,0.001,0.002,0.004};
float Ae_1S_pythia_rap2014[nRapBins_2014]={0.267,0.268,0.270,0.255,0.203,0.072};
float Ae_1S_pythia_rap2014e[nRapBins_2014]={0.0011,0.001,0.0012,0.0011,0.0010,0.0006};
 //B. 2S pythia with binning of 2S. with raa_pt-y 2S in mind. //  //tables 16 and 17 of note apr20 /// should be wrong because of the trigger.
///// July 15th - Prashant
float Ae_2S_pythia_pt2010[nPtBins_2010] = {0.199,0.183,0.310};
float Ae_2S_pythia_pt2010e[nPtBins_2010] = {0.0005,0.0012,0.0027};
float Ae_2S_pythia_rap2010[nRapBins_2010] = {0.237,0.161};
float Ae_2S_pythia_rap2010e[nRapBins_2010] = {0.0007,0.0006};
/* float Ae_2S_pythia_pt2010[nPtBins_2010] = {0.231,0.229,0.368}; */
/* float Ae_2S_pythia_pt2010e[nPtBins_2010] = {0.0004,0.0010,0.0022}; */
/* float Ae_2S_pythia_rap2010[nRapBins_2010] = {0.247,0.224}; */
/* float Ae_2S_pythia_rap2010e[nRapBins_2010] = {0.0005,0.0006}; */

//B.2 2S pythia with binning even, and binning of 1S. //tables 16 and 17 of note may 6th
// total acc*eff 2S = 0.375 +/- 0.006
///// July 15th - Prashant
float Ae_2S_pythia_pt2013[nPtBins_2013]={0.269,0.157,0.160,0.229,0.351};
float Ae_2S_pythia_pt2013e[nPtBins_2013]={0.0010,0.006,0.0008,0.0018,0.0043};
float Ae_2S_pythia_rap2014[nRapBins_2014]={0.237,0.239,0.236,0.219,0.175,0.060};//bug was here.
float Ae_2S_pythia_rap2014e[nRapBins_2014]={0.0012,0.0012,0.0012,0.0012,0.0011,0.0005};
/* float Ae_2S_pythia_pt2013[nPtBins_2013]={0.469,0.299,0.323,0.423,0.559}; */
/* float Ae_2S_pythia_pt2013e[nPtBins_2013]={0.0013,0.009,0.0013,0.0024,0.0047}; */
/* float Ae_2S_pythia_rap2014[nRapBins_2014]={0.421,0.392,0.377,0.361,0.347,0.261};//bug was here. */
/* float Ae_2S_pythia_rap2014e[nRapBins_2014]={0.0015,0.0014,0.0014,0.0014,0.0016,0.0021}; */

//B.3 3S pythia with binning even, and binning of 1S. //tables 16 and 17 of note may 6th
//total 3S acc*eff = 0.750 +/- 0.0013
///// July 15th, Prashant
float Ae_3S_pythia_pt2013[nPtBins_2013]={0.277,0.160,0.165,0.238,0.365};
float Ae_3S_pythia_pt2013e[nPtBins_2013]={0.0012,0.0006,0.0008,0.0015,0.0035};
float Ae_3S_pythia_rap2014[nRapBins_2014]={0.248,0.248,0.243,0.227,0.181,0.0062};//,0.1825};//bug is here
float Ae_3S_pythia_rap2014e[nRapBins_2014]={0.0011,0.0011,0.0012,0.0011,0.0010,0.0005};//,0.005149};
/* float Ae_3S_pythia_pt2013[nPtBins_2013]={0.538,0.351,0.357,0.455,0.588}; */
/* float Ae_3S_pythia_pt2013e[nPtBins_2013]={0.0018,0.001,0.0012,0.0018,0.0028}; */
/* float Ae_3S_pythia_rap2014[nRapBins_2014]={0.475,0.450,0.436,0.421,0.394,0.280};//,0.1825};//bug is here */
/* float Ae_3S_pythia_rap2014e[nRapBins_2014]={0.0015,0.0014,0.0014,0.0015,0.0017,0.0022};//,0.005149}; */

//C. Pyquen Sample 1S. //table 18 // updated June13

 // last bin is  good today.  but values are screwed
//float Ae_1S_pyquen_pt[nPtBins_2013] = {0.348,0.228,0.239,0.342,0.486};
//float Ae_1S_pyquen_pte[nPtBins_2013] = {0.0031,0.0020,0.0023,0.027,0.0031}; // 
///// July 15th, Prashant
float Ae_1S_pyquen_pt[nPtBins_2013] = {0.262,0.169,0.174,0.250,0.370};
float Ae_1S_pyquen_pte[nPtBins_2013] = {0.0026,0.0017,0.0020,0.0023,0.0032}; // 
/* float Ae_1S_pyquen_rap[nRapBins_2013] = {0.227,0.252,0.291,0.323,0.299};// table 21 */
float Ae_1S_pyquen_rap2014[nRapBins_2014]={0.223,0.235,0.249,0.243,0.197,0.073};
float Ae_1S_pyquen_rap2014e[nRapBins_2014]={0.0025,0.0026,0.0028,0.0029,0.0027,0.0018};//,0.003386};
//float Ae_1S_pyquen_rap2014[nRapBins_2014]={0.236,0.267,0.314,0.332,0.324,0.220}; // bug is here.
//float Ae_1S_pyquen_rap2014e[nRapBins_2014]={0.0023,0.0026,0.0031,0.0035,0.0041,0.0049};//,0.003386};

//float Ae_1S_pyquen_rape[nRapBins_2013] = {0.0024,0.0031,0.0036,0.0031,0.0031};
/* float Ae_1S_pyquen_cent[bin1] = {0.2,0.280,0.286,0.285,0.279,0.274,0.267}; //starting with peripheral bin!//table 21  */
//float Ae_1S_pyquen_cent2014[nCentBins_2014]={0.299,0.299,0.297,0.295,0.293,0.287,0.280,0.272};//starting with peripheral bin!//table 20 from may 6 (havent changed much.)
float Ae_1S_pyquen_cent2014[nCentBins_2014]={0.222,0.221,0.218,0.217,0.216,0.212,0.207,0.200};//starting with peripheral bin!//table 20 from may 6 (havent changed much.)
/* float Ae_1S_pyquen_cente[bin1] = {0.0017,0.0027,0.0027,0.0027,0.0026,0.0036,0.0031}; */
//float Ae_1S_pyquen_cent2014e[nCentBins_2014]={0.0018,0.0020,0.0026,0.0025,0.0026,0.0025,0.0035,0.0030};
float Ae_1S_pyquen_cent2014e[nCentBins_2014]={0.0015,0.0016,0.0021,0.0021,0.0021,0.020,0.0029,0.0024};

//I assume these are with 3p5 ? //table 15
float Ae_1S_pyquen_DD020[nRapBins_2010]={0.266,0.305};
float Ae_1S_pyquen_DD20100[nRapBins_2010]={0.278,0.328};
float Ae_1S_pyquen_DD020e[nRapBins_2010]={0.002,0.0032};
float Ae_1S_pyquen_DD20100e[nRapBins_2010]={0.0017,0.0027};

float Aet_1S_pyquen_DDR020[nRapBins_2010]={0.249,0.196};
float Aet_1S_pyquen_DDR020e[nRapBins_2010]={0.003,0.003};
float Aet_1S_pyquen_DDR20100[nRapBins_2010]={0.259,0.211};
float Aet_1S_pyquen_DDR20100e[nRapBins_2010]={0.003,0.003};

float Aet_1S_pyquen_DDP020[nPtBins_2010]={0.215,0.218,0.341};
float Aet_1S_pyquen_DDP020e[nPtBins_2010]={0.003,0.004,0.005};
float Aet_1S_pyquen_DDP20100[nPtBins_2010]={0.227,0.228,0.356};
float Aet_1S_pyquen_DDP20100e[nPtBins_2010]={0.002,0.004,0.004};


//D. Pyquen sample 2S. with binning of 2S. // with 
// total 2S acc*eff = 0.254 +/- 0.0022
float Ae_2S_pyquen_pt[nPtBins_2010] = {0.185,0.182,0.3051}; //table 25
float Ae_2S_pyquen_pte[nPtBins_2010] = {0.001,0.002,0.004};
float Ae_2S_pyquen_rap[nRapBins_2010] = {0.219,0.159};  //table 26
float Ae_2S_pyquen_rape[nRapBins_2010] = {0.002,0.002};
float Ae_2S_pyquen_cent2014[bin1]={0.200,0.198,0.197,0.195,0.193,0.192,0.185};
float Ae_2S_pyquen_cent2014e[bin1] = {0.0016,0.0022,0.0022,0.0022,0.0021,0.0030,0.0025};
float Ae_1S_pyquen_rap2014[nRapBins_2014]= {0.211,0.220,0.225,0.216,0.177,0.062};
float Ae_1S_pyquen_rap2014e[nRapBins_2014]={0.003,0.003,0.003,0.003,0.0030,0.002};//,0.003386};
float Ae_2S_pyquen_rap2014[nRapBins_2014]= {0.211,0.220,0.225,0.216,0.177,0.062};
float Ae_2S_pyquen_rap2014e[nRapBins_2014]={0.003,0.003,0.003,0.003,0.003,0.002};
//float Ae_2S_pyquen_cent2014[nCentBins_2014] ={0.266,0.264,0.265,0.263,0.263,0.255,0.254,0.245};
//float Ae_2S_pyquen_cent2014e[nCentBins_2014]={0.0017,0.0018,0.0024,0.0024,0.0024,0.0023,0.0033,0.0028};

float Ae_1S_pythia_tot= 0.2310;
float Ae_1S_pythia_tote=0.0006;
float Aet_1S_pythia_tot=0.251;
float Aet_1S_pythia_tote = 0.001;

float Ae_2S_pythia_tot= 0.201;
float Ae_2S_pythia_tote=0.0006;

float Ae_3S_pythia_tot= 0.210;
float Ae_3S_pythia_tote=0.0007;

float Ae_1S_pyquen_tot=0.210;
float Ae_1S_pyquen_tote=0.0009;

float Ae_2S_pyquen_tot=0.1918;
float Ae_2S_pyquen_tote=0.0010;

///for FUN
float Ae_3S_pyquen_tot=0.2003;
float Ae_3S_pyquen_tote=0.02;

/*

0,1918/0,201=0,954
0,210*0,1918/0,201=0.2003

max. 10%uncertainty on this from varying ratio eff vs pt or rapidity
-> effPyquen3S tot = 0.2003+/-0.02

 */
//June 18th
//fiducial purposes...
/* float A_1S_pythia_pp_pt3p5[nPtBins_2013]={0.56,0.552,0.598,0.637,0.677}; */
/* float A_1S_pythia_pp_pt3p5e[nPtBins_2013]={0.0109,0.0107,0.0142,0.243}; */
//rap in bins called 2014, pt in bins called 2013
//eff tot 1S = 0.656 +/- 0.0011
// acc tot 1S =0.3520 +/- 0.0007
// tnp tot 1S = 1.085 +/- 0.002

float A_1S_pythia_pt3p5[nPtBins_2013]={0.46,0.30,0.28,0.37,0.51};
float A_1S_pythia_pt3p5e[nPtBins_2013]={0.0001,0.0001,0.0001,0.001,0.002};
float e_1S_pythia_pt3p5[nPtBins_2013]={0.646,0.634,0.669,0.697,0.738};
float e_1S_pythia_pt3p5e[nPtBins_2013]={0.0019,0.0019,0.0027,0.0036,0.0051};
float A_1S_pythia_rap3p5[nRapBins_2014]={0.39,0.39,0.39,0.39,0.34,0.15};
float A_1S_pythia_rap3p5e[nRapBins_2014]={0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
float e_1S_pythia_rap3p5[nRapBins_2014]={0.686,0.689,0.692,0.654,0.598,0.484};
float e_1S_pythia_rap3p5e[nRapBins_2014]={0.0025,0.0025,0.0025,0.0025,0.0027,0.0038};
float t_1S_pythia_pt3p5[nPtBins_2013]={1.108,1.1,1.077,1.057,1.037};
float t_1S_pythia_pt3p5e[nPtBins_2013]={0.004,0.004,0.005,0.006,0.008};
float t_1S_pythia_rap3p5[nRapBins_2014]={1.069,1.069,1.076,1.099,1.138,1.160};
float t_1S_pythia_rap3p5e[nRapBins_2014]={0.004,0.004,0.004,0.005,0.006,0.011};
//eff tot 2S = 0.718 +/- 0.0013
//acc tot 2S =  0.28 +/- 0.000X
//tnp tot 2S = 1.067 +/- 0.002
float A_2S_pythia_pt2013[nPtBins_2013]=   {0.38,0.22,0.22,0.31,0.46};
float A_2S_pythia_pt2013e[nPtBins_2013]=  {0.0001,0.0001,0.0001,0.0001,0.0001};
float e_2S_pythia_pt2013[nPtBins_2013]=   {0.709,0.716,0.728,0.740,0.765};
float e_2S_pythia_pt2013e[nPtBins_2013]=  {0.0022,0.0024,0.0034,0.0047,0.0069};
float A_2S_pythia_rap2014[nRapBins_2014]= {0.31,0.31,0.31,0.31,0.28,0.12};
float A_2S_pythia_rap2014e[nRapBins_2014]={0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
float e_2S_pythia_rap2014[nRapBins_2014]= {0.766,0.722,0.761,0.708,0.625,0.503};
float e_2S_pythia_rap2014e[nRapBins_2014]={0.0031,0.0032,0.0032,0.0032,0.0032,0.0045};
float t_2S_pythia_pt4[nPtBins_2013]={1.08,1.072,1.057,1.045,1.031};
float t_2S_pythia_pt4e[nPtBins_2013]={0.004,0.004,0.005,0.007,0.010};
float t_2S_pythia_rap4[nRapBins_2014]={1.049,1.049,1.053,1.077,1.12,1.146};
float t_2S_pythia_rap4e[nRapBins_2014]={0.005,0.005,0.005,0.005,0.006,0.012};

float t_3S_pythia_pt4[nPtBins_2013]={1.072,1.067,1.054,1.043,1.03};
float t_3S_pythia_pt4e[nPtBins_2013]={0.004,0.004,0.005,0.007,0.010};
float t_3S_pythia_rap4[nRapBins_2014]={1.037,1.037,1.043,1.068,1.111,1.139};
float t_3S_pythia_rap4e[nRapBins_2014]={0.005,0.005,0.005,0.005,0.006,0.012};
//eff tot 3S = 0.750 +/- 0.0013
//acc tot 3S = 0.33 +/- 0.004
float A_3S_pythia_pt2013[nPtBins_2013]=   {0.448,0.260,0.259,0.366,0.542};// ooh...
float A_3S_pythia_pt2013e[nPtBins_2013]=  {0.002,0.0011,0.0013,0.002,0.0031};///oooooh
float e_3S_pythia_pt2013[nPtBins_2013]=   {0.731,0.731,0.751,0.769,0.795};
float e_3S_pythia_pt2013e[nPtBins_2013]=  {0.0026,0.0024,0.0029,0.0033,0.0041};
float A_3S_pythia_rap2014[nRapBins_2014]= {0.366,0.366,0.336,0.335,0.330,0.141}; // ooh...
float A_3S_pythia_rap2014e[nRapBins_2014]={0.0017,0.0016,0.0016,0.0017,0.002,0.0028}; // oooh
float e_3S_pythia_rap2014[nRapBins_2014]= {0.802,0.803,0.786,0.733,0.647,0.519};
float e_3S_pythia_rap2014e[nRapBins_2014]={0.0029,0.0029,0.003,0.0029,0.0031,0.0044};
//2010 bins ... :/
float A_2S_pythia_rap2010[nRapBins_2010]= {0.31,0.25};
float A_2S_pythia_rap2010e[nRapBins_2010]={0.0005,0.0008};
float e_2S_pythia_rap2010[nRapBins_2010]= {0.767,0.647};
float e_2S_pythia_rap2010e[nRapBins_2010]={0.0016,0.0017};
float t_2S_pythia_rap2010[nRapBins_2010]={1.050,1.101};
float t_2S_pythia_rap2010e[nRapBins_2010]={0.001,0.001};
float A_2S_pythia_pt2010[nPtBins_2010]= {0.28,0.25,0.41};
float A_2S_pythia_pt2010e[nPtBins_2010]={0.0002,0.0002,0.0041};
float e_2S_pythia_pt2010[nPtBins_2010]= {0.714,0.735,0.757};
float e_2S_pythia_pt2010e[nPtBins_2010]={0.0015,0.004,0.0051};
float t_2S_pythia_pt2010[nPtBins_2010] = {1.074,1.050,1.035};
float t_2S_pythia_pt2010e[nPtBins_2010] = {0.001,0.002,0.005};


//Pyquen 1S
float e_1S_pyquen_pt[nPtBins_2013]={0.571,0.565,0.623,0.677,0.726};
float e_1S_pyquen_pte[nPtBins_2013]={0.0055,0.0058,0.007,0.006,0.0049};
float A_1S_pyquen_pt[nPtBins_2013]={0.46,0.30,0.28,0.37,0.51};
float A_1S_pyquen_pte[nPtBins_2013]={0.0009,0.0001,0.0006,0.0012,0.0027};
float e_1S_pyquen_rap2014[nRapBins_2014]={0.572,0.603,0.639,0.625,0.580,0.491};
float e_1S_pyquen_rap2014e[nRapBins_2014]={0.0062,0.0066,0.0071,0.0073,0.0080,0.0121};
float A_1S_pyquen_rap2014[nRapBins_2014]={0.39,0.39,0.39,0.39,0.34,0.15};
float A_1S_pyquen_rap2014e[nRapBins_2014]={0.0009,0.0009,0.0009,0.0009,0.0009,0.0004};
float e_1S_pyquen_cent2014[nCentBins_2014]={0.634,0.632,0.625,0.621,0.618,0.606,0.594,0.572};
float e_1S_pyquen_cent2014e[nCentBins_2014]={0.0043,0.0046,0.0061,0.0060,0.0060,0.0058,0.0082,0.0070};
float A_1S_pyquen_cent2014[nCentBins_2014]={0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35};
float A_1S_pyquen_cent2014e[nCentBins_2014]={0.0043,0.0046,0.0061,0.0060,0.0060,0.0058,0.0082,0.0070}; // starting with the peripheral value!!!
float t_1S_pyquen_pt3p5[nPtBins_2013]={1.106,1.098,1.077,1.056,1.038};
float t_1S_pyquen_pt3p5e[nPtBins_2013]={0.012,0.013,0.014,0.010,0.008};
float t_1S_pyquen_rap3p5[nRapBins_2014]={1.064,1.064,1.072,1.097,1.143,1.168};
float t_1S_pyquen_rap3p5e[nRapBins_2014]={0.013,0.013,0.013,0.014,0.018,0.033};
float t_1S_pyquen_cent2014[nCentBins_2014]={1.090,1.088,1.087,1.087,1.087,1.087,1.087,1.087};
float t_1S_pyquen_cent2014e[nCentBins_2014]={0.008,0.009,0.012,0.012,0.012,0.012,0.017,0.015};
/* float e_1S_pyquen_rap[nRapBins_2010]={}; */
/* float e_1S_pyquen_rape[nRapBins_2010]={}; */
/* float A_1S_pyquen_rap[nRapBins_2010]={}; */
/* float A_1S_pyquen_rape[nRapBins_2010]={}; */


float e_1S_pyquen_DD020[nRapBins_2010]={0.596,0.577};
float e_1S_pyquen_DD20100[nRapBins_2010]={0.621,0.622};
float e_1S_pyquen_DD020e[nRapBins_2010]={0.0051,0.0066};
float e_1S_pyquen_DD20100e[nRapBins_2010]={0.0043,0.0057};
float A_1S_pyquen_DD020[nRapBins_2010]={0.447,0.530};
float A_1S_pyquen_DD20100[nRapBins_2010]={0.449,0.528};
float A_1S_pyquen_DD020e[nRapBins_2010]={0.0024,0.0039};
float A_1S_pyquen_DD20100e[nRapBins_2010]={0.0020,0.0032};
float t_1S_pyquen_DDR020[nRapBins_2010]= {1.067,1.122};
float t_1S_pyquen_DDR20100[nRapBins_2010]= {1.067,1.122};
float t_1S_pyquen_DDP020[nPtBins_2010]=   {1.099,1.066,1.043};
float t_1S_pyquen_DDP20100[nPtBins_2010]=   {1.1,1.066,1.043};

float t_2S_pyquen_pt[nPtBins_2013]={1.078,1.071,1.057,1.045,1.032};
float t_2S_pyquen_pte[nPtBins_2013]={0.012,0.014,0.015,0.011,0.007};
float t_2S_pyquen_rap[nRapBins_2014]={1.044,1.044,1.050,1.077,1.125,1.154};
float t_2S_pyquen_rape[nRapBins_2014]={0.013,0.013,0.013,0.014,0.018,0.033};
float t_2S_pyquen_cent2014[bin1]={1.067,1.067,1.067,1.0670,1.067,1.066,1.066};
float t_2S_pyquen_cent2014e[bin1]={0.010,0.013,0.013,0.013,0.013,0.018,0.016};
float e_2S_pyquen_rap2010[nRapBins_2010]= {0.710,0.647};
float e_2S_pyquen_rap2010e[nRapBins_2010]={0.005,0.006};
float A_2S_pyquen_rap2010[nRapBins_2010]= {0.343,0.424};
float A_2S_pyquen_rap2010e[nRapBins_2010]={0.0016,0.0025};
float e_2S_pyquen_pt2010[nPtBins_2010]=   {0.674,0.718,0.749};
float e_2S_pyquen_pt2010e[nPtBins_2010]=  {0.0047,0.0074,0.0049};
float A_2S_pyquen_pt2010[nPtBins_2010]=   {0.28,0.25,0.41};
float A_2S_pyquen_pt2010e[nPtBins_2010]=  {0.0004,0.0009,0.0024};
float t_2S_pyquen_rap2010[nRapBins_2010]= {1.046,1.103};
float t_2S_pyquen_rap2010e[nRapBins_2010]={0.008,0.012};
float t_2S_pyquen_pt2010[nPtBins_2010]=   {1.073,1.05,1.036};
float t_2S_pyquen_pt2010e[nPtBins_2010]=  {0.008,0.012,0.007};



//28th april 2014 
//integrated values!
/* pp Y(1S) */
/* Upsilon eff Total 0.560949 error 0.000954518 */
/* Upsilon Acc*eff Total 0.26641 error 0.000408323 */

/* pp Y(2S) */
/* Upsilon eff Total 0.640482 error 0.00123349 */
/* Upsilon Acc*eff Total 0.239197 error 0.000400378 */

/* PbPb Y(1S) */
/* Upsilon eff Total 0.587255 error 0.0073492 */
/* Upsilon Acc*eff Total 0.274313 error 0.0025024 */

/* PbPb Y(2S) */
/* Upsilon eff Total 0.686029 error 0.00928091 */
/* Upsilon Acc*eff Total 0.254253 error 0.00226777 */

///Tag and probe correction of corrections
float Aet_1S_pythia_pt[nPtBins_2013]={0.328,0.209,0.199,0.274,0.390};
float Aet_1S_pythia_pte[nPtBins_2013]={0.002,0.001,0.002,0.003,0.005};
float Aet_1S_pythia_rap2014[nRapBins_2014]={0.287,0.287,0.292,0.279,0.232,0.082};
float Aet_1S_pythia_rap2014e[nRapBins_2014]={0.002,0.002,0.002,0.002,0.002,0.001};

float Aet_2S_pythia_pt2013[nPtBins_2013]={0.295,0.168,0.167,0.236,0.363};
float Aet_2S_pythia_pt2013e[nPtBins_2013]={0.002,0.001,0.002,0.003,0.007};
float Aet_2S_pythia_rap2014[nRapBins_2014]={0.248,0.249,0.246,0.235,0.197,0.071};
float Aet_2S_pythia_rap2014e[nRapBins_2014]={0.002,0.002,0.002,0.002,0.002,0.001};

float Aet_3S_pythia_pt2013[nPtBins_2013]={0.301,0.171,0.171,0.245,0.378};
float Aet_3S_pythia_pt2013e[nPtBins_2013]={0.002,0.001,0.002,0.003,0.006};
float Aet_3S_pythia_rap2014[nRapBins_2014]={0.257,0.256,0.251,0.241,0.202,0.073};
float Aet_3S_pythia_rap2014e[nRapBins_2014]={0.002,0.002,0.002,0.002,0.002,0.001};

float Aet_1S_pyquen_pt[nPtBins_2013]={0.289,0.186,0.185,0.266,0.384};
float Aet_1S_pyquen_pte[nPtBins_2013]={0.004,0.003,0.003,0.004,0.005};
float Aet_1S_pyquen_rap2014[nRapBins_2014]={0.238,0.250,0.269,0.266,0.226,0.084};
float Aet_1S_pyquen_rap2014e[nRapBins_2014]={0.004,0.004,0.005,0.005,0.005,0.003};
float Aet_1S_pyquen_cent2014[nCentBins_2014]={0.243,0.242,0.239,0.238,0.236,0.232,0.227,0.219};
float Aet_1S_pyquen_cent2014e[nCentBins_2014]={0.003,0.003,0.004,0.003,0.004,0.003,0.005,0.004};

float Aet_2S_pythia_pt2013Large[nPtBins_2010]={0.212,0.196,0.321};
float Aet_2S_pythia_pt2013Largee[nPtBins_2010]={0.001,0.002,0.005};
float Aet_2S_pythia_rap2014Large[nRapBins_2010]={0.249,0.175};
float Aet_2S_pythia_rap2014Largee[nRapBins_2010]={0.001,0.001};

float Aet_2S_pyquen_pt2013Large[nPtBins_2010]={0.200,0.191,0.316};
float Aet_2S_pyquen_pt2013Largee[nPtBins_2010]={0.002,0.003,0.005};
float Aet_2S_pyquen_rap2014Large[nRapBins_2010]={0.230,0.177};
float Aet_2S_pyquen_rap2014Largee[nRapBins_2010]={0.003,0.003};

float Aet_2S_pyquen_pt2013[nPtBins_2013]={0.275,0.160,0.161,0.232,0.361};
float Aet_2S_pyquen_pt2013e[nPtBins_2013]={0.004,0.003,0.003,0.004,0.006};
float Aet_2S_pyquen_rap2014[nRapBins_2014]={0.220,0.230,0.236,0.233,0.199,0.071};
float Aet_2S_pyquen_rap20104e[nRapBins_2014]={0.004,0.004,0.005,0.005,0.005,0.003};
float Aet_2S_pyquen_cent2014[nCentBins_2014]={0.214,0.213,0.211,0.209,0.206,0.205,0.200,0.198};
float Aet_2S_pyquen_cent2014[nCentBins_2014]={0.003,0.004,0.004,0.004,0.003,0.005,0.004,0.004};








////yields in the note
/* float N1S_aa_pt3p5[nPtBins_2013] = {863,929,572,346,184}; // ,53 *///prior to Oct 20!
/* float N1S_aa_pt3p5e[nPtBins_2013] = {92,75,66,32,21}; //,8.9 */
/* float N1S_aa_pt3p5[nPtBins_2013] = {691,686,507,343,175}; // ,53 */
/* float N1S_aa_pt3p5e[nPtBins_2013] = {42,41,66,32,21}; //,8.9 */
float N1S_aa_pt3p5[nPtBins_2013]={725.578,700.694,515.951,349.402,182.037};
float N1S_aa_pt3p5e[nPtBins_2013]={43.222,41.8636,36.5595,40.0363,16.6537};
float N1S_aa_pt3p5Large[nPtBins_2010] = {1842,466,311}; // ,53
float N1S_aa_pt3p5eLarge[nPtBins_2010] = {110,67,24}; //,8.9   // wowowowow 2nd bin was 107, tricked it :) for unimportant reasons.

float N1S_aa_rap3p5[nRapBins_2013] = {492,388,403,694,677};
float N1S_aa_rap3p5e[nRapBins_2013] = {57,43,48,57,79};
//float N1S_aa_cent3p5[bin1] = {269,241,375,520,620,452,477};
// new values from july 2014, with 70-100 bin included!
/* float N1S_aa_cent3p5[nCentBins_2014] = {57,176,228,330,439,636,441,499}; // 70-100 bin, with pol2 */
/* float N1S_aa_cent3p5e[nCentBins_2014] ={10, 20, 31, 41, 45, 60, 53, 74}; */
float N1S_aa_cent3p5[nCentBins_2014] = {48.2829,173.952,172.083,288.917,438.602,646.191,411.335,412.229};
float N1S_aa_cent3p5e[nCentBins_2014] = {9.00598,17.2041,18.5096,24.7382,31.9025,40.348,33.3496,36.3436};
float N1S_aa_pt3p5U[nPtBins_2013] = {862,945,574,329,180};//,54};
float N1S_aa_pt3p5eU[nPtBins_2013] = {94.4,78.9,69.8,33.6,21.9}; //,8.9 9.3

float N2S_aa_pt3p5[bin1] = {65,53,9}; //2010 binning
float N2S_aa_pt3p5e[bin1] = {59,24,13}; //and their errors.
float N2S_aa_rap3p5[bin1] = {112,61}; //2010 binning
float N2S_aa_rap3p5e[bin1]= {46,38}; //and their errors.
//float N2S_aa_cent3p5[bin1] = {38,40,55,16,28,28,76}; //2011 binning
//float N2S_aa_cent3p5e[bin1]= {15,17,21,23,30,29,34}; //and their errors.
/* float N2S_aa_cent3p5[bin1] = {10.8,16,34,37,1.,22,12,74}; //2014 binning ! */
/* float N2S_aa_cent3p5e[bin1]= {6.9,11,16,20,22,30,26,36}; //and their errors. */
/* float N2S_aa_cent4[bin1] = {26,30,48,39,43,0.9,47}; //2011 binning ! */
/* float N2S_aa_cent4e[bin1]= {9,13,15,19,28,23,27}; //and their errors. */
float N2S_aa_cent4[bin1] = {26,21.7397,36.0665,30.8685,24.8739,2.59681,29.2158}; //
float N2S_aa_cent4e[bin1] = {9,9.76781,13.8025,15.7822,20.1582,15.5757,19.4574};// careful first bin crazy


float N1S_aa_tot4=1857;
float N1S_aa_tot4e=102;
float N1S_aa_tot3p5=2759;
float N1S_aa_tot3p5e=135;
float N2S_aa_tot4=171;
float N2S_aa_tot4e=43;
float N2S_aa_tot3p5=227;
float N2S_aa_tot3p5e=66;
float N3S_aa_tot4=23;
float N3S_aa_tot4e=43;
float N3S_aa_tot3p5=37;
float N3S_aa_tot3p5e=59;

//upper limit:
float N3S_aa_UL4=97.8;

///test bench
/* float N1S_aa_tot4=2042; */
/* float N1S_aa_tot4e=105; */
/* float N1S_aa_tot3p5=2762; */
/* float N1S_aa_tot3p5e=143; */
/* float N2S_aa_tot4=233; */
/* float N2S_aa_tot4e=53; */
/* float N2S_aa_tot3p5=225; */
/* float N2S_aa_tot3p5e=72; */
/* float N3S_aa_tot4=58; */
/* float N3S_aa_tot4e=46; */
/* float N3S_aa_tot3p5=39; */
/* float N3S_aa_tot3p5e=63; */

//LEAD LEAD
// yields from francois May 5th
/* float N1S_aa_rap3p5_2014[nRapBins_2014] = {486,520,623,532,464,113};// large bin=540} */
/* float N1S_aa_rap3p5_2014e[nRapBins_2014]= {40.6,51.7,69.3,44.3,76.2,30.3}; //largebin=78.6 */
float N1S_aa_rap3p5_2014[nRapBins_2014]={490.817,506.78,570.804,512.765,375.383,118.497};
float N1S_aa_rap3p5_2014e[nRapBins_2014]={31.3657,35.0903,38.568,37.6318,31.452,18.1268};
float N1S_aa_rap3p5_2014Large[nRapBins_2010] = {1610,970};
float N1S_aa_rap3p5_2014Largee[nRapBins_2010]= {91.3,64.7};

float N1S_aa_rap4_2014[nRapBins_2013]={312,320,397,388,426};//,334,59.6};
float N1S_aa_rap4_2014e[nRapBins_2013]={26.9,26.9,30.9,33.6,55.3};//,53.7,16.3};
float N1S_aa_rap4_2014Large[nRapBins_2010]={1080,665};
float N1S_aa_rap4_2014Largee[nRapBins_2010]={69.7,47.2};

float N2S_aa_rap3p5_2014Large[nRapBins_2010]={111,61.7};
float N2S_aa_rap3p5_2014Largee[nRapBins_2010]={45.4,38.2};
/* float N2S_aa_rap4_2014Large[nRapBins_2010]={84.3,58.8}; */
/* float N2S_aa_rap4_2014Largee[nRapBins_2010]={32.6,28.2}; */
float N2S_aa_rap4_2014Large[nRapBins_2010]={83.0694,52.3694};
float N2S_aa_rap4_2014Largee[nRapBins_2010]={29.9647,27.1147};
//from me

float N2S_aa_pt4_2013Large[nPtBins_2010]={71.,43.,20};
float N2S_aa_pt4_2013Largee[nPtBins_2010]={39.,26.,13};

//PROTON PROTON
/* float N1S_pp_rap3p5_2014[nRapBins_2014] = {1110,1130,1010,988,830,359}; //largebin=1289 */
/* float N1S_pp_rap3p5_2014e[nRapBins_2014]= {44.7,49.7,47.6,43.6,41.2,23.8}; // largebin=51.7 */
/* float N1S_pp_rap3p5_2014[nRapBins_2014]={1043.5,1124.69,1007.1,967.16,808.658,372.962}; *///Open trigger
/* float N1S_pp_rap3p5_2014e[nRapBins_2014]={36.6026,39.9872,39.123,39.7542,37.6477,23.2556}; */
float N1S_pp_rap3p5_2014[nRapBins_2014] = {1045.8,1114.99,971.794,892.208,655.377,249.099}; //highQ
float N1S_pp_rap3p5_2014e[nRapBins_2014] = {36.732,39.9483,38.4188,37.9055,32.4348,19.3658};

float N1S_pp_rap4_2014[nRapBins_2014]={785,836,699,660,601,265};  ///largebin 830};//
float N1S_pp_rap4_2014e[nRapBins_2014]={35.6,38,41.4,33.9,70.9,23}; //largebin = 36.9

/* float N2S_pp_rap3p5_2014[nRapBins_2014]= {360,347,315,341,276,63.9}; //open */
/* float N2S_pp_rap3p5_2014e[nRapBins_2014]={27.1,27.9,29.3,30,28.6,14.3}; */
float N2S_pp_rap3p5_2014[nRapBins_2014] = {316.234,338.689,308.931,305.157,229.726,50.7125}; //highQ
float N2S_pp_rap3p5_2014e[nRapBins_2014] = {24.0715,25.9207,26.8278,27.3601,23.4638,12.9031};
/* float N2S_pp_rap4_2014[nRapBins_2014]={265,277,239,241,254,62.3}; //open */
/* float N2S_pp_rap4_2014e[nRapBins_2014]={22.9,24.3,24,24.9,47.4,13.9}; */
float N2S_pp_rap4_2014[nRapBins_2014] = {241.004,257.264,229.192,218.071,179.066,45.9916}; //highQ
float N2S_pp_rap4_2014e[nRapBins_2014] = {19.9981,21.644,21.2744,22.3079,20.8062,11.7197};
/* float N2S_pp_rap4_2014Large[nRapBins_2010]={2290,1470}; */
/* float N2S_pp_rap4_2014Largee[nRapBins_2010]={67.9,63.6}; */
/* float N2S_pp_rap4_2014Large[nRapBins_2010]={714.061,497.992}; */ //open
/* float N2S_pp_rap4_2014Largee[nRapBins_2010]={36,37.8979}; */
float N2S_pp_rap4_2014Large[nRapBins_2010] = {706.144,435.857}; //highQ
float N2S_pp_rap4_2014Largee[nRapBins_2010] = {35.7708,36.5077};

/* float N3S_pp_rap3p5_2014[nRapBins_2013]={170,175,179,121,189};//,111,74.9}; //open */
/* float N3S_pp_rap3p5_2014e[nRapBins_2013]={21.2,23.3,25.7,24.9,28.2};//;,24.1,14.3}; */
float N3S_pp_rap3p5_2014[nRapBins_2014] = {140.577,167.214,165.89,113.125,87.8148,56.0211}; //highQ
float N3S_pp_rap3p5_2014e[nRapBins_2014] = {19.0652,21.8788,23.189,22.6344,19.3827,12.5086};
/* float N3S_pp_rap4_2014[nRapBins_2014]={134,138,138,96,108,67.3}; //175} */
/* float N3S_pp_rap4_2014e[nRapBins_2014]={18.1,20.7,20.5,20.7,37,14.5}; // 23.4 */
/* float N3S_pp_rap4_2014[nRapBins_2014]={117.8,126.87,131.51,95.9781,92.4825,67.7779}; //open */
/* float N3S_pp_rap4_2014e[nRapBins_2014]={15.9977,17.7523,18.5524,19.0335,31.3553,13.4012}; */
float N3S_pp_rap4_2014[nRapBins_2014] = {115.373,128.482,135.999,85.7138,71.7167,46.4221}; //highQ
float N3S_pp_rap4_2014e[nRapBins_2014] = {16.0236,17.9384,18.4957,18.3729,17.3947,11.4625};
float N3S_pp_rap4_2014Large[nRapBins_2010] = {360.713,196.843};
float N3S_pp_rap4_2014Largee[nRapBins_2010] = {29.9762,31.3335};

float N2S_pp_pt3p5_2013[nPtBins_2013] = {433,529,378,266,144};//31.2
float N2S_pp_pt3p5_2013e[nPtBins_2013] = {45,41,28,23,15};//7.6
/* float N2S_pp_pt4_2013[nPtBins_2013]={372.739,250.364,251.718,200.114,124.341}; */ //open
/* float N2S_pp_pt4_2013e[nPtBins_2013]={28.082,23.514,21.9586,17.523,13.319}; */
float N2S_pp_pt4_2013[nPtBins_2013] = {301.791,257.089,243.196,188.89,113.466};//highQ
float N2S_pp_pt4_2013e[nPtBins_2013] = {28.9903,23.9729,20.3855,17.0605,12.7171};
/* float N2S_pp_pt4_2013[nPtBins_2013] = {318,302,274,230,131};//31.7 */
/* float N2S_pp_pt4_2013e[nPtBins_2013] = {34,30,22,21,15};//8.0 */
/* float N2S_pp_pt4_2013Large[nPtBins_2010]={827,287,200}; */
/* float N2S_pp_pt4_2013Largee[nPtBins_2010]={168,23,31}; */
/* float N2S_pp_pt4_2013Large[nPtBins_2010]={725.214,248.916,189.512};//open */
/* float N2S_pp_pt4_2013Largee[nPtBins_2010]={40.9265,20.2315,16.7928}; */
float N2S_pp_pt4_2013Large[nPtBins_2010] = {695.703,234.325,173.997};//highQ
float N2S_pp_pt4_2013Largee[nPtBins_2010] = {39.8716,19.732,16.0895};

/* float N3S_pp_pt4_2013[nPtBins_2013] = {141,199,134,108,75};//23.4 */
/* float N3S_pp_pt4_2013e[nPtBins_2013] = {27,28,18,16,12};//7.2 */
/* float N3S_pp_pt4_2013[nPtBins_2013]={186.882,165.105,118.457,89.0954,71.8379}; //open */
/* float N3S_pp_pt4_2013e[nPtBins_2013]={23.4312,21.1747,17.7216,13.491,10.8839}; */
float N3S_pp_pt4_2013[nPtBins_2013] = {128.482,168.485,112.074,84.2136,69.3291}; //highQ
float N3S_pp_pt4_2013e[nPtBins_2013] = {24.0888,21.6864,16.513,13.1443,10.6023};

float N3S_pp_pt3p5_2013[nPtBins_2013] = {203,282,184,130,83};//21
float N3S_pp_pt3p5_2013e[nPtBins_2013] = {36,34,23,17,13};//6.5
float N3S_pp_pt4_2013Large[nPtBins_2010] = {358.951,91.1538,99.0478};
float N3S_pp_pt4_2013Largee[nPtBins_2010] = {34.5618,15.0331,13.0713};


////yields in the note
//trying now with yields from plots, waiting for the Arleo ones. 25th apr. 2014
/* float N1S_pp_pt3p5[nPtBins_2013] = {1717,1550,1107,691,352}; // ,62.5 */
/* float N1S_pp_pt3p5e[nPtBins_2013] = {80,64,43,43,23}; // ,9.9 */
/* float N1S_pp_pt3p5[nPtBins_2013]={1623.31,1520.25,1023.34,633.098,344.566}; */ //open
/* float N1S_pp_pt3p5e[nPtBins_2013]={50.737,47.6386,36.7225,28.3172,20.396}; */
float N1S_pp_pt3p5[nPtBins_2013] = {1499.27,1418.6,972.688,602.121,337.012}; //highQ
float N1S_pp_pt3p5e[nPtBins_2013] = {49.2967,45.7898,28.4194,27.578,20.1226};
float N1S_pp_pt4[nPtBins_2013] = {954.099,806.429,719.354,515.821,304.747};
float N1S_pp_pt4e[nPtBins_2013] = {39.8429,34.6068,30.325,25.261,19.0581};



float N1S_pp_pt3p5U[nPtBins_2013] = {1718,1546.18,1109.26,692.5,350}; // ,62.8
float N1S_pp_pt3p5eU[nPtBins_2013]= {81,65,44,44,23}; // ,9.9,10
float N1S_pp_pt3p5Large[nPtBins_2010]={3856,910,586};
float N1S_pp_pt3p5eLarge[nPtBins_2010]={97,37,29};
float N1S_pp_rap3p5[nRapBins_2013] = {1102,844,795,1198,1427}; 
float N1S_pp_rap3p5e[nRapBins_2013] = {47,42,37,47,49};


float N2S_pp_pt3p5Large[nPtBins_2010] = {1166,335,227}; // for CS2S_pp_ptLarge
float N2S_pp_pt3p5Largee[nPtBins_2010]= {64,26,35};
float N2S_pp_pt3p5[nPtBins_2013] = {372.615,491.722,318.957,218.694,123.526};
float N2S_pp_pt3p5e[nPtBins_2013] = {32.2694,32.2999,19.49,18.7975,13.4272};


float N2S_pp_rap3p5[nRapBins_2010] = {1033,697}; // for CS2S_pp_rapLarge
float N2S_pp_rap3p5e[nRapBins_2010]= {51,41};


/* float N1S_pp_tot3p5=5512; */
/* float N1S_pp_tot3p5e=116; */

/* float N2S_pp_tot3p5=1770; */
/* float N2S_pp_tot3p5e=72; */
/* float N2S_pp_tot4=1404; */
/* float N2S_pp_tot4e=65; */

/* float N3S_pp_tot3p5=906; */
/* float N3S_pp_tot3p5e=60; */
/* float N3S_pp_tot4= 755; */
/* float N3S_pp_tot4e=53; */


// test bench
float N1S_pp_tot3p5=4977;
float N1S_pp_tot3p5e=87;
float N1S_pp_tot4=3485;
float N1S_pp_tot4e=72;
float N2S_pp_tot3p5=1569;
float N2S_pp_tot3p5e=59;
float N2S_pp_tot4=1183;
float N2S_pp_tot4e=49;
float N3S_pp_tot3p5=757;
float N3S_pp_tot3p5e=50;
float N3S_pp_tot4= 601;
float N3S_pp_tot4e=41;


///p-Pb (1st part of the run)
float N1S_pa_pt3p5[nPtBins_2014]={1182,1010,940,659,381,77};
float N1S_pa_pt3p5e[nPtBins_2014]={136,37,49,40,23,10};
float N2S_pa_pt3p5[nPtBins_2014]={181,218,220,191,77,37.1};
float N2S_pa_pt3p5e[nPtBins_2014]={62,30,29,26,13,7.9};
float N3S_pa_pt3p5[nPtBins_2014]={37,71,48,46,74,27.2};
float N3S_pa_pt3p5e[nPtBins_2014]={40,27,23,19,12,7.2};

float N1S_pp_m146p239pt3p5[nPtBins_2014]={1766,1518,1091,664,347,62.2};
float N1S_pp_m146p239pt3p5e[nPtBins_2014]={57,48,42,26,23,9.6};
float N2S_pp_m146p239pt3p5[nPtBins_2014]={432,512,370,251,139,31.4};
float N2S_pp_m146p239pt3p5e[nPtBins_2014]={35,37,27,16,15,7.6};
float N3S_pp_m146p239pt3p5[nPtBins_2014]={192,266,176,118,79,21.4};
float N3S_pp_m146p239pt3p5e[nPtBins_2014]={28,37,22,13,13,6.4};


/* float N1S_pa_pt4[nPtBins_2014]={647,}; */
/* float N1S_pa_pt4e[nPtBins_2014]={50}; */
/* float N2S_pa_pt4[nPtBins_2014]={647,}; */
/* float N2S_pa_pt4e[nPtBins_2014]={50}; */
/* float N3S_pa_pt4[nPtBins_2014]={647,}; */
/* float N3S_pa_pt4e[nPtBins_2014]={50}; */

/* float N1S_pa_rap3p5[nRapBins_pPb]={647,}; */
/* float N1S_pa_rap3p5e[nRapBins_pPb]={50}; */
/* float N2S_pa_rap3p5[nRapBins_pPb]={647,}; */
/* float N2S_pa_rap3p5e[nRapBins_pPb]={50}; */
/* float N3S_pa_rap3p5[nRapBins_pPb]={647,}; */
/* float N3S_pa_rap3p5e[nRapBins_pPb]={50}; */

/* float N1S_pa_rap4[nRapBins_pPb]={647,}; */
/* float N1S_pa_rap4e[nRapBins_pPb]={50}; */
/* float N2S_pa_rap4[nRapBins_pPb]={647,}; */
/* float N2S_pa_rap4e[nRapBins_pPb]={50}; */
/* float N3S_pa_rap4[nRapBins_pPb]={647,}; */
/* float N3S_pa_rap4e[nRapBins_pPb]={50}; */



//// systematics

float N1S_pp_pt3p5s[nPtBins_2013]={0.134,0.0827,0.1236,0.0811,0.0345};
float N1S_aa_pt3p5s[nPtBins_2013]={0.316,0.1742,0.3831,0.0555,0.0410};

float N1S_pp_rap3p5s[nRapBins_2014]={0.0873,0.0630,0.0522,0.0962,0.0629,0.0436};
float N1S_aa_rap3p5s[nRapBins_2014]={0.1122,0.0746,0.1883,0.2029,0.4291,0.3936};

// from fitting.
float N1S_pp_tot3p5s=0.0673;
float N1S_aa_tot3p5s=0.1361;
 // from fitting.
float N1S_pp_tot4s=0.0791;
float N1S_aa_tot4s=0.1211;

//fit+bkg vars
float N1S_pp_pt3p5s_2p5[nfitvars] = {1550.15,1497.98,1591.3,1614.17,1574.54};
float N1S_pp_pt3p5s_5[nfitvars] = {1410.97,1412.42,1455.5,1460.92,1454.6};
float N1S_pp_pt3p5s_8[nfitvars] = {976.736,986.815,988.917,990.409,989.89};
float N1S_pp_pt3p5s_12[nfitvars] = {625.377,623.401,625.244,629.017,621.985};
float N1S_pp_pt3p5s_20[nfitvars] = {336.517,338.115,340.596,338.126,341.151};
float N2S_pp_pt4s_2p5[nfitvars] = {287.297,298.987,305.301,302.169,304.997};
float N2S_pp_pt4s_5[nfitvars] = {245.22,255.444,258.758,257.947,260.807};
float N2S_pp_pt4s_8[nfitvars] = {247.881,249.625,250.597,249.841,250.459};
float N2S_pp_pt4s_12[nfitvars] = {198.904,199.68,199.002,198.033,197.859};
float N2S_pp_pt4s_20[nfitvars] = {113.168,113.473,113.767,113.126,113.73};
float N3S_pp_pt4s_2p5[nfitvars] = {118.876,126.71,130.525,128.681,130.324};
float N3S_pp_pt4s_5[nfitvars] = {158.647,167.073,169.625,168.936,170.851};
float N3S_pp_pt4s_8[nfitvars] = {116.325,117.365,116.398,114.615,115.738};
float N3S_pp_pt4s_12[nfitvars] = {93.8629,92.1802,89.1508,86.8052,88.7491};
float N3S_pp_pt4s_20[nfitvars] = {68.4465,69.3293,69.5748,69.1804,69.4976};
float N1S_pp_rap3p5s_0p4[nfitvars] = {1083.6,1073.14,1091.47,1100.48,1079.09};
float N1S_pp_rap3p5s_0p8[nfitvars] = {1139.36,1097.09,1147.98,1118.5,1178.18};
float N1S_pp_rap3p5s_1p2[nfitvars] = {956.711,969.94,984.275,984.072,983.597};
float N1S_pp_rap3p5s_1p6[nfitvars] = {924.696,894.274,931.587,948.78,924.45};
float N1S_pp_rap3p5s_2p0[nfitvars] = {641.296,646.593,710.447,656.026,704.772};
float N1S_pp_rap3p5s_2p4[nfitvars] = {241,250.623,242.982,277.994,243.311};
float N2S_pp_rap4s_0p4[nfitvars] = {247.591,241.184,257.802,256.581,254.617};
float N2S_pp_rap4s_0p8[nfitvars] = {265.286,251.534,268.305,258.135,275.316};
float N2S_pp_rap4s_1p2[nfitvars] = {225.683,232.86,236.048,234.933,235.659};
float N2S_pp_rap4s_1p6[nfitvars] = {229.714,218.801,231.935,233.463,229.413};
float N2S_pp_rap4s_2p0[nfitvars] = {172.425,174.556,194.284,179.005,191.163};
float N2S_pp_rap4s_2p4[nfitvars] = {44.9216,47.3969,45.7954,45.8028,46.0241};
float N3S_pp_rap4s_0p4[nfitvars] = {120.145,115.556,126.331,125.206,123.859};
float N3S_pp_rap4s_0p8[nfitvars] = {134.834,124.428,137.599,128.89,139.942};
float N3S_pp_rap4s_1p2[nfitvars] = {133.156,138.711,139.785,138.495,139.226};
float N3S_pp_rap4s_1p6[nfitvars] = {93.4876,86.1959,91.1151,91.5756,90.268};
float N3S_pp_rap4s_2p0[nfitvars] = {67.7696,69.008,71.1997,71.6343,73.8208};
float N3S_pp_rap4s_2p4[nfitvars] = {44.5887,48.1206,46.2762,46.1297,46.4448};
float N1S_aa_pt3p5s_2p5[nfitvars] = {864.322,919.319,779.047,733.57,741.833};
float N1S_aa_pt3p5s_5[nfitvars] = {775.593,741.434,780.944,779.249,772.128};
float N1S_aa_pt3p5s_8[nfitvars] = {541.362,720.636,497.992,505.485,494.108};
float N1S_aa_pt3p5s_12[nfitvars] = {338.81,349.307,364.487,364.736,364.633};
float N1S_aa_pt3p5s_20[nfitvars] = {338.81,349.307,364.487,364.736,364.633};
float N2S_aa_pt4s_6p5[nfitvars] = {65.4384,63.553,59.7498,64.3404,63.8061};
float N2S_aa_pt4s_10[nfitvars] = {26.0016,35.9108,20.9523,21.2272,21.1887};
float N2S_aa_pt4s_20[nfitvars] = {21.4984,21.1819,22.1259,22.1392,22.1022};
float N1S_aa_rap3p5s_0p4[nfitvars] = {502.025,502.361,493.918,497.558,496.988};
float N1S_aa_rap3p5s_0p8[nfitvars] = {540.836,530.507,495.908,498.547,495.499};
float N1S_aa_rap3p5s_1p2[nfitvars] = {601.225,619.093,633.024,636.619,629.831};
float N1S_aa_rap3p5s_1p6[nfitvars] = {575.222,511.592,564.022,564.621,564.412};
float N1S_aa_rap3p5s_2p0[nfitvars] = {491.91,399.837,415.726,414.338,381.992};
float N1S_aa_rap3p5s_2p4[nfitvars] = {114.502,143.863,132.046,132.221,132.108};
float N2S_aa_rap4s_1p2[nfitvars] = {82.2576,85.7512,88.0356,86.1159,78.1231};
float N2S_aa_rap4s_2p4[nfitvars] = {89.3027,74.3289,59.6,62.39,63.4262};
float N1S_aa_cents_5[nfitvars] = {534.884,473.674,417.655,422.485,423.043};
float N1S_aa_cents_10[nfitvars] = {490.145,470.022,434.087,420.725,433.939};
float N1S_aa_cents_20[nfitvars] = {668.75,736.922,649.96,696.419,653.381};
float N1S_aa_cents_30[nfitvars] = {442.688,435.957,441.18,439.612,442.168};
float N1S_aa_cents_40[nfitvars] = {331.487,328.472,318.416,305.767,318.83};
float N1S_aa_cents_50[nfitvars] = {216.306,197.265,186.187,182.28,184.922};
float N1S_aa_cents_70[nfitvars] = {170.127,167.533,183.742,174.841,183.234};
float N1S_aa_cents_100[nfitvars] = {};
float N2S_aa_cents_5[nfitvars] = {48.4975,43.0919,33.6128,33.8641,32.8193};
float N2S_aa_cents_10[nfitvars] = {5.94285,8.63158,-0.318017,0.914677,0.190994};
float N2S_aa_cents_20[nfitvars] = {43.2279,43.5586,35.9874,35.2762,35.3462};
float N2S_aa_cents_30[nfitvars] = {33.9329,33.216,37.4242,33.6887,35.0151};
float N2S_aa_cents_40[nfitvars] = {43.3892,40.9809,43.8794,44.6109,44.818};
float N2S_aa_cents_50[nfitvars] = {23.5914,27.0468,23.9,24.106,23.8406};
float N2S_aa_cents_100[nfitvars] = {};
float N1B_pp_pt3p5s_2p5[nbkgdvars] = {1496.56,1496.98};
float N1B_pp_pt3p5s_5[nbkgdvars] = {1411.8,1414.98};
float N1B_pp_pt3p5s_8[nbkgdvars] = {960.313,962.119};
float N1B_pp_pt3p5s_12[nbkgdvars] = {602.519,602.295};
float N1B_pp_pt3p5s_20[nbkgdvars] = {335.663,336.39};
float N2B_pp_pt4s_2p5[nbkgdvars] = {359.307,353.335};
float N2B_pp_pt4s_5[nbkgdvars] = {266.619,262.225};
float N2B_pp_pt4s_8[nbkgdvars] = {237.328,237.889};
float N2B_pp_pt4s_12[nbkgdvars] = {188.656,188.804};
float N2B_pp_pt4s_20[nbkgdvars] = {112.203,112.903};
float N3B_pp_pt4s_2p5[nbkgdvars] = {160.93,169.281};
float N3B_pp_pt4s_5[nbkgdvars] = {177.804,179.161};
float N3B_pp_pt4s_8[nbkgdvars] = {109.725,110.124};
float N3B_pp_pt4s_12[nbkgdvars] = {83.8668,84.1256};
float N3B_pp_pt4s_20[nbkgdvars] = {67.9136,68.656};
float N1B_pp_rap3p5s_0p4[nbkgdvars] = {1031.25,1037.63};
float N1B_pp_rap3p5s_0p8[nbkgdvars] = {1117.05,1114.86};
float N1B_pp_rap3p5s_1p2[nbkgdvars] = {953.639,956.322};
float N1B_pp_rap3p5s_1p6[nbkgdvars] = {883.605,881.93};
float N1B_pp_rap3p5s_2p0[nbkgdvars] = {655.382,655.039};
float N1B_pp_rap3p5s_2p4[nbkgdvars] = {237.666,249.423};
float N2B_pp_rap4s_0p4[nbkgdvars] = {239.672,240.333};
float N2B_pp_rap4s_0p8[nbkgdvars] = {257.227,257.248};
float N2B_pp_rap4s_1p2[nbkgdvars] = {215.534,223.55};
float N2B_pp_rap4s_1p6[nbkgdvars] = {214.694,209.404};
float N2B_pp_rap4s_2p0[nbkgdvars] = {173.098,176.027};
float N2B_pp_rap4s_2p4[nbkgdvars] = {45.9987,43.8763};
float N3B_pp_rap4s_0p4[nbkgdvars] = {115.068,115.21};
float N3B_pp_rap4s_0p8[nbkgdvars] = {128.432,128.471};
float N3B_pp_rap4s_1p2[nbkgdvars] = {129.472,136.97};
float N3B_pp_rap4s_1p6[nbkgdvars] = {85.9007,83.6593};
float N3B_pp_rap4s_2p0[nbkgdvars] = {69.6524,71.8443};
float N3B_pp_rap4s_2p4[nbkgdvars] = {46.4128,45.8699};
float N1B_aa_pt3p5s_2p5[nbkgdvars] = {717.243,712.649};
float N1B_aa_pt3p5s_5[nbkgdvars] = {699.164,703.814};
float N1B_aa_pt3p5s_8[nbkgdvars] = {511.1,517.071};
float N1B_aa_pt3p5s_12[nbkgdvars] = {350.612,350.879};
float N1B_aa_pt3p5s_20[nbkgdvars] = {350.612,350.879};
float N2B_aa_pt4s_6p5[nbkgdvars] = {};
float N2B_aa_pt4s_10[nbkgdvars] = {};
float N2B_aa_pt4s_20[nbkgdvars] = {};
float N1B_aa_rap3p5s_0p4[nbkgdvars] = {482.485,482.742};
float N1B_aa_rap3p5s_0p8[nbkgdvars] = {497.149,502.336};
float N1B_aa_rap3p5s_1p2[nbkgdvars] = {570.466,570.447};
float N1B_aa_rap3p5s_1p6[nbkgdvars] = {545.56,549.237};
float N1B_aa_rap3p5s_2p0[nbkgdvars] = {485.341,360.041};
float N1B_aa_rap3p5s_2p4[nbkgdvars] = {122.955,117.653};
float N2B_aa_rap4s_1p2[nbkgdvars] = {};
float N2B_aa_rap4s_2p4[nbkgdvars] = {};
float N1B_aa_cents_5[nbkgdvars] = {};
float N1B_aa_cents_10[nbkgdvars] = {};
float N1B_aa_cents_20[nbkgdvars] = {};
float N1B_aa_cents_30[nbkgdvars] = {};
float N1B_aa_cents_40[nbkgdvars] = {};
float N1B_aa_cents_50[nbkgdvars] = {};
float N1B_aa_cents_70[nbkgdvars] = {};
float N1B_aa_cents_100[nbkgdvars] = {};
float N2B_aa_cents_5[nbkgdvars] = {};
float N2B_aa_cents_10[nbkgdvars] = {};
float N2B_aa_cents_20[nbkgdvars] = {};
float N2B_aa_cents_30[nbkgdvars] = {};
float N2B_aa_cents_40[nbkgdvars] = {};
float N2B_aa_cents_50[nbkgdvars] = {};
float N2B_aa_cents_100[nbkgdvars] = {};
