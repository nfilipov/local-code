
  //integers for binning.
  //centrality binning in 11-011


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
#define N_MB_uncorr 1126653312	
#define N_MB_corr   1138033648
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
  
  float raaMB1S[1]={0.388};
  float raaMB1SstatErr[1]={0.018};
  float raaMB1SsystErr[1]={0.05};
  float raaMB2S[1]={0.106};
  float raaMB2SstatErr[1]={0.026};
  float raaMB2SsystErr[1]={0.021};
  float raaMB3S[1]={0.033};
  float raaMB3SstatErr[1]={0.035};//{0.035};
  float raaMB3SsystErr[1]={0.006};
  
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
 
 ////EFFICIENCY and acceptance 24th apr. 2014
////so wrong..
 //A. Pythia sample. 1S //table 14 and 15
/* float Ae_1S_pythia_pt[nPtBins_2013] = {0.325,0.217,0.221,0.316,0.445}; */
/* float Ae_1S_pythia_pte[nPtBins_2013] = {0.001,0.0006,0.0008,0.0015,0.0029}; */
 /* float Ae_1S_pythia_rap[nRapBins_2013] = {0.246,0.265,0.299,0.315,0.221}; */
 /* float Ae_1S_pythia_rape[nRapBins_2013] = {0.0008,0.001,0.0011,0.0010,0.0008};  */
/* float Ae_1S_pythia_rap2014[nRapBins_2013]={0.2551,0.2799,0.3196,0.3174,0.2053};//0.2053};//bug is here, new table 15 of note may 6th. */
/* float Ae_1S_pythia_rap2014e[nRapBins_2013]={0.00495,0.00548,0.006429,0.006941,0.005202};//,0.005202}; */


////JUNE 18th
// total acc*eff 1S = 0.449 +/- 0.0007
/* float Ae_1S_pythia_pt[nPtBins_2013] = {0.531,0.388,0.390,0.485,0.602}; */
/* float Ae_1S_pythia_pte[nPtBins_2013] = {0.00151,0.0010,0.0014,0.0023,0.0040}; */
///// July 15th - Prashant
float Ae_1S_pythia_pt[nPtBins_2013] = {0.296,0.190,0.185,0.259,0.376};
float Ae_1S_pythia_pte[nPtBins_2013] = {0.001,0.001,0.001,0.002,0.004};
/* float Ae_1S_pythia_rap2014[nRapBins_2014]={0.526,0.473,0.452,0.433,0.411,0.294}; */
/* float Ae_1S_pythia_rap2014e[nRapBins_2014]={0.0018,0.0016,0.0015,0.0015,0.0017,0.0022}; */
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

//I assume these are with 3p5 ? //table 21
float Ae_1S_pyquen_DD020[nRapBins_2010]={0.266,0.305};
float Ae_1S_pyquen_DD20100[nRapBins_2010]={0.278,0.328};
float Ae_1S_pyquen_DD020e[nRapBins_2010]={0.002,0.0032};
float Ae_1S_pyquen_DD20100e[nRapBins_2010]={0.0017,0.0027};


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
float t_1S_pythia_pt3p5[nPtBins_2013]={1.102,1.095,1.075,1.057,1.039};
float t_1S_pythia_pt3p5e[nPtBins_2013]={0.004,0.004,0.005,0.006,0.008};
float t_1S_pythia_rap3p5[nRapBins_2014]={1.069,1.069,1.073,1.090,1.128,1.146};
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
float t_2S_pythia_pt3p5[nPtBins_2013]={1.078,1.071,1.057,1.046,1.031};
float t_2S_pythia_pt3p5e[nPtBins_2013]={0.004,0.004,0.005,0.007,0.010};
float t_2S_pythia_rap3p5[nRapBins_2014]={1.049,1.048,1.052,1.074,1.119,1.145};
float t_2S_pythia_rap3p5e[nRapBins_2014]={0.005,0.005,0.005,0.005,0.006,0.012};
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
float t_2S_pythia_rap2010[nRapBins_2010]={0.248,0.175};
float t_2S_pythia_rap2010e[nRapBins_2010]={0.001,0.001};
float A_2S_pythia_pt2010[nPtBins_2010]= {0.28,0.25,0.41};
float A_2S_pythia_pt2010e[nPtBins_2010]={0.0002,0.0002,0.0041};
float e_2S_pythia_pt2010[nPtBins_2010]= {0.714,0.735,0.757};
float e_2S_pythia_pt2010e[nPtBins_2010]={0.0015,0.004,0.0051};
float t_2S_pythia_pt2010[nPtBins_2010] = {0.211,0.195,0.320};
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
float t_1S_pyquen_pt[nPtBins_2013]={1.111,1.102,1.078,1.056,1.034};
float t_1S_pyquen_pte[nPtBins_2013]={0.012,0.013,0.014,0.010,0.008};
float t_1S_pyquen_rap[nRapBins_2014]={1.064,1.064,1.073,1.101,1.150,1.177};
float t_1S_pyquen_rape[nRapBins_2014]={0.013,0.013,0.013,0.014,0.018,0.033};
float t_1S_pyquen_cent2014[nCentBins_2014]={1.093,1.090,1.090,1.090,1.090,1.090,1.090,1.090};
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
/* float t_1S_pyquen_DD020[nRapBins_2010]= {}; */
/* float t_1S_pyquen_rap2010e[nRapBins_2010]={025}; */
/* float t_1S_pyquen_pt2010[nPtBins_2010]=   {0.41}; */
/* float t_1S_pyquen_pt2010e[nPtBins_2010]=  {009,0.0024}; */

float t_2S_pyquen_pt[nPtBins_2013]={1.083,1.074,1.058,1.044,1.028};
float t_2S_pyquen_pte[nPtBins_2013]={0.012,0.014,0.015,0.011,0.007};
float t_2S_pyquen_rap[nRapBins_2014]={1.044,1.044,1.050,1.081,1.133,1.164};
float t_2S_pyquen_rape[nRapBins_2014]={0.013,0.013,0.013,0.014,0.018,0.033};
float t_2S_pyquen_cent2014[nCentBins_2014]={1.070,1.069,1.070,1.070,1.070,1.069,1.068,1.068};
float t_2S_pyquen_cent2014e[nCentBins_2014]={0.009,0.010,0.013,0.013,0.013,0.013,0.018,0.016};
float e_2S_pyquen_rap2010[nRapBins_2010]= {0.710,0.647};
float e_2S_pyquen_rap2010e[nRapBins_2010]={0.005,0.006};
float A_2S_pyquen_rap2010[nRapBins_2010]= {0.343,0.424};
float A_2S_pyquen_rap2010e[nRapBins_2010]={0.0016,0.0025};
float e_2S_pyquen_pt2010[nPtBins_2010]=   {0.674,0.718,0.749};
float e_2S_pyquen_pt2010e[nPtBins_2010]=  {0.0047,0.0074,0.0049};
float A_2S_pyquen_pt2010[nPtBins_2010]=   {0.28,0.25,0.41};
float A_2S_pyquen_pt2010e[nPtBins_2010]=  {0.0004,0.0009,0.0024};
float t_2S_pyquen_rap2010[nRapBins_2010]= {1.046,1.109};
float t_2S_pyquen_rap2010e[nRapBins_2010]={0.008,0.012};
float t_2S_pyquen_pt2010[nPtBins_2010]=   {1.076,1.051,1.033};
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

///Tg and probe correction of corrections
float Aet_1S_pythia_pt[nPtBins_2013]={0.326,0.208,0.199,0.274,0.391};
float Aet_1S_pythia_pte[nPtBins_2013]={0.002,0.001,0.002,0.003,0.005};
float Aet_1S_pythia_rap2014[nRapBins_2014]={0.287,0.287,0.291,0.277,0.230,0.081};
float Aet_1S_pythia_rap2014e[nRapBins_2014]={0.002,0.002,0.002,0.002,0.002,0.001};

float Aet_2S_pythia_pt2013[nPtBins_2013]={0.294,0.168,0.167,0.236,0.365};
float Aet_2S_pythia_pt2013e[nPtBins_2013]={0.002,0.001,0.002,0.003,0.007};
float Aet_2S_pythia_rap2014[nRapBins_2014]={0.248,0.249,0.246,0.235,0.197,0.071};
float Aet_2S_pythia_rap2014e[nRapBins_2014]={0.002,0.002,0.002,0.002,0.002,0.001};

float Aet_3S_pythia_pt2013[nPtBins_2013]={0.301,0.171,0.172,0.245,0.379};
float Aet_3S_pythia_pt2013e[nPtBins_2013]={0.002,0.001,0.002,0.003,0.006};
float Aet_3S_pythia_rap2014[nRapBins_2014]={0.257,0.256,0.251,0.241,0.203,0.073};
float Aet_3S_pythia_rap2014e[nRapBins_2014]={0.002,0.002,0.002,0.002,0.002,0.001};

float Aet_1S_pyquen_pt[nPtBins_2013]={0.290,0.186,0.186,0.266,0.383};
float Aet_1S_pyquen_pte[nPtBins_2013]={0.004,0.003,0.003,0.004,0.005};
float Aet_1S_pyquen_rap2014[nRapBins_2014]={0.238,0.250,0.269,0.267,0.227,0.085};
float Aet_1S_pyquen_rap2014e[nRapBins_2014]={0.004,0.004,0.005,0.005,0.005,0.003};
float Aet_1S_pyquen_cent2014[nCentBins_2014]={0.244,0.243,0.240,0.238,0.237,0.233,0.228,0.220};
float Aet_1S_pyquen_cent2014e[nCentBins_2014]={0.003,0.003,0.004,0.003,0.004,0.003,0.005,0.004};

float Aet_2S_pyquen_pt3p5[nPtBins_2013]={0.276,0.160,0.161,0.231,0.359};
float Aet_2S_pyquen_pt3p5e[nPtBins_2013]={0.004,0.003,0.003,0.004,0.006};
float Aet_2S_pyquen_rap3p5[nRapBins_2014]={0.220,0.230,0.237,0.234,0.201,0.072};
float Aet_2S_pyquen_rap3p5e[nRapBins_2014]={0.004,0.004,0.005,0.005,0.005,0.003};
float Aet_2S_pyquen_cent2014[nCentBins_2014]={0.215,0.214,0.213,0.212,0.210,0.207,0.205,0.199};
float Aet_2S_pyquen_cent2014[nCentBins_2014]={0.003,0.003,0.004,0.004,0.004,0.003,0.005,0.004};








////yields in the note
/* float N1S_aa_pt3p5[nPtBins_2013] = {863,929,572,346,184}; // ,53 *///prior to Oct 20!
/* float N1S_aa_pt3p5e[nPtBins_2013] = {92,75,66,32,21}; //,8.9 */
float N1S_aa_pt3p5[nPtBins_2013] = {691,686,507,343,175}; // ,53
float N1S_aa_pt3p5e[nPtBins_2013] = {42,41,66,32,21}; //,8.9
float N1S_aa_pt3p5Large[nPtBins_2010] = {1842,466,311}; // ,53
float N1S_aa_pt3p5eLarge[nPtBins_2010] = {110,67,24}; //,8.9   // wowowowow 2nd bin was 107, tricked it :) for unimportant reasons.

float N1S_aa_rap3p5[nRapBins_2013] = {492,388,403,694,677};
float N1S_aa_rap3p5e[nRapBins_2013] = {57,43,48,57,79};
//float N1S_aa_cent3p5[bin1] = {269,241,375,520,620,452,477};
// new values from july 2014, with 70-100 bin included!
float N1S_aa_cent3p5[nCentBins_2014] = {57,176,228,330,439,636,441,499}; // 70-100 bin, with pol2
float N1S_aa_cent3p5e[nCentBins_2014] ={10, 20, 31, 41, 45, 60, 53, 74};
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
float N2S_aa_cent4[bin1] = {26,30,48,39,43,0.9,47}; //2011 binning !
float N2S_aa_cent4e[bin1]= {9,13,15,19,28,23,27}; //and their errors.

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
float N1S_aa_rap3p5_2014[nRapBins_2014] = {486,520,623,532,464,113};// large bin=540}
float N1S_aa_rap3p5_2014e[nRapBins_2014]= {40.6,51.7,69.3,44.3,76.2,30.3}; //largebin=78.6
float N1S_aa_rap3p5_2014Large[nRapBins_2010] = {1610,970};
float N1S_aa_rap3p5_2014Largee[nRapBins_2010]= {91.3,64.7};

float N1S_aa_rap4_2014[nRapBins_2013]={312,320,397,388,426};//,334,59.6};
float N1S_aa_rap4_2014e[nRapBins_2013]={26.9,26.9,30.9,33.6,55.3};//,53.7,16.3};
float N1S_aa_rap4_2014Large[nRapBins_2010]={1080,665};
float N1S_aa_rap4_2014Largee[nRapBins_2010]={69.7,47.2};

float N2S_aa_rap3p5_2014Large[nRapBins_2010]={111,61.7};
float N2S_aa_rap3p5_2014Largee[nRapBins_2010]={45.4,38.2};
float N2S_aa_rap4_2014Large[nRapBins_2010]={84.3,58.8};
float N2S_aa_rap4_2014Largee[nRapBins_2010]={32.6,28.2};
//from me

float N2S_aa_pt4_2013Large[nPtBins_2010]={71.,43.,20};
float N2S_aa_pt4_2013Largee[nPtBins_2010]={39.,26.,13};

//PROTON PROTON
float N1S_pp_rap3p5_2014[nRapBins_2014] = {1110,1130,1010,988,830,359}; //largebin=1289
float N1S_pp_rap3p5_2014e[nRapBins_2014]= {44.7,49.7,47.6,43.6,41.2,23.8}; // largebin=51.7
float N1S_pp_rap4_2014[nRapBins_2014]={785,836,699,660,601,265};  ///largebin 830};//
float N1S_pp_rap4_2014e[nRapBins_2014]={35.6,38,41.4,33.9,70.9,23}; //largebin = 36.9

float N2S_pp_rap3p5_2014[nRapBins_2014]= {360,347,315,341,276,63.9};
float N2S_pp_rap3p5_2014e[nRapBins_2014]={27.1,27.9,29.3,30,28.6,14.3};
float N2S_pp_rap4_2014[nRapBins_2014]={265,277,239,241,254,62.3};
float N2S_pp_rap4_2014e[nRapBins_2014]={22.9,24.3,24,24.9,47.4,13.9};
float N2S_pp_rap4_2014Large[nRapBins_2010]={2290,1470};
float N2S_pp_rap4_2014Largee[nRapBins_2010]={67.9,63.6};

float N3S_pp_rap3p5_2014[nRapBins_2013]={170,175,179,121,189};//,111,74.9};
float N3S_pp_rap3p5_2014e[nRapBins_2013]={21.2,23.3,25.7,24.9,28.2};//;,24.1,14.3};
float N3S_pp_rap4_2014[nRapBins_2014]={134,138,138,96,108,67.3}; //175}
float N3S_pp_rap4_2014e[nRapBins_2014]={18.1,20.7,20.5,20.7,37,14.5}; // 23.4

//may 6th from my fits
float N2S_pp_pt3p5_2013[nPtBins_2013] = {433,529,378,266,144};//31.2
float N2S_pp_pt3p5_2013e[nPtBins_2013] = {45,41,28,23,15};//7.6
float N2S_pp_pt4_2013[nPtBins_2013] = {318,302,274,230,131};//31.7
float N2S_pp_pt4_2013e[nPtBins_2013] = {34,30,22,21,15};//8.0
float N2S_pp_pt4_2013Large[nPtBins_2010]={827,287,200};
float N2S_pp_pt4_2013Largee[nPtBins_2010]={168,23,31};

float N3S_pp_pt4_2013[nPtBins_2013] = {141,199,134,108,75};//23.4
float N3S_pp_pt4_2013e[nPtBins_2013] = {27,28,18,16,12};//7.2
float N3S_pp_pt3p5_2013[nPtBins_2013] = {203,282,184,130,83};//21
float N3S_pp_pt3p5_2013e[nPtBins_2013] = {36,34,23,17,13};//6.5
//could add the large ones to 3S data too.

////yields in the note
//trying now with yields from plots, waiting for the Arleo ones. 25th apr. 2014
float N1S_pp_pt3p5[nPtBins_2013] = {1717,1550,1107,691,352}; // ,62.5
float N1S_pp_pt3p5e[nPtBins_2013] = {80,64,43,43,23}; // ,9.9
float N1S_pp_pt3p5U[nPtBins_2013] = {1718,1546.18,1109.26,692.5,350}; // ,62.8
float N1S_pp_pt3p5eU[nPtBins_2013]= {81,65,44,44,23}; // ,9.9,10
float N1S_pp_pt3p5Large[nPtBins_2010]={3856,910,586};
float N1S_pp_pt3p5eLarge[nPtBins_2010]={97,37,29};
float N1S_pp_rap3p5[nRapBins_2013] = {1102,844,795,1198,1427}; 
float N1S_pp_rap3p5e[nRapBins_2013] = {47,42,37,47,49};

float N2S_pp_pt3p5[nPtBins_2010] = {1166,335,227}; // for CS2S_pp_ptLarge
float N2S_pp_pt3p5e[nPtBins_2010]= {64,26,35};
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
float N1S_pp_tot3p5=5476;
float N1S_pp_tot3p5e=123;
float N2S_pp_tot3p5=1698;
float N2S_pp_tot3p5e=76;
float N2S_pp_tot4=1343;
float N2S_pp_tot4e=67;

float N3S_pp_tot3p5=841;
float N3S_pp_tot3p5e=64;

float N3S_pp_tot4= 700;
float N3S_pp_tot4e=55;


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
