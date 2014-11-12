#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>

#include <TH1.h>
#include <TH2D.h>

#include <TBranch.h>
#include <TCanvas.h>
#include "TClonesArray.h"
#include <TDirectory.h>
#include <TFile.h>
#include "TF1.h"
#include "TH1F.h"
#include <TLatex.h>
#include <TLegend.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "TMatrixT.h"
#include "TRandom.h"
#include <TStyle.h>
#include <TSystem.h>
#include "TTree.h"
#include "TString.h"
#include "TGraphErrors.h"

// miscellaneous  
#include <fstream>
#include <map>
#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#endif


void MatrixInversion(const char* inputOniaTree = "../../upsiMiniTree_Pyquen1S_QQtrigbit1_Trig_Unfolding_postCut_deltaRmatched_withCentrality.root") //upsiMiniTree_pythia1SwithCut.root //upsiMiniTree_Pyquen1S_QQtrigbit1_Trig_Unfolding_postCut_deltaRmatched_withCentrality.root
{
  gROOT->Macro("./logon.C+");
  double mass_min = 7.0; 
  double mass_max = 14.0;
  gStyle->SetOptStat(1);
  TFile *f = TFile::Open(Form("%s",inputOniaTree));

  TTree *t = (TTree*)f->Get("UpsilonTree");
  Long64_t nentries = t->GetEntries();


 
  const int nMCBins = 7;

  double filterEfficiency[nMCBins]={0.0554474,0.0612928,0.0238773,0.00887936,0.00368926,0.00303737,0.000169349};//to reweight the pt spectrum.
  double N_i[nMCBins]={172208,171048,172165,162777,117270,107560,107280};
  double NRec_i[nMCBins]={84185,83654,90231,97508,77935,79229,85703};//nb of Reco events for Gen pt bins (reco-events in bin [0,3] were generated in [0,3]).

  float GenUpsPt;
  float RecoUpsPt;
  float GenUpsM;
  float RecoUpsM;
  TBranch *b_GenUpsPt;
  TBranch *b_RecoUpsPt;
  TBranch *b_GenUpsM;
  TBranch *b_RecoUpsM;
  t->SetBranchAddress("GenUpsPt",      &GenUpsPt,      &b_GenUpsPt);
  t->SetBranchAddress("RecoUpsPt",      &RecoUpsPt,      &b_RecoUpsPt);
  t->SetBranchAddress("GenUpsM",      &GenUpsM,      &b_GenUpsM);
  t->SetBranchAddress("RecoUpsM",      &RecoUpsM,      &b_RecoUpsM);

  float GenUpsRap;
  float RecoUpsRap;
  float GenUpsPhi;
  float RecoUpsPhi;
  int centrality;
  TBranch *b_GenUpsRap;
  TBranch *b_RecoUpsRap;
  TBranch *b_GenUpsPhi;
  TBranch *b_RecoUpsPhi;
  TBranch *b_centrality;
  t->SetBranchAddress("GenUpsRap",      &GenUpsRap,      &b_GenUpsRap);
  t->SetBranchAddress("RecoUpsRap",      &RecoUpsRap,      &b_RecoUpsRap);
  t->SetBranchAddress("GenUpsPhi",      &GenUpsPhi,      &b_GenUpsPhi);
  t->SetBranchAddress("RecoUpsPhi",      &RecoUpsPhi,      &b_RecoUpsPhi);
  // t->SetBranchAddress("centrality",      &centrality,      &b_centrality);


  Int_t QQsize;
  double SumWeightsGen[7]={};
  double SumWeightsReco[7]={}; //sum d'une collone de matrice (sum des poids des events dans un bin de 0<GenUpsPt<3)
  double SumWeightsReco2[7]={};
  double SumWeightsReco2oui[7]={}; 
  double SumWeightsReco2non[7]={}; 
  double SumColM=0;
  double SumColI=0;


  double deltaM;
  double deltaR;
 
  TBranch *b_QQsize; //!
  t->SetBranchAddress("QQsize",  &QQsize, &b_QQsize); // 3eme argument : lieu d'origine de l'information dans le tree
 
  const int nBy=7;
  const int nBx=7;
  double binsx[nBx]={0,2.5,5,8,12,20,50};
  double binsy[nBy]={0,2.5,5,8,12,20,50};


  double jj;
  double ii;
  double IJ;
 

  const int nBX =15;
  const double wideX = 50;
  const int nBY = 200;
  const double wideY = 2;
  double nnBY = nBY;
  double nnBX = nBX;
  double binsY[nBY+1]={};
  for (int i=0; i<=nBY; i++)
    {
      ii=i;
      binsY[i]=-wideY/2 + i*wideY/nBY;
    }
 

  double distribSigma[nBx-1];
  double distribSigmaError[nBx-1];
  double distribMean[nBx-1];
  double distribMeanError[nBx-1];
  double fitSigma[nBx-1];
  double fitSigmaError[nBx-1];
  double fitMean[nBx-1];
  double fitMeanError[nBx-1];
  double xx[nBx-1];
  double yy[nBY];
  double xxError[nBx-1];
  double distribRMS[nBx-1];
  double SumWeightsHisto1DY[nBx-1]={};
 
  TMatrixT<double> N1S_obs_PbPb(1,nBy-1,1,1);
  N1S_obs_PbPb(1,1)=863;
  N1S_obs_PbPb(2,1)=929;
  N1S_obs_PbPb(3,1)=572;
  N1S_obs_PbPb(4,1)=346;
  N1S_obs_PbPb(5,1)=184;
  N1S_obs_PbPb(6,1)=53;
  TMatrixT<double> N1S_obs_PbPb_e(1,nBy-1,1,1);
  N1S_obs_PbPb_e(1,1)=92;
  N1S_obs_PbPb_e(2,1)=75;
  N1S_obs_PbPb_e(3,1)=66;
  N1S_obs_PbPb_e(4,1)=32;
  N1S_obs_PbPb_e(5,1)=21;
  N1S_obs_PbPb_e(6,1)=8.9;
  TMatrixT<double> N1S_unfolded(1,nBy-1,1,1);
  TMatrixT<double> N1S_unfolded_e(1,nBy-1,1,1);
  TMatrixT<double> N1S_obs_pp(1,nBy-1,1,1);
  N1S_obs_pp(1,1)=1717;
  N1S_obs_pp(2,1)=1550;
  N1S_obs_pp(3,1)=1107;
  N1S_obs_pp(4,1)=691;
  N1S_obs_pp(5,1)=352;
  N1S_obs_pp(6,1)=62.5;
  TMatrixT<double> N1S_obs_pp_e(1,nBy-1,1,1);
  N1S_obs_pp_e(1,1)=80;
  N1S_obs_pp_e(2,1)=64;
  N1S_obs_pp_e(3,1)=43;
  N1S_obs_pp_e(4,1)=43;
  N1S_obs_pp_e(5,1)=23;
  N1S_obs_pp_e(6,1)=9.9;
 
  TH1D *hPtGen_w = new TH1D("hPtGen_w",";p_{T} (GeV)", nBx-1,binsx); //nBx-1,binsx
  TH1D *hPtReco_w = new TH1D("hPtReco_w","", nBx-1,binsx);
  TH2D *hPtGenReco_w = new TH2D("hPtGenReco_w",";p_{T} GEN (GeV)  ;p_{T} RECO (GeV) ",nBx-1,binsx,nBy-1,binsy); // nBX,0,wideX,nBY,-wideY/2,wideY/2
  TMatrixT<double> TheMatrix(1,nBx-1,1,nBx-1);
  TMatrixT<double> TheInverse(1,nBx-1,1,nBx-1);
  TMatrixT<double> Id(1,nBx-1,1,nBx-1);
  TH2D *Unfolding = new TH2D("Unfolding Matrix",";p_{T}^{RECO} (GeV) ;p_{T}^{GEN} (GeV) ;fraction of p_{T}^{RECO}  ",nBx-1,binsx,nBy-1,binsy); //nBX,0,wideX,nBY,-wideY/2,wideY/2 //nBx-1,binsx,nBy-1,binsy

  TRandom3 r(0);
  double smear[404670];
   for(int i=0; i<nentries; i++)
     {
       smear[i]=0; //r.Gaus(0,0.1)
     }
   //cout<<smear[nentries-1]<<endl;


  cout << "--------------------------------------------------" << endl;
  cout << "filling the SumWeightsReco array for normalisation" << endl;
  cout << "--------------------------------------------------" << endl;
  for(int i=0; i<nentries ; i++)
    {
      //  if (i==2*((i)/2)) continue; // keep only events with odd number
      t->GetEntry(i);
      if (( RecoUpsPt+smear[i]>=50) || ( GenUpsPt>=50)) continue;
      if  (RecoUpsPt+smear[i]<0) continue;
      // deltaR=sqrt((RecoUpsRap-GenUpsRap)*(RecoUpsRap-GenUpsRap)+(RecoUpsPhi-GenUpsPhi)*(RecoUpsPhi-GenUpsPhi));
      // if ((deltaR>=0.1) ) continue; // deltaR cut 
      if (( RecoUpsM>=10.1) ||( RecoUpsM<=8.5) ) continue; // mass cut to keep the 'good' events
      // if (centrality<30) continue;

      if(GenUpsPt <2.5)
      	{
      	  SumWeightsReco[0]+= filterEfficiency[0]/ N_i[0];
  	  SumWeightsReco2[0]+= (( filterEfficiency[0]/ N_i[0]) * ( filterEfficiency[0]/ N_i[0]));
      	} 
      else if(GenUpsPt>=2.5 && GenUpsPt<3)
      	{
      	  SumWeightsReco[1]+= filterEfficiency[0]/ N_i[0];
  	  SumWeightsReco2[1]+=( filterEfficiency[0]/ N_i[0]) * ( filterEfficiency[0]/ N_i[0]);
      	} 
      else if(GenUpsPt >=3 && GenUpsPt< 5)
      	{
      	  SumWeightsReco[1]+= filterEfficiency[1]/ N_i[1];
  	  SumWeightsReco2[1]+=( filterEfficiency[1]/ N_i[1]) * ( filterEfficiency[1]/ N_i[1]);
      	}
      else if(GenUpsPt >=5 && GenUpsPt< 6)
      	{
      	  SumWeightsReco[2]+= filterEfficiency[1]/ N_i[1];
  	  SumWeightsReco2[2]+=( filterEfficiency[1]/ N_i[1]) * ( filterEfficiency[1]/ N_i[1]);
      	}
      else if(GenUpsPt >=6 && GenUpsPt< 8)
      	{
      	  SumWeightsReco[2]+= filterEfficiency[2]/ N_i[2];
  	  SumWeightsReco2[2]+=( filterEfficiency[2]/ N_i[2]) * ( filterEfficiency[2]/ N_i[2]);
      	}
      else if(GenUpsPt >=8 && GenUpsPt< 9)
      	{
      	  SumWeightsReco[3]+= filterEfficiency[2]/ N_i[2];
  	  SumWeightsReco2[3]+=( filterEfficiency[2]/ N_i[2]) * ( filterEfficiency[2]/ N_i[2]);
      	}
      else if(GenUpsPt >=9 && GenUpsPt< 12)
      	{
      	  SumWeightsReco[3]+= filterEfficiency[3]/ N_i[3];
  	  SumWeightsReco2[3]+=( filterEfficiency[3]/ N_i[3]) * ( filterEfficiency[3]/ N_i[3]);
      	}
      else if(GenUpsPt >=12 && GenUpsPt< 15)
      	{
      	  SumWeightsReco[4]+= filterEfficiency[4]/ N_i[4];
  	  SumWeightsReco2[4]+=( filterEfficiency[4]/ N_i[4]) * ( filterEfficiency[4]/ N_i[4]);
      	}
      else if(GenUpsPt >=15 && GenUpsPt< 20)  
      	{
      	  SumWeightsReco[4]+= filterEfficiency[5]/ N_i[5];
  	  SumWeightsReco2[4]+=( filterEfficiency[5]/ N_i[5]) * ( filterEfficiency[5]/ N_i[5]);
      	}
      else if(GenUpsPt >=20 && GenUpsPt< 30)
      	{
      	  SumWeightsReco[5]+= filterEfficiency[5]/ N_i[5];
  	  SumWeightsReco2[5]+=( filterEfficiency[5]/ N_i[5]) * ( filterEfficiency[5]/ N_i[5]);
      	}
      else if(GenUpsPt >=30)
      	{
      	  SumWeightsReco[5]+= filterEfficiency[6]/ N_i[6];
  	  SumWeightsReco2[5]+=( filterEfficiency[6]/ N_i[6]) * ( filterEfficiency[6]/ N_i[6]);
      	}
    }

  for (int i=1; i<=nBx-1; i++)
     {
  cout <<  SumWeightsReco2[i-1]<<"   "<<SumWeightsReco[i-1]<<"   "<<endl;
     }
  cout << "--------------------------------------------------" << endl;
  cout << "           Filling the Reco/Gen TH2D" << endl;
  cout << "--------------------------------------------------" << endl;
  //version matrice (avec renormalisation) PbPb
  for(int i=0; i<nentries; i++)
    {
      // if (i==2*((i)/2)) continue;
      t->GetEntry(i);
      if (( RecoUpsPt+smear[i]>=50) ||( GenUpsPt>=50) ) continue;
      if  (RecoUpsPt+smear[i]<0) continue;
        if (( RecoUpsM>=10.1) ||( RecoUpsM<=8.5) ) continue;
      //  if (centrality<30)  continue;
      // hMReco_w->Fill(RecoUpsM,1/404670.0);
      // deltaR=sqrt((RecoUpsRap-GenUpsRap)*(RecoUpsRap-GenUpsRap)+(RecoUpsPhi-GenUpsPhi)*(RecoUpsPhi-GenUpsPhi));
      //  if ((deltaR>=0.1) ) continue;

	//comment: Generation bin differ from analysis bin. So the structure of the loop is this:
	// depending on the Gen Pt of the pair:
	// 1. if we're the i-th Gen bin, fill the 1D histos for PtGen and PtReco with the i-th filter eff.
	// 2. if we're also in the i-th Reco bin, fill the 2D histo with PtGenReco weighted with the i-th filter eff and the i-th SumWeightsReco value.
	// 3. or if we're not in the i-th Reco bin, fill the 2D hist with PtGenReco weighted with the i-th filter eff and the j-th SumWeigthsReco value.
	// As a result, the rows are normalised (gen filter Eff) and the columns too (SumWeightsReco).
      if(GenUpsPt <3)
  	{
  	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[0]/(N_i[0]));
    	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[0]/N_i[0]);
  	  if (GenUpsPt<2.5)
  	    {
  	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[0]/( N_i[0]*SumWeightsReco[0]));
  	    }
  	  else
  	    {
  	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[0]/( N_i[0]*SumWeightsReco[1]));
  	    }
  	}
      else if(GenUpsPt>=3 && GenUpsPt< 6)
  	{
  	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[1]/(N_i[1]));
    	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[1]/N_i[1]);
	
  	  if (GenUpsPt<5)
  	    {
  	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[1]/( N_i[1]*SumWeightsReco[1]));
  	    }
  	  else
  	    {
  	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[1]/( N_i[1]*SumWeightsReco[2]));
  	    }
  	}
      else if(GenUpsPt >=6 && GenUpsPt< 9)
      	{
      	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[2]/(N_i[2]));
      	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[2]/N_i[2]);
      	  if (GenUpsPt<8)
      	    {
      	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[2]/( N_i[2]*SumWeightsReco[2]));
      	    }
      	  else
      	    {
      	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[2]/( N_i[2]*SumWeightsReco[3]));
      	    }
      	}
      else if(GenUpsPt >=9 && GenUpsPt< 12)
      	{
      	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[3]/(N_i[3]));
      	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[3]/N_i[3]);
      	  hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[3]/( N_i[3]*SumWeightsReco[3]));
      	}
         else if(GenUpsPt >=12 && GenUpsPt< 15)
  	{
  	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[4]/(N_i[4]));
  	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[4]/N_i[4]);
  	  hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[4]/( N_i[4]*SumWeightsReco[4]));

  	}
      else if(GenUpsPt >=15 && GenUpsPt< 30)
  	{
  	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[5]/(N_i[5]));
  	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[5]/N_i[5]);

  	  if (GenUpsPt<20)
  	    {
  	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[5]/( N_i[5]*SumWeightsReco[4]));
  	    }
  	  else
  	    {
  	      hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[5]/( N_i[5]*SumWeightsReco[5]));
  	    }
  	}
      else if(GenUpsPt >=30)
  	{
  	  hPtGen_w->Fill(GenUpsPt,filterEfficiency[6]/N_i[6]);
  	  hPtReco_w->Fill(RecoUpsPt+smear[i],filterEfficiency[6]/N_i[6]);
  	  hPtGenReco_w->Fill(GenUpsPt,RecoUpsPt+smear[i],  filterEfficiency[6]/( N_i[6]*SumWeightsReco[5]));
  	}    
    }
  // Drawing the 2D histo as a result.
  TCanvas *c02 = new TCanvas("c02","c02",1000,1000);
  c02->cd(); 
  TPad *p2 = new TPad("p2","p2",0,0,0.95,1);
  p2->SetRightMargin(0.2);
  p2->Draw();
  p2->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  gStyle->SetPaintTextFormat(".3f ");
  hPtGenReco_w->SetMarkerSize(0.9);
  hPtGenReco_w->GetZaxis()->SetRangeUser(0,0.1);
  // hPtGenReco_w->GetXaxis()->SetRangeUser(0,wideX);
  hPtGenReco_w->GetZaxis()->SetLabelSize(0.03);
  hPtGenReco_w->GetZaxis()->SetTitleOffset(1.25);
  hPtGenReco_w->Draw("TEXT00 colz");
  c02->Draw();
  //
  cout << "--------------------------------------------------" << endl;
  cout << "               Filling the TMatrixT" << endl;
  cout << "--------------------------------------------------" << endl;
  //contenu de l'histo2D => Matrice
  for (int i=1; i<=nBx-1; i++)
    {
      for (int j=1; j<=nBy-1; j++)
  	{
  	   TheMatrix(i,j)=hPtGenReco_w->GetBinContent(j,i); //GetBinContent takes the matrix with index (columns,rows) instead of (rows, columns) as usual (and starting from lower left-hand corner)
  	  // if (TheMatrix(i,j) < 1e-4)
  	  //   {
  	  //     TheMatrix(i,j)=0;
  	  //   }
  	}
    }

  for (int i=1; i<=nBx-1; i++)
    {
      cout<<TheMatrix(i,1)<< "    "<<TheMatrix(i,2)<< "     "<<TheMatrix(i,3)<< "    "<<TheMatrix(i,4)<< "    " <<TheMatrix(i,5)<< "    "<<TheMatrix(i,6)<< "    "<<endl;
    }
  cout << "--------------------------------------------------" << endl;
  cout << "                 Invert matrix:" << endl;
  cout << "--------------------------------------------------" << endl;

   //Inversion de la matrice
   TheInverse = TheMatrix;
   TheInverse = TheInverse.Invert();
   for (int i=1; i<=nBx-1; i++)
     {
       cout<<TheInverse(i,1)<< "    "<<TheInverse(i,2)<< "     "<<TheInverse(i,3)<< "    "<<TheInverse(i,4)<< "    " <<TheInverse(i,5)<< "    "<<TheInverse(i,6)<< "    "<<endl;
     }


   //inverse matrix content => histo2D
   for (int i=1; i<=5; i++)
     {
       for (int j=1; j<=5; j++)
   	{
   	  Unfolding->Fill(3*(i-1)+1e-4,3*(j-1)+1e-4,TheInverse(j,i)); //don't forget index are in reverse : (j,i) instead of (i,j)
   	}
     }
   for (int j=1; j<=5; j++)
     {
       Unfolding->Fill(3*(j-1)+1e-4,25,TheInverse(6,j));
       Unfolding->Fill(25,3*(j-1)+1e-4,TheInverse(j,6));
     }
   Unfolding->Fill(25,25,TheInverse(6,6));

  cout << "--------------------------------------------------" << endl;
  cout << "                 Inversion Check:" << endl;
  cout << "--------------------------------------------------" << endl;

   //Inversion check
  Id.Mult(TheMatrix,TheInverse);
   for (int i=1; i<=nBx-1; i++)
     {
       cout<<Id(i,1)<< "    "<<Id(i,2)<< "     "<<Id(i,3)<< "    "<<Id(i,4)<< "    " <<Id(i,5)<< "    "<<Id(i,6)<< "    "<<endl;
     }
   //Normalisation check
  cout << "--------------------------------------------------" << endl;
  cout << "               Normalisation Check:" << endl;
  cout << "--------------------------------------------------" << endl;
  for (int j=1; j<=nBx-1; j++)
     {
       SumColM=0;
       SumColI=0;
       for (int i=1; i<=nBx-1; i++)
  	{
  	  SumColM+=TheMatrix(i,j);
  	  SumColI+=TheInverse(i,j);
  	}
       cout<<"Matrix Column # " << j << " | " << SumColM << " " << endl;
       cout<<"Invert Column # " << j << " | " << SumColI << " " << endl;
     }

  cout << "--------------------------------------------------" << endl;
  cout << "          Apply the unfolding correction:" << endl;
  cout << "--------------------------------------------------" << endl;

  N1S_unfolded.Mult(TheInverse,N1S_obs_PbPb);
  //calculate the error 
  for (int i=1; i<=nBx-1; i++)
    {
      for (int j=1; j<=nBx-1; j++)
        {
  	 N1S_unfolded_e(i,1)+=TheInverse(i,j)*N1S_obs_PbPb_e(j,1)*TheInverse(i,j)*N1S_obs_PbPb_e(j,1);
        }
      N1S_unfolded_e(i,1)=sqrt(N1S_unfolded_e(i,1));
    }

  for (int i=1; i<=nBx-1; i++)
     {
       cout<<N1S_unfolded(i,1)<<" +/- "<<N1S_unfolded_e(i,1)<<" , obs:   "<<N1S_obs_PbPb(i,1)<<" +/- "<<N1S_obs_PbPb_e(i,1)<< endl;
     }

  TCanvas *c03 = new TCanvas("c03","c03",1000,1000);
  c03->cd(); 
  TPad *p3 = new TPad("p3","p3",0,0,1,1);
  p3->SetRightMargin(0.20);
  p3->Draw();
  p3->cd();
  gPad->SetLogx();
  gPad->SetLogy();
  gStyle->SetPaintTextFormat(".3f ");
  Unfolding->SetMarkerSize(0.9);
  Unfolding->GetZaxis()->SetRangeUser(-0.05,0.2);
  // Unfolding->GetXaxis()->SetRangeUser(0,wideX);
  Unfolding->GetZaxis()->SetLabelSize(0.03);
  Unfolding->GetZaxis()->SetTitleOffset(1.25);
  Unfolding->Draw("TEXT00 colz");
}
