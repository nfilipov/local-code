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

const bool plotCS =true; //unused for now.
const bool plotRAA=true;// centrality, transverse momentum, rapidity.
const bool plotScaledToMinBias=true;// need to scale the error on x to its value at x'.
const bool plot2010=true; //comparison with published 2010 results.
float computeRatio(float x, float y) ;
float computeRatioError(float x, float y, float xerr, float yerr);
float scaleToMB(float RAA_MB, float RAA_new, float RAA_new_err);
void plotRaa2014ForChad()
{
  gROOT->Macro("logon.C+");//it all looks much nicer with this.
  double RapBinWidth = 4.8;
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
  float CS1S_pp_tnp_rap2014[nRapBins_2014] = {};
  float CS1S_pp_tnp_rap2014e[nRapBins_2014] = {};
  float CS1S_pp_tnp_rap2014s[nRapBins_2014] = {};
  float CS1S_aa_tnp_rap2014[nRapBins_2014] = {};
  float CS1S_aa_tnp_rap2014e[nRapBins_2014] = {};
  float CS1S_aa_tnp_rap2014s[nRapBins_2014] = {};
  float CS1S_pp_tnp_pt[nPtBins_2013] = {};
  float CS1S_pp_tnp_pte[nPtBins_2013] = {};
  float CS1S_pp_tnp_pts[nPtBins_2013] = {};
  float CS1S_aa_tnp_pt[nPtBins_2013] = {};
  float CS1S_aa_tnp_pts[nPtBins_2013] = {};
  float CS1S_aa_tnp_pte[nPtBins_2013] = {};

  float RAA_1S_tnp_pt[nPtBins_2013]={};
  float RAA_1S_tnp_rap[nRapBins_2014]={};
  float RAA_1S_tnp_pte[nPtBins_2013]={};
  float RAA_1S_tnp_rape[nRapBins_2014]={};
  
  for(int i = 0 ; i<nPtBins_2013 ; i++)
    {
      CS1S_pp_tnp_pt[i]= computeRatio( N1S_pp_pt3p5[i] , Aet_1S_pythia_pt[i] ); 
      CS1S_aa_tnp_pt[i]= computeRatio( N1S_aa_pt3p5[i] , Aet_1S_pyquen_pt[i] );
      CS1S_pp_tnp_pte[i] = computeRatioError( N1S_pp_pt3p5[i] , Aet_1S_pythia_pt[i], N1S_pp_pt3p5e[i] , Aet_1S_pythia_pte[i]);
      CS1S_aa_tnp_pte[i] = computeRatioError(  N1S_aa_pt3p5[i] , Aet_1S_pyquen_pt[i] ,  N1S_aa_pt3p5e[i] , Aet_1S_pyquen_pte[i] );
      CS1S_pp_tnp_pt[i]=CS1S_pp_tnp_pt[i]/(L_pp_invNb*RapBinWidth*deltaPt[i]);    
      CS1S_aa_tnp_pt[i]=CS1S_aa_tnp_pt[i]/(N_MB_corr * T_AA_b*(RapBinWidth*deltaPt[i]));
      CS1S_pp_tnp_pte[i]=CS1S_pp_tnp_pte[i]/(L_pp_invNb*RapBinWidth*deltaPt[i]);
      CS1S_aa_tnp_pte[i]=CS1S_aa_tnp_pte[i]/(N_MB_corr * T_AA_b *(RapBinWidth*deltaPt[i]));
      RAA_1S_tnp_pt[i] = computeRatio( CS1S_aa_tnp_pt[i] , CS1S_pp_tnp_pt[i]);
      RAA_1S_tnp_pte[i] = computeRatioError( CS1S_aa_tnp_pt[i] , CS1S_pp_tnp_pt[i], CS1S_aa_tnp_pte[i] , CS1S_pp_tnp_pte[i]);
      if(plotScaledToMinBias){
	RAA_1S_tnp_pt[i] = 0.56;
	RAA_1S_tnp_pte[i]=scaleToMB(0.56,RAA_1S_tnp_pt[i],RAA_1S_tnp_pte[i]);
      }
    }
 
  for(int i=0; i < nRapBins_2014; i++)
    {
      CS1S_pp_tnp_rap2014[i]= computeRatio( N1S_pp_rap3p5_2014[i] , Aet_1S_pythia_rap2014[i] ); 
      CS1S_aa_tnp_rap2014[i]= computeRatio( N1S_aa_rap3p5_2014[i] , Aet_1S_pyquen_rap2014[i] );
      CS1S_pp_tnp_rap2014e[i] = computeRatioError(  N1S_pp_rap3p5_2014[i] , Aet_1S_pythia_rap2014[i] ,  N1S_pp_rap3p5_2014e[i] , Aet_1S_pythia_rap2014e[i] );
      CS1S_aa_tnp_rap2014e[i] = computeRatioError(N1S_aa_rap3p5_2014[i] , Aet_1S_pyquen_rap2014[i] ,  N1S_aa_rap3p5_2014e[i] , Aet_1S_pyquen_rap2014e[i]); 
      CS1S_pp_tnp_rap2014[i]=CS1S_pp_tnp_rap2014[i]/L_pp_invNb;
      CS1S_aa_tnp_rap2014[i]=CS1S_aa_tnp_rap2014[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_tnp_rap2014e[i]=CS1S_pp_tnp_rap2014e[i]/L_pp_invNb;
      CS1S_aa_tnp_rap2014e[i]=CS1S_aa_tnp_rap2014e[i]/(N_MB_corr * T_AA_b);
      CS1S_pp_tnp_rap2014[i]=CS1S_pp_tnp_rap2014[i]/(deltaRapEven[i]);
      CS1S_aa_tnp_rap2014[i]=CS1S_aa_tnp_rap2014[i]/deltaRapEven[i];
      CS1S_pp_tnp_rap2014e[i]=CS1S_pp_tnp_rap2014e[i]/(deltaRapEven[i]);
      CS1S_aa_tnp_rap2014e[i]=CS1S_aa_tnp_rap2014e[i]/deltaRapEven[i];
      RAA_1S_tnp_rap[i]= computeRatio( CS1S_aa_tnp_rap2014[i] , CS1S_pp_tnp_rap2014[i]);
      RAA_1S_tnp_rape[i]= computeRatioError( CS1S_aa_tnp_rap2014[i] , CS1S_pp_tnp_rap2014[i],  CS1S_aa_tnp_rap2014e[i] , CS1S_pp_tnp_rap2014e[i]);
      if(plotScaledToMinBias){
	RAA_1S_tnp_rap[i] = 0.56;
	RAA_1S_tnp_rape[i]=scaleToMB(0.56,RAA_1S_tnp_rap[i],RAA_1S_tnp_rape[i]);
      }
    }

cout << "  --- 1S Cross section in pp vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_pp = "<< CS1S_pp_tnp_rap2014[j] <<" +/- "<<CS1S_pp_tnp_rap2014e[j]  <<" nb" << endl;
   }

cout << "  --- 1S Cross section in pp vs. pt ---" << endl;
 for(int j =0 ; j<nRapBins_2013 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_pp = "<< CS1S_pp_tnp_pt[j] <<" +/- "<<CS1S_pp_tnp_pte[j]  <<" nb" << endl;
   }


cout << "  --- 1S Cross section in PbPb vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_PbPb = "<< CS1S_aa_tnp_rap2014[j] <<" +/- "<<CS1S_aa_tnp_rap2014e[j] <<" nb" << endl;
   }

cout << "  --- 1S Cross section in PbPb vs. pt ---" << endl;
 for(int j =0 ; j<nPtBins_2013 ; j++)
   {
     cout <<"j="<< j << "' ,sigma(1S)_PbPb = "<< CS1S_aa_tnp_pt[j] <<" +/- "<<CS1S_aa_tnp_pte[j]  <<" nb" << endl;
   }

cout << "  --- 1S RAA vs. p_{T} ---" << endl;
 for(int j =0 ; j<nPtBins_2013 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_1S_tnp_pt[j] <<" +/- "<< RAA_1S_tnp_pte[j] <<  endl;
   }

cout << "  --- 1S RAA vs. y ---" << endl;
 for(int j =0 ; j<nRapBins_2014 ; j++)
   {
     cout <<"j="<< j << "' , Raa = "<< RAA_1S_tnp_rap[j] <<" +/- "<< RAA_1S_tnp_rape[j]<< endl;
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
   //one pad to draw Raa vs Rap!
   TF1 *f4RaaRap = new TF1("f4RaaRap","1",0,2.45);
   f4RaaRap->SetLineWidth(0);
   f4RaaRap->GetXaxis()->SetTitle("|y^{#Upsilon}|");
   f4RaaRap->GetYaxis()->SetTitle("R_{AA}");
   f4RaaRap->GetYaxis()->SetTitleOffset(1.5);
   f4RaaRap->GetYaxis()->SetTitleSize(0.04);
   f4RaaRap->GetYaxis()->SetRangeUser(0.,1.6);
   f4RaaRap->GetXaxis()->CenterTitle(kTRUE);
   f4RaaRap->Draw();
   TLegend *legend = new TLegend(0.2,0.5,0.4,0.6);
   if(plot2010){
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
     legend->AddEntry(grap2010,"#varUpsilon(1S) JHEP 05 (2012) 063","lp");
     f4RaaRap->Draw("same");
   }
   // TGraphErrors *gRaaRap1syst = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_tnp_rap,rap2014e,RAA_1S_tnp_raps);
   // gRaaRap1syst->SetLineColor(kOrange+1);
   // gRaaRap1syst->SetFillStyle(0);
   // gRaaRap1syst->SetLineWidth(2);
   // gRaaRap1syst->SetMarkerSize(0);
   // gRaaRap1syst->Draw("2");
   TGraphErrors *gRaaRap1TNP = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_tnp_rap,0,RAA_1S_tnp_rape);
   gRaaRap1TNP->SetMarkerColor(kOrange+1);
   gRaaRap1TNP->SetMarkerStyle(21);
   gRaaRap1TNP->SetMarkerSize(1);
   TGraphErrors *gRaaRap1TNPcircle = new TGraphErrors(nRapBins_2014,rap2014,RAA_1S_tnp_rap,0,RAA_1S_tnp_rape);
   gRaaRap1TNPcircle->SetMarkerStyle(25);
   gRaaRap1TNPcircle->SetMarkerSize(1);
   gRaaRap1TNPcircle->SetLineColor(kBlack);
   gRaaRap1TNP->Draw("pe");
   gRaaRap1TNPcircle->Draw("p");
   f4RaaRap->Draw("same");
   gPad->RedrawAxis();

   // TGraphErrors *gRaaRap2TNP = new TGraphErrors(nRapBins_2010,rap2010,RAA_2S_tnp_rap,0,RAA_2S_tnp_rape);
   // gRaaRap2TNP->SetMarkerColor(kOrange+4);
   // gRaaRap2TNP->SetMarkerStyle(20);
   // gRaaRap2TNP->SetMarkerSize(1);
   // TGraphErrors *gRaaRap2TNPcircle = new TGraphErrors(nRapBins_2010,rap2010,RAA_2S_tnp_rap,0,RAA_2S_tnp_rape);
   // gRaaRap2TNPcircle->SetMarkerStyle(24);
   // gRaaRap2TNPcircle->SetMarkerSize(1);
   // gRaaRap2TNPcircle->SetLineColor(kBlack);
   // TGraphErrors *gRaaRap2syst = new TGraphErrors(nRapBins_2010,rap2010,RAA_2S_tnp_rap,rap2010e,RAA_2S_tnp_raps);
   // gRaaRap2syst->SetLineColor(kOrange+4);
   // gRaaRap2syst->SetFillStyle(0);
   // gRaaRap2syst->SetLineWidth(2);
   // gRaaRap2syst->SetMarkerSize(0);
   // gRaaRap2syst->Draw("2");
   // gRaaRap2TNP->Draw("pe");
   // gRaaRap2TNPcircle->Draw("p");
   // f4RaaRap->Draw("same");
   // gPad->RedrawAxis();
  
   TLatex *l1CMSrap = new TLatex(0.2,1.45, "CMS Internal #sqrt{s_{NN}} = 2.76 TeV");
   l1CMSrap->SetTextFont(42);
   l1CMSrap->SetTextSize(0.038);
   l1CMSrap->Draw();
   TLatex *lyLRAP= new TLatex(0.2,1.3,"2011 L_{int}^{PbPb} = 150 #mub^{-1}; 2013 L_{int}^{pp} = 5.4 pb^{-1};");
   lyLRAP->SetTextFont(42);
   lyLRAP->SetTextSize(0.027);
   lyLRAP->Draw();
   if(plot2010) {lyLRAP->DrawLatex(0.2,1.15,"2010 L_{int}^{PbPb} = 7.28 #mub^{-1}; L_{int}^{pp} = 225 nb^{-1};");}
   prap2->Update();
  
   legend->SetTextSize(0.029);
   legend->SetFillStyle(0);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->SetTextFont(42);
  
   gPad->RedrawAxis(); 

 
   //  legend->AddEntry(gRaaRap1,"#varUpsilon(1S)","lp");
   legend->AddEntry(gRaaRap1TNP,"#varUpsilon(1S), R_{AA}(2011) = 0.56","lp");
   // legend->AddEntry(gRaaRap2TNP,"#varUpsilon(2S)","lp");
   legend->Draw();

   gPad->RedrawAxis();
   cRaarap->SaveAs("~/Desktop/RAA_Rap.png");

  ///raa vs pt

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
 f4RaaPt->GetYaxis()->SetTitleOffset(1.5);
 f4RaaPt->GetYaxis()->SetTitleSize(0.04);
 f4RaaPt->GetYaxis()->SetRangeUser(0.,2);
 f4RaaPt->GetXaxis()->CenterTitle(kTRUE);
 f4RaaPt->Draw();
 TLegend *legend = new TLegend(0.2,0.65,0.4,0.75);
 if(plot2010){
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
   
   legend->AddEntry(gpt2010,"#varUpsilon(1S) JHEP 05 (2012) 063","lp");
 }
 TGraphErrors *gRaaPt1TNP = new TGraphErrors(nPtBins_2013,pt,RAA_1S_tnp_pt,pte,RAA_1S_tnp_pte);
 gRaaPt1TNP->SetMarkerColor(kOrange+1);
 gRaaPt1TNP->SetMarkerStyle(21);
 gRaaPt1TNP->SetMarkerSize(1);
 TGraphErrors *gRaaPt1TNPcircle = new TGraphErrors(nPtBins_2013,pt,RAA_1S_tnp_pt,pte,RAA_1S_tnp_pte);
 gRaaPt1TNPcircle->SetMarkerStyle(25);
 gRaaPt1TNPcircle->SetMarkerSize(1);
 gRaaPt1TNPcircle->SetLineColor(kBlack);
 gRaaPt1TNP->Draw("pe");
 gRaaPt1TNPcircle->Draw("p");
 f4RaaPt->Draw("same");
 
 // TGraphErrors *gRaaPt1syst = new TGraphErrors(nPtBins_2013,pt,RAA_1S_tnp_pt,pte,RAA_1S_tnp_pts);
 // gRaaPt1syst->SetLineColor(kOrange+1);
 // gRaaPt1syst->SetFillStyle(0);
 // gRaaPt1syst->SetLineWidth(2);
 // gRaaPt1syst->SetMarkerSize(0);
 // gRaaPt1syst->Draw("2");
 
 
 // TGraphErrors *gRaaPt2TNP = new TGraphErrors(nPtBins_2010,pt2010,RAA_2S_tnp_pt,pt2010e,RAA_2S_tnp_pte);
 // gRaaPt2TNP->SetMarkerColor(kOrange+4);
 // gRaaPt2TNP->SetMarkerStyle(22);
 // gRaaPt2TNP->SetMarkerSize(1);
 // TGraphErrors *gRaaPt2TNPcircle = new TGraphErrors(nPtBins_2010,pt2010,RAA_2S_tnp_pt,pt2010e,RAA_2S_tnp_pte);
 // gRaaPt2TNPcircle->SetMarkerStyle(26);
 // gRaaPt2TNPcircle->SetMarkerSize(1);
 // gRaaPt2TNPcircle->SetLineColor(kBlack);
 // gRaaPt2TNP->Draw("pe");
 // gRaaPt2TNPcircle->Draw("p");
 // f4RaaPt->Draw("same");
 // TGraphErrors *gRaaPt2syst = new TGraphErrors(nRapBins_2010,pt2010,RAA_2S_tnp_pt,pt2010e,RAA_2S_tnp_pts);
 // gRaaPt2syst->SetLineColor(kOrange+4);
 // gRaaPt2syst->SetFillStyle(0);
 // gRaaPt2syst->SetLineWidth(2);
 // gRaaPt2syst->SetMarkerSize(0);
 // gRaaPt2syst->Draw("2");

 TLatex *l1CMSpt = new TLatex(1,1.8, "CMS Internal  #sqrt{s_{NN}} = 2.76 TeV");
 l1CMSpt->SetTextFont(42);
 l1CMSpt->SetTextSize(0.04);
 l1CMSpt->Draw();
 TLatex *lyLPT= new TLatex(1,1.62,"2011 L_{int}^{PbPb} = 150 #mub^{-1}; 2013 L_{int}^{pp} = 5.4 pb^{-1};");
 lyLPT->SetTextFont(42);
 lyLPT->SetTextSize(0.027);
 lyLPT->Draw();
 if(plot2010){ lyLPT->DrawLatex(1,1.52,"2010 L_{int}^{PbPb} = 7.28 #mub^{-1}; L_{int}^{pp} = 225 nb^{-1};");}
 ppt2->Update();
 legend->SetTextSize(0.029);
 legend->SetFillStyle(0);
 legend->SetFillColor(0);
 legend->SetBorderSize(0);
 legend->SetTextFont(42);
 

 legend->AddEntry(gRaaPt1TNP,"#varUpsilon(1S), R_{AA}(2011) = 0.56","lp");
 //legend->AddEntry(gRaaPt2TNP,"#varUpsilon(2S)","lp");
 legend->Draw();
 gPad->RedrawAxis();
 cRaapt->SaveAs("~/Desktop/RAA_Pt.png");
 }

}

float computeRatio(float x, float y) 
{
  // pass the yield (x), and Acc*eff (y), and computes the corrected yield. then divide by lumi and delta rapidity to get the cross section. in case of pbpb, divide by taa*nMB to get the nColl scaled invariant yield.
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
float scaleToMB(float RAA_MB, float RAA_new, float RAA_new_err)
{
  return RAA_new_err*computeRatio(RAA_new,RAA_MB);
}
