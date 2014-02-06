// #include <TROOT.h>
// #include <TH1.h>
// #include <TH2D.h>
// #include <TBranch.h>
// #include <TCanvas.h>
// #include "TClonesArray.h"
// #include <TDirectory.h>
// #include <TFile.h>
// #include "TH1.h"
// #include <TLatex.h>
// #include <TLegend.h>
// #include "TLorentzVector.h"
// #include <TMath.h>
// #include "TRandom.h"
// #include <TStyle.h>
// #include "TPaveLabel.h"
// #include "TVirtualFitter.h"
// #include <TSystem.h>
// #include "TTree.h"
// #include "TString.h"
// #include "TChain.h"
// #include "TCut.h"
// #include "DataParameters.h"

// using namespace std;
// void Parameters()
{
  bool do3p5=true;
  bool do4 = false;
  bool doPt =true;
  bool doRap =true;
  gROOT->Macro("../code/cm/logon.C+");//it all looks much nicer with this.
  int nPtBins_2013 = 5;
  int nPtBins_2010 = 3;
  int nCentBins_2013 =8;
  int nRapBins_2013= 5;
  int nRapBins_2010= 2;
  double RapBinWidth = 4.8;
  float pt [nPtBins_2013] = {1.25, 3.75, 6.5, 10., 16.};
  float pte[nPtBins_2013] = {1.25, 1.25, 1.5, 2., 4.};
  float deltaPt[nPtBins_2013]   = {2.5,2.5,3,4,8};
  float pt_2010 [nPtBins_2010] = {3.25,8.25,15.};
  float pte_2010[nPtBins_2010] = {3.25,1.75,5};
  float deltaPt_2010[nPtBins_2010]   = {6.5,3.5,10};
  float deltaRap[nPtBins_2013]  = {0.8,0.6,0.6,1,1.8};
  float deltaRap2010[nRapBins_2010]={2.4,2.4};
  float nMB =1161498260; // 1.1E09 / 0.97 selection MB trigger eff
  float raapt[nPtBins_2013]={,,,,};
  float raarap[nRapBins_2013]={,,,,};
  float raapte[nPtBins_2013]={,,,,};
  float raarape[nRapBins_2013]={,,,,};
  float signal1S_pp_scaled_pt[nPtBins_2013];
  float signal1S_pp_scaled_rap[nRapBins_2013];
  float signal1S_aa_scaled_pt[nPtBins_2013];
  float signal1S_aa_scaled_rap[nRapBins_2013];
  float signal1S_pp_scaled_pte[nPtBins_2013];
  float signal1S_pp_scaled_rape[nRapBins_2013];
  float signal1S_aa_scaled_pte[nPtBins_2013];
  float signal1S_aa_scaled_rape[nRapBins_2013];
  float signal1S_aa_scaled_npart[nCentBins_2013];
  float signal1S_aa_scaled_nparte[nCentBins_2013];
  float lumi_pp=5300000000000;
  float taa=5660;
  float raapt3p5[nPtBins_2013];
  float raapte3p5[nPtBins_2013];
  float raarap3p5[nRapBins_2013];
  float raarape3p5[nRapBins_2013];
  
  float raapt4[nPtBins_2013];
  float raapte4[nPtBins_2013];
  float raarap4[nRapBins_2013];
  float raarape4[nRapBins_2013];
  
  //funny variables to scale raa3p5 errors to raa4 values!
  float raapte3p5_scaled[nPtBins_2013];
  float raarape3p5_scaled[nRapBins_2013];
  float RapForScaling[nRapBins_2013];
  float PtForScaling[nPtBins_2013];
  float taa2013[8]={25.90148438,
			20.47,
			14.47769531,
			8.782976563,
			5.089226563,
			2.748355469,
			0.982564453,
			0.125328724};
  float fMB[8]={0.05,0.05,0.1,0.1,0.1,0.1,0.2,0.3};
 
  for(int l =0;l<nPtBins_2013;l++)
    {
      if(do3p5){
	gROOT->Macro("DataParameters.h");
	raapt[l]=(signal1S_aa_pt[l]/(nMB*taa))/(signal1S_pp_pt[l]/lumi_pp);
	raapte[l]=raapt[l]*sqrt((signal1S_aa_pte[l]/signal1S_aa_pt[l])*(signal1S_aa_pte[l]/signal1S_aa_pt[l]) + (signal1S_aa_pte[l]/signal1S_aa_pt[l])*(signal1S_pp_pte[l]/signal1S_pp_pt[l]));
	raapt3p5[l]=raapt[l];
	raapte3p5[l]=raapte[l];
     
	raarap[l]=(signal1S_aa_rap[l]/(nMB*taa))/(signal1S_pp_rap[l]/lumi_pp);
	raarape[l]=raarap[l]*sqrt((signal1S_aa_rape[l]/signal1S_aa_rap[l])*(signal1S_aa_rape[l]/signal1S_aa_rap[l])+ (signal1S_aa_rape[l]/signal1S_aa_rap[l])*(signal1S_pp_rape[l]/signal1S_pp_rap[l]));
	raarap3p5[l]=raarap[l];
	raarape3p5[l]=raarape[l];
     
	signal1S_pp_scaled_pt[l]=signal1S_pp_pt[l]/(lumi_pp*deltaPt[l]);
	signal1S_pp_scaled_rap[l]=signal1S_pp_rap[l]/(lumi_pp*deltaRap[l]);    
	signal1S_aa_scaled_pt[l]=signal1S_aa_pt[l]/(taa*nMB*deltaPt[l]);
	signal1S_aa_scaled_rap[l]=signal1S_aa_rap[l]/(taa*nMB*deltaRap[l]);

	signal1S_pp_scaled_pte[l]=signal1S_pp_pte[l]/(lumi_pp*deltaPt[l]);
	signal1S_pp_scaled_rape[l]=signal1S_pp_rape[l]/(lumi_pp*deltaRap[l]);    
	signal1S_aa_scaled_pte[l]=signal1S_aa_pte[l]/(taa*nMB*deltaPt[l]);
	signal1S_aa_scaled_rape[l]=signal1S_aa_rape[l]/(taa*nMB*deltaRap[l]);
	cout << "3p5: "<< raapt[l] << " \pm " << raapte[l] << endl; 
	cout << "3p5: "<< raarap[l] << " \pm " << raarape[l] << endl; 



	do4 =true;
      }    
      do3p5=false;
      if(do4){
	  gROOT->Macro("DataParameters.h");
	raapt[l]=(signal1S_aa_pt[l]/(nMB*taa))/(signal1S_pp_pt[l]/lumi_pp);
	raapte[l]=raapt[l]*sqrt((signal1S_aa_pte[l]/signal1S_aa_pt[l])*(signal1S_aa_pte[l]/signal1S_aa_pt[l]) + (signal1S_aa_pte[l]/signal1S_aa_pt[l])*(signal1S_pp_pte[l]/signal1S_pp_pt[l]));
	raapt4[l]=raapt[l];
	raapte4[l]=raapte[l];
     
	raarap[l]=(signal1S_aa_rap[l]/(nMB*taa))/(signal1S_pp_rap[l]/lumi_pp);
	raarape[l]=raarap[l]*sqrt((signal1S_aa_rape[l]/signal1S_aa_rap[l])*(signal1S_aa_rape[l]/signal1S_aa_rap[l])+ (signal1S_aa_rape[l]/signal1S_aa_rap[l])*(signal1S_pp_rape[l]/signal1S_pp_rap[l]));
	raarap4[l]=raarap[l];
	raarape4[l]=raarape[l];
     
	signal1S_pp_scaled_pt[l]=signal1S_pp_pt[l]/(lumi_pp*deltaPt[l]);
	signal1S_pp_scaled_rap[l]=signal1S_pp_rap[l]/(lumi_pp*deltaRap[l]);    
	signal1S_aa_scaled_pt[l]=signal1S_aa_pt[l]/(taa*nMB*deltaPt[l]);
	signal1S_aa_scaled_rap[l]=signal1S_aa_rap[l]/(taa*nMB*deltaRap[l]);

	signal1S_pp_scaled_pte[l]=signal1S_pp_pte[l]/(lumi_pp*deltaPt[l]);
	signal1S_pp_scaled_rape[l]=signal1S_pp_rape[l]/(lumi_pp*deltaRap[l]);    
	signal1S_aa_scaled_pte[l]=signal1S_aa_pte[l]/(taa*nMB*deltaPt[l]);
	signal1S_aa_scaled_rape[l]=signal1S_aa_rape[l]/(taa*nMB*deltaRap[l]);
	cout << "4: "<< raapt[l] << " \pm " << raapte[l] << endl; 
	cout << "4: "<< raarap[l] << " \pm " << raarape[l] << endl; 

	///skip this if not playing with 3p5 and 4...
	raapte3p5_scaled[l]=raapte3p5[l]*(raapt4[l]/raapt3p5[l]);
	raarape3p5_scaled[l]=raarape3p5[l]*(raarap4[l]/raarap3p5[l]);
	PtForScaling[l]=pt[l]+0.3;
	RapForScaling[l]=rap[l]+0.05;
	do3p5=true;
      }
      do4=false;
    }

      TCanvas *cpt = new TCanvas("cpt","cpt"); 
      TCanvas *cpt_raa = new TCanvas("cpt_raa","cpt_raa"); 
      TCanvas *cpt2 = new TCanvas("cpt2","cpt2"); 
      TCanvas *cpt_raa2 = new TCanvas("cpt_raa2","cpt_raa2"); 
      TCanvas *crap = new TCanvas("crap","crap"); 
      TCanvas *crap_raa2 = new TCanvas("crap_raa2","crap_raa2"); 
      TCanvas *crap2 = new TCanvas("crap2","crap2"); 
      TCanvas *crap_raa = new TCanvas("crap_raa","crap_raa"); 
      TCanvas *cCBPt1= new TCanvas("cCBPt1","cCBPt1");
      TCanvas *cCBPt2= new TCanvas("cCBPt2","cCBPt2");
      TCanvas *cCBRap1 = new TCanvas("cCBRap1","cCBRap1");
      TCanvas *cCBRap2 = new TCanvas("cCBRap2","cCBRap2");
  
      
      TPad *ppt1 = new TPad("ppt1","ppt1",0.0,0.0,1.0,1.0);
      ppt1->SetBottomMargin(0.12);
      ppt1->SetTopMargin(0.08);
      ppt1->SetRightMargin(0.06);
      ppt1->SetLeftMargin(0.11);

      TF1 *f4Pt = new TF1("f4Pt","0",0,20.5);
      f4Pt->SetLineWidth(0);
      f4Pt->GetXaxis()->SetTitle("p_{T}^{#Upsilon_{cand.}} (GeV/c)");
      f4Pt->GetYaxis()->SetTitle("Mass peak resolution (GeV/c)");
      f4Pt->GetYaxis()->SetTitleOffset(1.6);
      f4Pt->GetYaxis()->SetLabelSize(0.032);
      f4Pt->GetYaxis()->SetTitleSize(0.03);
      f4Pt->GetXaxis()->CenterTitle(kTRUE);


      TPad *prap1 = new TPad("prap1","prap1",0.0,0.0,1.0,1.0);
      prap1->SetBottomMargin(0.12);
      prap1->SetTopMargin(0.08);
      prap1->SetRightMargin(0.06);
      prap1->SetLeftMargin(0.11);
      TF1 *f4Rap = new TF1("f4Rap","0",0,2.5);
      f4Rap->SetLineWidth(0);
      f4Rap->GetXaxis()->SetTitle("|y|^{#Upsilon_{cand.}}");
      f4Rap->GetYaxis()->SetTitle("Mass peak resolution (GeV/c)");
      f4Rap->GetYaxis()->SetTitleOffset(1.6);
      f4Rap->GetYaxis()->SetLabelSize(0.032);
      f4Rap->GetYaxis()->SetTitleSize(0.03);
      f4Rap->GetYaxis()->SetRangeUser(0.0,0.22);
      f4Rap->GetXaxis()->CenterTitle(kTRUE);
   
      //one pad to draw 
      if(do3p5)
	{
	 gROOT->Macro("DataParameters.h");
	  if(doPt){
	    cpt->cd();
	    ppt1->Draw();
	    ppt1->cd();
	      f4Pt->Draw();
	      f4Pt->GetYaxis()->SetRangeUser(0.0,0.22);
	    TGraphErrors *gMassResPtMC = new TGraphErrors(nPtBins_2013,pt,massRes1_MC_pt,pte,massRes1_MC_pte);
	    gMassResPtMC->SetMarkerColor(kGreen+2);
	    gMassResPtMC->SetMarkerStyle(33);
	    gMassResPtMC->SetMarkerSize(2);
	    TGraphErrors *gMassResPtMCcircle = new TGraphErrors(nPtBins_2013,pt,massRes1_MC_pt,pte,massRes1_MC_pte);
	    gMassResPtMCcircle->SetMarkerStyle(27);
	    gMassResPtMCcircle->SetMarkerSize(2);
	    gMassResPtMCcircle->SetLineColor(kBlack);
	    TGraphErrors *gMassResPtMC_2 = new TGraphErrors(nPtBins_2013,pt,massRes2_MC_pt,pte,massRes2_MC_pte);
	    gMassResPtMC_2->SetMarkerStyle(27);
	    gMassResPtMC_2->SetMarkerSize(2);
	    gMassResPtMC_2->SetLineColor(kBlack);

	    gMassResPtMCcircle->Draw("pe");
	    gMassResPtMC->Draw("p");
	    gMassResPtMC_2->Draw("p");
	    f4Pt->Draw("same");
	    TGraphErrors *gNoLabel = new TGraphErrors(bin,centnoErr,centnoErr,centnoErr,centnoErr);
	    gNoLabel->SetMarkerStyle(27);
	    gNoLabel->SetMarkerSize(2);
	    gNoLabel->SetLineColor(kWhite);
	    gNoLabel->SetFillColor(kWhite);
	  
	    TLegend *legendMuonCuts = new TLegend(0.1,0.8,0.5,0.9);
	    legendMuonCuts->SetTextSize(0.029);
	    legendMuonCuts->SetFillStyle(0);
	    legendMuonCuts->SetFillColor(0);
	    legendMuonCuts->SetBorderSize(0);
	    legendMuonCuts->SetTextFont(42);
	    legendMuonCuts->AddEntry(gMassResPtMC,"p_{T}^{#mu(1)} > 3.5, p_{T}^{#mu(2)} > 4","p");
	    legendMuonCuts->AddEntry(gNoLabel,"Width of additional CB shape (same cuts)","p");
	    legendMuonCuts->Draw();
	    TLatex *lCMSpt = new TLatex(8,0.225,"CMS 2.76 TeV Pythia Simulation (internal)");
	    lCMSpt->SetTextFont(4100);//42
	    lCMSpt->SetTextSize(0.03);
	    lCMSpt->Draw();
	    ppt1->Update(); 
	    TLatex *lMCPt= new TLatex(10,0.025,"#varUpsilon(1S), |y| < 2.4");
	    lMCPt->SetTextSize(0.029);
	    lMCPt->Draw();
	    //  ppt1->Update(); 
	    TLatex *sigmaFraction = new TLatex();
	    sigmaFraction->SetTextSize(0.03);    
	    sigmaFraction->DrawLatex(1.05,0.057,"68.3%");
	    sigmaFraction->DrawLatex(3.3,0.057,"68.2%");
	    sigmaFraction->DrawLatex(6,0.057,"70.8%");
	    sigmaFraction->DrawLatex(9,0.057,"70.2%");
	    sigmaFraction->DrawLatex(15,0.057,"73.2%");
	    
	    ppt1->Update(); 
	    cpt->SaveAs("MassResPtMC_3p5.pdf");
	    cpt->SaveAs("MassResPtMC_3p5.png");

	    cpt_raa->cd();
	    ppt1->Draw();
	    ppt1->cd();
	    f4Pt->Draw();
	    f4Pt->GetYaxis()->SetRangeUser(0.0,1.3);
	    f4Pt->GetYaxis()->SetTitle("almost R_{AA}");
	    TGraphErrors *gRaaPt = new TGraphErrors(nPtBins_2013,pt,raapt,pte,raapte);
	    gRaaPt->SetMarkerColor(kGreen+2);
	    gRaaPt->SetMarkerStyle(33);
	    gRaaPt->SetMarkerSize(2);
	    gRaaPt->Draw("pe");
	    TGraphErrors *gRaaPt3p5 = new TGraphErrors(nPtBins_2013,pt,raapt,pte,raapte);
	    gRaaPt3p5->SetMarkerColor(kGreen+2);
	    gRaaPt3p5->SetMarkerStyle(33);
	    gRaaPt3p5->SetMarkerSize(2);
	    gRaaPt3p5->Draw("pe");
	    cpt_raa->SaveAs("raatemporary_pt3p5.png");
	    // TGraphErrors *gYieldPtaa = new TGraphErrors(nPtBins_2013,pt,signal1S_aa_pt,pte,signal1S_aa_pte);
	    // gYieldPtaa->SetMarkerColor(kGreen+2);
	    // gYieldPtaa->SetMarkerStyle(33);
	    // gYieldPtaa->SetMarkerSize(2);
	    //npow and alpha now.
	    cCBPt1->cd();
	    ppt1->Draw();
	    ppt1->cd();
	    f4Pt->GetYaxis()->SetTitle("1S uncorr. Yield");
	    f4Pt->GetYaxis()->SetRangeUser(0.0,3.5);  
	    f4Pt->Draw();
	    TGraphErrors *gNpowPtMC = new TGraphErrors(nPtBins_2013,pt,npowMC_pt,pte,npowMC_pte);
	    gNpowPtMC->SetMarkerColor(kRed);
	    gNpowPtMC->SetMarkerStyle(33);
	    gNpowPtMC->SetMarkerSize(2);
	    TGraphErrors *gAlphaPtMC = new TGraphErrors(nPtBins_2013,pt,alphaMC_pt,pte,alphaMC_pte);
	    gAlphaPtMC->SetMarkerStyle(33);
	    gAlphaPtMC->SetMarkerSize(2);
	    gAlphaPtMC->SetMarkerColor(kBlue);
	    gNpowPtMC->Draw("pe");
	    gAlphaPtMC->Draw("p");
	    TLegend *legendCBpars = new TLegend(0.1,0.75,0.5,0.9);
	    legendCBpars->SetTextSize(0.032);
	    legendCBpars->SetFillStyle(0);
	    legendCBpars->SetFillColor(0);
	    legendCBpars->SetBorderSize(0);
	    legendCBpars->SetTextFont(42);
	    legendCBpars->AddEntry(gNpowPtMC,"n_{CB}","p");
	    legendCBpars->AddEntry(gAlphaPtMC,"#alpha_{CB}","p");
	    if(do3p5){ legendCBpars->AddEntry(gNoLabel,"single #mu asymetric p_{T} cuts","");}
	    
	    if(do4){ legendCBpars->AddEntry(gNoLabel," #mu p_{T} > 4 GeV/c ","");}
	    legendCBpars->Draw();
	    // lCMSpt->Draw();
	    TLatex *lMCPt= new TLatex();
    lMCPt->SetTextSize(0.029);
	    lMCPt->DrawLatex(10,0.5,"Crystal Ball parameters");
	    lMCPt->DrawLatex(10,0.65,"#varUpsilon(1S), |y| < 2.4");
	
	    lMCPt->Draw();
	    ppt1->Update(); 
	    cCBPt1->SaveAs("CBPtMC_3p5.pdf");
	    cCBPt1->SaveAs("CBPtMC_3p5.png");

	  }
	  if(doRap){
	    crap->cd();
	    prap1->Draw();
	    prap1->cd();
	    f4Rap->Draw();
	    TGraphErrors *gMassResRapMC = new TGraphErrors(nRapBins_2013,rap,massRes1_MC_rap,rape,massRes1_MC_rape);
	    gMassResRapMC->SetMarkerColor(kGreen+2);
	    gMassResRapMC->SetMarkerStyle(33);
	    gMassResRapMC->SetMarkerSize(2);
	    TGraphErrors *gMassResRapMCcircle = new TGraphErrors(nRapBins_2013,rap,massRes1_MC_rap,rape,massRes1_MC_rape);
	    gMassResRapMCcircle->SetMarkerStyle(27);
	    gMassResRapMCcircle->SetMarkerSize(2);
	    gMassResRapMCcircle->SetLineColor(kBlack);
	    TGraphErrors *gMassResRapMC_2 = new TGraphErrors(nRapBins_2013,rap,massRes2_MC_rap,rape,massRes2_MC_rape);
	    gMassResRapMC_2->SetMarkerStyle(27);
	    gMassResRapMC_2->SetMarkerSize(2);
	    gMassResRapMC_2->SetLineColor(kBlack);

	    gMassResRapMCcircle->Draw("pe");
	    gMassResRapMC->Draw("p");
	    gMassResRapMC_2->Draw("p");
	    f4Rap->Draw("same");
	    TGraphErrors *gNoLabel = new TGraphErrors(bin,centnoErr,centnoErr,centnoErr,centnoErr);
	    gNoLabel->SetMarkerStyle(27);
	    gNoLabel->SetMarkerSize(2);
	    gNoLabel->SetLineColor(kWhite);
	    gNoLabel->SetFillColor(kWhite);
	  
	    TLegend *legendMuonCuts = new TLegend(0.1,0.8,0.5,0.9);
	    legendMuonCuts->SetTextSize(0.029);
	    legendMuonCuts->SetFillStyle(0);
	    legendMuonCuts->SetFillColor(0);
	    legendMuonCuts->SetBorderSize(0);
	    legendMuonCuts->SetTextFont(42);
	    legendMuonCuts->AddEntry(gMassResRapMC,"p_{T}^{#mu(1)} > 3.5, p_{T}^{#mu(2)} > 4","p");
	    legendMuonCuts->AddEntry(gNoLabel,"Width of additional CB shape (same cuts)","p");
	    legendMuonCuts->Draw();
	    TLatex *lCMSrap = new TLatex(1.,0.225,"CMS 2.76 TeV Pythia Simulation (internal)");
	    lCMSrap->SetTextFont(4100);//42
	    lCMSrap->SetTextSize(0.03);
	    lCMSrap->Draw();
	    prap1->Update(); 
	    TLatex *lMCRap= new TLatex(1,0.025,"#varUpsilon(1S), |y| < 2.4");
	    lMCRap->SetTextSize(0.029);
	    lMCRap->Draw();
	    prap1->Update(); 
	    TLatex *sigmaFraction = new TLatex();
	    sigmaFraction->SetTextSize(0.03);    
	    sigmaFraction->DrawLatex(0.12,0.04,"89%");
	    sigmaFraction->DrawLatex(0.43,0.04,"87.8%");
	    sigmaFraction->DrawLatex(0.75,0.06,"86%");
	    sigmaFraction->DrawLatex(1.15,0.079,"84.8%");
	    sigmaFraction->DrawLatex(1.85,0.095,"71.8%");
	  
	    prap1->Update(); 
	    crap->SaveAs("MassResRapMC_3p5.pdf");
	    crap->SaveAs("MassResRapMC_3p5.png");
	    
	       crap_raa->cd();
	    prap1->Draw();
	    prap1->cd();
	    f4Rap->Draw();
	    f4Rap->GetYaxis()->SetRangeUser(0.0,1.3);
	    
	    TGraphErrors *gRaaRap = new TGraphErrors(nRapBins_2013,rap,raarap3p5,rape,raarape3p5);
	    gRaaRap->SetMarkerColor(kGreen+2);
	    gRaaRap->SetMarkerStyle(33);
	    gRaaRap->SetMarkerSize(2);
	    gRaaRap->Draw("pe");
	    crap_raa->SaveAs("raatemporary_rap3p5.png");
	    	    //npow and alpha now.
	    cCBRap1->cd();
	    prap1->Draw();
	    prap1->cd();
	    f4Rap->GetYaxis()->SetTitle("Crystal Ball tail parameters (arb. units)");
	    f4Rap->GetYaxis()->SetRangeUser(0.0,3.5);  
	    f4Rap->Draw();
	    TGraphErrors *gNpowRapMC = new TGraphErrors(nRapBins_2013,rap,npowMC_rap,rape,npowMC_rape);
	    gNpowRapMC->SetMarkerColor(kRed);
	    gNpowRapMC->SetMarkerStyle(33);
	    gNpowRapMC->SetMarkerSize(2);
	    TGraphErrors *gAlphaRapMC = new TGraphErrors(nRapBins_2013,rap,alphaMC_rap,rape,alphaMC_rape);
	    gAlphaRapMC->SetMarkerStyle(33);
	    gAlphaRapMC->SetMarkerSize(2);
	    gAlphaRapMC->SetMarkerColor(kBlue);
	    gNpowRapMC->Draw("pe");
	    gAlphaRapMC->Draw("p");
	    TLegend *legendCBpars = new TLegend(0.1,0.75,0.5,0.9);
	    legendCBpars->SetTextSize(0.032);
	    legendCBpars->SetFillStyle(0);
	    legendCBpars->SetFillColor(0);
	    legendCBpars->SetBorderSize(0);
	    legendCBpars->SetTextFont(42);
	    legendCBpars->AddEntry(gNpowRapMC,"n_{CB}","p");
	    legendCBpars->AddEntry(gAlphaRapMC,"#alpha_{CB}","p");
	    if(do3p5){ legendCBpars->AddEntry(gNoLabel,"single #mu asymetric p_{T} cuts","");}
	    
	    if(do4){ legendCBpars->AddEntry(gNoLabel," #mu p_{T} > 4 GeV/c ","");}
	    legendCBpars->Draw();
	    // lCMSrap->Draw();
	    TLatex *lMCRap= new TLatex();
	    lMCRap->SetTextSize(0.029);
	    lMCRap->DrawLatex(10,0.5,"Crystal Ball parameters");
	    lMCRap->DrawLatex(10,0.65,"#varUpsilon(1S), |y| < 2.4");
	
	    lMCRap->Draw();
	    prap1->Update(); 
	    cCBRap1->SaveAs("CBRapMC_3p5.pdf");
	    cCBRap1->SaveAs("CBRapMC_3p5.png");
	  }
	  
	do4=true;
	}


      //	raapte3p5_scaled[l]=  
      //	raarapt3p5_scaled[l]


      if(do4)
      	{    
	  if(doPt)
	    {
	      cpt2->cd();
	      ppt1->Draw();
	      ppt1->cd();
	      f4Pt->GetYaxis()->SetTitle("Mass peak resolution (GeV/c)");
	      f4Pt->GetYaxis()->SetRangeUser(0.0,0.22);
	      f4Pt->Draw();
	      do3p5=false;
	      gROOT->Macro("DataParameters.h");
	      TGraphErrors *gMassResPtMC = new TGraphErrors(nPtBins_2013,pt,massRes1_MC_pt,pte,massRes1_MC_pte);
	      gMassResPtMC->SetMarkerColor(kBlue+2);
	      gMassResPtMC->SetMarkerStyle(33);
	      gMassResPtMC->SetMarkerSize(2);
	      TGraphErrors *gMassResPtMCcircle = new TGraphErrors(nPtBins_2013,pt,massRes1_MC_pt,pte,massRes1_MC_pte);
	      gMassResPtMCcircle->SetMarkerColor(kBlack);
	      gMassResPtMCcircle->SetMarkerStyle(27);
	      gMassResPtMCcircle->SetMarkerSize(2);
	      TGraphErrors *gMassResPtMC_2 = new TGraphErrors(nPtBins_2013,pt,massRes2_MC_pt,pte,massRes2_MC_pte);
	      gMassResPtMC_2->SetMarkerStyle(27);
	      gMassResPtMC_2->SetMarkerSize(2);
	      gMassResPtMC_2->SetLineColor(kBlack);

	      gMassResPtMCcircle->Draw("pe");
	      gMassResPtMC->Draw("p");
	      gMassResPtMC_2->Draw("p");
	      f4Pt->Draw("same");
	      TGraphErrors *gNoLabel = new TGraphErrors(bin,centnoErr,centnoErr,centnoErr,centnoErr);
	      gNoLabel->SetMarkerStyle(27);
	      gNoLabel->SetMarkerSize(2);
	      gNoLabel->SetLineColor(kWhite);
	      gNoLabel->SetFillColor(kWhite);
	  
	      TLegend *legendMuonCuts = new TLegend(0.1,0.8,0.5,0.9);
	      legendMuonCuts->SetTextSize(0.029);
	      legendMuonCuts->SetFillStyle(0);
	      legendMuonCuts->SetFillColor(0);
	      legendMuonCuts->SetBorderSize(0);
	      legendMuonCuts->SetTextFont(42);
	      legendMuonCuts->AddEntry(gMassResPtMC,"p_{T}(#mu) > 4","p");
	      legendMuonCuts->AddEntry(gNoLabel,"Width of additional CB shape (same cuts)");
	      legendMuonCuts->Draw();
	      lCMSpt->Draw();
	      ppt1->Update(); 
	      lMCPt->Draw();
	      ppt1->Update(); 
	      TLatex *sigmaFraction_pt4 = new TLatex();
	      sigmaFraction_pt4->SetTextSize(0.03);    
	      sigmaFraction_pt4->DrawLatex(1.05,0.057,"66.3%");
	      sigmaFraction_pt4->DrawLatex(3.3,0.057,"67.9%");
	      sigmaFraction_pt4->DrawLatex(6,0.057,"72.3%");
	      sigmaFraction_pt4->DrawLatex(9,0.057,"69.5%");
	      sigmaFraction_pt4->DrawLatex(15,0.057,"71.8%");
   
	      ppt1->Update(); 
	      cpt2->SaveAs("MassResPtMC_4.pdf");
	      cpt2->SaveAs("MassResPtMC_4.png");

	      cpt_raa2->cd();
	      ppt1->Draw();
	      ppt1->cd();
	      f4Pt->Draw();
	      f4Pt->GetYaxis()->SetRangeUser(0.0,1.);
	      f4Pt->GetYaxis()->SetTitle("almost R_{AA}");
	      TGraphErrors *gRaaPt = new TGraphErrors(nPtBins_2013,pt,raapt4,pte,raapte4);
	      gRaaPt->SetMarkerColor(kBlue+2);
	      gRaaPt->SetMarkerStyle(33);
	      gRaaPt->SetMarkerSize(2);
	      gRaaPt->Draw("pe");
	    
	      // TGraphErrors *gRaaPtScaled = new TGraphErrors(nPtBins_2013,PtForScaling,raapt4,pte,raapte3p5_scaled);
	      // gRaaPtScaled->SetMarkerColor(kGreen+2);
	      // gRaaPtScaled->SetMarkerStyle(33);
	      // gRaaPtScaled->SetMarkerSize(2);
	      // gRaaPtScaled->Draw("pe");
	      TGraphErrors *gRaaPtNonScaled = new TGraphErrors(nPtBins_2013,PtForScaling,raapt3p5,pte,raapte3p5);
	      gRaaPtNonScaled->SetMarkerColor(kGreen+2);
	      gRaaPtNonScaled->SetMarkerStyle(33);
	      gRaaPtNonScaled->SetMarkerSize(2);
	      gRaaPtNonScaled->Draw("pe");
	      TLegend *legendRaa = new TLegend(0.1,0.8,0.5,0.9);
	      legendRaa->SetTextSize(0.029);
	      legendRaa->SetFillStyle(0);
	      legendRaa->SetFillColor(0);
	      legendRaa->SetBorderSize(0);
	      legendRaa->SetTextFont(42);
	      //	      legendRaa->AddEntry(gRaaPtScaled,"p_{T}(#mu) > 3.5 OR 4 (scaled to 4GeV measurement)","p");
	      legendRaa->AddEntry(gRaaPtNonScaled,"p_{T}(#mu) > 3.5 OR 4","p");
	      legendRaa->AddEntry(gRaaPt,"p_{T}(#mu) > 4","p");
	      legendRaa->Draw();
	      cpt_raa2->SaveAs("raatemporary_pt4_withoutscaling.png");


	      //npow and alpha now.
	      cCBPt2->cd();
	      ppt1->Draw();
	      ppt1->cd();
	      f4Pt->GetYaxis()->SetTitle("Crystal Ball tail parameters (arb. units)");
	      f4Pt->GetYaxis()->SetRangeUser(0.0,3.5);  
	      f4Pt->Draw();
	      TGraphErrors *gNpowPtMC = new TGraphErrors(nPtBins_2013,pt,npowMC_pt,pte,npowMC_pte);
	      gNpowPtMC->SetMarkerColor(kRed);
	      gNpowPtMC->SetMarkerStyle(33);
	      gNpowPtMC->SetMarkerSize(2);
	      TGraphErrors *gAlphaPtMC = new TGraphErrors(nPtBins_2013,pt,alphaMC_pt,pte,alphaMC_pte);
	      gAlphaPtMC->SetMarkerStyle(33);
	      gAlphaPtMC->SetMarkerSize(2);
	      gAlphaPtMC->SetMarkerColor(kBlue);
	      gNpowPtMC->Draw("pe");
	      gAlphaPtMC->Draw("p");
	      TLegend *legendCBpars = new TLegend(0.1,0.75,0.5,0.9);
	      legendCBpars->SetTextSize(0.032);
	      legendCBpars->SetFillStyle(0);
	      legendCBpars->SetFillColor(0);
	      legendCBpars->SetBorderSize(0);
	      legendCBpars->SetTextFont(42);
	      legendCBpars->AddEntry(gNpowPtMC,"n_{CB}","p");
	      legendCBpars->AddEntry(gAlphaPtMC,"#alpha_{CB}","p");
	      if(do3p5){ legendCBpars->AddEntry(gNoLabel,"single #mu asymetric p_{T} cuts","");}
	    
	      if(do4){ legendCBpars->AddEntry(gNoLabel," #mu p_{T} > 4 GeV/c ","");}
	      legendCBpars->Draw();
	      TLatex *lMCPt= new TLatex();
	      lMCPt->SetTextSize(0.029);
	      lMCPt->DrawLatex(10,0.5,"Crystal Ball parameters");
	      lMCPt->DrawLatex(10,0.65,"#varUpsilon(1S), |y| < 2.4");
	
	      lMCPt->Draw();
	      ppt1->Update(); 
	      cCBPt2->SaveAs("CBPtMC_4.pdf");
	      cCBPt2->SaveAs("CBPtMC_4.png");

	    }
	  if(doRap){
	    crap2->cd();
	    prap1->Draw();
	    prap1->cd();
	 
	    f4Rap->GetYaxis()->SetRangeUser(0.0,0.22);
	    f4Rap->Draw();
	    TGraphErrors *gMassResRapMC = new TGraphErrors(nRapBins_2013,rap,massRes1_MC_rap,rape,massRes1_MC_rape);
	    gMassResRapMC->SetMarkerColor(kBlue+2);
	    gMassResRapMC->SetMarkerStyle(33);
	    gMassResRapMC->SetMarkerSize(2);
	    TGraphErrors *gMassResRapMCcircle = new TGraphErrors(nRapBins_2013,rap,massRes1_MC_rap,rape,massRes1_MC_rape);
	    gMassResRapMCcircle->SetMarkerStyle(27);
	    gMassResRapMCcircle->SetMarkerSize(2);
	    gMassResRapMCcircle->SetLineColor(kBlack);
	    TGraphErrors *gMassResRapMC_2 = new TGraphErrors(nRapBins_2013,rap,massRes2_MC_rap,rape,massRes2_MC_rape);
	    gMassResRapMC_2->SetMarkerStyle(27);
	    gMassResRapMC_2->SetMarkerSize(2);
	    gMassResRapMC_2->SetLineColor(kBlack);

	    gMassResRapMCcircle->Draw("pe");
	    gMassResRapMC->Draw("p");
	    gMassResRapMC_2->Draw("p");
	    f4Rap->Draw("same");
	    TGraphErrors *gNoLabel = new TGraphErrors(bin,centnoErr,centnoErr,centnoErr,centnoErr);
	    gNoLabel->SetMarkerStyle(27);
	    gNoLabel->SetMarkerSize(2);
	    gNoLabel->SetLineColor(kWhite);
	    gNoLabel->SetFillColor(kWhite);
	  
	    TLegend *legendMuonCuts = new TLegend(0.1,0.8,0.5,0.9);
	    legendMuonCuts->SetTextSize(0.029);
	    legendMuonCuts->SetFillStyle(0);
	    legendMuonCuts->SetFillColor(0);
	    legendMuonCuts->SetBorderSize(0);
	    legendMuonCuts->SetTextFont(42);
	    legendMuonCuts->AddEntry(gMassResPtMC,"p_{T}(#mu) > 4","p");
	    legendMuonCuts->AddEntry(gNoLabel,"Width of additional CB shape (same cuts)","p");
	    legendMuonCuts->Draw();
	    TLatex *lCMSrap = new TLatex(1.,0.225,"CMS 2.76 TeV Pythia Simulation (internal)");
	    lCMSrap->SetTextFont(4100);//42
	    lCMSrap->SetTextSize(0.03);
	    lCMSrap->Draw();
	    prap1->Update(); 
	    TLatex *lMCRap= new TLatex(1,0.025,"#varUpsilon(1S), |y| < 2.4");
	    lMCRap->SetTextSize(0.029);
	    lMCRap->Draw();
	    prap1->Update(); 
	    TLatex *sigmaFraction = new TLatex();
	    sigmaFraction->SetTextSize(0.03);    
	    sigmaFraction->DrawLatex(0.12,0.04,"92.3%");
	    sigmaFraction->DrawLatex(0.43,0.04,"90.2%");
	    sigmaFraction->DrawLatex(0.75,0.06,"83.1%");
	    sigmaFraction->DrawLatex(1.15,0.079,"87.6%");
	    sigmaFraction->DrawLatex(1.85,0.095,"72.7%");
	  
	    prap1->Update(); 
	    crap2->SaveAs("MassResRapMC_4.pdf");
	    crap2->SaveAs("MassResRapMC_4.png");
	    
	    crap_raa2->cd();
	    prap1->Draw();
	    prap1->cd();
	    f4Rap->Draw();
	    f4Rap->GetYaxis()->SetRangeUser(0.0,1.);
	    f4Rap->GetYaxis()->SetTitle("almost R_{AA}");
	    TGraphErrors *gRaaRap = new TGraphErrors(nRapBins_2013,rap,raarap4,rape,raarape4);
	    gRaaRap->SetMarkerColor(kBlue+2);
	    gRaaRap->SetMarkerStyle(33);
	    gRaaRap->SetMarkerSize(2);
	    gRaaRap->Draw("pe");
	    
	    // TGraphErrors *gRaaRapScaled = new TGraphErrors(nRapBins_2013,RapForScaling,raarap4,rape,raarape3p5_scaled);
	    // gRaaRapScaled->SetMarkerColor(kGreen+2);
	    // gRaaRapScaled->SetMarkerStyle(33);
	    // gRaaRapScaled->SetMarkerSize(2);
	    // gRaaRapScaled->Draw("pe");
	    TGraphErrors *gRaaRapNonScaled = new TGraphErrors(nRapBins_2013,RapForScaling,raarap3p5,rape,raarape3p5);
	    gRaaRapNonScaled->SetMarkerColor(kGreen+2);
	    gRaaRapNonScaled->SetMarkerStyle(33);
	    gRaaRapNonScaled->SetMarkerSize(2);
	    gRaaRapNonScaled->Draw("pe");
	    TLegend *legendRaa = new TLegend(0.1,0.8,0.5,0.9);
	    legendRaa->SetTextSize(0.029);
	    legendRaa->SetFillStyle(0);
	    legendRaa->SetFillColor(0);
	    legendRaa->SetBorderSize(0);
	    legendRaa->SetTextFont(42);
	    //	    legendRaa->AddEntry(gRaaRapScaled,"p_{T}(#mu) > 3.5 OR 4 (scaled to 4GeV measurement)","p");
	    legendRaa->AddEntry(gRaaRapNonScaled,"p_{T}(#mu) > 3.5 OR 4","p");
	    legendRaa->AddEntry(gRaaRap,"p_{T}(#mu) > 4","p");
	    legendRaa->Draw();
	    crap_raa2->SaveAs("raatemporary_rap4_withoutscaling.png");

	    //npow and alpha now.
	    cCBRap2->cd();
	    prap1->Draw();
	    prap1->cd();
	    f4Rap->GetYaxis()->SetTitle("Crystal Ball tail parameters (arb. units)");
	    f4Rap->GetYaxis()->SetRangeUser(0.0,3.5);  
	    f4Rap->Draw();
	    TGraphErrors *gNpowRapMC = new TGraphErrors(nRapBins_2013,rap,npowMC_rap,rape,npowMC_rape);
	    gNpowRapMC->SetMarkerColor(kRed);
	    gNpowRapMC->SetMarkerStyle(33);
	    gNpowRapMC->SetMarkerSize(2);
	    TGraphErrors *gAlphaRapMC = new TGraphErrors(nRapBins_2013,rap,alphaMC_rap,rape,alphaMC_rape);
	    gAlphaRapMC->SetMarkerStyle(33);
	    gAlphaRapMC->SetMarkerSize(2);
	    gAlphaRapMC->SetMarkerColor(kBlue);
	    gNpowRapMC->Draw("pe");
	    gAlphaRapMC->Draw("p");
	    TLegend *legendCBpars = new TLegend(0.1,0.75,0.5,0.9);
	    legendCBpars->SetTextSize(0.032);
	    legendCBpars->SetFillStyle(0);
	    legendCBpars->SetFillColor(0);
	    legendCBpars->SetBorderSize(0);
	    legendCBpars->SetTextFont(42);
	    legendCBpars->AddEntry(gNpowRapMC,"n_{CB}","p");
	    legendCBpars->AddEntry(gAlphaRapMC,"#alpha_{CB}","p");
	    if(do3p5){ legendCBpars->AddEntry(gNoLabel,"single #mu asymetric p_{T} cuts","");}
	    
	    if(do4){ legendCBpars->AddEntry(gNoLabel," #mu p_{T} > 4 GeV/c ","");}
	    legendCBpars->Draw();
	    // lCMSrap->Draw();
	    TLatex *lMCRap= new TLatex();
	    lMCRap->SetTextSize(0.029);
	    lMCRap->DrawLatex(10,0.5,"Crystal Ball parameters");
	    lMCRap->DrawLatex(10,0.65,"#varUpsilon(1S), |y| < 2.4");
	    
	    lMCRap->Draw();
	    prap1->Update(); 
	    cCBRap2->SaveAs("CBRapMC_4.pdf");
	    cCBRap2->SaveAs("CBRapMC_4.png");   
	  }
	}
      




      TCanvas *cNpart = new TCanvas("cNpart","cNpart");
      float centReverse2013[8]={381.2,330.3,261.4,187.3,130.1,86.3,42.02,8.75};
      float centErr[8]={6,6,6,6,6,6,6};
      float centnoErr[8]={0,0,0,0,0,0,0,0};
      float scaled4GeV[8];
  for(int l =0;l<nCentBins_2013;l++)
    {
      do3p5=true;
      if(do3p5){
	gROOT->Macro("DataParameters.h");
	signal1S_aa_npart[l]=(signal1S_aa_npart[l]*fMB[l]);
	signal1S_aa_nparte[l]=(signal1S_aa_nparte[l]*fMB[l]);

	cout << signal1S_aa_npart[l] << " " << signal1S_aa_nparte[l] << endl;
	do4=true;
      }
      if(do4){
	do3p5=false;
	gROOT->Macro("DataParameters.h");

	signal1S_aa_scaled_npart[l]=signal1S_aa_npart[l]/(taa2013[l]/fMB[l]);
	signal1S_aa_scaled_nparte[l]=signal1S_aa_nparte[l]/(taa2013[l]/fMB[l]);
	
	cout << signal1S_aa_scaled_npart[l] << " " << signal1S_aa_scaled_nparte[l] << endl;
      }
    }
      TPad *pnpart1 = new TPad("pnpart1","pnpart1",0.0,0.0,1.0,1.0);
      pnpart1->SetBottomMargin(0.12);
      pnpart1->SetTopMargin(0.08);
      pnpart1->SetRightMargin(0.06);
      pnpart1->SetLeftMargin(0.11);

      TF1 *f4Npart = new TF1("f4Npart","0",0,416);
      f4Npart->SetLineWidth(0);
      f4Npart->GetXaxis()->SetTitle("N_{Part}");
      f4Npart->GetYaxis()->SetTitle("Yield / T_{AA}");
      f4Npart->GetYaxis()->SetTitleOffset(1.6);
      f4Npart->GetYaxis()->SetLabelSize(0.032);
      f4Npart->GetYaxis()->SetTitleSize(0.03);
      f4Npart->GetXaxis()->CenterTitle(kTRUE);
      cNpart->cd();
      pnpart1->Draw();
      pnpart1->cd();
        f4Npart->GetYaxis()->SetRangeUser(0.0,400);
      f4Npart->Draw();

      TGraphErrors *gRaaCent = new TGraphErrors(nCentBins_2013,centReverse2013,signal1S_aa_npart,centnoErr,signal1S_aa_nparte);
      gRaaCent->SetMarkerColor(kBlue+2);
      gRaaCent->SetMarkerStyle(33);
      gRaaCent->SetMarkerSize(2);
      gRaaCent->Draw("pe");
      

      // TGraphErrors *gRaaCentScaled = new TGraphErrors(nPtBins_2013,centReverse2013ForScaling,raapt4,centnoErr,raapte3p5_scaled);
      // gRaaCentScaled->SetMarkerColor(kGreen+2);
      // gRaaCentScaled->SetMarkerStyle(33);
      // gRaaCentScaled->SetMarkerSize(2);
      // gRaaCentScaled->Draw("pe");
      TLegend *legendYaa = new TLegend(0.1,0.8,0.5,0.9);
      legendYaa->SetTextSize(0.029);
      legendYaa->SetFillStyle(0);
      legendYaa->SetFillColor(0);
      legendYaa->SetBorderSize(0);
      legendYaa->SetTextFont(42);
      legendYaa->AddEntry(gRaaCent,"p_{T}(#mu) > 3.5 OR 4 (scaled to 4GeV measurement)","p");
      //    legendYaa->AddEntry(gRaaCent,"p_{T}(#mu) > 4","p");
      legendYaa->Draw();
      cNpart->SaveAs("YieldTaaNpartscaling.png");
      
}
