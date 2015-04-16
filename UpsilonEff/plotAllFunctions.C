#include "headerForTNP.h"
 double myFunc(double *x, double *par){
   Float_t xx =x[0];
   double f; 
   // par[0] is the normalisation;
   // par[1] is the mu_data;
   // par[2] is the sigma_data;
   // par[3] is the index in the array : 0 is the nominal case, 1-101 are the 100 variations.
   // par[4] is the Abs(eta) region : 0 is barrel, 1 is endcap muons.
   // par[5] is the collision type : 0 is pp, 1 is pbpb.

   ///let's go:
   if(par[5]>0){ //doing pbpb first,
     if(par[4]<1) { // barrel,
       if(par[3]<1) { // nominal 
	 f = (par[0]*TMath::Erf((xx-par[1])/par[2]))/(TMath::Erf((xx-1.7960)/2.5160));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2])/(0.9576*TMath::Erf((xx-1.7883)/2.6583)));  /// 
       }
     }
     if(par[4]>0){// endcap,
       if(par[3]<1) { // nominal case
     	 f = (par[0]*TMath::Erf((xx-par[1])/par[2]))/(0.7810*TMath::Erf((xx-1.3609)/2.1231));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2])/(0.7810*TMath::Erf((xx-1.3609)/2.1231)));  /// 
       }
     }
   }
   if(par[5]<1){ //doing pp now,
     if(par[4]<1) { // barrel,
       if(par[3]<1) { // nominal 
	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2]))/(TMath::Erf((xx-1.7960)/2.5160));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2])/ (0.9604*TMath::Erf((xx-2.0583)/2.1573)));  /// 
       }
     }
     if(par[4]>0){// endcap,
       if(par[3]<1) { // nominal case
     	 f = (par[0]*TMath::Erf((xx-par[1])/par[2]))/(0.7364*TMath::Erf((xx-1.2538)/2.2530));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2])/ (0.7364*TMath::Erf((xx-1.2149)/2.3352)));  /// 
       }
     }
   }
   return f;

   //pbpb midrap
   //
   //
   
   //pbpb fwdrap
   //j>0  (0.7810*TMath::Erf((x-1.3609)/2.1231))
   //j=0  (0.8299*TMath::Erf((pt-1.2785)/1.8833))/(0.7810*TMath::Erf((pt-1.3609)/2.1231)  
   
   //pp midrap
   //j>0  (0.9604*TMath::Erf((x-2.0583)/2.1573))
   //j=0  (0.998474*TMath::Erf((pt-1.63555)/(0.797933))/TMath::Erf((pt-0.222866)/(2.95593)))

   //pp fwd
   //j>0 (0.7364*TMath::Erf((x-1.2149)/2.3352))
   //j=0 (0.7788*TMath::Erf((pt-1.1903)/1.9880))/(0.7364*TMath::Erf((pt-1.2538)/2.2530))
 }

double myNumerator( double *x, double *par){
  double f;
  double xx = x[0];
    ///let's go:
   if(par[5]>0){ //doing pbpb first,
     if(par[4]<1) { // barrel,
       if(par[3]<1) { // nominal 
	 f = (par[0]*TMath::Erf((xx-par[1])/par[2]));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2]));  /// 
       }
     }
     if(par[4]>0){// endcap,
       if(par[3]<1) { // nominal case
     	 f = (par[0]*TMath::Erf((xx-par[1])/par[2]));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2]));  /// 
       }
     }
   }
   if(par[5]<1){ //doing pp now,
     if(par[4]<1) { // barrel,
       if(par[3]<1) { // nominal 
	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2]));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2]));  /// 
       }
     }
     if(par[4]>0){// endcap,
       if(par[3]<1) { // nominal case
     	 f = (par[0]*TMath::Erf((xx-par[1])/par[2]));
       }
       if(par[3]>0){ // variations
       	 f =   (par[0]*TMath::Erf((xx-par[1])/par[2]));  /// 
       }
     }
   }
  return f;
}
void plotAllFunctions()
{
  double pbpb=1;  // if pbpb=0, you're doing pp.
  double endcap=0;//if endcap=0, you're doing barrel ;)
  int j=0;     // cool increment.
  bool doOnlyNumerator= false;
  // /// ...the classic ones... // ///
  TCanvas *cBarrel_pbpb = new TCanvas ("cBarrel_pbpb","cBarrel_pbpb",1400,700);
  cBarrel_pbpb->Divide(2);
  cBarrel_pbpb->cd(1);
  TF1 *f1 = new TF1("myFunc",myFunc,0,20,6);
  if (doOnlyNumerator){ f1= new TF1("myFunk",myNumerator,0,20,6);}
  cout << j << endl;
  cout << endcap << endl;
  f1->SetParameters(alpha_pbpb_midrap[j],mu_data_pbpb_midrap[j],sigma_data_pbpb_midrap[j],j,endcap,pbpb);
  f1->SetRange(0,20);
  if(!doOnlyNumerator){ f1->GetYaxis()->SetRangeUser(0.8,1.5); }
  else if(doOnlyNumerator){   f1->GetYaxis()->SetRangeUser(0.0,1.5); }
  f1->Draw();
  std::cout << alpha_pbpb_midrap[j] <<"*TMath::Erf((x-"<< mu_data_pbpb_midrap[j] <<"/"<<sigma_data_pbpb_midrap[j]<<")) / (1.0000*TMath::Erf((x-1.7960)/2.5160))"<<std::endl;
  
  for(j=1;j< 101;j++){
    //reminder:
    // p[0] is the normalisation;
    // p[1] is the mu_data;
    // p[2] is the sigma_data;
    // p[3] is the index in the array : 0 is the nominal case, 1-101 are the 100 variations.
    // p[4] is the Abs(eta) region : 0 is barrel, 1 is endcap muons.
    // p[5] is the collision type : 0 is pp, 1 is pbpb.
    fb = new TF1("myfunc",myFunc,0,20,6);
    if (doOnlyNumerator){ fb= new TF1("myFunk",myNumerator,0,20,6);}
    //  fb->SetRange(0,20);
    std::cout << j << alpha_pbpb_midrap[j] <<"*TMath::Erf((x-"<< mu_data_pbpb_midrap[j] <<"/"<<sigma_data_pbpb_midrap[j]<<")) / (0.9576*TMath::Erf((x-1.7883)/2.6583))"<<std::endl;
    fb->SetParameters(alpha_pbpb_midrap[j],mu_data_pbpb_midrap[j],sigma_data_pbpb_midrap[j],j,endcap,pbpb);
    fb->SetLineColor(kAzure+1);
    fb->Draw("same");
   }
  f1->Draw("same");
  TLegend *legend = new TLegend(0.2,0.75,0.9,0.9);
  legend->SetTextSize(0.038);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  if(!doOnlyNumerator){ 
  legend->AddEntry(f1,"PbPb ","");
  legend->AddEntry(f1,"Nominal T&P scaling: |#eta^{#mu}| < 1.6 ","lp");
  legend->AddEntry(fb,"100 T&P stat. variations: |#eta^{#mu}| < 1.6 ","lp");}
  else{
  legend->AddEntry(f1,"PbPb ","");
  legend->AddEntry(f1," #varepsilon^{DATA}_{#mu}, |#eta^{#mu}| < 1.6 ","lp");
  legend->AddEntry(fb,"100 T&P stat. variations, |#eta^{#mu}| < 1.6 ","lp");
  }
  legend->Draw();

  //endcap pbpb:  (0.8299*TMath::Erf((pt-1.2785)/1.8833))/(0.7810*TMath::Erf((pt-1.3609)/2.1231))
  cBarrel_pbpb->cd(2);
  endcap=1;
  cout << "MOFO DOING ENDCAP NOW" << endl;
  TF1 *f2 = new TF1("myFunc",myFunc,0,20,6);
  if (doOnlyNumerator){ f2= new TF1("myFunk",myNumerator,0,20,6);}
  j=0; // reset, mofo!
  f2->SetParameters(alpha_pbpb_fwdrap[j],mu_data_pbpb_fwdrap[j],sigma_data_pbpb_fwdrap[j],j,endcap,pbpb);
  f2->SetRange(0,20);
  if(!doOnlyNumerator){   f2->GetYaxis()->SetRangeUser(0.8,1.5); }
  else if(doOnlyNumerator){   f2->GetYaxis()->SetRangeUser(0.0,1.5); }
  f2->Draw();
  // std::cout << alpha_pbpb_fwdrap[j] <<"*TMath::Erf((x-"<< mu_data_pbpb_fwdrap[j] <<"/"<<sigma_data_pbpb_fwdrap[j]<<")) / (1.0000*TMath::Erf((x-1.7960)/2.5160))"<<std::endl;
  
  for(j=1;j< 101;j++){
    //reminder:
    // p[0] is the normalisation;
    // p[1] is the mu_data;
    // p[2] is the sigma_data;
    // p[3] is the index in the array : 0 is the nominal case, 1-101 are the 100 variations.
    // p[4] is the Abs(eta) region : 0 is barrel, 1 is endcap muons.
    // p[5] is the collision type : 0 is pp, 1 is pbpb.
    fe = new TF1("myfunc",myFunc,0,20,6);
    if (doOnlyNumerator){ fe= new TF1("myFunk",myNumerator,0,20,6);}
    fe->SetRange(0,20);
    //   std::cout  << alpha_pbpb_fwdrap[j] <<"*TMath::Erf((x-"<< mu_data_pbpb_fwdrap[j] <<"/"<<sigma_data_pbpb_fwdrap[j]<<")) / (0.9576*TMath::Erf((x-1.7883)/2.6583))"<<std::endl;
    // cout << f << endl;
    fe->SetParameters(alpha_pbpb_fwdrap[j],mu_data_pbpb_fwdrap[j],sigma_data_pbpb_fwdrap[j],j,endcap,pbpb);
    fe->SetLineColor(kAzure+1);
    fe->Draw("same");
  }
  f2->Draw("same");
  TLegend *leg2 = new TLegend(0.2,0.75,0.9,0.9);
  leg2->SetTextSize(0.038);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  if(!doOnlyNumerator){ 
    leg2->AddEntry(f2,"PbPb ","");
    leg2->AddEntry(f2,"Nominal T&P scaling: 1.6 < |#eta^{#mu}| < 2.4 ","lp");
    leg2->AddEntry(fe,"100 T&P stat. variations:1.6 < |#eta^{#mu}| < 2.4 ","lp");
  }
  else{
    leg2->AddEntry(f2,"PbPb ","");
    leg2->AddEntry(f2,"#varepsilon^{DATA}_{#mu}, 1.6 < |#eta^{#mu}| < 2.4 ","lp");
    leg2->AddEntry(fe,"100 T&P stat. variations, 1.6 < |#eta^{#mu}| < 2.4 ","lp");
  }
  leg2->Draw();

  if(!doOnlyNumerator){  cBarrel_pbpb->SaveAs("~/Desktop/Grenelle/TNP_pbpb_100variations.pdf");}
  else{ cBarrel_pbpb->SaveAs("~/Desktop/Grenelle/numerator_variations.pdf");}

  // /// ...the crunchy ones... // ///
  cout << "OMG THIS IS SO PERVERTED!" << endl;
  pbpb=0;
  endcap=0;
  j=0;
 
  TCanvas *cBarrel_pp = new TCanvas ("cBarrel_pp","cBarrel_pp",1400,700);
  cBarrel_pp->Divide(2);
  cBarrel_pp->cd(1);
  TF1 *g1 = new TF1("myFunc",myFunc,0,20,6);
 if (doOnlyNumerator){ g1= new TF1("myFunk",myNumerator,0,20,6);}
  cout << j << endl;
  cout << endcap << endl;
  g1->SetParameters(alpha_pp_midrap[j],mu_data_pp_midrap[j],sigma_data_pp_midrap[j],j,endcap,pbpb);
  g1->SetRange(0,20);
  if(!doOnlyNumerator){   g1->GetYaxis()->SetRangeUser(0.8,1.5); }
  else if(doOnlyNumerator){   g1->GetYaxis()->SetRangeUser(0.0,1.5); }
  g1->Draw();
  // std::cout << alpha_pp_midrap[j] <<"*TMath::Erf((x-"<< mu_data_pp_midrap[j] <<"/"<<sigma_data_pp_midrap[j]<<")) / (1.0000*TMath::Erf((x-1.7960)/2.5160))"<<std::endl;
  
  for(j=1;j< 101;j++){
    //reminder:
    // p[0] is the normalisation;
    // p[1] is the mu_data;
    // p[2] is the sigma_data;
    // p[3] is the index in the array : 0 is the nominal case, 1-101 are the 100 variations.
    // p[4] is the Abs(eta) region : 0 is barrel, 1 is endcap muons.
    // p[5] is the collision type : 0 is pp, 1 is pp.
    gb = new TF1("myfunc",myFunc,0,20,6);
    if (doOnlyNumerator){ gb= new TF1("myFunk",myNumerator,0,20,6);}
    //  gb->SetRange(0,20);
    std::cout  << alpha_pp_midrap[j] <<"*TMath::Erf((x-"<< mu_data_pp_midrap[j] <<"/"<<sigma_data_pp_midrap[j]<<")) / (0.9576*TMath::Erf((x-1.7883)/2.6583))"<<std::endl;
    gb->SetParameters(alpha_pp_midrap[j],mu_data_pp_midrap[j],sigma_data_pp_midrap[j],j,endcap,pbpb);
    gb->SetLineColor(kAzure+1);
    gb->Draw("same");
   }
  g1->Draw("same");
  TLegend *legend = new TLegend(0.2,0.75,0.9,0.9);
  legend->SetTextSize(0.038);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  if(!doOnlyNumerator){ 
    legend->AddEntry(g1,"pp ","");
    legend->AddEntry(g1,"Nominal T&P scaling: |#eta^{#mu}| < 1.6 ","lp");
    legend->AddEntry(gb,"100 T&P stat. variations: |#eta^{#mu}| < 1.6 ","lp");
  }
  else{
    legend->AddEntry(g1,"pp ","");
    legend->AddEntry(g1,"#varepsilon^{DATA}_{#mu}, |#eta^{#mu}| < 1.6 ","lp");
    legend->AddEntry(gb,"100 T&P stat. variations: |#eta^{#mu}| < 1.6 ","lp");
  } 
  legend->Draw();

  //endcap pp:  (0.8299*TMath::Erf((pt-1.2785)/1.8833))/(0.7810*TMath::Erf((pt-1.3609)/2.1231))
  cBarrel_pp->cd(2);
  endcap=1;
  cout << "MOFO DOING ENDCAP NOW" << endl;
  TF1 *g2 = new TF1("myFunc",myFunc,0,20,6);
  if (doOnlyNumerator){ g2= new TF1("myFunk",myNumerator,0,20,6);}
  j=0; // reset, mofo!
  g2->SetParameters(alpha_pp_fwdrap[j],mu_data_pp_fwdrap[j],sigma_data_pp_fwdrap[j],j,endcap,pbpb);
  g2->SetRange(0,20);
  if(!doOnlyNumerator){   g2->GetYaxis()->SetRangeUser(0.8,1.5); }
  else if(doOnlyNumerator){   g2->GetYaxis()->SetRangeUser(0.0,1.5); }
  g2->Draw();
  //  std::cout << alpha_pp_fwdrap[j] <<"*TMath::Erf((x-"<< mu_data_pp_fwdrap[j] <<"/"<<sigma_data_pp_fwdrap[j]<<")) / (1.0000*TMath::Erf((x-1.7960)/2.5160))"<<std::endl;
  
  for(j=1;j< 101;j++){
    //reminder:
    // p[0] is the normalisation;
    // p[1] is the mu_data;
    // p[2] is the sigma_data;
    // p[3] is the index in the array : 0 is the nominal case, 1-101 are the 100 variations.
    // p[4] is the Abs(eta) region : 0 is barrel, 1 is endcap muons.
    // p[5] is the collision type : 0 is pp, 1 is pp.
    ge = new TF1("myfunc",myFunc,0,20,6);
    if (doOnlyNumerator){ ge= new TF1("myFunk",myNumerator,0,20,6);}
    ge->SetRange(0,20);
    std::cout  << alpha_pp_fwdrap[j] <<"*TMath::Erf((x-"<< mu_data_pp_fwdrap[j] <<"/"<<sigma_data_pp_fwdrap[j]<<")) / (0.9576*TMath::Erf((x-1.7883)/2.6583))"<<std::endl;
    // cout << f << endl;
    ge->SetParameters(alpha_pp_fwdrap[j],mu_data_pp_fwdrap[j],sigma_data_pp_fwdrap[j],j,endcap,pbpb);
    ge->SetLineColor(kAzure+1);
    ge->Draw("same");
  }
  g2->Draw("same");
  TLegend *leg2 = new TLegend(0.2,0.75,0.9,0.9);
  leg2->SetTextSize(0.038);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
    if(!doOnlyNumerator)
      { 
	leg2->AddEntry(g2,"pp ","");
	leg2->AddEntry(g2,"Nominal T&P scaling: 1.6 < |#eta^{#mu}| < 2.4 ","lp");
	leg2->AddEntry(ge,"100 T&P stat. variations:1.6 < |#eta^{#mu}| < 2.4 ","lp");
      }
    else{ 
      leg2->AddEntry(g2,"pp ","");
      leg2->AddEntry(g2,"#varepsilon^{DATA}_{#mu}, 1.6 < |#eta^{#mu}| < 2.4 ","lp");
      leg2->AddEntry(ge,"100 T&P stat. variations, 1.6 < |#eta^{#mu}| < 2.4 ","lp");
    }
    leg2->Draw();
    
     if(!doOnlyNumerator)
       {    cBarrel_pp->SaveAs("~/Desktop/Grenelle/TNP_pp_100variations.pdf");}
     else{	 cBarrel_pp->SaveAs("~/Desktop/Grenelle/numerator_variations_pp.pdf");}
    cout << "OMG THIS IS SO PERVERTED!" << endl;
 
}

