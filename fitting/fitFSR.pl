#!/opt/star/bin/perl -w
#
use strict;
#naming the directories: 
my $workDir      = "/Users/nicolas/Project/ups2013/code/";
# after $workdir, the structure is
# pdfOutput                        |txtOutput
# date (I use /MMDD)
# vP_0p01       |    vP_0p05       |vP_0p01     |   vP_0p05
#pt_3 | pt_3p5 | pt_4 | pt_4p5 (in each vP subdirectory)

my $pdfPrefix = "pdfOutput/MC_FSR/";
my $txtPrefix = "txtOutput/MC_FSR/";
my $outFigDirStep = $workDir.$pdfPrefix;
my $outDatDirStep = $workDir.$txtPrefix;
# my $outFigDir;
# my $outDatDir;	
my $outFigDir=$outFigDirStep;
my $outDatDir=$outDatDirStep;	

if(! -e $outFigDirStep){ system("mkdir $outFigDirStep");}
if(! -e $outDatDirStep){ system("mkdir $outDatDirStep");}
print "Welcome to SuperFitter! ";
###
###
#------------------------ OUTPUT DIRECTORIES per samples
my @sampleIndex  = ("","1"); #Input data sample. 1: pp@7TeV-d0 data ; 2: pp@7TeV-d3 data; 3: PbPb@276 regit 4: pp reco TT 5: pp reco GG 6: zhen tree, GG 7: pp@2.76TeV data 8: MC, 9:pPb5tevfirst part of the run.
my $samstart = 1;
my $samend   = 1;
my @choseSamples = ("","Pyquen");##pp2p76tev
# ------------------------------------------------------
my $doNarrowMass = 0; 
my $choseFitParams = 0; #what you want 0: (1s, 2s, 3s) 1: (1s, 2s/1s; 3s/1s); 2: (1S, (2s+3s)/1s); 3:(1S, 2S/1S, 3S/1S, 3S/2S)
my @choseWhat2Fit = ("MCpars");
my $prefix          = $choseWhat2Fit[$choseFitParams];

my $ptMuStart1 = 2;#single muon pt > 4GeV/c cut, that no-one wants to change... that makes me sad.
my $ptMuEnd1 = 2;
my $ptMuStart2 = 3;
my $ptMuEnd2 = 3;
my @chooseptMu = ("","3","3.5","4","4.5");
my @outFigPrefix_ptMu = ("","pt_3","pt_3p5","pt_4","pt_4p5");
my $ptMuPrefix;

my $vProbStart = 1;
my $vProbEnd = 1; #1
my @choosevProb = ("","0.01","0.05");
my @outFigPrefix_vProb = ("","vP_0p01","vP_0p05");
my $vProbPrefix;
# chose if fix sigma1:
my @sigma1  = ("0","1"); # fix or not sigma1
my $sigstart = 0;
my $sigend   = 0;

my @useRef   = ("0","1","2"); # 0 free; 1:data-driven estimation, 2:MC study, 3:old subscript.
my $refstart  = 0;
my $refend    = 0;

my $centstart =1;
my $centend   =1;
#my @centrality_min = ("0","4", "8", "0", "50", "20","0");
#my @centrality_max = ("4","8", "24", "8","100", "40","0");
# ----------------------0----1---2---3----4----5----6----7---8---9---10---11--12---13--; 0->7 1S binning, 8->11 2S binning, n>11 for trials
 my @centrality_min = ("0", "0","2","4", "8","12","16","20","28","0","4", "8","0","20"); 
 my @centrality_max = ("40","2","4","8","12","16","20","28","40","4","8","40","8","40");

my $isHI=1; ### turn to one ONLY when fitting Pyquen samples!
# 0-5 | 5-10 | 10-20 | 20-30 | 30-40 | 40-50(60) | 50(60)-100 | 0-20 | 20-90(alice comparisons)

my $modelstart =0;
my $modelend   =0;
#Background Model.  1: LS erf*exp + pol2; 2: LS RookeyPdf + pol2; 3: erf*exp; 4: pol2; 5: erf*exp+pol2 6: poly 3
my @bkgModel        = ("0","1","2","3","4","5","6");

my $sigModelstart =4;
my $sigModelend   =4;
##### only for MC!!!!
#signal Model. 1: CB; 2:Gaus; 3:CB+Gaus; 4: CB+CB 5:Gaus+Gaus
my @sigModel        = ("0","1","2","3","4","5");

# chose FSR param
#0: free; 1: both fixed 2: alpha fixed 3: npow fixed 
#changing the landscape here...
# fixing parameters one by one, for fun (considering both widths as one )
#0: all free, 1: all fixed, 2: alpha, 3: npow, 4: widths fixed, 5: sigmaFraction fixed.
my @fsr     = ("0","1","2","3","4","5");
my $fsrstart = 0;
my $fsrend   = 0;
### if($sigModelstart!=4 || $sigModelend!=4){$fsrend=3;} ###big danger here!!!! ---- CAREFUL -----
# fitting settings and plotting
my $paramsOnFigure = 0; #1: plot parameters;   0: plot CMS label
my $plotBkg        = 1; #0: hide LS or trkRot; 1: plot LS or trkRot data points and fit lines;
my $doTrkRot       = 0; #0: use LS;   1: use track rotation
my $doConstrainFit = 0; # 1: use constrain method

my $muonEtaMin     = -2.4;
my $muonEtaMax     = 2.4; 
my $dimuYMin       = -2.4; 
my $dimuYMax       = 2.4; 
# my $dimuYMin       = -1.93; 
# my $dimuYMax       = 1.93; #used when not binning in y_ups
my $upsPtCutMin    = 0;
my $upsPtCutMax    = 100; #used when not binning in pT
my $isam;

# my @rapBinMin = ("-2.4","-1.6","-0.8","0","0.8","1.6");
# my @rapBinMax = ("-1.6","-0.8","0","0.8","1.6","2.4");

# ----------------0----1-------2------3----4-----5-----6-----7--# rapidity binning suited for PbPb 1S (0->4), 2S (5->7);
# my @rapBinMin = ("0.","0.4" ,"0.7","1.0","1.5","1.2" ,"0.","0.");
# my @rapBinMax = ("0.4","0.7","1.0","1.5","2.4","2.4","1.2","2.4");
# new rap bins
my @rapBinMin = ("0.","0.4" ,"0.8","1.2","1.6","2.","0" ,"1.2.","0.");
my @rapBinMax = ("0.4","0.8","1.2","1.6","2.","2.4","1.2","2.4","2.4");
# my @rapBinMin = ("-2.465","-2.065","-1.665","-1.265","-0.865","-0.465","-0.065","0.335","0.735","1.135","1.535");       #### pPb
# my @rapBinMax = ("-2.065","-1.665","-1.265","-0.865","-0.465","-0.065","0.335","0.735","1.135","1.535","1.935");      #### pPb
#my $whatBin=0;
my  $rapStart=0;
my  $rapEnd=8;
my  $ptStart=8;
my  $ptEnd=9;
# pT binning for 2S:
#my @upsPtBinMin = ("0","6.5","10","0");	 
#my @upsPtBinMax = ("6.5","10","20","20");
# # ------------------0-----1-----2---3----4----5-----6-----7------8---9## pT binning for 1S (0->5), 2S (6->9);
# my @upsPtBinMin = ("0","2.5","5","8" ,"12","20" ,"0"  ,"6.5","10","0");	 
# my @upsPtBinMax = ("2.5","5" , "8","12","20","50","6.5","10","20","50");
# ------------------0-----1-----2---3----4----5-----6-----7------8---9## pT binning for 1S (0->X), 2S (Y->Z);
my @upsPtBinMin = ("0","2.5","5","8" ,"12","20" ,"0"  ,"5","12","0");
my @upsPtBinMax = ("2.5","5" , "8","12","20","50","5","12","20","20");
## (in theory, I could also do pt 20-50 GeV...)
my $dontDoRapNow =1;
my $dontDoPtNow =0;
#loop for mkdir purposes
 for (my $iPtMu1=$ptMuStart1; $iPtMu1<=$ptMuEnd1; $iPtMu1++)
    {
	$ptMuPrefix = $outFigPrefix_ptMu[$iPtMu1];
	$outFigDir = $outFigDirStep.$ptMuPrefix."/";
	$outDatDir = $outDatDirStep.$ptMuPrefix."/";

	if(! -e $outFigDir){ system("mkdir $outFigDir");}
	if(! -e $outDatDir){ system("mkdir $outDatDir");}
    }
my $firstloopStart;
my $firstloopEnd;
my $secondloopStart;
my $secondloopEnd;
if($dontDoRapNow==1 && $dontDoPtNow==1){
  $firstloopStart = $centstart;
  $firstloopEnd = $centend;
  $secondloopStart =0;
  $secondloopEnd = 0;
}

if($dontDoRapNow==1 && $dontDoPtNow==0){
   $firstloopStart = $ptStart;
   $firstloopEnd = $ptEnd;
   $secondloopStart = 0;
   $secondloopEnd = 0;
}

if($dontDoRapNow==0 && $dontDoPtNow==1){
    $secondloopStart = $rapStart;
    $secondloopEnd = $rapEnd;
    $firstloopStart =0;
    $firstloopEnd   =0;
}

if($dontDoRapNow==0 && $dontDoPtNow==0){
     $firstloopStart=  $ptStart;        
     $firstloopEnd  =  $ptEnd;	         
     $secondloopStart=  $rapStart;
     $secondloopEnd  =  $rapEnd;	
}
for( my $icent=$centstart; $icent <=$centend; $icent++)  
  {
    for ( $isam=$samstart; $isam<=$samend; $isam++)
      {
#	$doNarrowMass = 0; 
	if ($isam==1) {$doNarrowMass=1;}
	#loop over vertex probabilities
	# these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
	# for (my $ivProb=$vProbStart; $ivProb<=$vProbEnd; $ivProb++)
	# {
	#     $vProbPrefix = $outFigPrefix_vProb[$ivProb];
	#     my $vProb = $choosevProb[$ivProb];
	#     # these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
	#     my $outFigDirProb = $outFigDirStep.$vProbPrefix."/";
	#     my $outDatDirProb = $outDatDirStep.$vProbPrefix."/";
    
	    for (my $iPtMu1=$ptMuStart1; $iPtMu1<=$ptMuEnd1; $iPtMu1++)
	    {
		my $muonPtCut1      = $chooseptMu[$iPtMu1]; #single muon pT cut, loose one,
		for (my $iPtMu=$ptMuStart2; $iPtMu<=$ptMuEnd2; $iPtMu++)
		{
		    my $muonPtCut2      = $chooseptMu[$iPtMu]; #single muon pT cut, tight one!
		   # these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
		    $ptMuPrefix = $outFigPrefix_ptMu[$iPtMu1];	
		  #  these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
		    $outFigDir = $outFigDirStep.$ptMuPrefix."/";
		    $outDatDir = $outDatDirStep.$ptMuPrefix."/";
		    #loop over fixed|unfixed parameters for systematics
		    for (my $isig=$sigstart; $isig<=$sigend; $isig++)
		    {
			for (my $ifsr=$fsrstart; $ifsr<=$fsrend; $ifsr++)
			{
			    for (my $iref=$refstart; $iref<=$refend; $iref++)
			    {
				# for systm. studies
				next if($iref==0 && ($isig!=0 && $ifsr!=0));
				#  next if($iref==1);# do not fix to MC 
				#	next if(($iref==1 || $iref==2) && ($isig==0 && $ifsr==0));
				
				for ( my $ibkg=$modelstart; $ibkg<=$modelend; $ibkg++)
				{
				    my $myfsr   = $fsr[$ifsr];
				    my $myref   = $useRef[$iref];
				    my $bkgType = $bkgModel[$ibkg];
				    
				    my $mysigma1= $sigma1[$isig];
				    
				    if ($ibkg == 1 || $ibkg==2) {$plotBkg=1;}
				    else {$plotBkg=0;}
				    
				    my $choseSample = $sampleIndex[$isam];
				    my $sample = $choseSamples[$choseSample];
				    # next: bin loops
				    if($dontDoRapNow==1 && $dontDoPtNow==1)
				    {
					for(my $i1stLoop=$firstloopStart; $i1stLoop<=$firstloopEnd; $i1stLoop++)
					{
					    my $centMin = $centrality_min[$i1stLoop];
					    my $centMax = $centrality_max[$i1stLoop];
					    $upsPtCutMin=0.0;
					    $upsPtCutMax=50.0;
					    $dimuYMin=0.0;
					    $dimuYMax=2.4;
					
					    print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
					    print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
					    print " Sample is: $sample\n";
					    # print  "vProb is: $vProb\n";
					    print   "single Muon 1 pT > $muonPtCut1\n";
					    print   "single Muon 2 pT > $muonPtCut2\n";
					    print    " $upsPtCutMin < pT_Ups < $upsPtCutMax \n";
					    print     " $dimuYMin < y_Ups < $dimuYMax \n";
					    
					    
					    for (my$iMod=$sigModelstart; $iMod<=$sigModelend ; $iMod++)
					    {
						system("root -l -b -q	'fitBinnedHistograms.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$upsPtCutMin,$upsPtCutMax,$muonPtCut1,$muonPtCut2,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$dimuYMin,$dimuYMax,$isHI,$iMod,$dontDoPtNow,$dontDoRapNow)'");   #$dimuYMin,$dimuYMax,
					    }
					    
					}
				    }
				    if($dontDoRapNow==1 && $dontDoPtNow==0){  ## pt loop
					my $centMin = $centrality_min[0];
					my $centMax = $centrality_max[0];
					$dimuYMin=0.0;
					$dimuYMax=2.4;
					 # $dimuYMin=-1.46;
					 # $dimuYMax=2.39;
					for(my $i1stLoop=$firstloopStart; $i1stLoop<=$firstloopEnd; $i1stLoop++)
					{
					    $upsPtCutMin=$upsPtBinMin[$i1stLoop];
					    $upsPtCutMax=$upsPtBinMax[$i1stLoop];
					    print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
					    print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
					    print " Sample is: $sample\n";
					    # print  "vProb is: $vProb\n";
					    print   "single Muon 1 pT > $muonPtCut1\n";
					    print   "single Muon 2 pT > $muonPtCut2\n";
					    print "$upsPtCutMin > pT_Ups > $upsPtCutMax";
					    print "$dimuYMin > y_Ups > $dimuYMax";
					 	    for (my$iMod=$sigModelstart; $iMod<=$sigModelend ; $iMod++)
					    {
						system("root -l -b -q	'fitBinnedHistograms.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$upsPtCutMin,$upsPtCutMax,$muonPtCut1,$muonPtCut2,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$dimuYMin,$dimuYMax,$isHI,$iMod,$dontDoPtNow,$dontDoRapNow)'");   #$dimuYMin,$dimuYMax,
					    }
					}
				    } 
				    if($dontDoRapNow==0 && $dontDoPtNow==1){  ## rap loop
					my $centMin = $centrality_min[0];
					my $centMax = $centrality_max[0];
					$upsPtCutMin=0.0;
					$upsPtCutMax=50.0;
					for(my $i1stLoop=$secondloopStart; $i1stLoop<=$secondloopEnd; $i1stLoop++)
					{
				
					    $dimuYMin=$rapBinMin[$i1stLoop];
					    $dimuYMax=$rapBinMax[$i1stLoop];
					    print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
					    print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
					    print " Sample is: $sample\n";
					    # print  "vProb is: $vProb\n";
					    print   "single Muon 1 pT > $muonPtCut1\n";
					    print   "single Muon 2 pT > $muonPtCut2\n";
					    print "$upsPtCutMin > pT_Ups > $upsPtCutMax";
					    print "$dimuYMin > y_Ups > $dimuYMax";
					   	    for (my$iMod=$sigModelstart; $iMod<=$sigModelend ; $iMod++)
					    {
						system("root -l -b -q	'fitBinnedHistograms.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$upsPtCutMin,$upsPtCutMax,$muonPtCut1,$muonPtCut2,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$dimuYMin,$dimuYMax,$isHI,$iMod,$dontDoPtNow,$dontDoRapNow)'");   #$dimuYMin,$dimuYMax,
					    }
					    
					}
				    }

				    if($dontDoRapNow==0 && $dontDoPtNow==0){
					my $centMin = $centrality_min[0];
					my $centMax = $centrality_max[0];
					for(my $i1stLoop=$firstloopStart; $i1stLoop<=$firstloopEnd; $i1stLoop++)
					{   ### pt loop
					    # $dimuYMin=-1.93;
					    # $dimuYMax=1.93;
					     $dimuYMin=0.0;
					     $dimuYMax=2.4;
					    $upsPtCutMin=$upsPtBinMin[$i1stLoop];
					    $upsPtCutMax=$upsPtBinMax[$i1stLoop];
					    print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
					    print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
					    print " Sample is: $sample\n";
					    # print  "vProb is: $vProb\n";
					    print   "single Muon 1 pT > $muonPtCut1\n";
					    print   "single Muon 2 pT > $muonPtCut2\n";
					    print "$upsPtCutMin > pT_Ups > $upsPtCutMax";
					    print "$dimuYMin > y_Ups > $dimuYMax";
					     for (my$iMod=$sigModelstart; $iMod<=$sigModelend ; $iMod++)
					     {
						 system("root -l -b -q	'fitBinnedHistograms.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$upsPtCutMin,$upsPtCutMax,$muonPtCut1,$muonPtCut2,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$dimuYMin,$dimuYMax,$isHI,$iMod,$dontDoPtNow,$dontDoRapNow)'");   #$dimuYMin,$dimuYMax,
					     }
					    
					}

					
					for(my $i2ndLoop=$secondloopStart; $i2ndLoop<=$secondloopEnd; $i2ndLoop++)
					{   ### pt loop
					    $dimuYMin=$rapBinMin[$i2ndLoop];
					    $dimuYMax=$rapBinMax[$i2ndLoop];
					    $upsPtCutMin=0.0;
					    $upsPtCutMax=50.0;
					    print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
					    print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
					    print " Sample is: $sample\n";
					    # print  "vProb is: $vProb\n";
					    print   "single Muon 1 pT > $muonPtCut1\n";
					    print   "single Muon 2 pT > $muonPtCut2\n";
					    print    " $upsPtCutMin < pT_Ups < $upsPtCutMax \n";
					    print     " $dimuYMin < y_Ups < $dimuYMax \n";
					    

					    for (my$iMod=$sigModelstart; $iMod<=$sigModelend ; $iMod++)
					    {
						system("root -l -b -q	'fitBinnedHistograms.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$upsPtCutMin,$upsPtCutMax,$muonPtCut1,$muonPtCut2,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$dimuYMin,$dimuYMax,$isHI,$iMod,$dontDoPtNow,$dontDoRapNow)'");   #$dimuYMin,$dimuYMax,
					    }
					}
				    }### funky bin conditions
				}	#bkg 
			    } # ref
			}#fsr
		    }# sig
		}#muon pt cut1
	    }#muon pt cut2
      }#sample
  }#centrality loop
