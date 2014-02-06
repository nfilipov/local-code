 {
   //integers for binning.
   //centrality binning in 11-011
   const int bin  = 8;
  float cent[8]  ={22.1, 86.3, 130.0, 187.1, 261.4, 329.4, 381.3}; // with 40-50 and 50-100
   int bin1 = 7;
   int binPt = 5;
   int binPt2010 = 3;
   int binRap = 5;
   int binRap2010 = 2;
   int bin2010 = 4;
  
   // float cent1[8] ={14.2, 69.9, 130.0, 187.1, 261.4, 329.4, 381.3}; // with 40-60 and 60-100
   float centReverse2013[8]={381.2,330.3,261.4,187.3,130.1,86.3,42.02,8.75}; //reverse order
   float cent1[9] ={8.75,42.02, 86.3, 130.1, 187.3, 261.4, 330.3, 381.2}; // with 70-100, 50-70, 40-50, 30-40, 20-30, 10-20, 5-10, 0-5
   float cent2[7] ={17.8, 69.9, 130.0, 187.1, 261.4, 355.4};//??

   float cent2010[5]={308.6,64.24,261.3,355.7};//0-20,20-100,10-20,0-10

   float pt [5] = {1.25, 3.75, 6.5, 10., 16.};
   float pte[5] = {1.25, 1.25, 1.5, 2., 4.};
   float pt2010 [3] = {3.25,8.25,15.};
   float pt2010e[3] = {3.25,1.75,5};

   float rap2010[2]={0.,1.8};
   float rap2010e[2]={0.6,0.6};
   float rap2010paper[2]={0.64,1.54};
   float rap2010paperel[2]={0.64,0.34};
   float rap2010papereh[2]={0.56,0.86};
   float rap[5] = {0.2 , 0.55, 0.85, 1.25, 1.95};
   float rape[5]= {0.2, 0.15, 0.15, 0.25, 0.45};
   float centErr[8]={6,6,6,6,6,6,6};
   float centnoErr[8]={0,0,0,0,0,0,0,0};


   float taa2013[8]={25.90148438,//reverse order
		     20.47,
		     14.47769531,
		     8.782976563,
		     5.089226563,
		     2.748355469,
		     0.982564453,
		     0.125328724};

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
   if(do3p5)
     {
       float massRes1_MC_pt[5]={0.0702278,
     			       0.0704769,
     			       0.0723508,
     			       0.072539,
     			       0.0741795};
       float massRes2_MC_rap[5]={0.0987689,
     			    	0.10944,
     			    	0.12396,
     			    	0.149925,
     			    	0.181642};
       float massRes1_MC_pte[5]={0.000516951,
     			    	0.000501152,
     			    	0.000750781,
     			    	0.000938769,
     			    	0.00117969};
       float massRes2_MC_rape[5]={0.00255286,
     			    	 0.0042901,
     			    	 0.00416989,
     			    	 0.00394299,
     			    	 0.00285552};
       float massRes2_MC_pt[5]={0.135703,
     			    	0.139709,
     			    	0.140769,
     			    	0.141848,
     			    	0.14674};
       float massRes1_MC_rap[5]={0.055119,
     			    	 0.0654371,
     			    	 0.0789018,
     			    	 0.0947856,
     			    	 0.113839};
       float massRes2_MC_pte[5]={0.00116711,
     			    	 0.00114575,
     			    	 0.00185437,
     			    	 0.00228131,
     			    	 0.00339206};
       float massRes1_MC_rape[5]={0.000377702,
     			    	  0.000710405,
     			    	  0.000854388,
     			    	  0.000854612,
     			    	  0.00108084};
       float npowMC_pt[5]={1.63403,
     			   1.52982,
     			   1.53018,
     			   1.3078,
     			   1.35339};
       float npowMC_rap[5]={1.13925,
     			    1.22514,
     			    1.20287,
     			    1.29962,
     			    1.23529};
       float npowMC_pte[5]={0.0437184,
     			    0.0468337,
     			    0.0614992,
     			    0.072707,
     			    0.10773};
       float npowMC_rape[5]={0.0259398,
     			     0.0369573,
     			     0.0451537,
     			     0.0478129,
     			     0.064187};
       float alphaMC_pt[5]={ 1.75943,
     			     1.87361,
     			     1.82563,
     			     1.93345,
     			     1.93388};
       float alphaMC_rap[5]={1.98463,
     			     1.97219,
     			     2.01558,
     			     1.9974,
     			     2.06177};
       float alphaMC_pte[5]={0.0181543,
     			     0.0204457,
     			     0.0270012,
     			     0.0363508,
     			     0.0522739};
       float alphaMC_rape[5]={0.0165051,
			      0.0208522,
			      0.0237572,
			      0.0217666,
			      0.0265274};

       float sigmaFraction_rap[5]={ 0.10950, 0.12203, 0.13997,  0.15214, 0.28186};
       float sigmaFraction_rape[5]={0.0133, 0.0248, 0.0278, 0.0229, 0.0217 };
       float sigmaFraction_pt[5]={0.68306,0.6819,0.70786,0.70281,0.73227};
       float sigmaFraction_pte[5]={0.0105,0.00966,0.0148,0.0182,0.0221};

       float signal1S_pp_pt[5]={1650.85,
				1538.88,
				1038.01,
				626.173,
				348.501};
       float signal1S_pp_rap[5]={1067.23,
				 820.513,
				 778.079,
				 1244.95,
				 1440.24};
       float signal1S_pp_pte[5]={56.8967,
				 51.8017,
				 41.577,
				 32.0563,
				 22.6628};
       float signal1S_pp_rape[5]={38.9423,
				  35.293,
				  35.8386,
				  49.3019,
				  59.3885};
       float signal1S_aa_pt[5]={708.057,
				720.23,
				551.19,
				372.908,
				176.375};
       float signal1S_aa_rap[5]={517.073,
				 404.956,
				 382.291,
				 693.895,
				 598.465};
       float signal1S_aa_pte[5]={46.0242,
				 46.883,
				 42.303,
				 34.6487,
				 18.6638};
       float signal1S_aa_rape[5]={34.2973,
				  31.8904,
				  34.3959,
				  48.881,
				  47.1226};
  float signal_MC_pt[5]={255563,
			     242205,
			     130575,
			     76180.9,
			     40568.1};
       float signal_MC_rap[5]={160331,
			       118193,
			       114100,
			       181607,
			       179870}; //error on MC signal fits is roughly sqrt{signal}!
    
          

       float signal1S_aa_npart[8]={392.678,
				   220.287,
				   459.188,
				   324.863,
				   198.786,
				   117.616,
				   129.648,
				   31.6466};
       float signal2S_aa_npart[8]={31.2578,
				   27.1748,
				   25.9331,
				   53.1985,
				   15.4242,
				   11.3462,
				   9.958,
				   5.48026};
       float signal2S_aa_nparte[8]={1.90106,	
				    0.036798,
				    0.782569,
				    0.055619,
				    1.89971,	
				    1.33171,	
				    0.68534,	
				    0.470168};
       float signal1S_aa_nparte[8]={38.3359,
				    24.2698,
				    35.4018,
				    27.8966,
				    21.341 ,
				    16.2892,
				    15.6146,
				    8.20734 };


 
 
     }
   if (do4)
     { float sigmaFraction_rap[5]={0.07687, 0.097598, 0.169 , 0.12379, 0.27256};
       float sigmaFraction_rape[5]={0.0104,0.0216,0.0361,0.025,0.0255}; 
       float sigmaFraction_pt[5]={0.66387,0.67881,0.723,0.69545,0.71832};
       float sigmaFraction_pte[5]={0.0122,0.01250,0.0162,0.0202,0.0246};

       float massRes1_MC_pt[5]={0.0683712,
				0.0694958,
				0.0717929,
				0.0713796,
				0.0734393};
       float massRes2_MC_rap[5]={0.107213,
				 0.112244,
				 0.12054,
				 0.155516,
				 0.18437};
       float massRes1_MC_pte[5]={0.000637534,
				 0.000654816,
				 0.000821996,
				 0.00103575,
				 0.00131407};
       float massRes2_MC_rape[5]={0.00341234,
				  0.00499705,
				  0.00426943,
				  0.00566898,
				  0.00356049};
       float massRes2_MC_pt[5]={0.137299,
				0.140186,
				0.142773,
				0.140888,
				0.144339};
       float massRes1_MC_rap[5]={0.0552849,
				 0.0651953,
				 0.077617,
				 0.0962206,
				 0.113659};
       float massRes2_MC_pte[5]={0.00137735,
				 0.00150647,
				 0.00222818,
				 0.00247697,
				 0.0035438};
       float massRes1_MC_rape[5]={0.000356612,
				  0.000672031,
				  0.0010924,
				  0.000972713,
				  0.00131744};
       float npowMC_pt[5]={2.07882,
			   1.65775,
			   1.73509,
			   1.44071,
			   1.37867};
       float npowMC_rap[5]={1.25553,
			    1.34141,
			    1.22328,
			    1.5064,
			    1.48334};
       float npowMC_pte[5]={0.0845129,
			    0.0703862,
			    0.0798167,
			    0.0851945,
			    0.112018};
       float npowMC_rape[5]={0.0350327,
			     0.0493107,
			     0.0604923,
			     0.0678428,
			     0.0948897};
       float alphaMC_pt[5]={1.77283,
			    1.8701,
			    1.75824,
			    1.87573,
			    1.92024};
       float alphaMC_rap[5]={1.96474,
			     1.96184,
			     2.06546,
			     1.97356,
			     2.05739};
       float alphaMC_pte[5]={0.0267216,
			     0.0284054,
			     0.0310081,
			     0.0394019,
			     0.0535848};
       float alphaMC_rape[5]={0.0208169,
			      0.0258252,
			      0.0311686,
			      0.0269537,
			      0.0333402};
       float signal1S_pp_pt[5]={1101.69,
				870.149,
				756.12,
				522.075,
				319.197};
       float signal1S_pp_pte[5]={44.0769,
				 38.4254,
				 34.8172,
				 28.9365,
				 21.4794};
       float signal1S_pp_rap[5]={770.861,
				 562.779,
				 539.611,
				 859.788,
				 935.563};
       float signal1S_pp_rape[5]={32.157,
				  28.6306,
				  29.1644,
				  39.232,
				  46.3424};
       float signal1S_aa_pt[5]={460.59,
				421.201,
				357.055,
				297.813,
				156.354};
       float signal1S_aa_pte[5]={34.868,
				 33.3118,
				 31.3543,
				 29.2692,
				 17.3782};
       float signal1S_aa_rap[5]={323.368,
				 284.924,
				 260.601,
				 511.865,
				 357.765};
       float signal1S_aa_rape[5]={26.0745,
				  24.4786,
				  26.0905,
				  37.5099,
				  35.7498};
   float signal_MC_pt[5]={156088,
			      133770,
			      90824,
			      61280.8,
			      36369};
       float signal_MC_rap[5]={104730,
			       76646.8,
			       73101.9,
			       114331,
			       118081}; //error on MC signal fits is roughly sqrt{signal}!

       float signal1S_aa_npart[8]={292.481,
				   157.972,
				   289.354,
				   216.81,
				   145.959,
				   59.4915,
				   84.64,
				   22.5316};
       float signal1S_aa_nparte[8]={29.16,
				    18.8157,
				    26.5306,
				    21.3647,
				    16.5658,
				    11.7635,
				    12.4018,
				    7.35998};
       float signal2S_aa_npart[8]={33.6168,
				   7.34061,
				   12.0699,
				   17.9355,
				   27.4959,
				   8.25104,
				   7.53161,
				   4.14925};
       float signal2S_aa_nparte[8]={22.9018,
				    12.444,
				    19.5554,
				    15.454,
				    12.2451,
				    8.80566,
				    7.99619,
				    5.23672};


     }

   // mean mass resolution in AA "MB", and in pp "pp"
   float massRes_MB[1]={83.7};
   float massRes_MBe[1]={4.0};
   float massRes_pp[1]={91.6};
   float massRes_ppe[1]={2.1};
 
   //alice shit pT>0 GeV/c 
   int binA=2;

   float centAlice[binA]={72,308}; //20-90, 0-20%
   float centAliceErr[binA]={6,6};
   float centAliceNoErr[binA]={0,0};

}


