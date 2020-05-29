using namespace RooFit;                                                                                                           
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

double PDGB0Mass = 5.27958;
double PDGJpsiMass = 3.096916;
double PDGPsiPrimeMass = 3.686109;

void Geteff(int year, int q2Bin, int binnum){
	  // year format: [6] for 2016
  //              [7] for 2017
  //              [8] for 2018
  // q2-bin format: [0-8] for one bin
    // flags to mark which q2 bins should be filled

cout << "year 201" << year << " q2Bin" << q2Bin << " binnum" << binnum << endl;
	bool runBin [nBins];
  string shortString [nBins];
  string longString  [nBins];
  for (int i=0; i<nBins; ++i) {
    runBin [i] = false;
    if ( q2Bin!=-1 && q2Bin!=i ) continue;
    runBin [i] = true;
    shortString [i] = Form("b%i",i);
    longString  [i] = Form("q2 bin %i",i);
  }
  // Load ntuples
  TChain* t_gen = new TChain();
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_gen->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/GEN_BFilter_B0MuMuKstar_p*.root/ntuple");
  if ( year==6 ) {
    // 2016
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/NtupleMay20/2016GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/ntuple_01_10_2019/2016MC_LMNR.root/ntuple");
  } else if ( year==7 ) {
    // 2017
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/2017GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/2017MC_LMNR.root/ntuple");
  } else if ( year==8 ) {
    // 2018
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/2018GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/2018MC_LMNR.root/ntuple");
  }
  int genEntries = t_gen->GetEntries();
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();


    // Import branches from ntuples:
  // angular variables
  double genCosThetaK, genCosThetaL, genPhi;
  double recoCosThetaK, recoCosThetaL, recoPhi;
  t_gen->SetBranchAddress( "cos_theta_k"     , &genCosThetaK  );
  t_gen->SetBranchAddress( "cos_theta_l"     , &genCosThetaL  );
  t_gen->SetBranchAddress( "phi_kst_mumu"    , &genPhi        );
  t_den->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK  );
  t_den->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL  );
  t_den->SetBranchAddress( "gen_phi_kst_mumu", &genPhi        );
  t_num->SetBranchAddress( "cos_theta_k"     , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l"     , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu"    , &recoPhi       );

  // variables for applying GEN-filter
  double genmupEta, genmumEta, genkstTrkpEta, genkstTrkmEta, genmupPt, genmumPt, genkstTrkpPt, genkstTrkmPt;
  t_gen->SetBranchAddress( "genmupEta", &genmupEta );
  t_gen->SetBranchAddress( "genmumEta", &genmumEta );
  t_gen->SetBranchAddress( "genkstTrkpEta", &genkstTrkpEta );
  t_gen->SetBranchAddress( "genkstTrkmEta", &genkstTrkmEta );
  t_gen->SetBranchAddress( "genmupPt", &genmupPt );
  t_gen->SetBranchAddress( "genmumPt", &genmumPt );
  t_gen->SetBranchAddress( "genkstTrkpPt", &genkstTrkpPt );
  t_gen->SetBranchAddress( "genkstTrkmPt", &genkstTrkmPt );

  // dimuon mass variables
  double genDimuMass2, recoDimuMass;
  t_gen->SetBranchAddress( "genQ2"   , &genDimuMass2 );
  t_den->SetBranchAddress( "genQ2"   , &genDimuMass2 );
  t_num->SetBranchAddress( "mumuMass", &recoDimuMass );


    // B0 mass variable
  double recoB0Mass;
  t_num->SetBranchAddress( "tagged_mass", &recoB0Mass );

  // B0-kinematic variables
  // double genB0pT, genB0eta;
  // double recoB0pT, recoB0eta;
  // t_gen->SetBranchAddress( "genbPt" , &genB0pT   );
  // t_gen->SetBranchAddress( "genbEta", &genB0eta  );
  // t_den->SetBranchAddress( "genbPt" , &genB0pT   );
  // t_den->SetBranchAddress( "genbEta", &genB0eta  );
  // t_num->SetBranchAddress( "bPt"    , &recoB0pT  );
  // t_num->SetBranchAddress( "bEta"   , &recoB0eta );

  // flavour tagging variables
  double genSignal, tagB0;
  t_num->SetBranchAddress( "genSignal", &genSignal );
  t_num->SetBranchAddress( "tagB0"    , &tagB0     );

  // event number for even/odd splitting
  double eventN_Dou;
  Long64_t eventN;
  t_gen->SetBranchAddress( "eventN", &eventN_Dou );
  t_den->SetBranchAddress( "eventN", &eventN     );
  t_num->SetBranchAddress( "eventN", &eventN     );

  // event pileup weight
  float PUweight = 1;
  t_den->SetBranchAddress( "weight", &PUweight );
  t_num->SetBranchAddress( "weight", &PUweight );

    // final state radiation flag
  double genSignHasFSR;
  t_gen->SetBranchAddress( "genSignHasFSR", &genSignHasFSR );

TFile *fout = new TFile(Form("1Deffb%i.root",q2Bin),"recreate");

TH1D f1s0("f1s","f1s",binnum,0.,1.);
TH1D m6s0("m6s","m6s",binnum,-1,1.);
TH1D m6c0("m6c","m6c",binnum,-1,1.);
TH1D f30("f3","f3",binnum,-1,1.);
TH1D f40("f4","f4",binnum,-1,1.);
TH1D f50("f5","f5",binnum,-1,1.);
TH1D f70("f7","f7",binnum,-1,1.);
TH1D f80("f8","f8",binnum,-1,1.);
TH1D f90("f9","f9",binnum,-1,1.);



TH1D *f1s[7];
TH1D *m6s[7];
TH1D *m6c[7];
TH1D *f3[7];
TH1D *f4[7];
TH1D *f5[7];
TH1D *f7[7];
TH1D *f8[7];
TH1D *f9[7];

//0 genDen
//1 genNum
//2 den
//3 CTreco
//4 WTreco
//5 CTeff
//6 WTeff


for (int i=0;i<7;i++){
	f1s[i] = new TH1D(f1s0);
	m6s[i] = new TH1D(m6s0);
	m6c[i] = new TH1D(m6c0);
	f3[i] = new TH1D(f30);
	f4[i] = new TH1D(f40);
	f5[i] = new TH1D(f50);
	f7[i] = new TH1D(f70);
	f8[i] = new TH1D(f80);
	f9[i] = new TH1D(f90);
}
//int counter;
double cos_theta_k,cos_theta_l,phi_kst_mumu;
double f_1s,M_6s,M_6c,f_3,f_4,f_5,f_7,f_8,f_9;
// Prepare GEN-level datasets
  cout<<"Starting GEN datasets filling..."<<endl;
  int xBin;
  for (int iCand=0; iCand<genEntries; ++iCand) {
    t_gen->GetEntry(iCand);
	 if (iCand%1000000==0){
	 cout << "entries" << iCand << endl;
	 }
	 cos_theta_k = genCosThetaK;
	 cos_theta_l = genCosThetaL;
	 phi_kst_mumu = genPhi;
	         f_1s = 1 - cos_theta_k * cos_theta_k;
            M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
    // find q2 bin of current candidate
	    if ( ( genDimuMass2 < binBorders[q2Bin+1] ) &&
        ( genDimuMass2 > binBorders[q2Bin]   ) ) {
			    if ( genSignHasFSR<0.5 ) {
			f1s[0]->Fill(f_1s);
         m6s[0]->Fill(M_6s);
         m6c[0]->Fill(M_6c);
         f3[0]->Fill(f_3);
         f4[0]->Fill(f_4);
         f5[0]->Fill(f_5);
         f7[0]->Fill(f_7);
         f8[0]->Fill(f_8);
         f9[0]->Fill(f_9);      
				 }

		    if ( fabs(genmupEta)<2.5 && fabs(genmumEta)<2.5 &&
    fabs(genkstTrkpEta)<2.5 && fabs(genkstTrkmEta)<2.5 &&
    genmupPt>2.5 && genmumPt>2.5 &&
    genkstTrkpPt>0.4 && genkstTrkmPt>0.4){
         f1s[1]->Fill(f_1s);
         m6s[1]->Fill(M_6s);
         m6c[1]->Fill(M_6c);
         f3[1]->Fill(f_3);
         f4[1]->Fill(f_4);
         f5[1]->Fill(f_5);
         f7[1]->Fill(f_7);
         f8[1]->Fill(f_8);
         f9[1]->Fill(f_9); 
			 }



		 }
  }


  // Prepare denominator dataset
  cout<<"Starting denominator dataset filling..."<<endl;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);

	if (iCand%100000==0){
	cout << "entries" << iCand << endl;
	}
	cos_theta_k = genCosThetaK;
    cos_theta_l = genCosThetaL;
    phi_kst_mumu = genPhi;
            f_1s = 1 - cos_theta_k * cos_theta_k;
            M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
   if ( ( genDimuMass2 < binBorders[q2Bin+1] ) &&
        ( genDimuMass2 > binBorders[q2Bin]   ) ){

		   f1s[2]->Fill(f_1s,PUweight);
         m6s[2]->Fill(M_6s,PUweight);
         m6c[2]->Fill(M_6c,PUweight);
         f3[2]->Fill(f_3,PUweight);
         f4[2]->Fill(f_4,PUweight);
         f5[2]->Fill(f_5,PUweight);
         f7[2]->Fill(f_7,PUweight);
         f8[2]->Fill(f_8,PUweight);
         f9[2]->Fill(f_9,PUweight);
	}


  }

  // Prepare numerator dataset
  cout<<"Starting numerator dataset filling..."<<endl;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
   if (iCand%10000==0){
	cout << "entries" << iCand <<endl;
	}
	// anti-radiation cut
    if ( recoDimuMass < PDGJpsiMass ) { // below Jpsi
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGJpsiMass ) < 0.18 ) continue;
    } else if ( recoDimuMass > PDGPsiPrimeMass ) { // above PsiPrime
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGPsiPrimeMass ) < 0.08 ) continue;
    } else { // between the resonances
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGJpsiMass ) < 0.08 ) continue;
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGPsiPrimeMass ) < 0.09 ) continue;
    }
        cos_theta_k = recoCosThetaK;
    cos_theta_l = recoCosThetaL;
    phi_kst_mumu = recoPhi;
            f_1s = 1 - cos_theta_k * cos_theta_k;
            M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);

   if ( ( pow(recoDimuMass,2) < binBorders[q2Bin+1] ) &&
        ( pow(recoDimuMass,2) > binBorders[q2Bin]   ) ){
		    if (genSignal != tagB0+1) { // correctly tagged events
			f1s[3]->Fill(f_1s,PUweight);
         m6s[3]->Fill(M_6s,PUweight);
         m6c[3]->Fill(M_6c,PUweight);
         f3[3]->Fill(f_3,PUweight);
         f4[3]->Fill(f_4,PUweight);
         f5[3]->Fill(f_5,PUweight);
         f7[3]->Fill(f_7,PUweight);
         f8[3]->Fill(f_8,PUweight);
         f9[3]->Fill(f_9,PUweight);
			 }
			 else{
				 M_6s = -M_6s;
				 M_6c = -M_6c;
				 f_5 = -f_5;
				 f_8 = -f_8;
				 f_9 = -f_9;

			f1s[4]->Fill(f_1s,PUweight);
         m6s[4]->Fill(M_6s,PUweight);
         m6c[4]->Fill(M_6c,PUweight);
         f3[4]->Fill(f_3,PUweight);
         f4[4]->Fill(f_4,PUweight);
         f5[4]->Fill(f_5,PUweight);
         f7[4]->Fill(f_7,PUweight);
         f8[4]->Fill(f_8,PUweight);
         f9[4]->Fill(f_9,PUweight);
			 }


	}
  }

double effC[9],effW[9];
f1s[5]->SetName("f1seffC");
f1s[5]->SetTitle("Correctly tagged efficiency of f1s");
m6s[5]->SetName("m6seffC");
m6s[5]->SetTitle("Correctly tagged efficiency of m6s");
m6c[5]->SetName("m6ceffC");
m6c[5]->SetTitle("Correctly tagged efficiency of m6c");
f3[5]->SetName("f3effC");
f3[5]->SetTitle("Correctly tagged efficiency of f3");
f4[5]->SetName("f4effC");
f4[5]->SetTitle("Correctly tagged efficiency of f4");
f5[5]->SetName("f5effC");
f5[5]->SetTitle("Correctly tagged efficiency of f5");
f7[5]->SetName("f7effC");
f7[5]->SetTitle("Correctly tagged efficiency of f7");
f8[5]->SetName("f8effC");
f8[5]->SetTitle("Correctly tagged efficiency of f8");
f9[5]->SetName("f9effC");
f9[5]->SetTitle("Correctly tagged efficiency of f9");


f1s[6]->SetName("f1seffW");
f1s[6]->SetTitle("Wrongly tagged efficiency of f1s");
m6s[6]->SetName("m6seffW");
m6s[6]->SetTitle("Wrongly tagged efficiency of m6s");
m6c[6]->SetName("m6ceffW");
m6c[6]->SetTitle("Wrongly tagged efficiency of m6c");
f3[6]->SetName("f3effW");
f3[6]->SetTitle("Wrongly tagged efficiency of f3");
f4[6]->SetName("f4effW");
f4[6]->SetTitle("Wrongly tagged efficiency of f4");
f5[6]->SetName("f5effW");
f5[6]->SetTitle("Wrongly tagged efficiency of f5");
f7[6]->SetName("f7effW");
f7[6]->SetTitle("Wrongly tagged efficiency of f7");
f8[6]->SetName("f8effW");
f8[6]->SetTitle("Wrongly tagged efficiency of f8");
f9[6]->SetName("f9effW");
f9[6]->SetTitle("Wrongly tagged efficiency of f9");


cout << "begin calculate efficiency" << endl;
double genDen[9],genNum[9],den[9],ctRECO[9],wtRECO[9],cteff[9],wteff[9];
for (int i=1;i<=binnum;i++){
	genDen[0] = f1s[0]->GetBinContent(i);
	genDen[1] = m6s[0]->GetBinContent(i);
	genDen[2] = m6c[0]->GetBinContent(i);
	genDen[3] = f3[0]->GetBinContent(i);
	genDen[4] = f4[0]->GetBinContent(i);
	genDen[5] = f5[0]->GetBinContent(i);
	genDen[6] = f7[0]->GetBinContent(i);
	genDen[7] = f8[0]->GetBinContent(i);
	genDen[8] = f9[0]->GetBinContent(i);


	genNum[0] = f1s[1]->GetBinContent(i);
   genNum[1] = m6s[1]->GetBinContent(i);
   genNum[2] = m6c[1]->GetBinContent(i);
   genNum[3] = f3[1]->GetBinContent(i);
   genNum[4] = f4[1]->GetBinContent(i);
   genNum[5] = f5[1]->GetBinContent(i);
   genNum[6] = f7[1]->GetBinContent(i);
   genNum[7] = f8[1]->GetBinContent(i);
   genNum[8] = f9[1]->GetBinContent(i);


	den[0] = f1s[2]->GetBinContent(i);
   den[1] = m6s[2]->GetBinContent(i);
   den[2] = m6c[2]->GetBinContent(i);
   den[3] = f3[2]->GetBinContent(i);
   den[4] = f4[2]->GetBinContent(i);
   den[5] = f5[2]->GetBinContent(i);
   den[6] = f7[2]->GetBinContent(i);
   den[7] = f8[2]->GetBinContent(i);
   den[8] = f9[2]->GetBinContent(i);

	ctRECO[0] = f1s[3]->GetBinContent(i);
   ctRECO[1] = m6s[3]->GetBinContent(i);
   ctRECO[2] = m6c[3]->GetBinContent(i);
   ctRECO[3] = f3[3]->GetBinContent(i);
   ctRECO[4] = f4[3]->GetBinContent(i);
   ctRECO[5] = f5[3]->GetBinContent(i);
   ctRECO[6] = f7[3]->GetBinContent(i);
   ctRECO[7] = f8[3]->GetBinContent(i);
   ctRECO[8] = f9[3]->GetBinContent(i);


	wtRECO[0] = f1s[4]->GetBinContent(i);
   wtRECO[1] = m6s[4]->GetBinContent(i);
   wtRECO[2] = m6c[4]->GetBinContent(i);
   wtRECO[3] = f3[4]->GetBinContent(i);
   wtRECO[4] = f4[4]->GetBinContent(i);
   wtRECO[5] = f5[4]->GetBinContent(i);
   wtRECO[6] = f7[4]->GetBinContent(i);
   wtRECO[7] = f8[4]->GetBinContent(i);
   wtRECO[8] = f9[4]->GetBinContent(i);


	for (int k=0;k<9;k++){
		cteff[k] = ctRECO[k]*genNum[k]/genDen[k]/den[k];
		wteff[k] = wtRECO[k]*genNum[k]/genDen[k]/den[k];
	}


	f1s[5]->SetBinContent(i,cteff[0]);
	m6s[5]->SetBinContent(i,cteff[1]);
	m6c[5]->SetBinContent(i,cteff[2]);
	f3[5]->SetBinContent(i,cteff[3]);
	f4[5]->SetBinContent(i,cteff[4]);
	f5[5]->SetBinContent(i,cteff[5]);
	f7[5]->SetBinContent(i,cteff[6]);
	f8[5]->SetBinContent(i,cteff[7]);
	f9[5]->SetBinContent(i,cteff[8]);


	f1s[6]->SetBinContent(i,wteff[0]);
   m6s[6]->SetBinContent(i,wteff[1]);
   m6c[6]->SetBinContent(i,wteff[2]);
   f3[6]->SetBinContent(i,wteff[3]);
   f4[6]->SetBinContent(i,wteff[4]);
   f5[6]->SetBinContent(i,wteff[5]);
   f7[6]->SetBinContent(i,wteff[6]);
   f8[6]->SetBinContent(i,wteff[7]);
   f9[6]->SetBinContent(i,wteff[8]);

}

f1s[5]->Write();
f1s[6]->Write();
m6s[5]->Write();
m6s[6]->Write();
m6c[5]->Write();
m6c[6]->Write();
f3[5]->Write();
f3[6]->Write();
f4[5]->Write();
f4[6]->Write();
f5[5]->Write();
f5[6]->Write();
f7[5]->Write();
f7[6]->Write();
f8[5]->Write();
f8[6]->Write();
f9[5]->Write();
f9[6]->Write();

cout << "end year " << year << " bin" << q2Bin << endl;



}
	



