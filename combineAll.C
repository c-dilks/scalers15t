// combines sums.root and rdat_i.root into one final tree, with relative luminosities
// and total scaler counts for each run
//
// output is "rtree.root"

void combineAll(const char * rdatfile="rdat_i.root", const char * sumfile="sums.root")
{
  // open trees
  TFile * rdatfile_tf = new TFile(rdatfile,"READ");
  TTree * rdat = (TTree*) rdatfile_tf->Get("rtr");
  TFile * sumfile_tf = new TFile(sumfile,"READ");
  TTree * sums = (TTree*) sumfile_tf->Get("sum");


  // check if trees have same size
  if(rdat->GetEntries() != sums->GetEntries())
  {
    fprintf(stderr,"ERROR: tree sizes unequal\n");
    return;
  };

  enum detector_enum {kBBC,kZDC,kVPD};
  enum combo_enum {kE,kW,kX};
  enum spinbit_enum {kNN,kNP,kPN,kPP,kALL};

  // trigger bit character strings (tbit)
  char tbit[3][4];
  sprintf(tbit[kBBC],"bbc");
  sprintf(tbit[kZDC],"zdc");
  sprintf(tbit[kVPD],"vpd");


  // combination bit character strings (cbit)
  char cbit[3][4];
  sprintf(cbit[kE],"e");
  sprintf(cbit[kW],"w");
  sprintf(cbit[kX],"x");

  Int_t i_rdat,runnum_rdat,fill_rdat;
  Double_t t_rdat;

  Float_t RR[3][3][10]; // [tbit (bbc,zdc,vpd)] [cbit (e,w,x)] [rellum]
  Float_t RR_rsc[3][10]; 
  Float_t RR_mean[3][10]; // mean R
  Float_t RR_mean_err[3][10]; // error of mean R
  Float_t RR_rsc_err[3][10]; // error of rsc R

  char br_RR[3][3][10][32];
  char br_RR_rsc[3][10][32]; 
  char br_RR_mean[3][10][32];
  char br_RR_mean_err[3][10][32];
  char br_RR_rsc_err[3][10][32];

  char ty_RR[3][3][10][32];
  char ty_RR_rsc[3][10][32]; 
  char ty_RR_rsc_err[3][10][32];
  char ty_RR_mean[3][10][32];
  char ty_RR_mean_err[3][10][32];

  Int_t i_sums,runnum_sums,fill_sums,fi_sums; // (fi not in rdat tree)
  Double_t t_sums,tau,freq;
  Int_t num_runs;
  Double_t bbce,bbcw,bbcx;
  Double_t zdce,zdcw,zdcx;
  Double_t vpde,vpdw,vpdx;
  Double_t tot_bx;
  Int_t pattern;
  Float_t d_vz[3]; // [cbit]
  Float_t d_xx[3][3]; // [xbit] [tbit]
  Bool_t isConsistent; // true for a run if diagnostics are "consistent"

  rdat->SetBranchAddress("i",&i_rdat);
  rdat->SetBranchAddress("runnum",&runnum_rdat);
  rdat->SetBranchAddress("fill",&fill_rdat);
  rdat->SetBranchAddress("t",&t_rdat);

  for(Int_t r=1; r<10; r++)
  {
    for(Int_t t=0; t<3; t++)
    {
      for(Int_t c=0; c<3; c++)
      {
        sprintf(br_RR[t][c][r],"R%d_%s%s",r,tbit[t],cbit[c]);
        sprintf(ty_RR[t][c][r],"%s/F",br_RR[t][c][r]);
        rdat->SetBranchAddress(br_RR[t][c][r],&(RR[t][c][r]));
      };
      sprintf(br_RR_rsc[t][r],"R%d_%srsc",r,tbit[t]);
      sprintf(br_RR_rsc_err[t][r],"R%d_%s_rsc_err",r,tbit[t]);
      sprintf(br_RR_mean[t][r],"R%d_%s_mean",r,tbit[t]);
      sprintf(br_RR_mean_err[t][r],"R%d_%s_mean_err",r,tbit[t]);
      sprintf(ty_RR_rsc[t][r],"%s/F",br_RR_rsc[t][r]);
      sprintf(ty_RR_rsc_err[t][r],"%s/F",br_RR_rsc_err[t][r]);
      sprintf(ty_RR_mean[t][r],"%s/F",br_RR_mean[t][r]);
      sprintf(ty_RR_mean_err[t][r],"%s/F",br_RR_mean_err[t][r]);
      rdat->SetBranchAddress(br_RR_rsc[t][r],&(RR_rsc[t][r]));
      rdat->SetBranchAddress(br_RR_rsc_err[t][r],&(RR_rsc_err[t][r]));
      rdat->SetBranchAddress(br_RR_mean[t][r],&(RR_mean[t][r]));
      rdat->SetBranchAddress(br_RR_mean_err[t][r],&(RR_mean_err[t][r]));
    };
  };


  rdat->SetBranchAddress("d_vz_e",&(d_vz[0]));
  rdat->SetBranchAddress("d_vz_w",&(d_vz[1]));
  rdat->SetBranchAddress("d_vz_x",&(d_vz[2]));
  rdat->SetBranchAddress("d_ew_bbc",&(d_xx[0][0]));
  rdat->SetBranchAddress("d_ew_zdc",&(d_xx[0][1]));
  rdat->SetBranchAddress("d_ew_vpd",&(d_xx[0][2]));
  rdat->SetBranchAddress("d_ex_bbc",&(d_xx[1][0]));
  rdat->SetBranchAddress("d_ex_zdc",&(d_xx[1][1]));
  rdat->SetBranchAddress("d_ex_vpd",&(d_xx[1][2]));
  rdat->SetBranchAddress("d_wx_bbc",&(d_xx[2][0]));
  rdat->SetBranchAddress("d_wx_zdc",&(d_xx[2][1]));
  rdat->SetBranchAddress("d_wx_vpd",&(d_xx[2][2]));

  sums->SetBranchAddress("i",&i_sums);
  sums->SetBranchAddress("runnum",&runnum_sums);
  sums->SetBranchAddress("fi",&fi_sums);
  sums->SetBranchAddress("fill",&fill_sums);
  sums->SetBranchAddress("t",&t_sums);
  sums->SetBranchAddress("freq",&freq);
  sums->SetBranchAddress("tau",&tau);
  sums->SetBranchAddress("num_runs",&num_runs);
  sums->SetBranchAddress("bbce",&bbce);
  sums->SetBranchAddress("bbcw",&bbcw);
  sums->SetBranchAddress("bbcx",&bbcx);
  sums->SetBranchAddress("zdce",&zdce);
  sums->SetBranchAddress("zdcw",&zdcw);
  sums->SetBranchAddress("zdcx",&zdcx);
  sums->SetBranchAddress("vpde",&vpde);
  sums->SetBranchAddress("vpdw",&vpdw);
  sums->SetBranchAddress("vpdx",&vpdx);
  sums->SetBranchAddress("tot_bx",&tot_bx);
  sums->SetBranchAddress("pattern",&pattern);


  
  TFile * outfile = new TFile("rtree.root","RECREATE");
  TTree * rellum = new TTree("rellum","rellum");
  rellum->Branch("i",&i_sums,"i/I"); // run index
  rellum->Branch("runnum",&runnum_sums,"runnum/I"); // run number
  rellum->Branch("fi",&fi_sums,"fi/I"); // fill index
  rellum->Branch("fill",&fill_sums,"fill/I"); // fill number
  rellum->Branch("t",&t_sums,"t/D"); // run time (seconds)
  rellum->Branch("freq",&freq,"freq/D"); // clock frequency
  rellum->Branch("tau",&tau,"tau/D"); // total bXings / bXing rate (nominally == run time)
  rellum->Branch("num_runs",&num_runs,"num_runs/I"); // total no. runs in a fill
  rellum->Branch("bbce",&bbce,"bbce/D"); // multiples & accidentals corrected scaler counts
  rellum->Branch("bbcw",&bbcw,"bbcw/D");
  rellum->Branch("bbcx",&bbcx,"bbcx/D");
  rellum->Branch("zdce",&zdce,"zdce/D");
  rellum->Branch("zdcw",&zdcw,"zdcw/D");
  rellum->Branch("zdcx",&zdcx,"zdcx/D");
  rellum->Branch("vpde",&vpde,"vpde/D");
  rellum->Branch("vpdw",&vpdw,"vpdw/D");
  rellum->Branch("vpdx",&vpdx,"vpdx/D");
  rellum->Branch("tot_bx",&tot_bx,"tot_bx/D"); // total no. bXings


  for(Int_t t=0; t<3; t++)
  {
    for(Int_t c=0; c<3; c++)
    {
      for(Int_t r=1; r<10; r++)
      {
        rellum->Branch(br_RR[t][c][r],&(RR[t][c][r]),ty_RR[t][c][r]);
      };
    };
  };

  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      rellum->Branch(br_RR_rsc[t][r],&(RR_rsc[t][r]),ty_RR_rsc[t][r]);
    };
  };

  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      rellum->Branch(br_RR_rsc_err[t][r],&(RR_rsc_err[t][r]),ty_RR_rsc_err[t][r]);
    };
  };

  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      rellum->Branch(br_RR_mean[t][r],&(RR_mean[t][r]),ty_RR_mean[t][r]);
    };
  };

  for(Int_t t=0; t<3; t++)
  {
    for(Int_t r=1; r<10; r++)
    {
      rellum->Branch(br_RR_mean_err[t][r],&(RR_mean_err[t][r]),ty_RR_mean_err[t][r]);
    };
  };

  rellum->Branch("d_vz_e",&(d_vz[0]),"d_vz_e/F"); // zdc-vpd diagnostic
  rellum->Branch("d_vz_w",&(d_vz[1]),"d_vz_w/F");
  rellum->Branch("d_vz_x",&(d_vz[2]),"d_vz_x/F");
  rellum->Branch("d_ew_bbc",&(d_xx[0][0]),"d_ew_bbc/F"); // east-west diagnostic
  rellum->Branch("d_ew_zdc",&(d_xx[0][1]),"d_ew_zdc/F");
  rellum->Branch("d_ew_vpd",&(d_xx[0][2]),"d_ew_vpd/F");
  rellum->Branch("d_ex_bbc",&(d_xx[1][0]),"d_ex_bbc/F"); // east-coin diagnostic
  rellum->Branch("d_ex_zdc",&(d_xx[1][1]),"d_ex_zdc/F");
  rellum->Branch("d_ex_vpd",&(d_xx[1][2]),"d_ex_vpd/F");
  rellum->Branch("d_wx_bbc",&(d_xx[2][0]),"d_wx_bbc/F"); // west-coin diagnostic
  rellum->Branch("d_wx_zdc",&(d_xx[2][1]),"d_wx_zdc/F");
  rellum->Branch("d_wx_vpd",&(d_xx[2][2]),"d_wx_vpd/F");
  rellum->Branch("isConsistent",&isConsistent,"isConsistent/O"); // true if diagnostics passed
  rellum->Branch("pattern",&pattern,"pattern/I"); // spin pattern no. (see sumTree.C)


  // diagnostic consistency bounds (see if statement in tree for loop below)
  // -- these bounds were determined by eye
  // -- these are so far only done using cdf acc+mul corrections; but not rate-safe corrections;
  //    this can be modified later, but I so far see no reason that it's needed
  // -- currently these diagnostics are only for R3
  Float_t d_vz_cc[3];
  Float_t d_xx_cc[3][3];
  Float_t t_tau;
  d_vz_cc[0] = 1.000; // VPDE - ZDCE // ZDC not reliable for run15? leave cuts open for now (only cut on VPD)
  d_vz_cc[1] = 1.000; // VPDW - ZDCW
  d_vz_cc[2] = 1.000; // VPDX - ZDCX
  d_xx_cc[0][1] = 1.000; // ZDCE - ZDCW
  d_xx_cc[0][2] = 0.002;  // VPDE - VPDW  // <---
  d_xx_cc[1][1] = 1.000;  // ZDCE - ZDCX
  d_xx_cc[1][2] = 0.002;  // VPDE - VPDX  // <---
  d_xx_cc[2][1] = 1.000;  // ZDCW - ZDCX
  d_xx_cc[2][2] = 0.003;  // VPDW - VPDX  // <---
  t_tau = 3; // not used any more (since t is accurate up to a second; tau is much more accurate)

  gSystem->RedirectOutput("consistency_cuts.txt","w");
  printf("\nR3 diagnostic cuts\n");
  printf("| VPDE - ZDCE | < %f\n",d_vz_cc[0]);
  printf("| VPDW - ZDCW | < %f\n",d_vz_cc[1]);
  printf("| VPDX - ZDCX | < %f\n",d_vz_cc[2]);
  printf("| ZDCE - ZDCW | < %f\n",d_xx_cc[0][1]);
  printf("| VPDE - VPDW | < %f\n",d_xx_cc[0][2]);
  printf("| ZDCE - ZDCX | < %f\n",d_xx_cc[1][1]);
  printf("| VPDE - VPDX | < %f\n",d_xx_cc[1][2]);
  printf("| ZDCW - ZDCX | < %f\n",d_xx_cc[2][1]);
  printf("| VPDW - VPDX | < %f\n",d_xx_cc[2][2]);
  printf("t / tau < %f\n",t_tau);
  gSystem->RedirectOutput(0);

  // fill tree
  Int_t ent = rdat->GetEntries();
  for(Int_t e=0; e<ent; e++)
  {
    rdat->GetEntry(e);
    sums->GetEntry(e);

    // diagnostic pass check
    if( fabs(d_vz[0]) < d_vz_cc[0] &&
        fabs(d_vz[1]) < d_vz_cc[1] &&
        fabs(d_vz[2]) < d_vz_cc[2] &&
        fabs(d_xx[0][1]) < d_xx_cc[0][1] &&
        fabs(d_xx[0][2]) < d_xx_cc[0][2] &&
        fabs(d_xx[1][1]) < d_xx_cc[1][1] &&
        fabs(d_xx[1][2]) < d_xx_cc[1][2] &&
        fabs(d_xx[2][1]) < d_xx_cc[2][1] &&
        fabs(d_xx[2][2]) < d_xx_cc[2][2] &&
        t_sums/tau < t_tau )
    {
      isConsistent=1;
    }
    else isConsistent=0;


    // rdat and sum tree consistency check
    if( (i_rdat != i_sums) ||
        (runnum_rdat != runnum_sums) ||
        (fill_rdat != fill_sums) ||
        (t_rdat != t_sums))
    {
      fprintf(stderr,"ERROR: trees not consistent (either i,runnum,fill,t");
      return;
    }
    else rellum->Fill();
  };

  rellum->Write();


  gSystem->RedirectOutput("consistency_cuts.txt","a");
  printf("%d / %d rellum tree entries pass consistency check\n",
    rellum->GetEntries("isConsistent"),
    rellum->GetEntries());
  gSystem->RedirectOutput(0);
  system("cat consistency_cuts.txt");
  printf("\nrellum tree written to rtree.root\n");
};
