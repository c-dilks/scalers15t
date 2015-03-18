// builds sums.root, which is scaler counts summed
// for each run number
//
// --determines spin pattern and adds 3 branches to the tree which
//   state the spin patterns collided for each run
//
// --zeroAborts will filter out counts in abort gaps and
//   bunches listed as empty by cdev; this causes slope for
//   total bXings / bXing rate vs. run time to be less than unity
//   and therefore should be left as zeroAborts=false

void sumTree(const char * filename="counts.root")
{
  Bool_t zeroAborts = false;

  TFile * infile = new TFile(filename,"READ");
  TTree * str = (TTree*) infile->Get("sca");

  system("touch pattern_log.txt; rm pattern_log.txt");

  // read counts.root tree
  Int_t i,runnum,fi,fill,bx,blue,yell;
  Double_t t,freq,bbce,bbcw,bbcx,zdce,zdcw,zdcx,vpde,vpdw,vpdx,tot_bx;
  str->SetBranchAddress("i",&i);
  str->SetBranchAddress("runnum",&runnum);
  str->SetBranchAddress("fi",&fi);
  str->SetBranchAddress("fill",&fill);
  str->SetBranchAddress("t",&t);
  str->SetBranchAddress("freq",&freq);
  str->SetBranchAddress("bx",&bx);
  str->SetBranchAddress("bbce",&bbce);
  str->SetBranchAddress("bbcw",&bbcw);
  str->SetBranchAddress("bbcx",&bbcx);
  str->SetBranchAddress("zdce",&zdce);
  str->SetBranchAddress("zdcw",&zdcw);
  str->SetBranchAddress("zdcx",&zdcx);
  str->SetBranchAddress("vpde",&vpde);
  str->SetBranchAddress("vpdw",&vpdw);
  str->SetBranchAddress("vpdx",&vpdx);
  str->SetBranchAddress("tot_bx",&tot_bx);
  str->SetBranchAddress("blue",&blue);
  str->SetBranchAddress("yell",&yell);


  // get number of runs in each fill (used for averaging things over runs in a fill)
  Int_t fi_max_tmp = str->GetMaximum("fi");
  Int_t i_max_tmp = str->GetMaximum("i");
  const Int_t fi_max = fi_max_tmp;
  const Int_t i_max = i_max_tmp;
  Int_t num_runs[fi_max];
  for(Int_t n=0; n<fi_max; n++) num_runs[n]=0;
  for(Int_t n=0; n<str->GetEntries(); n++)
  {
    str->GetEntry(n);
    if(bx==0) num_runs[fi-1]+=1;
  };
  //for(Int_t n=0; n<fi_max; n++) printf("%d - %d\n",n,num_runs[n]);



  // spin pattern definitions
  /* (overall pattern no.) = 6 digits:
   *  (bsp)(bsp)(bsp)(ysp)(ysp)(ysp)
   *  where bsp = blue subpattern and
   *        ysp = yell subpattern
   *  and "000" indicates pattern recognition failure
   *
   * 2015 spin patterns
   * ------------------
   *  4 8-bXing subpatterns:
   *   [1] = + - + - - + - +
   *   [2] = - + - + + - + -
   *   [3] = + + - - + + - -
   *   [4] = - - + + - - + + 
   *  combine 3 subpatterns to form 4 24-bXing patterns
   *   P1 = [1] [2] [2]
   *   P2 = [2] [1] [1]
   *   P3 = [3] [3] [3]
   *   P4 = [4] [4] [4]
   *
   *   note that below, "pattern_no" is the pattern seen at STAR
   *   bXings, where blue bunch 0 collides with yell bunch 80 and
   *   there is an odd number of spin flips from the proton source;
   *   "pattern_no_cdev" is the pattern number seen at CDEV
   *   (i.e. without yellow cogging and spin flip)
  */

  Int_t subpattern[4][8]; // [subpattern no.] [bXing]

  subpattern[0][0] =  1; // subpattern 1
  subpattern[0][1] = -1;
  subpattern[0][2] =  1;
  subpattern[0][3] = -1;
  subpattern[0][4] = -1;
  subpattern[0][5] =  1;
  subpattern[0][6] = -1;
  subpattern[0][7] =  1;

  subpattern[2][0] =  1; // subpattern 3
  subpattern[2][1] =  1;
  subpattern[2][2] = -1;
  subpattern[2][3] = -1;
  subpattern[2][4] =  1;
  subpattern[2][5] =  1;
  subpattern[2][6] = -1;
  subpattern[2][7] = -1;

  for(Int_t b=0; b<8; b++)
  {
    subpattern[1][b] = -1 * subpattern[0][b]; // subpattern 2
    subpattern[3][b] = -1 * subpattern[2][b]; // subpattern 4
  };


  // fill sum tree
  Int_t num_runs_cur;
  Double_t bbce_tot,bbcw_tot,bbcx_tot,zdce_tot,zdcw_tot;
  Double_t zdcx_tot,vpde_tot,vpdw_tot,vpdx_tot,tot_bx_tot;
  Double_t tau; // total bXings / bXing rate
  Int_t pattern_no,pattern_no_cdev;
  Int_t cnt_blue_pos,cnt_blue_neg;
  Int_t cnt_yell_pos,cnt_yell_neg;
  Int_t cnt_0,cnt_1,cnt_2,cnt_3;
  cnt_blue_pos=cnt_blue_neg=cnt_yell_pos=cnt_yell_neg=0;
  cnt_0=cnt_1=cnt_2=cnt_3=0;
  TFile * outfile = new TFile("sums.root","RECREATE");
  TTree * sum = new TTree("sum","sum");
  sum->Branch("i",&i,"i/I");
  sum->Branch("runnum",&runnum,"runnum/I");
  sum->Branch("fi",&fi,"fi/I");
  sum->Branch("fill",&fill,"fill/I");
  sum->Branch("t",&t,"t/D");
  sum->Branch("freq",&freq,"freq/D");
  sum->Branch("tau",&tau,"tau/D");
  sum->Branch("num_runs",&num_runs_cur,"num_runs/I"); // total no. runs in a fill
  sum->Branch("bbce",&bbce_tot,"bbce/D");
  sum->Branch("bbcw",&bbcw_tot,"bbcw/D");
  sum->Branch("bbcx",&bbcx_tot,"bbcx/D");
  sum->Branch("zdce",&zdce_tot,"zdce/D");
  sum->Branch("zdcw",&zdcw_tot,"zdcw/D");
  sum->Branch("zdcx",&zdcx_tot,"zdcx/D");
  sum->Branch("vpde",&vpde_tot,"vpde/D");
  sum->Branch("vpdw",&vpdw_tot,"vpdw/D");
  sum->Branch("vpdx",&vpdx_tot,"vpdx/D");
  sum->Branch("tot_bx",&tot_bx_tot,"tot_bx/D");
  sum->Branch("pattern",&pattern_no,"pattern/I"); // overall pattern no at STAR bXings
  sum->Branch("patternCDEV",&pattern_no_cdev,"patternCDEV/I"); // overall pattern no from CDEV
  sum->Branch("cnt_blue_pos",&cnt_blue_pos,"cnt_blue_pos/I"); // no. blue=1
  sum->Branch("cnt_blue_neg",&cnt_blue_neg,"cnt_blue_neg/I"); // no. blue=-1
  sum->Branch("cnt_yell_pos",&cnt_yell_pos,"cnt_yell_pos/I"); // no. yell=1
  sum->Branch("cnt_yell_neg",&cnt_yell_neg,"cnt_yell_neg/I"); // no. yell=-1
  sum->Branch("cnt_0",&cnt_0,"cnt_0/I"); // no. - - collisions
  sum->Branch("cnt_1",&cnt_1,"cnt_1/I"); // no. - + collisions
  sum->Branch("cnt_2",&cnt_2,"cnt_2/I"); // no. + - collisions
  sum->Branch("cnt_3",&cnt_3,"cnt_3/I"); // no. + + collisions

  // tree loop; this assumes tree from counts.root is ordered by
  // run number and further by bXing number
  Int_t blue_pattern_cnt[4]; // patern check counters
  Int_t yell_pattern_cnt[4];
  Int_t blue_mul,yell_mul; // multipliers to determine pattern no.
  Int_t blue_stat,yell_stat; // error catching
  Int_t blue_subp[15];
  Int_t yell_subp[15];
  Int_t blue_subp_sum,yell_subp_sum;
  Bool_t blue_ok,yell_ok;
  Int_t fill_tmp=0;
  for(Int_t n=0; n<str->GetEntries(); n++)
  {
    str->GetEntry(n);
    num_runs_cur = num_runs[fi-1];

    if(bx%8==0)
    {
      // reset pattern check counters
      for(Int_t p=0; p<4; p++)
      {
        blue_pattern_cnt[p] = 0;
        yell_pattern_cnt[p] = 0;
      };
    };

    for(Int_t p=0; p<4; p++)
    {
      if(blue==subpattern[p][bx%8] || blue==0) blue_pattern_cnt[p]++;
      if(yell==subpattern[p][bx%8] || yell==0) yell_pattern_cnt[p]++;
    };

    if(bx%8==7)
    {
      blue_stat=yell_stat=0;
      for(Int_t p=0; p<4; p++)
      {
        if(blue_pattern_cnt[p]==8)
        {
          blue_subp[bx/8]=p+1;
          blue_stat++;
        }
        if(yell_pattern_cnt[p]==8)
        {
          yell_subp[bx/8]=p+1;
          yell_stat++;
        }
      };
      if(blue_stat!=1) blue_subp[bx/8]=0;
      if(yell_stat!=1) yell_subp[bx/8]=0;
    };

    if(bx==119 && fill!=fill_tmp)
    {
      fill_tmp=fill;
      blue_subp_sum=yell_subp_sum=0;


      // check if pattern is consistent for all bXings
      blue_ok=yell_ok=true;
      for(Int_t c=0; c<4; c++)
      {
        for(Int_t cc=0; cc<3; cc++)
        {
          if(blue_subp[3*c+cc]!=blue_subp[3*(c+1)+cc] && 
             blue_subp[3*c+cc]*blue_subp[3*(c+1)+cc]!=0) blue_ok=false;
          if(yell_subp[3*c+cc]!=yell_subp[3*(c+1)+cc] && 
             yell_subp[3*c+cc]*yell_subp[3*(c+1)+cc]!=0) yell_ok=false;
        };
      };

      // determine 6-digit pattern number
      pattern_no=pattern_no_cdev=0;
      for(Int_t c=0; c<3; c++)
      {
        pattern_no += ((Int_t)blue_ok) * blue_subp[c] * pow(10,5-c);
        pattern_no += ((Int_t)yell_ok) * yell_subp[c] * pow(10,2-c);
        pattern_no_cdev += ((Int_t)blue_ok) * Flip(blue_subp[c]) * pow(10,5-c);
        pattern_no_cdev += ((Int_t)yell_ok) * Flip(yell_subp[c+5]) * pow(10,2-c);
      };

      // print out pattern_log.txt
      gSystem->RedirectOutput("pattern_log.txt");
      printf("fill %d\n",fill);
      for(Int_t z=0; z<15; z++) 
      {
        printf("B%d x Y%d   bbx:%d-%d  ybx:%d-%d\n",blue_subp[z],yell_subp[z],
          z*8,z*8+7,(z*8+80)%120,(z*8+7+80)%120);
        blue_subp_sum += blue_subp[z];
        yell_subp_sum += yell_subp[z];
      };
      printf("sum: %d %d\n",blue_subp_sum,yell_subp_sum);
      printf("%d:  B[%d%d%d] x Y[%d%d%d] <%d> (@ STAR)\n",
        fill,blue_subp[0],blue_subp[1],blue_subp[2],
        yell_subp[0],yell_subp[1],yell_subp[2],pattern_no);
      printf("%d:  B[%d%d%d] x Y[%d%d%d] <%d> (@ CDEV)\n\n",
        fill,Flip(blue_subp[0]),Flip(blue_subp[1]),Flip(blue_subp[2]),
        Flip(yell_subp[5]),Flip(yell_subp[6]),Flip(yell_subp[7]),pattern_no_cdev);
      gSystem->RedirectOutput(0);
    };
    

    // increment counters
    if(blue==-1 && yell==-1)
    {
      cnt_blue_neg++;
      cnt_yell_neg++;
      cnt_0++;
    }
    else if(blue==-1 && yell==1)
    {
      cnt_blue_neg++;
      cnt_yell_pos++;
      cnt_1++;
    }
    else if(blue==1 && yell==-1)
    {
      cnt_blue_pos++;
      cnt_yell_neg++;
      cnt_2++;
    }
    else if(blue==1 && yell==1)
    {
      cnt_blue_pos++;
      cnt_yell_pos++;
      cnt_3++;
    };


    // increment scalers
    if(!zeroAborts || (blue!=0 && yell!=0))
    {
      bbce_tot += bbce;
      bbcw_tot += bbcw;
      bbcx_tot += bbcx;
      zdce_tot += zdce;
      zdcw_tot += zdcw;
      zdcx_tot += zdcx;
      vpde_tot += vpde;
      vpdw_tot += vpdw;
      vpdx_tot += vpdx;
      tot_bx_tot += tot_bx;
    };

    // fill sum tree
    if(bx==119)
    {
      tau = tot_bx_tot/(freq*pow(10,6));
      sum->Fill();
      bbce_tot = 0;
      bbcw_tot = 0;
      bbcx_tot = 0;
      zdce_tot = 0;
      zdcw_tot = 0;
      zdcx_tot = 0;
      vpde_tot = 0;
      vpdw_tot = 0;
      vpdx_tot = 0;
      tot_bx_tot = 0;
      cnt_blue_pos=cnt_blue_neg=cnt_yell_pos=cnt_yell_neg=0;
      cnt_0=cnt_1=cnt_2=cnt_3=0;
    };
  };

  sum->Write();
  //sum->Scan("fi:i:num_runs");
  sum->Print();
}

// Flip : 1|->2, 2|->1, 3|->4, 4|->3
Int_t Flip(Int_t x)
{
  return x+2*(x%2)-1;
};

