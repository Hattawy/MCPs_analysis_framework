#define analysis_theta_dep_cxx
#include "analysis_theta_dep.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
 
int Find_runNumber(const char* DIRNAME);

void analysis_theta_dep::Loop()
{

   // some setting for the histograms --------------------------
   gROOT->Reset();
   gStyle->SetOptStat("nemr");
   gStyle->SetPalette(55);   
   gStyle->SetLabelSize(0.03,"xyz"); // size of axis value font 
   gStyle->SetTitleSize(0.035,"xyz"); // size of axis title font 
   gStyle->SetTitleFont(22,"xyz"); // font option 
   gStyle->SetLabelFont(22,"xyz");
   gStyle->SetTitleOffset(1.2,"y");
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetPadBottomMargin(0.15); //margins... 
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.1);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetPaperSize(20,24);
   gStyle->SetLabelSize(0.05,"xy");
   gStyle->SetTitleSize(0.06,"xy");
   gStyle->SetStatX(0.9);
   gStyle->SetStatY(0.9);

   // output txt file
   ofstream output_txtfile;
   output_txtfile.open("final_theta_dep_txt_output.txt");

   // output root file 
   TFile *output_rootfile = new  TFile("final_theta_dep_output.root", "Recreate");
          TTree *tree = new TTree("TT","the root file");

   // define the ntuple root varaibles 
   int RunNumber;
   int HV_val;
   int B_val;
   double V_LED;
   double Theta;
   double Amplitude;
   double Gain;
   double Time_res;
   double Pos_res;
 
   tree->Branch("RunNumber",&RunNumber, "RunNumber/I");
   tree->Branch("HV_val",&HV_val, "HV_val/I");
   tree->Branch("B_val",&B_val, "B_val/I");
   tree->Branch("V_LED",&V_LED, "V_LED/D");
   tree->Branch("Theta",&Theta, "Theta/D");
   tree->Branch("Amplitude",&Amplitude, "Amplitude/D");
   tree->Branch("Gain",&Gain, "Gain/D");
   tree->Branch("Time_res",&Time_res, "Time_res/D");
   tree->Branch("Pos_res",&Pos_res, "Pos_res/D");
 
   // reading the runs sheet for the tube69-HV-diff-B-2nd-set
   const int N_runs = 34;   // total rumber of runs  
   const int SS = 2;     // two sets of fixed B field 
   const int NN = 17;    // 17 runs at each B field value 
   const int N_val = 5;   // run number, B field, HV, V_LED, theta
   double reading_data[N_runs][N_val];
   ifstream inputfile;
   inputfile.open("runs_sheet_theta_dep.txt");
   for(int i=0; i<N_runs; i++)
     for(int j=0; j<N_val; j++)
     inputfile>>reading_data[i][j];
   inputfile.close();

   // define the histograms here ----------------------
   TH1F *h_Gains[SS][NN]    ; 
   TH1F *h_amp[SS][NN]      ; 
   TH1F *h_timing[SS][NN]   ;  
   TH1F *h_pos_res[SS][NN]  ; 
   TH1F *h_rise_time[SS][NN]; 
   TH1F *h_fall_time[SS][NN];
   bool run_flag[SS][NN]    ; 

   for(int i=0; i<SS; i++){
   for(int j=0; j<NN; j++){
      h_Gains[i][j]     = new TH1F(Form("h_Gains[%d][%d]",i,j),"", 250, 1000,60000000);   
      h_amp[i][j]       = new TH1F(Form("h_amp[%d][%d]",i,j),"", 250, 0, 80);      
      h_timing[i][j]    = new TH1F(Form("h_timing[%d][%d]",i,j),"", 250, -12000,-8000);    
      h_pos_res[i][j]   = new TH1F(Form("h_pos_res[%d][%d]",i,j),"", 250, -100,100);      
      h_rise_time[i][j] = new TH1F(Form("h_rise_time[%d][%d]",i,j),"", 250, 145,165);      
      h_fall_time[i][j] = new TH1F(Form("h_fall_time[%d][%d]",i,j),"", 250, 145,165); 
      run_flag[i][j]  = false;     
    }}

 
 
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   // Here we do loop over the events to show them ------
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
 
      const char *input_file_name;
      input_file_name = fChain->GetCurrentFile()->GetName();
      int runNum = Find_runNumber( input_file_name); 

      // 0.25 T runs 
      if ( runNum < 4000){

          int var_rnum = runNum-2200;
          if (run_flag[0][var_rnum] == false ){
             cout<< var_rnum+2200<<endl;
             run_flag[0][var_rnum] = true;
           }  
       
          if( fwav_time[1]>0.0  && fwav_time[2]>0.0 && fwav_amp[1] >1.0  && fwav_amp[2] >1.0) {
       
               h_Gains[0][var_rnum]  -> Fill( fwav_gain[1]+fwav_gain[2]);
               h_amp[0][var_rnum]    -> Fill( (fwav_amp[1]+fwav_amp[2])/2.0);
               h_timing[0][var_rnum] -> Fill( (fwav_time[1]+fwav_time[2])/2.0 - fwav_time[0]);
               h_pos_res[0][var_rnum]-> Fill( ((fwav_time[1]-fwav_time[2])/2.0)*0.2);  //0.2 is the time bin size 
             }
        }

      // 0.5 T runs 
      if ( runNum >= 4000){

          int var_rnum = runNum-4000;
          if (run_flag[1][var_rnum] == false ){
             cout<< var_rnum+4000<<endl;
             run_flag[1][var_rnum] = true;
           }  
       
          if( fwav_time[1]>0.0  && fwav_time[2]>0.0 && fwav_amp[1] >1.0  && fwav_amp[2] >1.0) {
       
               h_Gains[1][var_rnum]  -> Fill( fwav_gain[1]+fwav_gain[2]);
               h_amp[1][var_rnum]    -> Fill( (fwav_amp[1]+fwav_amp[2])/2.0);
               h_timing[1][var_rnum] -> Fill( (fwav_time[1]+fwav_time[2])/2.0 - fwav_time[0]);
               h_pos_res[1][var_rnum]-> Fill( ((fwav_time[1]-fwav_time[2])/2.0)*0.2);  //0.2 is the time bin size 
             }
        }

   } // end the loop over the events 

   gSystem->mkdir("../main_plots");

   int   theta_val[SS][NN]    ; 
   float amp_val[SS][NN]      ; 
   float gain_val[SS][NN]     ; 
   float time_res_val[SS][NN] ; 
   float pos_res_val[SS][NN]  ; 
   float point_Err[SS][NN]          ; 

   for(int i=0; i<SS; i++){
       for(int j=0; j<NN; j++){
         amp_val[i][j]      =  0.0;
         gain_val[i][j]     =  0.0;
         time_res_val[i][j] =  0.0;
         pos_res_val[i][j]  =  0.0;
         point_Err[i][j]    =  0.0;
       }
    }
   
  output_txtfile<<"RunNum   B         HV       V_LED   Theta   Gain      Amp      time_res      pos_res"<<endl;
   
   // write the results to root and tx files  
   for(int i=0; i<SS; i++){
       for(int j=0; j<NN; j++){

            int kk = 0;
            if (i == 0) kk = j;
            else if ( i == 1) kk = j+17;
 
            RunNumber = reading_data[kk][0];
            HV_val    = reading_data[kk][2];
            B_val     = reading_data[kk][1];
            V_LED     = reading_data[kk][3];
            Theta     = reading_data[kk][4];

            Amplitude = 0.0;
            Gain      = 0.0;
            Time_res  = 0.0;
            Pos_res   = 0.0;

            theta_val[i][j]    = Theta ;
            
            if( (h_amp[i][j]->GetEntries()) > 200.0){ 
            
               h_amp[i][j]->Fit("gaus");
               TF1 *f1 = h_amp[i][j]->GetFunction("gaus");
               Amplitude  = f1->GetParameter(1);
               
               h_Gains[i][j]->Fit("gaus");
               TF1 *f2 = h_Gains[i][j]->GetFunction("gaus");
               Gain  = f2->GetParameter(1);
             
               h_timing[i][j]->Fit("gaus");
               TF1 *f3 = h_timing[i][j]->GetFunction("gaus");
               Time_res  = f3->GetParameter(2);
             
               h_pos_res[i][j]->Fit("gaus");
               TF1 *f4 = h_pos_res[i][j]->GetFunction("gaus");
               Pos_res  = f4->GetParameter(2);
             
               // variables for histograms
            
                  amp_val[i][j]      = Amplitude ;    
                  gain_val[i][j]     = Gain      ;
                  time_res_val[i][j] = Time_res  ;
                  pos_res_val[i][j]  = Pos_res   ; 
             } 
                                        
            // write the txt output file         
            output_txtfile<<RunNumber<<"     "<<B_val<<"     "<<
                            HV_val<<"     "<<V_LED<<"     "<<
                            Theta<<"     "<<Gain<<"     "<<
                            Amplitude<<"     "<<Time_res<<"     "<<
                            Pos_res<<endl;
           
            tree->Fill();
           
          }
        }


 tree->AutoSave();
 output_rootfile->Write();
 output_txtfile.close();

 // plot the amplitude, gain, time resolution as a function of B field at fixed HV

 TMultiGraph *mg_amp      = new TMultiGraph();
 TMultiGraph *mg_gain     = new TMultiGraph();
 TMultiGraph *mg_time_res = new TMultiGraph();
 TMultiGraph *mg_pos_res  = new TMultiGraph();

 TLegend* leg_B = new TLegend(0.75,0.75,0.9,0.9);
          leg_B->SetNColumns(1);          

 
 for( int i=0; i<SS; i++){

    float F_theta[NN]    ;
    float F_amp[NN]      ;
    float F_gain[NN]     ;
    float F_time_res[NN] ;
    float F_pos_res[NN]  ;
    float F_Err[NN]      ; 

    for(int j=0; j<NN; j++){
      
        F_theta[j]   =  theta_val[i][j]   ; 
        F_amp[j]     =  amp_val[i][j]     ;                          
        F_gain[j]    =  gain_val[i][j]    ;
        F_time_res[j]=  time_res_val[i][j];
        F_pos_res[j] =  pos_res_val[i][j] ;
        F_Err[j]     =  point_Err[i][j]   ;
     }                    
  
    TGraphErrors *tg_amp      = new TGraphErrors(NN, F_theta, F_amp, F_Err, F_Err);
    TGraphErrors *tg_gain     = new TGraphErrors(NN, F_theta, F_gain, F_Err, F_Err);
    TGraphErrors *tg_time_res = new TGraphErrors(NN, F_theta, F_time_res, F_Err, F_Err);
    TGraphErrors *tg_pos_res  = new TGraphErrors(NN, F_theta, F_pos_res, F_Err, F_Err);
    
    int kk;
    if(i == 0) kk =1;
    if(i == 1) kk =2;
    tg_amp      ->SetMarkerColor(kk);
    tg_gain     ->SetMarkerColor(kk);
    tg_time_res ->SetMarkerColor(kk);
    tg_pos_res  ->SetMarkerColor(kk);

    tg_amp      ->SetLineColor(kk);
    tg_gain     ->SetLineColor(kk);
    tg_time_res ->SetLineColor(kk);
    tg_pos_res  ->SetLineColor(kk);

    tg_amp      ->SetLineWidth(2);
    tg_gain     ->SetLineWidth(2);
    tg_time_res ->SetLineWidth(2);
    tg_pos_res  ->SetLineWidth(2);

    tg_amp      ->SetMarkerStyle(i+20);
    tg_gain     ->SetMarkerStyle(i+20);
    tg_time_res ->SetMarkerStyle(i+20);
    tg_pos_res  ->SetMarkerStyle(i+20);

    mg_amp     ->Add(tg_amp )     ; 
    mg_gain    ->Add(tg_gain )    ; 
    mg_time_res->Add(tg_time_res ); 
    mg_pos_res ->Add(tg_pos_res ) ;

    if(i==0) leg_B->AddEntry(tg_amp," 0.25 Tesla","PL"); 
    if(i==1) leg_B->AddEntry(tg_amp," 0.50 Tesla","PL"); 
  }

  // print the output histograms here ----------------------------------
   TCanvas *c1=new TCanvas("c1","c1",750,500);
   c1->SetBottomMargin(0.18);
   c1->SetLeftMargin(0.18);
   c1->cd();


   //c1->SetLogy();
   mg_amp-> Draw("APL");
   mg_amp-> GetXaxis()->SetTitle("#theta [deg.]");
   mg_amp-> GetYaxis()->SetTitle("Signal Amplitude [mV]");
   mg_amp-> GetYaxis()->CenterTitle(true);
   mg_amp-> GetXaxis()->CenterTitle(true);
   gPad->Modified();
   //mg_amp->GetYaxis()->SetRangeUser(1.0, 100.0);
   mg_amp-> Draw("APL");
   leg_B->Draw();
   c1->Print("../main_plots/amp_theta.png");

   mg_gain-> Draw("APL");
   mg_amp-> GetXaxis()->SetTitle("#theta [deg.]");
   mg_gain-> GetYaxis()->SetTitle("Gain");
   mg_gain-> GetYaxis()->CenterTitle(true);
   mg_gain-> GetXaxis()->CenterTitle(true);
   gPad->Modified();
   //mg_gain->GetYaxis()->SetRangeUser(1.0, 100.0);
   mg_gain-> Draw("APL");
   leg_B->Draw();
   c1->Print("../main_plots/gain_theta.png");

   mg_time_res-> Draw("APL");
   mg_amp-> GetXaxis()->SetTitle("#theta [deg.]");
   mg_time_res-> GetYaxis()->SetTitle("Time Resolution [ps]");
   mg_time_res-> GetYaxis()->CenterTitle(true);
   mg_time_res-> GetXaxis()->CenterTitle(true);
   gPad->Modified();
   //mg_time_res->GetYaxis()->SetRangeUser(1.0, 100.0);
   mg_time_res-> Draw("APL");
   leg_B->Draw();
   c1->Print("../main_plots/time_res_theta.png");

   mg_pos_res-> Draw("APL");
   mg_amp-> GetXaxis()->SetTitle("#theta [deg.]");
   mg_pos_res-> GetYaxis()->SetTitle("Position Resolution [mm]");
   mg_pos_res-> GetYaxis()->CenterTitle(true);
   mg_pos_res-> GetXaxis()->CenterTitle(true);
   gPad->Modified();
   //mg_pos_res->GetYaxis()->SetRangeUser(1.0, 100.0);
   mg_pos_res-> Draw("APL");
   leg_B->Draw();
   c1->Print("../main_plots/pos_res_theta.png");



} // end of the main loop ----------------------------------


 int Find_runNumber(const char* DIRNAME)
 {
    int RunNumber= 0;
    int fullNameLength = strlen(DIRNAME);
    string StringRunNumber;
    int  FNP = fullNameLength - 9;   // first number position
 
    for (int ii = FNP; ii<FNP+4; ii++)  StringRunNumber += DIRNAME[ii];
 
    std::istringstream iss(StringRunNumber);
    iss >> RunNumber;
   return RunNumber;
  } 


