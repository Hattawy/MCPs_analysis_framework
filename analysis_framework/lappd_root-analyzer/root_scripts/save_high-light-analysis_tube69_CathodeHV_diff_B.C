#define analysis_tube69_CathodeHV_diff_B_cxx
#include "analysis_tube69_CathodeHV_diff_B.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
 
int Find_runNumber(const char* DIRNAME);

void analysis_tube69_CathodeHV_diff_B::Loop()
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
   output_txtfile.open("final_analysis_tube69_CathodeHV_diff_B_txt_output.txt");

   // output root file 
   TFile *output_rootfile = new  TFile("final_analysis_tube69_CathodeHV_diff_B_output.root", "Recreate");
          TTree *tree = new TTree("TT","the root file");

   // define the ntuple root varaibles 
   int RunNumber;
   int HV_val;
   double B_val;
   double V_LED;
   double Theta;
   double Amplitude;
   double Gain;
   double Time_res;
   double Pos_res;
 
   tree->Branch("RunNumber",&RunNumber, "RunNumber/I");
   tree->Branch("HV_val",&HV_val, "HV_val/I");
   tree->Branch("B_val",&B_val, "B_val/D");
   tree->Branch("V_LED",&V_LED, "V_LED/D");
   tree->Branch("Theta",&Theta, "Theta/D");
   tree->Branch("Amplitude",&Amplitude, "Amplitude/D");
   tree->Branch("Gain",&Gain, "Gain/D");
   tree->Branch("Time_res",&Time_res, "Time_res/D");
   tree->Branch("Pos_res",&Pos_res, "Pos_res/D");
 
   // reading the runs sheet for the tube69-HV-diff-B-2nd-set
   const int N_runs = 48;
   const int N_val = 5;   // run number, B field, HV, V_LED, theta
   double reading_data[N_runs][N_val];
   ifstream inputfile;
   inputfile.open("runs_sheet_tube69-CathodeHV-diff-B.txt");
   for(int i=0; i<N_runs; i++)
     for(int j=0; j<N_val; j++)
     inputfile>>reading_data[i][j];
   inputfile.close();


   // define the histograms here ----------------------
   TH1F *h_Gains[100]    ; 
   TH1F *h_amp[100]      ; 
   TH1F *h_timing[100]   ;  
   TH1F *h_pos_res[100]  ; 
   TH1F *h_rise_time[100]; 
   TH1F *h_fall_time[100];
   bool run_flag[100]    ; 
 

   for(int i=0; i<100; i++){
      h_Gains[i]     = new TH1F(Form("h_Gains[%d]",i+3000),"", 250, 100, 300000000);   
      h_amp[i]       = new TH1F(Form("h_amp[%d]",i+3000),"", 250, 0, 500);      
      h_timing[i]    = new TH1F(Form("h_timing[%d]",i+3000),"", 250, 1000,40000);    
      h_pos_res[i]   = new TH1F(Form("h_pos_res[%d]",i+3000),"", 250, -100,100);      
      h_rise_time[i] = new TH1F(Form("h_rise_time[%d]",i+3000),"", 250, 145,165);      
      h_fall_time[i] = new TH1F(Form("h_fall_time[%d]",i+3000),"", 250, 145,165); 
      run_flag[i]  = false;     
    }
 
 
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

      // arrange the runs in an array
      int var_rnum;
      if(9299 <runNum && runNum< 9316) 
        var_rnum = runNum-9300;
      else if (9399 <runNum && runNum< 9416)
         var_rnum = runNum-9400 + 16;
      else var_rnum = runNum-9500 + 32;


      if (run_flag[var_rnum] == false ){
         cout<< var_rnum<<endl;
         run_flag[var_rnum] = true;
       }  
 
      if( fwav_time[1]>0.0  && fwav_time[2]>0.0 && fwav_amp[1] >1.0  && fwav_amp[2] >1.0)
        {

           h_Gains[var_rnum]  -> Fill( fwav_gain[1]+fwav_gain[2]);
           h_amp[var_rnum]    -> Fill( (fwav_amp[1]+fwav_amp[2])/2.0);
           h_timing[var_rnum] -> Fill( (fwav_time[1] + fwav_time[2])/2.0 - fwav_time[0]);
           h_pos_res[var_rnum]-> Fill( ((fwav_time[1] - fwav_time[2])/2.0) * 0.2);   // 0.2 is the time bin size
         }

   } // end the loop over the events 

   gSystem->mkdir("../main_plots");

   output_txtfile<<"RunNum   B         HV       V_LED   Theta   Gain      Amp      time_res      pos_res"<<endl;
   
// for(int i=0; i<N_runs; i++){
//    for(int j=0; j<17; j++){
//       if(reading_data[i][0] == j+9000){

//          RunNumber = reading_data[i][0];
//          HV_val    = reading_data[i][2];
//          B_val     = reading_data[i][1];
//          V_LED     = reading_data[i][3];
//          Theta     = reading_data[i][4];
//          
//          Amplitude = 0.0;
//          Gain      = 0.0;
//          Time_res  = 0.0;
//          Pos_res   = 0.0;

//          if( (h_amp[j]->GetEntries()) > 200.0) {
//             h_amp[j]->Fit("gaus");
//             TF1 *f1 = h_amp[j]->GetFunction("gaus");
//             Amplitude  = f1->GetParameter(1);
//             
//             h_Gains[j]->Fit("gaus");
//             TF1 *f2 = h_Gains[j]->GetFunction("gaus");
//             Gain  = f2->GetParameter(1);
//           
//             h_timing[j]->Fit("gaus");
//             TF1 *f3 = h_timing[j]->GetFunction("gaus");
//             Time_res  = f3->GetParameter(2);
//           
//             h_pos_res[j]->Fit("gaus");
//             TF1 *f4 = h_pos_res[j]->GetFunction("gaus");
//             Pos_res  = f4->GetParameter(2);
//           }


//          output_txtfile<<RunNumber<<"     "<<B_val<<"     "
//                        <<HV_val<<"     "<<   V_LED<<"     "
//                        <<Theta<<"     "<< Gain<<"     "
//                        <<Amplitude<<"     "<<Time_res<<"     "<<Pos_res<<endl;
//         
//          tree->Fill();
//        }
//      }
//   }  

 tree->AutoSave();
 output_rootfile->Write();
 output_txtfile.close();


   // plot the amplitude, gain, time resolution as a function of B field at fixed HV
  
   TMultiGraph *mg_amp      = new TMultiGraph();
   TMultiGraph *mg_gain     = new TMultiGraph();
   TMultiGraph *mg_time_res = new TMultiGraph();
   TMultiGraph *mg_pos_res  = new TMultiGraph();
  
   //1. list the HV sets
   const int N_HV_sets = 16;
   int HV[N_HV_sets] = { 0, 100, 200, 300, 400, 500, 600, 700, 
                         800, 900, 1000, 1100, 1200, 1300, 1400, 1500};
  
   //2. number of B field points
    const int N_B_field = 3;
  
   //3. find the B field, runNum, amp, gain, time_res, pos_res values at each fixed HV 
  
   vector<vector<int> > HV_RN_val         ; HV_RN_val.resize(N_B_field)       ;
   vector<vector<float> > HV_B_val        ; HV_B_val.resize(N_B_field)        ;
   vector<vector<float> > HV_amp_val      ; HV_amp_val.resize(N_B_field)      ;
   vector<vector<float> > HV_gain_val     ; HV_gain_val.resize(N_B_field)     ;
   vector<vector<float> > HV_time_res_val ; HV_time_res_val.resize(N_B_field) ;
   vector<vector<float> > HV_pos_res_val  ; HV_pos_res_val.resize(N_B_field)  ;
   vector<vector<float> > point_Err       ; point_Err.resize(N_B_field)       ;
   
   for( int i=0; i<N_B_field; i++){ 
      
      HV_RN_val[i].resize(N_HV_sets);
      HV_B_val[i].resize(N_HV_sets);
      HV_amp_val[i].resize(N_HV_sets);
      HV_gain_val[i].resize(N_HV_sets);
      HV_time_res_val[i].resize(N_HV_sets);
      HV_pos_res_val[i].resize(N_HV_sets);
      point_Err[i].resize(N_HV_sets);
    }
  
   for( int i=0; i<N_B_field; i++){ 
     for(int j=0; j<N_HV_sets; j++){ 
        HV_RN_val[i][j]        = 0 ;
        HV_B_val[i][j]         = 0.0; 
        HV_amp_val[i][j]       = 0.0;
        HV_gain_val[i][j]      = 0.0;
        HV_time_res_val[i][j]  = 0.0;
        HV_pos_res_val[i][j]   = 0.0;
        point_Err[i][j]        = 0.0;
      }
    }
  
   
   for( int i=0; i<N_B_field; i++){
      for(int j=0; j<N_HV_sets; j++){
           
           int nnum;
           if ( i== 0) { nnum = j; }
           else if(i== 1) { nnum = j+16; }
           else if(i== 2) { nnum = j+32; }        

           HV_RN_val[i][j]       = reading_data[nnum][0]; 
           HV_B_val[i][j]        = reading_data[nnum][1];
  
           if( (h_amp[nnum]->GetEntries()) > 100.0) {
   
              cout<<" nnum = "<< nnum<<endl;
              h_amp[nnum]->Fit("gaus");
              TF1 *f1 = h_amp[nnum]->GetFunction("gaus");
              HV_amp_val[i][j]  = f1->GetParameter(1);
              
              h_Gains[nnum]->Fit("gaus");
              TF1 *f2 = h_Gains[nnum]->GetFunction("gaus");
              HV_gain_val[i][j]  = f2->GetParameter(1);
            
              h_timing[nnum]->Fit("gaus");
              TF1 *f3 = h_timing[nnum]->GetFunction("gaus");
              HV_time_res_val[i][j]  = f3->GetParameter(2);
            
              h_pos_res[nnum]->Fit("gaus");
              TF1 *f4 = h_pos_res[nnum]->GetFunction("gaus");
              HV_pos_res_val[i][j]  = f4->GetParameter(2);
            
              point_Err[i][j]       = 0.0 ;
  
              cout<<"first loop    "<<  HV_RN_val[i][j] <<"      "<< HV_amp_val[i][j]<<endl;
           }

           else {
               HV_amp_val[i][j]  = 0.0           ;
               HV_gain_val[i][j] = 0.0           ;
               HV_time_res_val[i][j]  = -9999.0  ;
               HV_pos_res_val[i][j]  =  -9999.0  ;
           }
  
       }
    }  
  
   //4. creating the histograms
  
   TLegend* leg_HV = new TLegend(0.75,0.75,0.9,0.9);
            leg_HV->SetNColumns(1);          
   float B_vall[3] = { 1.0, 0.0, 0.5 };

   for( int i=0; i<3; i++){
  
      const int   F_N  = 16;;
      float F_Cathode_HV[F_N];
      float F_B[F_N]        ;
      float F_amp[F_N]      ;
      float F_gain[F_N]     ;
      float F_time_res[F_N] ;
      float F_pos_res[F_N]  ;
      float F_Err[F_N]      ; 
  
      for(int j=0; j<F_N; j++){
          
          F_Cathode_HV[j] = HV[j]               ;
          F_B[j]       =  HV_B_val[i][j]        ; 
          F_amp[j]     =  HV_amp_val[i][j]      ;                          
          F_gain[j]    =  HV_gain_val[i][j]     ;
          F_time_res[j]=  HV_time_res_val[i][j] ;
          F_pos_res[j] =  HV_pos_res_val[i][j]  ;
          F_Err[j]     =  point_Err[i][j]       ;
  
          cout<< F_Cathode_HV[j] <<"       "<<F_amp[j]<<endl;
       }                    
    
      TGraphErrors *tg_amp      = new TGraphErrors(F_N, F_Cathode_HV, F_amp, F_Err, F_Err);
      TGraphErrors *tg_gain     = new TGraphErrors(F_N, F_Cathode_HV, F_gain, F_Err, F_Err);
      TGraphErrors *tg_time_res = new TGraphErrors(F_N, F_Cathode_HV, F_time_res, F_Err, F_Err);
      TGraphErrors *tg_pos_res  = new TGraphErrors(F_N, F_Cathode_HV, F_pos_res, F_Err, F_Err);
  
      int kk = i;
      if ( i == 0) kk=4;
      tg_amp      ->SetLineColor(kk);
      tg_gain     ->SetLineColor(kk);
      tg_time_res ->SetLineColor(kk);
      tg_pos_res  ->SetLineColor(kk);

      tg_amp      ->SetMarkerColor(kk);
      tg_gain     ->SetMarkerColor(kk);
      tg_time_res ->SetMarkerColor(kk);
      tg_pos_res  ->SetMarkerColor(kk);

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
  
      leg_HV->AddEntry(tg_amp,Form(" %.2f T",B_vall[i]),"PL"); 
  }
  
    // print the output histograms here ----------------------------------
     TCanvas *c1=new TCanvas("c1","c1",750,500);
     c1->SetBottomMargin(0.18);
     c1->SetLeftMargin(0.18);
     c1->cd();
  
  
     //c1->SetLogy();
     mg_amp-> Draw("APL");
     mg_amp-> GetXaxis()->SetTitle("#Delta V (Cathode-MCP)[V]");
     //mg_amp-> GetXaxis()->SetTitle("Magnetic Field Strength [Tesla]");
     mg_amp-> GetYaxis()->SetTitle("Signal Amplitude [mV]");
     mg_amp-> GetYaxis()->CenterTitle(true);
     mg_amp-> GetXaxis()->CenterTitle(true);
     gPad->Modified();
     //mg_amp->GetYaxis()->SetRangeUser(1.0, 100.0);
     mg_amp-> Draw("APL");
     leg_HV->Draw();
     c1->Print("../main_plots/Cathode_amp_B.png");
  
     mg_gain-> Draw("APL");
     mg_gain-> GetXaxis()->SetTitle("#Delta V (Cathode-MCP)[V]");
     mg_gain-> GetYaxis()->SetTitle("Gain");
     mg_gain-> GetYaxis()->CenterTitle(true);
     mg_gain-> GetXaxis()->CenterTitle(true);
     gPad->Modified();
     //mg_gain->GetYaxis()->SetRangeUser(1.0, 100.0);
     mg_gain-> Draw("APL");
     leg_HV->Draw();
     c1->Print("../main_plots/Cathode_gain_B.png");
  
     mg_time_res-> Draw("APL");
     mg_time_res-> GetXaxis()->SetTitle("#Delta V (Cathode-MCP)[V]");
     mg_time_res-> GetYaxis()->SetTitle("Time Resolution [ps]");
     mg_time_res-> GetYaxis()->CenterTitle(true);
     mg_time_res-> GetXaxis()->CenterTitle(true);
     gPad->Modified();
     //mg_time_res->GetYaxis()->SetRangeUser(1.0, 100.0);
     mg_time_res-> Draw("APL");
     leg_HV->Draw();
     c1->Print("../main_plots/Cathode_time_res__B.png");
  
     mg_pos_res-> Draw("APL");
     mg_pos_res-> GetXaxis()->SetTitle("#Delta V (Cathode-MCP)[V]");
     mg_pos_res-> GetYaxis()->SetTitle("Position Resolution [mm]");
     mg_pos_res-> GetYaxis()->CenterTitle(true);
     mg_pos_res-> GetXaxis()->CenterTitle(true);
     gPad->Modified();
     //mg_pos_res->GetYaxis()->SetRangeUser(1.0, 100.0);
     mg_pos_res-> Draw("APL");
     leg_HV->Draw();
     c1->Print("../main_plots/Cathode_pos_res__B.png");
  

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


