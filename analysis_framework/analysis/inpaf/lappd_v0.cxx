//----------APS Analysis 
//----------Author's Name:Jingbo WANG
//----------Copyright:Those valid for ANL
//----------Modified:22/07/2014
////  Run this Code by command :
//		(1) ./lappd inputfile XXX.root parameterfile
//		(2) ./lappd inputfile XXX.root parameterfile UVeventnumber
//		(3) ./lappd inputfile XXX.root parameterfile StartUVeventNumber EndUVeventNumber
//
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <sys/types.h>
#include "math.h"
#include <string.h>
#include <string>
#ifndef __CINT__
#include "TROOT.h"
#include "TFrame.h"
#include "TFile.h" 
#include "TChain.h"
#include <Riostream.h>
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TF1.h"
#include "TPostScript.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "Calibration.h"
#include "TSplineFit.h"
#include "UVevent.h"
using namespace std;
using std::cout;
using std::endl;
#endif


char*	treename="lappd"; 
int		UsersEntries = 100;						// anslyze a user defined number of UVevents
int		StartEntry = 0;							// if(argc==4) Start analysis from this entry
int		EndEntry = 100;							// if(argc==4) End analysis at this entry

//------------------Parameter file format-------------------------------------------// 
int				InputFileFormat = 1;
int				Nacquisition = 10;						// Number of acquisitions
int				NfastFrame = 1000;						// 999 UVevents per acquisition
int				NChannel = 60;								// number of channels

float			WavSamplingRate = 1.e10;				//Sampling rate
float			WavDeltat = 100;
int				WavDimSize = 256;							// 1032 samples per waveform

int				NtrigChannel = 2;
vector<int>		WhichTrig;

int				IfSelectChannel = 0;
int				NChanSelected = 0;
vector<int>		WhichChannel;

int				IfDynamicWindow = 1;
float			FixTrigTwindow_min = 0, FixTrigTwindow_max = 25600;
float			FixMcpTwindow_min = 0, FixMcpTwindow_max = 25600;

int				DoFFT = 0;
float			TrigCutoffFreq = 1.e9;
float			McpCutoffFreq = 5.e8;

float			TrigFraction_CFD = 0.5;
float			McpFraction_CFD = 0.5;
float			TrigDelay_CFD = 800;
float			McpDelay_CFD = 1000;
int				TrigPointsPerSpline = 1;
int				McpPointsPerSpline = 5;

int				DoAmpCut = 1;
float			AmpThreshold = 50;

int				DoTotCut = 1;
float			TrigTotThreshold = 50;
float			McpTotThreshold = 50;
float			TrigMinimumTot = 800;
float			McpMinimumTot = 1500;

int				DoCalibration = 0;
void SetParameters(char *par);

//Slewing correction
int		ProfileMode = 0;
int		BinSizeMode = 0;
vector<double> CalibY[20];							// Amplitude
vector<double> CalibX[20];							// Time
void GetAnaInfo(UVevent* fevent);
void Analyze();

//===========================================================================================================
//	Start main() function
//-----------------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {
	
	cout<<"-------------------------------------------------------------------"<<endl;
	cout<<"-------------------------------------------------------------------"<<endl;

//------------------Set input and output filenames-------------------------------------------// 
	
	TStopwatch clock;								//runtime clock
	clock.Start();									//clock started	  
	char *FileInput = "test.fff";			//default input file
	char *FileOutput = "otp.root";			//default output file
	char *FileParameters = "par.txt";
	switch(argc){
		case 1: {cout<<"Error: File name required as the 2nd statement"<<endl; break;} 
		case 2: {FileInput = argv[1]; break;}
		case 3: {FileInput = argv[1]; FileOutput = argv[2]; break;} 
		case 4: {FileInput = argv[1]; FileOutput = argv[2]; FileParameters = argv[3]; break;}
		case 5: {FileInput = argv[1]; FileOutput = argv[2]; FileParameters = argv[3]; UsersEntries = atoi(argv[4]); break;}
		case 6: {FileInput = argv[1]; FileOutput = argv[2]; FileParameters = argv[3]; StartEntry = atoi(argv[4]); EndEntry = atoi(argv[5]); break;}
   }
	cout<<"Input file:"<<"\t"<< FileInput << endl;
	cout<< "Output file:"<<"\t" << FileOutput << endl;
	cout<<"Parameter file:"<<"\t"<< FileParameters << endl;
	

//--------------------------Set parameters----------------------------------------------------//
	SetParameters(FileParameters);
	vector<int> polarity;
	for(int i=0;i<NChannel;i++) {polarity.push_back(-1);}	//Set all polarity to -1;
	for(int i=0;i<NtrigChannel;i++) {
		int k = WhichTrig.at(i);
		if(k<=NChannel) polarity.at(k) = 1;					//Set trigger polarity to 1;
		else cout<<"Error: Trigger Channel out of range "<<NChannel<<endl;
	} 
	vector<int> selection;
	for(int i=0;i<NChannel;i++) {selection.push_back(0);}	//Set all channel-flag to 0;
	for(int i=0;i<NChanSelected;i++) {
		int k = WhichChannel.at(i);
		if(k<=NChannel) selection.at(k) = 1;					//Set selected channel-flag to 1;
		else cout<<"Error: Selected Channel out of range "<<NChannel<<endl;
	}
   
//------------------Create a root file--------------------------------------------------------//
	TFile f(FileOutput, "RECREATE");				//open a file
	TTree lappd(treename, treename);				//create a tree
	UVevent *event = new UVevent(WavSamplingRate, WavDimSize);				//point to class UVevent
	lappd.Branch("event", &event, 16000, 99);	//Set branches with "split" option
	f.mkdir("wav");							//histograms of waveform
	f.mkdir("cfd");									//constant fraction discrimination

//---------------------Fill a tree------------------------------------------------------------//
	ifstream isin;
	isin.open(FileInput, ios::in);					//open input .txt files
	float WavSample[NChannel][WavDimSize];
	int Nentries = Nacquisition*NfastFrame;			//number of entries
	
	cout<<"Total # of UVevents: "<<Nentries<<endl;

	for (int i=0;i<Nentries;i++) {
		if(argc==5) {if(i>UsersEntries) break;}
		if(argc==6) {if(i<StartEntry) continue; if(i>EndEntry) break;}
		clock.Stop();
		if(i%100 == 0) cout<<"processing "<<i<<"th UVevent..."<<"\t"<<"CpuTime = "<<clock.CpuTime()<<" s"<<endl;
		clock.Continue();		
		event->Setevtno(i);						//Set UVevent number
		event->SetWavDimSize(WavDimSize);
		event->SetWavSampleRate(WavSamplingRate);
		event->Initialize();

//------------------------------------read scope data file---------------------------------------------
		if(InputFileFormat == 0) {
			for(int m=0;m<NChannel;m++) {
				for(int n=0;n<WavDimSize;n++) {				
				isin>>WavSample[m][n];
				WavSample[m][n] *= polarity.at(m);
				}
			}
		}
//------------------------------------read psec data file---------------------------------------------
		if(InputFileFormat == 1) {
			int nsa = 0;					//read waveforms
			float meta = 0;
			for(int n=0;n<WavDimSize;n++) {
				isin>>nsa;
				for(int m=0;m<NChannel/2;m++) {
					isin>>WavSample[m][n];
					WavSample[m][n] *= -polarity.at(m);
				}
				isin>>meta;
				isin>>nsa;
				for(int m=NChannel-1;m>=NChannel/2;m--) {
					isin>>WavSample[m][n];
					WavSample[m][n] *= -polarity.at(m);
				}
				isin>>meta;
			}
		}
//------------------------------------------------------------------------------------------------------		
				
		int fnwav = 0;
		for(int m=0;m<NChannel;m++) {				//loop over channels								
			if(IfSelectChannel==1) {if(selection.at(m)==0) continue;}
			Waveform *waveform = (Waveform*)(event->fwav->ConstructedAt(event->fnwav++));
			waveform->Setevtno(i);
			waveform->Setchno(m);
			waveform->SetSamplingrate(WavSamplingRate);
			waveform->SetDeltat(WavDeltat);
			waveform->SetDimSize(WavDimSize);
			waveform->Initialize(&WavSample[m][0]);
			waveform->SetAmpThreshold(AmpThreshold);
			
			if(polarity.at(m)==1) {					//trigger
				waveform->SetMinimumTot(TrigMinimumTot);
				waveform->SetCutoffFrequency(TrigCutoffFreq);
				waveform->SetPointsPerSpline(TrigPointsPerSpline);
				waveform->SetFraction_CFD(TrigFraction_CFD);
				waveform->SetDelay_CFD(TrigDelay_CFD);
				waveform->SetFitWindow(FixTrigTwindow_min, FixTrigTwindow_max);
			}
			else {									//mcp
				waveform->SetMinimumTot(McpMinimumTot);
				waveform->SetCutoffFrequency(McpCutoffFreq);
				waveform->SetPointsPerSpline(McpPointsPerSpline);
				waveform->SetFraction_CFD(McpFraction_CFD);
				waveform->SetDelay_CFD(McpDelay_CFD);
				waveform->SetFitWindow(FixMcpTwindow_min, FixMcpTwindow_max);
			}
			if(DoFFT==1) waveform->EnableFFT();
			if(IfDynamicWindow == 1) waveform->EnableDynamicWindow();
			waveform->Analyze();					//analyze waveform
			f.cd("wav");								//go to folder "/waveform"
			if(i<10) waveform->hwav->Write();					//write histograms for waveforms
			f.cd("cfd");
			if(i<10) waveform->hcfd->Write();
			
		}
//		GetAnaInfo(event);
		lappd.Fill();								//Fill tree
		
		event->Clear();							//free memory (! to be carefully checked)	
	}

	f.cd();											//home folder
	lappd.Write();									//write tree
//	Analyze();
	
	isin.close();
	cout<<"-------------------------------------------------------------------"<<endl;
	cout<<"-------------------------------------------------------------------"<<endl;
	
	f.Close();

	return 0;
}
//-----------------------------------------------------------------------------------------------------------
//	End main() function
//===========================================================================================================

void SetParameters(char* par) {
	string name;
	int fileopt, k;
	ifstream parin;
	parin.open(par, ios::in);					//open input .txt files
	parin>>name>>InputFileFormat>>name>>Nacquisition>>name>>NfastFrame>>name>>NChannel;
	parin>>name>>WavSamplingRate>>name>>WavDeltat>>name>>WavDimSize;
	parin>>name>>NtrigChannel>>name;
	for(int i=0;i<NtrigChannel;i++) {
		parin>>k;
		WhichTrig.push_back(k);
	}
	parin>>name>>IfSelectChannel>>name>>NChanSelected>>name;
	for(int i=0;i<NChanSelected;i++) {
		parin>>k;
		WhichChannel.push_back(k);
	}
	parin>>name>>IfDynamicWindow>>name>>FixTrigTwindow_min>>FixTrigTwindow_max>>name>>FixMcpTwindow_min>>FixMcpTwindow_max;
	parin>>name>>DoFFT>>name>>TrigCutoffFreq>>name>>McpCutoffFreq;
	parin>>name>>TrigFraction_CFD>>name>>McpFraction_CFD>>name>>TrigDelay_CFD>>name>>McpDelay_CFD>>name>>TrigPointsPerSpline>>name>>McpPointsPerSpline;
	parin>>name>>DoAmpCut>>name>>AmpThreshold;
	parin>>name>>DoTotCut>>name>>TrigTotThreshold>>name>>McpTotThreshold>>name>>TrigMinimumTot>>name>>McpMinimumTot;
	parin>>name>>DoCalibration;
	cout<<"!!!!!!"<<DoCalibration<<endl;
}

void GetAnaInfo(UVevent* fevent) {
	float t[10], q[10];
	t[0] = ((Waveform*)(fevent->fwav->ConstructedAt(0)))->time;
	t[1] = ((Waveform*)(fevent->fwav->ConstructedAt(1)))->time;
	t[2] = ((Waveform*)(fevent->fwav->ConstructedAt(2)))->time;
	t[3] = ((Waveform*)(fevent->fwav->ConstructedAt(3)))->time;
	q[0] = ((Waveform*)(fevent->fwav->ConstructedAt(0)))->qfast;
	q[1] = ((Waveform*)(fevent->fwav->ConstructedAt(1)))->qfast;
	q[2] = ((Waveform*)(fevent->fwav->ConstructedAt(2)))->qfast;
	q[3] = ((Waveform*)(fevent->fwav->ConstructedAt(3)))->qfast;
	if(t[2]==0 || t[3]==0) return;
	CalibY[0].push_back((t[2]+t[3])/2-(t[0]+t[1])/2);
	CalibX[0].push_back(q[0]);
	CalibX[1].push_back(q[1]);
	CalibX[2].push_back(q[2]);
	CalibX[3].push_back(q[3]);
	cout<<t[3]<<endl;
}
void Analyze() {
	cout<<"come to ana"<<endl;
	float mean = 8627;
	TH1D *tabs = new TH1D("tabs","tabs",200,-1500,1500);
	int len = CalibY[0].size();
	for(int i=0;i<len;i++) {tabs->Fill(CalibY[0].at(i)-mean);}
	float rms = tabs->GetRMS();
	TF1 *f1 = new TF1("f1","gaus");
	tabs->Fit("f1");
	tabs->Write();
}

