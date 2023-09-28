#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "TRotation.h"
#include "TEfficiency.h"

//#include "muonhits_canvas.C"
//#include "muonhits_canvas_new.C"


using namespace ROOT::Math;

void Calohits_macro(){


  TChain* c_hits = new TChain("c_hits");
  c_hits->Add("ntuple_CaloHit_onlyBIB.root/MyLCTuple"); 
  
 // ofstream myfile;
  //myfile.open ("hits_muchambes.txt");
  Int_t *ca_ori;
  Int_t *ca_ci0;
  Float_t *ca_tim;
  Float_t *ca_pox;
  Float_t *ca_poy;
  Float_t *ca_poz;
  Float_t *ca_ene;
  
  gStyle->SetOptStat("nemrou");
  
  int num_ca=1000000;
  Int_t n_cah;
  ca_ori = (int*) malloc(sizeof(int)*num_ca);
  //ca_typ = (int*) malloc(sizeof(int)*num_ca);
  ca_ci0 = (int*) malloc(sizeof(int)*num_ca);
	ca_tim = (float*) malloc(sizeof(float)*num_ca);
	ca_pox = (float*) malloc(sizeof(float)*num_ca);
	ca_poy = (float*) malloc(sizeof(float)*num_ca);
	ca_poz = (float*) malloc(sizeof(float)*num_ca);
	ca_tim = (float*) malloc(sizeof(float)*num_ca);
	ca_ene = (float*) malloc(sizeof(float)*num_ca);
 
   
   
  
  
  c_hits->SetBranchAddress("ncah", &n_cah);
  c_hits->SetBranchAddress("caori", ca_ori);
  //c_hits->SetBranchAddress("catyp", ca_typ);
  c_hits->SetBranchAddress("catim", ca_tim);
  c_hits->SetBranchAddress("capox", ca_pox);
  c_hits->SetBranchAddress("capoy", ca_poy);
  c_hits->SetBranchAddress("capoz", ca_poz);
  c_hits->SetBranchAddress("caci0", ca_ci0);
  c_hits->SetBranchAddress("catim", ca_tim);
  c_hits->SetBranchAddress("caene", ca_ene);
  
   
  int system, layer;
  
  
  //TH2F* barrel_HCalHits_histo= new TH2F("Barrel Hcal", "Barrel Hcal", 224, -3360, 3360, 224, -3360,3360);
  //TH2F* endcap_HCalHits_histo= new TH2F("Endcap Hcal", "Endcap Hcal", 218, -3300, 3300, 218, -3300,3300);
  //TH2F* ring_HCalHits_histo= new TH2F("Ring Hcal", "Ring Hcal", 218, -3300,3300, 218, -3300, 3300);
  
  TH3D* barrel_HCalHits_histo= new TH3D("Barrel Hcal", "Barrel Hcal", 224, -3360, 3360, 224, -3360,3360, 61,-0.5,60.5);
  TH3D* endcap_HCalHits_histo= new TH3D("Endcap Hcal", "Endcap Hcal", 218, -3300, 3300, 218, -3300,3300, 61,-0.5,60.5);
  TH3D* ring_HCalHits_histo= new TH3D("Ring Hcal", "Ring Hcal", 218, -3300,3300, 218, -3300, 3300, 10, -0.5, 9.5);
  TH1F* time_Barrel= new TH1F("hit time in HCal Barrel","hit time in HCal Barrel",100,-0.2,0.4);
  TH1F* time_Endcap= new TH1F("hit time in HCal Endcap","hit time in HCal Endcap",100,-0.2,0.4);
  TH1F* energy_Endcap= new TH1F("hit energy in HCal Endcap","hit energy in HCal Endcap",80,0,0.4);
  TH1F* energy_Barrel= new TH1F("hit energy in HCal Barrel","hit energy in HCal Barrel",80,0,0.4);
  
 float hits_layer_barrel[60]={0};
 float hits_layer_endcap[60]={0};
 //float hits_layer_other[60]={0};
 
  
  for(unsigned int ientry=0; ientry<c_hits->GetEntries(); ++ientry){
  	c_hits->GetEntry(ientry);
 		for(unsigned int ihit=0; ihit<n_cah; ++ihit){
 			system= ca_ci0[ihit] & 31;
 			layer= ca_ci0[ihit]>>19 & 511;
 			if(ca_ori[ihit]==10){ // Barrel
 				barrel_HCalHits_histo->Fill(ca_pox[ihit],ca_poy[ihit],layer);
 				time_Barrel->Fill(ca_tim[ihit]);
 				energy_Barrel->Fill(ca_ene[ihit]);
 				hits_layer_barrel[layer]+=1.;
 				//cout<<"Barrel"<<endl;
 				//cout<<"system: "<<system<<"    layer: "<<layer<<endl;
 			}
 			
 			if(ca_ori[ihit]==11){ // Endcap
 				endcap_HCalHits_histo->Fill(ca_pox[ihit],ca_poy[ihit],layer);
 				time_Endcap->Fill(ca_tim[ihit]);
 				energy_Endcap->Fill(ca_ene[ihit]);
 				hits_layer_endcap[layer]+=1.;
 				//cout<<"Endcap"<<endl;
 				//cout<<"system: "<<system<<"    layer: "<<layer<<endl;
 			}
 			
 			if(ca_ori[ihit]==12){ // Other
 				ring_HCalHits_histo->Fill(ca_pox[ihit],ca_poy[ihit],layer);
 				//cout<<"Other"<<endl;
 				//cout<<"system: "<<system<<"    layer: "<<layer<<endl;
 			}
 		}
  	
  }
  
  for(unsigned int i=0; i<60;i++){
  
  	cout<<"hits_layer_barrel["<<i<<"]: "<<hits_layer_barrel[i]<<endl;
  	cout<<"hits_layer_endcap["<<i<<"]: "<<hits_layer_endcap[i]<<endl;
  	cout<<"after dividing by 50"<<endl;
  	hits_layer_barrel[i]=hits_layer_barrel[i]/50;
  	hits_layer_endcap[i]=hits_layer_endcap[i]/50;
  	
  	cout<<"hits_layer_barrel["<<i<<"]: "<<hits_layer_barrel[i]<<endl;
  	cout<<"hits_layer_endcap["<<i<<"]: "<<hits_layer_endcap[i]<<endl;
  	
  	
  }
  
  time_Barrel->Scale(1./time_Barrel->Integral());
  time_Endcap->Scale(1./time_Endcap->Integral());
  energy_Barrel->Scale(1./energy_Barrel->Integral());
  energy_Endcap->Scale(1./energy_Endcap->Integral());
  
  gStyle->SetOptStat(0);
  
 /// Project a 3-d histogram into 1 or 2-d histograms depending on the
 /// option parameter, which may contain a combination of the characters x,y,z,e
 ///  - option = "x" return the x projection into a TH1D histogram
 ///  - option = "y" return the y projection into a TH1D histogram
 ///  - option = "z" return the z projection into a TH1D histogram
 ///  - option = "xy" return the x versus y projection into a TH2D histogram
 ///  - option = "yx" return the y versus x projection into a TH2D histogram
 ///  - option = "xz" return the x versus z projection into a TH2D histogram
 ///  - option = "zx" return the z versus x projection into a TH2D histogram
 ///  - option = "yz" return the y versus z projection into a TH2D histogram
 ///  - option = "zy" return the z versus y projection into a TH2D histogram
 
 // !!! axes are inverted --> in order to have the 2D projection of TH3 you should use Project3D("yx")
  
  TH2D* barrel_HCalHits[12];
  TH2D* endcap_HCalHits[12];
  TH3D *histo;
  TCanvas* c2= new TCanvas("c2","c2",1000,1000);
  gStyle->SetOptStat(0);
  c2->Divide(6,2);
  int j=1; 
  int l=0;
  
  float aver_hits_barrel[12]={0}; // barrel -> average number of hits per event for each group of 5 layers
  float aver_hits_endcap[12]={0}; // endcap -> average number of hits per event for each group of 5 layers
  for(int i=0; i<60;i+=5){
  	for(int k=0; k<5;k++){
  		cout<<i<<" "<<k<<endl;
  		aver_hits_barrel[l]= aver_hits_barrel[l]+hits_layer_barrel[i+k];
  		aver_hits_endcap[l]= aver_hits_endcap[l]+hits_layer_endcap[i+k];
  		cout<<hits_layer_barrel[i+k]<<endl;
  	}
  	
  	cout<<"Barrel -> averge number of hits per event, layers "<<i<<"-"<<i+5<<":  "<<aver_hits_barrel[l]<<endl;
  	cout<<"Endcap -> averge number of hits per event, layers "<<i<<"-"<<i+5<<":  "<<aver_hits_endcap[l]<<endl;
  	l++;
  }
  
  //barrel
  for(int i=0; i<60; i+=5){
  	c2->cd(j);
  	barrel_HCalHits_histo->GetZaxis()->SetRange(i+1,i+5);
  	//cout<<"i: "<<i<<endl;
  	//cout<<"j-1: "<<j-1<<endl;
  	barrel_HCalHits[j-1]= (TH2D*) barrel_HCalHits_histo->Project3D("yx")->Clone();
  	barrel_HCalHits[j-1]->SetTitle(Form("Hcal - Barrel Layers:%d-%d",i,i+5));
  	barrel_HCalHits[j-1]->GetXaxis()->SetTitle("x [mm]");
  	barrel_HCalHits[j-1]->GetYaxis()->SetTitle("y [mm]");
  	barrel_HCalHits[j-1]->GetYaxis()->SetTitleOffset(1.3);
  	barrel_HCalHits[j-1]->Draw("COLZ");
  	j++;
  }
  
  //endcap
  j=1;
  TCanvas* c2_1= new TCanvas("c2_1","c2_1",1000,1000);
  c2_1->Divide(6,2);
  for(int i=0; i<60; i+=5){
  	c2_1->cd(j);
  	endcap_HCalHits_histo->GetZaxis()->SetRange(i+1,i+5);
  	//cout<<"i: "<<i<<endl;
  	//cout<<"j-1: "<<j-1<<endl;
  	endcap_HCalHits[j-1]= (TH2D*) endcap_HCalHits_histo->Project3D("yx")->Clone();
  	endcap_HCalHits[j-1]->SetTitle(Form("Hcal - Endcap Layers:%d-%d",i,i+5));
  	endcap_HCalHits[j-1]->GetXaxis()->SetTitle("x [mm]");
  	endcap_HCalHits[j-1]->GetYaxis()->SetTitle("y [mm]");
  	endcap_HCalHits[j-1]->GetYaxis()->SetTitleOffset(1.4);
  	endcap_HCalHits[j-1]->Draw("COLZ");
  	j++;
  }
  
  
	
	TCanvas* c4 = new TCanvas("c4","c4", 800,600);
  gStyle->SetOptStat(0);
	barrel_HCalHits[0]->Draw("COLZ");
	
	TCanvas* c5 = new TCanvas("c5","c5", 800,600);
  gStyle->SetOptStat(0);
	barrel_HCalHits[2]->Draw("COLZ");
	
	TCanvas* c6 = new TCanvas("c6","c6", 800,600);
  gStyle->SetOptStat(0);
	barrel_HCalHits[4]->Draw("COLZ");
	
	TCanvas* c7 = new TCanvas("c7","c7", 800,600);
  gStyle->SetOptStat(0);
	barrel_HCalHits[6]->Draw("COLZ");
	
	TCanvas* c8 = new TCanvas("c8","c8", 800,600);
  gStyle->SetOptStat(0);
	barrel_HCalHits[8]->Draw("COLZ");
	
	TCanvas* c9 = new TCanvas("c9","c9", 800,600);
  gStyle->SetOptStat(0);
	barrel_HCalHits[10]->Draw("COLZ");
	
	TCanvas* c4_1 = new TCanvas("c4_1","c4_1", 800,600);
  gStyle->SetOptStat(0);
	endcap_HCalHits[0]->Draw("COLZ");
	
	TCanvas* c5_1 = new TCanvas("c5_1","c5_1", 800,600);
  gStyle->SetOptStat(0);
	endcap_HCalHits[2]->Draw("COLZ");
	
	TCanvas* c6_1 = new TCanvas("c6_1","c6_1", 800,600);
  gStyle->SetOptStat(0);
	endcap_HCalHits[4]->Draw("COLZ");
	
	TCanvas* c7_1 = new TCanvas("c7_1","c7_1", 800,600);
  gStyle->SetOptStat(0);
	endcap_HCalHits[6]->Draw("COLZ");
	
	TCanvas* c8_1 = new TCanvas("c8_1","c8_1", 800,600);
  gStyle->SetOptStat(0);
	endcap_HCalHits[8]->Draw("COLZ");
	
	TCanvas* c9_1 = new TCanvas("c9_1","c9_1", 800,600);
  gStyle->SetOptStat(0);
	endcap_HCalHits[10]->Draw("COLZ");
	
	barrel_HCalHits_histo->GetZaxis()->SetRange(0,-1);
	TH2D* barrel_HCalHits_all = (TH2D*) barrel_HCalHits_histo->Project3D("yx")->Clone();
  
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	barrel_HCalHits_all->SetTitle("HCal - Barrel");
	barrel_HCalHits_all->GetXaxis()->SetTitle("x [mm]");
	barrel_HCalHits_all->GetYaxis()->SetTitle("y [mm]");
	barrel_HCalHits_all->Draw("COLZ");
	
	/*TCanvas *c2 = new TCanvas("c2","c2",800,600);
	endcap_HCalHits_histo->SetTitle("HCal - Endcap");
	endcap_HCalHits_histo->GetXaxis()->SetTitle("x [mm]");
	endcap_HCalHits_histo->GetYaxis()->SetTitle("y [mm]");
	endcap_HCalHits_histo->Draw("COLZ");
	
	TCanvas *c3 = new TCanvas("c3","c3",800,600);
	ring_HCalHits_histo->SetTitle("HCal - Other");
	ring_HCalHits_histo->GetXaxis()->SetTitle("x [mm]");
	ring_HCalHits_histo->GetYaxis()->SetTitle("y [mm]");
	ring_HCalHits_histo->Draw("COLZ");*/
   
  TCanvas *c10 = new TCanvas("c10","c10",800,600);
  time_Barrel->GetXaxis()->SetTitle("time [ns]");
  time_Barrel->GetYaxis()->SetTitle("Normalized to unity");
  time_Barrel->SetLineColor(kBlue);
  time_Barrel->SetLineWidth(2);
  time_Barrel->Draw("hist");
  
  TCanvas *c11 = new TCanvas("c11","c11",800,600);
  time_Endcap->GetXaxis()->SetTitle("time [ns]");
  time_Endcap->GetYaxis()->SetTitle("Normalized to unity");
  time_Endcap->SetLineColor(kBlue);
  time_Endcap->SetLineWidth(2);
  time_Endcap->Draw("hist");
  
  TCanvas *c12 = new TCanvas("c12","c12",800,600);
  energy_Barrel->GetXaxis()->SetTitle("energy [GeV]");
  energy_Barrel->GetYaxis()->SetTitle("Normalized to unity");
  energy_Barrel->SetLineColor(kBlue);
  energy_Barrel->SetLineWidth(2);
  energy_Barrel->Draw("hist");
  
  TCanvas *c13 = new TCanvas("c13","c13",800,600);
  energy_Endcap->GetXaxis()->SetTitle("energy [GeV]");
  energy_Endcap->GetYaxis()->SetTitle("Normalized to unity");
  energy_Endcap->SetLineColor(kBlue);
  energy_Endcap->SetLineWidth(2);
  energy_Endcap->Draw("hist");
  
  
  TFile *rootFile = new TFile("CaloHit_onlyBIB.root","RECREATE");
	time_Barrel->Write();
	time_Endcap->Write();
	energy_Barrel->Write();
	energy_Endcap->Write();
	rootFile->Close();
  
  free(ca_ori);  free(ca_pox); free(ca_poy); free(ca_poz); free(ca_ci0); free(ca_tim); 
  
}
