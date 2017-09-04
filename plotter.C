#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TFile.h"
#include "tdrstyle.C"

TCanvas* compareHistos(TString title, TH1F* histo1, TH1F* histo2, TString leg1, TString leg2)
{

  TCanvas *canvas = new TCanvas(title, title, 500, 500);
  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(1);
  histo1->SetMarkerColor(kRed);
  histo1->SetLineColor(kRed);
  histo1->Scale(1./histo1->Integral());
  histo2->SetMarkerStyle(21);
  histo2->SetMarkerSize(1);
  histo2->SetMarkerColor(kBlue);
  histo2->SetLineColor(kBlue);
  histo2->Scale(1./histo2->Integral());
  histo2->Draw();
  histo1->Draw("same");
  histo2->SetMaximum(1.01);
  TLegend *legend = new TLegend(0.65,0.75,0.9,0.9);   
  legend->SetBorderSize(1);
  legend->AddEntry(histo1, leg1, "lep");
  legend->AddEntry(histo2, leg2, "lep");    
  legend->Draw();

  return canvas;
  
}

TCanvas* compareEfficiencies(TString title, TEfficiency* eff1, TEfficiency* eff2, TString leg1, TString leg2)
{
  
  TCanvas *canvas = new TCanvas(title, title, 500, 500);
  eff1->SetMarkerStyle(20);
  eff2->SetMarkerStyle(21);
  eff1->SetMarkerSize(1);
  eff2->SetMarkerSize(1);
  eff1->SetMarkerColor(kRed);
  eff2->SetMarkerColor(kBlue);
  eff1->SetLineColor(kRed);
  eff2->SetLineColor(kBlue);
  eff1->Draw(); 
  canvas->Update(); 
  auto graph = eff1->GetPaintedGraph(); 
  graph->SetMinimum(0.);
  graph->SetMaximum(1.05); 
  canvas->Update();
  eff2->Draw("same");   
  TLegend *legend = new TLegend(0.75,0.5,0.9,0.6);
  legend->SetBorderSize(1);
  legend->AddEntry(eff1, leg1, "lep");
  legend->AddEntry(eff2, leg2, "lep");    
  legend->Draw();

  return canvas;
  
}


void plotter()
{

  bool save = false;
  TString inputFile = "analysis_results.root";
	
  TFile *file = new TFile(inputFile,"open");
  
  if (!(file->IsOpen())) {
    std::cout<<("File cannot be opened\n");
    return;
  }
	
  ///////////////////////////////////////////////
  // Recommended style macro used as strating /// 
  // point for many CMS plots                 /// 
  ///////////////////////////////////////////////

  setTDRStyle();
  

  //////////////////////////
  // Read the histograms ///
  //////////////////////////
  
  // From kinematics folder
  
  TH1F *h_Zmumu_mass = (TH1F*)gDirectory->Get("/kinematics/h_Zmumu_mass");
  TH1F *h_mu_eta = (TH1F*)gDirectory->Get("/kinematics/h_mu_eta");
  TH1F *h_mu_phi = (TH1F*)gDirectory->Get("/kinematics/h_mu_phi");
  TH1F *h_mu_pt = (TH1F*)gDirectory->Get("/kinematics/h_mu_pt");
  
  // From distance folder (matchings etc ...)
      
  // From efficiency folder  
  
  // From trigger folder (TwinMux quality and BX)  
    
  /////////////////////////
  ////  Now plot the  /////
  //// the histograms /////
  /////////////////////////
  
  TCanvas *c_kinematcis = new TCanvas("Muon (pair) kinematcis", "Muon (pair) kinematcis",750,750);
  
  c_kinematcis->Divide(2,2);
  
  c_kinematcis->cd(1);
  h_Zmumu_mass->Draw();
  
  c_kinematcis->cd(2);
  h_mu_eta->Draw();	
  
  c_kinematcis->cd(3);
  h_mu_pt->Draw();
  
  c_kinematcis->cd(4);
  h_mu_phi->SetMinimum(0.);	
  h_mu_phi->Draw();	


}
