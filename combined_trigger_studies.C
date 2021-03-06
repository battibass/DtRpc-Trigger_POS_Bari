﻿#define combined_trigger_studies_cxx
#include "combined_trigger_studies.h"

#include "TMath.h"
#include "TFile.h"

#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"

#include <iostream>

void combined_trigger_studies::Loop()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //   In a ROOT session, you can do:
  //      root> .L combined_trigger_studies.C
  //      root> combined_trigger_studies t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries


  TFile *outputFile = new TFile(output_file_name,"recreate");
 
  if (!(outputFile->IsOpen())) {
    std::cout<<("Output File cannot be opened\n");
    return;
  }

  outputFile->mkdir("kinematics");
  outputFile->mkdir("efficiencies");
  outputFile->mkdir("distance");
  outputFile->mkdir("trigger");

  // muon kinematics

  outputFile->cd("/kinematics");

  TH1F *h_Zmumu_mass = new TH1F("h_Zmumu_mass", "Z boson candidate mass; m(tag_{#mu},probe_{#mu}) (GeV);# entries", 100, 50., 150.);
  TH1F *h_mu_eta = new TH1F("h_mu_eta", "muon #eta;muon #eta;# entries", 48, -2.4, 2.4);
  TH1F *h_mu_phi = new TH1F("h_mu_phi", "muon #phi;muon #phi (rad);# entries", 48, -pig,  pig);
  TH1F *h_mu_pt  = new TH1F("h_mu_pt" , "muon p_{T};muon p_{T} (GeV/c);# entries", 150, 0., 150.);
  
  // distance

  outputFile->cd("/distance");

  // TH1F *h_dR_seg_muon_MB1 = new TH1F("h_dR_seg_muon_MB1", 
  // 				     "Distance in MB1 between segment and muon;#DeltaR (cm);# entires", 
  // 				     50, 0., 50.);

  // TH1F *h_dR_seg_muon_MB2 = new TH1F("h_dR_seg_muon_MB2", 
  // 				     "Distance in MB1 between segment and muon;#DeltaR (cm);# entires", 
  // 				     50, 0., 50.);

  // TH2F *h_dX_dY_seg_muon_MB1 = new TH2F("h_dX_dY_seg_muon_MB1", 
  // 					"Distance in MB1 between segment and muon;#Deltax (cm);#Deltay (cm)", 
  // 					100, -50., 50., 100, -50., 50.);

  // TH2F *h_dX_dY_seg_muon_MB2 = new TH2F("h_dX_dY_seg_muon_MB2", 
  // 					"Distance DT MB2 between segment and muon;#Deltax (cm);#Deltay (cm)", 
  // 					100, -50., 50., 100, -50., 50.);

  // TH1F *h_dPhi_seg_TrigIn_MB1 = new TH1F("h_dPhi_seg_TrigIn_MB1", 
  // 					 "Distance in MB1 between segment and TwinMux In;#Delta#phi (rad);# entires", 
  // 					 100, 0., .5);

  // TH1F *h_dPhi_seg_TrigIn_MB2 = new TH1F("h_dPhi_seg_TrigIn_MB2", 
  // 					 "Distance in MB2 between segment and TwinMux In;#Delta#phi (rad);# entires", 
  // 					 100, 0., .5);  

  // TH1F *h_dX_MB1_layer_1 = new TH1F("h_dX_MB1_layer_1", "Distance on MB1 layer 1;#Deltax (cm); # entries", 100, 0, 10);
  // TH1F *h_dX_MB1_layer_2 = new TH1F("h_dX_MB1_layer_2", "Distance on MB1 layer 2;#Deltax (cm); # entries", 100, 0, 10);
  
  // TH1F *h_dX_MB2_layer_1 = new TH1F("h_dX_MB2_layer_1", "Distance on MB2 layer 1;#Deltax (cm); # entries", 100, 0, 10);
  // TH1F *h_dX_MB2_layer_2 = new TH1F("h_dX_MB2_layer_2", "Distance on MB2 layer 2;#Deltax (cm); # entires", 100, 0, 10);

  

  
  // efficiency

  outputFile->cd("/efficiencies");

  // TEfficiency* h_eff_twinmux_in_pt_MB1 = new TEfficiency("h_eff_twinmux_in_pt_MB1", 
  // 							 "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);
  // TEfficiency* h_eff_twinmux_in_pt_MB2 = new TEfficiency("h_eff_twinmux_in_pt_MB2", 
  // 							 "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);    
  
  // TEfficiency* h_eff_twinmux_in_eta_MB1 = new TEfficiency("h_eff_twinmux_in_eta_MB1",
  // 							  "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  // TEfficiency* h_eff_twinmux_in_eta_MB2 = new TEfficiency("h_eff_twinmux_in_eta_MB2",
  // 							  "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  // TEfficiency* h_eff_twinmux_in_phi_MB1 = new TEfficiency("h_eff_twinmux_in_phi_MB1",
  // 							  "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  // TEfficiency* h_eff_twinmux_in_phi_MB2 = new TEfficiency("h_eff_twinmux_in_phi_MB2",
  // 							  "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  
  // TEfficiency* h_eff_rpc_pt_MB1 = new TEfficiency("h_eff_rpc_pt_MB1", 
  // 						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);
  // TEfficiency* h_eff_rpc_pt_MB2 = new TEfficiency("h_eff_rpc_pt_MB2", 
  // 						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);    
  
  // TEfficiency* h_eff_rpc_eta_MB1 = new TEfficiency("h_eff_rpc_eta_MB1",
  // 						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  // TEfficiency* h_eff_rpc_eta_MB2 = new TEfficiency("h_eff_rpc_eta_MB2",
  // 						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  // TEfficiency* h_eff_rpc_phi_MB1 = new TEfficiency("h_eff_rpc_phi_MB1",
  // 						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  // TEfficiency* h_eff_rpc_phi_MB2 = new TEfficiency("h_eff_rpc_phi_MB2",
  // 						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);

  // TEfficiency* h_eff_combined_toy_pt_MB1 = new TEfficiency("h_eff_combined_toy_pt_MB1", 
  // 						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);
  // TEfficiency* h_eff_combined_toy_pt_MB2 = new TEfficiency("h_eff_combined_toy_pt_MB2", 
  // 						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);    
  
  // TEfficiency* h_eff_combined_toy_eta_MB1 = new TEfficiency("h_eff_combined_toy_eta_MB1",
  // 						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  // TEfficiency* h_eff_combined_toy_eta_MB2 = new TEfficiency("h_eff_combined_toy_eta_MB2",
  // 						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  // TEfficiency* h_eff_combined_toy_phi_MB1 = new TEfficiency("h_eff_combined_toy_phi_MB1",
  // 						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  // TEfficiency* h_eff_combined_toy_phi_MB2 = new TEfficiency("h_eff_combined_toy_phi_MB2",
  // 						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);

  // TEfficiency* h_eff_twinmux_out_pt_MB1 = new TEfficiency("h_eff_twinmux_out_pt_MB1", 
  // 							 "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);
  // TEfficiency* h_eff_twinmux_out_pt_MB2 = new TEfficiency("h_eff_twinmux_out_pt_MB2", 
  // 							 "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 20,0.,100.);    
  
  // TEfficiency* h_eff_twinmux_out_eta_MB1 = new TEfficiency("h_eff_twinmux_out_eta_MB1",
  // 							  "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  // TEfficiency* h_eff_twinmux_out_eta_MB2 = new TEfficiency("h_eff_twinmux_out_eta_MB2",
  // 							  "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  // TEfficiency* h_eff_twinmux_out_phi_MB1 = new TEfficiency("h_eff_twinmux_out_phi_MB1",
  // 							  "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  // TEfficiency* h_eff_twinmux_out_phi_MB2 = new TEfficiency("h_eff_twinmux_out_phi_MB2",
  // 							  "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);


  // TwinMux quality and BX

  outputFile->cd("/trigger");

  // TH1F *h_BX_twinmux_in_MB1 = new TH1F("h_BX_twinmux_in_MB1", "TwinMux In BX for MB1; BX; # entries", 7, -3.5, 3.5);
  // TH1F *h_BX_twinmux_in_MB2 = new TH1F("h_BX_twinmux_in_MB2", "TwinMux In BX for MB2; BX; # entries", 7, -3.5, 3.5);

  // TH1F *h_BX_twinmux_out_MB1 = new TH1F("h_BX_twinmux_out_MB1", "TwinMux Out BX for MB1; BX; # entries", 7, -3.5, 3.5);
  // TH1F *h_BX_twinmux_out_MB2 = new TH1F("h_BX_twinmux_out_MB2", "TwinMux Out BX for MB2; BX; # entries", 7, -3.5, 3.5);

  // TH1F *h_qual_twinmux_in_MB1 = new TH1F("h_qual_twinmux_in_MB1", "TwinMux In quality for MB1; quality; # entries", 7, -0.5, 6.5);
  // TH1F *h_qual_twinmux_in_MB2 = new TH1F("h_qual_twinmux_in_MB2", "TwinMux In quality for MB2; quality; # entries", 7, -0.5, 6.5);
 


  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast() < n_events ?
                      fChain->GetEntriesFast() : n_events ; 

  Long64_t nbytes = 0, nb = 0;
   
  Int_t count_muon = 0;
  
  std::cout << "[combined_trigger_studies::Loop] processing : " 
	    << nentries << " entries " << std::endl;
 
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if(jentry % 10000 == 0) 
      std::cout << "[combined_trigger_studies::Loop] processed : " << jentry << " entries\r" << std::flush;
    
    auto tnpPairs = combined_trigger_studies::TnPSelection(min_TnP_mass,max_TnP_mass);

    for(const auto & pair : tnpPairs) {

      count_muon++;

      // the first  element of the pair is the index of the tag muon
      // the second element of the pair is the index of the probe muon
      // we want to study the probe

      Int_t iTag   = pair.first;
      Int_t iProbe = pair.second;
      
      TLorentzVector tag_vec;
      tag_vec.SetXYZM(Mu_px->at(iTag),
		      Mu_py->at(iTag),
		      Mu_pz->at(iTag),
		      0.106);
      
      TLorentzVector probe_vec;
      probe_vec.SetXYZM(Mu_px->at(iProbe),
			Mu_py->at(iProbe),
			Mu_pz->at(iProbe),
			0.106);
      
      h_mu_pt->Fill(probe_vec.Pt());
      h_mu_phi->Fill(probe_vec.Phi());
      h_mu_eta->Fill(probe_vec.Eta());
      
      Float_t mass = (probe_vec + tag_vec).M();
      h_Zmumu_mass->Fill(mass);
      
      
      // ************************************
      // We want to look for segmets close to 
      // the reconstructed mu in MB1 and MB2
      // ************************************

      // **************************************
      // 1st part step1:
      // select the segment closest to the
      // propagated muon track and tune the cut 
      // between segments and reconstructed 
      // muons to a reasonable value
      // **************************************

      bool has_match_DT_MB1_muon = false;
      bool has_match_DT_MB2_muon = false;
      Int_t dtsegment_index[2] = {-999, -999};

      Float_t dr_seg_muon_cut = 999.; // parameter to be tuned as part of the exercise

      Float_t minDrMB1 = 999.;
      Float_t minDrMB2 = 999.;      

      /// Loop on DT 4D segments ///
      for(Int_t iDtSegm = 0; iDtSegm < Ndtsegments; ++iDtSegm) {
	
	if(dtsegm4D_station->at(iDtSegm) > 2)   continue; // MB3 and MB4 do not matter
	
	if(dtsegm4D_hasZed->at(iDtSegm)   <= 0) continue; // Miminal quality cuts:
	if(dtsegm4D_phinhits->at(iDtSegm) <  4) continue; // ask for z view and >= 4 hits in phi  			

	///// station 1 /////
	
	if(Mu_sector_MB1->at(iProbe)      > 0 && // non valid matches with MB1 are -999
	   dtsegm4D_station->at(iDtSegm) == 1 &&
	   dtsegm4D_sector->at(iDtSegm) == Mu_sector_MB1->at(iProbe) &&
	   dtsegm4D_wheel->at(iDtSegm)  == Mu_wheel_MB1->at(iProbe)  ) {

	  // ...

	  // h_dX_dY_seg_muon_MB1->Fill(dXMB1,dYMB1);

	  // ...
	  
	}
	
	///// station 2 /////
	
 	if(Mu_sector_MB2->at(iProbe)      > 0  && // non valid matches with MB2 are -999
	   dtsegm4D_station->at(iDtSegm) == 2 &&
	   dtsegm4D_sector->at(iDtSegm) == Mu_sector_MB2->at(iProbe) &&
	   dtsegm4D_wheel->at(iDtSegm)  == Mu_wheel_MB2->at(iProbe)  ) {
	  
	  // ...
	  
	  // h_dX_dY_seg_muon_MB2->Fill(dXMB2,dYMB2);
	  
	  
	}
	
      }

      // ...

      // ************************************
      // We want to look for *all* DT trigger
      // primitives (TwinMux input) close
      // to the segments found above
      // ************************************

      // ************************************
      // 1st part step2:
      // tune the matching between segments
      // and trigger (can do the same for
      // TwinMux Out later on!)
      // ************************************
     
      std::vector<Int_t> twinmux_in_MB1;
      std::vector<Int_t> twinmux_in_MB2;
      Float_t dphi_twinumx_in_seg_cut = 999.; // parameter to be tuned as part of the exercise

      for(Int_t iTrig = 0; iTrig < NdtltTwinMuxIn; ++iTrig) {
	
	if(has_match_DT_MB1_muon) {

	  // ...
	  
	}
	  
	if(has_match_DT_MB2_muon) {

	  // ...
	    
	}

      } // loop on TwinMux In


      // ************************************
      // We want to look at the quality, BX
      // assignment and efficiency of trigger
      // primitives from "good" muons
      // ************************************

      // ************************************
      // 1st part step3:
      // plot BX and quality for all the
      // matched trigger primitives and 
      // compute primitive efficiency (w.r.t. 
      // segments) to trigger at BX == 0
      // ************************************

      bool has_twinmux_in_BX0_MB1 = false;
      bool has_twinmux_in_BX0_MB2 = false;

      for (const auto & iTrig : twinmux_in_MB1)
	{

	  // ...
	  
	}	 

      for (const auto & iTrig : twinmux_in_MB2)
	{
	  
	  // ...

	}	 

      // ...

      // ************************************
      // We want to study RPC trigger 
      // quantities
      // ************************************

      // ************************************
      // 2nd part step1:
      // extrapolate DT segment on inner and
      // outer layer of the RPC and match
      // with the closest RPC recHit
      // ************************************

      bool has_dt_extrapolation_InOut_MB1 = false;
      bool has_dt_extrapolation_InOut_MB2 = false;

      bool has_rpc_match_MB1 = false;
      bool has_rpc_match_MB2 = false;

      // indexes in arrays are inner and outer layer
      bool has_dt_extrapolation_MB1[2] = {false, false};
      
      Float_t minDxMB1[2] = { 999., 999.};
	  
      Float_t rpcCluSizeMB1[2] = {0., 0.}; 
      Float_t stripWidthMB1[2] = {0., 0.};	    

      bool has_dt_extrapolation_MB2[2] = {false, false};
      
      Float_t minDxMB2[2] = { 999., 999.};
	  
      Float_t rpcCluSizeMB2[2] = {0., 0.}; 
      Float_t stripWidthMB2[2] = {0., 0.};	    

      // MB1

      if( has_match_DT_MB1_muon ) {
	for(Int_t iOnRpc = 0; iOnRpc < NDTsegmentonRPC->at(dtsegment_index[0]); ++iOnRpc) {	  

	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[0],iOnRpc) == 1) { // inner ring
	    
	    // ...
	      
	  }
	  
	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[0],iOnRpc) == 2) { // outer ring
	    
	    // ...
	      
	  }
	  
	}
      }
    	
      // MB2
      
      if( has_match_DT_MB2_muon ) {
	
	for(Int_t iOnRpc = 0; iOnRpc < NDTsegmentonRPC->at(dtsegment_index[1]); ++iOnRpc) {
	  
	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[1],iOnRpc) == 1) { // inner ring

	    // ...
	  
	  }
	
	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[1],iOnRpc) == 2) { // outer ring
	  
	    // ...
	    	    
	  }	  
	}
      }
    
      // ************************************
      // 2nd part step2:
      // Tune matching between extrapolated
      // recHit and RPC recHit: our goal is 
      // to evaluate muon efficiency of RPC
      // pseudo-segments with 2 hits
      // ************************************      
      
      // RPC matching cuts
      Float_t cluster_size_cut = 99.; // parameter to be tuned as part of the exercise
      Float_t range_strips = 4.;

      if (has_dt_extrapolation_MB1[0] && has_dt_extrapolation_MB1[1] ) {
	
	// ...

      }
      
      if (has_dt_extrapolation_MB2[0] && has_dt_extrapolation_MB2[1] ) {
	
	// ...

      }
      
      // ...

      
      // ************************************
      // 2nd part step3:
      // Look how RPC improve DT primitive
      // BX assignment looking at the TWinMux
      // Out BX distribution and efficiency 
      // ************************************	
	
      std::vector<Int_t> twinmux_out_MB1;
      std::vector<Int_t> twinmux_out_MB2;
      Float_t dphi_twinumx_out_seg_cut = 999.; // parameter to be tuned as part of the exercise

      for(Int_t iTrig = 0; iTrig < NdtltTwinMuxOut; ++iTrig) {
	
	if(has_match_DT_MB1_muon) {

	  // ...
	  
	}
	  
	if(has_match_DT_MB2_muon) {	    

	  // ...
	    
	}

      } // loop on TwinMux Out

      bool has_twinmux_out_BX0_MB1 = false;
      bool has_twinmux_out_BX0_MB2 = false;

      for (const auto & iTrig : twinmux_out_MB1)
	{

	  // ...
	  
	}	 

      for (const auto & iTrig : twinmux_out_MB2)
	{

	  // ...

	}	 

      // ...

    } //loop on muons			
    
  } // loop on event
  
  std::cout << std::endl;

  outputFile->Write();
  outputFile->Close();

}


vector<std::pair<Int_t,Int_t>> combined_trigger_studies::TnPSelection(Float_t minMass,
								      Float_t maxMass) const
{

  vector<std::pair<Int_t,Int_t>> pairs;

  for(Int_t iTag = 0; iTag < Nmuons; ++iTag) 
    {
    
      TLorentzVector tagVec;
      tagVec.SetXYZM(Mu_px->at(iTag),
		     Mu_py->at(iTag),
		     Mu_pz->at(iTag),
		     0.106);
      
      bool tagQuality = 
	Mu_isMuGlobal->at(iTag)  == 1 &&
	Mu_isMuTracker->at(iTag) == 1 &&
	std::abs(Mu_dxy_glb->at(iTag)) < Dxy_cut &&
	Mu_normchi2_glb->at(iTag)      < muchi2_cut &&
	Mu_numberOfHits_sta->at(iTag)        > N_hits_cut     && 
	Mu_numberOfPixelHits_glb->at(iTag)   >= npix_cut      &&
	Mu_numberOfTrackerHits_glb->at(iTag) >= ntkr_cut      &&
 	Mu_tkIsoR03_glb->at(iTag) / tagVec.Pt() < tkr_iso_cut &&
	tagVec.Pt() > min_pt_tag ;
      
      if(tagQuality && HasTrigger(tagVec,muon_hlt_dR)) 
	{
	  
	  for(Int_t iProbe = 0; iProbe < Nmuons; ++iProbe) 
	    {
	      
	      if (iTag == iProbe) 
		continue;

	      TLorentzVector probeVec;
	      probeVec.SetXYZM(Mu_px->at(iProbe),
			       Mu_py->at(iProbe),
			       Mu_pz->at(iProbe),
			       0.106);
	      
	      bool probeQuality =
		Mu_isMuTracker->at(iProbe) == 1 &&
		Mu_isMuGlobal->at(iProbe)  == 1 &&
		std::abs(Mu_dz_glb->at(iProbe))  < Dz_cut  &&
		Mu_numberOfPixelHits_glb->at(iProbe)   >= npix_cut        &&
		Mu_numberOfTrackerHits_glb->at(iProbe) >= ntkr_cut        &&
		Mu_tkIsoR03_glb->at(iProbe) / probeVec.Pt() < tkr_iso_cut &&
		std::abs(probeVec.Eta()) < max_eta_probe &&
		probeVec.Pt() > min_pt_probe;

	      if (probeQuality && std::abs(Mu_dz_glb->at(iTag) - Mu_dz_glb->at(iProbe)) < Dz_cut)
		{
		  Float_t mass = (tagVec + probeVec).M();

		  if (mass > minMass && mass < maxMass)
		    {
		      pairs.push_back(std::make_pair(iTag,iProbe));
		      break; // just one probe per tag
		    }
		} 
	    }
	}
    }
  
  return pairs;

}

bool combined_trigger_studies::HasTrigger(const TLorentzVector & muon,
					  Float_t deltaR) const 
{

  for(Int_t iTrig = 0; iTrig < NhltFilters; ++iTrig)
    {

      if (hlt_filter->at(iTrig).Contains(trigger_filter_name))
	{
	  TLorentzVector trigVec;
	  trigVec.SetPtEtaPhiM(hlt_filter_pt->at(iTrig),
			       hlt_filter_eta->at(iTrig),
			       hlt_filter_phi->at(iTrig),
			       0.106);
	  
	  if (trigVec.DeltaR(muon) < deltaR)
	    return true;
  	}
    }
  
  return false;
  
}

Float_t combined_trigger_studies::PhiConversion(Int_t phi_In, Int_t sector) const
{	
  Float_t locphi = phi_In / 4096.0;
  Float_t newphi = locphi + ((sector-1) * (pig/6.));
  if (newphi > pig) newphi -= 2 * pig;
  
  return newphi;
}
