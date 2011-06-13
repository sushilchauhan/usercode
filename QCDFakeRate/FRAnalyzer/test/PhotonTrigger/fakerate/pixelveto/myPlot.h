//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 10 22:14:15 2011 by ROOT version 5.26/00e
// from TTree myEvent/a tree with histograms
// found on file: FR_Data_AOD.root
//////////////////////////////////////////////////////////

#ifndef myPlot_h
#define myPlot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h> 
#include <TGraphAsymmErrors.h>
#include <map>
//#include <map>
#include <vector>
#ifdef __MAKECINT__ 
#pragma link C++ class vector<bool>+;
#endif 
 
using namespace std;        
using namespace ROOT;


class myPlot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //My Variables

   bool keepThisPhoton;
   bool isPassed;

   // Declaration of leaf types
   Int_t           nevents;
   UInt_t          run;
   UInt_t          event;
   UInt_t          luminosityBlock;
   UInt_t          beamCrossing;
   UInt_t          totalIntensityBeam1;
   UInt_t          totalIntensityBeam2;
   Float_t         avgInsDelLumi;
   Float_t         avgInsDelLumiErr;
   Float_t         avgInsRecLumi;
   Float_t         avgInsRecLumiErr;
   Int_t           ntriggers;
   vector<string>  *triggernames;
   vector<int>     *triggerprescales;
   vector<bool>    *ifTriggerpassed;
   Float_t         trobjpt[100][100][10];   //[ntriggers]
   Float_t         trobjeta[100][100][10];   //[ntriggers]
   Float_t         trobjphi[100][100][10];   //[ntriggers]
   Int_t           Vertex_n;
   Float_t         Vertex_x[200];   //[Vertex_n]
   Float_t         Vertex_y[200];   //[Vertex_n]
   Float_t         Vertex_z[200];   //[Vertex_n]
   Int_t           Vertex_tracksize[200];   //[Vertex_n]
   Int_t           Vertex_ndof[200];   //[Vertex_n]
   Float_t         Vertex_chi2[200];   //[Vertex_n]
   Float_t         Vertex_d0[200];   //[Vertex_n]
   Bool_t          Vertex_isFake[200];   //[Vertex_n]
   Bool_t          Scraping_isScrapingEvent;
   Int_t           Scraping_numOfTracks;
   Float_t         Scraping_fractionOfGoodTracks;
   Int_t           Jet_n;
   Float_t         Jet_px[200];   //[Jet_n]
   Float_t         Jet_py[200];   //[Jet_n]
   Float_t         Jet_E[200];   //[Jet_n]
   Float_t         Jet_pz[200];   //[Jet_n]
   Float_t         Jet_vx[200];   //[Jet_n]
   Float_t         Jet_vy[200];   //[Jet_n]
   Float_t         Jet_vz[200];   //[Jet_n]
   Float_t         Jet_pt[200];   //[Jet_n]
   Float_t         Jet_eta[200];   //[Jet_n]
   Float_t         Jet_phi[200];   //[Jet_n]
   Float_t         Jet_emEnergyFraction[200];   //[Jet_n]
   Float_t         Jet_energyFractionHadronic[200];   //[Jet_n]
   Int_t           Jet_hitsInN90[200];   //[Jet_n]
   Int_t           Jet_n90Hits[200];   //[Jet_n]
   Int_t           Jet_nTowers[200];   //[Jet_n]
   Float_t         Jet_fHPD[200];   //[Jet_n]
   Float_t         Jet_fRBX[200];   //[Jet_n]
   Float_t         Jet_RHF[200];   //[Jet_n]
   Float_t         Jet_jecCorr[200];   //[Jet_n]
   Float_t         ucJet_px[200];   //[Jet_n]
   Float_t         ucJet_py[200];   //[Jet_n]
   Float_t         ucJet_E[200];   //[Jet_n]
   Float_t         ucJet_pz[200];   //[Jet_n]
   Float_t         ucJet_pt[200];   //[Jet_n]
   Float_t         ucJet_eta[200];   //[Jet_n]
   Float_t         ucJet_phi[200];   //[Jet_n]
   Int_t           pfJet_n;
   Float_t         pfJet_px[200];   //[pfJet_n]
   Float_t         pfJet_py[200];   //[pfJet_n]
   Float_t         pfJet_E[200];   //[pfJet_n]
   Float_t         pfJet_pz[200];   //[pfJet_n]
   Float_t         pfJet_vx[200];   //[pfJet_n]
   Float_t         pfJet_vy[200];   //[pfJet_n]
   Float_t         pfJet_vz[200];   //[pfJet_n]
   Float_t         pfJet_pt[200];   //[pfJet_n]
   Float_t         pfJet_eta[200];   //[pfJet_n]
   Float_t         pfJet_phi[200];   //[pfJet_n]
   Float_t         pfJet_jecCorr[200];   //[pfJet_n]
   Float_t         ucpfJet_px[200];   //[pfJet_n]
   Float_t         ucpfJet_py[200];   //[pfJet_n]
   Float_t         ucpfJet_E[200];   //[pfJet_n]
   Float_t         ucpfJet_pz[200];   //[pfJet_n]
   Float_t         ucpfJet_pt[200];   //[pfJet_n]
   Float_t         ucpfJet_eta[200];   //[pfJet_n]
   Float_t         ucpfJet_phi[200];   //[pfJet_n]
   Int_t           Electron_n;
   Float_t         Electron_px[200];   //[Electron_n]
   Float_t         Electron_py[200];   //[Electron_n]
   Float_t         Electron_pz[200];   //[Electron_n]
   Float_t         Electron_vx[200];   //[Electron_n]
   Float_t         Electron_vy[200];   //[Electron_n]
   Float_t         Electron_vz[200];   //[Electron_n]
   Float_t         Electron_pt[200];   //[Electron_n]
   Float_t         Electron_eta[200];   //[Electron_n]
   Float_t         Electron_phi[200];   //[Electron_n]
   Float_t         Electron_energy[200];   //[Electron_n]
   Float_t         Electron_charge[200];   //[Electron_n]
   Float_t         Electron_trkIso[200];   //[Electron_n]
   Float_t         Electron_ecalIso[200];   //[Electron_n]
   Float_t         Electron_hcalIso[200];   //[Electron_n]
   Float_t         Electron_SigmaIetaIeta[200];   //[Electron_n]
   Float_t         Electron_dEtaIn[200];   //[Electron_n]
   Float_t         Electron_dPhiIn[200];   //[Electron_n]
   Float_t         Electron_HoE[200];   //[Electron_n]
   Float_t         Electron_sc_energy[200];   //[Electron_n]
   Float_t         Electron_sc_eta[200];   //[Electron_n]
   Float_t         Electron_sc_phi[200];   //[Electron_n]
   Int_t           Photon_n;
   Float_t         Photon_E[200];   //[Photon_n]
   Float_t         Photon_pt[200];   //[Photon_n]
   Float_t         Photon_eta[200];   //[Photon_n]
   Float_t         Photon_phi[200];   //[Photon_n]
   Float_t         Photon_theta[200];   //[Photon_n]
   Float_t         Photon_et[200];   //[Photon_n]
   Float_t         Photon_swissCross[200];   //[Photon_n]
   Float_t         Photon_e6e2[200];   //[Photon_n]
   Float_t         Photon_e4e1[200];   //[Photon_n]
   Float_t         Photonr9[200];   //[Photon_n]
   Float_t         Photon_e1x5[200];   //[Photon_n]
   Float_t         Photon_e2x5[200];   //[Photon_n]
   Float_t         Photon_e3x3[200];   //[Photon_n]
   Float_t         Photon_e5x5[200];   //[Photon_n]
   Float_t         Photon_r1x5[200];   //[Photon_n]
   Float_t         Photon_r2x5[200];   //[Photon_n]
   Float_t         Photon_maxEnergyXtal[200];   //[Photon_n]
   Float_t         Photon_SigmaEtaEta[200];   //[Photon_n]
   Float_t         Photon_SigmaIetaIeta[200];   //[Photon_n]
   Float_t         Photon_SigmaEtaPhi[200];   //[Photon_n]
   Float_t         Photon_SigmaIetaIphi[200];   //[Photon_n]
   Float_t         Photon_SigmaPhiPhi[200];   //[Photon_n]
   Float_t         Photon_SigmaIphiIphi[200];   //[Photon_n]
   Float_t         Photon_Roundness[200];   //[Photon_n]
   Float_t         Photon_Angle[200];   //[Photon_n]
   Float_t         Photon_ecalRecHitSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_hcalTowerSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_trkSumPtSolidConeDR03[200];   //[Photon_n]
   Float_t         Photon_trkSumPtHollowConeDR03[200];   //[Photon_n]
   Int_t           Photon_nTrkSolidConeDR03[200];   //[Photon_n]
   Int_t           Photon_nTrkHollowConeDR03[200];   //[Photon_n]
   Float_t         Photon_hcalDepth1TowerSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_hcalDepth2TowerSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_ecalRecHitSumEtConeDR04[200];   //[Photon_n]
   Float_t         Photon_hcalTowerSumEtConeDR04[200];   //[Photon_n]
   Float_t         Photon_trkSumPtSolidConeDR04[200];   //[Photon_n]
   Float_t         Photon_trkSumPtHollowConeDR04[200];   //[Photon_n]
   Int_t           Photon_nTrkSolidConeDR04[200];   //[Photon_n]
   Int_t           Photon_nTrkHollowConeDR04[200];   //[Photon_n]
   Float_t         Photon_hcalDepth1TowerSumEtConeDR04[200];   //[Photon_n]
   Float_t         Photon_hcalDepth2TowerSumEtConeDR04[200];   //[Photon_n]
   Bool_t          Photon_hasPixelSeed[200];   //[Photon_n]
   Bool_t          Photon_isEB[200];   //[Photon_n]
   Bool_t          Photon_isEE[200];   //[Photon_n]
   Bool_t          Photon_isEBGap[200];   //[Photon_n]
   Bool_t          Photon_isEEGap[200];   //[Photon_n]
   Bool_t          Photon_isEBEEGap[200];   //[Photon_n]
   Float_t         Photon_e2e9[200];   //[Photon_n]
   Float_t         Photon_HoE[200];   //[Photon_n]
   Float_t         Photon_px[200];   //[Photon_n]
   Float_t         Photon_py[200];   //[Photon_n]
   Float_t         Photon_pz[200];   //[Photon_n]
   Float_t         Photon_vx[200];   //[Photon_n]
   Float_t         Photon_vy[200];   //[Photon_n]
   Float_t         Photon_vz[200];   //[Photon_n]
   Int_t           Photon_no_of_basic_clusters[200];   //[Photon_n]
   Float_t         Photon_sc_energy[200];   //[Photon_n]
   Float_t         Photon_sc_eta[200];   //[Photon_n]
   Float_t         Photon_sc_phi[200];   //[Photon_n]
   Float_t         Photon_sc_x[200];   //[Photon_n]
   Float_t         Photon_sc_y[200];   //[Photon_n]
   Float_t         Photon_sc_z[200];   //[Photon_n]
   Float_t         Photon_etaWidth[200];   //[Photon_n]
   Float_t         Photon_phiWidth[200];   //[Photon_n]
   Float_t         Photon_sc_et[200];   //[Photon_n]
   Float_t         matchphotonE[200];   //[Photon_n]
   Float_t         matchphotonpt[200];   //[Photon_n]
   Float_t         matchphotoneta[200];   //[Photon_n]
   Float_t         matchphotonphi[200];   //[Photon_n]
   Float_t         matchphotonpx[200];   //[Photon_n]
   Float_t         matchphotonpy[200];   //[Photon_n]
   Float_t         matchphotonpz[200];   //[Photon_n]
   Bool_t          ismatchedphoton[200];   //[Photon_n]
   Bool_t          Photon_hasConvTrk[200];   //[Photon_n]
   Int_t           Photon_ntracks[200];   //[Photon_n]
   Bool_t          Photon_isconverted[200];   //[Photon_n]
   Float_t         Photon_pairInvmass[200];   //[Photon_n]
   Float_t         Photon_pairCotThetaSeperation[200];   //[Photon_n]
   Float_t         Photon_pairmomentumX[200];   //[Photon_n]
   Float_t         Photon_pairmomentumY[200];   //[Photon_n]
   Float_t         Photon_pairmomentumZ[200];   //[Photon_n]
   Float_t         Photon_EoverP[200];   //[Photon_n]
   Float_t         Photon_ConvVx[200];   //[Photon_n]
   Float_t         Photon_ConvVy[200];   //[Photon_n]
   Float_t         Photon_ConvVz[200];   //[Photon_n]
   Float_t         Photon_ZOfPrimaryVertex[200];   //[Photon_n]
   Float_t         Photon_distOfMinimumApproach[200];   //[Photon_n]
   Float_t         Photon_dPhiTracksAtVtx[200];   //[Photon_n]
   Float_t         Photon_dPhiTracksAtEcal[200];   //[Photon_n]
   Float_t         Photon_dEtaTracksAtEcal[200];   //[Photon_n]
   Int_t           Photon_ncrys[200];   //[Photon_n]
   Float_t         Photon_timing_xtal[200][100];   //[Photon_n]
   Float_t         Photon_timingavg_xtal[200];   //[Photon_n]
   Float_t         Photon_energy_xtal[200][100];   //[Photon_n]
   Int_t           Photon_ieta_xtalEB[200][100];   //[Photon_n]
   Int_t           Photon_iphi_xtalEB[200][100];   //[Photon_n]
   Int_t           Photon_recoFlag_xtalEB[200][100];   //[Photon_n]
   Float_t         Photon_timeError_xtal[200][100];   //[Photon_n]
   Float_t         Photon_s9[200];   //[Photon_n]
   Bool_t          isBeamHaloGlobalLoosePass;
   Bool_t          isBeamHaloGlobalTightPass;
   Bool_t          isBeamHaloHcalLoosePass;
   Bool_t          isBeamHaloHcalTightPass;
   Bool_t          isBeamHaloCSCLoosePass;
   Bool_t          isBeamHaloCSCTightPass;
   Bool_t          isBeamHaloEcalLoosePass;
   Bool_t          isBeamHaloEcalTightPass;
   Bool_t          isBeamHaloIDTightPass;
   Bool_t          isBeamHaloIDLoosePass;
   Bool_t          isSmellsLikeHalo_Tag;
   Bool_t          isLooseHalo_Tag;
   Bool_t          isTightHalo_Tag;
   Bool_t          isExtremeTightHalo_Tag;
   Float_t         CaloMetSigma;
   Float_t         CaloMetEz;
   Float_t         CaloEtFractionHadronic;
   Float_t         CaloEmEtFraction;
   Float_t         CaloHadEtInHB;
   Float_t         CaloHadEtInHE;
   Float_t         CaloHadEtInHO;
   Float_t         CaloHadEtInHF;
   Float_t         CaloEmEtInEB;
   Float_t         CaloEmEtInEE;
   Float_t         CaloEmEtInHF;
   Float_t         CaloMaxEtInEmTowers;
   Float_t         CaloMaxEtInHadTowers;
   Float_t         CaloMetPt[6];
   Float_t         CaloMetPx[6];
   Float_t         CaloMetPy[6];
   Float_t         CaloMetPhi[6];
   Float_t         CaloMetSumEt[6];
   Float_t         Delta_phi;
   Float_t         PFMetPt[6];
   Float_t         PFMetPx[6];
   Float_t         PFMetPy[6];
   Float_t         PFMetPhi[6];
   Float_t         PFMetSumEt[6];
   Float_t         Delta_phiPF;

   // List of branches
   TBranch        *b_nevents;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_BXNumber;   //!
   TBranch        *b_totalIntensityBeam1;   //!
   TBranch        *b_totalIntensityBeam2;   //!
   TBranch        *b_avgInsDelLumi;   //!
   TBranch        *b_avgInsDelLumiErr;   //!
   TBranch        *b_avgInsRecLumi;   //!
   TBranch        *b_avgInsRecLumiErr;   //!
   TBranch        *b_ntriggers;   //!
   TBranch        *b_triggernames;   //!
   TBranch        *b_triggerprescales;   //!
   TBranch        *b_ifTriggerpassed;   //!
   TBranch        *b_trobjpt;   //!
   TBranch        *b_trobjeta;   //!
   TBranch        *b_trobjphi;   //!
   TBranch        *b_Vertex_n;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_tracksize;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_d0;   //!
   TBranch        *b_Vertex_isFake;   //!
   TBranch        *b_Scraping_isScrapingEvent;   //!
   TBranch        *b_Scraping_numOfTracks;   //!
   TBranch        *b_Scraping_fractionOfGoodTracks;   //!
   TBranch        *b_Jet_n;   //!
   TBranch        *b_Jet_px;   //!
   TBranch        *b_Jet_py;   //!
   TBranch        *b_Jet_E;   //!
   TBranch        *b_Jet_pz;   //!
   TBranch        *b_Jet_vx;   //!
   TBranch        *b_Jet_vy;   //!
   TBranch        *b_Jet_vz;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_emEnergyFraction;   //!
   TBranch        *b_Jet_energyFractionHadronic;   //!
   TBranch        *b_Jet_hitsInN90;   //!
   TBranch        *b_Jet_n90Hits;   //!
   TBranch        *b_Jet_nTowers;   //!
   TBranch        *b_Jet_fHPD;   //!
   TBranch        *b_Jet_fRBX;   //!
   TBranch        *b_Jet_RHF;   //!
   TBranch        *b_Jet_jecCorr;   //!
   TBranch        *b_ucJet_px;   //!
   TBranch        *b_ucJet_py;   //!
   TBranch        *b_ucJet_E;   //!
   TBranch        *b_ucJet_pz;   //!
   TBranch        *b_ucJet_pt;   //!
   TBranch        *b_ucJet_eta;   //!
   TBranch        *b_ucJet_phi;   //!
   TBranch        *b_pfJet_n;   //!
   TBranch        *b_pfJet_px;   //!
   TBranch        *b_pfJet_py;   //!
   TBranch        *b_pfJet_E;   //!
   TBranch        *b_pfJet_pz;   //!
   TBranch        *b_pfJet_vx;   //!
   TBranch        *b_pfJet_vy;   //!
   TBranch        *b_pfJet_vz;   //!
   TBranch        *b_pfJet_pt;   //!
   TBranch        *b_pfJet_eta;   //!
   TBranch        *b_pfJet_phi;   //!
   TBranch        *b_pfJet_jecCorr;   //!
   TBranch        *b_ucpfJet_px;   //!
   TBranch        *b_ucpfJet_py;   //!
   TBranch        *b_ucpfJet_E;   //!
   TBranch        *b_ucpfJet_pz;   //!
   TBranch        *b_ucpfJet_pt;   //!
   TBranch        *b_ucpfJet_eta;   //!
   TBranch        *b_ucpfJet_phi;   //!
   TBranch        *b_Electron_n;   //!
   TBranch        *b_Electron_px;   //!
   TBranch        *b_Electron_py;   //!
   TBranch        *b_Electron_pz;   //!
   TBranch        *b_Electron_vx;   //!
   TBranch        *b_Electron_vy;   //!
   TBranch        *b_Electron_vz;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_trkIso;   //!
   TBranch        *b_Electron_ecalIso;   //!
   TBranch        *b_Electron_hcalIso;   //!
   TBranch        *b_Electron_SigmaIetaIeta;   //!
   TBranch        *b_Electron_dEtaIn;   //!
   TBranch        *b_Electron_dPhiIn;   //!
   TBranch        *b_Electron_HoE;   //!
   TBranch        *b_Electron_sc_energy;   //!
   TBranch        *b_Electron_sc_eta;   //!
   TBranch        *b_Electron_sc_phi;   //!
   TBranch        *b_Photon_n;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_theta;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_swissCross;   //!
   TBranch        *b_Photon_e6e2;   //!
   TBranch        *b_Photon_e4e1;   //!
   TBranch        *b_Photonr9;   //!
   TBranch        *b_Photon_e1x5;   //!
   TBranch        *b_Photon_e2x5;   //!
   TBranch        *b_Photon_e3x3;   //!
   TBranch        *b_Photon_e5x5;   //!
   TBranch        *b_Photon_r1x5;   //!
   TBranch        *b_Photon_r2x5;   //!
   TBranch        *b_Photon_maxEnergyXtal;   //!
   TBranch        *b_Photon_SigmaEtaEta;   //!
   TBranch        *b_Photon_SigmaIetaIeta;   //!
   TBranch        *b_Photon_SigmaEtaPhi;   //!
   TBranch        *b_Photon_SigmaIetaIphi;   //!
   TBranch        *b_Photon_SigmaPhiPhi;   //!
   TBranch        *b_Photon_SigmaIphiIphi;   //!
   TBranch        *b_Photon_Roundness;   //!
   TBranch        *b_Photon_Angle;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR03;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_nTrkSolidConeDR03;   //!
   TBranch        *b_Photon_nTrkHollowConeDR03;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_Photon_nTrkSolidConeDR04;   //!
   TBranch        *b_Photon_nTrkHollowConeDR04;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!
   TBranch        *b_Photon_isEB;   //!
   TBranch        *b_Photon_isEE;   //!
   TBranch        *b_Photon_isEBGap;   //!
   TBranch        *b_Photon_isEEGap;   //!
   TBranch        *b_Photon_isEBEEGap;   //!
   TBranch        *b_Photon_e2e9;   //!
   TBranch        *b_Photon_HoE;   //!
   TBranch        *b_Photon_px;   //!
   TBranch        *b_Photon_py;   //!
   TBranch        *b_Photon_pz;   //!
   TBranch        *b_Photon_vx;   //!
   TBranch        *b_Photon_vy;   //!
   TBranch        *b_Photon_vz;   //!
   TBranch        *b_Photon_no_of_basic_clusters;   //!
   TBranch        *b_Photon_sc_energy;   //!
   TBranch        *b_Photon_sc_eta;   //!
   TBranch        *b_Photon_sc_phi;   //!
   TBranch        *b_Photon_sc_x;   //!
   TBranch        *b_Photon_sc_y;   //!
   TBranch        *b_Photon_sc_z;   //!
   TBranch        *b_Photon_etaWidth;   //!
   TBranch        *b_Photon_phiWidth;   //!
   TBranch        *b_Photon_sc_et;   //!
   TBranch        *b_matchphotonE;   //!
   TBranch        *b_matchphotonpt;   //!
   TBranch        *b_matchphotoneta;   //!
   TBranch        *b_matchphotonphi;   //!
   TBranch        *b_matchphotonpx;   //!
   TBranch        *b_matchphotonpy;   //!
   TBranch        *b_matchphotonpz;   //!
   TBranch        *b_ismatchedphoton;   //!
   TBranch        *b_Photon_hasConvTrk;   //!
   TBranch        *b_Photon_ntracks;   //!
   TBranch        *b_Photon_isconverted;   //!
   TBranch        *b_Photon_pairInvmass;   //!
   TBranch        *b_Photon_pairCotThetaSeperation;   //!
   TBranch        *b_Photon_pairmomentumX;   //!
   TBranch        *b_Photon_pairmomentumY;   //!
   TBranch        *b_Photon_pairmomentumZ;   //!
   TBranch        *b_Photon_EoverP;   //!
   TBranch        *b_Photon_ConvVx;   //!
   TBranch        *b_Photon_ConvVy;   //!
   TBranch        *b_Photon_ConvVz;   //!
   TBranch        *b_Photon_ZOfPrimaryVertex;   //!
   TBranch        *b_Photon_distOfMinimumApproach;   //!
   TBranch        *b_Photon_dPhiTracksAtVtx;   //!
   TBranch        *b_Photon_dPhiTracksAtEcal;   //!
   TBranch        *b_Photon_dEtaTracksAtEcal;   //!
   TBranch        *b_Photon_ncrys;   //!
   TBranch        *b_Photon_timing_xtal;   //!
   TBranch        *b_Photon_timingavg_xtal;   //!
   TBranch        *b_Photon_energy_xtal;   //!
   TBranch        *b_Photon_ieta_xtalEB;   //!
   TBranch        *b_Photon_iphi_xtalEB;   //!
   TBranch        *b_Photon_recoFlag_xtalEB;   //!
   TBranch        *b_Photon_timeError_xtal;   //!
   TBranch        *b_Photon_s9;   //!
   TBranch        *b_isBeamHaloGlobalLoosePass;   //!
   TBranch        *b_isBeamHaloGloablTightPass;   //!
   TBranch        *b_isBeamHaloHcalLoosePass;   //!
   TBranch        *b_isBeamHaloHcalTightPass;   //!
   TBranch        *b_isBeamHaloCSCLoosePass;   //!
   TBranch        *b_isBeamHaloCSCTightPass;   //!
   TBranch        *b_isBeamHaloEcalLoosePass;   //!
   TBranch        *b_isBeamHaloEcalTightPass;   //!
   TBranch        *b_isBeamHaloIDTightPass;   //!
   TBranch        *b_isBeamHaloIDLoosePass;   //!
   TBranch        *b_isSmellsLikeHalo_Tag;   //!
   TBranch        *b_isLooseHalo_Tag;   //!
   TBranch        *b_isTightHalo_Tag;   //!
   TBranch        *b_isExtremeTightHalo_Tag;   //!
   TBranch        *b_CaloMetSig;   //!
   TBranch        *b_CaloMetEz;   //!
   TBranch        *b_CaloEtFractionHadronic;   //!
   TBranch        *b_CaloEmEtFraction;   //!
   TBranch        *b_CaloHadEtInHB;   //!
   TBranch        *b_CaloHadEtInHE;   //!
   TBranch        *b_CaloHadEtInHO;   //!
   TBranch        *b_CaloHadEtInHF;   //!
   TBranch        *b_CaloEmEtInEB;   //!
   TBranch        *b_CaloEmEtInEE;   //!
   TBranch        *b_CaloEmEtInHF;   //!
   TBranch        *b_CaloMaxEtInEmTowers;   //!
   TBranch        *b_CaloMaxEtInHadTowers;   //!
   TBranch        *b_CaloMetPt;   //!
   TBranch        *b_CaloMetPx;   //!
   TBranch        *b_CaloMetPy;   //!
   TBranch        *b_CaloMetPhi;   //!
   TBranch        *b_CaloMetSumEt;   //!
   TBranch        *b_Delta_phi;   //!
   TBranch        *b_PFMetPt;   //!
   TBranch        *b_PFMetPx;   //!
   TBranch        *b_PFMetPy;   //!
   TBranch        *b_PFMetPhi;   //!
   TBranch        *b_PFMetSumEt;   //!
   TBranch        *b_Delta_phiPF;   //!

   myPlot();
   virtual ~myPlot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool isHLTMatch(std::string HLTPhotonTrigger[], int NTriggers, int ipho);
   virtual bool isTriggerPassed(std::string HLTPhotonTrigger[], int NTriggers, float &Weight, bool unprescale, bool usePSfactor);   
   virtual bool  CheckThisTrigger(std::string TriggerName); 

};

#endif



#ifdef myPlot_cxx   
myPlot::myPlot()    
{                   

TChain *chain = new TChain("myEvent");


TString inputFolder[1]={
                       "Data/",
                        };



for(int j=0;j<1;j++){

  TString location = "dcap://cmsgridftp.fnal.gov:24125";
  TString main_path = "/pnfs/cms/WAX/11/store/user/sushil/Summer11/QCDFakeRate_PhotonTrigger/";  //<<-----Change the in input File path
  TString main_short_path ="/pnfs/fnal.gov/usr/cms/WAX/11/store/user/sushil/Summer11/QCDFakeRate_PhotonTrigger/";

  TString subdirectory = inputFolder[j];

  TSystemDirectory sourceDir("hi",main_path+subdirectory);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter next(fileList);
  TSystemFile* fileName;

  int fileNumber = 1;
  int maxFiles   =-1;
  int fileCount=0;

  while ((fileName = (TSystemFile*)next()) && fileNumber >  maxFiles){
  if(TString(fileName->GetName()) == "." || TString(fileName->GetName()) == ".."  ){continue;}

   TString  FullPathInputFile = (location+main_short_path+subdirectory+fileName->GetName());
   cout<<j<<"   Root File Name: "<<subdirectory<<"   "<<(fileName->GetName())<<endl;
   chain->Add(FullPathInputFile);
   fileCount++;
   std::cout<<chain->GetEntries()<<"   For file nmber = "<<fileCount<<std::endl;
  }//loop over while

}//loop over j

  Init(chain);
}





myPlot::~myPlot()
{
   cout<<"destructor"<<endl;
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myPlot::GetEntry(Long64_t entry)
{

   cout<<"GetEntry"<<endl;
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


Long64_t myPlot::LoadTree(Long64_t entry)
{
   cout<<"LoadTree"<<endl;
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   cout<<"In the middle"<<endl;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myPlot::Init(TTree *tree)
{

    cout<<"Init"<<endl;
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggernames = 0;
   triggerprescales = 0;
   ifTriggerpassed = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nevents", &nevents, &b_nevents);
   fChain->SetBranchAddress("run", &run, &b_RunNumber);
   fChain->SetBranchAddress("event", &event, &b_EventNumber);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_LumiNumber);
   fChain->SetBranchAddress("beamCrossing", &beamCrossing, &b_BXNumber);
   fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1, &b_totalIntensityBeam1);
   fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2, &b_totalIntensityBeam2);
   fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi, &b_avgInsDelLumi);
   fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr, &b_avgInsDelLumiErr);
   fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_avgInsRecLumi);
   fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr, &b_avgInsRecLumiErr);
   fChain->SetBranchAddress("ntriggers", &ntriggers, &b_ntriggers);
   fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
   fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
   fChain->SetBranchAddress("ifTriggerpassed", &ifTriggerpassed, &b_ifTriggerpassed);
   fChain->SetBranchAddress("trobjpt", trobjpt, &b_trobjpt);
   fChain->SetBranchAddress("trobjeta", trobjeta, &b_trobjeta);
   fChain->SetBranchAddress("trobjphi", trobjphi, &b_trobjphi);
   fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
   fChain->SetBranchAddress("Vertex_x", Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex_tracksize", Vertex_tracksize, &b_Vertex_tracksize);
   fChain->SetBranchAddress("Vertex_ndof", Vertex_ndof, &b_Vertex_ndof);
   fChain->SetBranchAddress("Vertex_chi2", Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_d0", Vertex_d0, &b_Vertex_d0);
   fChain->SetBranchAddress("Vertex_isFake", Vertex_isFake, &b_Vertex_isFake);
   fChain->SetBranchAddress("Scraping_isScrapingEvent", &Scraping_isScrapingEvent, &b_Scraping_isScrapingEvent);
   fChain->SetBranchAddress("Scraping_numOfTracks", &Scraping_numOfTracks, &b_Scraping_numOfTracks);
   fChain->SetBranchAddress("Scraping_fractionOfGoodTracks", &Scraping_fractionOfGoodTracks, &b_Scraping_fractionOfGoodTracks);
   fChain->SetBranchAddress("Jet_n", &Jet_n, &b_Jet_n);
   fChain->SetBranchAddress("Jet_px", Jet_px, &b_Jet_px);
   fChain->SetBranchAddress("Jet_py", Jet_py, &b_Jet_py);
   fChain->SetBranchAddress("Jet_E", Jet_E, &b_Jet_E);
   fChain->SetBranchAddress("Jet_pz", Jet_pz, &b_Jet_pz);
   fChain->SetBranchAddress("Jet_vx", Jet_vx, &b_Jet_vx);
   fChain->SetBranchAddress("Jet_vy", Jet_vy, &b_Jet_vy);
   fChain->SetBranchAddress("Jet_vz", Jet_vz, &b_Jet_vz);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_emEnergyFraction", Jet_emEnergyFraction, &b_Jet_emEnergyFraction);
   fChain->SetBranchAddress("Jet_energyFractionHadronic", Jet_energyFractionHadronic, &b_Jet_energyFractionHadronic);
   fChain->SetBranchAddress("Jet_hitsInN90", Jet_hitsInN90, &b_Jet_hitsInN90);
   fChain->SetBranchAddress("Jet_n90Hits", Jet_n90Hits, &b_Jet_n90Hits);
   fChain->SetBranchAddress("Jet_nTowers", Jet_nTowers, &b_Jet_nTowers);
   fChain->SetBranchAddress("Jet_fHPD", Jet_fHPD, &b_Jet_fHPD);
   fChain->SetBranchAddress("Jet_fRBX", Jet_fRBX, &b_Jet_fRBX);
   fChain->SetBranchAddress("Jet_RHF", Jet_RHF, &b_Jet_RHF);
   fChain->SetBranchAddress("Jet_jecCorr", Jet_jecCorr, &b_Jet_jecCorr);
   fChain->SetBranchAddress("ucJet_px", ucJet_px, &b_ucJet_px);
   fChain->SetBranchAddress("ucJet_py", ucJet_py, &b_ucJet_py);
   fChain->SetBranchAddress("ucJet_E", ucJet_E, &b_ucJet_E);
   fChain->SetBranchAddress("ucJet_pz", ucJet_pz, &b_ucJet_pz);
   fChain->SetBranchAddress("ucJet_pt", ucJet_pt, &b_ucJet_pt);
   fChain->SetBranchAddress("ucJet_eta", ucJet_eta, &b_ucJet_eta);
   fChain->SetBranchAddress("ucJet_phi", ucJet_phi, &b_ucJet_phi);
   fChain->SetBranchAddress("pfJet_n", &pfJet_n, &b_pfJet_n);
   fChain->SetBranchAddress("pfJet_px", pfJet_px, &b_pfJet_px);
   fChain->SetBranchAddress("pfJet_py", pfJet_py, &b_pfJet_py);
   fChain->SetBranchAddress("pfJet_E", pfJet_E, &b_pfJet_E);
   fChain->SetBranchAddress("pfJet_pz", pfJet_pz, &b_pfJet_pz);
   fChain->SetBranchAddress("pfJet_vx", pfJet_vx, &b_pfJet_vx);
   fChain->SetBranchAddress("pfJet_vy", pfJet_vy, &b_pfJet_vy);
   fChain->SetBranchAddress("pfJet_vz", pfJet_vz, &b_pfJet_vz);
   fChain->SetBranchAddress("pfJet_pt", pfJet_pt, &b_pfJet_pt);
   fChain->SetBranchAddress("pfJet_eta", pfJet_eta, &b_pfJet_eta);
   fChain->SetBranchAddress("pfJet_phi", pfJet_phi, &b_pfJet_phi);
   fChain->SetBranchAddress("pfJet_jecCorr", pfJet_jecCorr, &b_pfJet_jecCorr);
   fChain->SetBranchAddress("ucpfJet_px", ucpfJet_px, &b_ucpfJet_px);
   fChain->SetBranchAddress("ucpfJet_py", ucpfJet_py, &b_ucpfJet_py);
   fChain->SetBranchAddress("ucpfJet_E", ucpfJet_E, &b_ucpfJet_E);
   fChain->SetBranchAddress("ucpfJet_pz", ucpfJet_pz, &b_ucpfJet_pz);
   fChain->SetBranchAddress("ucpfJet_pt", ucpfJet_pt, &b_ucpfJet_pt);
   fChain->SetBranchAddress("ucpfJet_eta", ucpfJet_eta, &b_ucpfJet_eta);
   fChain->SetBranchAddress("ucpfJet_phi", ucpfJet_phi, &b_ucpfJet_phi);
   fChain->SetBranchAddress("Electron_n", &Electron_n, &b_Electron_n);
   fChain->SetBranchAddress("Electron_px", Electron_px, &b_Electron_px);
   fChain->SetBranchAddress("Electron_py", Electron_py, &b_Electron_py);
   fChain->SetBranchAddress("Electron_pz", Electron_pz, &b_Electron_pz);
   fChain->SetBranchAddress("Electron_vx", Electron_vx, &b_Electron_vx);
   fChain->SetBranchAddress("Electron_vy", Electron_vy, &b_Electron_vy);
   fChain->SetBranchAddress("Electron_vz", Electron_vz, &b_Electron_vz);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_energy", Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_trkIso", Electron_trkIso, &b_Electron_trkIso);
   fChain->SetBranchAddress("Electron_ecalIso", Electron_ecalIso, &b_Electron_ecalIso);
   fChain->SetBranchAddress("Electron_hcalIso", Electron_hcalIso, &b_Electron_hcalIso);
   fChain->SetBranchAddress("Electron_SigmaIetaIeta", Electron_SigmaIetaIeta, &b_Electron_SigmaIetaIeta);
   fChain->SetBranchAddress("Electron_dEtaIn", Electron_dEtaIn, &b_Electron_dEtaIn);
   fChain->SetBranchAddress("Electron_dPhiIn", Electron_dPhiIn, &b_Electron_dPhiIn);
   fChain->SetBranchAddress("Electron_HoE", Electron_HoE, &b_Electron_HoE);
   fChain->SetBranchAddress("Electron_sc_energy", Electron_sc_energy, &b_Electron_sc_energy);
   fChain->SetBranchAddress("Electron_sc_eta", Electron_sc_eta, &b_Electron_sc_eta);
   fChain->SetBranchAddress("Electron_sc_phi", Electron_sc_phi, &b_Electron_sc_phi);
   fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
   fChain->SetBranchAddress("Photon_E", Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_theta", Photon_theta, &b_Photon_theta);
   fChain->SetBranchAddress("Photon_et", Photon_et, &b_Photon_et);
   fChain->SetBranchAddress("Photon_swissCross", Photon_swissCross, &b_Photon_swissCross);
   fChain->SetBranchAddress("Photon_e6e2", Photon_e6e2, &b_Photon_e6e2);
   fChain->SetBranchAddress("Photon_e4e1", Photon_e4e1, &b_Photon_e4e1);
   fChain->SetBranchAddress("Photonr9", Photonr9, &b_Photonr9);
   fChain->SetBranchAddress("Photon_e1x5", Photon_e1x5, &b_Photon_e1x5);
   fChain->SetBranchAddress("Photon_e2x5", Photon_e2x5, &b_Photon_e2x5);
   fChain->SetBranchAddress("Photon_e3x3", Photon_e3x3, &b_Photon_e3x3);
   fChain->SetBranchAddress("Photon_e5x5", Photon_e5x5, &b_Photon_e5x5);
   fChain->SetBranchAddress("Photon_r1x5", Photon_r1x5, &b_Photon_r1x5);
   fChain->SetBranchAddress("Photon_r2x5", Photon_r2x5, &b_Photon_r2x5);
   fChain->SetBranchAddress("Photon_maxEnergyXtal", Photon_maxEnergyXtal, &b_Photon_maxEnergyXtal);
   fChain->SetBranchAddress("Photon_SigmaEtaEta", Photon_SigmaEtaEta, &b_Photon_SigmaEtaEta);
   fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta, &b_Photon_SigmaIetaIeta);
   fChain->SetBranchAddress("Photon_SigmaEtaPhi", Photon_SigmaEtaPhi, &b_Photon_SigmaEtaPhi);
   fChain->SetBranchAddress("Photon_SigmaIetaIphi", Photon_SigmaIetaIphi, &b_Photon_SigmaIetaIphi);
   fChain->SetBranchAddress("Photon_SigmaPhiPhi", Photon_SigmaPhiPhi, &b_Photon_SigmaPhiPhi);
   fChain->SetBranchAddress("Photon_SigmaIphiIphi", Photon_SigmaIphiIphi, &b_Photon_SigmaIphiIphi);
   fChain->SetBranchAddress("Photon_Roundness", Photon_Roundness, &b_Photon_Roundness);
   fChain->SetBranchAddress("Photon_Angle", Photon_Angle, &b_Photon_Angle);
   fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", Photon_ecalRecHitSumEtConeDR03, &b_Photon_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", Photon_hcalTowerSumEtConeDR03, &b_Photon_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", Photon_trkSumPtSolidConeDR03, &b_Photon_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", Photon_nTrkSolidConeDR03, &b_Photon_nTrkSolidConeDR03);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", Photon_nTrkHollowConeDR03, &b_Photon_nTrkHollowConeDR03);
   fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", Photon_hcalDepth1TowerSumEtConeDR03, &b_Photon_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", Photon_hcalDepth2TowerSumEtConeDR03, &b_Photon_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", Photon_ecalRecHitSumEtConeDR04, &b_Photon_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
   fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", Photon_hcalDepth1TowerSumEtConeDR04, &b_Photon_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", Photon_hcalDepth2TowerSumEtConeDR04, &b_Photon_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   fChain->SetBranchAddress("Photon_isEB", Photon_isEB, &b_Photon_isEB);
   fChain->SetBranchAddress("Photon_isEE", Photon_isEE, &b_Photon_isEE);
   fChain->SetBranchAddress("Photon_isEBGap", Photon_isEBGap, &b_Photon_isEBGap);
   fChain->SetBranchAddress("Photon_isEEGap", Photon_isEEGap, &b_Photon_isEEGap);
   fChain->SetBranchAddress("Photon_isEBEEGap", Photon_isEBEEGap, &b_Photon_isEBEEGap);
   fChain->SetBranchAddress("Photon_e2e9", Photon_e2e9, &b_Photon_e2e9);
   fChain->SetBranchAddress("Photon_HoE", Photon_HoE, &b_Photon_HoE);
   fChain->SetBranchAddress("Photon_px", Photon_px, &b_Photon_px);
   fChain->SetBranchAddress("Photon_py", Photon_py, &b_Photon_py);
   fChain->SetBranchAddress("Photon_pz", Photon_pz, &b_Photon_pz);
   fChain->SetBranchAddress("Photon_vx", Photon_vx, &b_Photon_vx);
   fChain->SetBranchAddress("Photon_vy", Photon_vy, &b_Photon_vy);
   fChain->SetBranchAddress("Photon_vz", Photon_vz, &b_Photon_vz);
   fChain->SetBranchAddress("Photon_no_of_basic_clusters", Photon_no_of_basic_clusters, &b_Photon_no_of_basic_clusters);
   fChain->SetBranchAddress("Photon_sc_energy", Photon_sc_energy, &b_Photon_sc_energy);
   fChain->SetBranchAddress("Photon_sc_eta", Photon_sc_eta, &b_Photon_sc_eta);
   fChain->SetBranchAddress("Photon_sc_phi", Photon_sc_phi, &b_Photon_sc_phi);
   fChain->SetBranchAddress("Photon_sc_x", Photon_sc_x, &b_Photon_sc_x);
   fChain->SetBranchAddress("Photon_sc_y", Photon_sc_y, &b_Photon_sc_y);
   fChain->SetBranchAddress("Photon_sc_z", Photon_sc_z, &b_Photon_sc_z);
   fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
   fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
   fChain->SetBranchAddress("Photon_sc_et", Photon_sc_et, &b_Photon_sc_et);
   fChain->SetBranchAddress("matchphotonE", matchphotonE, &b_matchphotonE);
   fChain->SetBranchAddress("matchphotonpt", matchphotonpt, &b_matchphotonpt);
   fChain->SetBranchAddress("matchphotoneta", matchphotoneta, &b_matchphotoneta);
   fChain->SetBranchAddress("matchphotonphi", matchphotonphi, &b_matchphotonphi);
   fChain->SetBranchAddress("matchphotonpx", matchphotonpx, &b_matchphotonpx);
   fChain->SetBranchAddress("matchphotonpy", matchphotonpy, &b_matchphotonpy);
   fChain->SetBranchAddress("matchphotonpz", matchphotonpz, &b_matchphotonpz);
   fChain->SetBranchAddress("ismatchedphoton", ismatchedphoton, &b_ismatchedphoton);
   fChain->SetBranchAddress("Photon_hasConvTrk", Photon_hasConvTrk, &b_Photon_hasConvTrk);
   fChain->SetBranchAddress("Photon_ntracks", Photon_ntracks, &b_Photon_ntracks);
   fChain->SetBranchAddress("Photon_isconverted", Photon_isconverted, &b_Photon_isconverted);
   fChain->SetBranchAddress("Photon_pairInvmass", Photon_pairInvmass, &b_Photon_pairInvmass);
   fChain->SetBranchAddress("Photon_pairCotThetaSeperation", Photon_pairCotThetaSeperation, &b_Photon_pairCotThetaSeperation);
   fChain->SetBranchAddress("Photon_pairmomentumX", Photon_pairmomentumX, &b_Photon_pairmomentumX);
   fChain->SetBranchAddress("Photon_pairmomentumY", Photon_pairmomentumY, &b_Photon_pairmomentumY);
   fChain->SetBranchAddress("Photon_pairmomentumZ", Photon_pairmomentumZ, &b_Photon_pairmomentumZ);
   fChain->SetBranchAddress("Photon_EoverP", Photon_EoverP, &b_Photon_EoverP);
   fChain->SetBranchAddress("Photon_ConvVx", Photon_ConvVx, &b_Photon_ConvVx);
   fChain->SetBranchAddress("Photon_ConvVy", Photon_ConvVy, &b_Photon_ConvVy);
   fChain->SetBranchAddress("Photon_ConvVz", Photon_ConvVz, &b_Photon_ConvVz);
   fChain->SetBranchAddress("Photon_ZOfPrimaryVertex", Photon_ZOfPrimaryVertex, &b_Photon_ZOfPrimaryVertex);
   fChain->SetBranchAddress("Photon_distOfMinimumApproach", Photon_distOfMinimumApproach, &b_Photon_distOfMinimumApproach);
   fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", Photon_dPhiTracksAtVtx, &b_Photon_dPhiTracksAtVtx);
   fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", Photon_dPhiTracksAtEcal, &b_Photon_dPhiTracksAtEcal);
   fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", Photon_dEtaTracksAtEcal, &b_Photon_dEtaTracksAtEcal);
   fChain->SetBranchAddress("Photon_ncrys", Photon_ncrys, &b_Photon_ncrys);
   fChain->SetBranchAddress("Photon_timing_xtal", Photon_timing_xtal, &b_Photon_timing_xtal);
   fChain->SetBranchAddress("Photon_timingavg_xtal", Photon_timingavg_xtal, &b_Photon_timingavg_xtal);
   fChain->SetBranchAddress("Photon_energy_xtal", Photon_energy_xtal, &b_Photon_energy_xtal);
   fChain->SetBranchAddress("Photon_ieta_xtalEB", Photon_ieta_xtalEB, &b_Photon_ieta_xtalEB);
   fChain->SetBranchAddress("Photon_iphi_xtalEB", Photon_iphi_xtalEB, &b_Photon_iphi_xtalEB);
   fChain->SetBranchAddress("Photon_recoFlag_xtalEB", Photon_recoFlag_xtalEB, &b_Photon_recoFlag_xtalEB);
   fChain->SetBranchAddress("Photon_timeError_xtal", Photon_timeError_xtal, &b_Photon_timeError_xtal);
   fChain->SetBranchAddress("Photon_s9", Photon_s9, &b_Photon_s9);
   fChain->SetBranchAddress("isBeamHaloGlobalLoosePass", &isBeamHaloGlobalLoosePass, &b_isBeamHaloGlobalLoosePass);
   fChain->SetBranchAddress("isBeamHaloGlobalTightPass", &isBeamHaloGlobalTightPass, &b_isBeamHaloGloablTightPass);
   fChain->SetBranchAddress("isBeamHaloHcalLoosePass", &isBeamHaloHcalLoosePass, &b_isBeamHaloHcalLoosePass);
   fChain->SetBranchAddress("isBeamHaloHcalTightPass", &isBeamHaloHcalTightPass, &b_isBeamHaloHcalTightPass);
   fChain->SetBranchAddress("isBeamHaloCSCLoosePass", &isBeamHaloCSCLoosePass, &b_isBeamHaloCSCLoosePass);
   fChain->SetBranchAddress("isBeamHaloCSCTightPass", &isBeamHaloCSCTightPass, &b_isBeamHaloCSCTightPass);
   fChain->SetBranchAddress("isBeamHaloEcalLoosePass", &isBeamHaloEcalLoosePass, &b_isBeamHaloEcalLoosePass);
   fChain->SetBranchAddress("isBeamHaloEcalTightPass", &isBeamHaloEcalTightPass, &b_isBeamHaloEcalTightPass);
   fChain->SetBranchAddress("isBeamHaloIDTightPass", &isBeamHaloIDTightPass, &b_isBeamHaloIDTightPass);
   fChain->SetBranchAddress("isBeamHaloIDLoosePass", &isBeamHaloIDLoosePass, &b_isBeamHaloIDLoosePass);
   fChain->SetBranchAddress("isSmellsLikeHalo_Tag", &isSmellsLikeHalo_Tag, &b_isSmellsLikeHalo_Tag);
   fChain->SetBranchAddress("isLooseHalo_Tag", &isLooseHalo_Tag, &b_isLooseHalo_Tag);
   fChain->SetBranchAddress("isTightHalo_Tag", &isTightHalo_Tag, &b_isTightHalo_Tag);
   fChain->SetBranchAddress("isExtremeTightHalo_Tag", &isExtremeTightHalo_Tag, &b_isExtremeTightHalo_Tag);
   fChain->SetBranchAddress("CaloMetSigma", &CaloMetSigma, &b_CaloMetSig);
   fChain->SetBranchAddress("CaloMetEz", &CaloMetEz, &b_CaloMetEz);
   fChain->SetBranchAddress("CaloEtFractionHadronic", &CaloEtFractionHadronic, &b_CaloEtFractionHadronic);
   fChain->SetBranchAddress("CaloEmEtFraction", &CaloEmEtFraction, &b_CaloEmEtFraction);
   fChain->SetBranchAddress("CaloHadEtInHB", &CaloHadEtInHB, &b_CaloHadEtInHB);
   fChain->SetBranchAddress("CaloHadEtInHE", &CaloHadEtInHE, &b_CaloHadEtInHE);
   fChain->SetBranchAddress("CaloHadEtInHO", &CaloHadEtInHO, &b_CaloHadEtInHO);
   fChain->SetBranchAddress("CaloHadEtInHF", &CaloHadEtInHF, &b_CaloHadEtInHF);
   fChain->SetBranchAddress("CaloEmEtInEB", &CaloEmEtInEB, &b_CaloEmEtInEB);
   fChain->SetBranchAddress("CaloEmEtInEE", &CaloEmEtInEE, &b_CaloEmEtInEE);
   fChain->SetBranchAddress("CaloEmEtInHF", &CaloEmEtInHF, &b_CaloEmEtInHF);
   fChain->SetBranchAddress("CaloMaxEtInEmTowers", &CaloMaxEtInEmTowers, &b_CaloMaxEtInEmTowers);
   fChain->SetBranchAddress("CaloMaxEtInHadTowers", &CaloMaxEtInHadTowers, &b_CaloMaxEtInHadTowers);
   fChain->SetBranchAddress("CaloMetPt", CaloMetPt, &b_CaloMetPt);
   fChain->SetBranchAddress("CaloMetPx", CaloMetPx, &b_CaloMetPx);
   fChain->SetBranchAddress("CaloMetPy", CaloMetPy, &b_CaloMetPy);
   fChain->SetBranchAddress("CaloMetPhi", CaloMetPhi, &b_CaloMetPhi);
   fChain->SetBranchAddress("CaloMetSumEt", CaloMetSumEt, &b_CaloMetSumEt);
   fChain->SetBranchAddress("Delta_phi", &Delta_phi, &b_Delta_phi);
   fChain->SetBranchAddress("PFMetPt", PFMetPt, &b_PFMetPt);
   fChain->SetBranchAddress("PFMetPx", PFMetPx, &b_PFMetPx);
   fChain->SetBranchAddress("PFMetPy", PFMetPy, &b_PFMetPy);
   fChain->SetBranchAddress("PFMetPhi", PFMetPhi, &b_PFMetPhi);
   fChain->SetBranchAddress("PFMetSumEt", PFMetSumEt, &b_PFMetSumEt);
   fChain->SetBranchAddress("Delta_phiPF", &Delta_phiPF, &b_Delta_phiPF);
   Notify();

  cout<<" Chain Init Is Done "<<endl;
}

Bool_t myPlot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   cout<<"Insdie a Notify "<<endl;
   return kTRUE;
}

void myPlot::Show(Long64_t entry)
{
   cout<<"Show"<<endl;
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myPlot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


bool myPlot::isTriggerPassed(std::string HLTPhotonTrigger[],
                             int NTriggers,
                              float & Weight, 
                               bool unprescale,
                                 bool usePSfactor){
 
   bool HLTPassed=false;
   int  HLTPreScale=1;
        Weight=-1;

 
   for(int p=0; p < (*triggernames).size();p++){
 
     for(int k=0; k < NTriggers; k++){
      
         string string_search (HLTPhotonTrigger[k]);
         size_t found = (*triggernames)[p].find(string_search);     
                                                                                           
         if(found != string::npos){
         //if( HLTPhotonTrigger[k] == (*triggernames)[p])
                     HLTPreScale =  (*triggerprescales)[p];
                     HLTPassed   = (*ifTriggerpassed)[p];
           //cout<<"Trigger = "<<((*triggernames)[p])<<"  Passed= "<<HLTPassed<<"  Pre-Scale= "<<HLTPreScale<<endl; 
                    if( HLTPassed){

                          if(!unprescale && usePSfactor)Weight=HLTPreScale;
                          if(!unprescale && !usePSfactor)Weight=1.0;
                          if(unprescale && HLTPreScale==1)Weight=1.0;

                          if(unprescale && HLTPreScale!=1.0)
                                { Weight=-1.0;
                                  HLTPassed=false;
                                }
                       }
 
              }
        }
  }
  

  return HLTPassed;
 
}                           


bool myPlot::isHLTMatch(std::string HLTPhotonTrigger[],
                             int NTriggers,
                              int ipho){

   bool thisPhotonOK=true;

  for(int p=0; p < (*triggernames).size();p++){
          
     for(int k=0; k < NTriggers;k++){

         string string_search (HLTPhotonTrigger[k]);
         size_t found = (*triggernames)[p].find(string_search);     
                                                                                           
         if(found != string::npos){

                   for(int f=0; f<100; f++)
                       {
                           for(int key=0; key< 10; key++)
                            { if(trobjpt[p][f][key]>0.)
                               {
                               double  dphi = Photon_phi[ipho]- trobjphi[p][f][key];
                               double  deta = Photon_eta[ipho]- trobjeta[p][f][key];
                               double  delta_r = pow( ((dphi*dphi)+ (deta*deta))  ,0.5);

                                if(delta_r < 0.2) thisPhotonOK = false;
                                }
                            }
                       }
    
              }//check Trigger exist or not
   }//k=0;
  }//p=0;

  //cout<<"keep this photon = "<<thisPhotonOK<<endl;
  return thisPhotonOK;
}

bool myPlot::CheckThisTrigger(std::string TriggerName){

 bool trigexist=false;
 bool HLTPassed=false;

 for(int p=0; p < (*triggernames).size();p++){
 
         string string_search (TriggerName);
         size_t found = (*triggernames)[p].find(string_search);     
                                                                                           
         if(found != string::npos){
                     HLTPassed   = (*ifTriggerpassed)[p];

                    if( HLTPassed)trigexist=true;

                   }//fount the trigger

        }//for loop

  return trigexist;

}


#endif // #ifdef myPlot_cxx
