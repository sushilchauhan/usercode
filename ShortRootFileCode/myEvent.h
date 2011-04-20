//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan  5 06:16:02 2011 by ROOT version 5.22/00d
// from TTree myEvent/a tree with histograms
// found on file: dcache:///pnfs/cms/WAX/11/store/user/sushil/MonoPhoton/387_Ntuples_V26/WENu_Gamma/Histo_Wenu_Gamma_10_1_fNU.root
//////////////////////////////////////////////////////////

#ifndef myEvent_h
#define myEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <list>
#include <vector>
#include <TCanvas.h>
#include <TSystem.h>
#include <TPostScript.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TRef.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
using namespace std;
#ifdef __MAKECINT__
#pragma link C++ class vector<std::bool>;
#pragma link C++ class vector<std::string>+;
#pragma link C++ class vector<int>+;
#pragma extra_include "vector";
#endif

class myEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   Int_t           Vertex_n;
   Float_t         Vertex_x[100];   //[Vertex_n]
   Float_t         Vertex_y[100];   //[Vertex_n]
   Float_t         Vertex_z[100];   //[Vertex_n]
   Int_t           Vertex_tracksize[100];   //[Vertex_n]
   Int_t           Vertex_ndof[100];   //[Vertex_n]
   Float_t         Vertex_chi2[100];   //[Vertex_n]
   Float_t         Vertex_d0[100];   //[Vertex_n]
   Bool_t          Vertex_isFake[100];   //[Vertex_n]
   Int_t           Track_n;
   Float_t         Track_px[100];   //[Track_n]
   Float_t         Track_py[100];   //[Track_n]
   Float_t         Track_pz[100];   //[Track_n]
   Float_t         Track_vx[100];   //[Track_n]
   Float_t         Track_vy[100];   //[Track_n]
   Float_t         Track_vz[100];   //[Track_n]
   Float_t         Track_pt[100];   //[Track_n]
   Float_t         Track_eta[100];   //[Track_n]
   Float_t         Track_phi[100];   //[Track_n]
   Int_t           Jet_n;
   Float_t         Jet_px[100];   //[Jet_n]
   Float_t         Jet_py[100];   //[Jet_n]
   Float_t         Jet_E[100];   //[Jet_n]
   Float_t         Jet_pz[100];   //[Jet_n]
   Float_t         Jet_vx[100];   //[Jet_n]
   Float_t         Jet_vy[100];   //[Jet_n]
   Float_t         Jet_vz[100];   //[Jet_n]
   Float_t         Jet_pt[100];   //[Jet_n]
   Float_t         Jet_eta[100];   //[Jet_n]
   Float_t         Jet_phi[100];   //[Jet_n]
   Float_t         Jet_emEnergyFraction[100];   //[Jet_n]
   Float_t         Jet_energyFractionHadronic[100];   //[Jet_n]
   Int_t           Jet_hitsInN90[100];   //[Jet_n]
   Int_t           Jet_n90Hits[100];   //[Jet_n]
   Int_t           Jet_nTowers[100];   //[Jet_n]
   Float_t         Jet_fHPD[100];   //[Jet_n]
   Float_t         Jet_fRBX[100];   //[Jet_n]
   Float_t         Jet_RHF[100];   //[Jet_n]

   Int_t            pfJet_n;
   Float_t          pfJet_pt[100];                                                                                                                                
   Float_t          pfJet_px[100];
   Float_t          pfJet_py[100];
   Float_t          pfJet_E[100];
   Float_t          pfJet_pz[100];
   Float_t          pfJet_vx[100];
   Float_t          pfJet_vy[100];
   Float_t          pfJet_vz[100];
   Float_t          pfJet_eta[100];
   Float_t          pfJet_phi[100];


   Int_t           Electron_n;
   Float_t         Electron_px[100];   //[Electron_n]
   Float_t         Electron_py[100];   //[Electron_n]
   Float_t         Electron_pz[100];   //[Electron_n]
   Float_t         Electron_vx[100];   //[Electron_n]
   Float_t         Electron_vy[100];   //[Electron_n]
   Float_t         Electron_vz[100];   //[Electron_n]
   Float_t         Electron_pt[100];   //[Electron_n]
   Float_t         Electron_eta[100];   //[Electron_n]
   Float_t         Electron_phi[100];   //[Electron_n]
   Float_t         Electron_energy[100];   //[Electron_n]
   Float_t         Electron_charge[100];   //[Electron_n]
   Float_t         Electron_trkIso[100];   //[Electron_n]
   Float_t         Electron_ecalIso[100];   //[Electron_n]
   Float_t         Electron_hcalIso[100];   //[Electron_n]
   Float_t         Electron_SigmaIetaIeta[100];   //[Electron_n]
   Float_t         Electron_dEtaIn[100];   //[Electron_n]
   Float_t         Electron_dPhiIn[100];   //[Electron_n]
   Float_t         Electron_HoE[100];   //[Electron_n]
   Float_t         Electron_sc_energy[100];   //[Electron_n]
   Float_t         Electron_sc_eta[100];   //[Electron_n]
   Float_t         Electron_sc_phi[100];   //[Electron_n]
   Int_t           Muon_n;
   Float_t         Muon_px[100];   //[Muon_n]
   Float_t         Muon_py[100];   //[Muon_n]
   Float_t         Muon_pz[100];   //[Muon_n]
   Float_t         Muon_vx[100];   //[Muon_n]
   Float_t         Muon_vy[100];   //[Muon_n]
   Float_t         Muon_vz[100];   //[Muon_n]
   Float_t         Muon_pt[100];   //[Muon_n]
   Float_t         Muon_eta[100];   //[Muon_n]
   Float_t         Muon_phi[100];   //[Muon_n]
   Float_t         Muon_energy[100];   //[Muon_n]
   Float_t         Muon_charge[100];   //[Muon_n]
   Bool_t          Muon_isGlobalMuon[100];   //[Muon_n]
   Bool_t          Muon_isTrackerMuon[100];   //[Muon_n]
   Bool_t          Muon_isStandAloneMuon[100];   //[Muon_n]
   Bool_t          Muon_InnerTrack_isNonnull[100];   //[Muon_n]
   Bool_t          Muon_OuterTrack_isNonnull[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_x[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_y[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_z[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_px[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_py[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_pz[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_x[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_y[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_z[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_px[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_py[100];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_pz[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_x[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_y[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_z[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_px[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_py[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_pz[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_x[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_y[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_z[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_px[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_py[100];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_pz[100];   //[Muon_n]
   Int_t           CosmicMuon_n;
   Float_t         CosmicMuon_px[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_py[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_pz[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_pt[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_eta[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_phi[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_energy[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_charge[100];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_isGlobalMuon[100];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_isTrackerMuon[100];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_isStandAloneMuon[100];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_InnerTrack_isNonnull[100];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_OuterTrack_isNonnull[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_x[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_y[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_z[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_px[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_py[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_pz[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_x[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_y[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_z[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_px[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_py[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_pz[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_x[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_y[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_z[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_px[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_py[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_pz[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_x[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_y[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_z[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_px[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_py[100];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_pz[100];   //[CosmicMuon_n]
   Int_t           Tau_n;
   Float_t         Tau_px[100];   //[Tau_n]
   Float_t         Tau_py[100];   //[Tau_n]
   Float_t         Tau_pz[100];   //[Tau_n]
   Float_t         Tau_vx[100];   //[Tau_n]
   Float_t         Tau_vy[100];   //[Tau_n]
   Float_t         Tau_vz[100];   //[Tau_n]
   Float_t         Tau_pt[100];   //[Tau_n]
   Float_t         Tau_eta[100];   //[Tau_n]
   Float_t         Tau_phi[100];   //[Tau_n]
   Float_t         Tau_energy[100];   //[Tau_n]
   Float_t         Tau_charge[100];   //[Tau_n]
   Float_t         gen_pthat;
   Int_t           ngenphotons;
   Float_t         gen_photonpt[1000];   //[ngenphotons]
   Float_t         gen_photoneta[1000];   //[ngenphotons]
   Float_t         gen_photonphi[1000];   //[ngenphotons]
   Float_t         gen_photonpx[1000];   //[ngenphotons]
   Float_t         gen_photonpy[1000];   //[ngenphotons]
   Float_t         gen_photonpz[1000];   //[ngenphotons]
   Float_t         gen_photonE[1000];   //[ngenphotons]
   Int_t           gen_photonstatus[1000];   //[ngenphotons]
   Int_t           gen_photonMotherID[1000];   //[ngenphotons]
   Float_t         gen_photonMotherPt[1000];   //[ngenphotons]
   Float_t         gen_photonMotherEta[1000];   //[ngenphotons]
   Float_t         gen_photonMotherPhi[1000];   //[ngenphotons]
   Int_t           gen_photonMotherStatus[1000];   //[ngenphotons]
   Int_t           gen_photonGrandmotherID[1000];   //[ngenphotons]
   Float_t         gen_photonGrandmotherPt[1000];   //[ngenphotons]
   Float_t         gen_photonGrandmotherEta[1000];   //[ngenphotons]
   Float_t         gen_photonGrandmotherPhi[1000];   //[ngenphotons]
   Int_t           gen_photonGrandmotherStatus[1000];   //[ngenphotons]
   Int_t           nhardphotons;

   Float_t         gen_hardphotonpt[100];   //[nhardphotons]
   Float_t         gen_hardphotoneta[100];   //[nhardphotons]
   Float_t         gen_hardphotonphi[100];   //[nhardphotons]
   Float_t         gen_hardphotonpx[100];   //[nhardphotons]
   Float_t         gen_hardphotonpy[100];   //[nhardphotons]
   Float_t         gen_hardphotonpz[100];   //[nhardphotons]
   Float_t         gen_hardphotonE[100];   //[nhardphotons]
   Float_t         gen_gravitonpt;
   Float_t         gen_gravitoneta;
   Float_t         gen_gravitonphi;
   Float_t         gen_gravitonpx;
   Float_t         gen_gravitonpy;
   Float_t         gen_gravitonpz;
   Float_t         gen_gravitonE;
   Float_t         gen_Wdaughterpt[100];
   Float_t         gen_Wdaughtereta[100];
   Float_t         gen_Wdaughterphi[100];
   Float_t         gen_Wdaughterpx[100];
   Float_t         gen_Wdaughterpy[100];
   Float_t         gen_Wdaughterpz[100];
   Float_t         gen_WdaughterE[100];
   Int_t           gen_Wdaughter_charge[100];
   Int_t           gen_WdaughterID[100];
   Float_t         gen_Wbosonpt;
   Float_t         gen_Wbosoneta;
   Float_t         gen_Wbosonphi;
   Float_t         gen_Wbosonpx;
   Float_t         gen_Wbosonpy;
   Float_t         gen_Wbosonpz;
   Float_t         gen_WbosonE;
   Int_t           gen_Wbosoncharge;
   Int_t           gen_WbosonID;
   Float_t         gen_Zdaughterpt[100];
   Float_t         gen_Zdaughtereta[100];
   Float_t         gen_Zdaughterphi[100];
   Float_t         gen_Zdaughterpx[100];
   Float_t         gen_Zdaughterpy[100];
   Float_t         gen_Zdaughterpz[100];
   Float_t         gen_ZdaughterE[100];
   Int_t           gen_Zdaughter_charge[100];
   Int_t           gen_ZdaughterID[100];
   Float_t         gen_Zbosonpt;
   Float_t         gen_Zbosoneta;
   Float_t         gen_Zbosonphi;
   Float_t         gen_Zbosonpx;
   Float_t         gen_Zbosonpy;
   Float_t         gen_Zbosonpz;
   Float_t         gen_ZbosonE;
   Bool_t          is_signal_event;
   Bool_t          is_Z_event;
   Bool_t          is_W_event;
   Bool_t          is_Znunu_event;
   Bool_t          is_Zelec_event;
   Bool_t          is_Zmu_event;
   Bool_t          is_Ztau_event;
   Bool_t          is_Welec_event;
   Bool_t          is_Wmu_event;
   Bool_t          is_Wtau_event;
   Bool_t          is_SingleHardPhoton_event;
   Bool_t          is_diphoton_event;
   Bool_t          is_isr_photon_event;
   Int_t           n_signal_events;
   Int_t           n_Z_events;
   Int_t           n_W_events;
   Int_t           n_Znunu_events;
   Int_t           n_Zelec_events;
   Int_t           n_Zmu_events;
   Int_t           n_Ztau_events;
   Int_t           n_Welec_events;
   Int_t           n_Wmu_events;
   Int_t           n_Wtau_events;
   Int_t           n_SingleHardPhoton_events;
   Int_t           n_diphoton_events;
   Float_t         gen_MuonID[100];
   Float_t         gen_MuonStatus[100];
   Float_t         gen_MuonPt[100];
   Float_t         gen_MuonDaughterpt[100];
   Float_t         gen_MuonDaughtereta[100];
   Float_t         gen_MuonDaughterphi[100];
   Float_t         gen_MuonDaughterpx[100];
   Float_t         gen_MuonDaughterpy[100];
   Float_t         gen_MuonDaughterpz[100];
   Float_t         gen_MuonDaughterE[100];
   Int_t           gen_MuonDaughterCharge[100];
   Int_t           gen_MuonDaughterStatus[100];
   Int_t           gen_MuonDaughterID[100];
   Float_t         gen_tauID[100];
   Float_t         gen_tauStatus[100];
   Float_t         gen_tauPt[100];
   Float_t         gen_tauDaughterpt[100];
   Float_t         gen_tauDaughtereta[100];
   Float_t         gen_tauDaughterphi[100];
   Float_t         gen_tauDaughterpx[100];
   Float_t         gen_tauDaughterpy[100];
   Float_t         gen_tauDaughterpz[100];
   Float_t         gen_tauDaughterE[100];
   Int_t           gen_tauDaughterCharge[100];
   Int_t           gen_tauDaughterStatus[100];
   Int_t           gen_tauDaughterID[100];

    Bool_t          Scraping_isScrapingEvent;
   Int_t           Photon_n;
   Float_t         Photon_E[100];   //[Photon_n]
   Float_t         Photon_pt[100];   //[Photon_n]
   Float_t         Photon_eta[100];   //[Photon_n]
   Float_t         Photon_phi[100];   //[Photon_n]
   Float_t         Photon_theta[100];   //[Photon_n]
   Float_t         Photon_et[100];   //[Photon_n]
   Float_t         Photon_swissCross[100];   //[Photon_n]
   Float_t         Photonr9[100];   //[Photon_n]
   Float_t         Photon_e1x5[100];   //[Photon_n]
   Float_t         Photon_e2x5[100];   //[Photon_n]
   Float_t         Photon_e3x3[100];   //[Photon_n]
   Float_t         Photon_e5x5[100];   //[Photon_n]
   Float_t         Photon_r1x5[100];   //[Photon_n]
   Float_t         Photon_r2x5[100];   //[Photon_n]
   Float_t         Photon_maxEnergyXtal[100];   //[Photon_n]
   Float_t         Photon_SigmaEtaEta[100];   //[Photon_n]
   Float_t         Photon_SigmaIetaIeta[100];   //[Photon_n]
   Float_t         Photon_SigmaEtaPhi[100];   //[Photon_n]
   Float_t         Photon_SigmaIetaIphi[100];   //[Photon_n]
   Float_t         Photon_SigmaPhiPhi[100];   //[Photon_n]
   Float_t         Photon_SigmaIphiIphi[100];   //[Photon_n]
   Float_t         Photon_Roundness[100];   //[Photon_n]
   Float_t         Photon_Angle[100];   //[Photon_n]
   Float_t         Photon_ecalRecHitSumEtConeDR03[100];   //[Photon_n]
   Float_t         Photon_hcalTowerSumEtConeDR03[100];   //[Photon_n]
   Float_t         Photon_trkSumPtSolidConeDR03[100];   //[Photon_n]
   Float_t         Photon_trkSumPtHollowConeDR03[100];   //[Photon_n]
   Int_t           Photon_nTrkSolidConeDR03[100];   //[Photon_n]
   Int_t           Photon_nTrkHollowConeDR03[100];   //[Photon_n]
   Float_t         Photon_hcalDepth1TowerSumEtConeDR03[100];   //[Photon_n]
   Float_t         Photon_hcalDepth2TowerSumEtConeDR03[100];   //[Photon_n]
   Float_t         Photon_ecalRecHitSumEtConeDR04[100];   //[Photon_n]
   Float_t         Photon_hcalTowerSumEtConeDR04[100];   //[Photon_n]
   Float_t         Photon_trkSumPtSolidConeDR04[100];   //[Photon_n]
   Float_t         Photon_trkSumPtHollowConeDR04[100];   //[Photon_n]
   Int_t           Photon_nTrkSolidConeDR04[100];   //[Photon_n]
   Int_t           Photon_nTrkHollowConeDR04[100];   //[Photon_n]
   Float_t         Photon_hcalDepth1TowerSumEtConeDR04[100];   //[Photon_n]
   Float_t         Photon_hcalDepth2TowerSumEtConeDR04[100];   //[Photon_n]
   Bool_t          Photon_hasPixelSeed[100];   //[Photon_n]
   Bool_t          Photon_isEB[100];   //[Photon_n]
   Bool_t          Photon_isEE[100];   //[Photon_n]
   Bool_t          Photon_isEBGap[100];   //[Photon_n]
   Bool_t          Photon_isEEGap[100];   //[Photon_n]
   Bool_t          Photon_isEBEEGap[100];   //[Photon_n]
   Float_t         Photon_HoE[100];   //[Photon_n]
   Float_t         Photon_px[100];   //[Photon_n]
   Float_t         Photon_py[100];   //[Photon_n]
   Float_t         Photon_pz[100];   //[Photon_n]
   Float_t         Photon_vx[100];   //[Photon_n]
   Float_t         Photon_vy[100];   //[Photon_n]
   Float_t         Photon_vz[100];   //[Photon_n]
   Int_t           Photon_no_of_basic_clusters[100];   //[Photon_n]
   Float_t         Photon_sc_energy[100];   //[Photon_n]
   Float_t         Photon_sc_eta[100];   //[Photon_n]
   Float_t         Photon_sc_phi[100];   //[Photon_n]
   Float_t         Photon_sc_x[100];   //[Photon_n]
   Float_t         Photon_sc_y[100];   //[Photon_n]
   Float_t         Photon_sc_z[100];   //[Photon_n]
   Float_t         Photon_etaWidth[100];   //[Photon_n]
   Float_t         Photon_phiWidth[100];   //[Photon_n]
   Float_t         Photon_sc_et[100];   //[Photon_n]
   Float_t         matchphotonE[100];   //[Photon_n]
   Float_t         matchphotonpt[100];   //[Photon_n]
   Float_t         matchphotoneta[100];   //[Photon_n]
   Float_t         matchphotonphi[100];   //[Photon_n]
   Float_t         matchphotonpx[100];   //[Photon_n]
   Float_t         matchphotonpy[100];   //[Photon_n]
   Float_t         matchphotonpz[100];   //[Photon_n]
   Bool_t          ismatchedphoton[100];   //[Photon_n]
   Bool_t          Photon_hasConvTrk[100];   //[Photon_n]
   Int_t           Photon_ntracks[100];   //[Photon_n]
   Bool_t          Photon_isconverted[100];   //[Photon_n]
   Float_t         Photon_pairInvmass[100];   //[Photon_n]
   Float_t         Photon_pairCotThetaSeperation[100];   //[Photon_n]
   Float_t         Photon_pairmomentumX[100];   //[Photon_n]
   Float_t         Photon_pairmomentumY[100];   //[Photon_n]
   Float_t         Photon_pairmomentumZ[100];   //[Photon_n]
   Float_t         Photon_EoverP[100];   //[Photon_n]
   Float_t         Photon_ConvVx[100];   //[Photon_n]
   Float_t         Photon_ConvVy[100];   //[Photon_n]
   Float_t         Photon_ConvVz[100];   //[Photon_n]
   Float_t         Photon_ZOfPrimaryVertex[100];   //[Photon_n]
   Float_t         Photon_distOfMinimumApproach[100];   //[Photon_n]
   Float_t         Photon_dPhiTracksAtVtx[100];   //[Photon_n]
   Float_t         Photon_dPhiTracksAtEcal[100];   //[Photon_n]
   Float_t         Photon_dEtaTracksAtEcal[100];   //[Photon_n]
   Int_t           Photon_ncrys[100];   //[Photon_n]
   Float_t         Photon_timing_xtal[100][100];   //[Photon_n]
   Float_t         Photon_timingavg_xtal[100];   //[Photon_n]
   Float_t         Photon_energy_xtal[100][100];   //[Photon_n]
   Int_t           Photon_ieta_xtalEB[100][100];   //[Photon_n]
   Int_t           Photon_iphi_xtalEB[100][100];   //[Photon_n]
   Float_t         Photon_s9[100];   //[Photon_n]
   Int_t           HERecHit_subset_n;
   UInt_t          HERecHit_subset_detid[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_energy[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_time[10000];   //[HERecHit_subset_n]
   Int_t           HERecHit_subset_depth[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_phi[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_eta[100000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_x[100000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_y[100000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_z[100000];   //[HERecHit_subset_n]
   Int_t           CSCseg_n;
   Float_t         CSCseg_time[10000];   //[CSCseg_n]
   Float_t         CSCseg_x[10000];   //[CSCseg_n]
   Float_t         CSCseg_y[10000];   //[CSCseg_n]
   Float_t         CSCseg_z[10000];   //[CSCseg_n]
   Float_t         CSCseg_phi[10000];   //[CSCseg_n]
   Float_t         CSCseg_DirectionX[10000];   //[CSCseg_n]
   Float_t         CSCseg_DirectionY[10000];   //[CSCseg_n]
   Float_t         CSCseg_DirectionZ[10000];   //[CSCseg_n]

  Bool_t isBeamHaloIDTightPass;
  Bool_t isBeamHaloIDLoosePass;
    
  Bool_t isBeamHaloEcalLoosePass;
  Bool_t isBeamHaloHcalLoosePass;
  Bool_t isBeamHaloCSCLoosePass;
  Bool_t isBeamHaloGlobalLoosePass;
    
  Bool_t isBeamHaloEcalTightPass;
  Bool_t isBeamHaloHcalTightPass;
  Bool_t isBeamHaloCSCTightPass;
  Bool_t isBeamHaloGlobalTightPass;
    
  Bool_t isSmellsLikeHalo_Tag;
  Bool_t isLooseHalo_Tag;
  Bool_t isTightHalo_Tag;
  Bool_t isExtremeTightHalo_Tag;

   Int_t           RPChit_n;
   Float_t         RPChit_x[10000];   //[RPChit_n]
   Float_t         RPChit_y[10000];   //[RPChit_n]
   Float_t         RPChit_z[10000];   //[RPChit_n]
   Int_t           RPChit_BunchX[10000];   //[RPChit_n]
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

   Float_t         genMetPt;
   Float_t         genMetPx;
   Float_t         genMetPy;
   Float_t         genMetPhi;
   Float_t         genMetSumEt;
   Float_t         Delta_phi;
   Float_t         Delta_phiGEN;
   Float_t         PFMetPt[6];
   Float_t         PFMetPx[6];
   Float_t         PFMetPy[6];
   Float_t         PFMetPhi[6];
   Float_t         PFMetSumEt[6];
   Float_t         Delta_phiPF;
   Float_t         TCMetPt[6];
   Float_t         TCMetPx[6];
   Float_t         TCMetPy[6];
   Float_t         TCMetPhi[6];
   Float_t         TCMetSumEt[6];
   Float_t         Delta_phiTC;

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
   TBranch        *b_Vertex_n;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_tracksize;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_d0;   //!
   TBranch        *b_Vertex_isFake;   //!
   TBranch        *b_Track_n;   //!
   TBranch        *b_Track_px;   //!
   TBranch        *b_Track_py;   //!
   TBranch        *b_Track_pz;   //!
   TBranch        *b_Track_vx;   //!
   TBranch        *b_Track_vy;   //!
   TBranch        *b_Track_vz;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_phi;   //!
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
   TBranch        *b_Muon_n;   //!
   TBranch        *b_Muon_px;   //!
   TBranch        *b_Muon_py;   //!
   TBranch        *b_Muon_pz;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_InnerTrack_isNonnull;   //!
   TBranch        *b_Muon_OuterTrack_isNonnull;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_x;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_y;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_z;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_px;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_py;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_pz;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_x;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_y;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_z;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_px;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_py;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_pz;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_x;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_y;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_z;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_px;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_py;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_pz;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_x;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_y;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_z;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_px;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_py;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_pz;   //!
   TBranch        *b_CosmicMuon_n;   //!
   TBranch        *b_CosmicMuon_px;   //!
   TBranch        *b_CosmicMuon_py;   //!
   TBranch        *b_CosmicMuon_pz;   //!
   TBranch        *b_CosmicMuon_pt;   //!
   TBranch        *b_CosmicMuon_eta;   //!
   TBranch        *b_CosmicMuon_phi;   //!
   TBranch        *b_CosmicMuon_energy;   //!
   TBranch        *b_CosmicMuon_charge;   //!
   TBranch        *b_CosmicMuon_isGlobalMuon;   //!
   TBranch        *b_CosmicMuon_isTrackerMuon;   //!
   TBranch        *b_CosmicMuon_isStandAloneMuon;   //!
   TBranch        *b_CosmicMuon_InnerTrack_isNonnull;   //!
   TBranch        *b_CosmicMuon_OuterTrack_isNonnull;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_x;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_y;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_z;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_px;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_py;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_pz;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_x;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_y;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_z;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_px;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_py;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_pz;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_x;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_y;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_z;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_px;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_py;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_pz;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_x;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_y;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_z;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_px;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_py;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_pz;   //!
   TBranch        *b_Tau_n;   //!
   TBranch        *b_Tau_px;   //!
   TBranch        *b_Tau_py;   //!
   TBranch        *b_Tau_pz;   //!
   TBranch        *b_Tau_vx;   //!
   TBranch        *b_Tau_vy;   //!
   TBranch        *b_Tau_vz;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!

   TBranch        *b_gen_pthat;   //!
   TBranch        *b_ngenphotons;   //!
   TBranch        *b_gen_photonpt;   //!
   TBranch        *b_gen_photoneta;   //!
   TBranch        *b_gen_photonphi;   //!
   TBranch        *b_gen_photonpx;   //!
   TBranch        *b_gen_photonpy;   //!
   TBranch        *b_gen_photonpz;   //!
   TBranch        *b_gen_photonE;   //!
   TBranch        *b_gen_photonstatus;   //!
   TBranch        *b_gen_photonMotherID;   //!
   TBranch        *b_gen_photonMotherPt;   //!
   TBranch        *b_gen_photonMotherEta;   //!
   TBranch        *b_gen_photonMotherPhi;   //!
   TBranch        *b_gen_photonMotherStatus;   //!
   TBranch        *b_gen_photonGrandmotherID;   //!
   TBranch        *b_gen_photonGrandmotherPt;   //!
   TBranch        *b_gen_photonGrandmotherEta;   //!
   TBranch        *b_gen_photonGrandmotherPhi;   //!
   TBranch        *b_gen_photonGrandmotherStatus;   //!
   TBranch        *b_nhardphotons;   //!
   TBranch        *b_gen_hardphotonpt;   //!
   TBranch        *b_gen_hardphotoneta;   //!
   TBranch        *b_gen_hardphotonphi;   //!
   TBranch        *b_gen_hardphotonpx;   //!
   TBranch        *b_gen_hardphotonpy;   //!
   TBranch        *b_gen_hardphotonpz;   //!
   TBranch        *b_gen_hardphotonE;   //!
   TBranch        *b_gen_graviton_pt;   //!
   TBranch        *b_gen_graviton_eta;   //!
   TBranch        *b_gen_graviton_phi;   //!
   TBranch        *b_gen_graviton_px;   //!
   TBranch        *b_gen_graviton_py;   //!
   TBranch        *b_gen_graviton_pz;   //!
   TBranch        *b_gen_graviton_E;   //!
   TBranch        *b_gen_Wdaughter_pt;   //!
   TBranch        *b_gen_Wdaughter_eta;   //!
   TBranch        *b_gen_Wdaughter_phi;   //!
   TBranch        *b_gen_Wdaughter_px;   //!
   TBranch        *b_gen_Wdaughter_py;   //!
   TBranch        *b_gen_Wdaughter_pz;   //!
   TBranch        *b_gen_Wdaughter_E;   //!
   TBranch        *b_gen_Wdaughter_charge;   //!
   TBranch        *b_gen_Wdaughter_ID;   //!
   TBranch        *b_gen_Wboson_pt;   //!
   TBranch        *b_gen_Wboson_eta;   //!
   TBranch        *b_gen_Wboson_phi;   //!
   TBranch        *b_gen_Wboson_px;   //!
   TBranch        *b_gen_Wboson_py;   //!
   TBranch        *b_gen_Wboson_pz;   //!
   TBranch        *b_gen_Wboson_E;   //!
   TBranch        *b_gen_Wboson_charge;   //!
   TBranch        *b_gen_Wboson_ID;   //!
   TBranch        *b_gen_Zdaughter_pt;   //!
   TBranch        *b_gen_Zdaughter_eta;   //!
   TBranch        *b_gen_Zdaughter_phi;   //!
   TBranch        *b_gen_Zdaughter_px;   //!
   TBranch        *b_gen_Zdaughter_py;   //!
   TBranch        *b_gen_Zdaughter_pz;   //!
   TBranch        *b_gen_Zdaughter_E;   //!
   TBranch        *b_gen_Zdaughter_charge;   //!
   TBranch        *b_gen_Zdaughter_ID;   //!
   TBranch        *b_gen_Zboson_pt;   //!
   TBranch        *b_gen_Zboson_eta;   //!
   TBranch        *b_gen_Zboson_phi;   //!
   TBranch        *b_gen_Zboson_px;   //!
   TBranch        *b_gen_Zboson_py;   //!
   TBranch        *b_gen_Zboson_pz;   //!
   TBranch        *b_gen_Zboson_E;   //!
   TBranch        *b_is_signal_event;   //!
   TBranch        *b_is_Z_event;   //!
   TBranch        *b_is_W_event;   //!
   TBranch        *b_is_Znunu_event;   //!
   TBranch        *b_is_Zelec_event;   //!
   TBranch        *b_is_Zmu_event;   //!
   TBranch        *b_is_Ztau_event;   //!
   TBranch        *b_is_Welec_event;   //!
   TBranch        *b_is_Wmu_event;   //!
   TBranch        *b_is_Wtau_event;   //!
   TBranch        *b_is_SingleHardPhoton_event;   //!
   TBranch        *b_is_diphoton_event;   //!
   TBranch        *b_is_isr_photon_event;   //!
   TBranch        *b_n_signal_events;   //!
   TBranch        *b_n_Z_events;   //!
   TBranch        *b_n_W_events;   //!
   TBranch        *b_n_Znunu_events;   //!
   TBranch        *b_n_Zelec_events;   //!
   TBranch        *b_n_Zmu_events;   //!
   TBranch        *b_n_Ztau_events;   //!
   TBranch        *b_n_Welec_events;   //!
   TBranch        *b_n_Wmu_events;   //!
   TBranch        *b_n_Wtau_events;   //!
   TBranch        *b_n_SingleHardPhoton_events;   //!
   TBranch        *b_n_diphoton_events;   //!
   TBranch        *b_gen_Muon_ID;   //!
   TBranch        *b_gen_Muon_Status;   //!
   TBranch        *b_gen_Muon_Pt;   //!
   TBranch        *b_gen_MuonDaughter_pt;   //!
   TBranch        *b_gen_MuonDaughter_eta;   //!
   TBranch        *b_gen_MuonDaughter_phi;   //!
   TBranch        *b_gen_MuonDaughter_px;   //!
   TBranch        *b_gen_MuonDaughter_py;   //!
   TBranch        *b_gen_MuonDaughter_pz;   //!
   TBranch        *b_gen_MuonDaughter_E;   //!
   TBranch        *b_gen_MuonDaughter_charge;   //!
   TBranch        *b_gen_MuonDaughter_status;   //!
   TBranch        *b_gen_MuonDaughter_ID;   //!
   TBranch        *b_gen_tau_ID;   //!
   TBranch        *b_gen_tau_Status;   //!
   TBranch        *b_gen_tau_Pt;   //!
   TBranch        *b_gen_tauDaughter_pt;   //!
   TBranch        *b_gen_tauDaughter_eta;   //!
   TBranch        *b_gen_tauDaughter_phi;   //!
   TBranch        *b_gen_tauDaughter_px;   //!
   TBranch        *b_gen_tauDaughter_py;   //!
   TBranch        *b_gen_tauDaughter_pz;   //!
   TBranch        *b_gen_tauDaughter_E;   //!
   TBranch        *b_gen_tauDaughter_charge;   //!
   TBranch        *b_gen_tauDaughter_status;   //!
   TBranch        *b_gen_tauDaughter_ID;   //!

   TBranch        *b_Scraping_isScrapingEvent;//!
   TBranch        *b_Photon_n;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_theta;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_swissCross;   //!
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
   TBranch        *b_Photon_s9;   //!
   TBranch        *b_HERecHit_subset_n;   //!
   TBranch        *b_HERecHit_subset_detid;   //!
   TBranch        *b_HERecHit_subset_energy;   //!
   TBranch        *b_HERecHit_subset_time;   //!
   TBranch        *b_HERecHit_subset_depth;   //!
   TBranch        *b_HERecHit_subset_phi;   //!
   TBranch        *b_HERecHit_subset_eta;   //!
   TBranch        *b_HERecHit_subset_x;   //!
   TBranch        *b_HERecHit_subset_y;   //!
   TBranch        *b_HERecHit_subset_z;   //!
   TBranch        *b_CSCseg_n;   //!
   TBranch        *b_CSCseg_time;   //!
   TBranch        *b_CSCseg_x;   //!
   TBranch        *b_CSCseg_y;   //!
   TBranch        *b_CSCseg_z;   //!
   TBranch        *b_CSCseg_phi;   //!
   TBranch        *b_CSCseg_DirectionX;   //!
   TBranch        *b_CSCseg_DirectionY;   //!
   TBranch        *b_CSCseg_DirectionZ;   //!

   TBranch       *b_isBeamHaloIDTightPass;
   TBranch       *b_isBeamHaloIDLoosePass;
   TBranch       *b_isBeamHaloEcalLoosePass;
   TBranch       *b_isBeamHaloHcalLoosePass;
   TBranch       *b_isBeamHaloCSCLoosePass;
   TBranch       *b_isBeamHaloGlobalLoosePass;
   TBranch       *b_isBeamHaloEcalTightPass;
   TBranch       *b_isBeamHaloHcalTightPass;
   TBranch       *b_isBeamHaloCSCTightPass;
   TBranch       *b_isBeamHaloGlobalTightPass;
   TBranch       *b_isSmellsLikeHalo_Tag;
   TBranch       *b_isLooseHalo_Tag;
   TBranch       *b_isTightHalo_Tag;
   TBranch       *b_isExtremeTightHalo_Tag;


   TBranch        *b_RPChit_n;   //!
   TBranch        *b_RPChit_x;   //!
   TBranch        *b_RPChit_y;   //!
   TBranch        *b_RPChit_z;   //!
   TBranch        *b_RPChit_BunchX;   //!
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

   TBranch        *b_genMetPt;   //!
   TBranch        *b_genMetPx;   //!
   TBranch        *b_genMetPy;   //!
   TBranch        *b_genMetPhi;   //!
   TBranch        *b_genMetSumEt;   //!
   TBranch        *b_Delta_phi;   //!
   TBranch        *b_Delta_phiGEN;   //!
   TBranch        *b_PFMetPt;   //!
   TBranch        *b_PFMetPx;   //!
   TBranch        *b_PFMetPy;   //!
   TBranch        *b_PFMetPhi;   //!
   TBranch        *b_PFMetSumEt;   //!
   TBranch        *b_Delta_phiPF;   //!
   TBranch        *b_TCMetPt;   //!
   TBranch        *b_TCMetPx;   //!
   TBranch        *b_TCMetPy;   //!
   TBranch        *b_TCMetPhi;   //!
   TBranch        *b_TCMetSumEt;   //!
   TBranch        *b_Delta_phiTC;   //!

   myEvent(TTree *tree=0);
   virtual ~myEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myEvent_cxx
myEvent::myEvent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcache:///pnfs/cms/WAX/11/store/user/sushil/MonoPhoton/387_Ntuples_V26/Data_noSkim/Data_A/Histo_Data_A_10_0_L0n.root");
      if (!f) {
        // f = new TFile("dcache:///pnfs/cms/WAX/11/store/user/sushil/MonoPhoton/387_Ntuples_V26/WENu_Gamma/Histo_Wenu_Gamma_10_1_fNU.root");
      }
      tree = (TTree*)gDirectory->Get("myEvent");

   }
   Init(tree);
}

myEvent::~myEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myEvent::Init(TTree *tree)
{
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
   fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
   fChain->SetBranchAddress("Vertex_x", Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex_tracksize", Vertex_tracksize, &b_Vertex_tracksize);
   fChain->SetBranchAddress("Vertex_ndof", Vertex_ndof, &b_Vertex_ndof);
   fChain->SetBranchAddress("Vertex_chi2", Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_d0", Vertex_d0, &b_Vertex_d0);
   fChain->SetBranchAddress("Vertex_isFake", Vertex_isFake, &b_Vertex_isFake);
   fChain->SetBranchAddress("Track_n", &Track_n, &b_Track_n);
   fChain->SetBranchAddress("Track_px", Track_px, &b_Track_px);
   fChain->SetBranchAddress("Track_py", Track_py, &b_Track_py);
   fChain->SetBranchAddress("Track_pz", Track_pz, &b_Track_pz);
   fChain->SetBranchAddress("Track_vx", Track_vx, &b_Track_vx);
   fChain->SetBranchAddress("Track_vy", Track_vy, &b_Track_vy);
   fChain->SetBranchAddress("Track_vz", Track_vz, &b_Track_vz);
   fChain->SetBranchAddress("Track_pt", Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_eta", Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_phi", Track_phi, &b_Track_phi);
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
   fChain->SetBranchAddress("Muon_n", &Muon_n, &b_Muon_n);
   fChain->SetBranchAddress("Muon_px", Muon_px, &b_Muon_px);
   fChain->SetBranchAddress("Muon_py", Muon_py, &b_Muon_py);
   fChain->SetBranchAddress("Muon_pz", Muon_pz, &b_Muon_pz);
   fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
   fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
   fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_energy", Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_InnerTrack_isNonnull", Muon_InnerTrack_isNonnull, &b_Muon_InnerTrack_isNonnull);
   fChain->SetBranchAddress("Muon_OuterTrack_isNonnull", Muon_OuterTrack_isNonnull, &b_Muon_OuterTrack_isNonnull);
   fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_x", Muon_OuterTrack_InnerPoint_x, &b_Muon_OuterTrack_InnerPoint_x);
   fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_y", Muon_OuterTrack_InnerPoint_y, &b_Muon_OuterTrack_InnerPoint_y);
   fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_z", Muon_OuterTrack_InnerPoint_z, &b_Muon_OuterTrack_InnerPoint_z);
   fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_px", Muon_OuterTrack_InnerPoint_px, &b_Muon_OuterTrack_InnerPoint_px);
   fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_py", Muon_OuterTrack_InnerPoint_py, &b_Muon_OuterTrack_InnerPoint_py);
   fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_pz", Muon_OuterTrack_InnerPoint_pz, &b_Muon_OuterTrack_InnerPoint_pz);
   fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_x", Muon_OuterTrack_OuterPoint_x, &b_Muon_OuterTrack_OuterPoint_x);
   fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_y", Muon_OuterTrack_OuterPoint_y, &b_Muon_OuterTrack_OuterPoint_y);
   fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_z", Muon_OuterTrack_OuterPoint_z, &b_Muon_OuterTrack_OuterPoint_z);
   fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_px", Muon_OuterTrack_OuterPoint_px, &b_Muon_OuterTrack_OuterPoint_px);
   fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_py", Muon_OuterTrack_OuterPoint_py, &b_Muon_OuterTrack_OuterPoint_py);
   fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_pz", Muon_OuterTrack_OuterPoint_pz, &b_Muon_OuterTrack_OuterPoint_pz);
   fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_x", Muon_InnerTrack_InnerPoint_x, &b_Muon_InnerTrack_InnerPoint_x);
   fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_y", Muon_InnerTrack_InnerPoint_y, &b_Muon_InnerTrack_InnerPoint_y);
   fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_z", Muon_InnerTrack_InnerPoint_z, &b_Muon_InnerTrack_InnerPoint_z);
   fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_px", Muon_InnerTrack_InnerPoint_px, &b_Muon_InnerTrack_InnerPoint_px);
   fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_py", Muon_InnerTrack_InnerPoint_py, &b_Muon_InnerTrack_InnerPoint_py);
   fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_pz", Muon_InnerTrack_InnerPoint_pz, &b_Muon_InnerTrack_InnerPoint_pz);
   fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_x", Muon_InnerTrack_OuterPoint_x, &b_Muon_InnerTrack_OuterPoint_x);
   fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_y", Muon_InnerTrack_OuterPoint_y, &b_Muon_InnerTrack_OuterPoint_y);
   fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_z", Muon_InnerTrack_OuterPoint_z, &b_Muon_InnerTrack_OuterPoint_z);
   fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_px", Muon_InnerTrack_OuterPoint_px, &b_Muon_InnerTrack_OuterPoint_px);
   fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_py", Muon_InnerTrack_OuterPoint_py, &b_Muon_InnerTrack_OuterPoint_py);
   fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_pz", Muon_InnerTrack_OuterPoint_pz, &b_Muon_InnerTrack_OuterPoint_pz);
   fChain->SetBranchAddress("CosmicMuon_n", &CosmicMuon_n, &b_CosmicMuon_n);
   fChain->SetBranchAddress("CosmicMuon_px", CosmicMuon_px, &b_CosmicMuon_px);
   fChain->SetBranchAddress("CosmicMuon_py", CosmicMuon_py, &b_CosmicMuon_py);
   fChain->SetBranchAddress("CosmicMuon_pz", CosmicMuon_pz, &b_CosmicMuon_pz);
   fChain->SetBranchAddress("CosmicMuon_pt", CosmicMuon_pt, &b_CosmicMuon_pt);
   fChain->SetBranchAddress("CosmicMuon_eta", CosmicMuon_eta, &b_CosmicMuon_eta);
   fChain->SetBranchAddress("CosmicMuon_phi", CosmicMuon_phi, &b_CosmicMuon_phi);
   fChain->SetBranchAddress("CosmicMuon_energy", CosmicMuon_energy, &b_CosmicMuon_energy);
   fChain->SetBranchAddress("CosmicMuon_charge", CosmicMuon_charge, &b_CosmicMuon_charge);
   fChain->SetBranchAddress("CosmicMuon_isGlobalMuon", CosmicMuon_isGlobalMuon, &b_CosmicMuon_isGlobalMuon);
   fChain->SetBranchAddress("CosmicMuon_isTrackerMuon", CosmicMuon_isTrackerMuon, &b_CosmicMuon_isTrackerMuon);
   fChain->SetBranchAddress("CosmicMuon_isStandAloneMuon", CosmicMuon_isStandAloneMuon, &b_CosmicMuon_isStandAloneMuon);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_isNonnull", CosmicMuon_InnerTrack_isNonnull, &b_CosmicMuon_InnerTrack_isNonnull);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_isNonnull", CosmicMuon_OuterTrack_isNonnull, &b_CosmicMuon_OuterTrack_isNonnull);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_x", CosmicMuon_OuterTrack_InnerPoint_x, &b_CosmicMuon_OuterTrack_InnerPoint_x);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_y", CosmicMuon_OuterTrack_InnerPoint_y, &b_CosmicMuon_OuterTrack_InnerPoint_y);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_z", CosmicMuon_OuterTrack_InnerPoint_z, &b_CosmicMuon_OuterTrack_InnerPoint_z);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_px", CosmicMuon_OuterTrack_InnerPoint_px, &b_CosmicMuon_OuterTrack_InnerPoint_px);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_py", CosmicMuon_OuterTrack_InnerPoint_py, &b_CosmicMuon_OuterTrack_InnerPoint_py);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_pz", CosmicMuon_OuterTrack_InnerPoint_pz, &b_CosmicMuon_OuterTrack_InnerPoint_pz);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_x", CosmicMuon_OuterTrack_OuterPoint_x, &b_CosmicMuon_OuterTrack_OuterPoint_x);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_y", CosmicMuon_OuterTrack_OuterPoint_y, &b_CosmicMuon_OuterTrack_OuterPoint_y);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_z", CosmicMuon_OuterTrack_OuterPoint_z, &b_CosmicMuon_OuterTrack_OuterPoint_z);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_px", CosmicMuon_OuterTrack_OuterPoint_px, &b_CosmicMuon_OuterTrack_OuterPoint_px);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_py", CosmicMuon_OuterTrack_OuterPoint_py, &b_CosmicMuon_OuterTrack_OuterPoint_py);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_pz", CosmicMuon_OuterTrack_OuterPoint_pz, &b_CosmicMuon_OuterTrack_OuterPoint_pz);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_x", CosmicMuon_InnerTrack_InnerPoint_x, &b_CosmicMuon_InnerTrack_InnerPoint_x);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_y", CosmicMuon_InnerTrack_InnerPoint_y, &b_CosmicMuon_InnerTrack_InnerPoint_y);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_z", CosmicMuon_InnerTrack_InnerPoint_z, &b_CosmicMuon_InnerTrack_InnerPoint_z);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_px", CosmicMuon_InnerTrack_InnerPoint_px, &b_CosmicMuon_InnerTrack_InnerPoint_px);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_py", CosmicMuon_InnerTrack_InnerPoint_py, &b_CosmicMuon_InnerTrack_InnerPoint_py);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_pz", CosmicMuon_InnerTrack_InnerPoint_pz, &b_CosmicMuon_InnerTrack_InnerPoint_pz);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_x", CosmicMuon_InnerTrack_OuterPoint_x, &b_CosmicMuon_InnerTrack_OuterPoint_x);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_y", CosmicMuon_InnerTrack_OuterPoint_y, &b_CosmicMuon_InnerTrack_OuterPoint_y);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_z", CosmicMuon_InnerTrack_OuterPoint_z, &b_CosmicMuon_InnerTrack_OuterPoint_z);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_px", CosmicMuon_InnerTrack_OuterPoint_px, &b_CosmicMuon_InnerTrack_OuterPoint_px);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_py", CosmicMuon_InnerTrack_OuterPoint_py, &b_CosmicMuon_InnerTrack_OuterPoint_py);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_pz", CosmicMuon_InnerTrack_OuterPoint_pz, &b_CosmicMuon_InnerTrack_OuterPoint_pz);
   fChain->SetBranchAddress("Tau_n", &Tau_n, &b_Tau_n);
   fChain->SetBranchAddress("Tau_px", Tau_px, &b_Tau_px);
   fChain->SetBranchAddress("Tau_py", Tau_py, &b_Tau_py);
   fChain->SetBranchAddress("Tau_pz", Tau_pz, &b_Tau_pz);
   fChain->SetBranchAddress("Tau_vx", Tau_vx, &b_Tau_vx);
   fChain->SetBranchAddress("Tau_vy", Tau_vy, &b_Tau_vy);
   fChain->SetBranchAddress("Tau_vz", Tau_vz, &b_Tau_vz);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_energy", Tau_energy, &b_Tau_energy);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("gen_pthat", &gen_pthat, &b_gen_pthat);
   fChain->SetBranchAddress("ngenphotons", &ngenphotons, &b_ngenphotons);
   fChain->SetBranchAddress("gen_photonpt", gen_photonpt, &b_gen_photonpt);
   fChain->SetBranchAddress("gen_photoneta", gen_photoneta, &b_gen_photoneta);
   fChain->SetBranchAddress("gen_photonphi", gen_photonphi, &b_gen_photonphi);
   fChain->SetBranchAddress("gen_photonpx", gen_photonpx, &b_gen_photonpx);
   fChain->SetBranchAddress("gen_photonpy", gen_photonpy, &b_gen_photonpy);
   fChain->SetBranchAddress("gen_photonpz", gen_photonpz, &b_gen_photonpz);
   fChain->SetBranchAddress("gen_photonE", gen_photonE, &b_gen_photonE);
   fChain->SetBranchAddress("gen_photonstatus", gen_photonstatus, &b_gen_photonstatus);
   fChain->SetBranchAddress("gen_photonMotherID", gen_photonMotherID, &b_gen_photonMotherID);
   fChain->SetBranchAddress("gen_photonMotherPt", gen_photonMotherPt, &b_gen_photonMotherPt);
   fChain->SetBranchAddress("gen_photonMotherEta", gen_photonMotherEta, &b_gen_photonMotherEta);
   fChain->SetBranchAddress("gen_photonMotherPhi", gen_photonMotherPhi, &b_gen_photonMotherPhi);
   fChain->SetBranchAddress("gen_photonMotherStatus", gen_photonMotherStatus, &b_gen_photonMotherStatus);
   fChain->SetBranchAddress("gen_photonGrandmotherID", gen_photonGrandmotherID, &b_gen_photonGrandmotherID);
   fChain->SetBranchAddress("gen_photonGrandmotherPt", gen_photonGrandmotherPt, &b_gen_photonGrandmotherPt);
   fChain->SetBranchAddress("gen_photonGrandmotherEta", gen_photonGrandmotherEta, &b_gen_photonGrandmotherEta);
   fChain->SetBranchAddress("gen_photonGrandmotherPhi", gen_photonGrandmotherPhi, &b_gen_photonGrandmotherPhi);
   fChain->SetBranchAddress("gen_photonGrandmotherStatus", gen_photonGrandmotherStatus, &b_gen_photonGrandmotherStatus);
   fChain->SetBranchAddress("nhardphotons", &nhardphotons, &b_nhardphotons);
   fChain->SetBranchAddress("gen_hardphotonpt", gen_hardphotonpt, &b_gen_hardphotonpt);
   fChain->SetBranchAddress("gen_hardphotoneta", gen_hardphotoneta, &b_gen_hardphotoneta);
   fChain->SetBranchAddress("gen_hardphotonphi", gen_hardphotonphi, &b_gen_hardphotonphi);
   fChain->SetBranchAddress("gen_hardphotonpx", gen_hardphotonpx, &b_gen_hardphotonpx);
   fChain->SetBranchAddress("gen_hardphotonpy", gen_hardphotonpy, &b_gen_hardphotonpy);
   fChain->SetBranchAddress("gen_hardphotonpz", gen_hardphotonpz, &b_gen_hardphotonpz);
   fChain->SetBranchAddress("gen_hardphotonE", gen_hardphotonE, &b_gen_hardphotonE);
   fChain->SetBranchAddress("gen_gravitonpt", &gen_gravitonpt, &b_gen_graviton_pt);
   fChain->SetBranchAddress("gen_gravitoneta", &gen_gravitoneta, &b_gen_graviton_eta);
   fChain->SetBranchAddress("gen_gravitonphi", &gen_gravitonphi, &b_gen_graviton_phi);
   fChain->SetBranchAddress("gen_gravitonpx", &gen_gravitonpx, &b_gen_graviton_px);
   fChain->SetBranchAddress("gen_gravitonpy", &gen_gravitonpy, &b_gen_graviton_py);
   fChain->SetBranchAddress("gen_gravitonpz", &gen_gravitonpz, &b_gen_graviton_pz);
   fChain->SetBranchAddress("gen_gravitonE", &gen_gravitonE, &b_gen_graviton_E);
   fChain->SetBranchAddress("gen_Wdaughterpt", gen_Wdaughterpt, &b_gen_Wdaughter_pt);
   fChain->SetBranchAddress("gen_Wdaughtereta", gen_Wdaughtereta, &b_gen_Wdaughter_eta);
   fChain->SetBranchAddress("gen_Wdaughterphi", gen_Wdaughterphi, &b_gen_Wdaughter_phi);
   fChain->SetBranchAddress("gen_Wdaughterpx", gen_Wdaughterpx, &b_gen_Wdaughter_px);
   fChain->SetBranchAddress("gen_Wdaughterpy", gen_Wdaughterpy, &b_gen_Wdaughter_py);
   fChain->SetBranchAddress("gen_Wdaughterpz", gen_Wdaughterpz, &b_gen_Wdaughter_pz);
   fChain->SetBranchAddress("gen_WdaughterE", gen_WdaughterE, &b_gen_Wdaughter_E);
   fChain->SetBranchAddress("gen_Wdaughter_charge", gen_Wdaughter_charge, &b_gen_Wdaughter_charge);
   fChain->SetBranchAddress("gen_WdaughterID", gen_WdaughterID, &b_gen_Wdaughter_ID);
   fChain->SetBranchAddress("gen_Wbosonpt", &gen_Wbosonpt, &b_gen_Wboson_pt);
   fChain->SetBranchAddress("gen_Wbosoneta", &gen_Wbosoneta, &b_gen_Wboson_eta);
   fChain->SetBranchAddress("gen_Wbosonphi", &gen_Wbosonphi, &b_gen_Wboson_phi);
   fChain->SetBranchAddress("gen_Wbosonpx", &gen_Wbosonpx, &b_gen_Wboson_px);
   fChain->SetBranchAddress("gen_Wbosonpy", &gen_Wbosonpy, &b_gen_Wboson_py);
   fChain->SetBranchAddress("gen_Wbosonpz", &gen_Wbosonpz, &b_gen_Wboson_pz);
   fChain->SetBranchAddress("gen_WbosonE", &gen_WbosonE, &b_gen_Wboson_E);
   fChain->SetBranchAddress("gen_Wbosoncharge", &gen_Wbosoncharge, &b_gen_Wboson_charge);
   fChain->SetBranchAddress("gen_WbosonID", &gen_WbosonID, &b_gen_Wboson_ID);
   fChain->SetBranchAddress("gen_Zdaughterpt", gen_Zdaughterpt, &b_gen_Zdaughter_pt);
   fChain->SetBranchAddress("gen_Zdaughtereta", gen_Zdaughtereta, &b_gen_Zdaughter_eta);
   fChain->SetBranchAddress("gen_Zdaughterphi", gen_Zdaughterphi, &b_gen_Zdaughter_phi);
   fChain->SetBranchAddress("gen_Zdaughterpx", gen_Zdaughterpx, &b_gen_Zdaughter_px);
   fChain->SetBranchAddress("gen_Zdaughterpy", gen_Zdaughterpy, &b_gen_Zdaughter_py);
   fChain->SetBranchAddress("gen_Zdaughterpz", gen_Zdaughterpz, &b_gen_Zdaughter_pz);
   fChain->SetBranchAddress("gen_ZdaughterE", gen_ZdaughterE, &b_gen_Zdaughter_E);
   fChain->SetBranchAddress("gen_Zdaughter_charge", gen_Zdaughter_charge, &b_gen_Zdaughter_charge);
   fChain->SetBranchAddress("gen_ZdaughterID", gen_ZdaughterID, &b_gen_Zdaughter_ID);
   fChain->SetBranchAddress("gen_Zbosonpt", &gen_Zbosonpt, &b_gen_Zboson_pt);
   fChain->SetBranchAddress("gen_Zbosoneta", &gen_Zbosoneta, &b_gen_Zboson_eta);
   fChain->SetBranchAddress("gen_Zbosonphi", &gen_Zbosonphi, &b_gen_Zboson_phi);
   fChain->SetBranchAddress("gen_Zbosonpx", &gen_Zbosonpx, &b_gen_Zboson_px);
   fChain->SetBranchAddress("gen_Zbosonpy", &gen_Zbosonpy, &b_gen_Zboson_py);
   fChain->SetBranchAddress("gen_Zbosonpz", &gen_Zbosonpz, &b_gen_Zboson_pz);
   fChain->SetBranchAddress("gen_ZbosonE", &gen_ZbosonE, &b_gen_Zboson_E);
   fChain->SetBranchAddress("is_signal_event", &is_signal_event, &b_is_signal_event);
   fChain->SetBranchAddress("is_Z_event", &is_Z_event, &b_is_Z_event);
   fChain->SetBranchAddress("is_W_event", &is_W_event, &b_is_W_event);
   fChain->SetBranchAddress("is_Znunu_event", &is_Znunu_event, &b_is_Znunu_event);
   fChain->SetBranchAddress("is_Zelec_event", &is_Zelec_event, &b_is_Zelec_event);
   fChain->SetBranchAddress("is_Zmu_event", &is_Zmu_event, &b_is_Zmu_event);
   fChain->SetBranchAddress("is_Ztau_event", &is_Ztau_event, &b_is_Ztau_event);
   fChain->SetBranchAddress("is_Welec_event", &is_Welec_event, &b_is_Welec_event);
   fChain->SetBranchAddress("is_Wmu_event", &is_Wmu_event, &b_is_Wmu_event);
   fChain->SetBranchAddress("is_Wtau_event", &is_Wtau_event, &b_is_Wtau_event);
   fChain->SetBranchAddress("is_SingleHardPhoton_event", &is_SingleHardPhoton_event, &b_is_SingleHardPhoton_event);
   fChain->SetBranchAddress("is_diphoton_event", &is_diphoton_event, &b_is_diphoton_event);
   fChain->SetBranchAddress("is_isr_photon_event", &is_isr_photon_event, &b_is_isr_photon_event);
   fChain->SetBranchAddress("Scraping_isScrapingEvent", &Scraping_isScrapingEvent, &b_Scraping_isScrapingEvent);
   fChain->SetBranchAddress("n_signal_events", &n_signal_events, &b_n_signal_events);
   fChain->SetBranchAddress("n_Z_events", &n_Z_events, &b_n_Z_events);
   fChain->SetBranchAddress("n_W_events", &n_W_events, &b_n_W_events);
   fChain->SetBranchAddress("n_Znunu_events", &n_Znunu_events, &b_n_Znunu_events);
   fChain->SetBranchAddress("n_Zelec_events", &n_Zelec_events, &b_n_Zelec_events);
   fChain->SetBranchAddress("n_Zmu_events", &n_Zmu_events, &b_n_Zmu_events);
   fChain->SetBranchAddress("n_Ztau_events", &n_Ztau_events, &b_n_Ztau_events);
   fChain->SetBranchAddress("n_Welec_events", &n_Welec_events, &b_n_Welec_events);
   fChain->SetBranchAddress("n_Wmu_events", &n_Wmu_events, &b_n_Wmu_events);
   fChain->SetBranchAddress("n_Wtau_events", &n_Wtau_events, &b_n_Wtau_events);
   fChain->SetBranchAddress("n_SingleHardPhoton_events", &n_SingleHardPhoton_events, &b_n_SingleHardPhoton_events);
   fChain->SetBranchAddress("n_diphoton_events", &n_diphoton_events, &b_n_diphoton_events);
   fChain->SetBranchAddress("gen_MuonID", gen_MuonID, &b_gen_Muon_ID);
   fChain->SetBranchAddress("gen_MuonStatus", gen_MuonStatus, &b_gen_Muon_Status);
   fChain->SetBranchAddress("gen_MuonPt", gen_MuonPt, &b_gen_Muon_Pt);
   fChain->SetBranchAddress("gen_MuonDaughterpt", gen_MuonDaughterpt, &b_gen_MuonDaughter_pt);
   fChain->SetBranchAddress("gen_MuonDaughtereta", gen_MuonDaughtereta, &b_gen_MuonDaughter_eta);
   fChain->SetBranchAddress("gen_MuonDaughterphi", gen_MuonDaughterphi, &b_gen_MuonDaughter_phi);
   fChain->SetBranchAddress("gen_MuonDaughterpx", gen_MuonDaughterpx, &b_gen_MuonDaughter_px);
   fChain->SetBranchAddress("gen_MuonDaughterpy", gen_MuonDaughterpy, &b_gen_MuonDaughter_py);
   fChain->SetBranchAddress("gen_MuonDaughterpz", gen_MuonDaughterpz, &b_gen_MuonDaughter_pz);
   fChain->SetBranchAddress("gen_MuonDaughterE", gen_MuonDaughterE, &b_gen_MuonDaughter_E);
   fChain->SetBranchAddress("gen_MuonDaughterCharge", gen_MuonDaughterCharge, &b_gen_MuonDaughter_charge);
   fChain->SetBranchAddress("gen_MuonDaughterStatus", gen_MuonDaughterStatus, &b_gen_MuonDaughter_status);
   fChain->SetBranchAddress("gen_MuonDaughterID", gen_MuonDaughterID, &b_gen_MuonDaughter_ID);
   fChain->SetBranchAddress("gen_tauID", gen_tauID, &b_gen_tau_ID);
   fChain->SetBranchAddress("gen_tauStatus", gen_tauStatus, &b_gen_tau_Status);
   fChain->SetBranchAddress("gen_tauPt", gen_tauPt, &b_gen_tau_Pt);
   fChain->SetBranchAddress("gen_tauDaughterpt", gen_tauDaughterpt, &b_gen_tauDaughter_pt);
   fChain->SetBranchAddress("gen_tauDaughtereta", gen_tauDaughtereta, &b_gen_tauDaughter_eta);
   fChain->SetBranchAddress("gen_tauDaughterphi", gen_tauDaughterphi, &b_gen_tauDaughter_phi);
   fChain->SetBranchAddress("gen_tauDaughterpx", gen_tauDaughterpx, &b_gen_tauDaughter_px);
   fChain->SetBranchAddress("gen_tauDaughterpy", gen_tauDaughterpy, &b_gen_tauDaughter_py);
   fChain->SetBranchAddress("gen_tauDaughterpz", gen_tauDaughterpz, &b_gen_tauDaughter_pz);
   fChain->SetBranchAddress("gen_tauDaughterE", gen_tauDaughterE, &b_gen_tauDaughter_E);
   fChain->SetBranchAddress("gen_tauDaughterCharge", gen_tauDaughterCharge, &b_gen_tauDaughter_charge);
   fChain->SetBranchAddress("gen_tauDaughterStatus", gen_tauDaughterStatus, &b_gen_tauDaughter_status);
   fChain->SetBranchAddress("gen_tauDaughterID", gen_tauDaughterID, &b_gen_tauDaughter_ID);

   fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
   fChain->SetBranchAddress("Photon_E", Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_theta", Photon_theta, &b_Photon_theta);
   fChain->SetBranchAddress("Photon_et", Photon_et, &b_Photon_et);
   fChain->SetBranchAddress("Photon_swissCross", Photon_swissCross, &b_Photon_swissCross);
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
   fChain->SetBranchAddress("Photon_s9", Photon_s9, &b_Photon_s9);
   fChain->SetBranchAddress("HERecHit_subset_n", &HERecHit_subset_n, &b_HERecHit_subset_n);
   fChain->SetBranchAddress("HERecHit_subset_detid", HERecHit_subset_detid, &b_HERecHit_subset_detid);
   fChain->SetBranchAddress("HERecHit_subset_energy", HERecHit_subset_energy, &b_HERecHit_subset_energy);
   fChain->SetBranchAddress("HERecHit_subset_time", HERecHit_subset_time, &b_HERecHit_subset_time);
   fChain->SetBranchAddress("HERecHit_subset_depth", HERecHit_subset_depth, &b_HERecHit_subset_depth);
   fChain->SetBranchAddress("HERecHit_subset_phi", HERecHit_subset_phi, &b_HERecHit_subset_phi);
   fChain->SetBranchAddress("HERecHit_subset_eta", HERecHit_subset_eta, &b_HERecHit_subset_eta);
   fChain->SetBranchAddress("HERecHit_subset_x", HERecHit_subset_x, &b_HERecHit_subset_x);
   fChain->SetBranchAddress("HERecHit_subset_y", HERecHit_subset_y, &b_HERecHit_subset_y);
   fChain->SetBranchAddress("HERecHit_subset_z", HERecHit_subset_z, &b_HERecHit_subset_z);
   fChain->SetBranchAddress("CSCseg_n", &CSCseg_n, &b_CSCseg_n);
   fChain->SetBranchAddress("CSCseg_time", CSCseg_time, &b_CSCseg_time);
   fChain->SetBranchAddress("CSCseg_x", CSCseg_x, &b_CSCseg_x);
   fChain->SetBranchAddress("CSCseg_y", CSCseg_y, &b_CSCseg_y);
   fChain->SetBranchAddress("CSCseg_z", CSCseg_z, &b_CSCseg_z);
   fChain->SetBranchAddress("CSCseg_phi", CSCseg_phi, &b_CSCseg_phi);
   fChain->SetBranchAddress("CSCseg_DirectionX", CSCseg_DirectionX, &b_CSCseg_DirectionX);
   fChain->SetBranchAddress("CSCseg_DirectionY", CSCseg_DirectionY, &b_CSCseg_DirectionY);
   fChain->SetBranchAddress("CSCseg_DirectionZ", CSCseg_DirectionZ, &b_CSCseg_DirectionZ);

   fChain->SetBranchAddress("isBeamHaloIDTightPass", &isBeamHaloIDTightPass, &b_isBeamHaloIDTightPass);
   fChain->SetBranchAddress("isBeamHaloIDLoosePass", &isBeamHaloIDLoosePass, &b_isBeamHaloIDLoosePass);
   fChain->SetBranchAddress("isBeamHaloEcalLoosePass",&isBeamHaloEcalLoosePass, &b_isBeamHaloEcalLoosePass);
   fChain->SetBranchAddress("isBeamHaloHcalLoosePass", &isBeamHaloHcalLoosePass, &b_isBeamHaloHcalLoosePass);
   fChain->SetBranchAddress("isBeamHaloCSCLoosePass", &isBeamHaloCSCLoosePass, &b_isBeamHaloCSCLoosePass);
   fChain->SetBranchAddress("isBeamHaloGlobalLoosePass", &isBeamHaloGlobalLoosePass, &b_isBeamHaloGlobalLoosePass);
    
   fChain->SetBranchAddress("isBeamHaloEcalTightPass", &isBeamHaloEcalTightPass, &b_isBeamHaloEcalTightPass);
   fChain->SetBranchAddress("isBeamHaloHcalTightPass", &isBeamHaloHcalTightPass, &b_isBeamHaloHcalTightPass);
   fChain->SetBranchAddress("isBeamHaloCSCTightPass", &isBeamHaloCSCTightPass, &b_isBeamHaloCSCTightPass);
   fChain->SetBranchAddress("isBeamHaloGlobalTightPass", &isBeamHaloGlobalTightPass, &b_isBeamHaloGlobalTightPass);
    
   fChain->SetBranchAddress("isSmellsLikeHalo_Tag", &isSmellsLikeHalo_Tag, &b_isSmellsLikeHalo_Tag);
   fChain->SetBranchAddress("isLooseHalo_Tag", &isLooseHalo_Tag, &b_isLooseHalo_Tag);
   fChain->SetBranchAddress("isTightHalo_Tag", &isTightHalo_Tag,  &b_isTightHalo_Tag);
   fChain->SetBranchAddress("isExtremeTightHalo_Tag", &isExtremeTightHalo_Tag, &b_isExtremeTightHalo_Tag);
    

   fChain->SetBranchAddress("RPChit_n", &RPChit_n, &b_RPChit_n);
   fChain->SetBranchAddress("RPChit_x", RPChit_x, &b_RPChit_x);
   fChain->SetBranchAddress("RPChit_y", RPChit_y, &b_RPChit_y);
   fChain->SetBranchAddress("RPChit_z", RPChit_z, &b_RPChit_z);
   fChain->SetBranchAddress("RPChit_BunchX", RPChit_BunchX, &b_RPChit_BunchX);
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
   fChain->SetBranchAddress("genMetPt", &genMetPt, &b_genMetPt);
   fChain->SetBranchAddress("genMetPx", &genMetPx, &b_genMetPx);
   fChain->SetBranchAddress("genMetPy", &genMetPy, &b_genMetPy);
   fChain->SetBranchAddress("genMetPhi", &genMetPhi, &b_genMetPhi);
   fChain->SetBranchAddress("genMetSumEt", &genMetSumEt, &b_genMetSumEt);

   fChain->SetBranchAddress("Delta_phi", &Delta_phi, &b_Delta_phi);
   fChain->SetBranchAddress("Delta_phiGEN", &Delta_phiGEN, &b_Delta_phiGEN);
   fChain->SetBranchAddress("PFMetPt", PFMetPt, &b_PFMetPt);
   fChain->SetBranchAddress("PFMetPx", PFMetPx, &b_PFMetPx);
   fChain->SetBranchAddress("PFMetPy", PFMetPy, &b_PFMetPy);
   fChain->SetBranchAddress("PFMetPhi", PFMetPhi, &b_PFMetPhi);
   fChain->SetBranchAddress("PFMetSumEt", PFMetSumEt, &b_PFMetSumEt);
   fChain->SetBranchAddress("Delta_phiPF", &Delta_phiPF, &b_Delta_phiPF);
   fChain->SetBranchAddress("TCMetPt", TCMetPt, &b_TCMetPt);
   fChain->SetBranchAddress("TCMetPx", TCMetPx, &b_TCMetPx);
   fChain->SetBranchAddress("TCMetPy", TCMetPy, &b_TCMetPy);
   fChain->SetBranchAddress("TCMetPhi", TCMetPhi, &b_TCMetPhi);
   fChain->SetBranchAddress("TCMetSumEt", TCMetSumEt, &b_TCMetSumEt);
   fChain->SetBranchAddress("Delta_phiTC", &Delta_phiTC, &b_Delta_phiTC);
   Notify();
}

Bool_t myEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myEvent_cxx
