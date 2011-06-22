#define myCuts_cxx
#define pow2(x) pow(x,2.)
#include "myEvent.C"
#include "NonCollisionBG.h"
#include <stdexcept> 
#include <assert.h>
#include "string.h"

class myCuts : public myEvent,public NonCollisionBG { 
public: 
  myCuts();
  ~myCuts();
  void define_Tree();
  void Loop();

 
private: 
  //NonCollisionBG *bg;
  //NonCollisionBG *bg=new NonCollisionBG();
  TChain *fchain; 
  TFile *outfile;
  TTree *tree;	
  
  Int_t    NFilled;
  Int_t    NJFilled;
  Int_t    Nnevents;
  UInt_t   Nrun;
  UInt_t   Nevent;
  UInt_t   NluminosityBlock;
  UInt_t   NbeamCrossing;
  UInt_t   NtotalIntensityBeam1;
  UInt_t   NtotalIntensityBeam2;
  Float_t  NavgInsDelLumi;
  Float_t  NavgInsDelLumiErr;
  Float_t  NavgInsRecLumi;
  Float_t  NavgInsRecLumiErr;
  
  Float_t pthat;
  
  //VERTEX
  Int_t   vertex_n;
  Float_t vertex_x[100];   
  Float_t vertex_y[100];   
  Float_t vertex_z[100];   
  Float_t vertex_tracksize[100];   
  Float_t vertex_ndof[100];   
  Float_t vertex_chi2[100];   
  Float_t vertex_d0[100];   
  Bool_t  vertex_isFake[100];   
  
  //TRACK
  Int_t  track_n;
  Float_t track_px[100];   
  Float_t track_py[100];   
  Float_t track_pz[100];   
  Float_t track_pt[100];   
  Float_t track_eta[100];   
  Float_t track_phi[100];   
  //PHOTON
  Int_t   photon_n;
  //Float_t photon_iLICTD[100];
  Float_t photon_pt[100];
  Float_t photon_et[100];
  Float_t photon_E[100];
  Float_t photon_eta[100];
  Float_t photon_phi[100];
  Float_t photon_px[100];
  Float_t photon_py[100];
  Float_t photon_pz[100];
  Float_t photon_theta[100];

  Float_t photon_e2e9[100];
  Float_t photon_e6e2[100];
  Float_t photon_e4e1[100];
  Float_t photon_swissCross[100];
  Float_t photon_r9[100];
  Float_t photon_e1x5[100];
  Float_t photon_e2x5[100];
  Float_t photon_e3x3[100];
  Float_t photon_e5x5[100];
  Float_t photon_r1x5[100];
  Float_t photon_r2x5[100];
  Float_t photon_maxEnergyXtal[100];

  Float_t photon_SigmaEtaEta[100];
  Float_t photon_SigmaPhiPhi[100];
  Float_t photon_SigmaIetaIeta[100];
  Float_t photon_SigmaIphiIphi[100];
  Float_t photon_Roundness[100];
  Float_t photon_Angle[100];

  Float_t photon_ecalRecHitSumEtConeDR03[100];
  Float_t photon_hcalTowerSumEtConeDR03[100];
  Float_t photon_hcalDepth1TowerSumEtConeDR03[100];
  Float_t photon_hcalDepth2TowerSumEtConeDR03[100];
  Float_t photon_trkSumPtSolidConeDR03[100];
  Float_t photon_trkSumPtHollowConeDR03[100];
  Int_t   photon_nTrkSolidConeDR03[100];
  Int_t   photon_nTrkHollowConeDR03[100];
  Float_t photon_ecalRecHitSumEtConeDR04[100];
  Float_t photon_hcalTowerSumEtConeDR04[100];
  Float_t photon_hcalDepth1TowerSumEtConeDR04[100];
  Float_t photon_hcalDepth2TowerSumEtConeDR04[100];
  Float_t photon_trkSumPtSolidConeDR04[100];
  Float_t photon_trkSumPtHollowConeDR04[100];
  Int_t   photon_nTrkSolidConeDR04[100];
  Int_t   photon_nTrkHollowConeDR04[100];
  Bool_t  photon_hasPixelSeed[100];
  Float_t photon_HoE[100]; 
  Bool_t  photon_isEB[100];   
  Bool_t  photon_isEE[100];   
  Bool_t  photon_isEBGap[100];
  Bool_t  photon_isEEGap[100];
  Bool_t  photon_isEBEEGap[100];   
  Float_t photon_sc_energy[100];  
  Float_t photon_sc_eta[100];   
  Float_t photon_sc_phi[100];  
  Float_t photon_sc_x[100];   
  Float_t photon_sc_y[100];   
  Float_t photon_sc_z[100];    
  Float_t photon_etaWidth[100];   
  Float_t photon_phiWidth[100];   
  Float_t photon_sc_et[100];   
  Bool_t   photon_isconverted[100];   
  Float_t  photon_pairInvmass[100];   
  Float_t  photon_pairCotThetaSeparation[100];   
  Float_t  photon_pairmomentumX[100];   
  Float_t  photon_pairmomentumY[100];   
  Float_t  photon_pairmomentumZ[100];   
  Float_t  photon_EoverP[100];   
  Float_t  photon_ConvVx[100];  
  Float_t  photon_ConvVy[100];  
  Float_t  photon_ConvVz[100];  
  Float_t  photon_ZOfPrimaryVertex[100];  
  Float_t  photon_distOfMinimumApproach[100];  
  Float_t  photon_dPhiTracksAtVtx[100];   
  Float_t  photon_dPhiTracksAtEcal[100];  
  Float_t  photon_dEtaTracksAtEcal[100];   
  Int_t    photon_ncrys[100];   
  Float_t  photon_timing_xtal[100][2];   
  //Float_t  photon_timingavg_xtal[100];   
  //Float_t  photon_energy_xtal[100][100]; 
  //Float_t  photon_ieta_xtalEB[100][100]; 
  //Float_t  photon_iphi_xtalEB[100][100]; 
  Float_t  photon_s9[100]; 
  Int_t photon_basiccluster_size[100];
  Int_t photon_ntracks[100];
  Float_t photon_deltaR_mc[100];

  //JET 
  Int_t   jet_n;
  Int_t   pfjet_n;

 //ELECTRON 
  Int_t elec_n;
  Float_t elec_px[100];
  Float_t elec_py[100];
  Float_t elec_pz[100];
  Float_t elec_pt[100];
  Float_t elec_eta[100];
  Float_t elec_phi[100];
  Float_t elec_energy[100];
  Float_t elec_charge[100];
  Float_t elec_trkIso[100];
  Float_t elec_ecalIso[100];
  Float_t elec_hcalIso[100];
  Float_t elec_HoE[100]; 
  Float_t elec_SigmaIetaIeta[100];
  Float_t elec_SigmaIphiIphi[100];
  Float_t elec_dEtaIn[100];
  Float_t elec_dPhiIn[100];
  Float_t elec_sc_energy[100];
  Float_t elec_sc_eta[100];
  Float_t elec_sc_phi[100];


  //MUON
  Int_t    muon_n;  
  Float_t  muon_px[100];   
  Float_t  muon_py[100];   
  Float_t  muon_pz[100];   
  Float_t  muon_pt[100];   
  Float_t  muon_eta[100];   
  Float_t  muon_phi[100];   
  Float_t  muon_energy[100];   
  Float_t  muon_charge[100];   
  Bool_t   muon_isGlobalMuon[100];   
  Bool_t   muon_isTrackerMuon[100];   
  Bool_t   muon_isStandAloneMuon[100];   
  Bool_t   muon_InnerTrack_isNonnull[100];   
  Bool_t   muon_OuterTrack_isNonnull[100];   
  Float_t  muon_OuterTrack_InnerPoint_x[100];   
  Float_t  muon_OuterTrack_InnerPoint_y[100];   
  Float_t  muon_OuterTrack_InnerPoint_z[100];   
  Float_t  muon_OuterTrack_InnerPoint_px[100];   
  Float_t  muon_OuterTrack_InnerPoint_py[100];   
  Float_t  muon_OuterTrack_InnerPoint_pz[100];   
  Float_t  muon_OuterTrack_OuterPoint_x[100];   
  Float_t  muon_OuterTrack_OuterPoint_y[100];   
  Float_t  muon_OuterTrack_OuterPoint_z[100];   
  Float_t  muon_OuterTrack_OuterPoint_px[100];   
  Float_t  muon_OuterTrack_OuterPoint_py[100];   
  Float_t  muon_OuterTrack_OuterPoint_pz[100];   
  Float_t  muon_InnerTrack_InnerPoint_x[100];   
  Float_t  muon_InnerTrack_InnerPoint_y[100];   
  Float_t  muon_InnerTrack_InnerPoint_z[100];   
  Float_t  muon_InnerTrack_InnerPoint_px[100];   
  Float_t  muon_InnerTrack_InnerPoint_py[100];   
  Float_t  muon_InnerTrack_InnerPoint_pz[100];   
  Float_t  muon_InnerTrack_OuterPoint_x[100];   
  Float_t  muon_InnerTrack_OuterPoint_y[100];   
  Float_t  muon_InnerTrack_OuterPoint_z[100];   
  Float_t  muon_InnerTrack_OuterPoint_px[100];   
  Float_t  muon_InnerTrack_OuterPoint_py[100];   
  Float_t  muon_InnerTrack_OuterPoint_pz[100];   
  Float_t  muon_OuterPoint_x[100];   
  Float_t  muon_OuterPoint_y[100]; 
  Float_t  muon_OuterPoint_z[100];   
  Float_t  muon_InnerPoint_x[100];   
  Float_t  muon_InnerPoint_y[100]; 
  Float_t  muon_InnerPoint_z[100];   
  Float_t  muon_vx[100];   
  Float_t  muon_vy[100]; 
  Float_t  muon_vz[100];


  //cosmicmuon
  Int_t    cosmicmuon_n;  
  Float_t  cosmicmuon_px[100];   
  Float_t  cosmicmuon_py[100];   
  Float_t  cosmicmuon_pz[100];   
  Float_t  cosmicmuon_pt[100];   
  Float_t  cosmicmuon_eta[100];   
  Float_t  cosmicmuon_phi[100];   
  Float_t  cosmicmuon_energy[100];   
  Float_t  cosmicmuon_charge[100];   
  Bool_t   cosmicmuon_isGlobalcosmicmuon[100];   
  Bool_t   cosmicmuon_isTrackercosmicmuon[100];   
  Bool_t   cosmicmuon_isStandAlonecosmicmuon[100];   
  Bool_t   cosmicmuon_InnerTrack_isNonnull[100];   
  Bool_t   cosmicmuon_OuterTrack_isNonnull[100];   
  Float_t  cosmicmuon_OuterTrack_InnerPoint_x[100];   
  Float_t  cosmicmuon_OuterTrack_InnerPoint_y[100];   
  Float_t  cosmicmuon_OuterTrack_InnerPoint_z[100];   
  Float_t  cosmicmuon_OuterTrack_InnerPoint_px[100];   
  Float_t  cosmicmuon_OuterTrack_InnerPoint_py[100];   
  Float_t  cosmicmuon_OuterTrack_InnerPoint_pz[100];   
  Float_t  cosmicmuon_OuterTrack_OuterPoint_x[100];   
  Float_t  cosmicmuon_OuterTrack_OuterPoint_y[100];   
  Float_t  cosmicmuon_OuterTrack_OuterPoint_z[100];   
  Float_t  cosmicmuon_OuterTrack_OuterPoint_px[100];   
  Float_t  cosmicmuon_OuterTrack_OuterPoint_py[100];   
  Float_t  cosmicmuon_OuterTrack_OuterPoint_pz[100];   
  Float_t  cosmicmuon_InnerTrack_InnerPoint_x[100];   
  Float_t  cosmicmuon_InnerTrack_InnerPoint_y[100];   
  Float_t  cosmicmuon_InnerTrack_InnerPoint_z[100];   
  Float_t  cosmicmuon_InnerTrack_InnerPoint_px[100];   
  Float_t  cosmicmuon_InnerTrack_InnerPoint_py[100];   
  Float_t  cosmicmuon_InnerTrack_InnerPoint_pz[100];   
  Float_t  cosmicmuon_InnerTrack_OuterPoint_x[100];   
  Float_t  cosmicmuon_InnerTrack_OuterPoint_y[100];   
  Float_t  cosmicmuon_InnerTrack_OuterPoint_z[100];   
  Float_t  cosmicmuon_InnerTrack_OuterPoint_px[100];   
  Float_t  cosmicmuon_InnerTrack_OuterPoint_py[100];   
  Float_t  cosmicmuon_InnerTrack_OuterPoint_pz[100];  
  Float_t  cosmicmuon_OuterPoint_x[100];   
  Float_t  cosmicmuon_OuterPoint_y[100]; 
  Float_t  cosmicmuon_OuterPoint_z[100];   

  //MC photon_
  Float_t MC_photon_pt[1000];
  Float_t MC_photon_eta[1000];
  Float_t MC_photon_phi[1000];
  Float_t MC_photon_px[1000]; 
  Float_t MC_photon_py[1000];
  Float_t MC_photon_pz[1000];
  Float_t MC_photon_E[1000];
  Int_t   MC_photon_MotherID[1000];
  Float_t MC_photon_MotherPt[1000];
  Float_t MC_photon_MotherEta[1000];
  Float_t MC_photon_MotherPhi[1000];
  Int_t   MC_photon_GrandMotherID[1000]; 

  ////////////////////////added by bhawna//////////////// 
  ///////////////////////////////////////////////////////
  //tau info
  Int_t   tau_n;
  Int_t   OneProng0Pi0[100];  
  Int_t   OneProng1Pi0[100];  
  Int_t   OneProng2Pi0[100];  
  Int_t   ThreeProng0Pi0[100];
  Int_t   ThreeProng1Pi0[100];
  Float_t  GenHadTauPt[100];  
  Float_t  GenHadTauEta[100]; 
  Float_t  GenHadTauPhi[100]; 
  Int_t     NPions[100];  
  Int_t     pionPdgId[100][100]; 
  Float_t  pionPt[100][100];  
  Float_t  pionEta[100][100]; 
  Float_t  pionPhi[100][100]; 
  Int_t     NPi0[100];   
  Int_t     pi0PdgId[100][100]; 
  Float_t  pi0Pt[100][100];  
  Float_t  pi0Eta[100][100]; 
  Float_t  pi0Phi[100][100]; 
  Int_t     NPhotons[100];  
  Float_t  photonPt[100][100]; //change these doubles to float after the next time when you will ran since this  is wrong
  Float_t  photonEta[100][100];
  Float_t  photonPhi[100][100];
  Int_t     photonPdgId[100][100];
  ////////////////ends here//////////////////////
  ////////////////////////////////////////////////


  //MET
  Float_t calo_MetSigma;
  //Float_t calo_MetCorr;
  Float_t calo_MetPt[6];
  Float_t calo_MetPx[6];
  Float_t calo_MetPy[6];
  Float_t calo_MetEz;
  Float_t calo_MetPhi[6];
  Float_t calo_MetSumEt[6];
  Float_t calo_EtFractionHadronic;
  Float_t calo_EmEtFraction;
  Float_t calo_HadEtInHB;
  Float_t calo_HadEtInHE;
  Float_t calo_HadEtInHO;
  Float_t calo_HadEtInHF;
  Float_t calo_EmEtInEB;
  Float_t calo_EmEtInEE;
  Float_t calo_EmEtInHF;
  Float_t calo_MaxEtInEmTowers;
  Float_t calo_MaxEtInHadTowers;
 
  Float_t delta_phi;
  Float_t delta_phiGEN;
  Float_t pf_MetPt[6];
  Float_t pf_MetPx[6];
  Float_t pf_MetPy[6];
  Float_t pf_MetPhi[6];
  Float_t pf_MetSumEt[6];   
  Float_t delta_phiPF;  
  Float_t tc_MetPt[6];
  Float_t tc_MetPx[6];
  Float_t tc_MetPy[6];
  Float_t tc_MetPhi[6];
  Float_t tc_MetSumEt[6];
  Float_t delta_phiTC;   

  //GEN MET
  Float_t gen_MetPt;
  Float_t gen_MetPx;
  Float_t gen_MetPy;
  Float_t gen_MetPhi;
  Float_t gen_MetSumEt;

  //TRIGGER
  vector<std::string>    HLTNames;
  vector<int>            HLTPrescales;	
  vector<bool>           HLTIfPassed;
  vector<std::string>    *HLT_Names;
  vector<int>            *HLT_Prescales;	
  vector<bool>           *HLT_IfPassed;
 

  //Graviton
  Float_t  graviton_pt;

  //MC EVENT TYPE
  bool isZ_event, isW_event;
  bool isZnunu_event, isZelec_event, isZmu_event, isZtau_event;
  bool isWelec_event, isWmu_event, isWtau_event;
  bool isSingleHardPhoton_event;
  bool isdiphoton_event;
  bool isisr_photon_event;

  //CSC variables
  Int_t cscseg_n;
  Float_t cscseg_x[100];
  Float_t cscseg_y[100];
  Float_t cscseg_z[100];
  Float_t cscseg_time[100];

 //BeamHaloSummary                                                                                                              
  bool is_BeamHaloIDTightPass;
  bool is_BeamHaloIDLoosePass;
         
  bool is_BeamHaloEcalLoosePass;
  bool is_BeamHaloHcalLoosePass;
  bool is_BeamHaloCSCLoosePass;
  bool is_BeamHaloGlobalLoosePass;
         
  bool is_BeamHaloEcalTightPass;
  bool is_BeamHaloHcalTightPass;
  bool is_BeamHaloCSCTightPass;
  bool is_BeamHaloGlobalTightPass;

  bool is_SmellsLikeHalo_Tag;
  bool is_LooseHalo_Tag;
  bool is_TightHalo_Tag;
  bool is_ExtremeTightHalo_Tag;

 
  //RPChit info



  float RPChit_x[10000];
  float RPChit_y[10000];


  //calculated here in this code:
  Float_t jet_pt[100];
  Float_t jet_px[100];
  Float_t jet_py[100];
  Float_t jet_pz[100];
  Float_t jet_E[100];
  Float_t jet_eta[100];
  Float_t jet_phi[100];
  Float_t jet_EMenergyFraction[100];
  Float_t jet_energyFractionHadronic[100];
  Int_t   jet_hitsInN90[100];
  Int_t   jet_n90Hits[100];
  Int_t   jet_nTowers[100];
  Float_t jet_fHPD[100];
  Float_t jet_fRBX[100];
  Float_t jet_RHF[100];
  Int_t   jet_n_PhotonRemoved;
  Int_t   jet_n_FromHT;


  //pfjets information
  Float_t pfjet_pt[100];
  Float_t pfjet_px[100];
  Float_t pfjet_py[100];
  Float_t pfjet_pz[100];
  Float_t pfjet_E[100];
  Float_t pfjet_eta[100];
  Float_t pfjet_phi[100];
  Int_t   pfjet_n_PhotonRemoved;
  Int_t   pfjet_n_FromHT;



  Float_t HT;
  Float_t HT_JetpT50;
  Float_t HT_JetpT100;
  Float_t HT_JetpT150;


  Float_t pfHT;
  Float_t HT_pfJetpT50;
  Float_t HT_pfJetpT100;
  Float_t HT_pfJetpT150;

  bool scraping_isScrapingEvent;
  bool IsHEHalo[100];
  bool IsCSCHalo[100];
  bool IsTrackHalo[100];
  bool IsE2E9Spike[100];
  bool IsTimeSpike[100];
  Float_t DeltaPhiCSCHalo[100][100];
  //For cosmic
  Float_t sigmaRC[100];
  Float_t fisherd2[100];
  Float_t fisherd3[100];
};

myCuts::myCuts(){
  fchain = new TChain("myEvent");
  TString location = "dcap://cmsgridftp.fnal.gov:24125";
  TString path = "/pnfs/cms/WAX/11/store/user/sushil/MonoPhoton/38X_Ntuples_V22/Data/Wgamma/"; 
  TString short_path ="/pnfs/fnal.gov/usr/cms/WAX/11/store/user/sushil/MonoPhoton/38X_Ntuples_V22/Data/Wgamma/"; 

  TSystemDirectory sourceDir("hi",path);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter next(fileList);
  TSystemFile* fileName;
  int fileNumber = 1;
  int maxFiles = -1;
   
  while ((fileName = (TSystemFile*)next()) && fileNumber >  maxFiles){
    if(TString(fileName->GetName()) == "." || TString(fileName->GetName()) == ".."  ){continue;}
    cout<<" File Name = "<<(fileName->GetName())<<endl;
    TString  FullPathInputFile = (location+short_path+fileName->GetName());
    cout<<FullPathInputFile<<endl;
    fchain->Add(FullPathInputFile); cout<<fchain->GetEntries()<<endl;
  }//loop over while
  Init(fchain);
  cout<<"ADDED files:"<<endl;
  define_Tree();
}


myCuts::~myCuts(){
}

void myCuts::define_Tree(){

  HLT_Names        = &HLTNames;
  HLT_Prescales    = &HLTPrescales;
  HLT_IfPassed     = &HLTIfPassed;	



  Long64_t maxsize = 500000000;                  //100GB
  maxsize *= 100;                                //to bypass some compiler limitations with big constants
  TTree::SetMaxTreeSize(maxsize);
  outfile = new TFile("/uscmst1b_scratch/lpc1/3DayLifetime/sushil/ShortRootFileCode/Wgamma.root","RECREATE");	
  outfile->cd();
  
  //TREE NAME
  tree   = new TTree("tree","a tree with histograms");
  tree->Branch("NFilled",&NFilled,"NFilled/I");
  tree->Branch("NJFilled",&NJFilled,"NJFilled/I");

  //Event infor
  tree->Branch("Nnevents",&Nnevents,"Nnevents/I");
  tree->Branch("Nrun",&Nrun,"Nrun/i");
  tree->Branch("Nevent",&Nevent,"Nevent/i");
  tree->Branch("NluminosityBlock",&NluminosityBlock,"NluminosityBlock/i");
  tree->Branch("NbeamCrossing",&NbeamCrossing,"NbeamCrossing/i");
  tree->Branch("NtotalIntensityBeam1",&NtotalIntensityBeam1,"NtotalIntensityBeam1/i");
  tree->Branch("NtotalIntensityBeam2",&NtotalIntensityBeam2,"NtotalIntensityBeam2/i");
  tree->Branch("NavgInsDelLumi",&NavgInsDelLumi,"NavgInsDelLumi/F");
  tree->Branch("NavgInsDelLumiErr",&NavgInsDelLumiErr,"NavgInsDelLumiErr/F");
  tree->Branch("NavgInsRecLumi",&NavgInsRecLumi,"NavgInsRecLumi/F");
  tree->Branch("NavgInsRecLumiErr",&NavgInsRecLumiErr,"NavgInsRecLumiErr/F");
  
  tree->Branch("pthat",&pthat,"pthat/F");

  //HLT
  tree->Branch("HLT_Names","vector<std::string>",&HLT_Names);
  tree->Branch("HLT_Prescales","vector<int>",&HLT_Prescales);
  tree->Branch("HLT_IfPassed","vector<bool>",&HLT_IfPassed);

  //VERTICES
  tree->Branch("vertex_n",&vertex_n,"vertex_n/I");
  tree->Branch("vertex_x",vertex_x,"vertex_x[NFilled]/F");
  tree->Branch("vertex_y",vertex_y,"vertex_y[NFilled]/F");
  tree->Branch("vertex_z",vertex_z,"vertex_z[NFilled]/F");
  tree->Branch("vertex_tracksize",vertex_tracksize,"vertex_tracksize[NFilled]/F");
  tree->Branch("vertex_ndof",vertex_ndof,"vertex_ndof[NFilled]/F");
  tree->Branch("vertex_chi2",vertex_chi2,"vertex_chi2[NFilled]/F");
  tree->Branch("vertex_d0",vertex_d0,"vertex_d0[NFilled]/F");
  tree->Branch("vertex_isFake",vertex_isFake,"vertex_isFake[NFilled]/O");
  
  //TRACK
  tree->Branch("track_n",&track_n,"track_n/I");
  tree->Branch("track_px",track_px,"track_px[NFilled]/F");
  tree->Branch("track_py",track_py,"track_py[NFilled]/F");
  tree->Branch("track_pz",track_pz,"track_pz[NFilled]/F");
  tree->Branch("track_pt",track_pt,"track_pt[NFilled]/F");
  tree->Branch("track_eta",track_eta,"track_eta[NFilled]/F");
  tree->Branch("track_phi",track_phi,"track_phi[NFilled]/F");
  
  //PHOTON
  tree->Branch("photon_n",&Photon_n,"photon_n/I");
  tree->Branch("photon_E",photon_E,"photon_E[NFilled]/F");
  //tree->Branch("photon_iLICTD",photon_iLICTD,"photon_iLICTD[NFilled]/F");
  tree->Branch("photon_pt",photon_pt,"photon_pt[NFilled]/F");
  tree->Branch("photon_eta",photon_eta,"photon_eta[NFilled]/F");
  tree->Branch("photon_phi",photon_phi,"photon_phi[NFilled]/F");
  tree->Branch("photon_theta",photon_theta,"photon_theta[NFilled]/F");
  tree->Branch("photon_et",photon_et,"photon_et[NFilled]/F");
  tree->Branch("photon_swissCross",photon_swissCross,"photon_swissCross[NFilled]/F");
  tree->Branch("photon_r9",photon_r9,"photon_r9[NFilled]/F");
  tree->Branch("photon_e1x5",photon_e1x5,"photon_e1x5[NFilled]/F");
  tree->Branch("photon_e2x5",photon_e2x5,"photon_e2x5[NFilled]/F");
  tree->Branch("photon_e3x3",photon_e3x3,"photon_e3x3[NFilled]/F");
  tree->Branch("photon_e5x5",photon_e5x5,"photon_e5x5[NFilled]/F");
  tree->Branch("photon_r1x5",photon_r1x5,"photon_erx5[NFilled]/F");
  tree->Branch("photon_r2x5",photon_r2x5,"photon_erx5[NFilled]/F");
  tree->Branch("photon_maxEnergyXtal",photon_maxEnergyXtal,"photon_maxEnergyXtal[NFilled]/F");
  tree->Branch("photon_SigmaEtaEta",photon_SigmaEtaEta,"photon_SigmaEtaEta[NFilled]/F");
  tree->Branch("photon_SigmaIetaIeta",photon_SigmaIetaIeta,"photon_SigmaIetaIeta[NFilled]/F");
  tree->Branch("photon_SigmaPhiPhi",photon_SigmaPhiPhi,"photon_SigmaPhiPhi[NFilled]/F");
  tree->Branch("photon_SigmaIphiIphi",photon_SigmaIphiIphi,"photon_SigmaIphiIphi[NFilled]/F");
  tree->Branch("photon_Roundness",photon_Roundness,"photon_Roundness[NFilled]/F");
  tree->Branch("photon_Angle",photon_Angle,"photon_Angle[NFilled]/F");
  tree->Branch("photon_ecalRecHitSumEtConeDR03",photon_ecalRecHitSumEtConeDR03,"photon_ecalRecHitSumEtConeDR03[NFilled]/F");
  tree->Branch("photon_hcalTowerSumEtConeDR03",photon_hcalTowerSumEtConeDR03,"photon_hcalTowerSumEtConeDR03[NFilled]/F");
  tree->Branch("photon_trkSumPtSolidConeDR03",photon_trkSumPtSolidConeDR03,"photon_trkSumPtSolidConeDR03[NFilled]/F");
  tree->Branch("photon_trkSumPtHollowConeDR03",photon_trkSumPtHollowConeDR03,"photon_trkSumPtHollowConeDR03[NFilled]/F");
  tree->Branch("photon_nTrkSolidConeDR03",photon_nTrkSolidConeDR03,"photon_nTrkSolidConeDR03[NFilled]/I");
  tree->Branch("photon_nTrkHollowConeDR03",photon_nTrkHollowConeDR03,"photon_nTrkHollowConeDR03[NFilled]/I");
  tree->Branch("photon_hcalDepth1TowerSumEtConeDR03",photon_hcalDepth1TowerSumEtConeDR03,"photon_hcalDepth1TowerSumEtConeDR03[NFilled]/F");
  tree->Branch("photon_hcalDepth2TowerSumEtConeDR03",photon_hcalDepth2TowerSumEtConeDR03,"photon_hcalDepth2TowerSumEtConeDR03[NFilled]/F");
  tree->Branch("photon_ecalRecHitSumEtConeDR04",photon_ecalRecHitSumEtConeDR04,"photon_ecalRecHitSumEtConeDR04[NFilled]/F");
  tree->Branch("photon_hcalTowerSumEtConeDR04",photon_hcalTowerSumEtConeDR04,"photon_hcalTowerSumEtConeDR04[NFilled]/F");
  tree->Branch("photon_trkSumPtSolidConeDR04",photon_trkSumPtSolidConeDR04,"photon_trkSumPtSolidConeDR04[NFilled]/F");
  tree->Branch("photon_trkSumPtHollowConeDR04",photon_trkSumPtHollowConeDR04,"photon_trkSumPtHollowConeDR04[NFilled]/F");
  tree->Branch("photon_nTrkSolidConeDR04",photon_nTrkSolidConeDR04,"photon_nTrkSolidConeDR04[NFilled]/I");
  tree->Branch("photon_nTrkHollowConeDR04",photon_nTrkHollowConeDR04,"photon_nTrkHollowConeDR04[NFilled]/I");
  tree->Branch("photon_hcalDepth1TowerSumEtConeDR04",photon_hcalDepth1TowerSumEtConeDR04,"photon_hcalDepth1TowerSumEtConeDR04[NFilled]/F");
  tree->Branch("photon_hcalDepth2TowerSumEtConeDR04",photon_hcalDepth2TowerSumEtConeDR04,"photon_hcalDepth2TowerSumEtConeDR04[NFilled]/F");
  tree->Branch("photon_hasPixelSeed",photon_hasPixelSeed,"photon_hasPixelSeed[NFilled]/O");
  tree->Branch("photon_isEB",photon_isEB,"photon_isEB[NFilled]/O");
  tree->Branch("photon_isEE",photon_isEE,"photon_isEE[NFilled]/O");
  tree->Branch("photon_isEBGap",photon_isEBGap,"photon_isEBGap[NFilled]/O");
  tree->Branch("photon_isEEGap",photon_isEEGap,"photon_isEEGap[NFilled]/O");
  tree->Branch("photon_isEBEEGap",photon_isEBEEGap,"photon_isEBEEGap[NFilled]/O");
  tree->Branch("photon_HoE",photon_HoE,"photon_HoE[NFilled]/F");
  tree->Branch("photon_px",photon_px,"photon_px[NFilled]/F");
  tree->Branch("photon_py",photon_py,"photon_py[NFilled]/F");
  tree->Branch("photon_pz",photon_pz,"photon_pz[NFilled]/F");
  tree->Branch("photon_basiccluster_size",photon_basiccluster_size,"photon_basiccluster_size[NFilled]/I");
  tree->Branch("photon_sc_energy",photon_sc_energy,"photon_sc_energy[NFilled]/F");
  tree->Branch("photon_sc_eta",photon_sc_eta,"photon_sc_eta[NFilled]/F");
  tree->Branch("photon_sc_phi",photon_sc_phi,"photon_sc_phi[NFilled]/F");
  tree->Branch("photon_sc_x",photon_sc_x,"photon_sc_x[NFilled]/F");
  tree->Branch("photon_sc_y",photon_sc_y,"photon_sc_y[NFilled]/F");
  tree->Branch("photon_sc_z",photon_sc_z,"photon_sc_z[NFilled]/F");
  tree->Branch("photon_etaWidth",photon_etaWidth,"photon_etaWidth[NFilled]/F");
  tree->Branch("photon_phiWidth",photon_phiWidth,"photon_phiWidth[NFilled]/F");
  tree->Branch("photon_sc_et",photon_sc_et,"photon_sc_et[NFilled]/F");
  tree->Branch("photon_ntracks",photon_ntracks,"photon_ntracks[NFilled]/I");
  tree->Branch("photon_isconverted",photon_isconverted,"photon_isconverted[NFilled]/O");
  tree->Branch("photon_pairInvarmass",photon_pairInvmass,"photon_pairInvmass[NFilled]/F");
  tree->Branch("photon_pairCotThetaSeparation",photon_pairCotThetaSeparation,"photon_pairCotThetaSeparation[NFilled]/F");
  tree->Branch("photon_pairmomentumX",photon_pairmomentumX,"photon_pairmomentumX[NFilled]/F");
  tree->Branch("photon_pairmomentumY",photon_pairmomentumY,"photon_pairmomentumY[NFilled]/F");
  tree->Branch("photon_pairmomentumZ",photon_pairmomentumZ,"photon_pairmomentumZ[NFilled]/F");
  tree->Branch("photon_EoverP",photon_EoverP,"photon_EoverP[NFilled]/F");
  tree->Branch("photon_ConvVx",photon_ConvVx,"photon_ConvVx[NFilled]/F");
  tree->Branch("photon_ConvVy",photon_ConvVy,"photon_ConvVy[NFilled]/F");
  tree->Branch("photon_ConvVz",photon_ConvVz,"photon_ConvVz[NFilled]/F");
  tree->Branch("photon_ZOfPrimaryVertex",photon_ZOfPrimaryVertex,"photon_ZOfPrimaryVertex[NFilled]/F");
  tree->Branch("photon_distOfMinimumApproach",photon_distOfMinimumApproach,"photon_distOfMinimumApproach[NFilled]/F");
  tree->Branch("photon_dPhiTracksAtVtx",photon_dPhiTracksAtVtx,"photon_dPhiTracksAtVtx[NFilled]/F");
  tree->Branch("photon_dPhiTracksAtEcal",photon_dPhiTracksAtEcal,"photon_dPhiTracksAtEcal[NFilled]/F");
  tree->Branch("photon_dEtaTracksAtEcal",photon_dEtaTracksAtEcal,"photon_dEtaTracksAtEcal[NFilled]/F");
  tree->Branch("photon_e2e9",photon_e2e9,"photon_e2e9[NFilled]/F");
  tree->Branch("photon_e4e1",photon_e4e1,"photon_e4e1[NFilled]/F");
  tree->Branch("photon_e6e2",photon_e6e2,"photon_e6e2[NFilled]/F");
  tree->Branch("photon_timing_xtal",photon_timing_xtal,"photon_timing_xtal[NFilled][2]/F");

  //
  tree->Branch("photon_deltaR_mc",photon_deltaR_mc,"photon_deltaR_mc[NFilled]/F"); 
 
  //JET
  tree->Branch("jet_n",&jet_n,"jet_n/I");
  tree->Branch("pfjet_n",&pfjet_n,"pfjet_n/I");

  //ELECTRON 
  tree->Branch("elec_n",&elec_n,"elec_n/I");
  tree->Branch("elec_px",elec_px,"elec_px[NFilled]/F");
  tree->Branch("elec_py",elec_py,"elec_py[NFilled]/F");
  tree->Branch("elec_pz",elec_pz,"elec_pz[NFilled]/F");
  tree->Branch("elec_pt",elec_pt,"elec_pt[NFilled]/F");
  tree->Branch("elec_eta",elec_eta,"elec_eta[NFilled]/F");
  tree->Branch("elec_phi",elec_phi,"elec_phi[NFilled]/F");
  tree->Branch("elec_energy",elec_energy,"elec_energy[NFilled]/F");
  tree->Branch("elec_charge",elec_charge,"elec_charge[NFilled]/F");
  tree->Branch("elec_trkIso",elec_trkIso,"elec_trkIso[NFilled]/F");
  tree->Branch("elec_ecalIso",elec_ecalIso,"elec_ecalIso[NFilled]/F");
  tree->Branch("elec_hcalIso",elec_hcalIso,"elec_hcalIso[NFilled]/F");
  tree->Branch("elec_HoE",elec_HoE,"elec_HoE[NFilled]/F");
  tree->Branch("elec_SigmaIetaIeta",elec_SigmaIetaIeta,"elec_SigmaIetaIeta[NFilled]/F");
  tree->Branch("elec_dEtaIn",elec_dEtaIn,"elec_dEtaIn[NFilled]/F");	
  tree->Branch("elec_dPhiIn",elec_dPhiIn,"elec_dPhiIn[NFilled]/F");
  tree->Branch("elec_sc_energy",elec_sc_energy,"elec_sc_energy[NFilled]/F");	
  tree->Branch("elec_sc_eta",elec_sc_eta,"elec_sc_eta[NFilled]/F");	
  tree->Branch("elec_sc_phi",elec_sc_phi,"elec_sc_phi[NFilled]/F");	
  
  //MUON 
  tree->Branch("muon_n",&muon_n,"muon_n/I"); 
  tree->Branch("muon_px",muon_px,"muon_px[NFilled]/F");
  tree->Branch("muon_py",muon_py,"muon_py[NFilled]/F");
  tree->Branch("muon_pz",muon_pz,"muon_pz[NFilled]/F");
  tree->Branch("muon_pt",muon_pt,"muon_pt[NFilled]/F");
  tree->Branch("muon_eta",muon_eta,"muon_eta[NFilled]/F");
  tree->Branch("muon_phi",muon_phi,"muon_phi[NFilled]/F");
  tree->Branch("muon_energy",muon_energy,"muon_energy[NFilled]/F");
  tree->Branch("muon_charge",muon_charge,"muon_charge[NFilled]/F");
  tree->Branch("muon_isGlobalMuon",muon_isGlobalMuon,"muon_isGlobalMuon[NFilled]/O");
  tree->Branch("muon_isTrackerMuon",muon_isTrackerMuon,"muon_isTrackerMuon[NFilled]/O");
  tree->Branch("muon_isStandAloneMuon",muon_isStandAloneMuon,"muon_isStandAloneMuon[NFilled]/O");
  tree->Branch("muon_InnerTrack_isNonnull",muon_InnerTrack_isNonnull,"muon_InnerTrack_isNonnull[NFilled]/O");
  tree->Branch("muon_OuterTrack_isNonnull",muon_OuterTrack_isNonnull,"muon_OuterTrack_isNonnull[NFilled]/O");
  tree->Branch("muon_OuterTrack_InnerPoint_x",muon_OuterTrack_InnerPoint_x,"muon_OuterTrack_InnerPoint_x[NFilled]/F");
  tree->Branch("muon_OuterTrack_InnerPoint_y",muon_OuterTrack_InnerPoint_y,"muon_OuterTrack_InnerPoint_y[NFilled]/F");
  tree->Branch("muon_OuterTrack_InnerPoint_z",muon_OuterTrack_InnerPoint_z,"muon_OuterTrack_InnerPoint_z[NFilled]/F");
  tree->Branch("muon_OuterTrack_InnerPoint_px",muon_OuterTrack_InnerPoint_px,"muon_OuterTrack_InnerPoint_px[NFilled]/F");
  tree->Branch("muon_OuterTrack_InnerPoint_py",muon_OuterTrack_InnerPoint_py,"muon_OuterTrack_InnerPoint_py[NFilled]/F");
  tree->Branch("muon_OuterTrack_InnerPoint_pz",muon_OuterTrack_InnerPoint_pz,"muon_OuterTrack_InnerPoint_pz[NFilled]/F");
  tree->Branch("muon_OuterTrack_OuterPoint_x",muon_OuterTrack_OuterPoint_x,"muon_OuterTrack_OuterPoint_x[NFilled]/F");
  tree->Branch("muon_OuterTrack_OuterPoint_y",muon_OuterTrack_OuterPoint_y,"muon_OuterTrack_OuterPoint_y[NFilled]/F");
  tree->Branch("muon_OuterTrack_OuterPoint_z",muon_OuterTrack_OuterPoint_z,"muon_OuterTrack_OuterPoint_z[NFilled]/F");
  tree->Branch("muon_OuterTrack_OuterPoint_px",muon_OuterTrack_OuterPoint_px,"muon_OuterTrack_OuterPoint_px[NFilled]/F");
  tree->Branch("muon_OuterTrack_OuterPoint_py",muon_OuterTrack_OuterPoint_py,"muon_OuterTrack_OuterPoint_py[NFilled]/F");
  tree->Branch("muon_OuterTrack_OuterPoint_pz",muon_OuterTrack_OuterPoint_pz,"muon_OuterTrack_OuterPoint_pz[NFilled]/F");
  tree->Branch("muon_InnerTrack_InnerPoint_x",muon_InnerTrack_InnerPoint_x,"muon_InnerTrack_InnerPoint_x[NFilled]/F");
  tree->Branch("muon_InnerTrack_InnerPoint_y",muon_InnerTrack_InnerPoint_y,"muon_InnerTrack_InnerPoint_y[NFilled]/F");
  tree->Branch("muon_InnerTrack_InnerPoint_z",muon_InnerTrack_InnerPoint_z,"muon_InnerTrack_InnerPoint_z[NFilled]/F");
  tree->Branch("muon_InnerTrack_InnerPoint_px",muon_InnerTrack_InnerPoint_px,"muon_InnerTrack_InnerPoint_px[NFilled]/F");
  tree->Branch("muon_InnerTrack_InnerPoint_py",muon_InnerTrack_InnerPoint_py,"muon_InnerTrack_InnerPoint_py[NFilled]/F");
  tree->Branch("muon_InnerTrack_InnerPoint_pz",muon_InnerTrack_InnerPoint_pz,"muon_InnerTrack_InnerPoint_pz[NFilled]/F");
  tree->Branch("muon_InnerTrack_OuterPoint_x",muon_InnerTrack_OuterPoint_x,"muon_InnerTrack_OuterPoint_x[NFilled]/F");
  tree->Branch("muon_InnerTrack_OuterPoint_y",muon_InnerTrack_OuterPoint_y,"muon_InnerTrack_OuterPoint_y[NFilled]/F");
  tree->Branch("muon_InnerTrack_OuterPoint_z",muon_InnerTrack_OuterPoint_z,"muon_InnerTrack_OuterPoint_z[NFilled]/F");
  tree->Branch("muon_InnerTrack_OuterPoint_px",muon_InnerTrack_OuterPoint_px,"muon_InnerTrack_OuterPoint_px[NFilled]/F");
  tree->Branch("muon_InnerTrack_OuterPoint_py",muon_InnerTrack_OuterPoint_py,"muon_InnerTrack_OuterPoint_py[NFilled]/F");
  tree->Branch("muon_InnerTrack_OuterPoint_pz",muon_InnerTrack_OuterPoint_pz,"muon_InnerTrack_OuterPoint_pz[NFilled]/F");
  tree->Branch("muon_OuterPoint_x",muon_OuterPoint_x,"muon_OuterPoint_x[NFilled]/F");
  tree->Branch("muon_OuterPoint_y",muon_OuterPoint_y,"muon_OuterPoint_y[NFilled]/F");
  tree->Branch("muon_OuterPoint_z",muon_OuterPoint_z,"muon_OuterPoint_z[NFilled]/F"); 
  tree->Branch("muon_InnerPoint_x",muon_InnerPoint_x,"muon_InnerPoint_x[NFilled]/F");
  tree->Branch("muon_InnerPoint_y",muon_InnerPoint_y,"muon_InnerPoint_y[NFilled]/F");
  tree->Branch("muon_InnerPoint_z",muon_InnerPoint_z,"muon_InnerPoint_z[NFilled]/F"); 
  tree->Branch("muon_vx",muon_vx,"muon_vx[NFilled]/F");
  tree->Branch("muon_vy",muon_vy,"muon_vy[NFilled]/F");
  tree->Branch("muon_vz",muon_vz,"muon_vz[NFilled]/F"); 

  //cosmicmuon 
  tree->Branch("cosmicmuon_n",&cosmicmuon_n,"cosmicmuon_n/I"); 
  tree->Branch("cosmicmuon_px",cosmicmuon_px,"cosmicmuon_px[NFilled]/F");
  tree->Branch("cosmicmuon_py",cosmicmuon_py,"cosmicmuon_py[NFilled]/F");
  tree->Branch("cosmicmuon_pz",cosmicmuon_pz,"cosmicmuon_pz[NFilled]/F");
  tree->Branch("cosmicmuon_pt",cosmicmuon_pt,"cosmicmuon_pt[NFilled]/F");
  tree->Branch("cosmicmuon_eta",cosmicmuon_eta,"cosmicmuon_eta[NFilled]/F");
  tree->Branch("cosmicmuon_phi",cosmicmuon_phi,"cosmicmuon_phi[NFilled]/F");
  tree->Branch("cosmicmuon_energy",cosmicmuon_energy,"cosmicmuon_energy[NFilled]/F");
  tree->Branch("cosmicmuon_charge",cosmicmuon_charge,"cosmicmuon_charge[NFilled]/F");
  tree->Branch("cosmicmuon_isGlobalcosmicmuon",cosmicmuon_isGlobalcosmicmuon,"cosmicmuon_isGlobalcosmicmuon[NFilled]/O");
  tree->Branch("cosmicmuon_isTrackercosmicmuon",cosmicmuon_isTrackercosmicmuon,"cosmicmuon_isTrackercosmicmuon[NFilled]/O");
  tree->Branch("cosmicmuon_isStandAlonecosmicmuon",cosmicmuon_isStandAlonecosmicmuon,"cosmicmuon_isStandAlonecosmicmuon[NFilled]/O");
  tree->Branch("cosmicmuon_InnerTrack_isNonnull",cosmicmuon_InnerTrack_isNonnull,"cosmicmuon_InnerTrack_isNonnull[NFilled]/O");
  tree->Branch("cosmicmuon_OuterTrack_isNonnull",cosmicmuon_OuterTrack_isNonnull,"cosmicmuon_OuterTrack_isNonnull[NFilled]/O");
  tree->Branch("cosmicmuon_OuterTrack_InnerPoint_x",cosmicmuon_OuterTrack_InnerPoint_x,"cosmicmuon_OuterTrack_InnerPoint_x[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_InnerPoint_y",cosmicmuon_OuterTrack_InnerPoint_y,"cosmicmuon_OuterTrack_InnerPoint_y[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_InnerPoint_z",cosmicmuon_OuterTrack_InnerPoint_z,"cosmicmuon_OuterTrack_InnerPoint_z[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_InnerPoint_px",cosmicmuon_OuterTrack_InnerPoint_px,"cosmicmuon_OuterTrack_InnerPoint_px[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_InnerPoint_py",cosmicmuon_OuterTrack_InnerPoint_py,"cosmicmuon_OuterTrack_InnerPoint_py[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_InnerPoint_pz",cosmicmuon_OuterTrack_InnerPoint_pz,"cosmicmuon_OuterTrack_InnerPoint_pz[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_OuterPoint_x",cosmicmuon_OuterTrack_OuterPoint_x,"cosmicmuon_OuterTrack_OuterPoint_x[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_OuterPoint_y",cosmicmuon_OuterTrack_OuterPoint_y,"cosmicmuon_OuterTrack_OuterPoint_y[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_OuterPoint_z",cosmicmuon_OuterTrack_OuterPoint_z,"cosmicmuon_OuterTrack_OuterPoint_z[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_OuterPoint_px",cosmicmuon_OuterTrack_OuterPoint_px,"cosmicmuon_OuterTrack_OuterPoint_px[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_OuterPoint_py",cosmicmuon_OuterTrack_OuterPoint_py,"cosmicmuon_OuterTrack_OuterPoint_py[NFilled]/F");
  tree->Branch("cosmicmuon_OuterTrack_OuterPoint_pz",cosmicmuon_OuterTrack_OuterPoint_pz,"cosmicmuon_OuterTrack_OuterPoint_pz[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_InnerPoint_x",cosmicmuon_InnerTrack_InnerPoint_x,"cosmicmuon_InnerTrack_InnerPoint_x[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_InnerPoint_y",cosmicmuon_InnerTrack_InnerPoint_y,"cosmicmuon_InnerTrack_InnerPoint_y[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_InnerPoint_z",cosmicmuon_InnerTrack_InnerPoint_z,"cosmicmuon_InnerTrack_InnerPoint_z[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_InnerPoint_px",cosmicmuon_InnerTrack_InnerPoint_px,"cosmicmuon_InnerTrack_InnerPoint_px[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_InnerPoint_py",cosmicmuon_InnerTrack_InnerPoint_py,"cosmicmuon_InnerTrack_InnerPoint_py[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_InnerPoint_pz",cosmicmuon_InnerTrack_InnerPoint_pz,"cosmicmuon_InnerTrack_InnerPoint_pz[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_OuterPoint_x",cosmicmuon_InnerTrack_OuterPoint_x,"cosmicmuon_InnerTrack_OuterPoint_x[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_OuterPoint_y",cosmicmuon_InnerTrack_OuterPoint_y,"cosmicmuon_InnerTrack_OuterPoint_y[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_OuterPoint_z",cosmicmuon_InnerTrack_OuterPoint_z,"cosmicmuon_InnerTrack_OuterPoint_z[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_OuterPoint_px",cosmicmuon_InnerTrack_OuterPoint_px,"cosmicmuon_InnerTrack_OuterPoint_px[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_OuterPoint_py",cosmicmuon_InnerTrack_OuterPoint_py,"cosmicmuon_InnerTrack_OuterPoint_py[NFilled]/F");
  tree->Branch("cosmicmuon_InnerTrack_OuterPoint_pz",cosmicmuon_InnerTrack_OuterPoint_pz,"cosmicmuon_InnerTrack_OuterPoint_pz[NFilled]/F");
  tree->Branch("cosmicmuon_OuterPoint_x",cosmicmuon_OuterPoint_x,"cosmicmuon_OuterPoint_x[NFilled]/F");
  tree->Branch("cosmicmuon_OuterPoint_y",cosmicmuon_OuterPoint_y,"cosmicmuon_OuterPoint_y[NFilled]/F");
  tree->Branch("cosmicmuon_OuterPoint_z",cosmicmuon_OuterPoint_z,"cosmicmuon_OuterPoint_z[NFilled]/F");
 
  //MC
  tree->Branch("MC_photon_pt",MC_photon_pt,"MC_photon_pt[NFilled]/F");
  tree->Branch("MC_photon_eta",MC_photon_eta,"MC_photon_eta[NFilled]/F");
  tree->Branch("MC_photon_phi",MC_photon_phi,"MC_photon_phi[NFilled]/F");
  tree->Branch("MC_photon_px",MC_photon_pz,"MC_photon_px[NFilled]/F");   
  tree->Branch("MC_photon_py",MC_photon_py,"MC_photon_py[NFilled]/F");
  tree->Branch("MC_photon_pz",MC_photon_pz,"MC_photon_pz[NFilled]/F");
  tree->Branch("MC_photon_E",MC_photon_E,"MC_photon_E[NFilled]/F"); 

  tree->Branch("MC_photon_MotherPt",MC_photon_MotherPt,"MC_photon_MotherPt[NFilled]/F");
  tree->Branch("MC_photon_MotherID",MC_photon_MotherID,"MC_photon_MotherID[NFilled]/I");
  tree->Branch("MC_photon_MotherEta",MC_photon_MotherEta,"MC_photon_MotherEta[NFilled]/F");
  tree->Branch("MC_photon_MotherPhi",MC_photon_MotherPhi,"MC_photon_MotherPhi[NFilled]/F");
  tree->Branch("MC_photon_GrandMotherID",MC_photon_GrandMotherID,"MC_photon_GrandMotherID[NFilled]/I");

  //GRAVITON
  tree->Branch("graviton_pt",&graviton_pt,"graviton_pt/F");

  //MET
  tree->Branch("calo_MetSigma",&calo_MetSigma,"calo_MetSigma/F");
  //tree->Branch("calo_MetCorr",&calo_MetCorr,"calo_MetCorr/F");
  tree->Branch("calo_MetPt",calo_MetPt,"calo_MetPt[6]/F");
  tree->Branch("calo_MetPx",calo_MetPx,"calo_MetPx[6]/F");
  tree->Branch("calo_MetPy",calo_MetPy,"calo_MetPy[6]/F");
  tree->Branch("calo_MetEz",&calo_MetEz,"calo_MetEz/F");
  tree->Branch("calo_MetPhi",calo_MetPhi,"calo_MetPhi[6]/F");
  tree->Branch("calo_MetSumEt",calo_MetSumEt,"calo_MetSumEt[6]/F");
  tree->Branch("calo_EtFractionHadronic",&calo_EtFractionHadronic,"calo_EtFractionHadronic/F");
  tree->Branch("calo_EmEtFraction",&calo_EmEtFraction,"calo_EmEtFraction/F");
  tree->Branch("calo_HadEtInHB",&calo_HadEtInHB,"calo_HadEtInHB/F");
  tree->Branch("calo_HadEtInHE",&calo_HadEtInHE,"calo_HadEtInHE/F");
  tree->Branch("calo_HadEtInHO",&calo_HadEtInHO,"calo_HadEtInHO/F");
  tree->Branch("calo_HadEtInHF",&calo_HadEtInHF,"calo_HadEtInHF/F");
  tree->Branch("calo_EmEtInEB",&calo_EmEtInEB,"calo_EmEtInEB/F");
  tree->Branch("calo_EmEtInEE",&calo_EmEtInEE,"calo_EmEtInEE/F");
  tree->Branch("calo_EmEtInHF",&calo_EmEtInHF,"calo_EmEtInHF/F");
  tree->Branch("calo_MaxEtInEmTowers",&calo_MaxEtInEmTowers,"calo_MaxEtInEmTowers/F");
  tree->Branch("calo_MaxEtInHadTowers",&calo_MaxEtInHadTowers,"calo_MaxEtInHadTowers/F");
  tree->Branch("gen_MetPt",&gen_MetPt,"gen_MetPt/F");
  tree->Branch("gen_MetPx",&gen_MetPx,"gen_MetPx/F");
  tree->Branch("gen_MetPy",&gen_MetPy,"gen_MetPy/F");
  tree->Branch("gen_MetPhi",&gen_MetPhi,"gen_MetPhi/F");
  tree->Branch("gen_MetSumEt",&gen_MetSumEt,"gen_MetSumEt/F");
  tree->Branch("delta_phi",&delta_phi,"delta_phi/F");
  tree->Branch("delta_phiGEN",&delta_phiGEN,"delta_phiGEN/F");
  tree->Branch("pf_MetPt",pf_MetPt,"pf_MetPt[6]/F");
  tree->Branch("pf_MetPx",pf_MetPx,"pf_MetPx[6]/F");
  tree->Branch("pf_MetPy",pf_MetPy,"pf_MetPy[6]/F");
  tree->Branch("pf_MetPhi",pf_MetPhi,"pf_MetPhi[6]/F");
  tree->Branch("pf_MetSumEt",pf_MetSumEt,"pf_MetSumEt[6]/F");
  tree->Branch("delta_phiPF",&delta_phiPF,"delta_phiPF/F");
  tree->Branch("tc_MetPt",tc_MetPt,"tc_MetPt[6]/F");
  tree->Branch("tc_MetPx",tc_MetPx,"tc_MetPx[6]/F");
  tree->Branch("tc_MetPy",tc_MetPy,"tc_MetPy[6]/F");
  tree->Branch("tc_MetPhi",tc_MetPhi,"tc_MetPhi[6]/F");
  tree->Branch("tc_MetSumEt",tc_MetSumEt,"tc_MetSumEt[6]/F");
  tree->Branch("delta_phiTC",&delta_phiTC,"delta_phiTC/F");
  
  //EVENT TYPE
  tree->Branch("isZ_event",&isZ_event,"isZ_event/O");
  tree->Branch("isW_event",&isW_event,"isW_event/O");
  tree->Branch("isZnunu_event",&isZnunu_event,"isZnunu_event/O");
  tree->Branch("isZelec_event",&isZelec_event,"isZelec_event/O");
  tree->Branch("isZmu_event",&isZmu_event,"isZmu_event/O");
  tree->Branch("isZtau_event",&isZtau_event,"isZtau_event/O");
  tree->Branch("isWelec_event",&isWelec_event,"isWelec_event/O");
  tree->Branch("isWmu_event",&isWmu_event,"isWmu_event/O");
  tree->Branch("isWtau_event",&isWtau_event,"isWtau_event/O");
  tree->Branch("isSingleHardPhoton_event",&isSingleHardPhoton_event,"isSingleHardPhoton_event/O");
  tree->Branch("isdiphoton_event",&isdiphoton_event,"isdiphoton_event/O");
  tree->Branch("isisr_photon_event",&isisr_photon_event,"isisr_photon_event/O");

  //csc
  tree->Branch("cscseg_n",&cscseg_n,"cscseg_n/I");
  tree->Branch("cscseg_x",cscseg_x,"cscseg_x[cscseg_n]/F");
  tree->Branch("cscseg_y",cscseg_y,"cscseg_y[cscseg_n]/F");
  tree->Branch("cscseg_z",cscseg_z,"cscseg_z[cscseg_n]/F");
  tree->Branch("cscseg_time",cscseg_time,"cscseg_time[cscseg_n]/F");

  //BeamHalo Sumary
   tree->Branch("is_BeamHaloGlobalLoosePass",&is_BeamHaloGlobalLoosePass,"is_BeamHaloGlobalLoosePass/O");
   tree->Branch("is_BeamHaloGlobalTightPass",&is_BeamHaloGlobalTightPass,"is_BeamHaloGloablTightPass/O");
   tree->Branch("is_BeamHaloHcalLoosePass",&is_BeamHaloHcalLoosePass,"is_BeamHaloHcalLoosePass/O");
   tree->Branch("is_BeamHaloHcalTightPass",&is_BeamHaloHcalTightPass,"is_BeamHaloHcalTightPass/O");
   tree->Branch("is_BeamHaloCSCLoosePass",&is_BeamHaloCSCLoosePass,"is_BeamHaloCSCLoosePass/O");
   tree->Branch("is_BeamHaloCSCTightPass",&is_BeamHaloCSCTightPass,"is_BeamHaloCSCTightPass/O");
   tree->Branch("is_BeamHaloEcalLoosePass",&is_BeamHaloEcalLoosePass,"is_BeamHaloEcalLoosePass/O");
   tree->Branch("is_BeamHaloEcalTightPass",&is_BeamHaloEcalTightPass,"is_BeamHaloEcalTightPass/O");
   tree->Branch("is_BeamHaloIDTightPass",&is_BeamHaloIDTightPass,"is_BeamHaloIDTightPass/O");
   tree->Branch("is_BeamHaloIDLoosePass",&is_BeamHaloIDLoosePass,"is_BeamHaloIDLoosePass/O");
   tree->Branch("is_SmellsLikeHalo_Tag",&is_SmellsLikeHalo_Tag, "is_SmellsLikeHalo_Tag/O");
   tree->Branch("is_LooseHalo_Tag",&is_LooseHalo_Tag, "is_LooseHalo_Tag/O");
   tree->Branch("is_TightHalo_Tag",&is_TightHalo_Tag, "is_TightHalo_Tag/O");
   tree->Branch("is_ExtremeTightHalo_Tag",&is_ExtremeTightHalo_Tag, "is_ExtremeTightHalo_Tag/O");


  ///information different separate from the original ntuple
  tree->Branch("jet_pt",jet_pt,"jet_pt[NJFilled]/F"); 
  tree->Branch("jet_px",jet_px,"jet_px[NJFilled]/F"); 
  tree->Branch("jet_py",jet_py,"jet_py[NJFilled]/F"); 
  tree->Branch("jet_pz",jet_pz,"jet_pz[NJFilled]/F"); 
  tree->Branch("jet_E",jet_E,"jet_E[NJFilled]/F");
  tree->Branch("jet_eta",jet_eta,"jet_eta[NJFilled]/F");
  tree->Branch("jet_phi",jet_phi,"jet_phi[NJFilled]/F");
  tree->Branch("jet_EMenergyFraction",jet_EMenergyFraction,"jet_EMenergyFraction[NJFilled]/F");
  tree->Branch("jet_energyFractionHadronic",jet_energyFractionHadronic,"jet_energyFractionHadronic[NJFilled]/F");
  tree->Branch("jet_hitsInN90",jet_hitsInN90,"jet_hitsInN90[NJFilled]/I");
  tree->Branch("jet_n90Hits",jet_n90Hits,"jet_n90Hits[NJFilled]/I");
  tree->Branch("jet_nTowers",jet_nTowers,"jet_nTowers[NJFilled]/I");
  tree->Branch("jet_fHPD",jet_fHPD,"jet_fHPD[NJFilled]/F");
  tree->Branch("jet_fRBX",jet_fRBX,"jet_fRBX[NJFilled]/F");
  tree->Branch("jet_RHF",jet_RHF,"jet_RHF[NJFilled]/F");
  tree->Branch("jet_n_PhotonRemoved",&jet_n_PhotonRemoved,"jet_n_PhotonRemoved/I");
  tree->Branch("jet_n_FromHT",&jet_n_FromHT,"jet_n_FromHT/I");
  tree->Branch("HT",&HT,"HT/F");
  tree->Branch("HT_JetpT50",&HT_JetpT50,"HT_JetpT50/F");  
  tree->Branch("HT_JetpT100",&HT_JetpT100,"HT_JetpT100/F");
  tree->Branch("HT_JetpT150",&HT_JetpT150,"HT_JetpT150/F");


  //PF jet informatoin  
  tree->Branch("pfjet_pt",pfjet_pt,"pfjet_pt[NJFilled]/F"); 
  tree->Branch("pfjet_px",pfjet_px,"pfjet_px[NJFilled]/F"); 
  tree->Branch("pfjet_py",pfjet_py,"pfjet_py[NJFilled]/F"); 
  tree->Branch("pfjet_pz",pfjet_pz,"pfjet_pz[NJFilled]/F"); 
  tree->Branch("pfjet_E",pfjet_E,"pfjet_E[NJFilled]/F");
  tree->Branch("pfjet_eta",pfjet_eta,"pfjet_eta[NJFilled]/F");
  tree->Branch("pfjet_phi",pfjet_phi,"pfjet_phi[NJFilled]/F");
<<<<<<< myCuts_halo.C
  tree->Branch("pfjet_n_PhotonRemoved",&pfjet_n_PhotonRemoved,"pfjet_n_PhotonRemoved/I");
  tree->Branch("pfjet_n_FromHT",&pfjet_n_FromHT,"pfjet_n_FromHT/I");
  tree->Branch("pfHT",&pfHT,"pfHT/F");
  tree->Branch("HT_pfJetpT50",&HT_pfJetpT50,"HT_pfJetpT50/F");  
  tree->Branch("HT_pfJetpT100",&HT_pfJetpT100,"HT_pfJetpT100/F");
  tree->Branch("HT_pfJetpT150",&HT_pfJetpT150,"HT_pfJetpT150/F");
                 


  ///////////////////////added by bhawna/////////////////////
  //////////////////////////////////////////////////////////
  //tau info
  tree->Branch("tau_n",&tau_n,"tau_n/I");
  tree->Branch("OneProng0Pi0",OneProng0Pi0,"OneProng0Pi0[NFilled]/I");
  tree->Branch("OneProng1Pi0",OneProng1Pi0,"OneProng1Pi0[NFilled]/I");
  tree->Branch("OneProng2Pi0",OneProng2Pi0,"OneProng2Pi0[NFilled]/I");
  tree->Branch("ThreeProng0Pi0",ThreeProng0Pi0,"ThreeProng0Pi0[NFilled]/I");
  tree->Branch("ThreeProng1Pi0",ThreeProng1Pi0,"ThreeProng1Pi0[NFilled]/I");
  tree->Branch("GenHadTauPt",GenHadTauPt,"GenHadTauPt[NFilled]/F");
  tree->Branch("GenHadTauEta",GenHadTauEta,"GenHadTauEta[NFilled]/F");
  tree->Branch("GenHadTauPhi",GenHadTauPhi,"GenHadTauPhi[NFilled]/F");
  tree->Branch("NPions",&NPions,"NPions/I");
  tree->Branch("NPi0",&NPi0,"NPi0/I");
  tree->Branch("NPhotons",&NPhotons,"NPhotons/I");
  tree->Branch("pionPdgId",pionPdgId,"pionPdgId[NFilled][NFilled]/I");
  tree->Branch("pionPt",pionPt,"pionPt[NFilled][NFilled]/F");
  tree->Branch("pionEta",pionEta,"pionEta[NFilled][NFilled]/F");
  tree->Branch("pionPhi",pionPhi,"pionPhi[NFilled][NFilled]/F");
  
  tree->Branch("pi0PdgId",pi0PdgId,"pi0PdgId[NFilled][NFilled]/I");
  tree->Branch("pi0Pt",pi0Pt,"pi0Pt[NFilled][NFilled]/F");
  tree->Branch("pi0Eta",pi0Eta,"pi0Eta[NFilled][NFilled]/F");
  tree->Branch("pi0Phi",pi0Phi,"pi0Phi[NFilled][NFilled]/F");


  tree->Branch("photonPdgId",photonPdgId,"photonPdgId[NFilled][NFilled]/I");
  tree->Branch("photonPt",photonPt,"photonPt[NFilled][NFilled]/F");
  tree->Branch("photonEta",photonEta,"photonEta[NFilled][NFilled]/F");
  tree->Branch("photonPhi",photonPhi,"photonPhi[NFilled][NFilled]/F");  
=======
  tree->Branch("pfjet_n_PhotonRemoved",&pfjet_n_PhotonRemoved,"pfjet_n_PhotonRemoved/I");
  tree->Branch("pfjet_n_FromHT",&pfjet_n_FromHT,"pfjet_n_FromHT/I");
  tree->Branch("pfHT",&pfHT,"pfHT/F");
  tree->Branch("HT_pfJetpT50",&HT_pfJetpT50,"HT_pfJetpT50/F");  
  tree->Branch("HT_pfJetpT100",&HT_pfJetpT100,"HT_pfJetpT100/F");
  tree->Branch("HT_pfJetpT150",&HT_pfJetpT150,"HT_pfJetpT150/F");
                 


  ///////////////////////added by bhawna/////////////////////
  //////////////////////////////////////////////////////////
  //tau info
  tree->Branch("tau_n",&tau_n,"tau_n/I");
  tree->Branch("OneProng0Pi0",OneProng0Pi0,"OneProng0Pi0[NFilled]/I");
  tree->Branch("OneProng1Pi0",OneProng1Pi0,"OneProng1Pi0[NFilled]/I");
  tree->Branch("OneProng2Pi0",OneProng2Pi0,"OneProng2Pi0[NFilled]/I");
  tree->Branch("ThreeProng0Pi0",ThreeProng0Pi0,"ThreeProng0Pi0[NFilled]/I");
  tree->Branch("ThreeProng1Pi0",ThreeProng1Pi0,"ThreeProng1Pi0[NFilled]/I");
  tree->Branch("GenHadTauPt",GenHadTauPt,"GenHadTauPt[NFilled]/F");
  tree->Branch("GenHadTauEta",GenHadTauEta,"GenHadTauEta[NFilled]/F");
  tree->Branch("GenHadTauPhi",GenHadTauPhi,"GenHadTauPhi[NFilled]/F");
  tree->Branch("NPions",&NPions,"NPions/I");
  tree->Branch("NPi0",&NPi0,"NPi0/I");
  tree->Branch("NPhotons",&NPhotons,"NPhotons/I");
  tree->Branch("pionPdgId",pionPdgId,"pionPdgId[NFilled][NFilled]/I");
  tree->Branch("pionPt",pionPt,"pionPt[NFilled][NFilled]/F");
  tree->Branch("pionEta",pionEta,"pionEta[NFilled][NFilled]/F");
  tree->Branch("pionPhi",pionPhi,"pionPhi[NFilled][NFilled]/F");
  
  tree->Branch("pi0PdgId",pi0PdgId,"pi0PdgId[NFilled][NFilled]/I");
  tree->Branch("pi0Pt",pi0Pt,"pi0Pt[NFilled][NFilled]/F");
  tree->Branch("pi0Eta",pi0Eta,"pi0Eta[NFilled][NFilled]/F");
  tree->Branch("pi0Phi",pi0Phi,"pi0Phi[NFilled][NFilled]/F");
>>>>>>> 1.2

  /////////////ends here///////////////////////////
  
  
  
  
  

  tree->Branch("photonPdgId",photonPdgId,"photonPdgId[NFilled][NFilled]/I");
  tree->Branch("photonPt",photonPt,"photonPt[NFilled][NFilled]/F");
  tree->Branch("photonEta",photonEta,"photonEta[NFilled][NFilled]/F");
  tree->Branch("photonPhi",photonPhi,"photonPhi[NFilled][NFilled]/F");  

  /////////////ends here///////////////////////////
  
  
  
  
  

  //scrapping Halo etc
  tree->Branch("scraping_isScrapingEvent",&scraping_isScrapingEvent,"scraping_isScrapingEvent/O");
  tree->Branch("IsHEHalo",IsHEHalo,"IsHEHalo[NFilled]/O");
  tree->Branch("IsTrackHalo",IsTrackHalo,"IsTrackHalo[NFilled]/O");
  tree->Branch("IsE2E9Spike",IsE2E9Spike,"IsE2E9Spike[NFilled]/O");
  tree->Branch("IsTimeSpike",IsTimeSpike,"IsTimeSpike[NFilled]/O");
  tree->Branch("IsCSCHalo",IsCSCHalo,"IsCSCHalo[NFilled]/O");
  tree->Branch("DeltaPhiCSCHalo",DeltaPhiCSCHalo,"DeltaPhiCSCHalo[NFilled][100]/F");
  //For comsmic rejection
  tree->Branch("sigmaRC",sigmaRC,"sigmaRC[NFilled]/F");
  tree->Branch("fisherd2",fisherd2,"fisherd2[NFilled]/F");
  tree->Branch("fisherd3",fisherd3,"fisherd3[NFilled]/F");

	
}
//======================================================================

void myCuts::Loop(){
  cout<<"starting event loop:"<<endl;
  bool debug_;
    debug_=false;
  Long64_t jentry, nentries,nb,nbytes;
  nentries = Int_t(fchain->GetEntries());
  cout << " nentries " << nentries << endl;
  nbytes = 0, nb = 0; NonCollisionBG *bg=new NonCollisionBG();


   //add branches in cache for faster run
    Int_t cachesize = 100000000; //100 MBytes
    fChain->SetCacheSize(cachesize); //<<<
    fChain->AddBranchToCache("*");


//--------Loop over events
  for(jentry =0; jentry < nentries; jentry++){
    Int_t ientry = LoadTree(jentry) ; 
    if(ientry == 0)cout<<"file name:"<<fchain->GetFile()->GetName()<<endl;
    //initializing variables here:
   if(debug_)cout<<"starting Initilalization of variable"<<endl;


    NFilled=3;
    NJFilled=5;

    for(int i = 0; i <NFilled; i++){
      vertex_x[i]                           =-99.;
      vertex_y[i]                           =-99.;
      vertex_z[i]                           =-99.;
      vertex_tracksize[i]                   =-99.;
      vertex_ndof[i]                        =-99.;
      vertex_chi2[i]                        =-99.;
      vertex_d0[i]                          =-99.;
      vertex_isFake[i]                      =false;
      
      //Photon
      photon_e2e9[i]                        =-99.; 
      photon_e6e2[i]                        =-99.; 
      photon_e4e1[i]                        =-99.; 
      photon_E[i]                           =-99.;
      photon_pt[i]                          =-99.;
      //photon_iLICTD[i]                      =-99.;
      photon_eta[i]                         =-99.;
      photon_phi[i]                         =-99.;
      photon_theta[i]                       =-99.;
      photon_et[i]                          =-99.;
      photon_swissCross[i]                  =-99.;
      photon_r9[i]                          =-99.;
      photon_e1x5[i]                        =-99.;
      photon_e2x5[i]                        =-99.;
      photon_e3x3[i]                        =-99.;
      photon_e5x5[i]                        =-99.;
      photon_r1x5[i]                        =-99.;
      photon_r2x5[i]                        =-99.;
      photon_maxEnergyXtal[i]               =-99.;
      photon_SigmaEtaEta[i]                 =-99.;
      photon_SigmaIetaIeta[i]               =-99.;
      photon_SigmaPhiPhi[i]                 =-99.;
      photon_SigmaIphiIphi[i]               =-99.;
      photon_Roundness[i]                   =-99.;
      photon_Angle[i]                       =-99.;
      photon_ecalRecHitSumEtConeDR03[i]     =-99.;
      photon_hcalTowerSumEtConeDR03[i]      =-99.;
      photon_trkSumPtSolidConeDR03[i]       =-99.;
      photon_trkSumPtHollowConeDR03[i]      =-99.;
      photon_nTrkSolidConeDR03[i]           =-99;
      photon_nTrkHollowConeDR03[i]          =-99;
      photon_hcalDepth1TowerSumEtConeDR03[i]=-99.;
      photon_hcalDepth2TowerSumEtConeDR03[i]=-99.;
      photon_ecalRecHitSumEtConeDR04[i]     =-99.;
      photon_hcalTowerSumEtConeDR04[i]      =-99.;
      photon_trkSumPtSolidConeDR04[i]       =-99.;
      photon_trkSumPtHollowConeDR04[i]      =-99.;
      photon_nTrkSolidConeDR04[i]           =-99;
      photon_nTrkHollowConeDR04[i]          =-99;
      photon_hcalDepth1TowerSumEtConeDR04[i]=-99.;
      photon_hcalDepth2TowerSumEtConeDR04[i]=-99.;
      photon_HoE[i]                         =-99.;
      photon_px[i]                          =-99.;
      photon_py[i]                          =-99.;
      photon_pz[i]                          =-99.;
      photon_basiccluster_size[i]           =-99;
      photon_sc_energy[i]                   =-99.;
      photon_sc_eta[i]                      =-99.;
      photon_sc_phi[i]                      =-99.;
      photon_sc_x[i]                        =-99.;
      photon_sc_y[i]                        =-99.;
      photon_sc_z[i]                        =-99.;
      photon_etaWidth[i]                    =-99.;
      photon_phiWidth[i]                    =-99.;
      photon_sc_et[i]                       =-99.;
      photon_ntracks[i]                     =-99;
      photon_pairInvmass[i]                 =-99.;
      photon_pairCotThetaSeparation[i]      =-99.;
      photon_pairmomentumX[i]               =-99.;
      photon_pairmomentumY[i]               =-99.;
      photon_pairmomentumZ[i]               =-99.;
      photon_EoverP[i]                      =-99.;
      photon_ConvVx[i]                      =-99.;
      photon_ConvVy[i]                      =-99.;
      photon_ConvVz[i]                      =-99.;
      photon_ZOfPrimaryVertex[i]            =-99.;
      photon_distOfMinimumApproach[i]       =-99.;
      photon_dPhiTracksAtVtx[i]             =-99.;
      photon_dPhiTracksAtEcal[i]            =-99.;
      photon_dEtaTracksAtEcal[i]            =-99.;
      photon_isconverted[i]  = false;
      photon_hasPixelSeed[i] = false;
      photon_isEB[i]         = false;
      photon_isEE[i]         = false;
      photon_isEBGap[i]      = false;
      photon_isEEGap[i]      = false;
      photon_isEBEEGap[i]    = false; 
      for(int j=0;j<2;j++)
	{	
	  photon_timing_xtal[i][j]         = -99.;
	  //photon_energy_xtal[i][j]       = -99.;
  	}
      
      //Tracks
      track_px[i]                          =-99.;   
      track_py[i]                          =-99.;   
      track_pz[i]                          =-99.;   
      track_pt[i]                          =-99.;   
      track_eta[i]                         =-99.;   
      track_phi[i]                         =-99.;   
      
      //MC photon
      MC_photon_pt[i]                      =-99.;
      MC_photon_eta[i]                     =-99.;
      MC_photon_phi[i]                     =-99.;

      MC_photon_px[i]                      =-99.;
      MC_photon_py[i]                      =-99.;
      MC_photon_pz[i]                      =-99.;
      MC_photon_E[i]                       =-99.;

      MC_photon_MotherID[i]                =-99;
      MC_photon_MotherPt[i]                =-99.;
      MC_photon_MotherEta[i]               =-99.;
      MC_photon_MotherPhi[i]               =-99.;
      MC_photon_GrandMotherID[i]           =-99;

      ////////////////////added by bhawna///////////////////
      //tau info
      
      OneProng0Pi0[i]                    =-99;
      OneProng1Pi0[i]                    =-99;
      OneProng2Pi0[i]                    =-99;
      ThreeProng0Pi0[i]                  =-99;
      ThreeProng1Pi0[i]                  =-99;
      GenHadTauPt[i]                     =-99;
      GenHadTauEta[i]                    =-99;
      GenHadTauPhi[i]                    =-99;
      NPions[i]                          =-99;
      NPi0[i]                            =-99;
      NPi0[i]                            =-99;


      for(int j = 0; j <2; j++){
      pionPdgId[i][j]                     =-99;
      pionPt[i][j]                        =-99;
      pionEta[i][j]                       =-99;
      pionPhi[i][j]                       =-99;

      pi0PdgId[i][j]                      =-99;
      pi0Pt[i][j]                         =-99;
      pi0Eta[i][j]                        =-99;
      pi0Phi[i][j]                        =-99;

      photonPt[i][j]                       =-99; 
      photonEta[i][j]                      =-99;
      photonPhi[i][j]                      =-99;
      photonPdgId[i][j]                    =-99;
      }
      /////////////////////////ends here///////////////////////

      //gsfElectron
      elec_px[i]                          =-99.;
      elec_py[i]                          =-99.;
      elec_pz[i]                          =-99.;
      elec_pt[i]                          =-99.;
      elec_eta[i]                         =-99.;
      elec_phi[i]                         =-99.;
      elec_energy[i]                      =-99.;
      elec_charge[i]                      =-99.;
      elec_trkIso[i]                      =-99.;    
      elec_ecalIso[i]                     =-99.;    
      elec_hcalIso[i]                     =-99.;    
      elec_HoE[i]                         =-99.;    
      elec_SigmaIetaIeta[i]               =-99.;
      elec_dEtaIn[i]                      =-99.;
      elec_dPhiIn[i]                      =-99.;
      elec_sc_energy[i]                   =-99.;
      elec_sc_eta[i]                      =-99.;
      elec_sc_phi[i]                      =-99.;
      
      //Muon
      muon_px[i]                           =-99.;
      muon_py[i]                           =-99.;
      muon_pz[i]                           =-99.;
      muon_pt[i]                           =-99.;
      muon_eta[i]                          =-99.;
      muon_phi[i]                          =-99.;
      muon_energy[i]                       =-99.;
      muon_charge[i]                       =-99.;
      muon_OuterTrack_InnerPoint_x[i]      =-99.;
      muon_OuterTrack_InnerPoint_y[i]      =-99.;
      muon_OuterTrack_InnerPoint_z[i]      =-99.;
      muon_OuterTrack_InnerPoint_px[i]     =-99.;
      muon_OuterTrack_InnerPoint_py[i]     =-99.;
      muon_OuterTrack_InnerPoint_pz[i]     =-99.;
      muon_OuterTrack_OuterPoint_x[i]      =-99.;
      muon_OuterTrack_OuterPoint_y[i]      =-99.;
      muon_OuterTrack_OuterPoint_z[i]      =-99.;
      muon_OuterTrack_OuterPoint_px[i]     =-99.;
      muon_OuterTrack_OuterPoint_py[i]     =-99.;
      muon_OuterTrack_OuterPoint_pz[i]     =-99.;
      muon_InnerTrack_InnerPoint_x[i]      =-99.;
      muon_InnerTrack_InnerPoint_y[i]      =-99.;
      muon_InnerTrack_InnerPoint_z[i]      =-99.;
      muon_InnerTrack_InnerPoint_px[i]     =-99.;
      muon_InnerTrack_InnerPoint_py[i]     =-99.;
      muon_InnerTrack_InnerPoint_pz[i]     =-99.;
      muon_InnerTrack_OuterPoint_x[i]      =-99.;
      muon_InnerTrack_OuterPoint_y[i]      =-99.;
      muon_InnerTrack_OuterPoint_z[i]      =-99.;
      muon_InnerTrack_OuterPoint_px[i]     =-99.;
      muon_InnerTrack_OuterPoint_py[i]     =-99.;
      muon_InnerTrack_OuterPoint_pz[i]     =-99.;  
      
      muon_isGlobalMuon[i]                 = false;
      muon_isTrackerMuon[i]                = false;
      muon_isStandAloneMuon[i]             = false;
      muon_InnerTrack_isNonnull[i]         = false;
      muon_OuterTrack_isNonnull[i]         = false;

      muon_OuterPoint_x[i]                   = -99.;   
      muon_OuterPoint_y[i]                   = -99.;
      muon_OuterPoint_z[i]                   = -99.;
      muon_InnerPoint_x[i]                   = -99.;   
      muon_InnerPoint_y[i]                   = -99.;
      muon_InnerPoint_z[i]                   = -99.;
      
      muon_vx[i]                            = -99.;   
      muon_vy[i]                            = -99.;
      muon_vz[i]                            = -99.;

     //cosmicmuon
      cosmicmuon_px[i]                           =-99.;
      cosmicmuon_py[i]                           =-99.;
      cosmicmuon_pz[i]                           =-99.;
      cosmicmuon_pt[i]                           =-99.;
      cosmicmuon_eta[i]                          =-99.;
      cosmicmuon_phi[i]                          =-99.;
      cosmicmuon_energy[i]                       =-99.;
      cosmicmuon_charge[i]                       =-99.;
      cosmicmuon_OuterTrack_InnerPoint_x[i]      =-99.;
      cosmicmuon_OuterTrack_InnerPoint_y[i]      =-99.;
      cosmicmuon_OuterTrack_InnerPoint_z[i]      =-99.;
      cosmicmuon_OuterTrack_InnerPoint_px[i]     =-99.;
      cosmicmuon_OuterTrack_InnerPoint_py[i]     =-99.;
      cosmicmuon_OuterTrack_InnerPoint_pz[i]     =-99.;
      cosmicmuon_OuterTrack_OuterPoint_x[i]      =-99.;
      cosmicmuon_OuterTrack_OuterPoint_y[i]      =-99.;
      cosmicmuon_OuterTrack_OuterPoint_z[i]      =-99.;
      cosmicmuon_OuterTrack_OuterPoint_px[i]     =-99.;
      cosmicmuon_OuterTrack_OuterPoint_py[i]     =-99.;
      cosmicmuon_OuterTrack_OuterPoint_pz[i]     =-99.;
      cosmicmuon_InnerTrack_InnerPoint_x[i]      =-99.;
      cosmicmuon_InnerTrack_InnerPoint_y[i]      =-99.;
      cosmicmuon_InnerTrack_InnerPoint_z[i]      =-99.;
      cosmicmuon_InnerTrack_InnerPoint_px[i]     =-99.;
      cosmicmuon_InnerTrack_InnerPoint_py[i]     =-99.;
      cosmicmuon_InnerTrack_InnerPoint_pz[i]     =-99.;
      cosmicmuon_InnerTrack_OuterPoint_x[i]      =-99.;
      cosmicmuon_InnerTrack_OuterPoint_y[i]      =-99.;
      cosmicmuon_InnerTrack_OuterPoint_z[i]      =-99.;
      cosmicmuon_InnerTrack_OuterPoint_px[i]     =-99.;
      cosmicmuon_InnerTrack_OuterPoint_py[i]     =-99.;
      cosmicmuon_InnerTrack_OuterPoint_pz[i]     =-99.;  
      
      cosmicmuon_isGlobalcosmicmuon[i]        = false;
      cosmicmuon_isTrackercosmicmuon[i]       = false;
      cosmicmuon_isStandAlonecosmicmuon[i]    = false;
      cosmicmuon_InnerTrack_isNonnull[i]      = false;
      cosmicmuon_OuterTrack_isNonnull[i]      = false;

       cosmicmuon_OuterPoint_x[i]                   = -99.;   
       cosmicmuon_OuterPoint_y[i]                   = -99.;
       cosmicmuon_OuterPoint_z[i]                   = -99.;


    }//loop over 100 entries
	
    for (int i =0; i< 6;i++)
      {		
	calo_MetPt[i]                      =-99.;
	calo_MetPx[i]                      =-99.;
	calo_MetPy[i]                      =-99.;
	calo_MetPhi[i]                     =-99.;
	calo_MetSumEt[i]                   =-99.;
	pf_MetPt[i]                        =-99.;
	pf_MetPx[i]                        =-99.;
	pf_MetPy[i]                        =-99.;
	pf_MetPhi[i]                       =-99.;
	pf_MetSumEt[i]                     =-99.;
	tc_MetPt[i]                        =-99.;
	tc_MetPx[i]                        =-99.;
	tc_MetPy[i]                        =-99.;
	tc_MetPhi[i]                       =-99.;
	tc_MetSumEt[i]                     =-99.;
      }

    calo_MetSigma                          =-99.;
    calo_MetEz                             =-99.;
    calo_EtFractionHadronic                =-99.;
    calo_EmEtFraction                      =-99.;
    calo_HadEtInHB                         =-99.;
    calo_HadEtInHE                         =-99.;
    calo_HadEtInHO                         =-99.;
    calo_HadEtInHF                         =-99.;
    calo_EmEtInEB                          =-99.;
    calo_EmEtInEE                          =-99.;
    calo_EmEtInHF                          =-99.;
    calo_MaxEtInEmTowers                   =-99.;
    calo_MaxEtInHadTowers                  =-99.;
    delta_phi                              =-99.;
    delta_phiGEN                           =-99.;
    delta_phiPF                            =-99.;
    delta_phiTC                            =-99.;
    gen_MetPt                              =-99.;
    gen_MetPx                              =-99.;
    gen_MetPy                              =-99.;
    gen_MetPhi                             =-99.;
    gen_MetSumEt                           =-99.;
	
    //Size of objects
    photon_n                              =-99;
    vertex_n                              =-99;
    jet_n                                 =-99;
    pfjet_n                               =-99;
    track_n                               =-99;
    muon_n                                =-99;
    elec_n                                =-99;
    cscseg_n                              =-99;
    gen_gravitonpt                        =-99.;
    tau_n                              =-99;

    //BeamHaloo Summary 
    
  is_BeamHaloIDTightPass     = false;
  is_BeamHaloIDLoosePass     = false; 
         
  is_BeamHaloEcalLoosePass  = false;
  is_BeamHaloHcalLoosePass  = false;
  is_BeamHaloCSCLoosePass   = false;
  is_BeamHaloGlobalLoosePass =false;
         
  is_BeamHaloEcalTightPass   =false;
  is_BeamHaloHcalTightPass   =false;
  is_BeamHaloCSCTightPass    =false;
  is_BeamHaloGlobalTightPass =false;

  is_SmellsLikeHalo_Tag      =false;
  is_LooseHalo_Tag           =false; 
  is_TightHalo_Tag           =false;
  is_ExtremeTightHalo_Tag    =false;


    //Event Type
    isZ_event                = false;
    isW_event                = false;
    isZnunu_event            = false;
    isZelec_event            = false;
    isZmu_event              = false;
    isZtau_event             = false;
    isWelec_event            = false;
    isWmu_event              = false;
    isWtau_event             = false;
    isSingleHardPhoton_event = false;
    isdiphoton_event         = false;
    isisr_photon_event       = false;
    scraping_isScrapingEvent = false;

    //csc
    for(int j=0;j<100;j++){
      cscseg_x[j]                   =-99.;
      cscseg_y[j]                   =-99.;
      cscseg_z[j]                   =-99.;
      cscseg_time[j]                =-99.;
    }
    
    //HLT
    ntriggers = 0;
    HLTNames.clear();	
    HLTPrescales.clear();	
    HLTIfPassed.clear();	
     
    //information different separate from the original ntuple
    for(int i = 0; i<NJFilled;i++){
      jet_pt[i]                         =-99.;
      jet_px[i]                         =-99.;
      jet_py[i]                         =-99.;
      jet_pz[i]                         =-99.;
      jet_E[i]                          =-99.;
      jet_eta[i]                        =-99.;
      jet_phi[i]                        =-99.;
      jet_EMenergyFraction[i]           =-99.;
      jet_energyFractionHadronic[i]     =-99.;
      jet_hitsInN90[i]                  =-99;
      jet_n90Hits[i]                    =-99;
      jet_nTowers[i]                    =-99;
      jet_fHPD[i]                       =-99.;
      jet_fRBX[i]                       =-99.;
      jet_RHF[i]                        =-99.;


      pfjet_pt[i]                         =-99.;
      pfjet_px[i]                         =-99.;
      pfjet_py[i]                         =-99.;
      pfjet_pz[i]                         =-99.;
      pfjet_E[i]                          =-99.;
      pfjet_eta[i]                        =-99.;
      pfjet_phi[i]                        =-99.;

      IsHEHalo[i]                 = false;
      IsCSCHalo[i]                = false;     
      IsTrackHalo[i]              = false;
      IsE2E9Spike[i]              = false;
      IsTimeSpike[i]              = false;

      for(int j=0;j<100;j++) DeltaPhiCSCHalo[i][j]=-99.;
 
      
      //For Cosmic
      fisherd2[i]                       = -99.;
      fisherd3[i]                       = -99.;
      sigmaRC[i]                        = -99.;
    }

     jet_n_PhotonRemoved = 0;
     jet_n_FromHT        = 0;

     HT          =0.0;
     HT_JetpT50  =0.0;
     HT_JetpT100 =0.0;
     HT_JetpT150 =0.0;

     pfjet_n_PhotonRemoved = 0;
     pfjet_n_FromHT        = 0;

     pfHT          =0.0;
     HT_pfJetpT50  =0.0;
     HT_pfJetpT100 =0.0;
     HT_pfJetpT150 =0.0;
  
    

    /////////////////////////////////////////////////////////////////////
    //*********************Initialization complete*********************//
    /////////////////////////////////////////////////////////////////////
    
    
    /////////////////////////////////////////////////////////////////////
    //***************NOW FILL SOME EVENT Variables*********************//
     ////////////////////////////////////////////////////////////////////

   if(debug_)cout<<"starting assigning values to new brancehs"<<endl;

  //Get this event
   nb = fChain->GetEntry(jentry);

   if(debug_)cout<<"nb ="<<nb<<endl;
    Nnevents             = nevents;
    Nrun                 = run;
    Nevent               = event;
    NluminosityBlock     = luminosityBlock ;
    NbeamCrossing        = beamCrossing ;
    NtotalIntensityBeam1 = totalIntensityBeam1;
    NtotalIntensityBeam2 = totalIntensityBeam2;
    NavgInsDelLumi       = avgInsDelLumi;
    NavgInsDelLumiErr    = avgInsDelLumiErr;
    NavgInsRecLumi       = avgInsRecLumi;
    NavgInsRecLumiErr    = avgInsRecLumiErr;
    
    /////////////////////////////////////////////////////////////////////
    //***************How  many Entries to Fill*************************//
    /////////////////////////////////////////////////////////////////////
    for(int i=0;i<ntriggers;i++){	
      HLTNames.push_back((*triggernames)[i]);
      HLTPrescales.push_back((*triggerprescales)[i]);
      HLTIfPassed.push_back((*ifTriggerpassed)[i]);
    }
    //////////////////////////////////////////////////////
    /////to be turned on for MC only SSC /////////////////
    //////////////////////////////////////////////////////


    isZ_event                = is_Z_event;
    isW_event                = is_W_event;
    isZnunu_event            = is_Znunu_event;
    isZelec_event            = is_Zelec_event;
    isZmu_event              = is_Zmu_event;
    isZtau_event             = is_Ztau_event;
    isWelec_event            = is_Welec_event;
    isWmu_event              = isWmu_event;
    isWtau_event             = is_Wtau_event;
    isSingleHardPhoton_event = is_SingleHardPhoton_event;
    isdiphoton_event         = is_diphoton_event;
    isisr_photon_event       = is_isr_photon_event;
    
    gen_MetPt                = genMetPt;
    gen_MetPx                = genMetPx;
    gen_MetPy                = genMetPy;
    gen_MetPhi               = genMetPhi;
    gen_MetSumEt             = genMetSumEt;
    delta_phiGEN             = Delta_phiGEN; 
   
    pthat                    = gen_pthat;
    graviton_pt              = gen_gravitonpt;	
    //////////////////////////////////////////////////////
    
    
    //////////////////////////////////////////////////////
    /////to be turned on for data only SSC ///////////////
    //////////////////////////////////////////////////////
    //scraping_isScrapingEvent = Scraping_isScrapingEvent;//
    /////////////////////////////////////////////////////

    
     photon_n                 = Photon_n; 
     jet_n                    = Jet_n;
     pfjet_n                  = pfJet_n;	
     elec_n                   = Electron_n;
     muon_n                   = Muon_n;
     cosmicmuon_n             = CosmicMuon_n;
     vertex_n                 = Vertex_n;
<<<<<<< myCuts_halo.C
     //cscseg_n                 = CSCseg_n;
     tau_n                    = Tau_n;



=======
     //cscseg_n                 = CSCseg_n;
     tau_n                    = Tau_n;

>>>>>>> 1.2
    if(debug_)cout<<"Starting Filling: vertex, photon, ele, muon  and tracks variable "<<endl;

     for(int i=0; i<NFilled;i++){
       vertex_x[i]                            = Vertex_x[i];
       vertex_y[i]                            = Vertex_y[i];
       vertex_z[i]                            = Vertex_z[i];
       vertex_tracksize[i]                    = Vertex_tracksize[i];
       vertex_ndof[i]                         = Vertex_ndof[i];
       vertex_chi2[i]                         = Vertex_chi2[i];
       vertex_d0[i]                           = Vertex_d0[i];
       vertex_isFake[i]                       = Vertex_isFake[i];
       //photon_iLICTD[i]                        = computeLICTD(Photon_ncrys[i], Photon_timing_xtal[i], Photon_energy_xtal[i]);
       photon_E[i]                            = Photon_E[i];
       photon_pt[i]                           = Photon_pt[i];
       photon_eta[i]                          = Photon_eta[i];
       photon_phi[i]                          = Photon_phi[i];
       photon_theta[i]                        = Photon_theta[i];
       photon_et[i]                           = Photon_et[i];
       photon_swissCross[i]                   = Photon_swissCross[i];
       photon_r9[i]                           = Photonr9[i];
       photon_e1x5[i]                         = Photon_e1x5[i];
       photon_e2x5[i]                         = Photon_e2x5[i];
       photon_e3x3[i]                         = Photon_e3x3[i];
       photon_e5x5[i]                         = Photon_e5x5[i];
       photon_r1x5[i]                         = Photon_r1x5[i];
       photon_r2x5[i]                         = Photon_r2x5[i];
       photon_maxEnergyXtal[i]                = Photon_maxEnergyXtal[i];
       photon_SigmaEtaEta[i]                  = Photon_SigmaEtaEta[i];
       photon_SigmaIetaIeta[i]                = Photon_SigmaIetaIeta[i];
       photon_SigmaPhiPhi[i]                  = Photon_SigmaPhiPhi[i];
       photon_SigmaIphiIphi[i]                = Photon_SigmaIphiIphi[i];
       photon_Roundness[i]                    = Photon_Roundness[i];
       photon_Angle[i]                        = Photon_Angle[i];


       photon_ecalRecHitSumEtConeDR03[i]      = Photon_ecalRecHitSumEtConeDR03[i];
       photon_hcalTowerSumEtConeDR03[i]       = Photon_hcalTowerSumEtConeDR03[i];
       photon_trkSumPtSolidConeDR03[i]        = Photon_trkSumPtSolidConeDR03[i];
       photon_trkSumPtHollowConeDR03[i]       = Photon_trkSumPtHollowConeDR03[i];
       photon_nTrkSolidConeDR03[i]            = Photon_nTrkSolidConeDR03[i];
       photon_nTrkHollowConeDR03[i]           = Photon_nTrkHollowConeDR03[i];

       photon_hcalDepth1TowerSumEtConeDR03[i] = Photon_hcalDepth1TowerSumEtConeDR03[i];
       photon_hcalDepth2TowerSumEtConeDR03[i] = Photon_hcalDepth2TowerSumEtConeDR03[i];
       photon_ecalRecHitSumEtConeDR04[i]      = Photon_ecalRecHitSumEtConeDR04[i];
       photon_hcalTowerSumEtConeDR04[i]       = Photon_hcalTowerSumEtConeDR04[i];
       photon_trkSumPtSolidConeDR04[i]        = Photon_trkSumPtSolidConeDR04[i];
       photon_trkSumPtHollowConeDR04[i]       = Photon_trkSumPtHollowConeDR04[i];

       photon_nTrkSolidConeDR04[i]            = Photon_nTrkSolidConeDR04[i];
       photon_nTrkHollowConeDR04[i]           = Photon_nTrkHollowConeDR04[i];
       photon_hcalDepth1TowerSumEtConeDR04[i] = Photon_hcalDepth1TowerSumEtConeDR04[i];
       photon_hcalDepth2TowerSumEtConeDR04[i] = Photon_hcalDepth2TowerSumEtConeDR04[i];
       photon_hasPixelSeed[i]                 = Photon_hasPixelSeed[i];
       photon_isEB[i]                         = Photon_isEB[i];
       photon_isEE[i]                         = Photon_isEE[i];
       photon_isEBGap[i]                      = Photon_isEBGap[i];
       photon_isEEGap[i]                      = Photon_isEEGap[i];
       photon_isEBEEGap[i]                    = Photon_isEBEEGap[i];
       photon_HoE[i]                          = Photon_HoE[i];
       photon_px[i]                           = Photon_px[i];
       photon_py[i]                           = Photon_py[i];
       photon_pz[i]                           = Photon_pz[i];

       TLorentzVector v1;
       v1.SetPxPyPzE(photon_px[i],
                     photon_py[i],
                     photon_pz[i],
                     photon_E[i]); 


       photon_basiccluster_size[i]            = Photon_no_of_basic_clusters[i];
       photon_sc_energy[i]                    = Photon_sc_energy[i];
       photon_sc_eta[i]                       = Photon_sc_eta[i];
       photon_sc_phi[i]                       = Photon_sc_phi[i];
       photon_sc_x[i]                         = Photon_sc_x[i];
       photon_sc_y[i]                         = Photon_sc_y[i];
       photon_sc_z[i]                         = Photon_sc_z[i];
       photon_etaWidth[i]                     = Photon_etaWidth[i];
       photon_phiWidth[i]                     = Photon_phiWidth[i];
       photon_sc_et[i]                        = Photon_sc_et[i];
       photon_ntracks[i]                      = Photon_ntracks[i];

       photon_isconverted[i]                  = Photon_isconverted[i];
       photon_pairInvmass[i]                  = Photon_pairInvmass[i];
       photon_pairCotThetaSeparation[i]       = Photon_pairCotThetaSeperation[i];
       photon_pairmomentumX[i]                = Photon_pairmomentumX[i];
       photon_pairmomentumY[i]                = Photon_pairmomentumY[i];
       photon_pairmomentumZ[i]                = Photon_pairmomentumZ[i];
       photon_EoverP[i]                       = Photon_EoverP[i];
       photon_ConvVx[i]                       = Photon_ConvVx[i];
       photon_ConvVy[i]                       = Photon_ConvVy[i];
       photon_ConvVz[i]                       = Photon_ConvVz[i];
       photon_ZOfPrimaryVertex[i]             = Photon_ZOfPrimaryVertex[i];
       photon_distOfMinimumApproach[i]        = Photon_distOfMinimumApproach[i];
       photon_dPhiTracksAtVtx[i]              = Photon_dPhiTracksAtVtx[i];
       photon_dPhiTracksAtEcal[i]             = Photon_dPhiTracksAtEcal[i];
       photon_dEtaTracksAtEcal[i]             = Photon_dEtaTracksAtEcal[i];

       //photon_ncrys[i]                      = Photon_ncrys[i];
      
       for(int j=0;j<2;j++){	
	 photon_timing_xtal[i][j]             = Photon_timing_xtal[i][j];
	//  photon_ieta_xtalEB[i][j]          = Photon_ieta_xtalEB[i][j];	
	//  photon_iphi_xtalEB[i][j]          = Photon_iphi_xtalEB[i][j];
	//  photon_energy_xtalEB[i][j]        = Photon_energy_xtalEB[i][j];
       }

       //////////////////////////////////////////////////////
       /////to be turned on for MC only SSC /////////////////
       //////////////////////////////////////////////////////


       //Gen Photon
       MC_photon_pt[i]                        = gen_photonpt[i];
       MC_photon_eta[i]                       = gen_photoneta[i];
       MC_photon_phi[i]                       = gen_photonphi[i];
       MC_photon_px[i]                        = gen_photonpx[i]; 
       MC_photon_py[i]                        = gen_photonpy[i];
       MC_photon_pz[i]                        = gen_photonpx[i];
       MC_photon_E[i]                         = gen_photonE[i];
       MC_photon_MotherID[i]                  = gen_photonMotherID[i];
       MC_photon_MotherPt[i]                  = gen_photonMotherPt[i];
       MC_photon_MotherEta[i]                 = gen_photonMotherEta[i];
       MC_photon_MotherPhi[i]                 = gen_photonMotherPhi[i];
       MC_photon_GrandMotherID[i]             = gen_photonGrandmotherID[i];

       TLorentzVetor v2;
       v2.clear();
       v2.SetPxPyPzE(MC_photon_px[i],
                     MC_photon_py[i],
                     MC_photon_pz[i],
                     MC_photon_E[i]); 
 
       photon_deltR_mc[i]=v2.DeltaR(v1); 
 
       //////////////////////////////////////////////////////


       track_px[i]                            = Track_px[i];   
       track_py[i]                            = Track_py[i];   
       track_pz[i]                            = Track_pz[i];   
       track_pt[i]                            = Track_pt[i];   
       track_eta[i]                           = Track_eta[i];   
       track_phi[i]                           = Track_phi[i];   
       
       elec_px[i]                             = Electron_px[i];
       elec_py[i]                             = Electron_py[i];
       elec_pz[i]                             = Electron_pz[i];
       elec_pt[i]                             = Electron_pt[i];
       elec_eta[i]                            = Electron_eta[i];
       elec_phi[i]                            = Electron_phi[i];
       elec_energy[i]                         = Electron_energy[i];
       elec_charge[i]                         = Electron_charge[i]; 
       elec_trkIso[i]                         = Electron_trkIso[i];
       elec_ecalIso[i]                        = Electron_ecalIso[i];
       elec_hcalIso[i]                        = Electron_hcalIso[i];	
       elec_HoE[i]                            = Electron_HoE[i];
       elec_SigmaIetaIeta[i]                  = Electron_SigmaIetaIeta[i];
       elec_dEtaIn[i]                         = Electron_dEtaIn[i];	
       elec_dPhiIn[i]                         = Electron_dPhiIn[i];	
       elec_sc_energy[i]                      = Electron_sc_energy[i];	
       elec_sc_eta[i]                         = Electron_sc_eta[i];	
       elec_sc_phi[i]                         = Electron_sc_phi[i];	

       ///////////////////added by bhawna//////////
       OneProng0Pi0[i]                     =oneProng0Pi0[i];
       OneProng1Pi0[i]                     =oneProng1Pi0[i];
       OneProng2Pi0[i]                     =oneProng2Pi0[i];
       ThreeProng0Pi0[i]                   =threeProng0Pi0[i];
       ThreeProng1Pi0[i]                   =threeProng1Pi0[i];
       GenHadTauPt[i]                      =genHadTauPt[i];
       GenHadTauEta[i]                     =genHadTauEta[i];
       GenHadTauPhi[i]                     =genHadTauPhi[i];
       NPions[i]                           =nPions[i]; 
       NPhotons[i]                         =nPhotons[i];
       NPi0[i]                             =nPi0[i];

       for(int j = 0; j <2; j++){
	 pionPdgId[i][j]                =PionPdgId[i][j];
	 pionPt[i][j]                   =PionPt[i][j];
	 pionEta[i][j]                  =PionEta[i][j];
	 pionPhi[i][j]                  =PionPhi[i][j];
	 pi0PdgId[i][j]                 =Pi0PdgId[i][j];
	 pi0Pt[i][j]                    =Pi0Pt[i][j];
	 pi0Eta[i][j]                   =Pi0Eta[i][j];
	 pi0Phi[i][j]                   =Pi0Phi[i][j];
	 photonPt[i][j]                 =PhotonPt[i][j]; 
	 photonEta[i][j]                =PhotonEta[i][j];
	 photonPhi[i][j]                =PhotonPhi[i][j];
	 photonPdgId[i][j]              =PhotonPdgId[i][j];
       }
       

       muon_px[i]                             = Muon_px[i];
       muon_py[i]                             = Muon_py[i];
       muon_pz[i]                             = Muon_pz[i];
       muon_pt[i]                             = Muon_pt[i];
       muon_eta[i]                            = Muon_eta[i];
       muon_phi[i]                            = Muon_phi[i];
       muon_energy[i]                         = Muon_energy[i];
       muon_charge[i]                         = Muon_charge[i];
       muon_isGlobalMuon[i]                   = Muon_isGlobalMuon[i];
       muon_isTrackerMuon[i]                  = Muon_isTrackerMuon[i];
       muon_isStandAloneMuon[i]               = Muon_isStandAloneMuon[i];
       muon_InnerTrack_isNonnull[i]           = Muon_InnerTrack_isNonnull[i];
       muon_OuterTrack_isNonnull[i]           = Muon_OuterTrack_isNonnull[i];
       muon_OuterTrack_InnerPoint_x[i]        = Muon_OuterTrack_InnerPoint_x[i];
       muon_OuterTrack_InnerPoint_y[i]        = Muon_OuterTrack_InnerPoint_y[i];
       muon_OuterTrack_InnerPoint_z[i]        = Muon_OuterTrack_InnerPoint_z[i];
       muon_OuterTrack_InnerPoint_px[i]       = Muon_OuterTrack_InnerPoint_px[i];
       muon_OuterTrack_InnerPoint_py[i]       = Muon_OuterTrack_InnerPoint_py[i];
       muon_OuterTrack_InnerPoint_pz[i]       = Muon_OuterTrack_InnerPoint_pz[i];
       muon_OuterTrack_OuterPoint_x[i]        = Muon_OuterTrack_OuterPoint_x[i];
       muon_OuterTrack_OuterPoint_y[i]        = Muon_OuterTrack_OuterPoint_y[i];
       muon_OuterTrack_OuterPoint_z[i]        = Muon_OuterTrack_OuterPoint_z[i];
       muon_OuterTrack_OuterPoint_px[i]       = Muon_OuterTrack_OuterPoint_px[i];
       muon_OuterTrack_OuterPoint_py[i]       = Muon_OuterTrack_OuterPoint_py[i];
       muon_OuterTrack_OuterPoint_pz[i]       = Muon_OuterTrack_OuterPoint_pz[i];
       muon_InnerTrack_InnerPoint_x[i]        = Muon_InnerTrack_InnerPoint_x[i];
       muon_InnerTrack_InnerPoint_y[i]        = Muon_InnerTrack_InnerPoint_y[i];
       muon_InnerTrack_InnerPoint_z[i]        = Muon_InnerTrack_InnerPoint_z[i];
       muon_InnerTrack_InnerPoint_px[i]       = Muon_InnerTrack_InnerPoint_px[i];
       muon_InnerTrack_InnerPoint_py[i]       = Muon_InnerTrack_InnerPoint_py[i];
       muon_InnerTrack_InnerPoint_pz[i]       = Muon_InnerTrack_InnerPoint_pz[i];
       muon_InnerTrack_OuterPoint_x[i]        = Muon_InnerTrack_OuterPoint_x[i];
       muon_InnerTrack_OuterPoint_y[i]        = Muon_InnerTrack_OuterPoint_y[i];
       muon_InnerTrack_OuterPoint_z[i]        = Muon_InnerTrack_OuterPoint_z[i];
       muon_InnerTrack_OuterPoint_px[i]       = Muon_InnerTrack_OuterPoint_px[i];
       muon_InnerTrack_OuterPoint_py[i]       = Muon_InnerTrack_OuterPoint_py[i];
       muon_InnerTrack_OuterPoint_pz[i]       = Muon_InnerTrack_OuterPoint_pz[i];  
       muon_InnerTrack_OuterPoint_pz[i]       = Muon_InnerTrack_OuterPoint_pz[i];   
       muon_OuterPoint_x[i]                   = Muon_OuterPoint_x[i];  
       muon_OuterPoint_y[i]                   = Muon_OuterPoint_y[i];
       muon_OuterPoint_z[i]                   = Muon_OuterPoint_z[i];
       muon_InnerPoint_x[i]                   = Muon_InnerPoint_x[i];  
       muon_InnerPoint_y[i]                   = Muon_InnerPoint_y[i];
       muon_InnerPoint_z[i]                   = Muon_InnerPoint_z[i];
       muon_vx[i]                             = Muon_vx[i];  
       muon_vy[i]                             = Muon_vy[i];
       muon_vz[i]                             = Muon_vz[i];
 

       //Cosmic Muon
       cosmicmuon_px[i]                             = CosmicMuon_px[i];
       cosmicmuon_py[i]                             = CosmicMuon_py[i];
       cosmicmuon_pz[i]                             = CosmicMuon_pz[i];
       cosmicmuon_pt[i]                             = CosmicMuon_pt[i];
       cosmicmuon_eta[i]                            = CosmicMuon_eta[i];
       cosmicmuon_phi[i]                            = CosmicMuon_phi[i];
       cosmicmuon_energy[i]                         = CosmicMuon_energy[i];
       cosmicmuon_charge[i]                         = CosmicMuon_charge[i];
       cosmicmuon_isGlobalcosmicmuon[i]             = CosmicMuon_isGlobalMuon[i];
       cosmicmuon_isTrackercosmicmuon[i]            = CosmicMuon_isTrackerMuon[i];
       cosmicmuon_isStandAlonecosmicmuon[i]         = CosmicMuon_isStandAloneMuon[i];
       cosmicmuon_InnerTrack_isNonnull[i]           = CosmicMuon_InnerTrack_isNonnull[i];
       cosmicmuon_OuterTrack_isNonnull[i]           = CosmicMuon_OuterTrack_isNonnull[i];
       cosmicmuon_OuterTrack_InnerPoint_x[i]        = CosmicMuon_OuterTrack_InnerPoint_x[i];
       cosmicmuon_OuterTrack_InnerPoint_y[i]        = CosmicMuon_OuterTrack_InnerPoint_y[i];
       cosmicmuon_OuterTrack_InnerPoint_z[i]        = CosmicMuon_OuterTrack_InnerPoint_z[i];
       cosmicmuon_OuterTrack_InnerPoint_px[i]       = CosmicMuon_OuterTrack_InnerPoint_px[i];
       cosmicmuon_OuterTrack_InnerPoint_py[i]       = CosmicMuon_OuterTrack_InnerPoint_py[i];
       cosmicmuon_OuterTrack_InnerPoint_pz[i]       = CosmicMuon_OuterTrack_InnerPoint_pz[i];
       cosmicmuon_OuterTrack_OuterPoint_x[i]        = CosmicMuon_OuterTrack_OuterPoint_x[i];
       cosmicmuon_OuterTrack_OuterPoint_y[i]        = CosmicMuon_OuterTrack_OuterPoint_y[i];
       cosmicmuon_OuterTrack_OuterPoint_z[i]        = CosmicMuon_OuterTrack_OuterPoint_z[i];
       cosmicmuon_OuterTrack_OuterPoint_px[i]       = CosmicMuon_OuterTrack_OuterPoint_px[i];
       cosmicmuon_OuterTrack_OuterPoint_py[i]       = CosmicMuon_OuterTrack_OuterPoint_py[i];
       cosmicmuon_OuterTrack_OuterPoint_pz[i]       = CosmicMuon_OuterTrack_OuterPoint_pz[i];
       cosmicmuon_InnerTrack_InnerPoint_x[i]        = CosmicMuon_InnerTrack_InnerPoint_x[i];
       cosmicmuon_InnerTrack_InnerPoint_y[i]        = CosmicMuon_InnerTrack_InnerPoint_y[i];
       cosmicmuon_InnerTrack_InnerPoint_z[i]        = CosmicMuon_InnerTrack_InnerPoint_z[i];
       cosmicmuon_InnerTrack_InnerPoint_px[i]       = CosmicMuon_InnerTrack_InnerPoint_px[i];
       cosmicmuon_InnerTrack_InnerPoint_py[i]       = CosmicMuon_InnerTrack_InnerPoint_py[i];
       cosmicmuon_InnerTrack_InnerPoint_pz[i]       = CosmicMuon_InnerTrack_InnerPoint_pz[i];
       cosmicmuon_InnerTrack_OuterPoint_x[i]        = CosmicMuon_InnerTrack_OuterPoint_x[i];
       cosmicmuon_InnerTrack_OuterPoint_y[i]        = CosmicMuon_InnerTrack_OuterPoint_y[i];
       cosmicmuon_InnerTrack_OuterPoint_z[i]        = CosmicMuon_InnerTrack_OuterPoint_z[i];
       cosmicmuon_InnerTrack_OuterPoint_px[i]       = CosmicMuon_InnerTrack_OuterPoint_px[i];
       cosmicmuon_InnerTrack_OuterPoint_py[i]       = CosmicMuon_InnerTrack_OuterPoint_py[i];
       cosmicmuon_InnerTrack_OuterPoint_pz[i]       = CosmicMuon_InnerTrack_OuterPoint_pz[i];  
       cosmicmuon_InnerTrack_OuterPoint_pz[i]       = CosmicMuon_InnerTrack_OuterPoint_pz[i];   
       cosmicmuon_OuterPoint_x[i]                   = CosmicMuon_OuterPoint_x[i];  
       cosmicmuon_OuterPoint_y[i]                   = CosmicMuon_OuterPoint_y[i];
       cosmicmuon_OuterPoint_z[i]                   = CosmicMuon_OuterPoint_z[i];
       


       if(Photon_isEB[i]){       
	 vector<int> thisPho_ietaRH;
	 vector<int> thisPho_iphiRH;
	 vector<float> thisPho_eRH;
	 vector<float> thisPho_tRH;
	 for(int j=0; j<Photon_ncrys[i];j++){
	   thisPho_ietaRH.push_back(Photon_ieta_xtalEB[i][j]);
	   thisPho_iphiRH.push_back(Photon_iphi_xtalEB[i][j]);
	   thisPho_eRH.push_back(Photon_energy_xtal[i][j]);
	   thisPho_tRH.push_back(Photon_timing_xtal[i][j]);
       }
	 
	 photon_e2e9[i]    = bg->simpleBarrelE2E9(Photon_ncrys[0],thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH);
	 IsE2E9Spike[i]    = bg->isE2E9Spike(Photon_ncrys[0], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH);
	 IsTimeSpike[i]    = bg->isTimeSpike(thisPho_tRH);
       
       }

       photon_e4e1[i] = Photon_e4e1[i];
       photon_e6e2[i] = Photon_e6e2[i];

       IsHEHalo[i]       = bg->isHEHalo(Photon_sc_phi[i], HERecHit_subset_n, HERecHit_subset_x, HERecHit_subset_y, HERecHit_subset_energy, HERecHit_subset_time);
       //IsCSCHalo[i]      = bg->isCSCHalo(Photon_sc_phi[i], CSCseg_n, CSCseg_x, CSCseg_y, CSCseg_time);
       IsTrackHalo[i]    = bg->isTrackHalo(Photon_sc_phi[i], CosmicMuon_n, CosmicMuon_OuterTrack_InnerPoint_x, CosmicMuon_OuterTrack_InnerPoint_y);
       /*   for(int j=0;j < cscseg_n;j++)
	 {
	   DeltaPhiCSCHalo[i][j]   = bg->deltaPhiCSCHalo(Photon_sc_phi[i], CSCseg_x[j], CSCseg_y[j]);
	   }*/
      // sigmaRC[i] = bg->sigmaRCosmic(Photon_sc_x[i],  Photon_sc_y[i], Photon_sc_z[i],CosmicMuon_OuterTrack_InnerPoint_x,CosmicMuon_OuterTrack_InnerPoint_y,CosmicMuon_OuterTrack_InnerPoint_z, CosmicMuon_OuterTrack_InnerPoint_px,CosmicMuon_OuterTrack_InnerPoint_py,CosmicMuon_OuterTrack_InnerPoint_pz, CosmicMuon_n);
       fisherd2[i] = bg->fisherCosmic2(Photon_Roundness[i],Photon_Angle[i]);
       fisherd3[i] = bg->fisherCosmic3(Photon_Roundness[i],Photon_Angle[i],Photon_hasConvTrk[i]);


     }

    
   //Fill BeamHaloSummary

  is_BeamHaloIDTightPass = isBeamHaloIDTightPass;
  is_BeamHaloIDLoosePass  = isBeamHaloIDLoosePass;
         
  is_BeamHaloEcalLoosePass =isBeamHaloEcalLoosePass;
  is_BeamHaloHcalLoosePass =isBeamHaloHcalLoosePass;
  is_BeamHaloCSCLoosePass  =isBeamHaloCSCLoosePass;
  is_BeamHaloGlobalLoosePass =isBeamHaloGlobalLoosePass;
         
  is_BeamHaloEcalTightPass = isBeamHaloEcalTightPass;
  is_BeamHaloHcalTightPass  =isBeamHaloHcalTightPass;
  is_BeamHaloCSCTightPass  = isBeamHaloCSCTightPass;
  is_BeamHaloGlobalTightPass =isBeamHaloGlobalTightPass;

  is_SmellsLikeHalo_Tag  = isSmellsLikeHalo_Tag;
  is_LooseHalo_Tag       = isLooseHalo_Tag;
  is_TightHalo_Tag       = isTightHalo_Tag;
  is_ExtremeTightHalo_Tag =isExtremeTightHalo_Tag;

 
  //RPChit info

 
     ///////////////////////////////////////////////////////////////////////////////// 
     //***Check photon_ Isoalitons : Want to match jets with only Isoalted photon***//
     /////////////////////////////////////////////////////////////////////////////////
     if(debug_)cout<<"Removing Photon from Jet Block"<<endl;

     if(photon_ecalRecHitSumEtConeDR04[0]<4.2+0.006*photon_pt[0]   && 
        photon_hcalTowerSumEtConeDR04[0]< 2.2+0.0025*photon_pt[0]  &&
        photon_trkSumPtHollowConeDR04[0] < 3.5+0.001*photon_pt[0]  &&
        photon_HoE[0]<0.05 && photon_n >0)
     {
       int tmpFilled=0;
       for (int j=0;j<Jet_n;++j){
	 Float_t dPhi = fabs(photon_phi[0]-Jet_phi[j]);
	 if (dPhi > TMath::Pi()) dPhi = TMath::Pi()*2. - dPhi;
	 Float_t dEta = fabs(photon_eta[0]-Jet_eta[j]);
	 Float_t dR = sqrt(dPhi*dPhi + dEta*dEta);

	 ///////////////////////////////////////////////////////////////////////////////// 
	 //**************Check Delta_R and Jet ID Cuts then fill jets*******************//
	 ///////////////////////////////////////////////////////////////////////////////// 
	 if (dR > 0.5 && fabs(Jet_eta[j])<2.6 && 
             Jet_emEnergyFraction[j]>0.01     &&
             Jet_n90Hits[j]>1 && Jet_fHPD[j]<0.98)
         {
	   jet_pt[tmpFilled]                    = Jet_pt[j];
	   jet_px[tmpFilled]                    = Jet_px[j];
	   jet_py[tmpFilled]                    = Jet_py[j];
	   jet_pz[tmpFilled]                    = Jet_pz[j];
	   jet_E[tmpFilled]                     = Jet_E[j];
	   jet_eta[tmpFilled]                   = Jet_eta[j];
	   jet_phi[tmpFilled]                   = Jet_phi[j];
	   jet_energyFractionHadronic[tmpFilled]= Jet_energyFractionHadronic[j];
	   jet_EMenergyFraction[tmpFilled]      = Jet_emEnergyFraction[j]; 
	   jet_hitsInN90[tmpFilled]             = Jet_hitsInN90[j]; 
	   jet_n90Hits[tmpFilled]               = Jet_n90Hits[j];
	   jet_nTowers[tmpFilled]               = Jet_nTowers[j];
	   jet_fHPD[tmpFilled]                  = Jet_fHPD[j];
	   jet_fRBX[tmpFilled]                  = Jet_fRBX[j];
	   jet_RHF[tmpFilled]                   = Jet_RHF[j];
	   tmpFilled++;
	  if(tmpFilled >= NJFilled ) break;
	 }

       }
     }

    if(debug_)cout<<"Again doing the same "<<endl;
    ///Same above loop but to calculate HT for different threshold of Jets pT
           if(photon_ecalRecHitSumEtConeDR04[0]<4.2+0.006*photon_pt[0]   && 
              photon_hcalTowerSumEtConeDR04[0]< 2.2+0.0025*photon_pt[0]  &&
              photon_trkSumPtHollowConeDR04[0]< 3.5+0.001*photon_pt[0]   &&
              photon_HoE[0]<0.05                                         &&
              photon_n >0)
     {
       int tmpFilled=0;
       for (int j=0;j<Jet_n;++j){
         Float_t dPhi = fabs(photon_phi[0]-Jet_phi[j]);
         if (dPhi > TMath::Pi()) dPhi = TMath::Pi()*2. - dPhi;
          
         Float_t dEta = fabs(photon_eta[0]-Jet_eta[j]);
         Float_t dR = sqrt(dPhi*dPhi + dEta*dEta);

         ///////////////////////////////////////////////////////////////////////////////// 
         //**************Check Delta_R and Jet ID Cuts then fill jets*******************//
         ///////////////////////////////////////////////////////////////////////////////// 
         if (dR > 0.5 && fabs(Jet_eta[j])<2.6 && 
             Jet_emEnergyFraction[j]>0.01     && 
             Jet_n90Hits[j]>1                 &&
             Jet_fHPD[j]<0.98)
            {
                             HT += Jet_pt[j];
             if(Jet_pt[j]> 50.0)HT_JetpT50 += Jet_pt[j];
             if(Jet_pt[j]> 100.)HT_JetpT100 += Jet_pt[j];
             if(Jet_pt[j]> 150.)HT_JetpT150 += Jet_pt[j];  

             jet_n_FromHT++;

        }
       }
     }

   if(Jet_n>0)jet_n_PhotonRemoved=Jet_n-jet_n_FromHT;


///-----------Same for pfJets
     if(debug_)cout<<"Removing Photon from pfJet Block"<<endl;

     if(photon_ecalRecHitSumEtConeDR04[0]<4.2+0.006*photon_pt[0]   &&
        photon_hcalTowerSumEtConeDR04[0]< 2.2+0.0025*photon_pt[0]  &&
        photon_trkSumPtHollowConeDR04[0] < 3.5+0.001*photon_pt[0]  &&
        photon_HoE[0]<0.05 && photon_n >0)
     {
       int tmpFilledPF=0;
       for (int j=0;j<pfJet_n;++j){
	 Float_t pfdPhi = fabs(photon_phi[0]-pfJet_phi[j]);
	 if (pfdPhi > TMath::Pi()) pfdPhi = TMath::Pi()*2. - pfdPhi;
	 Float_t pfdEta = fabs(photon_eta[0]-pfJet_eta[j]);
	 Float_t pfdR = sqrt(pfdPhi*pfdPhi + pfdEta*pfdEta);

	 if (pfdR > 0.5 && fabs(pfJet_eta[j])<2.6){
	   pfjet_pt[tmpFilledPF]                    = pfJet_pt[j];
	   pfjet_px[tmpFilledPF]                    = pfJet_px[j];
	   pfjet_py[tmpFilledPF]                    = pfJet_py[j];
	   pfjet_pz[tmpFilledPF]                    = pfJet_pz[j];
	   pfjet_E[tmpFilledPF]                     = pfJet_E[j];
	   pfjet_eta[tmpFilledPF]                   = pfJet_eta[j];
	   pfjet_phi[tmpFilledPF]                   = pfJet_phi[j];
	   tmpFilledPF++;
	  if(tmpFilledPF >= NJFilled ) break;
	 }

       }
     }//check if photon is isolated or not

       ///Same above loop but to calculate HT for different threshold of Jets pT
       if(photon_ecalRecHitSumEtConeDR04[0]<4.2+0.006*photon_pt[0]  && 
          photon_hcalTowerSumEtConeDR04[0]< 2.2+0.0025*photon_pt[0] && 
          photon_trkSumPtHollowConeDR04[0] < 3.5+0.001*photon_pt[0] &&
          photon_HoE[0]<0.05 && photon_n >0)
         {
          int tmpFilled=0; 
           for(int j=0;j<pfJet_n;++j)
              {  Float_t dPhi = fabs(photon_phi[0]-pfJet_phi[j]);
                 if (dPhi > TMath::Pi()) dPhi = TMath::Pi()*2. - dPhi;
                     Float_t dEta = fabs(photon_eta[0]-pfJet_eta[j]);
                     Float_t dR = sqrt(dPhi*dPhi + dEta*dEta);
                    
         /////////////////////////////////////////////////////////////////////////////////
          //**************Check Delta_R and Jet ID Cuts then fill jets*******************//
         /////////////////////////////////////////////////////////////////////////////////
          if (dR > 0.5 && fabs(pfJet_eta[j])<5.0)
             {
                  pfHT += pfJet_pt[j];
                 if(pfJet_pt[j]> 50.0)HT_pfJetpT50 += pfJet_pt[j];
                 if(pfJet_pt[j]> 100.)HT_pfJetpT100 += pfJet_pt[j];
                 if(pfJet_pt[j]> 150.)HT_pfJetpT150 += pfJet_pt[j]; 
                  pfjet_n_FromHT++;
            }
          }//loop over j
       }//check if an isolated photon           

    if(pfJet_n > 0)pfjet_n_PhotonRemoved = pfJet_n - pfjet_n_FromHT;
                    
  //-------------------


     if(debug_)cout<<"starting CSC infor filling"<<endl;
          
    //csc 
     /*
     for(int j=0;j < cscseg_n;j++)
       {
	 cscseg_x[j] = CSCseg_x[j];
	 cscseg_y[j] = CSCseg_y[j];
	 cscseg_z[j] = CSCseg_z[j];
	 cscseg_time[j] = CSCseg_time[j];
       }
     */
     //lets calculate Ht based on the 
     
    if(debug_)cout<<"Starting filling MET"<<endl;
    for(int i =0; i<6;i++)
    {
      calo_MetPt[i]                = CaloMetPt[i];
      calo_MetPx[i]                = CaloMetPx[i];
      calo_MetPy[i]                = CaloMetPy[i];
      calo_MetPhi[i]               = CaloMetPhi[i];
      calo_MetSumEt[i]             = CaloMetSumEt[i];
      pf_MetPt[i]                  = PFMetPt[i];
      pf_MetPx[i]                  = PFMetPx[i];
      pf_MetPy[i]                  = PFMetPy[i];
      pf_MetPhi[i]                 = PFMetPhi[i];
      pf_MetSumEt[i]               = PFMetSumEt[i]; 
      tc_MetPt[i]                  = TCMetPt[i];
      tc_MetPx[i]                  = TCMetPx[i];
      tc_MetPy[i]                  = TCMetPy[i];
      tc_MetPhi[i]                 = TCMetPhi[i];
      tc_MetSumEt[i]               = TCMetSumEt[i];
    }
    calo_MetSigma                  = CaloMetSigma;
    //calo_MetCorr                 = CaloMetCorr;
    calo_MetEz                     = CaloMetEz;
    calo_EtFractionHadronic        = CaloEtFractionHadronic;
    calo_EmEtFraction              = CaloEmEtFraction;
    calo_HadEtInHB                 = CaloHadEtInHB;
    calo_HadEtInHE                 = CaloHadEtInHE;
    calo_HadEtInHO                 = CaloHadEtInHO;
    calo_HadEtInHF                 = CaloHadEtInHF;
    calo_EmEtInEB                  = CaloEmEtInEB;
    calo_EmEtInEE                  = CaloEmEtInEE;
    calo_EmEtInHF                  = CaloEmEtInHF;
    calo_MaxEtInEmTowers           = CaloMaxEtInEmTowers;
    calo_MaxEtInHadTowers          = CaloMaxEtInHadTowers;
    delta_phi                      = Delta_phi;
    delta_phiPF                    = Delta_phiPF; 
    delta_phiTC                    = Delta_phiTC;
    
    if(debug_)cout<<" Fill the tree: Entry ="<<jentry<<"/"<<nentries<<endl;
    outfile->cd();
    if(debug_)cout<<"now Fill"<<endl;
//    if(photon_n > 0 && (photon_pt[0]> 30.0 && photon_HoE[0]< 0.05)){

    tree->Fill();
                                                                // }
     
  }//Event Loop
  
  //Write the Tree
  if(debug_)cout<<"Write root files now"<<endl;
  tree->Write();
  cout<<"written the tree"<<endl;
  outfile->Write() ;	
  outfile->Close() ;
  
}//Loop()  


/*
float makeHistogramsW_7::computeLICTD(int nCry, float cryTime[], float cryEnergy[]){
    int seed = GetSeed(nCry, cryEnergy);
    if(seed<0) return -99;   
    Float_t LICTD =0;        
    Int_t crysCrys=-1;       
    Int_t crysThresh=0;      
    for(int k=0;k<nCry&&k<100;++k){
        if (seed==k) continue;
        Float_t crysE = cryEnergy[k];
        if(crysE > 1.){      
            crysThresh++;    
            Float_t tdiff = cryTime[seed] - cryTime[k];
            if(fabs(tdiff) > fabs(LICTD)){
                LICTD = tdiff;
                crysCrys=k;  
            }                
        }                    
    }                        
    return LICTD;            
}                            

//  seed finding function    
//useage makeHistogramsW_7::GetSeed(Photon_ncrys[x], Photon_energy_xtal[x])
int makeHistogramsW_7::GetSeed(int Ncry, float cryE[]){
    int seedIdx=-1;          
    float seedE=-99;         
    for(int k=0;k<Ncry&&k<100;++k){
        if(cryE[k] > seedE){ 
            seedIdx = k;     
            seedE=cryE[seedIdx];
        }                    
    }                        
    return seedIdx;          
}                            

*/

