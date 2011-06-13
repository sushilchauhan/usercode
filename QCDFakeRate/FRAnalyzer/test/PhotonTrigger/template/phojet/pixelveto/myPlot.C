#define myPlot_cxx
#include "myPlot.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include "TMath.h"
#include <fstream>
#include <vector>
#include <map>
#include "TColor.h"
#include "TLegend.h"
#include "stdio.h"
#include <string>
#include "TRFIOFile.h"


void myPlot::Loop()
{
  using namespace std;

   if (fChain == 0) return;

   TFile *f = new TFile("/uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_4_2_3/src/QCDFakeRate/FRAnalyzer/test/PhotonTrigger/template/phojet/pixelveto_monophoton/pixelphojet_e2e9.root","RECREATE");

   //EB                                                                                              
   TH1D *eb20to30 = new TH1D("sigmaietaEB_eb20to30","sigmaietaieta  for pt:20to30 : EB",80,0.0,0.04);
   TH1D *eb30to40 = new TH1D("sigmaietaEB_eb30to40","sigmaietaieta  for pt:30to40 : EB",80,0.0,0.04);
   TH1D *eb40to50 = new TH1D("sigmaietaEB_eb40to50","sigmaietaieta  for pt:40to50 : EB",80,0.0,0.04);
   TH1D *eb50to60 = new TH1D("sigmaietaEB_eb50to60","sigmaietaieta  for pt:50to60 : EB",80,0.0,0.04);
   TH1D *eb60to70 = new TH1D("sigmaietaEB_eb60to70","sigmaietaieta  for pt:60to70 : EB",80,0.0,0.04);
   TH1D *eb70to80 = new TH1D("sigmaietaEB_eb70to80","sigmaietaieta  for pt:70to80 : EB",80,0.0,0.04);
   TH1D *eb80to120 = new TH1D("sigmaietaEB_eb80to120","sigmaietaieta  for pt:80to120 : EB",80,0.0,0.04);
   TH1D *eb120to160 = new TH1D("sigmaietaEB_eb120to160","sigmaietaieta  for pt:120to160 : EB",80,0.0,0.04);
   TH1D *eb160to300 = new TH1D("sigmaietaEB_eb160to300","sigmaietaieta  for pt:160to300 : EB",80,0.0,0.04);
   TH1D *eb300to1000 = new TH1D("sigmaietaEB_eb300to1000","sigmaietaieta  for pt:300to1000 : EB",80,0.0,0.04);

   //EE                                                                                               
   TH1D *ee20to30 = new TH1D("sigmaietaEE_ee20to30","sigmaietaieta  for pt:30to35 : EE",120,0.0,0.06);
   TH1D *ee30to40 = new TH1D("sigmaietaEE_ee30to40","sigmaietaieta  for pt:30to40 : EE",120,0.0,0.06);
   TH1D *ee40to50 = new TH1D("sigmaietaEE_ee40to50","sigmaietaieta  for pt:40to50 : EE",120,0.0,0.06);
   TH1D *ee50to60 = new TH1D("sigmaietaEE_ee50to60","sigmaietaieta  for pt:50to60 : EE",120,0.0,0.06);
   TH1D *ee60to70 = new TH1D("sigmaietaEE_ee60to70","sigmaietaieta  for pt:60to70 : EE",120,0.0,0.06);
   TH1D *ee70to80 = new TH1D("sigmaietaEE_ee70to80","sigmaietaieta  for pt:70to80 : EE",120,0.0,0.06);
   TH1D *ee80to120 = new TH1D("sigmaietaEE_ee80to120","sigmaietaieta  for pt:80to120 : EE",120,0.0,0.06);
   TH1D *ee120to160 = new TH1D("sigmaietaEE_ee120to160","sigmaietaieta  for pt:120to160 : EE",120,0.0,0.06);
   TH1D *ee160to300 = new TH1D("sigmaietaEE_ee160to300","sigmaietaieta  for pt:160to300 : EE",120,0.0,0.06);
   TH1D *ee300to1000 = new TH1D("sigmaietaEE_ee300to1000","sigmaietaieta  for pt:300to1000 : EE",120,0.0,0.06);


   //FLIPPED TRK FOR DATA
   //EB
   //
   TH1D *teb20to30 = new TH1D("tsigmaietaEB_eb20to30","sigmaietaieta  for pt:20to30 : EB",80,0.0,0.04);
   TH1D *teb30to40 = new TH1D("tsigmaietaEB_eb30to40","sigmaietaieta  for pt:30to40 : EB",80,0.0,0.04);
   TH1D *teb40to50 = new TH1D("tsigmaietaEB_eb40to50","sigmaietaieta  for pt:40to50 : EB",80,0.0,0.04);
   TH1D *teb50to60 = new TH1D("tsigmaietaEB_eb50to60","sigmaietaieta  for pt:50to60 : EB",80,0.0,0.04);
   TH1D *teb60to70 = new TH1D("tsigmaietaEB_eb60to70","sigmaietaieta  for pt:60to70 : EB",80,0.0,0.04);
   TH1D *teb70to80 = new TH1D("tsigmaietaEB_eb70to80","sigmaietaieta  for pt:70to80 : EB",80,0.0,0.04);
   TH1D *teb80to120 = new TH1D("tsigmaietaEB_eb80to120","sigmaietaieta  for pt:80to120 : EB",80,0.0,0.04);
   TH1D *teb120to160 = new TH1D("tsigmaietaEB_eb120to160","sigmaietaieta  for pt:120to160 : EB",80,0.0,0.04);
   TH1D *teb160to300 = new TH1D("tsigmaietaEB_eb160to300","sigmaietaieta  for pt:160to300 : EB",80,0.0,0.04);
   TH1D *teb300to1000 = new TH1D("tsigmaietaEB_eb300to1000","sigmaietaieta  for pt:300to1000 : EB",80,0.0,0.04);

   //EE                                                                       
   TH1D *tee20to30 = new TH1D("tsigmaietaEE_ee20to30","sigmaietaieta  for pt:30to35 : EE",120,0.0,0.06);
   TH1D *tee30to40 = new TH1D("tsigmaietaEE_ee30to40","sigmaietaieta  for pt:30to40 : EE",120,0.0,0.06);
   TH1D *tee40to50 = new TH1D("tsigmaietaEE_ee40to50","sigmaietaieta  for pt:40to50 : EE",120,0.0,0.06);
   TH1D *tee50to60 = new TH1D("tsigmaietaEE_ee50to60","sigmaietaieta  for pt:50to60 : EE",120,0.0,0.06);
   TH1D *tee60to70 = new TH1D("tsigmaietaEE_ee60to70","sigmaietaieta  for pt:60to70 : EE",120,0.0,0.06);
   TH1D *tee70to80 = new TH1D("tsigmaietaEE_ee70to80","sigmaietaieta  for pt:70to80 : EE",120,0.0,0.06);
   TH1D *tee80to120 = new TH1D("tsigmaietaEE_ee80to120","sigmaietaieta  for pt:80to120 : EE",120,0.0,0.06);
   TH1D *tee120to160 = new TH1D("tsigmaietaEE_ee120to160","sigmaietaieta  for pt:120to160 : EE",120,0.0,0.06);
   TH1D *tee160to300 = new TH1D("tsigmaietaEE_ee160to300","sigmaietaieta  for pt:160to300 : EE",120,0.0,0.06);
   TH1D *tee300to1000 = new TH1D("tsigmaietaEE_ee300to1000","sigmaietaieta  for pt:300to1000 : EE",120,0.0,0.06);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(ientry == 0)
	{
	  string filename = fChain->GetCurrentFile()->GetName();
	  std::cout<<"filename = "<<filename<<std::endl;
	}
      
      
      for(int ipho=0; ipho<Photon_n && ipho< 200; ipho++)
	{
	  int mypatphoisEB = 0;
	  int mypatphoisEE = 0;
	  
	  mypatphoisEB = Photon_isEB[ipho];
	  mypatphoisEE = Photon_isEE[ipho];
	  
	  //int spikeveto = (patphoseedtime[ipho])<3.5 && Photon_swissCross[ipho]<0.95;
          int spikeveto = Photon_e2e9[ipho]<0.95 && Photon_swissCross[ipho]<0.95; 	
		  
	  if( fabs(Photon_eta[ipho])<=2.5 )
	    {
	      //----- tight cuts                                                                                                 
	      double tecalcut = 4.2+0.006*Photon_pt[ipho];
	      double thcalcut = 2.2+0.0025*Photon_pt[ipho];
	      double ttrkcut  = 2.0+0.001*Photon_pt[ipho];
	      double tsigietacut = 0;
	      if(mypatphoisEB)
		tsigietacut = 0.013;
	      
	      if(mypatphoisEE)
		tsigietacut = 0.03;
	      
	      double thoecut = 0.05;
	      int cutpix  = Photon_hasPixelSeed[ipho];                                                                                       
	      //int cutpix  = 0; //not using for the time being so this will pass my numerator always                                           
	      
	      //----- loose cuts                                                                                                 
	      double lecalcut = 4.2+0.006*Photon_pt[ipho];
	      double lhcalcut = 2.2+0.0025*Photon_pt[ipho];
	      double ltrkcut  = 3.5+0.001*Photon_pt[ipho];
	      double lsigmaietacut;
	      if(mypatphoisEB)
		lsigmaietacut = 0.013;
	      
	      if(mypatphoisEE)
		lsigmaietacut = 0.03;
	      
	      double lminecalcut = TMath::Min( 5*lecalcut,0.2*Photon_pt[ipho] );
	      double lminhcalcut = TMath::Min( 5*lhcalcut,0.2*Photon_pt[ipho] );
	      double lmintrkcut  = TMath::Min( 5*ltrkcut,0.2*Photon_pt[ipho] ) ;
	      double lminhoecut  = 0.05;
	      
	      int ftrkflip = (Photon_trkSumPtHollowConeDR04[ipho]>ltrkcut) &&
                             (Photon_trkSumPtHollowConeDR04[ipho]<(5*ltrkcut));
	      
	      
	      int ftightphoid = Photon_HoE[ipho]<thoecut && 
                                Photon_ecalRecHitSumEtConeDR04[ipho]<(tecalcut) && 
                                Photon_hcalTowerSumEtConeDR04[ipho]<(thcalcut) &&
                                Photon_trkSumPtHollowConeDR04[ipho]<(ttrkcut) && (cutpix==0);

             //Make sure the photon's mother is not a jet	      

              bool truePhoton=true;
 
	      if(truePhoton )
		{
		  if(ftightphoid)
		    {
		      if(mypatphoisEB)
			{
			  if(spikeveto)
			    {
			      if(Photon_pt[ipho]>=20. && Photon_pt[ipho]<30.)
				eb20to30->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=30. && Photon_pt[ipho]<40.)
				eb30to40->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      
			      if(Photon_pt[ipho]>=40. && Photon_pt[ipho]<50.)
				eb40to50->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=50. && Photon_pt[ipho]<60.)
				eb50to60->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=60. && Photon_pt[ipho]<70.)
				eb60to70->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=70. && Photon_pt[ipho]<80.)
				eb70to80->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			      if(Photon_pt[ipho]>=80. && Photon_pt[ipho]<120.)
				eb80to120->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=120. && Photon_pt[ipho]<160.)
				eb120to160->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=160. && Photon_pt[ipho]<300.)
				eb160to300->Fill(Photon_SigmaIetaIeta[ipho]);
			      
			      if(Photon_pt[ipho]>=300. && Photon_pt[ipho]<1000.)
				eb300to1000->Fill(Photon_SigmaIetaIeta[ipho]);

			      
			    }//if(spikeveto)
			}//if(mypatphoisEB)
		      
		      if(mypatphoisEE)
			{
			  if(Photon_pt[ipho]>=20. && Photon_pt[ipho]<30.)
			    ee20to30->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  
			  if(Photon_pt[ipho]>=30. && Photon_pt[ipho]<40.)
			    ee30to40->Fill(Photon_SigmaIetaIeta[ipho]);
			  
		      
			  if(Photon_pt[ipho]>=40. && Photon_pt[ipho]<50.)
			    ee40to50->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=50. && Photon_pt[ipho]<60.)
			    ee50to60->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=60. && Photon_pt[ipho]<70.)
			    ee60to70->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=70. && Photon_pt[ipho]<80.)
			    ee70to80->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=80. && Photon_pt[ipho]<120.)
			    ee80to120->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=120. && Photon_pt[ipho]<160.)
			    ee120to160->Fill(Photon_SigmaIetaIeta[ipho]);
		      
			  if(Photon_pt[ipho]>=160. && Photon_pt[ipho]<300.)
			    ee160to300->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=300. && Photon_pt[ipho]<1000.)
			    ee300to1000->Fill(Photon_SigmaIetaIeta[ipho]);
		      
			  
			}//if(mypatphoisEE) 
		    }//if(ftightphoid)
		}//if( patphomGenpdgId[ipho] == 22 && (patphomGenmompdgId[ipho]!=21) && !(patphonGenmompdgId[ipho]<10) )
			      
	      int tight = Photon_HoE[ipho]<thoecut && 
                          Photon_ecalRecHitSumEtConeDR04[ipho]<(tecalcut) && 
                          Photon_hcalTowerSumEtConeDR04[ipho]<(thcalcut) && 
                          (cutpix==0);
	      
	      if( ftrkflip  && tight)
		{
		  
		  if(mypatphoisEB)
		    {
		      if(spikeveto)
			{
			  
			  if(Photon_pt[ipho]>=20. && Photon_pt[ipho]<30.)
			    teb20to30->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  
			  if(Photon_pt[ipho]>=30. && Photon_pt[ipho]<40.)
			    teb30to40->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  
			  if(Photon_pt[ipho]>=40. && Photon_pt[ipho]<50.)
			    teb40to50->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=50. && Photon_pt[ipho]<60.)
			    teb50to60->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=60. && Photon_pt[ipho]<70.)
			    teb60to70->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=70. && Photon_pt[ipho]<80.)
			    teb70to80->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=80. && Photon_pt[ipho]<120.)
			    teb80to120->Fill(Photon_SigmaIetaIeta[ipho]);

			  if(Photon_pt[ipho]>=120. && Photon_pt[ipho]<160.)
			    teb120to160->Fill(Photon_SigmaIetaIeta[ipho]);
			  
			  if(Photon_pt[ipho]>=160. && Photon_pt[ipho]<300.)
			    teb160to300->Fill(Photon_SigmaIetaIeta[ipho]);

			  if(Photon_pt[ipho]>=300. && Photon_pt[ipho]<1000.)
			    teb300to1000->Fill(Photon_SigmaIetaIeta[ipho]);



			  
			}//if(spikeveto)
		    }//if(mypatphoisEB) 
		  
		  if(mypatphoisEE)
		    {
		      if(Photon_pt[ipho]>=20. && Photon_pt[ipho]<30.)
			tee20to30->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      
		      if(Photon_pt[ipho]>=30. && Photon_pt[ipho]<40.)
			tee30to40->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      
		      if(Photon_pt[ipho]>=40. && Photon_pt[ipho]<50.)
			tee40to50->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=50. && Photon_pt[ipho]<60.)
			tee50to60->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=60. && Photon_pt[ipho]<70.)
			tee60to70->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=70. && Photon_pt[ipho]<80.)
			tee70to80->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=80. && Photon_pt[ipho]<120.)
			tee80to120->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=120. && Photon_pt[ipho]<160.)
			tee120to160->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=160. && Photon_pt[ipho]<300.)
			tee160to300->Fill(Photon_SigmaIetaIeta[ipho]);
		      
		      if(Photon_pt[ipho]>=300. && Photon_pt[ipho]<1000.)
			tee300to1000->Fill(Photon_SigmaIetaIeta[ipho]);
		      


		    }//if(mypatphoisEE) 
		}//if( ftrkflip  && floosephoid )
	      
	    }//if( fabs(Photon_eta[ipho])<=2.5 ).
	  
	}//for(int ipho=0; ipho<patphosize; ipho++)
   }//for (Long64_t jentry=0; jentry<nentries;jentry++)

   f->cd();
   eb20to30->Write();
   eb30to40->Write();
   eb40to50->Write();
   eb50to60->Write();
   eb60to70->Write();
   eb70to80->Write();
   eb80to120->Write();
   eb120to160->Write();
   eb160to300->Write();
   eb300to1000->Write();



   ee20to30->Write();
   ee30to40->Write();
   ee40to50->Write();
   ee50to60->Write();
   ee60to70->Write();
   ee70to80->Write();
   ee80to120->Write();
   ee120to160->Write();
   ee160to300->Write();
   ee300to1000->Write();


   teb20to30->Write();
   teb30to40->Write();
   teb40to50->Write();
   teb50to60->Write();
   teb60to70->Write();
   teb70to80->Write();
   teb80to120->Write();
   teb120to160->Write();
   teb160to300->Write();
   teb300to1000->Write();


   tee20to30->Write();
   tee30to40->Write();
   tee40to50->Write();
   tee50to60->Write();
   tee60to70->Write();
   tee70to80->Write();
   tee80to120->Write();
   tee120to160->Write();
   tee160to300->Write();
   tee300to1000->Write();


   f->Write();
   f->Close();
}//void myPlot::Loop()
