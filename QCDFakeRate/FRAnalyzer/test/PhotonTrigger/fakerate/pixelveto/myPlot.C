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


#ifdef __MAKECINT__
#pragma link C++ class map<TString,TH1D*>+;
#endif



void myPlot::Loop()
{

  using namespace std;
  

 //-------------Various Options------------
  bool removeHLTMatch = false;  


  //Triggers to use for selection
  string SelectedTrigger[3]={  "HLT_Photon30_CaloIdVL_v",
                               "HLT_Photon75_CaloIdVL_v", 
                               "HLT_DoublePhoton33_v" };

 float TriggThrMin[3] ={  0.0,
                         80.0,  
                          0.0 
                        };

 float TriggThrMax[3] ={ 80.0,
                         10000., 
                         80.0
                        };

  int NTriggers    = 3;

  bool UseUnPrescales = false;
  bool UseScaleFactor = false;
  double weight = 1; 


  //Output Files
  TString region   = "fakerate_pixel_photonTrigger"; 
 string myregion   = (string)region;
  Long64_t maxsize = 500000000;                  //100GB
  maxsize *= 100;                                //to bypass some compiler limitations with big constants
  TTree::SetMaxTreeSize(maxsize);
  TFile *rfile = new TFile(("/uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_4_2_3/src/QCDFakeRate/FRAnalyzer/test/PhotonTrigger/fakerate/pixelveto/"+myregion+".root").c_str(),"RECREATE");
  rfile->cd();


   //---------------------------------







 //Thinks which remain same 
  double pi = 4.0*atan(1.0);
  TString i_file   = "data_";
  string filename = "data";

  //define histos here
  std::map<TString, TH1D*> histmap;
  std::map<TString, TH1D*> histmapratio;
  std::map<TString, TH2D*> histmap2d;



  //various plotting ranges  
  int etarange;
  TString xvar;
  TString title;
  int nbins;
  double xmin, xmax;
  string xtitle;
  string ytitle;


  //Get the Variable to plot
  ifstream infile_xvar;
  char cxvariable[200]; 
  std::vector<string> xvariable;
  infile_xvar.open("xvariables.list", ifstream::in );

  while(!infile_xvar.eof()){
    infile_xvar >> cxvariable;
    if(strncmp(cxvariable,"#",1)==0)
      {
	continue;
      }
    xvariable.push_back((string)cxvariable);
  }//while(!infile_xvar.eof())

  infile_xvar.close();
  



      
 //----------Define all pt bins and histo labels----------     

 
  double lower[12];
  for(int ivar=0; ivar<xvariable.size(); ivar++)
     {
	  xvar = (TString)xvariable[ivar];
	
        if(xvar == "phopt" || xvar == "frpt" || xvar == "phoptnocuts" )
	  {
	    title  = "pt of photons" ;
	    xtitle = "pt";
	    nbins  = 10;
	    lower[0] = 20.0; 
            lower[1] = 30.0; 
            lower[2] = 40.0; 
            lower[3] = 50.0; 
            lower[4] = 60.0; 
            lower[5] = 70.0; 
            lower[6] = 80.0; 
            lower[7] = 120.0; 
            lower[8] = 160.0; 
	    lower[9] = 300.0;
            lower[10] = 1000.0; 
	    ytitle = "No of entries";
	  }


	  if(xvar == "phoeta" || xvar == "phoetaisEB" || xvar == "phoetaisEE" || xvar == "freta")
	    {
	      title  = "#eta of photons" ;
	      xtitle = "#eta";
	      nbins  = 24;
	      xmin   = -3.0;
	      xmax   = 3.0;
	      Double_t binwidth = (xmax-xmin)/nbins;
	      char bin_width[10];
	      sprintf(bin_width,"%f",binwidth);
	      string binw = bin_width;
	      ytitle = "No of entries/"+binw;
	  }
	
	  
	
	if(xvar == "ecaliso")
          {
            title  = "ecaliso" ;
            xtitle = "ecaliso(GeV)";
            nbins  = 100;
            xmin   = 0.0;
            xmax   = 50.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }

	if(xvar == "hcaliso")
          {
            title  = "hcaliso" ;
            xtitle = "hcaliso(GeV)";
            nbins  = 100;
            xmin   = 0.0;
            xmax   = 50.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }

	if(xvar == "trkiso")
          {
            title  = "trkiso" ;
            xtitle = "trkiso(GeV)";
            nbins  = 100;
            xmin   = 0.0;
            xmax   = 50.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }

	if(xvar == "hovere")
          {
            title  = "hovere" ;
            xtitle = "hovere";
            nbins  = 100;
            xmin   = 0.0;
            xmax   = 50.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }

	if(xvar == "spikevar")
          {
            title  = "spikevar" ;
            xtitle = "spikevar";
            nbins  = 100;
            xmin   = 10.0;
            xmax   = 2.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }
	
	if(xvar == "seedtime")
          {
            title  = "seedtime" ;
            xtitle = "seedtime";
            nbins  = 100;
            xmin   = 10.0;
            xmax   = 2.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }

	if(  xvar == "sigmaieta" || xvar == "sigietaisEB" || xvar == "sigietaisEE" ||
             xvar == "sigietaisEB_afterspikeRem" || xvar == "sigietaisEE_afterspikeRem" ||
             xvar == "sigietaisEB_onlyspikeVar" || xvar == "sigietaisEE_onlyspikeVar"     )
          {
            title  = "sigmaieta" ;
            xtitle = "sigmaieta";
            nbins  = 100;
            xmin   = 0.0;
            xmax   = 0.08;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }
	
	
	if(xvar == "deltar")
          {
            title  = "deltar" ;
            xtitle = "deltar";
            nbins  = 100;
            xmin   = 0.0;
            xmax   = 10.0;
            Double_t binwidth = (xmax-xmin)/nbins;
            char bin_width[10];
            sprintf(bin_width,"%f",binwidth);
            string binw = bin_width;
            ytitle = "No of entries/"+binw;
          }
	
	if(xvar!="freta" && xvar!="frpt" && xvar!="phopt" && xvar!="phoptnocuts" )
	  histmap[i_file+xvariable[ivar]]     = new TH1D((filename+"_"+xvariable[ivar]).c_str(),(xvariable[ivar]).c_str(),nbins,xmin,xmax);
	
	
	if(xvar == "freta")
	  {
	    histmap[i_file+"etanum"]             = new TH1D((filename+"_"+"eta_num").c_str(),"numerator  ",nbins,xmin,xmax);
	    histmap[i_file+"etanum"]->Sumw2();
	    histmap[i_file+"etadeno"]            = new TH1D((filename+"_"+"eta_deno").c_str(),"denominator  ",nbins,xmin,xmax);
	    histmap[i_file+"etadeno"]->Sumw2();
	    histmapratio[i_file+"eta"]           = new TH1D((filename+"_"+"eta_ratio").c_str(),"fake rate of Vs #eta of #gamma  ",nbins,xmin,xmax);
	    histmapratio[i_file+"eta"]->Sumw2();
	  }


	if(xvar == "phopt" || xvar == "phoptnocuts")
	  histmap[i_file+xvariable[ivar]]        = new TH1D((filename+"_"+xvariable[ivar]).c_str(),(xvariable[ivar]).c_str(),nbins,lower);

	
	if(xvar == "frpt")
          {
            histmap[i_file+"ptnum"]              = new TH1D((filename+"_"+"pt_num").c_str(),"numerator  ",nbins,lower);
            histmap[i_file+"ptnum"]->Sumw2();
            histmap[i_file+"ptdeno"]             = new TH1D((filename+"_"+"pt_deno").c_str(),"denominator  ",nbins,lower);
            histmap[i_file+"ptdeno"]->Sumw2();
	    histmap[i_file+"ptnumnoweight"]      = new TH1D((filename+"_"+"ptnoweight_num").c_str(),"numerator(no weight) ",nbins,lower);
	    histmap[i_file+"ptdenonoweight"]     = new TH1D((filename+"_"+"ptnoweight_deno").c_str(),"denominator(no weight) ",nbins,lower);
            histmapratio[i_file+"pt"]            = new TH1D((filename+"_"+"pt_ratio").c_str(),"fake rate of Vs #pt of #gamma  ",nbins,lower);
            histmapratio[i_file+"pt"]->Sumw2();
	    
	    //EB
	    histmap[i_file+"ptEBnum"]            = new TH1D((filename+"_"+"ptEB_num").c_str(),"numerator : EB ",nbins,lower);
            histmap[i_file+"ptEBnum"]->Sumw2();
            histmap[i_file+"ptEBdeno"]           = new TH1D((filename+"_"+"ptEB_deno").c_str(),"denominator : EB ",nbins,lower);
            histmap[i_file+"ptEBdeno"]->Sumw2();
            histmapratio[i_file+"ptEB"]          = new TH1D((filename+"_"+"ptEB_ratio").c_str(),"fake rate of Vs #pt of #gamma : EB ",nbins,lower);
            histmapratio[i_file+"ptEB"]->Sumw2();
	    histmap[i_file+"ptEBnumnoweight"]    = new TH1D((filename+"_"+"ptEBnoweight_num").c_str(),"numerator(no weight) ",nbins,lower);
	    histmap[i_file+"ptEBdenonoweight"]   = new TH1D((filename+"_"+"ptEBnoweight_deno").c_str(),"denominator(no weight) : EB ",nbins,lower);

	    //EE
	    histmap[i_file+"ptEEnum"]            = new TH1D((filename+"_"+"ptEE_num").c_str(),"numerator : EE ",nbins,lower);
            histmap[i_file+"ptEEnum"]->Sumw2();
            histmap[i_file+"ptEEdeno"]           = new TH1D((filename+"_"+"ptEE_deno").c_str(),"denominator : EE ",nbins,lower);
            histmap[i_file+"ptEEdeno"]->Sumw2();
            histmapratio[i_file+"ptEE"]          = new TH1D((filename+"_"+"ptEE_ratio").c_str(),"fake rate of Vs #pt of #gamma : EE ",nbins,lower);
            histmapratio[i_file+"ptEE"]->Sumw2();
	    histmap[i_file+"ptEEnumnoweight"]    = new TH1D((filename+"_"+"ptEEnoweight_num").c_str(),"numerator(no weight) ",nbins,lower);
	    histmap[i_file+"ptEEdenonoweight"]   = new TH1D((filename+"_"+"ptEEnoweight_deno").c_str(),"denominator(no weight) : EE ",nbins,lower);
		    
          }
	
	
	if(xvar!="freta" && xvar!="frpt"){
	    histmap[i_file+xvariable[ivar]]->GetXaxis()->SetTitle(xtitle.c_str());
	    histmap[i_file+xvariable[ivar]]->GetYaxis()->SetTitle(ytitle.c_str());
	    histmap[i_file+xvariable[ivar]]->GetYaxis()->SetTitleOffset(1.3);
	  }//if(xvar !="freta" && xvar!="frpt")
	
  }//for(int ivar=0; ivar<xvariable.size; ivar++)


//--------Histos for Templates
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




      int mypatphoisEB;
      int mypatphoisEE;
     

//-----------------Event Loop sstart here---------------------------------
 if (fChain == 0) return;
  
    int bugnum=0;  //for tracking deno<num entries

      cout<<"before entrie loop in C"<<endl;
      Long64_t nentries = fChain->GetEntriesFast();

      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {

        cout<<"Loading Tree "<<endl;

	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;
        cout<<"Get Entry "<<endl;
        
        	
       int  goodVertex=0;
       for(int v=0;v< Vertex_n; v++){
          if(  (fabs(Vertex_z[v]) < 24.0) && 
               (Vertex_ndof[v] > 4)       &&
               (!Vertex_isFake[v])        && 
               (fabs(Vertex_d0[0])< 2.0 )
            ){
              goodVertex++;
              }

         }

     

       //check if any of the desired trigger is fired or not
        float  ScaleFactor=1;
               isPassed = isTriggerPassed(SelectedTrigger, NTriggers , ScaleFactor, UseUnPrescales, UseScaleFactor);
               if(UseScaleFactor)weight= ScaleFactor;
	
  
     for(int ivar=0; ivar<xvariable.size(); ivar++)
     {
	      xvar = (TString)xvariable[ivar];
	      double wei = 1  ;
	      
	      if(xvar == "phoptnocuts")
		{
                 for(int ii=0; ii<Photon_n; ii++)
		 histmap[i_file+"phoptnocuts"]->Fill(Photon_pt[ii],wei);
		}

     if(isPassed && goodVertex >=1){

        for(int ipho=0; ipho<Photon_n && ipho< 200; ipho++)
 	   {  
             bool keepThisPhoton=true;
             if(removeHLTMatch)
                     {keepThisPhoton=isHLTMatch(SelectedTrigger, NTriggers ,ipho);}

                  
           //this will remove those photon which matches with hltbjects to remove the trigger bias
        if(keepThisPhoton)
           {
               //cout<<"C: Scale Factor is = "<<weight<<endl;  

		if(Photon_pt[ipho]>20.0)
                   {
				      mypatphoisEB = 0;
				      mypatphoisEE = 0;
				      
				      mypatphoisEB = Photon_isEB[ipho];
				      mypatphoisEE = Photon_isEE[ipho];
				      
				      
				      //spike removal cut
				      
				      int spikeveto = Photon_e2e9[ipho]>0.90 && Photon_swissCross[ipho]<0.95;
				      
				      if( xvar == "phoetaisEB" && mypatphoisEB)
					histmap[i_file+xvariable[ivar]]->Fill(Photon_eta[ipho],weight);
				      if( xvar == "sigietaisEB" && mypatphoisEB)
					histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
				      if( xvar == "phoetaisEE" && mypatphoisEE)
					histmap[i_file+xvariable[ivar]]->Fill(Photon_eta[ipho],weight);
				      if( xvar == "sigietaisEE" && mypatphoisEE)
					histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
				      
				      
                           if( fabs(Photon_eta[ipho])<=2.5 )
			      {
					  if( fabs(Photon_eta[ipho])>2.5 )
					    {
					      cout<<"WARNING!eta out of required range"<<endl;
					      cout<<"i_file = "<<i_file<<" & eta = "<<Photon_eta[ipho]<< endl;
					    }
					 
                          /*                
                                         //use Triggers for specific pt-bin of photon
                                          bool TrigExistPassed=false;

                                          for(int th=0; th < NTriggers; th++){
                                             if(Photon_pt[ipho] >= TriggThrMin[th] && Photon_pt[ipho]< TriggThrMax[th] )
                                               { TrigExistPassed = CheckThisTrigger(SelectedTrigger[th]);}
                
                                               if(TrigExistPassed) continue;
                                           }//for loop

                                   if(!TrigExistPassed) continue; // come out of eta < 2.5 cut 
                           */
					  //no cuts
					  if(xvar == "ecaliso")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_ecalRecHitSumEtConeDR04[ipho],weight);
					  if(xvar == "hcaliso")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_hcalTowerSumEtConeDR04[ipho],weight);
					  if(xvar == "trkiso")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_trkSumPtHollowConeDR04[ipho],weight);
					  if(xvar == "phopt")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_pt[ipho],weight);
					  if(xvar == "phoeta")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_eta[ipho],weight);
					  if(xvar == "spikevar")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_swissCross[ipho],weight);
					  if(xvar == "seedtime")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_timing_xtal[ipho][0],weight);
					  if(xvar == "hovere")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_HoE[ipho],weight);
					  if(xvar == "sigmaieta")
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
					  
					  // tight cuts
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
					  
					  //after subsequent cuts 
					  int fecal = Photon_ecalRecHitSumEtConeDR04[ipho]<(tecalcut);
					  int fhcal = Photon_hcalTowerSumEtConeDR04[ipho]<(thcalcut);
					  int ftrk  = Photon_trkSumPtHollowConeDR04[ipho]<(ttrkcut);
					  int fhoe  = Photon_HoE[ipho]<thoecut ;
					  int fsigieta = 0;
					  
					  if(mypatphoisEB)
					    fsigieta = Photon_SigmaIetaIeta[ipho]<tsigietacut;
					  
					  if(mypatphoisEE)
					    fsigieta = Photon_SigmaIetaIeta[ipho]<tsigietacut;
					  
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
					  
					  
					  //sigmaieta-ieta after spike removal
					  if(xvar == "sigietaisEB_afterspikeRem" && mypatphoisEB && spikeveto)
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
					  if(xvar == "sigietaisEE_afterspikeRem" && mypatphoisEE && spikeveto)
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
					  
					  //sigmaieta after spike var and no time cut       
					  if(mypatphoisEB && xvar == "sigietaisEB_onlyspikeVar" && Photon_swissCross[ipho]<0.95)
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
					  if(mypatphoisEE && xvar == "sigietaisEE_onlyspikeVar" && Photon_swissCross[ipho]<0.95)
					    histmap[i_file+xvariable[ivar]]->Fill(Photon_SigmaIetaIeta[ipho],weight);
					  
					 
					  int ftrkflip = ((Photon_ecalRecHitSumEtConeDR04[ipho]>lecalcut) || 
                                                          (Photon_hcalTowerSumEtConeDR04[ipho]>lhcalcut)  ||
                                                          (Photon_trkSumPtHollowConeDR04[ipho]>ltrkcut)   || 
                                                           Photon_SigmaIetaIeta[ipho]>lsigmaietacut         );
					  
					  
					  
					  //loose photonid
					  int floosephoid = 0;
					  int ftightphoid = 0;
					  
					  floosephoid =(Photon_HoE[ipho]<lminhoecut            &&
                                                        Photon_ecalRecHitSumEtConeDR04[ipho]<lminecalcut  &&
                                                        Photon_hcalTowerSumEtConeDR04[ipho]<lminhcalcut   && 
                                                        Photon_trkSumPtHollowConeDR04[ipho]<lmintrkcut );
					  
					  ftightphoid = ( Photon_HoE[ipho]<thoecut             && 
                                                          Photon_ecalRecHitSumEtConeDR04[ipho]<(tecalcut) &&
                                                          Photon_hcalTowerSumEtConeDR04[ipho]<(thcalcut)  && 
                                                          Photon_trkSumPtHollowConeDR04[ipho]<(ttrkcut)   && 
                                                          Photon_SigmaIetaIeta[ipho]<tsigietacut          &&
                                                         (cutpix==0)                          );
					  
					  
					  
					  
					//Fill denominator  
					  if( ftrkflip  && floosephoid )
					    {
					      if(mypatphoisEB)
						{
						  if(spikeveto)
						    {
						      if(xvar == "freta")
							histmap[i_file+"etadeno"]->Fill(Photon_eta[ipho],weight);
						      
						      if(xvar == "frpt")
							{
							  histmap[i_file+"ptdeno"]->Fill(Photon_pt[ipho],weight);
							  histmap[i_file+"ptdenonoweight"]->Fill(Photon_pt[ipho]);
							}
						      
						      //EB
						      if(xvar == "frpt")
							{
							  histmap[i_file+"ptEBdeno"]->Fill(Photon_pt[ipho],weight);
							  histmap[i_file+"ptEBdenonoweight"]->Fill(Photon_pt[ipho]);
							}
						    }//if(spikeveto)
						}//if(mypatphoisEB)
					      


					      if(mypatphoisEE)
						{
						  if(xvar == "freta")
						    histmap[i_file+"etadeno"]->Fill(Photon_eta[ipho],weight);
						  if(xvar == "frpt")
						    {
						      histmap[i_file+"ptdeno"]->Fill(Photon_pt[ipho],weight);
						      histmap[i_file+"ptdenonoweight"]->Fill(Photon_pt[ipho]);
						    } 
						  //EE
						  if(xvar == "frpt")
						    {
						      histmap[i_file+"ptEEdeno"]->Fill(Photon_pt[ipho],weight);
						      histmap[i_file+"ptEEdenonoweight"]->Fill(Photon_pt[ipho]);
						    }
						}//if(mypatphoisEE)
					    }//if( floosephoid... )

					  

					 //Fill Numerator 
					  if( ftightphoid )
					    {
					      if(mypatphoisEB)
						{
						  if(spikeveto)
						    {
						      if(xvar == "freta")
							histmap[i_file+"etanum"]->Fill(Photon_eta[ipho],weight);
						      if(xvar == "frpt")
							{
							  histmap[i_file+"ptnum"]->Fill(Photon_pt[ipho],weight);
							  histmap[i_file+"ptnumnoweight"]->Fill(Photon_pt[ipho]);
							}
						      
						      //EB
						      if(xvar == "frpt")
							{
							  histmap[i_file+"ptEBnum"]->Fill(Photon_pt[ipho],weight);
							  histmap[i_file+"ptEBnumnoweight"]->Fill(Photon_pt[ipho]);
							}
						    }//if(spikeveto)
						}//if(mypatphoisEB)
					      
					      if(mypatphoisEE)
						{
						  if(xvar == "freta")
						    histmap[i_file+"etanum"]->Fill(Photon_eta[ipho],weight);
						  if(xvar == "frpt")
						    {
						      histmap[i_file+"ptnum"]->Fill(Photon_pt[ipho],weight);
						      histmap[i_file+"ptnumnoweight"]->Fill(Photon_pt[ipho]);
						    }
						  
						  //EE
						  if(xvar == "frpt")
						    {
						      histmap[i_file+"ptEEnum"]->Fill(Photon_pt[ipho],weight);
						      histmap[i_file+"ptEEnumnoweight"]->Fill(Photon_pt[ipho]);
						    }
						}//if(mypatphoisEE)
					      
					    }//if( ftightphoid )



        //-------------Add here the template Methods------------------

         if(xvar==0)//<--this runs once
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
	      
			      
	      int tight = Photon_HoE[ipho]<thoecut && 
                          Photon_ecalRecHitSumEtConeDR04[ipho]<(tecalcut) && 
                          Photon_hcalTowerSumEtConeDR04[ipho]<(thcalcut) && 
                          (cutpix==0);

	     int tmpftrkflip = (Photon_trkSumPtHollowConeDR04[ipho]>ltrkcut) && 
                                (Photon_trkSumPtHollowConeDR04[ipho]<(5*ltrkcut));
 
	      if( tmpftrkflip  && tight)
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
            }//---if(xvar==0)  Fills only once


                    }//if( fabs(Photon_eta[ipho])<=2.5 )

		}//if(Photon_pt[ipho]>20.0)


           }//Keep this photo


      }//for(int ipho=0; ipho<Photon_n; ipho++)
                      			    

    }//ifPassed
	
  }//for(int ivar=0; ivar<xvariable.size; ivar++)


 if((jentry%100000)==0) cout<<"Event Analyzerd Till Now  =  "<<(jentry+1)<<endl;

}//for (Long64_t jentry=0; jentry<nentries;jentry++)
  
      double integ_deno = 0;
      double integ_num  = 0;
      double eff = 0;

      double integ_denoEB = 0;
      double integ_numEB  = 0;
      double effEB = 0;
      
      double integ_denoEE = 0;
      double integ_numEE  = 0;
      double effEE = 0;

 

//--------------Set Bin content for  raw fake rate-----------------


      //----------->eta
      int nbins_eta = histmap[i_file+"etanum"]->GetNbinsX();
      //cout<<"nbins_eta = "<<nbins_eta<<endl;
      for(int ibin=1; ibin<=nbins_eta; ibin++)
	{
	  //cout<<"ibin = "<<ibin<<endl;
	  
	  integ_deno  = histmap[i_file+"etadeno"]->GetBinContent(ibin);

	  integ_num = histmap[i_file+"etanum"]->GetBinContent(ibin);
	 
	  
	  if(integ_deno < integ_num)
	    {
	      std::cout<<"WARNING!!! deno < num in eta and bin = "<<ibin<<endl;
	      cout<<"integ_deno =  "<<integ_deno<<endl;
	      cout<<"integ_num =  "<<integ_num<<endl;
	      cout<<"bugnum = "<<bugnum<<endl;
	    }
	  
	  eff   = integ_num/TMath::Max(0.001,integ_deno);
	  
	  double binerror = 0; 
	  //binerror  = ( (integ_num+1)*(integ_num+2) )/( (integ_deno+2)*(integ_deno+3)) - pow( (integ_num+1)/(integ_deno+2),2 );
	  binerror  = sqrt( (eff*(1-eff))/TMath::Max(0.001,integ_deno) );
	  
	  histmapratio[i_file+"eta"]->SetBinContent(ibin,eff);
	  histmapratio[i_file+"eta"]->SetBinError(ibin,binerror); 
	}//for(int ibin=1; ibin<=nbins; ibin++)


      //--------->pt
      int nbins_pt = histmap[i_file+"ptnum"]->GetNbinsX();
	for(int ibin=1; ibin<=nbins_pt; ibin++)
	  {
	    //cout<<"ibin = "<<ibin<<endl;
	    
	    integ_deno  = histmap[i_file+"ptdeno"]->GetBinContent(ibin);

	    integ_num = histmap[i_file+"ptnum"]->GetBinContent(ibin);
	    
	    if(integ_deno < integ_num)
	      {
		std::cout<<"WARNING!!! deno < num in pt and bin = "<<ibin<<endl;
		cout<<"integ_deno =  "<<integ_deno<<endl;
		cout<<"integ_num =  "<<integ_num<<endl;
		cout<<"bugnum = "<<bugnum<<endl;
	      }

	    
	    eff   = integ_num/TMath::Max(0.001,integ_deno);
	    
	    double binerror = 0; 
	    //binerror  = ( (integ_num+1)*(integ_num+2) )/( (integ_deno+2)*(integ_deno+3)) - pow( (integ_num+1)/(integ_deno+2),2 );
	    binerror  = sqrt( (eff*(1-eff))/TMath::Max(0.001,integ_deno) );
	    
	    histmapratio[i_file+"pt"]->SetBinContent(ibin,eff);
	    histmapratio[i_file+"pt"]->SetBinError(ibin,binerror); //
	    
	    //EB
	    integ_denoEB  = histmap[i_file+"ptEBdeno"]->GetBinContent(ibin);
	    integ_numEB   = histmap[i_file+"ptEBnum"]->GetBinContent(ibin);
	    if(integ_denoEB < integ_numEB)
	      {
		std::cout<<"WARNING!!! deno < num in pt and bin = "<<ibin<<endl;
		cout<<"integ_deno =  "<<integ_denoEB<<endl;
		cout<<"integ_num =  "<<integ_numEB<<endl;
		cout<<"bugnum = "<<bugnum<<endl;
	      }
	    
	    effEB   = integ_numEB/TMath::Max(0.001,integ_denoEB);
	    
	    double binerrorEB = 0; 
	    binerrorEB  = sqrt( (effEB*(1-effEB))/TMath::Max(0.001,integ_denoEB) );
	    histmapratio[i_file+"ptEB"]->SetBinContent(ibin,effEB);
	    histmapratio[i_file+"ptEB"]->SetBinError(ibin,binerrorEB); //
	    
	    //EE
	    integ_denoEE  = histmap[i_file+"ptEEdeno"]->GetBinContent(ibin);
	    integ_numEE   = histmap[i_file+"ptEEnum"]->GetBinContent(ibin);
	    if(integ_denoEE < integ_numEE)
	      {
		std::cout<<"WARNING!!! deno < num in pt and bin = "<<ibin<<endl;
		cout<<"integ_deno =  "<<integ_denoEE<<endl;
		cout<<"integ_num =  "<<integ_numEE<<endl;
		cout<<"bugnum = "<<bugnum<<endl;
	      }
	    
	    effEE   = integ_numEE/TMath::Max(0.001,integ_denoEE);
	    
	    double binerrorEE = 0; 
	    binerrorEE  = sqrt( (effEE*(1-effEE))/TMath::Max(0.001,integ_denoEE) );
	    histmapratio[i_file+"ptEE"]->SetBinContent(ibin,effEE);
	    histmapratio[i_file+"ptEE"]->SetBinError(ibin,binerrorEE); //
	    
	  }//for(int ibin=1; ibin<=nbins; ibin++)
	

//--------Now write the raw fake rate histos here----------
	rfile->cd();
	
	for(int ivar=0; ivar<xvariable.size(); ivar++)
	    {
	      xvar = (TString)xvariable[ivar];
	      
	      if(xvar!="freta" && xvar!="frpt")
		histmap[i_file+xvariable[ivar]]->Write();
	      
	      if(xvar=="freta")
		{
		  histmap[i_file+"etanum"]->Write();
		  histmap[i_file+"etadeno"]->Write();
		  histmapratio[i_file+"eta"]->Write();
		}
	      
	      if(xvar=="frpt")
		{
		  histmap[i_file+"ptnum"]->Write();
		  histmap[i_file+"ptdeno"]->Write();
		  histmapratio[i_file+"pt"]->Write();
		  histmap[i_file+"ptnumnoweight"]->Write();
		  histmap[i_file+"ptdenonoweight"]->Write();

		  //EB
		  histmap[i_file+"ptEBnum"]->Write();
		  histmap[i_file+"ptEBdeno"]->Write();
		  histmapratio[i_file+"ptEB"]->Write();
		  histmap[i_file+"ptEBnumnoweight"]->Write();
		  histmap[i_file+"ptEBdenonoweight"]->Write();
		  //EE
		  histmap[i_file+"ptEEnum"]->Write();
		  histmap[i_file+"ptEEdeno"]->Write();
		  histmapratio[i_file+"ptEE"]->Write();
		  histmap[i_file+"ptEEnumnoweight"]->Write();
		  histmap[i_file+"ptEEdenonoweight"]->Write();
		}
	    }//for(int ivar=0; ivar<xvariable.size; ivar++)
	  

  //-----------Now Write the Template histos	
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
	 


  rfile->Write();
  rfile->Close();
  std::cout<<"done ! "<<std::endl;


}//void loop

