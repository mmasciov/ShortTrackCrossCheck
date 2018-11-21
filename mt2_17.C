#define mt2_17_cxx
#include "mt2_17.h"
#include "goodrun.cc"
#include "dorky.cc"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;
using namespace duplicate_removal;

float maxmt2=200;
float isoSTC = 6;
float qualSTC = 3;

TFile VetoFile("/nfs-6/userdata/dpgilber/MT2_Inputs/VetoHists.root");
TH2F* veto_bar = (TH2F*) VetoFile.Get("h_VetoEtaPhi_bar");
TH2F* veto_ecp = (TH2F*) VetoFile.Get("h_VetoEtaPhi_ecp");
TH2F* veto_ecn = (TH2F*) VetoFile.Get("h_VetoEtaPhi_ecn");

int InEtaPhiVetoRegion(float eta, float phi) {
  if (fabs(eta) > 2.4) return -1;
  else if (eta > 1.4) {
    float bc = veto_ecp->GetBinContent(veto_ecp->FindBin(eta,phi));
    if (bc > 0) return 3;
  } else if (eta < -1.4) {
    float bc = veto_ecn->GetBinContent(veto_ecn->FindBin(eta,phi));
    if (bc > 0) return 2;
  } else {
    float bc = veto_bar->GetBinContent(veto_bar->FindBin(eta,phi));
    if (bc > 0) return 1;
  }
  if (eta < -1.05 && eta > -1.15 && phi < -1.8 && phi > -2.1) return 4;  
  return 0;
}

void getYield()
{

  set_goodrun_file("Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1_snt.txt");
  
  TChain* chain =new TChain("mt2");
  chain->Add("/nfs-6/userdata/dpgilber/mt2babies/data_2017_loose/data_Run2017*.root");

  mt2_17 t(chain);
  
  long int nentries = chain->GetEntries();
  std::cout << "nentries: " << nentries << std::endl;
  
  int nDuplicates;
  int ne_stF23=0;
  int ne_stcF23=0;
  int ne_stF4=0;
  int ne_stcF4=0;
  int ne_stV23lHt=0;
  int ne_stcV23lHt=0;
  int ne_stV4lHt=0;
  int ne_stcV4lHt=0;
  int ne_stV23hHt=0;
  int ne_stcV23hHt=0;
  int ne_stV4hHt=0;
  int ne_stcV4hHt=0;

  int ne_PstF23=0;
  int ne_PstcF23=0;
  int ne_PstF4=0;
  int ne_PstcF4=0;
  int ne_PstV23lHt=0;
  int ne_PstcV23lHt=0;
  int ne_PstV4lHt=0;
  int ne_PstcV4lHt=0;
  int ne_PstV23hHt=0;
  int ne_PstcV23hHt=0;
  int ne_PstV4hHt=0;
  int ne_PstcV4hHt=0;

  int ne_MstF23=0;
  int ne_MstcF23=0;
  int ne_MstF4=0;
  int ne_MstcF4=0;
  int ne_MstV23lHt=0;
  int ne_MstcV23lHt=0;
  int ne_MstV4lHt=0;
  int ne_MstcV4lHt=0;
  int ne_MstV23hHt=0;
  int ne_MstcV23hHt=0;
  int ne_MstV4hHt=0;
  int ne_MstcV4hHt=0;

  int ne_LstF23=0;
  int ne_LstcF23=0;
  int ne_LstF4=0;
  int ne_LstcF4=0;
  int ne_LstV23lHt=0;
  int ne_LstcV23lHt=0;
  int ne_LstV4lHt=0;
  int ne_LstcV4lHt=0;
  int ne_LstV23hHt=0;
  int ne_LstcV23hHt=0;
  int ne_LstV4hHt=0;
  int ne_LstcV4hHt=0;

  int nt_PstF23=0;
  int nt_PstcF23=0;
  int nt_PstF4=0;
  int nt_PstcF4=0;
  int nt_PstV23lHt=0;
  int nt_PstcV23lHt=0;
  int nt_PstV4lHt=0;
  int nt_PstcV4lHt=0;
  int nt_PstV23hHt=0;
  int nt_PstcV23hHt=0;
  int nt_PstV4hHt=0;
  int nt_PstcV4hHt=0;

  int nt_MstF23=0;
  int nt_MstcF23=0;
  int nt_MstF4=0;
  int nt_MstcF4=0;
  int nt_MstV23lHt=0;
  int nt_MstcV23lHt=0;
  int nt_MstV4lHt=0;
  int nt_MstcV4lHt=0;
  int nt_MstV23hHt=0;
  int nt_MstcV23hHt=0;
  int nt_MstV4hHt=0;
  int nt_MstcV4hHt=0;

  int nt_LstF23=0;
  int nt_LstcF23=0;
  int nt_LstF4=0;
  int nt_LstcF4=0;
  int nt_LstV23lHt=0;
  int nt_LstcV23lHt=0;
  int nt_LstV4lHt=0;
  int nt_LstcV4lHt=0;
  int nt_LstV23hHt=0;
  int nt_LstcV23hHt=0;
  int nt_LstV4hHt=0;
  int nt_LstcV4hHt=0;
//  int ne_st=0;
//  int ne_stc=0;
  for(long int e = 0; e<nentries; ++e){
    t.GetEntry(e);
    
    if(e%100000==0)
      std::cout<<"Entry "<<e<<" / "<<nentries<<std::endl;
    
    DorkyEventIdentifier id(t.run, t.evt, t.lumi);
    if (is_duplicate(id) ){
      ++nDuplicates;
      continue;
    }

    //Golden JSON
    if(!goodrun(t.run, t.lumi)) continue;

    // MET filters
    if(!(t.Flag_eeBadScFilter && t.Flag_globalSuperTightHalo2016Filter && t.Flag_goodVertices && t.Flag_HBHENoiseFilter && t.Flag_HBHENoiseIsoFilter && t.Flag_EcalDeadCellTriggerPrimitiveFilter && t.Flag_ecalBadCalibFilter && t.Flag_badMuonFilter && t.Flag_badChargedCandidateFilter)) continue;    

    // All triggers
    //    if(!((t.HLT_PFHT1050 || t.HLT_PFHT890_Prescale || t.HLT_PFHT780_Prescale || t.HLT_PFHT680_Prescale || t.HLT_PFHT590_Prescale || t.HLT_PFHT510_Prescale || t.HLT_PFHT430_Prescale || t.HLT_PFHT370_Prescale || t.HLT_PFHT250_Prescale || t.HLT_PFHT180_Prescale) || (t.HLT_PFHT500_PFMET100_PFMHT100 || t.HLT_PFMET120_PFMHT120 || t.HLT_PFMETNoMu120_PFMHTNoMu120 || t.HLT_PFMETNoMu120_PFMHTNoMu120_PFHT60 || t.HLT_PFHT800_PFMET75_PFMHT75))) continue;
    // Only prescaled triggers
    if(!((t.HLT_PFHT1050 || t.HLT_PFHT890_Prescale || t.HLT_PFHT780_Prescale || t.HLT_PFHT680_Prescale || t.HLT_PFHT590_Prescale || t.HLT_PFHT510_Prescale || t.HLT_PFHT430_Prescale || t.HLT_PFHT370_Prescale || t.HLT_PFHT250_Prescale || t.HLT_PFHT180_Prescale))) continue;
    // Only un-prescaled triggers
    //if(!((t.HLT_PFHT1050) || (t.HLT_PFHT500_PFMET100_PFMHT100 || t.HLT_PFMET120_PFMHT120 || t.HLT_PFMETNoMu120_PFMHTNoMu120 || t.HLT_PFMETNoMu120_PFMHTNoMu120_PFHT60 || t.HLT_PFHT800_PFMET75_PFMHT75))) continue;
    
    // Event pre-selection
    if(t.nVert<=0) continue;
    if(t.nJet30<=1) continue;
    if(t.nJet30FailId>=1) continue;
    if(t.ht<250) continue;
    if(t.met_pt<30) continue;
    if(t.mt2<60) continue;
    if(t.deltaPhiMin<0.3) continue;
    if(t.diffMetMht/t.met_pt>0.5) continue;
    if(t.met_miniaodPt/t.met_caloPt>=5.0 || t.nJet200MuFrac50DphiMet>0) continue;
    if(t.mt2>maxmt2) continue;

    // Track selection
    int nst=0;
    int nstc=0;

    int ntpe_PstF23=0;
    int ntpe_PstcF23=0;
    int ntpe_PstF4=0;
    int ntpe_PstcF4=0;
    int ntpe_PstV23lHt=0;
    int ntpe_PstcV23lHt=0;
    int ntpe_PstV4lHt=0;
    int ntpe_PstcV4lHt=0;
    int ntpe_PstV23hHt=0;
    int ntpe_PstcV23hHt=0;
    int ntpe_PstV4hHt=0;
    int ntpe_PstcV4hHt=0;
    
    int ntpe_MstF23=0;
    int ntpe_MstcF23=0;
    int ntpe_MstF4=0;
    int ntpe_MstcF4=0;
    int ntpe_MstV23lHt=0;
    int ntpe_MstcV23lHt=0;
    int ntpe_MstV4lHt=0;
    int ntpe_MstcV4lHt=0;
    int ntpe_MstV23hHt=0;
    int ntpe_MstcV23hHt=0;
    int ntpe_MstV4hHt=0;
    int ntpe_MstcV4hHt=0;
    
    int ntpe_LstF23=0;
    int ntpe_LstcF23=0;
    int ntpe_LstF4=0;
    int ntpe_LstcF4=0;
    int ntpe_LstV23lHt=0;
    int ntpe_LstcV23lHt=0;
    int ntpe_LstV4lHt=0;
    int ntpe_LstcV4lHt=0;
    int ntpe_LstV23hHt=0;
    int ntpe_LstcV23hHt=0;
    int ntpe_LstV4hHt=0;
    int ntpe_LstcV4hHt=0;

    for (int it = 0; it<t.ntracks; ++it){
      
      if(t.track_pt[it]<15.0) continue;
      if(std::fabs(t.track_eta[it])>2.4) continue;

      bool CaloSel = !(t.track_DeadECAL[it] || t.track_DeadHCAL[it]) && InEtaPhiVetoRegion(t.track_eta[it],t.track_phi[it]) == 0;
      if(!CaloSel) continue;
      bool extraCaloSel = (!(t.track_eta[it] < -0.7 && t.track_eta[it] > -0.9 && t.track_phi[it] > 1.5 && t.track_phi[it] < 1.7 )
			   && !(t.track_eta[it] < 0.30 && t.track_eta[it] > 0.10 && t.track_phi[it] > 2.2 && t.track_phi[it] < 2.5 )
			   && !(t.track_eta[it] < 0.50 && t.track_eta[it] > 0.40 && t.track_phi[it] > -0.7 && t.track_phi[it] < -0.5  )
			   && !(t.track_eta[it] < 0.70 && t.track_eta[it] > 0.60 && t.track_phi[it] > -1.1 && t.track_phi[it] < -0.9  ));
      if(!extraCaloSel) continue;
      
      if(t.track_recoveto[it]!=0) continue;

      if(t.track_nLostOuterHits[it] < 2) continue;
      if(t.track_nLostInnerPixelHits[it] > 0) continue;
      if(!t.track_isHighPurity[it]) continue;

      // Isolation
      float niso = t.track_neuIso0p05[it];
      bool nisosel = niso<10.0;
      bool nisoselSTC = niso<10.0*isoSTC;
      
      float nreliso = t.track_neuRelIso0p05[it];
      bool nrelisosel = nreliso<0.1;
      bool nrelisoselSTC = nreliso<0.1*isoSTC;
      
      float iso = t.track_iso[it];
      bool isosel = iso<10.0;
      bool isoselSTC = iso<10.0*isoSTC;
      
      float reliso = t.track_reliso[it];
      bool relisosel = reliso<0.2;
      bool relisoselSTC = reliso<0.2*isoSTC;
      
      bool allisosel = nisosel && nrelisosel && isosel && relisosel;
      bool allisoselSTC = nisoselSTC && nrelisoselSTC && isoselSTC && relisoselSTC;
      if(!allisoselSTC) continue;

      // Track categorization
      bool isP = (t.track_nLayersWithMeasurement[it] == t.track_nPixelLayersWithMeasurement[it]) && (t.track_nLayersWithMeasurement[it]>=3);
      bool isL = t.track_nLayersWithMeasurement[it]>=7;
      bool isM = !isP && !isL;
      bool isShort = isP || isM || isL;
      
      if(!isShort) continue;
      
      if(t.track_nLayersWithMeasurement[it]>4 && t.track_nPixelLayersWithMeasurement[it]<2) continue;
      if(t.track_nLayersWithMeasurement[it]<=4 && t.track_nPixelLayersWithMeasurement[it]<3) continue;

      // Quality
      float pterr = t.track_ptErr[it];
      float pterrOPt2 = pterr/((t.track_pt[it])*(t.track_pt[it]));
      bool pteselL = pterrOPt2<0.2;
      bool pteselLSTC = pterrOPt2<0.2*qualSTC;
      bool pteselM = pterrOPt2<0.02;
      bool pteselMSTC = pterrOPt2<0.02*qualSTC;
      bool pteselT = pterrOPt2<0.005;
      bool pteselTSTC = pterrOPt2<0.005*qualSTC;
      
      bool ptesel = pteselT;
      bool pteselSTC = pteselTSTC;
      if(isP){
	ptesel = pteselL;
	pteselSTC = pteselLSTC;
      }
      if(isM){
	ptesel = pteselM;
	pteselSTC = pteselMSTC;
      }
      
      float dz = std::fabs(t.track_dz[it]);
      bool dzsel = dz<0.05;
      bool dzselSTC = dz<0.05*qualSTC;

      float dxy = std::fabs(t.track_dxy[it]);
      bool dxyselL = dxy<0.02;
      bool dxyselLSTC = dxy<0.02*qualSTC;
      bool dxyselT = dxy<0.01;
      bool dxyselTSTC = dxy<0.01*qualSTC;
      
      bool dxysel = dxyselT;
      bool dxyselSTC = dxyselTSTC;
      if(isP){
	dxysel = dxyselL;
	dxyselSTC = dxyselLSTC;
      }
      
      bool qualsel = ptesel && dxysel && dzsel;
      bool qualselSTC = pteselSTC && dxyselSTC && dzselSTC;
      
      if(!qualselSTC) continue;

      if(isL){
       
	float mt = TMath::Sqrt(2*(t.track_pt[it])*(t.met_pt)*(1-TMath::Cos(t.track_phi[it]-t.met_phi)));
	if (mt<100. && t.track_pt[it]<150.) continue;

      }

      bool isST = false;
      bool isSTC = false;

      if(allisosel && qualsel) isST=true;
      if(!isST) isSTC=true;

      if(isST){ 
	++nst;
	
	if(isP){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++nt_PstF23;
	    else if(t.nJet30>=4) ++nt_PstF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++nt_PstV23lHt;
	      else if(t.nJet30>=4) ++nt_PstV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++nt_PstV23hHt;
	      else if(t.nJet30>=4)++nt_PstV4hHt;
	    }
	  }
	}
	else if(isM){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++nt_MstF23;
	    else if(t.nJet30>=4) ++nt_MstF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++nt_MstV23lHt;
	      else if(t.nJet30>=4) ++nt_MstV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++nt_MstV23hHt;
	      else if(t.nJet30>=4)++nt_MstV4hHt;
	    }
	  }
	}
	else if(isL){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++nt_LstF23;
	    else if(t.nJet30>=4) ++nt_LstF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++nt_LstV23lHt;
	      else if(t.nJet30>=4) ++nt_LstV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++nt_LstV23hHt;
	      else if(t.nJet30>=4)++nt_LstV4hHt;
	    }
	  }
	}

	if(isP){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++ntpe_PstF23;
	    else if(t.nJet30>=4) ++ntpe_PstF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++ntpe_PstV23lHt;
	      else if(t.nJet30>=4) ++ntpe_PstV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++ntpe_PstV23hHt;
	      else if(t.nJet30>=4)++ntpe_PstV4hHt;
	    }
	  }
	}
	else if(isM){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++ntpe_MstF23;
	    else if(t.nJet30>=4) ++ntpe_MstF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++ntpe_MstV23lHt;
	      else if(t.nJet30>=4) ++ntpe_MstV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++ntpe_MstV23hHt;
	      else if(t.nJet30>=4)++ntpe_MstV4hHt;
	    }
	  }
	}
	else if(isL){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++ntpe_LstF23;
	    else if(t.nJet30>=4) ++ntpe_LstF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++ntpe_LstV23lHt;
	      else if(t.nJet30>=4) ++ntpe_LstV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++ntpe_LstV23hHt;
	      else if(t.nJet30>=4)++ntpe_LstV4hHt;
	    }
	  }
	}
	
      }
      if(isSTC){
	++nstc;
	
	if(isP){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++nt_PstcF23;
	    else if(t.nJet30>=4) ++nt_PstcF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++nt_PstcV23lHt;
	      else if(t.nJet30>=4) ++nt_PstcV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++nt_PstcV23hHt;
	      else if(t.nJet30>=4)++nt_PstcV4hHt;
	    }
	  }
	}
	else if(isM){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++nt_MstcF23;
	    else if(t.nJet30>=4) ++nt_MstcF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++nt_MstcV23lHt;
	      else if(t.nJet30>=4) ++nt_MstcV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++nt_MstcV23hHt;
	      else if(t.nJet30>=4)++nt_MstcV4hHt;
	    }
	  }
	}
	else if(isL){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++nt_LstcF23;
	    else if(t.nJet30>=4) ++nt_LstcF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++nt_LstcV23lHt;
	      else if(t.nJet30>=4) ++nt_LstcV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++nt_LstcV23hHt;
	      else if(t.nJet30>=4)++nt_LstcV4hHt;
	    }
	  }
	}

	if(isP){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++ntpe_PstcF23;
	    else if(t.nJet30>=4) ++ntpe_PstcF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++ntpe_PstcV23lHt;
	      else if(t.nJet30>=4) ++ntpe_PstcV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++ntpe_PstcV23hHt;
	      else if(t.nJet30>=4)++ntpe_PstcV4hHt;
	    }
	  }
	}
	else if(isM){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++ntpe_MstcF23;
	    else if(t.nJet30>=4) ++ntpe_MstcF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++ntpe_MstcV23lHt;
	      else if(t.nJet30>=4) ++ntpe_MstcV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++ntpe_MstcV23hHt;
	      else if(t.nJet30>=4)++ntpe_MstcV4hHt;
	    }
	  }
	}
	else if(isL){
	  if(t.mt2<100){
	    if(t.nJet30<4) ++ntpe_LstcF23;
	    else if(t.nJet30>=4) ++ntpe_LstcF4;
	  }
	  else if(t.mt2>100 && t.mt2<200){
	    if(t.ht<1000){
	      if(t.nJet30<4) ++ntpe_LstcV23lHt;
	      else if(t.nJet30>=4) ++ntpe_LstcV4lHt;
	    }
	    else if(t.ht>1000){
	      if(t.nJet30<4) ++ntpe_LstcV23hHt;
	      else if(t.nJet30>=4)++ntpe_LstcV4hHt;
	    }
	  }
	}

      }
    }

    
    if(nstc>0){
      if(t.mt2<100){
	if(t.nJet30<4) ++ne_stcF23;
	else if(t.nJet30>=4) ++ne_stcF4;
      }
      else if(t.mt2>100 && t.mt2<200){
	if(t.ht<1000){
	  if(t.nJet30<4) ++ne_stcV23lHt;
	  else if(t.nJet30>=4) ++ne_stcV4lHt;
	}
	else if(t.ht>1000){
	  if(t.nJet30<4) ++ne_stcV23hHt;
	  else if(t.nJet30>=4)++ne_stcV4hHt;
	}
      }
    }

    if(nst>0){
      if(t.mt2<100){
	if(t.nJet30<4) ++ne_stF23;
	else if(t.nJet30>=4) ++ne_stF4;
      }
      else if(t.mt2>100 && t.mt2<200){
	if(t.ht<1000){
	  if(t.nJet30<4) ++ne_stV23lHt;
	  else if(t.nJet30>=4) ++ne_stV4lHt;
	}
	else if(t.ht>1000){
	  if(t.nJet30<4) ++ne_stV23hHt;
	  else if(t.nJet30>=4)++ne_stV4hHt;
	}
      }
    }
   
    if(ntpe_PstF23>0) ++ne_PstF23;
    if(ntpe_PstcF23>0) ++ne_PstcF23;
    if(ntpe_PstF4>0) ++ne_PstF4;
    if(ntpe_PstcF4>0) ++ne_PstcF4;
    if(ntpe_PstV23lHt>0) ++ne_PstV23lHt;
    if(ntpe_PstcV23lHt>0) ++ne_PstcV23lHt;
    if(ntpe_PstV4lHt>0) ++ne_PstV4lHt;
    if(ntpe_PstcV4lHt>0) ++ne_PstcV4lHt;
    if(ntpe_PstV23hHt>0) ++ne_PstV23hHt;
    if(ntpe_PstcV23hHt>0) ++ne_PstcV23hHt;
    if(ntpe_PstV4hHt>0) ++ne_PstV4hHt;
    if(ntpe_PstcV4hHt>0) ++ne_PstcV4hHt;
    
    if(ntpe_MstF23>0) ++ne_MstF23;
    if(ntpe_MstcF23>0) ++ne_MstcF23;
    if(ntpe_MstF4>0) ++ne_MstF4;
    if(ntpe_MstcF4>0) ++ne_MstcF4;
    if(ntpe_MstV23lHt>0) ++ne_MstV23lHt;
    if(ntpe_MstcV23lHt>0) ++ne_MstcV23lHt;
    if(ntpe_MstV4lHt>0) ++ne_MstV4lHt;
    if(ntpe_MstcV4lHt>0) ++ne_MstcV4lHt;
    if(ntpe_MstV23hHt>0) ++ne_MstV23hHt;
    if(ntpe_MstcV23hHt>0) ++ne_MstcV23hHt;
    if(ntpe_MstV4hHt>0) ++ne_MstV4hHt;
    if(ntpe_MstcV4hHt>0) ++ne_MstcV4hHt;
    
    if(ntpe_LstF23>0) ++ne_LstF23;
    if(ntpe_LstcF23>0) ++ne_LstcF23;
    if(ntpe_LstF4>0) ++ne_LstF4;
    if(ntpe_LstcF4>0) ++ne_LstcF4;
    if(ntpe_LstV23lHt>0) ++ne_LstV23lHt;
    if(ntpe_LstcV23lHt>0) ++ne_LstcV23lHt;
    if(ntpe_LstV4lHt>0) ++ne_LstV4lHt;
    if(ntpe_LstcV4lHt>0) ++ne_LstcV4lHt;
    if(ntpe_LstV23hHt>0) ++ne_LstV23hHt;
    if(ntpe_LstcV23hHt>0) ++ne_LstcV23hHt;
    if(ntpe_LstV4hHt>0) ++ne_LstV4hHt;
    if(ntpe_LstcV4hHt>0) ++ne_LstcV4hHt;
    
  }
  
  std::cout<<std::endl;
  std::cout<<"Event counts:"<<std::endl;
  std::cout<<std::endl;
  std::cout << "FSR 2-3j (STC | ST): " << ne_stcF23 << " " << ne_stF23 << std::endl;
  std::cout << "isP: " << ne_PstcF23 << " " << ne_PstF23 << std::endl;
  std::cout << "isM: " << ne_MstcF23 << " " << ne_MstF23 << std::endl;
  std::cout << "isL: " << ne_LstcF23 << " " << ne_LstF23 << std::endl;
  
  std::cout<<std::endl;
  std::cout << "FSR  +4j (STC | ST): " << ne_stcF4  << " " << ne_stF4  << std::endl;
  std::cout << "isP: " << ne_PstcF4  << " " << ne_PstF4  << std::endl;
  std::cout << "isM: " << ne_MstcF4  << " " << ne_MstF4  << std::endl;
  std::cout << "isL: " << ne_LstcF4  << " " << ne_LstF4  << std::endl;
  
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<std::endl;
  std::cout << "VSR lHT 2-3j (STC | ST): " << ne_stcV23lHt << " " << ne_stV23lHt << std::endl;
  std::cout << "isP: " << ne_PstcV23lHt << " " << ne_PstV23lHt << std::endl;
  std::cout << "isM: " << ne_MstcV23lHt << " " << ne_MstV23lHt << std::endl;
  std::cout << "isL: " << ne_LstcV23lHt << " " << ne_LstV23lHt << std::endl;

  std::cout<<std::endl;
  std::cout << "VSR lHT  +4j (STC | ST): " << ne_stcV4lHt  << " " << ne_stV4lHt   << std::endl;
  std::cout << "isP: " << ne_PstcV4lHt  << " " << ne_PstV4lHt   << std::endl;
  std::cout << "isM: " << ne_MstcV4lHt  << " " << ne_MstV4lHt   << std::endl;
  std::cout << "isL: " << ne_LstcV4lHt  << " " << ne_LstV4lHt   << std::endl;
  
  std::cout<<std::endl;

  std::cout<<std::endl;
  std::cout << "VSR hHT 2-3j (STC | ST): " << ne_stcV23hHt << " " << ne_stV23hHt << std::endl;
  std::cout << "isP: " << ne_PstcV23hHt << " " << ne_PstV23hHt << std::endl;
  std::cout << "isM: " << ne_MstcV23hHt << " " << ne_MstV23hHt << std::endl;
  std::cout << "isL: " << ne_LstcV23hHt << " " << ne_LstV23hHt << std::endl;

  std::cout<<std::endl;
  std::cout << "VSR hHT  +4j (STC | ST): " << ne_stcV4hHt  << " " << ne_stV4hHt   << std::endl;
  std::cout << "isP: " << ne_PstcV4hHt  << " " << ne_PstV4hHt   << std::endl;
  std::cout << "isM: " << ne_MstcV4hHt  << " " << ne_MstV4hHt   << std::endl;
  std::cout << "isL: " << ne_LstcV4hHt  << " " << ne_LstV4hHt   << std::endl;
  

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<"Track counts:"<<std::endl;
  std::cout<<std::endl;
  std::cout << "FSR 2-3j (STC | ST): "<< std::endl;
  std::cout << "isP: " << nt_PstcF23 << " " << nt_PstF23 << std::endl;
  std::cout << "isM: " << nt_MstcF23 << " " << nt_MstF23 << std::endl;
  std::cout << "isL: " << nt_LstcF23 << " " << nt_LstF23 << std::endl;
  
  std::cout<<std::endl;
  std::cout << "FSR  +4j (STC | ST): " << std::endl;
  std::cout << "isP: " << nt_PstcF4  << " " << nt_PstF4  << std::endl;
  std::cout << "isM: " << nt_MstcF4  << " " << nt_MstF4  << std::endl;
  std::cout << "isL: " << nt_LstcF4  << " " << nt_LstF4  << std::endl;
  
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<std::endl;
  std::cout << "VSR lHT 2-3j (STC | ST): " << std::endl;
  std::cout << "isP: " << nt_PstcV23lHt << " " << nt_PstV23lHt << std::endl;
  std::cout << "isM: " << nt_MstcV23lHt << " " << nt_MstV23lHt << std::endl;
  std::cout << "isL: " << nt_LstcV23lHt << " " << nt_LstV23lHt << std::endl;

  std::cout<<std::endl;
  std::cout << "VSR lHT  +4j (STC | ST): " << std::endl;
  std::cout << "isP: " << nt_PstcV4lHt  << " " << nt_PstV4lHt   << std::endl;
  std::cout << "isM: " << nt_MstcV4lHt  << " " << nt_MstV4lHt   << std::endl;
  std::cout << "isL: " << nt_LstcV4lHt  << " " << nt_LstV4lHt   << std::endl;
  
  std::cout<<std::endl;

  std::cout<<std::endl;
  std::cout << "VSR hHT 2-3j (STC | ST): " << std::endl;
  std::cout << "isP: " << nt_PstcV23hHt << " " << nt_PstV23hHt << std::endl;
  std::cout << "isM: " << nt_MstcV23hHt << " " << nt_MstV23hHt << std::endl;
  std::cout << "isL: " << nt_LstcV23hHt << " " << nt_LstV23hHt << std::endl;

  std::cout<<std::endl;
  std::cout << "VSR hHT  +4j (STC | ST): " << std::endl;
  std::cout << "isP: " << nt_PstcV4hHt  << " " << nt_PstV4hHt   << std::endl;
  std::cout << "isM: " << nt_MstcV4hHt  << " " << nt_MstV4hHt   << std::endl;
  std::cout << "isL: " << nt_LstcV4hHt  << " " << nt_LstV4hHt   << std::endl;
  
}
