#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EventModel/Event.h"
//#include "TrigEvent/TrigEvent.h"
//#include "TrigEvent/TrigData.h"
#include "McTruth/McParticle.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"


#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector; 
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
	typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "KSKPiPiAlg/KSKPiPi.h"

#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/WTrackParameter.h"
#include "ParticleID/ParticleID.h"
#include <vector>
const double mpi = 0.13957;
const double mk = 0.493677;
const double  xmass[7] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272, 0.497614, 1.115683};//e, muon, pion, kaon, proton,  k_s, Lamda 
const double velc = 299.792458; //tof_path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,Ncut_ks,Ncut_svtxfit;

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

KSKPiPi::KSKPiPi(const std::string& name, ISvcLocator* pSvcLocator):
 Algorithm(name,pSvcLocator)
{	
	declareProperty("ChisqCut",        m_chisq_cut=200);
	declareProperty("KsVr0cut",	   m_ks_vr0cut =20);
	declareProperty("KsVz0cut",	   m_ks_vz0cut =25);
	declareProperty("KsCut",	   m_kshort_cut=0.030);
	declareProperty("Vr0cut",          m_vr0cut=1.0);
	declareProperty("Vz0cut",          m_vz0cut=5.0);
	declareProperty("EnergyThreshold", m_energyThreshold=0.04);
	declareProperty("GammaPhiCut",     m_gammaPhiCut=20.0);
	declareProperty("GammaThetaCut",   m_gammaThetaCut=20.0);
	declareProperty("GammaAngleCut",   m_gammaAngleCut=20.0);
	declareProperty("Test4C",          m_test4C = 1);
	declareProperty("Test5C", 	   m_test5C = 1);
	declareProperty("CheckDedx", 	   m_checkDedx = 1);
 	declareProperty("CheckTof",  	   m_checkTof = 1);
 	declareProperty("Dorecoil",  	   m_Dorecoil = 1);
 	declareProperty("McTruth",  	   m_McTruth = 1);
}

//*******************************************************
StatusCode KSKPiPi::initialize()
{
	MsgStream log(msgSvc(), name());

	log << MSG::INFO << "in initialize()" <<endmsg;

	StatusCode status;
	NTuplePtr nt1(ntupleSvc(),"FILE1/vxyz");
	if (nt1) m_tuple1 = nt1;
	else 
	{
		m_tuple1 = ntupleSvc()->book("FILE1/vxyz", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if (m_tuple1)
		{
			status = m_tuple1->addItem ("vx0",       m_vx0);
			status = m_tuple1->addItem ("vy0",       m_vy0);
			status = m_tuple1->addItem ("vz0",       m_vz0);
			status = m_tuple1->addItem ("vr0",       m_vr0);
			status = m_tuple1->addItem ("rvxy0",   m_rvxy0);
			status = m_tuple1->addItem ("rvx0",     m_rvz0);
			status = m_tuple1->addItem ("rvphio", m_rvphi0);
		}
		else 
		{	
			log << MSG::ERROR <<"              Cannot book N-tuple1:" << long (m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	

	NTuplePtr nt2(ntupleSvc(), "FILE1/photon");
	if ( nt2 ) m_tuple2 = nt2;
	else 
	{
		m_tuple2 = ntupleSvc()->book ("FILE1/photon", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if(m_tuple2)
		{
			status = m_tuple2->addItem("dthe", m_dthe);
			status = m_tuple2->addItem("dphi", m_dphi);
			status = m_tuple2->addItem("dang", m_dang);
			status = m_tuple2->addItem("eraw", m_eraw);
		}	
		else 
		{
			log << MSG::ERROR << "        Cannot book N-tuple2:" << long (m_tuple2) << endmsg;
			return StatusCode::FAILURE;
		}
	}


	NTuplePtr nt3(ntupleSvc(), "FILE1/etot");
	if ( nt3 ) m_tuple3 =nt3;
	else 
	{
		m_tuple3 = ntupleSvc()->book ("FILE1/etot", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if (m_tuple3)
		{
			status = m_tuple3->addItem ("m2gg", m_m2gg);
			status = m_tuple3->addItem ("etot", m_etot);
			status = m_tuple3->addItem ("momtot", m_mom_tot);
			status = m_tuple3->addItem ("mkp",   m_mkp);
			status = m_tuple3->addItem ("mpim", m_mpim);
			status = m_tuple3->addItem ("mks",   m_mks);
		}
		else
		{
			log << MSG::ERROR << "          Cannot book N-tuple3:" << long(m_tuple3) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	NTuplePtr nt_chi(ntupleSvc(), "FILE1/mchi2");
	if (nt_chi) m_tuple_chi = nt_chi;
	else
	{
		m_tuple_chi = ntupleSvc()->book("FILE1/mchi2", CLID_ColumnWiseTuple, "check 4c fit chi2");
		if (m_tuple_chi)
		{
			status = m_tuple_chi->addItem("chi2_4c", m_chi_4c);
		}
		else
		{
			log << MSG::ERROR << "    Cannot book N-Tuple_chi:" << long(m_tuple_chi) << endmsg;
			return StatusCode::FAILURE;

		}
		
	}
	if(m_McTruth == 1)
	{
		NTuplePtr nt_McTruth(ntupleSvc(), "FILE1/mctruth");
		if (nt_McTruth) m_tuple_mctruth = nt_McTruth;
		else
		{
			m_tuple_mctruth = ntupleSvc()->book("FILE1/mctruth",CLID_ColumnWiseTuple, "check mctruth");
			if(m_tuple_mctruth)
			{
				status = m_tuple_mctruth->addItem("pi0_mc_recoil_mass",		       m_pi0_mc_recoil);
				status = m_tuple_mctruth->addItem("pi0_mc_mass",		       m_pi0_mc_mass);
				status = m_tuple_mctruth->addItem("mctrth_pionp_bef",		       m_pionp_ks_bef);
				status = m_tuple_mctruth->addItem("mctrth_pionm_bef",		       m_pionm_bef);
				status = m_tuple_mctruth->addItem("mctrth_pionm_ks_bef",	       m_pionm_ks_bef);
				status = m_tuple_mctruth->addItem("mctrth_kshort_bef",		      m_kshort_bef);
				status = m_tuple_mctruth->addItem("mctrth_kaonp_bef",		       m_kaonp_bef);
				status = m_tuple_mctruth->addItem("mctrth_pion0_bef",   	       m_pion0_bef);
				status = m_tuple_mctruth->addItem("mctrth_pionp_costh_bef",	 m_pionp_ks_costh_bef);
				status = m_tuple_mctruth->addItem("mctrth_pionm_costh_bef",	 m_pionm_costh_bef);
				status = m_tuple_mctruth->addItem("mctrth_pionm_ks_costh_bef",	 m_pionm_ks_costh_bef);
				status = m_tuple_mctruth->addItem("mctrth_kshort_costh_bef",	m_kshort_costh_bef);
				status = m_tuple_mctruth->addItem("mctrth_kaonp_costh_bef",	 m_kaonp_costh_bef);
				status = m_tuple_mctruth->addItem("mctrth_pion0_costh_bef",	 m_pion0_costh_bef);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-Tuple_mctrth:" << long(m_tuple_mctruth) << endmsg;
			}
		}
	}
	
	NTuplePtr nt_vf(ntupleSvc(), "FILE1/vtxfit");
	if (nt_vf) m_tuple_vf = nt_vf;
	else
	{
		m_tuple_vf = ntupleSvc()->book("FILE1/vtxfit", CLID_ColumnWiseTuple, "after vertexfit");
		if ( m_tuple_vf )
		{
			status = m_tuple_vf->addItem ("mks_pip",		m_mpip_ks);
			status = m_tuple_vf->addItem ("mks_pim",		m_mpim_ks);
			status = m_tuple_vf->addItem ("mpim",			m_mpim_vf);
			status = m_tuple_vf->addItem ("mkp",			 m_mkp_vf);
			status = m_tuple_vf->addItem ("mks",			 m_mks_vf);		
			status = m_tuple_vf->addItem ("etot",			m_etot_vf);
		}
		else
		{
			log << MSG::ERROR << "	Cannot book	N-Tuple_vf:"	<< long(m_tuple_vf) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	if (m_Dorecoil == 1)
	{
		NTuplePtr nt_recoil(ntupleSvc(), "FILE1/recoil");
		if (nt_recoil) m_tuple_recoil = nt_recoil;
		else
		{
			m_tuple_recoil = ntupleSvc()->book("FILE1/recoil", CLID_ColumnWiseTuple, "for recoil");
			if (m_tuple_recoil)
			{
				status = m_tuple_recoil->addItem ("4momentum_index", 	m_4momentum_index_recoil, 0,10);
				status = m_tuple_recoil->addIndexedItem ("4momentum_matrix", 	m_4momentum_index_recoil, 4,m_4momentum_recoil);
			}	
			else
			{
				log <<MSG::ERROR << "Cannot book   N-Tuple_recoil" 	<< long(m_tuple_recoil) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	if (m_test4C==1)
	{
		NTuplePtr nt4(ntupleSvc(),"FILE1/fit4c");
		if (nt4) m_tuple4 = nt4;
		else
		{
			m_tuple4 = ntupleSvc()->book ("FILE1/fit4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple4)
			{	
				status = m_tuple4->addItem ("chi2",   m_chi1);
				status = m_tuple4->addItem ("mpi0",   m_mpi0_4c);
				status = m_tuple4->addItem ("mks",  m_mks_4c);
				status = m_tuple4->addItem ("etot",  m_etot_4c);
				status = m_tuple4->addItem ("egam1",  m_egam1_4c);
				status = m_tuple4->addItem ("egam2",  m_egam2_4c);
				status = m_tuple4->addItem ("4momentum_index",	m_4momentum_index_4c, 0, 10);
				status = m_tuple4->addIndexedItem("4momentum_matrix",	m_4momentum_index_4c, 4,	m_4momentum_4c);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-Tuple4:" << long(m_tuple4) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	if (m_test5C==1)
	{
		NTuplePtr nt5(ntupleSvc(),"FILE1/fit5c");
		if (nt5) m_tuple5 = nt5;
		else
		{
			m_tuple5 = ntupleSvc()->book ("FILE1/fit5c",CLID_ColumnWiseTuple, "ks N-Tuple example");
			if (m_tuple5)
			{
				status = m_tuple5->addItem ("chi2", m_chi2);
				status = m_tuple5->addItem ("mpi0", m_mpi0);
				status = m_tuple5->addItem ("mks", m_mks_5c);
				status = m_tuple5->addItem ("etot", m_etot_5c);
                                status = m_tuple5->addItem ("4momentum_index",  m_4momentum_index_5c, 0, 10);
                                status = m_tuple5->addIndexedItem("4momentum_matrix",  m_4momentum_index_5c, 4,   m_4momentum_5c);
			}
			else 
			{
				log << MSG::ERROR << "          Cannot book N-Tuple5:" << long(m_tuple5) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
		
/*		
		NTuplePtr nt6(ntupleSvc(),"FILE1/geff");
		if (nt6) m_tuple6 = nt6;
		else
		{
			m_tuple6 = ntupleSvc()->book ("FILE1/geff",CLID_ColumnWiseTuple, "ks N-Tuple example");
			if(m_tuple6)	
			{
				status = m_tuple6->addItem ("fcos", m_fcos);
				status = m_tuple6->addItem ("elow", m_elow);
			}
			else
			{
				log << MSG::ERROR << "    Cannot book N-Tuple6:" << long(m_tuple6) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	

	if(m_checkDedx == 1)
	{
		NTuplePtr nt7(ntupleSvc(),"FILE1/dedx");
		if ( nt7 ) m_tuple7 = nt7;
		else
		{
			m_tuple7 = ntupleSvc()->book ("FILE1/dedx", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if ( m_tuple7 )
			{
				status = m_tuple7->addItem ("ptrk",      m_ptrk);
				status = m_tuple7->addItem ("chie",      m_chie);
				status = m_tuple7->addItem ("chimu",    m_chimu);
				status = m_tuple7->addItem ("chipi",    m_chipi);
				status = m_tuple7->addItem ("chik",      m_chik);
				status = m_tuple7->addItem ("chip",      m_chip);
				status = m_tuple7->addItem ("proPH",   m_probPH);
				status = m_tuple7->addItem ("normPH",  m_normPH);
				status = m_tuple7->addItem ("ghit",      m_ghit);
				status = m_tuple7->addItem ("thit",      m_thit);
			}
			else
			{
				log << MSG::ERROR << "         Cannot book N-Tuple7:" << long(m_tuple7) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	

	if(m_checkTof ==1)
	{
		NTuplePtr nt8(ntupleSvc(),"FILE1/tofe");
		if ( nt8 ) m_tuple8 = nt8;
		else
		{
			m_tuple8 = ntupleSvc()->book ("FILE1/tofe", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if ( m_tuple8 )
			{
				status = m_tuple8->addItem ("ptrk",  m_ptot_etof);
				status = m_tuple8->addItem ("cntr",  m_cntr_etof);
				status = m_tuple8->addItem ("ph",      m_ph_etof);
				status = m_tuple8->addItem ("rhit",  m_rhit_etof);
				status = m_tuple8->addItem ("qual",  m_qual_etof);
				status = m_tuple8->addItem ("te",      m_te_etof);
				status = m_tuple8->addItem ("tmu",    m_tmu_etof);
				status = m_tuple8->addItem ("tpi",    m_tpi_etof);
				status = m_tuple8->addItem ("tk",      m_tk_etof);
				status = m_tuple8->addItem ("tp",      m_tp_etof);
			}
			else
			{
				log << MSG::ERROR << "           Cannot book N-Tuple8:" <<long( m_tuple8) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	if( m_checkTof == 1 )
	{
		NTuplePtr nt9(ntupleSvc(),"FILE1/tof1");
		if ( nt9 ) m_tuple9 = nt9 ;
		else
		{
			m_tuple9 = ntupleSvc()->book ("FILE1/tof1", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if ( m_tuple9 )
			{
				status = m_tuple9->addItem("ptrk",  m_ptot_btof1);
				status = m_tuple9->addItem("cntr",  m_cntr_btof1);
				status = m_tuple9->addItem("ph",      m_ph_btof1);
				status = m_tuple9->addItem("zhit",  m_zhit_btof1);
				status = m_tuple9->addItem("qual",  m_qual_btof1);
				status = m_tuple9->addItem("te",      m_te_btof1);
                                status = m_tuple9->addItem("tmu",    m_tmu_btof1);
                                status = m_tuple9->addItem("tpi",    m_tpi_btof1);
                                status = m_tuple9->addItem("tk",      m_tk_btof1);
                                status = m_tuple9->addItem("tp",      m_tp_btof1);
			}
			else  
			{
				log << MSG::ERROR <<"   Cannot book N-Tuple9:" << long(m_tuple9) << endmsg;
				return StatusCode::FAILURE;
			}
		}
	}
	

	if (m_checkTof == 1)
	{
		NTuplePtr nt10(ntupleSvc(),"FILE1/tof2");
		if (nt10)  m_tuple10 = nt10;
		else
		{
			m_tuple10 = ntupleSvc()->book ("FILE1/tof2",CLID_ColumnWiseTuple, "ks N-Tuple example");
			if ( m_tuple10 )
			{
				status = m_tuple10->addItem("ptrk", m_ptrk_btof2);
				status = m_tuple10->addItem("cntr", m_cntr_btof2);
				status = m_tuple10->addItem("ph",     m_ph_btof2);
				status = m_tuple10->addItem("zhit", m_zhit_btof2);
				status = m_tuple10->addItem("qual", m_qual_btof2);
				status = m_tuple10->addItem("te",     m_te_btof2);
				status = m_tuple10->addItem("tmu",   m_tmu_btof2);
				status = m_tuple10->addItem("tpi",   m_tpi_btof2);
				status = m_tuple10->addItem("tk",     m_tk_btof2);
				status = m_tuple10->addItem("tp",     m_tp_btof2);
			}
			else
			{
				log << MSG::ERROR << "      Cannot  book N-Tuple10:" << long(m_tuple10) << endmsg;
				return StatusCode::FAILURE;
			}
		}	
	}


	NTuplePtr nt11(ntupleSvc(),"FILE1/pid");
	if (nt11)  m_tuple11 = nt11;
	else
	{
		m_tuple11 = ntupleSvc()->book ("FILE1/pid",CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple11 )  
		{
			status = m_tuple11->addItem ("ptrk", m_ptrk_pid);
			status = m_tuple11->addItem ("cost", m_cost_pid);
			status = m_tuple11->addItem ("dedx", m_dedx_pid);
			status = m_tuple11->addItem ("tof1", m_tof1_pid);
			status = m_tuple11->addItem ("tof2", m_tof2_pid);
			status = m_tuple11->addItem ("prob", m_prob_pid);			
		}
		else
		{
			log << MSG::ERROR << " 		Cannot book N-Tuple11:" << long(m_tuple11) << endmsg;
			return StatusCode::FAILURE;
		}
	}
*/
//----------------------end of book----------------------------
	log << MSG::INFO << "successfully return from initialize()" << endmsg;
	return StatusCode::SUCCESS;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode KSKPiPi::execute()
{
	//std::cout << "execute()" << std::endl;
	
	MsgStream log(msgSvc(),name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	//log << MSG::DEBUG <<"run, evtnum = "
	  //  << runNo << " , "
	    //<< event << endreq;
	//cout<<"event "<<event<<endl;
	//cout <<"  " << endl;
	Ncut0++;
	


/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MC TRUTH   ////////////////////////////////////////////////////////
	if(m_McTruth ==2)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),                                                    
				"/Event/MC/McParticleCol");                                    
		if(runNo<0)
		{                                                                                   
			if(!mcParticleCol)
			{	                                                                                           
				//log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;                                          
				return StatusCode::FAILURE;                                                                                 
			}                                                                                                             
			else	
			{                                                                                                         
				bool PsiDecay(false);                                                                                      
				Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();                                            
				HepLorentzVector pY4260, ppi0_mc, pks_mc, pkp_mc, ppim_mc, ppim_ks_mc, ppip_ks_mc; 
				HepLorentzVector ecms( 0.04, 0, 0, 4.26 );
//				cout <<"ecms"<<ecms.boost(-ecms.boostVector())<<endl;
				int i=0;
				for (; iter_mc != mcParticleCol->end(); iter_mc++)
				{      ;                                                   
					if ((*iter_mc)->primaryParticle() ) continue;                                                             
					if (!(*iter_mc)->decayFromGenerator() ) continue;                                                         
					//if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;                                                  
					if ((*iter_mc)->particleProperty()==9030443) 
					{
						pY4260 = (*iter_mc)->initialFourMomentum();
						PsiDecay = true; 
					}                                            
					if (!PsiDecay) continue;                                                                                 
					//if(!(*iter_mc)->leafParticle()) continue;                                                               
					if((*iter_mc)->particleProperty() == 211)
					{                                                                
						if ( ((*iter_mc)->mother()).particleProperty() == 310 )
						{	++i;
							ppip_ks_mc = (*iter_mc)->initialFourMomentum();
							m_pionp_ks_bef = (*iter_mc)->initialFourMomentum().vect().mag();
							m_pionp_ks_costh_bef = (*iter_mc)->initialFourMomentum().vect().cosTheta();   
						}                                    
						else
						{
						cout << "                         strange pi+"	<< endl;
						}
					}                                                                                                         
					if((*iter_mc)->particleProperty() == -211) 	//pi-
					{                                          
						if ( ((*iter_mc)->mother()).particleProperty() == 310 )
						{	++i;
							ppim_ks_mc   = (*iter_mc)->initialFourMomentum();
							m_pionm_ks_bef   = (*iter_mc)->initialFourMomentum().vect().mag();
							m_pionm_ks_costh_bef = (*iter_mc)->initialFourMomentum().vect().cosTheta();        
						}
						if ( ((*iter_mc)->mother()).particleProperty() == 9030443 )
						{	++i;
							ppim_mc   = (*iter_mc)->initialFourMomentum();                                     
							m_pionm_bef   = (*iter_mc)->initialFourMomentum().vect().mag();
							m_pionm_costh_bef = (*iter_mc)->initialFourMomentum().vect().cosTheta();        
						}
						else
						{
							cout << "                       strange pi-"	<< endl;
						}

					}
					if((*iter_mc)->particleProperty() == 310)	//kshort
					{	++i;
						pks_mc    = (*iter_mc)->initialFourMomentum();
						m_kshort_bef    = (*iter_mc)->initialFourMomentum().vect().mag();
						m_kshort_costh_bef = (*iter_mc)->initialFourMomentum().vect().cosTheta();
					}                                     
					if((*iter_mc)->particleProperty() == 321)     	//k+                       
					{	++i;
						pkp_mc    = (*iter_mc)->initialFourMomentum();
						m_kaonp_bef    = (*iter_mc)->initialFourMomentum().vect().mag();
						m_kaonp_costh_bef = (*iter_mc)->initialFourMomentum().vect().cosTheta();
					}
					if((*iter_mc)->particleProperty() == 111)    	//pi0                            
					{	++i;
						ppi0_mc    = (*iter_mc)->initialFourMomentum();
						m_pion0_bef    = (*iter_mc)->initialFourMomentum().vect().mag();
						m_pion0_costh_bef = (*iter_mc)->initialFourMomentum().vect().cosTheta();
					}                           
					else
					{
						cout << "                       strange particle "  << endl;
					}
				}       
				
//				HepLorentzVector a1,a2,a3,a4,a5,a6,a7;
				HepLorentzVector a1(ecms - pks_mc           - pkp_mc - ppim_mc );
				HepLorentzVector a2(ecms - pks_mc - ppi0_mc          - ppim_mc );
				HepLorentzVector a3(ecms - pks_mc - ppi0_mc - pkp_mc           );
				HepLorentzVector a4(ecms         - ppi0_mc - pkp_mc - ppim_mc );
				  
cout <<"===============================================" << i << endl;
				           
cout <<"++++++++++++++++++++++++++++++++++++++" <<endl;
cout <<"  ppi0_mc =    " <<ppi0_mc    << endl;//cout for test
cout <<"  pks_pip_mc = " <<ppip_ks_mc << endl;
cout <<"  pks_pim_mc = " <<ppim_ks_mc << endl;
cout <<"  pkp_mc   =   " <<pkp_mc     << endl;
cout <<"  ppim_mc  =   " <<ppim_mc    << endl;
cout <<"  pks_mc   =   " <<pks_mc     << endl;
cout <<"++++++++++++++++++++++++++++++++++++++" <<endl;
cout << " a2 = " << a2 << endl;
cout <<"  ptot =    " << ( pks_mc + ppi0_mc + pkp_mc + ppim_mc -ecms)  << endl;//cout for test
cout <<" a1 - a1 = " << 
   a1 + a2 + a3 + a4 
 - (ecms - pks_mc           - pkp_mc - ppim_mc )
 - (ecms - pks_mc - ppi0_mc          - ppim_mc ) 
 - (ecms - pks_mc - ppi0_mc - pkp_mc           )
 - (ecms         - ppi0_mc - pkp_mc - ppim_mc )
 << endl;
cout <<" a1 - a1 = " << 
   a1 + a2 + a3 + a4 -  4 * ecms 
  + pks_mc + pkp_mc + ppim_mc 
  + pks_mc + ppi0_mc + ppim_mc 
  + pks_mc + ppi0_mc + pkp_mc           
  + ppi0_mc + pkp_mc + ppim_mc 
 << endl;
cout <<"  a1+a2+a3+a4 = " << a1 + a2 + a3 + a4 - ecms + 3 * (pks_mc + ppi0_mc + pkp_mc + ppim_mc -ecms)<< endl;
cout <<"  ppi0_recoil =    " << a1  << "mpi0_recoil"<< endl;//cout for test
cout <<"  pkp_recoil   =   " << a2  << "mkp_recoil" << endl;
cout <<"  ppim_recoil  =   " << a3  << "mpim_recoil"<<endl;
cout <<"  pks_recoil   =   " << a4  << "mks_recoil" <<endl;

cout << " a1 - ppi0_mc" << a1 - ppi0_mc << endl;
cout <<"  mpi0_mc =    " <<ppi0_mc.m()    << endl;
cout <<"  mpi0_recoil =    " << a1.m()  << endl;

				m_pi0_mc_recoil=(ecms - ppim_mc - pkp_mc - pks_mc).m();
				m_pi0_mc_mass=ppi0_mc.m2();
			}                                                                                                             
		}                                                 
		m_tuple_mctruth->write();
	}





        SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
        log << MSG::DEBUG <<"ncharge, nneu, tottks = "
            << evtRecEvent->totalCharged() << " , "
            << evtRecEvent->totalNeutral() << " , "
            << evtRecEvent->totalTracks()  << endreq;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////





	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	if(evtRecEvent->totalCharged() !=4 || evtRecTrkCol->size()<6 ||evtRecTrkCol->size()>20)
	return StatusCode::SUCCESS;
	
	int m_trk_index[4]={-1, -1, -1, -1};
	Vint m_ptrk_list, m_ntrk_list;
	m_ptrk_list.clear();
	m_ntrk_list.clear();
	EvtRecTrackIterator m_itTrk_begin = evtRecTrkCol->begin();
	for( int i = 0; i < evtRecEvent->totalCharged(); i++ )
	{
		EvtRecTrackIterator itTrk = m_itTrk_begin + i ;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;
		RecMdcKalTrack *mdcTrk = (*itTrk)->mdcKalTrack();
		if(fabs(mdcTrk->z()) > m_ks_vz0cut) continue;
		if(fabs(mdcTrk->r()) > m_ks_vr0cut) continue;
		if( mdcTrk->charge() >0) m_ptrk_list.push_back(i);
		else m_ntrk_list.push_back(i);
	}
	



	if(m_ptrk_list.size() < 1 || m_ntrk_list.size() < 1) 
		return StatusCode::SUCCESS;


	HepLorentzVector m_lv_lab(0.04, 0, 0, 4.260);
	HepLorentzVector pkshort;
	WTrackParameter wvks_pipTrk, wvks_pimTrk, vtxfit0_track;
	VertexParameter	vtxfit0_vpar;
	int temp_ptrk_index(-1), temp_ntrk_index(-1);
	double temp_chisq = 9999.0;
	for( int i = 0; i<m_ptrk_list.size(); i++)
	{
		for( int j = 0; j < m_ntrk_list.size(); j++)
		{
			EvtRecTrackIterator itTrk = m_itTrk_begin+m_ptrk_list[i];
			RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();
			HepLorentzVector m_lv_ptrk(mdcTrk->p4(xmass[2]));
			itTrk = m_itTrk_begin + m_ntrk_list[j];
			mdcTrk = (*itTrk)->mdcKalTrack();
			HepLorentzVector m_lv_ntrk(mdcTrk->p4(xmass[2]));
			if(fabs((m_lv_ptrk + m_lv_ntrk).m()-xmass[5]) > m_kshort_cut ) continue;
			// do second vertex fit for k_s
			HepPoint3D vx(0., 0., 0. );
			HepSymMatrix Evx(3, 0);
			double bx = 1E+6;
			double by = 1E+6;
			double bz = 1E+6;
			Evx[0][0] = bx*bx;
			Evx[1][1] = by*by;
			Evx[2][2] = bz*bz;
	
			VertexParameter vxpar;
			vxpar.setVx(vx);
			vxpar.setEvx(Evx);
			
			VertexFit *vtxfit0 = VertexFit::instance();
			SecondVertexFit *vtxfit = SecondVertexFit::instance();
	
			itTrk = m_itTrk_begin + m_ptrk_list[i];
			RecMdcKalTrack *mdcTrk_1 = (*itTrk)->mdcKalTrack();
			mdcTrk_1->setPidType(RecMdcKalTrack::pion);
			itTrk = m_itTrk_begin + m_ntrk_list[j];
			RecMdcKalTrack *mdcTrk_2 = (*itTrk)->mdcKalTrack();
			mdcTrk_2->setPidType(RecMdcKalTrack::pion);
			WTrackParameter wpip(xmass[2], mdcTrk_1->getZHelix(), mdcTrk_1->getZError());
			WTrackParameter wpim(xmass[2], mdcTrk_2->getZHelix(), mdcTrk_2->getZError());
			
			vtxfit0->init();
			vtxfit0->AddTrack(0, wpip);
			vtxfit0->AddTrack(1, wpim);
			vtxfit0->AddVertex(0, vxpar, 0, 1);
			if(!(vtxfit0->Fit(0))) continue;
			vtxfit0->BuildVirtualParticle(0);
	
			vtxfit->init();
			vtxfit->AddTrack(0, vtxfit0->wVirtualTrack(0));
			vtxfit->setVpar(vtxfit0->vpar(0));
			if(!(vtxfit->Fit())) continue;
			
			HepLorentzVector m_lv_kshort(vtxfit->p4par());
			if(fabs(m_lv_kshort.m()-xmass[5]) > m_kshort_cut) continue;
			if( vtxfit->chisq() > m_chisq_cut) continue;
			if( vtxfit->chisq() < temp_chisq )
			{
				temp_chisq = vtxfit->chisq();
				temp_ptrk_index = i;
				temp_ntrk_index = j;
				//m_declen = vtxfit->decayLength();
				//m_decerr = vtxfit->decayLengthError();
				//m_four_mom
			pkshort = m_lv_kshort;
			wvks_pipTrk = vtxfit0->wtrk(0);
			wvks_pimTrk = vtxfit0->wtrk(1);
			vtxfit0_track = vtxfit0->wVirtualTrack(0);
			vtxfit0_vpar = vtxfit0->vpar(0);
			}
			
		}		
	}
	if(temp_ptrk_index < 0 || temp_ntrk_index < 0)   //why not < 0?
		return StatusCode::SUCCESS;
	Ncut_ks++;



/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////// check x0, y0, z0,r0,
        //suggest cut: |z0|<5 && r0<1
	Vint iGood, ikp, ikm, ipip, ipim;
	iGood.clear();
	ikp.clear();
	ikm.clear();
	ipip.clear();
	ipim.clear();
	Vp4 pkp, pkm, ppip, ppim;
	pkp.clear();
	pkm.clear();
	ppip.clear();
	ppim.clear();
	
	int nCharge = 0;

	Hep3Vector xorigin(0,0,0);

	IVertexDbSvc* vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if (vtxsvc->isVertexValid())
	{
		double *dbv = vtxsvc->PrimaryVertex();
		double *vv = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	
	for (int i =0; i<evtRecEvent->totalCharged();i++)
	{
	
		if( i == m_ptrk_list[temp_ptrk_index] ) continue;
		if( i == m_ntrk_list[temp_ntrk_index] ) continue;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		double pch=mdcTrk->p();
		double x0=mdcTrk->x();
		double y0=mdcTrk->y();
		double z0=mdcTrk->z();
		double phi0=mdcTrk->helix(1);
		double xv=xorigin.x();
		double yv=xorigin.y();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		m_vx0 = x0;
		m_vy0 = y0;
		m_vz0 = z0;
		m_vr0 = Rxy;

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);    //the initial point for MDC reconstruction
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double Rvxy0=fabs(vecipa[0]);   //the nearest distance to IP in xy plane
		double Rvz0=vecipa[3];          //the nearest distance to IP in z direction
		double Rvphi0=vecipa[1];
		m_rvxy0=Rvxy0;
		m_rvz0=Rvz0;
		m_rvphi0=Rvphi0;
	
		m_tuple1->write();
		
		if(fabs(Rvz0) > 10.0 || fabs(Rvxy0) > 1.0 )	continue;
		iGood.push_back(i);
		nCharge += mdcTrk->charge();
	}

	//Finish GOOD Charged Track Selection

	int nGood = iGood.size();
	//log << MSG::DEBUG << "ngood1, totcharge = " << nGood << " , " <<nCharge << endmsg;
	if((nGood !=2) || (nCharge!=0))
	{
		return StatusCode::SUCCESS;
	}
	Ncut1++;
	

			


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////



	
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nGood; i++)   
	{		
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		
		pid->init();
		pid->setMethod(pid->methodProbability() );
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE() );//use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon() );  //seperater  Pion/Kaon
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
//		if(pid->probPion() < 0.001 && pid->probKaon() < 0.001) 
//			return StatusCode::SUCCESS;
	

		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
		if(pid->probPion() > pid->probKaon())
		{	
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron ,muon, pion, kaon and proton; The default setting is pion
		
			if(mdcKalTrk->charge() > 0 )
			{
				ipip.push_back(iGood[i]);
				HepLorentzVector ptrk;
				ptrk.setPx(mdcKalTrk->px() );
				ptrk.setPy(mdcKalTrk->py() );
				ptrk.setPz(mdcKalTrk->pz() );
				double p3 = ptrk.mag();
				ptrk.setE(sqrt(p3*p3+mpi*mpi) );
			
				ppip.push_back(ptrk);
			}
			else
			{
				ipim.push_back(iGood[i]);
				HepLorentzVector ptrk;
				ptrk.setPx(mdcKalTrk->px() );
				ptrk.setPy(mdcKalTrk->py() );
				ptrk.setPz(mdcKalTrk->pz() );
				double p3 = ptrk.mag();
				ptrk.setE(sqrt(p3*p3+mpi*mpi));
				ppim.push_back(ptrk);
			}
		}
		else
		{
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::kaon);//PID can set to electron ,muon, pion, kaon and proton; The default setting is pion
			if(mdcKalTrk->charge() > 0 )
			{
				ikp.push_back(iGood[i]);
				HepLorentzVector ptrk;
				ptrk.setPx(mdcKalTrk->px() );
				ptrk.setPy(mdcKalTrk->py() );
				ptrk.setPz(mdcKalTrk->pz() );
				double p3 = ptrk.mag();
				ptrk.setE(sqrt(p3*p3+mk*mk) );
				pkp.push_back(ptrk);
			}
			else
			{
				ikm.push_back(iGood[i]);
				HepLorentzVector ptrk;
				ptrk.setPx(mdcKalTrk->px() );
				ptrk.setPy(mdcKalTrk->py() );
				ptrk.setPz(mdcKalTrk->pz() );
				double p3 = ptrk.mag();
				ptrk.setE(sqrt(p3*p3+mk*mk) );
				pkm.push_back(ptrk);
			}
		}
	}
	int npip = ppip.size();
	int npim = ppim.size();
	int nkp  = pkp.size();
	int nkm  = pkm.size();
	if( npim*nkp != 1) return SUCCESS;
	Ncut2++;
	

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


	Vint iGam;
        iGam.clear();
        for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++)
        {
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		//find the nearest charged track
		double dthe = 200.;
		double dphi = 200.;
		double dang = 200.;
		for(int j = 0; j < evtRecEvent->totalCharged(); j++)
		{
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();

			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);
			thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			if(angd < dang)
			{
				dang = angd;
				dthe = thed;
				dphi = phid;
			}
		}
		if(dang >= 200) continue;
		double eraw = emcTrk->energy();
		dthe = dthe * 180 / (CLHEP::pi);
		dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);
		m_dthe = dthe;
		m_dphi = dphi;
		m_dang = dang;
		m_eraw = eraw;
		m_tuple2->write();
		if(eraw < m_energyThreshold) continue;
		if(fabs(dang) < m_gammaAngleCut) continue;
		iGam.push_back(i);
        }
	// finish Good Photon Selection
        int nGam = iGam.size();

        //log << MSG::DEBUG << "num Good Photon " << nGam << " , " << evtRecEvent->totalNeutral() << endreq;
        if(nGam < 2 )
        {
                return StatusCode::SUCCESS;
        }
	Ncut3++;



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
	
  	Vp4 pGam;
 	pGam.clear();
 	for(int i = 0; i < nGam; i++) 
	{
    		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
    		RecEmcShower* emcTrk = (*itTrk)->emcShower();
    		double eraw = emcTrk->energy();
    		double phi = emcTrk->phi();
    		double the = emcTrk->theta();
    		HepLorentzVector ptrk;
    		ptrk.setPx(eraw*sin(the)*cos(phi));
    		ptrk.setPy(eraw*sin(the)*sin(phi));
    		ptrk.setPz(eraw*cos(the));
    		ptrk.setE(eraw);


   		pGam.push_back(ptrk);
  	}

	
	
	//Loop each gamma pair, check ppi0 and pTot
	HepLorentzVector pTot;
	for(int i = 0; i < nGam - 1;i++)
	{
		for(int j = i+1; j < nGam; j++)
		{
			HepLorentzVector p2g = pGam[i] + pGam[j];
			pTot   = pkp[0] + ppim[0]+pkshort;
			pTot  += p2g;
			m_mkp  = pkp[0].m();
			m_mpim = ppim[0].m();
			m_mks  = pkshort.m();
			m_m2gg = p2g.m();
			m_mom_tot = pTot.vect().mag();
			m_tuple3-> write();  
		}
		
	}




///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
	//RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+ipip[0]))->mdcKalTrack();
//**************vertexfit********************//////
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();
	RecMdcKalTrack *kpTrk = (*(evtRecTrkCol->begin() + ikp[0]))->mdcKalTrack();
	RecMdcKalTrack *ks_pipTrk = (*(evtRecTrkCol->begin() + m_ptrk_list[temp_ptrk_index]))->mdcKalTrack();
	RecMdcKalTrack *ks_pimTrk = (*(evtRecTrkCol->begin() + m_ntrk_list[temp_ntrk_index]))->mdcKalTrack();	

	//pimTrk->setPidType(RecMdcKalTrack::pion);	
	//kpTrk->setPidType(RecMdcKalTrack::kaon);
	//ks_pipTrk->setPidType(RecMdcKalTrack::pion);
	//ks_pimTrk->setPidType(RecMdckalTrack::pion);

	WTrackParameter wvkpTrk, wvpimTrk;// wvks_pipTrk, wvks_pimTrk;
	wvkpTrk = WTrackParameter (mk, kpTrk->getZHelixK(), kpTrk->getZErrorK() );
	wvpimTrk = WTrackParameter (mpi, pimTrk->getZHelix(), pimTrk->getZError() );
	//wvks_pipTrk = WTrackParameter (mpi, ks_pipTrk->getZHelix(), pimTrk->getZError() );
	//wvks_pimTrk = WTrackParameter (mpi, ks_pimTrk->getZHelix(), pimTrk->getZError() );
	
			
	HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3.0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx*bx;
	Evx[1][1] = by*by;
	Evx[2][2] = bz*bz;
	
	VertexParameter vxpar,vxpar1;
	vxpar.setVx(vx);
	vxpar.setEvx(Evx);
	vxpar1=vxpar;
/*	
	VertexFit* vtxfit1 = VertexFit::instance();
	vtxfit1->init();
	vtxfit1->AddTrack(0, wvks_pipTrk);
	vtxfit1->AddTrack(1, wvks_pimTrk);
	vtxfit1->AddVertex(0,vxpar,0,1);
	if(!vtxfit1->Fit(0)) return SUCCESS;
	if(!vtxfit1->Fit()) return SUCCESS;

	vtxfit1->BuildVirtualParticle(0);
	vtxfit1->Swim(0);
	
*/
	VertexFit *vtxfit2 = VertexFit::instance();
	vtxfit2->init();
	vtxfit2->AddTrack(0, wvkpTrk);
	vtxfit2->AddTrack(1, wvpimTrk);
	vtxfit2->AddVertex(0, vxpar1, 0, 1);
	if(!vtxfit2->Fit(0)) return SUCCESS;
	if(!vtxfit2->Fit()) return SUCCESS;
	vtxfit2->Swim(0);
	


	SecondVertexFit *svtxfit = SecondVertexFit::instance();
	svtxfit->init();
	svtxfit->AddTrack(0, vtxfit0_track);  	
	//svtxfit->AddTrack(0,vtxfit1->wVirtualTrack(0));
	svtxfit->setVpar(vtxfit0_vpar); 	
	//svtxfit->setVpar(vtxfit1->vpar(0));
	svtxfit->setPrimaryVertex(vtxfit2->vpar(0));
	if (!svtxfit->Fit())	return SUCCESS;
	Ncut_svtxfit++;
///////////// finished  vertexfit

	WTrackParameter    	wks_pip 	= wvks_pipTrk;		
	//WTrackParameter    	wks_pip 	= vtxfit1->wtrk(0);
	HepLorentzVector 	lv_ks_pip 	= wks_pip.p();
	//lv_ks_pip.boost(-0.094, 0, 0);
	
	WTrackParameter    	wks_pim 	= wvks_pimTrk;	
	//WTrackParameter    	wks_pim 	= vtxfit1->wtrk(1);
	HepLorentzVector 	lv_ks_pim 	= wks_pim.p();
	//lv_ks_pim.boost(-0.094, 0, 0);

	WTrackParameter		wks		= vtxfit0_track;
	//WTrackParameter	wks		= vtxfit1->wVirtualTrack(0);
	HepLorentzVector	lv_ks		= wks.p();
	//lv_ks.boost(-0.094, 0, 0);
	
	WTrackParameter       	wpim 		= vtxfit2->wtrk(1);
	HepLorentzVector    	lv_pim 		= wpim.p();
	//lv_pim.boost(-0.094, 0, 0);
	
	WTrackParameter        	wkp 		= vtxfit2->wtrk(0);
	HepLorentzVector     	lv_kp 		= wkp.p();	
	//lv_kp.boost(-0.094, 0, 0);
/*

	if(m_Dorecoil == 1)
	{	
		lv_ks_pip.boost(-0.094, 0, 0);
		lv_ks_pim.boost(-0.094, 0, 0);
		//lv_ks.boost(-0.094, 0, 0);
		lv_pim.boost(-0.094, 0, 0);
		lv_kp.boost(-0.094, 0, 0);
	
		HepLorentzVector lv_gam;
		for( int i=0; i<nGam; i++)
		{
			lv_gam +=pGam[i];
		}

		HepLorentzVector ecms(0., 0., 0., 4.260);
		m_4momentum_index_recoil=0;
		
cout <<"=============="<< endl;		
		m_4momentum_recoil[0][0] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pim).e(); //for ks_pi+`s recoil_fourmomentum
cout << m_4momentum_recoil[0][0] <<endl;
		m_4momentum_recoil[0][1] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pim).px(); 
cout << m_4momentum_recoil[0][1] <<endl;
		m_4momentum_recoil[0][2] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pim).py(); 
cout << m_4momentum_recoil[0][2] <<endl;
		m_4momentum_recoil[0][3] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pim).pz(); 
cout << m_4momentum_recoil[0][3] <<endl;
		++m_4momentum_index_recoil;
cout <<"=============="<< endl;		
		m_4momentum_recoil[1][0] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pip).e(); //for ks_pi-`s recoil_fourmomentum
		m_4momentum_recoil[1][1] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pip).px(); 
		m_4momentum_recoil[1][2] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pip).py(); 
		m_4momentum_recoil[1][3] = (ecms - lv_gam - lv_kp - lv_pim - lv_ks_pip).pz(); 
		++m_4momentum_index_recoil;

		m_4momentum_recoil[2][0] = (ecms - lv_gam - lv_pim - lv_ks_pip - lv_ks_pim).e(); //for k+`s recoil_fourmomentum
		m_4momentum_recoil[2][1] = (ecms - lv_gam - lv_pim - lv_ks_pip - lv_ks_pim).px(); 
		m_4momentum_recoil[2][2] = (ecms - lv_gam - lv_pim - lv_ks_pip - lv_ks_pim).py(); 
		m_4momentum_recoil[2][3] = (ecms - lv_gam - lv_pim - lv_ks_pip - lv_ks_pim).pz(); 
		++m_4momentum_index_recoil;
		
		m_4momentum_recoil[3][0] = (ecms - lv_gam - lv_kp - lv_ks_pip - lv_ks_pim).e(); //for pi-`s recoil_fourmomentum
		m_4momentum_recoil[3][1] = (ecms - lv_gam - lv_kp - lv_ks_pip - lv_ks_pim).px(); 
		m_4momentum_recoil[3][2] = (ecms - lv_gam - lv_kp - lv_ks_pip - lv_ks_pim).py(); 
		m_4momentum_recoil[3][3] = (ecms - lv_gam - lv_kp - lv_ks_pip - lv_ks_pim).pz(); 
		++m_4momentum_index_recoil;
		
		m_4momentum_recoil[4][0] = (ecms - lv_kp - lv_pim - lv_ks_pip - lv_ks_pim).e(); //for pi0`s recoil fourmomentum
		m_4momentum_recoil[4][1] = (ecms - lv_kp - lv_pim - lv_ks_pip - lv_ks_pim).px(); 
		m_4momentum_recoil[4][2] = (ecms - lv_kp - lv_pim - lv_ks_pip - lv_ks_pim).py(); 
		m_4momentum_recoil[4][3] = (ecms - lv_kp - lv_pim - lv_ks_pip - lv_ks_pim).pz(); 
		++m_4momentum_index_recoil;
		
		m_4momentum_recoil[5][0] = (ecms - lv_gam - lv_kp - lv_pim).e();             //for kshort`s recoil_fourmomentum
		m_4momentum_recoil[5][1] = (ecms - lv_gam - lv_kp - lv_pim).px();             
		m_4momentum_recoil[5][2] = (ecms - lv_gam - lv_kp - lv_pim).py();             
		m_4momentum_recoil[5][3] = (ecms - lv_gam - lv_kp - lv_pim).pz();             
		++m_4momentum_index_recoil;
		
		m_tuple_recoil->write();
	}*/
/*
	m_mpip_ks = lv_ks_pip.vect().mag();
	m_mpim_ks = lv_ks_pim.vect().mag();
	m_mpim_vf = lv_pim.vect().mag();
	m_mkp_vf  = lv_kp.vect().mag();
	
	//m_mks_vf  = lv_ks.vect().mag();
	m_etot_vf  = (lv_ks_pip+lv_ks_pim+lv_pim+lv_kp).e()+m_m2gg;
	
	m_tuple_vf->write();
*/


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	
	//Apply Kinematic 4C fit
	
	if(m_test4C==1)
	{
		HepLorentzVector ecms(0.0, 0, 0, 4.260);
		
		double chisq = 9999.;
		int ig1 = -1;
		int ig2 = -1;
		for(int i =0; i < nGam-1; i++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
			for(int j = i+1; j < nGam; j++)
			{
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
				kmfit->init();
				kmfit->AddTrack(0, wks_pip);
				kmfit->AddTrack(1, wks_pim);
				//kmfit->AddTrack(0, wks);
				kmfit->AddTrack(2, wkp);
				kmfit->AddTrack(3, wpim);
				kmfit->AddTrack(4, 0.0, g1Trk);
				kmfit->AddTrack(5, 0.0, g2Trk);
				kmfit->AddFourMomentum(0, ecms);
				if(!kmfit->Fit(0)) continue;
				bool oksq = kmfit->Fit();
				if(oksq)
				{
					double chi2 = kmfit->chisq();
					if(chi2 < chisq)
					{
						chisq = chi2;
						ig1 = iGam[i];
						ig2 = iGam[j];
					}
				}
			}
		}
		//m_chi_4c = kmfit->chisq();
		//m_tuple_chi->write();
		
		if(chisq < 200)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
			kmfit->init();
			kmfit->AddTrack(0, wks_pip);
			kmfit->AddTrack(1, wks_pim);
			kmfit->AddTrack(2, wkp);
			kmfit->AddTrack(3, wpim);
			kmfit->AddTrack(4, 0.0, g1Trk);
			kmfit->AddTrack(5, 0.0, g2Trk);
			kmfit->AddFourMomentum(0, ecms);
			bool oksq = kmfit->Fit();
			if(oksq)
			{
				HepLorentzVector ppi0 = kmfit->pfit(4) + kmfit->pfit(5);
				HepLorentzVector pks_pip = kmfit->pfit(0);
				HepLorentzVector pks_pim = kmfit->pfit(1);
				HepLorentzVector pkp = kmfit->pfit(2);
				HepLorentzVector ppim = kmfit->pfit(3);
				HepLorentzVector pgam1 = kmfit->pfit(4);
				HepLorentzVector pgam2 = kmfit->pfit(5);
				HepLorentzVector pks  = kmfit->pfit(0) + kmfit->pfit(1);
				HepLorentzVector ptot = kmfit->pfit(0) + kmfit->pfit(1)+ kmfit->pfit(2)+
							 kmfit->pfit(3)+ kmfit->pfit(4)+ kmfit->pfit(5);
				m_mpi0_4c    = ppi0.m();
				m_egam1_4c   = pgam1.e();
				m_egam2_4c   = pgam2.e();
				m_mks_4c  = pks.m();
				m_etot_4c = ptot.e();
				m_chi1 = kmfit->chisq();
			
				for ( int k = 0; k < 6; k++ )
				{
					kmfit->pfit(k).boost(-ecms.boostVector());
				}
				m_4momentum_index_4c=0;
				for ( int k = 0; k < 6; k++ )
				{
					m_4momentum_4c[k][0] = kmfit->pfit(k).e();
					m_4momentum_4c[k][1] = kmfit->pfit(k).px();
					m_4momentum_4c[k][2] = kmfit->pfit(k).py();
					m_4momentum_4c[k][3] = kmfit->pfit(k).pz();
				++m_4momentum_index_4c;
				}
				m_4momentum_4c[6][0] = (kmfit->pfit(4)+kmfit->pfit(5)).e();
				m_4momentum_4c[6][1] = (kmfit->pfit(4)+kmfit->pfit(5)).px();
				m_4momentum_4c[6][2] = (kmfit->pfit(4)+kmfit->pfit(5)).py();
				m_4momentum_4c[6][3] = (kmfit->pfit(4)+kmfit->pfit(5)).pz();
				++m_4momentum_index_4c;
				
				m_4momentum_4c[7][0] = (kmfit->pfit(0)+kmfit->pfit(1)).e();
				m_4momentum_4c[7][1] = (kmfit->pfit(0)+kmfit->pfit(1)).px();
				m_4momentum_4c[7][2] = (kmfit->pfit(0)+kmfit->pfit(1)).py();
				m_4momentum_4c[7][3] = (kmfit->pfit(0)+kmfit->pfit(1)).pz();
				++m_4momentum_index_4c;


				m_tuple4->write();
				Ncut4++;
			}
		}
		else
		{
			return StatusCode::SUCCESS;
		}
	}







	if(m_Dorecoil == 1)
	{
		HepLorentzVector ecms(0.00, 0, 0, 4.260);
		lv_ks_pip.boost(-ecms.boostVector());
		lv_ks_pim.boost(-ecms.boostVector());
		//lv_ks.boost(-ecms.boostVector());
		lv_pim.boost(-ecms.boostVector());
		lv_kp.boost(-ecms.boostVector());

		HepLorentzVector lv_gam;
		for( int i=0; i<nGam; i++)
		{
			lv_gam +=pGam[i];
		}

		m_4momentum_index_recoil=0;
		for(int i=0; i<6;  i++)
		{
			m_4momentum_recoil[i][0] = (ecms - kmfit->pfit(i)).e(); //for ks_pi+ ks_pi- k+ pi- gamma1 gamma2 recoil_fourmomentum
			m_4momentum_recoil[i][1] = (ecms - kmfit->pfit(i)).px();
 			m_4momentum_recoil[i][2] = (ecms - kmfit->pfit(i)).py();
			m_4momentum_recoil[i][3] = (ecms - kmfit->pfit(i)).pz();
			++m_4momentum_index_recoil;
		}
		m_4momentum_recoil[6][0] = (ecms - kmfit->pfit(0) - kmfit->pfit(1)).e(); //for ks recoil_fourmomentum
		m_4momentum_recoil[6][1] = (ecms - kmfit->pfit(0) - kmfit->pfit(1)).px();
		m_4momentum_recoil[6][2] = (ecms - kmfit->pfit(0) - kmfit->pfit(1)).py();
		m_4momentum_recoil[6][3] = (ecms - kmfit->pfit(0) - kmfit->pfit(1)).pz();
		++m_4momentum_index_recoil;

		m_4momentum_recoil[7][0] = (ecms - kmfit->pfit(4) - kmfit->pfit(5)).e(); //for pi0 recoil_fourmomentum
		m_4momentum_recoil[7][1] = (ecms - kmfit->pfit(4) - kmfit->pfit(5)).px();
		m_4momentum_recoil[7][2] = (ecms - kmfit->pfit(4) - kmfit->pfit(5)).py();
		m_4momentum_recoil[7][3] = (ecms - kmfit->pfit(4) - kmfit->pfit(5)).pz();
		++m_4momentum_index_recoil;
		m_tuple_recoil->write();
	}








	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	//Apply Kinematic 5C Fit

	//find the best combination over all possible pi+ pi- gamma gamma pair
	if(m_test5C==1)
	{
		HepLorentzVector ecms(0.034, 0, 0, 4.260);
		double chisq = 9999.;
		int ig1 = -1;
		int ig2 = -1;
		for(int i = 0; i < nGam-1; i++)
		{
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
			for(int j = i+1; j < nGam; j++)
			{
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
				kmfit->init();
				kmfit->AddTrack(0, wks_pip);
				kmfit->AddTrack(1, wks_pim);
				kmfit->AddTrack(2, wkp);
				kmfit->AddTrack(3, wpim);
				kmfit->AddTrack(4, 0.0, g1Trk);
				kmfit->AddTrack(5, 0.0, g2Trk);
				kmfit->AddResonance(0, 0.135, 4, 5);
				kmfit->AddFourMomentum(1, ecms);
				if(!kmfit->Fit(0)) continue;
				if(!kmfit->Fit(1)) continue;
				bool oksq = kmfit->Fit();
				if(oksq)
				{
					double chi2 = kmfit->chisq();
					if(chi2 < chisq)
					{
						chisq = chi2;
						ig1 = iGam[i];
						ig2 = iGam[j];
					}
				}
			}
		}
		//log << MSG::INFO << " chisq = " << chisq << endreq;

		if(chisq <200)
		{
			RecEmcShower* g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
			RecEmcShower* g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
			
			kmfit->init();
			kmfit->AddTrack(0, wks_pip);
			kmfit->AddTrack(1, wks_pim);
			kmfit->AddTrack(2, wkp);
			kmfit->AddTrack(3, wpim);			
			kmfit->AddTrack(4, 0.0, g1Trk);
			kmfit->AddTrack(5, 0.0, g2Trk);
			kmfit->AddResonance(0, 0.135, 4, 5);
			kmfit->AddFourMomentum(1, ecms);
			bool oksq = kmfit->Fit();
			if(oksq)
			{
				HepLorentzVector pks_pip  = kmfit->pfit(0) ;
				HepLorentzVector pks_pim  = kmfit->pfit(1) ;
				HepLorentzVector pkp  = kmfit->pfit(2) ;
				HepLorentzVector ppim  = kmfit->pfit(3) ;
				HepLorentzVector pgam1  = kmfit->pfit(4);
				HepLorentzVector pgam2  = kmfit->pfit(5) ;
				HepLorentzVector ppi0  = kmfit->pfit(4) + kmfit->pfit(5);
				HepLorentzVector pks = kmfit->pfit(0) + kmfit->pfit(1);
				HepLorentzVector ptot = kmfit->pfit(0) + kmfit->pfit(1)+ kmfit->pfit(2)+
							kmfit->pfit(3) + kmfit->pfit(4)+ kmfit->pfit(5);
	
				m_chi2    = kmfit->chisq();
				m_mks_5c  = pks.m();
				m_mpi0    = ppi0.m();
				m_etot_5c = ptot.e();
				//m_mrhm = prhom.m();
				//double eg1 = (kimfit->pfit(2)).e();
				//double eg2 = (kimfit->pfit(3)).e();
				//double fcos = abs(eg1-eg2)/ppi0.rho();	//meaning?

				for ( int k = 0; k < 6; k++ )
                                {
                                        kmfit->pfit(k).boost(-ecms.boostVector());
                                }
				m_4momentum_index_5c=0;
                                for ( int k = 0; k < 6; k++ )
                                {
                                        m_4momentum_5c[k][0] = kmfit->pfit(k).e();
                                        m_4momentum_5c[k][1] = kmfit->pfit(k).px();
                                        m_4momentum_5c[k][2] = kmfit->pfit(k).py();
                                        m_4momentum_5c[k][3] = kmfit->pfit(k).pz();
				++m_4momentum_index_5c;
                                }
                                m_4momentum_5c[6][0] = (kmfit->pfit(4)+kmfit->pfit(5)).e();
                                m_4momentum_5c[6][1] = (kmfit->pfit(4)+kmfit->pfit(5)).px();
                                m_4momentum_5c[6][2] = (kmfit->pfit(4)+kmfit->pfit(5)).py();
                                m_4momentum_5c[6][3] = (kmfit->pfit(4)+kmfit->pfit(5)).pz();
				++m_4momentum_index_5c;

                                m_4momentum_5c[7][0] = (kmfit->pfit(0)+kmfit->pfit(1)).e();
                                m_4momentum_5c[7][1] = (kmfit->pfit(0)+kmfit->pfit(1)).px();
                                m_4momentum_5c[7][2] = (kmfit->pfit(0)+kmfit->pfit(1)).py();
                                m_4momentum_5c[7][3] = (kmfit->pfit(0)+kmfit->pfit(1)).pz();      
				++m_4momentum_index_5c;



				m_tuple5->write();
				Ncut5++;
	
				// Measure the photon detection efficiences via  J/psi -> rho0 pi0
				
				/*if(fabs(prho0.m()-0.770) < 0.150)
				{
					if(fabs(fcos)<0.99)
					{
						m_fcos = (eg1-eg2)/ppi0.rho();
						m_elow = (eg1 < eg2) ? eg1 :eg2;
						m_tuple6->write();
						Ncut6++;
					}
				}*/
			}
		}
	}
	return StatusCode::SUCCESS;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode KSKPiPi::finalize()
{
	cout << "total number:			" << Ncut0 << endl;
   	//cout << "Pass k_short_cut:		" << Ncut_ks << endl;	
	//cout << "nGood == 2, nCharge == 0:	" << Ncut1 << endl;
	//cout << "Pass Pid:			" << Ncut2 << endl;
	//cout << "nGam >= 2:			" << Ncut3 << endl;
	//cout << "Pass svtxfit			" << Ncut_svtxfit << endl;
	cout << "Pass 4C:			" << Ncut4 << endl;
	cout << "Pass 5C: 			" << Ncut5 << endl;
	//cout << "J/psi->rho0 pi0:		" << Ncut6 << endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << " in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
