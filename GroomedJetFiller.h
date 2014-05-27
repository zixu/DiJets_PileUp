/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VplusJetSubstructure
 *
 *
 * Authors:
 *   Zijun Xu, PKU - zijun.xu@cern.ch
 *   Nhan V Tran, Fermilab - kalanand@fnal.gov
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   To fill groomed jet related quantities into a specified TTree
 *   Works with jets in PAT data format.
 * History:
 *   
 *
 * Copyright (C) 2012 FNAL 
 *****************************************************************************/

#ifndef GroomedJetFiller_h
#define GroomedJetFiller_h

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLorentzVector.h"


#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <iostream>
#include "fastjet/tools/Filter.hh"
// the following includes are only needed when combining filtering with subtraction
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"

#include <cstdio>   // needed for io

using namespace fastjet;
using namespace std;



// system include files
#include "tools.h"
#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"
#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h"


#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
//#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include <fastjet/PseudoJet.hh>

#include "../../fjcontrib-1.009/include/fastjet/contrib/GenericSubtractor.hh"
#include "../../fjcontrib-1.009/include/fastjet/contrib/JetCleanser.hh"
#include "../../fjcontrib-1.009/GenericSubtractor/ExampleShapes.hh"
using namespace fastjet::contrib;
#include "NSubjettinessRatio.hh"
#include "GroomTool.hh"



//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
//#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

// class decleration
namespace ewk
{

	// typedefs
	typedef std::auto_ptr<fastjet::ClusterSequence>        ClusterSequencePtr;
	typedef std::auto_ptr<fastjet::GhostedAreaSpec>        ActiveAreaSpecPtr;
	typedef std::auto_ptr<fastjet::RangeDefinition>        RangeDefPtr;


	class GroomedJetFiller {
		public:
			/// specify the name of the TTree, and the configuration for it
			GroomedJetFiller(const char *name, 
						TTree* tree, 
						const std::string jetAlgorithmLabel,
						const std::string additionalLabel,
						//const edm::ParameterSet& iConfig,
						bool isGen = 0);

			/// default constructor
			GroomedJetFiller() {};
			/// Destructor, does nothing 
			~GroomedJetFiller(){  
				//if(jec_) delete jec_; 
				//if(jecL1_) delete jecL1_; 
				//if(jecL2_) delete jecL2_; 
				//if(jecL3_) delete jecL3_; 
				//if(jecUnc_) delete jecUnc_; 
				if(mJetDef) delete mJetDef;
				if(mAreaDefinition) delete mAreaDefinition;
			};
			//~GroomedJetFiller(){ };

			/// To be called once per event to fill the values for groomed jets
			void fill(
						//	const edm::Event& iEvent, 
						std::vector<fastjet::PseudoJet> FJparticles, bool doJetCleansing=0, std::vector<fastjet::PseudoJet> FJparticles_hardcharge=std::vector<fastjet::PseudoJet>(), std::vector<fastjet::PseudoJet> FJparticles_pileupcharge=std::vector<fastjet::PseudoJet>(), std::vector<fastjet::PseudoJet> FJparticles_fullneutral=std::vector<fastjet::PseudoJet>() ); //full event, bool, hard charge, pile up charge, full neutral

			void Init();//init all the branch value

			// ----------member data ---------------------------
			static const int NUM_JET_MAX = 4;
			static const int NUM_JETCLEANSING_MODE_MAX = 200;

		protected:
			// 'mutable' because we will fill it from a 'const' method
			mutable std::vector<std::string> bnames;
			/// Helper function for main constructor
			void SetBranch( float* x, std::string name);
			void SetBranch( int* x, std::string name);
			void SetBranchSingle( double* x, std::string name);
			void SetBranchSingle( float* x, std::string name);
			void SetBranchSingle( int* x, std::string name);
			double getJEC(double curJetEta, double curJetPt, double curJetE, double curJetArea); 
			double getJECL1(double curJetEta, double curJetPt, double curJetE, double curJetArea); 
			double getJECL2(double curJetEta, double curJetPt, double curJetE, double curJetArea); 
			double getJECL3(double curJetEta, double curJetPt, double curJetE, double curJetArea); 
			TLorentzVector getCorrectedJet(fastjet::PseudoJet& jet, Int_t jeclevel=100, bool debug=0);
			fastjet::PseudoJet getScaledJet(fastjet::PseudoJet& jet, double scale);
			void computeCore( std::vector<fastjet::PseudoJet> constits, double Rval, float &m_core, float &pt_core );
			void computePlanarflow(std::vector<fastjet::PseudoJet> constits,double Rval,fastjet::PseudoJet jet,std::string mJetAlgo,float &planarflow);
			float computeJetCharge( std::vector<fastjet::PseudoJet> constits, std::vector<float> pdgIds, float Ejet );        
			float getPdgIdCharge( float fid );        
			Double_t getnPV( );//const edm::Event& iEvent );        
			Double_t getrho( );//const edm::Event& iEvent );        
			Double_t getrho_Hand(std::vector<fastjet::PseudoJet>  FJparticles);        
			Double_t getrho_Hand2(std::vector<fastjet::PseudoJet>  FJparticles, fastjet::Subtractor** subtractor);
			Double_t getrho_Grid(std::vector<fastjet::PseudoJet>  FJparticles, fastjet::Subtractor** subtractor);
			fastjet::PseudoJet do_rhoA_correction(fastjet::PseudoJet jet_origin, double rho, double area);
			void do_GenericShapeSubtract_correction(fastjet::PseudoJet jet_origin, fastjet::BackgroundEstimatorBase* bge_rho, float& jetpt_new, float& jetmass_new, float& tau2tau1_shapesubtraction, Double_t& rhom);
			void get_nsubjettiness(fastjet::PseudoJet jet_origin, float &tau1, float &tau2, float &tau3, float &tau4, float & tau2tau1);

			//Jet Cleansing
			JetCleanser makeJVFCleanser(fastjet::JetDefinition subjet_def, std::string projectmode="CMS", double fcut=-1.0, int nsj=-1 );//projectmode: CMS or ATLAS
			JetCleanser makeLinearCleanser(fastjet::JetDefinition subjet_def, double linear_para0,std::string projectmode="CMS", double fcut=-1, int nsj=-1 );//projectmode: CMS or ATLAS
			JetCleanser makeGausCleanser(fastjet::JetDefinition subjet_def, double gaus_para0, double gaus_para1, double gaus_para2, double gaus_para3, std::string projectmode="CMS", double fcut=-1, int nsj=-1 );//projectmode: CMS or ATLAS
			void DoJetCleansing(fastjet::JetDefinition jetDef, std::vector<fastjet::PseudoJet> FJparticles, std::vector<fastjet::PseudoJet> FJparticles_hardcharge, std::vector<fastjet::PseudoJet> FJparticles_pileupcharge, std::vector<fastjet::PseudoJet> FJparticles_fullneutral, fastjet::PseudoJet recoJet  );

			// ----------member data ---------------------------
			TTree* tree_;
			bool runningOverMC_;
			bool applyJECToGroomedJets_;

			//FactorizedJetCorrector* jec_;
			//FactorizedJetCorrector* jecL1_;
			//FactorizedJetCorrector* jecL2_;
			//FactorizedJetCorrector* jecL3_;
			//JetCorrectionUncertainty* jecUnc_;

			std::string jetAlgorithmLabel_;
			std::string jetAlgorithmAdditonalLabel_;
			std::string JetsFor_rho;
			std::string JEC_GlobalTag_forGroomedJet;
			//edm::InputTag mPrimaryVertex;
			//std::string mGroomedJet;

			// specific parameters
			std::string lableGen;        

			//Jet Def
			double                  mJetRadius;
			std::string             mJetAlgo;        
			fastjet::JetDefinition *mJetDef;
			//Jet Area
			fastjet::AreaDefinition *mAreaDefinition;
			//Jet background estimate
			std::auto_ptr<fastjet::JetMedianBackgroundEstimator>  mBgeMedi;
			std::auto_ptr<fastjet::GridMedianBackgroundEstimator> mBgeGrid;

			double mJetChargeKappa;
			bool   mDoQJets; 
			int    mQJetsPreclustering;
			int    mQJetsN;
			double mNsubjettinessKappa;
			bool   mSaveConstituents;

			std::vector<GroomTool> vect_groomtools;

		private:

			// ----------member data ---------------------------
			// names of modules, producing object collections    
			float jetmass_uncorr[NUM_JET_MAX];
			float jetmass_tr_uncorr[NUM_JET_MAX];
			float jetmass_tr1_uncorr[NUM_JET_MAX]; //added
			float jetmass_tr2_uncorr[NUM_JET_MAX]; //added
			float jetmass_tr3_uncorr[NUM_JET_MAX]; //added


			float jetmass_ft_uncorr[NUM_JET_MAX];
			float jetmass_pr_uncorr[NUM_JET_MAX];
			float jetmass_pr1_uncorr[NUM_JET_MAX];// added 
			float jetmass_pr2_uncorr[NUM_JET_MAX];// added
			float jetmass_pr3_uncorr[NUM_JET_MAX];// added
			float tau2tau1[NUM_JET_MAX];
			float tau2tau1_shapesubtraction[NUM_JET_MAX];
			float tau2tau1_trimmingshapesubtraction[NUM_JET_MAX];
			float tau2tau1_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float tau1[NUM_JET_MAX];
			float tau2[NUM_JET_MAX];
			float tau3[NUM_JET_MAX];
			float tau4[NUM_JET_MAX];
			float tau1_shapesubtraction[NUM_JET_MAX];
			float tau2_shapesubtraction[NUM_JET_MAX];
			float tau3_shapesubtraction[NUM_JET_MAX];
			float tau4_shapesubtraction[NUM_JET_MAX];
			float tau1_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float tau2_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float tau3_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float tau4_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float massdrop_pr_uncorr[NUM_JET_MAX];
			float massdrop_pr1_uncorr[NUM_JET_MAX];// added
			float massdrop_pr2_uncorr[NUM_JET_MAX]; // added
			float massdrop_pr3_uncorr[NUM_JET_MAX]; // added

			int   number_jet_central;
			float jetpt_uncorr[NUM_JET_MAX];
			float jetpt[NUM_JET_MAX];
			float jeteta[NUM_JET_MAX];
			float jetphi[NUM_JET_MAX];
			float jete[NUM_JET_MAX];

			float jetpt_JECL1[NUM_JET_MAX];
			float jeteta_JECL1[NUM_JET_MAX];
			float jetphi_JECL1[NUM_JET_MAX];
			float jete_JECL1[NUM_JET_MAX];

			float jetpt_L1_rhoSW[NUM_JET_MAX];//rho from kt6PF SW
			float jetpt_L1_rhoHand[NUM_JET_MAX];//rho from kt6PF Hand
			float jetpt_L1_rhoHand2[NUM_JET_MAX];//rho from kt6PF Hand2
			float jetpt_L1_rhoGrid[NUM_JET_MAX];//rho from kt6PF Grid

			float jetpt_rho4A[NUM_JET_MAX];
			float jetpt_shapesubtraction[NUM_JET_MAX];
			float jetpt_trimmingshapesubtraction[NUM_JET_MAX];
			float jetpt_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float jeteta_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float jetphi_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];

			float jetpt_tr_uncorr[NUM_JET_MAX];
			float jetpt_tr[NUM_JET_MAX];
			float jeteta_tr[NUM_JET_MAX];
			float jetphi_tr[NUM_JET_MAX];
			float jete_tr[NUM_JET_MAX];

			float jetpt_tr1_uncorr[NUM_JET_MAX];//added 
			float jetpt_tr1[NUM_JET_MAX]; // added
			float jeteta_tr1[NUM_JET_MAX];//added
			float jetphi_tr1[NUM_JET_MAX];//added
			float jete_tr1[NUM_JET_MAX]; //added

			float jetpt_tr2_uncorr[NUM_JET_MAX];//added 
			float jetpt_tr2[NUM_JET_MAX]; // added
			float jeteta_tr2[NUM_JET_MAX];//added
			float jetphi_tr2[NUM_JET_MAX];//added
			float jete_tr2[NUM_JET_MAX]; //added

			float jetpt_tr3_uncorr[NUM_JET_MAX];//added 
			float jetpt_tr3[NUM_JET_MAX]; // added
			float jeteta_tr3[NUM_JET_MAX];//added
			float jetphi_tr3[NUM_JET_MAX];//added
			float jete_tr3[NUM_JET_MAX]; //added





			float jetpt_ft_uncorr[NUM_JET_MAX];
			float jetpt_ft[NUM_JET_MAX];
			float jeteta_ft[NUM_JET_MAX];
			float jetphi_ft[NUM_JET_MAX];
			float jete_ft[NUM_JET_MAX];

			float jetpt_pr_uncorr[NUM_JET_MAX];
			float jetpt_pr[NUM_JET_MAX];
			float jeteta_pr[NUM_JET_MAX];
			float jetphi_pr[NUM_JET_MAX];
			float jete_pr[NUM_JET_MAX];

			float jetpt_pr1_uncorr[NUM_JET_MAX];// added
			float jetpt_pr1[NUM_JET_MAX];
			float jeteta_pr1[NUM_JET_MAX];
			float jetphi_pr1[NUM_JET_MAX];
			float jete_pr1[NUM_JET_MAX];

			float jetpt_pr2_uncorr[NUM_JET_MAX];// added
			float jetpt_pr2[NUM_JET_MAX];
			float jeteta_pr2[NUM_JET_MAX];
			float jetphi_pr2[NUM_JET_MAX];
			float jete_pr2[NUM_JET_MAX];

			float jetpt_pr3_uncorr[NUM_JET_MAX];// added
			float jetpt_pr3[NUM_JET_MAX];
			float jeteta_pr3[NUM_JET_MAX];
			float jetphi_pr3[NUM_JET_MAX];
			float jete_pr3[NUM_JET_MAX];




			float prsubjet1_px[NUM_JET_MAX];
			float prsubjet1_py[NUM_JET_MAX];
			float prsubjet1_pz[NUM_JET_MAX];
			float prsubjet1_e[NUM_JET_MAX];        
			float prsubjet2_px[NUM_JET_MAX];
			float prsubjet2_py[NUM_JET_MAX];
			float prsubjet2_pz[NUM_JET_MAX];
			float prsubjet2_e[NUM_JET_MAX];        


			float jetmass[NUM_JET_MAX];
			float jetmass_JECL1[NUM_JET_MAX];
			float jetmass_rhoArea[NUM_JET_MAX];
			float jetmass_rhoGArea[NUM_JET_MAX];
			float jetmass_rho4A[NUM_JET_MAX];
			float jetmass_rhoG4A[NUM_JET_MAX];
			float jetmass_rhom4Am[NUM_JET_MAX];
			float jetmass_shapesubtraction[NUM_JET_MAX];
			float jetmass_trimmingshapesubtraction[NUM_JET_MAX];
			float jetmass_JetCleansing_DiffMode[NUM_JETCLEANSING_MODE_MAX];
			float jetmass_tr[NUM_JET_MAX];//main
			float jetmass_tr1[NUM_JET_MAX]; // added
			float jetmass_tr2[NUM_JET_MAX]; // added
			float jetmass_tr3[NUM_JET_MAX]; // added

			float jetmass_ft[NUM_JET_MAX];
			float jetmass_pr[NUM_JET_MAX];
			float jetmass_pr1[NUM_JET_MAX];// added
			float jetmass_pr2[NUM_JET_MAX];// added
			float jetmass_pr3[NUM_JET_MAX];// added

			float jetarea[NUM_JET_MAX];
			float jetarea_tr[NUM_JET_MAX];
			float jetarea_tr1[NUM_JET_MAX];// added
			float jetarea_tr2[NUM_JET_MAX];// added
			float jetarea_tr3[NUM_JET_MAX];// added 

			float jetarea_ft[NUM_JET_MAX];
			float jetarea_pr[NUM_JET_MAX]; 
			float jetarea_pr1[NUM_JET_MAX];// added
			float jetarea_pr2[NUM_JET_MAX];// added
			float jetarea_pr3[NUM_JET_MAX];// added          
			float massdrop_pr[NUM_JET_MAX];
			float massdrop_pr1[NUM_JET_MAX];// added
			float massdrop_pr2[NUM_JET_MAX];//added
			float massdrop_pr3[NUM_JET_MAX];//added

			float jetconstituents[NUM_JET_MAX];   
			float jetcharge[NUM_JET_MAX];           

			float rcores[11][NUM_JET_MAX];
			float ptcores[11][NUM_JET_MAX];

			float planarflow[11][NUM_JET_MAX];

			float qjetmass[50];
			float qjetmassdrop[50];

			float constituents0_eta[100];
			float constituents0_phi[100];        
			float constituents0_e[100];        
			int nconstituents0;    

			float constituents0pr_eta[100];
			float constituents0pr_phi[100];        
			float constituents0pr_e[100];        
			int nconstituents0pr;   

			std::vector<int> neutrals;
			std::vector<int> positives;
			std::vector<int> negatives;        
			std::vector<float>  charge_handle_Gen;
			std::vector<float>  PF_id_handle;


			double rhoVal_;//rho kt6PF SW
			double rhoVal_hand;//rho kt6PF Hand
			double rhoVal_hand2;//rho kt6PF Hand2 
			double rhoVal_grid;//rho Grid
			double nPV_;
			bool isGenJ;
	};
}
#endif
