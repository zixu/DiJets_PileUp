/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VplusJetSubstructure
 *
 *
 * Authors:
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


// user include files
#include "GroomedJetFiller.h" 

ewk::GroomedJetFiller::GroomedJetFiller(const char *name, 
			TTree* tree, 
			const std::string jetAlgorithmLabel,
			const std::string additionalLabel,
			//const edm::ParameterSet& iConfig,
			bool isGen)
{
	tree_     = tree;

	//Gen Jet
	isGenJ = isGen;
	if(isGen){ lableGen = "Gen"; } else{ lableGen = ""; }

	// get algo and radius
	jetAlgorithmLabel_  = jetAlgorithmLabel;
	jetAlgorithmAdditonalLabel_  = additionalLabel;
	unsigned int labelSize = jetAlgorithmLabel.size();
	mJetAlgo = "";
	mJetAlgo.push_back( jetAlgorithmLabel.at(0) ); mJetAlgo.push_back( jetAlgorithmLabel.at(1) );
	if (labelSize == 3 || labelSize == 4 ){
		const char* tmp1 = &jetAlgorithmLabel.at(2);
		mJetRadius = atof( tmp1 )/10.;
	} else{
		std::cout << "problem in defining jet type!" << std::endl;
	}
	mJetDef= new fastjet::JetDefinition(fastjet::cambridge_algorithm, mJetRadius);
	if (mJetAlgo == "AK") mJetDef->set_jet_algorithm( fastjet::antikt_algorithm );
	else if (mJetAlgo == "CA") mJetDef->set_jet_algorithm( fastjet::cambridge_algorithm );
	else if (mJetAlgo == "KT") mJetDef->set_jet_algorithm( fastjet::kt_algorithm );
	//else throw cms::Exception("GroomedJetFiller") << " unknown jet algorithm " << std::endl;
	else cout<<" In GroomedJetFiller " << " unknown jet algorithm " << std::endl;
	cout<<"mJetRadius="<<mJetRadius<<endl;

	//Jet Area
	double ghostEtaMax = 4.4;
	int activeAreaRepeats = 1;
	double ghostArea = 0.01;
	fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
	fjActiveArea.set_fj2_placement(true);
	mAreaDefinition =new fastjet::AreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );

	// Declare all the branches of the tree
	SetBranch( jetpt_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_uncorr");
	SetBranchSingle( &number_jet_central, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_number_jet_central");
	SetBranch( jetmass_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_uncorr");
	SetBranch( jetmass_tr_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr_uncorr");
	SetBranch( jetmass_tr1_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr1_uncorr"); ///Added
	SetBranch( jetmass_tr2_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr2_uncorr"); // added
	SetBranch( jetmass_tr3_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr3_uncorr"); // added

	SetBranch( jetmass_ft_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_ft_uncorr");

	SetBranch( jetmass_pr_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr_uncorr");
	SetBranch( jetmass_pr1_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr1_uncorr");// added
	SetBranch( jetmass_pr2_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr2_uncorr");// added
	SetBranch( jetmass_pr3_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr3_uncorr");// added



	SetBranch( tau2tau1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2tau1");
	SetBranch( tau2tau1_shapesubtraction, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2tau1_shapesubtraction");
	SetBranch( tau2tau1_trimmingshapesubtraction, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2tau1_trimmingshapesubtraction");
	SetBranch( tau1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau1");
	SetBranch( tau2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2");
	SetBranch( tau3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau3");
	SetBranch( tau4, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau4");
	SetBranch( massdrop_pr_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_massdrop_pr_uncorr");

	SetBranch( jetpt, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt");
	SetBranch( jeteta, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta");
	SetBranch( jetphi, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi");
	SetBranch( jete, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e");

	SetBranch( jetpt_JECL1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_JECL1");
	SetBranch( jeteta_JECL1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_JECL1");
	SetBranch( jetphi_JECL1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_JECL1");
	SetBranch( jete_JECL1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_JECL1");

	SetBranch( jetpt_L1_rhoSW, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_L1_rhoSW");
	SetBranch( jetpt_L1_rhoHand, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_L1_rhoHand");
	SetBranch( jetpt_L1_rhoHand2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_L1_rhoHand2");
	SetBranch( jetpt_L1_rhoGrid, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_L1_rhoGrid");

	SetBranch( jetpt_rho4A, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_rho4A");
	SetBranch( jetpt_shapesubtraction, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_shapesubtraction");
	SetBranch( jetpt_trimmingshapesubtraction, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_trimmingshapesubtraction");

	// trimming
	SetBranch( jetpt_tr_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr_uncorr");
	SetBranch( jetpt_tr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr");
	SetBranch( jeteta_tr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_tr");
	SetBranch( jetphi_tr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_tr");
	SetBranch( jete_tr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_tr");

	SetBranch( jetpt_tr1_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr1_uncorr"); // Added
	SetBranch( jetpt_tr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr1");
	SetBranch( jeteta_tr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_tr1");
	SetBranch( jetphi_tr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_tr1");
	SetBranch( jete_tr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_tr1");

	SetBranch( jetpt_tr2_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr2_uncorr"); // Added
	SetBranch( jetpt_tr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr2");
	SetBranch( jeteta_tr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_tr2");
	SetBranch( jetphi_tr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_tr2");
	SetBranch( jete_tr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_tr2");

	SetBranch( jetpt_tr3_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr3_uncorr"); // Added
	SetBranch( jetpt_tr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_tr3");
	SetBranch( jeteta_tr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_tr3");
	SetBranch( jetphi_tr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_tr3");
	SetBranch( jete_tr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_tr3");

	//filtering
	SetBranch( jetpt_ft_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_ft_uncorr");
	SetBranch( jetpt_ft, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_ft");
	SetBranch( jeteta_ft, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_ft");
	SetBranch( jetphi_ft, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_ft");
	SetBranch( jete_ft, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_ft");

	//pruning
	SetBranch( jetpt_pr_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr_uncorr");
	SetBranch( jetpt_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr");
	SetBranch( jeteta_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_pr");
	SetBranch( jetphi_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_pr");
	SetBranch( jete_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_pr");

	SetBranch( jetpt_pr1_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr1_uncorr");// added
	SetBranch( jetpt_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr1");
	SetBranch( jeteta_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_pr1");
	SetBranch( jetphi_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_pr1");
	SetBranch( jete_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_pr1");

	SetBranch( jetpt_pr2_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr2_uncorr");// added
	SetBranch( jetpt_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr2");
	SetBranch( jeteta_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_pr2");
	SetBranch( jetphi_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_pr2");
	SetBranch( jete_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_pr2");

	SetBranch( jetpt_pr3_uncorr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr3_uncorr");// added
	SetBranch( jetpt_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_pr3");
	SetBranch( jeteta_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_pr3");
	SetBranch( jetphi_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_pr3");
	SetBranch( jete_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_e_pr3");



	SetBranch( prsubjet1_px, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet1_px");
	SetBranch( prsubjet1_py, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet1_py");
	SetBranch( prsubjet1_pz, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet1_pz");
	SetBranch( prsubjet1_e, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet1_e");
	SetBranch( prsubjet2_px, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet2_px");
	SetBranch( prsubjet2_py, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet2_py");
	SetBranch( prsubjet2_pz, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet2_pz");
	SetBranch( prsubjet2_e, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_prsubjet2_e");

	SetBranch( jetmass, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass");
	SetBranch( jetmass_JECL1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_JECL1");
	SetBranch( jetmass_rhoArea, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_rhoArea");
	SetBranch( jetmass_rhoGArea, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_rhoGArea");
	SetBranch( jetmass_rho4A, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_rho4A");
	SetBranch( jetmass_rhoG4A, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_rhoG4A");
	SetBranch( jetmass_rhom4Am, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_rhom4Am");
	SetBranch( jetmass_shapesubtraction, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_shapesubtraction");
	SetBranch( jetmass_trimmingshapesubtraction, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_trimmingshapesubtraction");

	SetBranch( jetmass_tr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr");
	SetBranch( jetmass_tr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr1");// added
	SetBranch( jetmass_tr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr2");// added
	SetBranch( jetmass_tr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_tr3");// added

	SetBranch( jetmass_ft, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_ft");
	SetBranch( jetmass_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr");
	SetBranch( jetmass_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr1");// added
	SetBranch( jetmass_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr2");// added
	SetBranch( jetmass_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_pr3");// added

	SetBranch( massdrop_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_massdrop_pr");
	SetBranch( massdrop_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_massdrop_pr1");// added
	SetBranch( massdrop_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_massdrop_pr2");// added
	SetBranch( massdrop_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_massdrop_pr3");// added

	SetBranch( jetarea, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area");

	SetBranch( jetarea_tr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_tr");
	SetBranch( jetarea_tr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_tr1");// added
	SetBranch( jetarea_tr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_tr2");// added
	SetBranch( jetarea_tr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_tr3");// added  

	SetBranch( jetarea_ft, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_ft");
	SetBranch( jetarea_pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_pr");
	SetBranch( jetarea_pr1, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_pr1");// added
	SetBranch( jetarea_pr2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_pr2");//added
	SetBranch( jetarea_pr3, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_area_pr3");// added

	SetBranch( jetconstituents, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_jetconstituents");
	SetBranch( jetcharge, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_jetcharge");

	// cores
	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rcores").c_str(), rcores, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rcores"+"[11][6]/F").c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rcores").c_str() );
	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_ptcores").c_str(), ptcores, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_ptcores"+"[11][6]/F").c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_ptcores").c_str() );

	//planarflow
	tree_->Branch((lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_planarflow").c_str(),planarflow, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_planarflow"+"[11][6]/F").c_str());
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_planarflow").c_str() );

	// qjets
	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_qjetmass").c_str(), qjetmass, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_qjetmass"+"[50]/F").c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_qjetmass").c_str() );
	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_qjetmassdrop").c_str(), qjetmassdrop, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_qjetmassdrop"+"[50]/F").c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_qjetmassdrop").c_str() );

	//Jet mass, pt after JetCleansing
	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_JetCleansing_DiffMode").c_str(), jetmass_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_mass_JetCleansing_DiffMode").c_str() );
	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_JetCleansing_DiffMode").c_str(), jetpt_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_pt_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_JetCleansing_DiffMode").c_str(), jeteta_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_eta_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_JetCleansing_DiffMode").c_str(), jetphi_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_phi_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2tau1_JetCleansing_DiffMode").c_str(), tau2tau1_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2tau1_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2tau1_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau1_JetCleansing_DiffMode").c_str(), tau1_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau1_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau1_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2_JetCleansing_DiffMode").c_str(), tau2_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau2_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau3_JetCleansing_DiffMode").c_str(), tau3_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau3_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau3_JetCleansing_DiffMode").c_str() );

	tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau4_JetCleansing_DiffMode").c_str(), tau4_JetCleansing_DiffMode, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau4_JetCleansing_DiffMode"+Form("[%i]/F",NUM_JETCLEANSING_MODE_MAX)).c_str() );
	bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_tau4_JetCleansing_DiffMode").c_str() );

	//if( iConfig.existsAs<bool>("GroomedJet_saveConstituents") ) 
	//  mSaveConstituents=iConfig.getParameter< bool >("GroomedJet_saveConstituents");
	//else mSaveConstituents = true;
	mSaveConstituents = true;

	if (mSaveConstituents){
		tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_eta").c_str(), constituents0_eta, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_eta"+"[100]/F").c_str() );
		bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_eta").c_str() );   
		tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_phi").c_str(), constituents0_phi, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_phi"+"[100]/F").c_str() );
		bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_phi").c_str() );   
		tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_e").c_str(), constituents0_e, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_e"+"[100]/F").c_str() );
		bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0_e").c_str() );   
		SetBranchSingle( &nconstituents0, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_nconstituents0" );

		tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_eta").c_str(), constituents0pr_eta, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_eta"+"[100]/F").c_str() );
		bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_eta").c_str() );   
		tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_phi").c_str(), constituents0pr_phi, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_phi"+"[100]/F").c_str() );
		bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_phi").c_str() );   
		tree_->Branch( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_e").c_str(), constituents0pr_e, (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_e"+"[100]/F").c_str() );
		bnames.push_back( (lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_constituents0pr_e").c_str() );   
		SetBranchSingle( &nconstituents0pr, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_nconstituents0pr" );
	}
	SetBranchSingle( &rhoVal_, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rhoSW" );
	SetBranchSingle( &rhoVal_hand, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rhohand" );
	SetBranchSingle( &rhoVal_hand2, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rhohand2" );
	SetBranchSingle( &rhoVal_grid, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel + "_rhogrid" );

	////////////////////////////////////
	// CORRECTIONS ON THE FLY
	////////////////////////////////////     
	//// --- groomed jet label -------
	//mGroomedJet =  srcGroomedJet;

	//// --- Are we running over Monte Carlo ? --- 
	//if( iConfig.existsAs<bool>("runningOverMC") ) 
	//  runningOverMC_=iConfig.getParameter< bool >("runningOverMC");
	//else runningOverMC_= false;
	runningOverMC_= false;

	// --- Are we applying AK7 JEC to our groomed jet ? --- 
	//if( iConfig.existsAs<bool>("applyJECToGroomedJets") ) 
	//  applyJECToGroomedJets_=iConfig.getParameter< bool >("applyJECToGroomedJets");
	//else applyJECToGroomedJets_ = false;
	applyJECToGroomedJets_ = false;

	//// --- fastjet rho label -------
	//JetsFor_rho =  iConfig.getParameter<std::string>("srcJetsforRho") ; 
	//if(applyJECToGroomedJets_)
	//  JEC_GlobalTag_forGroomedJet = iConfig.getParameter<std::string>("JEC_GlobalTag_forGroomedJet") ; 
	JEC_GlobalTag_forGroomedJet="START53_V26";

	//// --- primary vertex -------
	//if(  iConfig.existsAs<edm::InputTag>("srcPrimaryVertex") )
	//  mPrimaryVertex = iConfig.getParameter<edm::InputTag>("srcPrimaryVertex"); 
	//else mPrimaryVertex =  edm::InputTag("offlinePrimaryVertices");

	// ---- setting up the jec on-the-fly from text files...    
	//    std::string fDir = "JEC/" + JEC_GlobalTag_forGroomedJet;   
	//std::string fDir = JEC_GlobalTag_forGroomedJet;   
	////std::string fDir = "JEC/"+JEC_GlobalTag_forGroomedJet;   
	////std::cout<<"JEC_GlobalTag_forGroomedJet="<<JEC_GlobalTag_forGroomedJet<<std::endl;
	//std::vector< JetCorrectorParameters > jecPars;
	//std::vector< JetCorrectorParameters > jecL1Pars;
	//std::vector< JetCorrectorParameters > jecL2Pars;
	//std::vector< JetCorrectorParameters > jecL3Pars;
	//std::vector< std::string > jecStr;

	//if(applyJECToGroomedJets_) {
	//	if (jetAlgorithmAdditonalLabel_=="_PF"){
	//		if(mJetAlgo == "AK" && fabs(mJetRadius-0.5)<0.001) {
	//			jecStr.push_back( fDir +  "_L1FastJet_AK5PF.txt" );
	//			jecStr.push_back( fDir + "_L2Relative_AK5PF.txt" );
	//			jecStr.push_back( fDir + "_L3Absolute_AK5PF.txt" );
	//			if (!runningOverMC_)
	//			  jecStr.push_back( fDir + "_L2L3Residual_AK5PF.txt" );
	//		}else{
	//			jecStr.push_back( fDir + "_L1FastJet_AK7PF.txt" );
	//			jecStr.push_back( fDir + "_L2Relative_AK7PF.txt" );
	//			jecStr.push_back( fDir + "_L3Absolute_AK7PF.txt" );
	//			if (!runningOverMC_)
	//			  jecStr.push_back( fDir + "_L2L3Residual_AK7PF.txt" );
	//		}
	//	}else //if (jetAlgorithmAdditonalLabel_=="_PFCHS")
	//	{
	//		if(mJetAlgo == "AK" && fabs(mJetRadius-0.5)<0.001) {
	//			jecStr.push_back( fDir +  "_L1FastJet_AK5PFchs.txt" );
	//			jecStr.push_back( fDir + "_L2Relative_AK5PFchs.txt" );
	//			jecStr.push_back( fDir + "_L3Absolute_AK5PFchs.txt" );
	//			if (!runningOverMC_)
	//			  jecStr.push_back( fDir + "_L2L3Residual_AK5PFchs.txt" );
	//		}else{
	//			jecStr.push_back( fDir + "_L1FastJet_AK7PFchs.txt" );
	//			jecStr.push_back( fDir + "_L2Relative_AK7PFchs.txt" );
	//			jecStr.push_back( fDir + "_L3Absolute_AK7PFchs.txt" );
	//			if (!runningOverMC_)
	//			  jecStr.push_back( fDir + "_L2L3Residual_AK7PFchs.txt" );
	//		}
	//	}

	//	for (unsigned int i = 0; i < jecStr.size(); ++i){
	//		//cout<<" JetCorrectorParameters( jecStr[i] )"<< jecStr[i]<<endl ; 
	//		JetCorrectorParameters* ijec = new JetCorrectorParameters( jecStr[i] );
	//		jecPars.push_back( *ijec );
	//		if( i==0 ){ jecL1Pars.push_back( *ijec );jecL2Pars.push_back( *ijec );jecL3Pars.push_back( *ijec ); }
	//		if( i==1 ){ jecL2Pars.push_back( *ijec );jecL3Pars.push_back( *ijec ); }
	//		if( i==2 ){ jecL3Pars.push_back( *ijec ); }
	//	}
	//	//jec_ = new FactorizedJetCorrector(jecPars);
	//	//jecL1_ = new FactorizedJetCorrector(jecL1Pars);
	//	//jecL2_ = new FactorizedJetCorrector(jecL2Pars);
	//	//jecL3_ = new FactorizedJetCorrector(jecL3Pars);

	//	//if (jetAlgorithmAdditonalLabel_=="_PF"){
	//	//	if(mJetAlgo == "AK" && fabs(mJetRadius-0.5)<0.001) {
	//	//		jecUnc_ = new JetCorrectionUncertainty( fDir + "_Uncertainty_AK5PF.txt" );
	//	//	}else{
	//	//		jecUnc_ = new JetCorrectionUncertainty( fDir + "_Uncertainty_AK7PF.txt" );
	//	//	}
	//	//}else //if (jetAlgorithmAdditonalLabel_=="_PFCHS")
	//	//{
	//	//	if(mJetAlgo == "AK" && fabs(mJetRadius-0.5)<0.001) {
	//	//		jecUnc_ = new JetCorrectionUncertainty( fDir + "_Uncertainty_AK5PFchs.txt" );
	//	//	}else{
	//	//		jecUnc_ = new JetCorrectionUncertainty( fDir + "_Uncertainty_AK7PFchs.txt" );
	//	//	}
	//	//}
	//}

	// specific configurations
	//if( iConfig.existsAs<double>("GroomedJet_JetChargeKappa") ) 
	//  mJetChargeKappa=iConfig.getParameter< double >("GroomedJet_JetChargeKappa");
	//else mJetChargeKappa = 0.3;
	mJetChargeKappa = 0.3;
	//if( iConfig.existsAs<int>("GroomedJet_QJetsPreclustering") ) 
	//  mQJetsPreclustering=iConfig.getParameter< int >("GroomedJet_QJetsPreclustering");
	//else mQJetsPreclustering = 30;
	mQJetsPreclustering = 30;
	//if( iConfig.existsAs<int>("GroomedJet_QJetsN") ) 
	//  mQJetsN=iConfig.getParameter< int >("GroomedJet_QJetsN");
	//else mQJetsN = 50;
	mQJetsN = 50;
	//if( iConfig.existsAs<double>("GroomedJet_NsubjettinessKappa") ) 
	//  mNsubjettinessKappa=iConfig.getParameter< double >("GroomedJet_NsubjettinessKappa");
	//else mNsubjettinessKappa = 1.0;
	mNsubjettinessKappa = 1.0;
	//if( iConfig.existsAs<bool>("GroomedJet_doQJets") ) 
	//  mDoQJets=iConfig.getParameter< bool >("GroomedJet_doQJets");
	//else mDoQJets = true;
	mDoQJets = false;

	// define charges of pdgIds
	neutrals.push_back( 22 ); neutrals.push_back( 130 ); neutrals.push_back( 310 ); neutrals.push_back( 311 ); neutrals.push_back( 111 ); 
	neutrals.push_back( 1 ); neutrals.push_back( 2 ); neutrals.push_back( 3 ); neutrals.push_back( 4 ); neutrals.push_back( 5 ); 
	neutrals.push_back( -1 ); neutrals.push_back( -2 ); neutrals.push_back( -3 ); neutrals.push_back( -4 ); neutrals.push_back( -5 ); 
	neutrals.push_back( 2112 ); neutrals.push_back( -2112 );
	neutrals.push_back( 3122 ); neutrals.push_back( -3122 );

	positives.push_back( 321 ); positives.push_back( 211 ); ; positives.push_back( -11 ); positives.push_back( -13); positives.push_back( 2212);
	negatives.push_back( -321 ); negatives.push_back( -211 ); negatives.push_back( 11 ); negatives.push_back( 13 );  negatives.push_back(-2212); 

	// define groomers
	//Trimming
	for(Int_t nradii=0;nradii<3;nradii++){
		Double_t tmp_radii=0.1+0.1*nradii;
		for(Int_t mpara=0;mpara<4;mpara++){
			Double_t tmp_para=0.03+0.01*mpara;
			//fastjet::Transformer* trimmer =new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, tmp_radii), fastjet::SelectorPtFractionMin(tmp_para)));
			fastjet::Transformer* trimmer =new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, tmp_radii), fastjet::SelectorPtFractionMin(tmp_para));
			GroomTool* gt0=new GroomTool("gttr", trimmer);
			vect_groomtools.push_back(*gt0);
		}
	}
	//Filtering
	for(Int_t nradii=0;nradii<3;nradii++){
		Double_t tmp_radii=0.1+0.1*nradii;
		for(Int_t mpara=0;mpara<3;mpara++){
			Int_t tmp_para=2+1*mpara;
			//fastjet::Transformer* filter  =new fastjet::Filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, tmp_radii), fastjet::SelectorNHardest(tmp_para)));
			fastjet::Transformer* filter  =new fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, tmp_radii), fastjet::SelectorNHardest(tmp_para));
			GroomTool* gt1=new GroomTool("gtfl", filter);
			vect_groomtools.push_back(*gt1);
		}
	}

	//Prunning
	for(Int_t n_zcut=0;n_zcut<3;n_zcut++){
		Double_t tmp_zcut=0.05+0.05*n_zcut;
		for(Int_t mpara=0;mpara<3;mpara++){
			Double_t tmp_para=0.4+0.1*mpara;
			fastjet::Transformer* pruner  =new fastjet::Pruner(fastjet::cambridge_algorithm, tmp_zcut, tmp_para);
			GroomTool* gt2=new GroomTool("gtpr", pruner);
			vect_groomtools.push_back(*gt2);
		}
	}


	for ( Int_t k=0; k< Int_t(vect_groomtools.size()); k++){
		vect_groomtools[k].SetBranchs(tree_, lableGen + "GroomedJet_" + jetAlgorithmLabel_ + additionalLabel);
	}

}


//////////////////////////////////////////////////////////////////
/////// Helper for above function ////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void ewk::GroomedJetFiller::SetBranchSingle( float* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"/F").c_str() );
	bnames.push_back( name );
}

void ewk::GroomedJetFiller::SetBranchSingle( double* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"/D").c_str() );
	bnames.push_back( name );
}

void ewk::GroomedJetFiller::SetBranchSingle( int* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"/I").c_str() );
	bnames.push_back( name );
}

void ewk::GroomedJetFiller::SetBranch( float* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"[6]/F").c_str() );
	bnames.push_back( name );
}


void ewk::GroomedJetFiller::SetBranch( int* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"[6]/I").c_str() );
	bnames.push_back( name );
}

//////////////////////////////////////////////////////////////////

void ewk::GroomedJetFiller::Init(){
	for (int j =0; j< NUM_JET_MAX; ++j) {
		jetmass_uncorr[j] = -1.;

		jetmass_tr_uncorr[j] = -1.;
		jetmass_tr1_uncorr[j] = -1.; // added
		jetmass_tr2_uncorr[j] = -1.; // added
		jetmass_tr3_uncorr[j] = -1.; // added

		jetmass_ft_uncorr[j] = -1.;
		jetmass_pr_uncorr[j] = -1.;
		jetmass_pr1_uncorr[j] = -1.;//added
		jetmass_pr2_uncorr[j] = -1.;// added
		jetmass_pr3_uncorr[j] = -1.;// added

		tau2tau1[j] = -1.;
		tau2tau1_shapesubtraction[j] = -1.;
		tau2tau1_trimmingshapesubtraction[j] = -1.;
		tau1[j] = -1.;
		tau2[j] = -1.;
		tau3[j] = -1.;
		tau4[j] = -1.;
		massdrop_pr_uncorr[j] = -1.; 
		number_jet_central=0.;
		jetpt_uncorr[j] = -1.;
		jetpt[j] = -1.;
		jeteta[j] = -10.;
		jetphi[j] = -10.;
		jete[j] = -1.;
		jetmass[j] = -1.;
		jetpt_JECL1[j] = -1.;
		jeteta_JECL1[j] = -10.;
		jetphi_JECL1[j] = -10.;
		jete_JECL1[j] = -1.;
		jetmass_JECL1[j] = -1.;
		jetmass_rhoArea[j] = -1.;
		jetmass_rhoGArea[j] = -1.;
		jetmass_rho4A[j] = -1.;
		jetmass_rhoG4A[j] = -1.;
		jetmass_rhom4Am[j] = -1.;
		jetmass_shapesubtraction[j] = -1.;
		jetmass_trimmingshapesubtraction[j] = -1.;

		jetpt_L1_rhoSW[j] = -1.;
		jetpt_L1_rhoHand[j] = -1.;
		jetpt_L1_rhoHand2[j] = -1.;
		jetpt_L1_rhoGrid[j] = -1.;

		jetpt_rho4A[j] = -1.;
		jetpt_shapesubtraction[j] = -1.;
		jetpt_trimmingshapesubtraction[j] = -1.;

		jetmass_tr[j] = -1.;
		jetmass_tr1[j] = -1.; //Added
		jetmass_tr2[j] = -1.; //Added
		jetmass_tr3[j] = -1.; //Added 

		jetmass_ft[j] = -1.;

		jetmass_pr[j] = -1.;
		jetmass_pr1[j] = -1.;// Added
		jetmass_pr2[j] = -1.;// added  
		jetmass_pr3[j] = -1.;// added 

		jetarea[j] = -1.;

		jetarea_tr[j] = -1.;
		jetarea_tr1[j] = -1.;// added 
		jetarea_tr2[j] = -1.;// added
		jetarea_tr3[j] = -1.;// added  

		jetarea_ft[j] = -1.;

		jetarea_pr[j] = -1.;
		jetarea_pr1[j] = -1.;//Added
		jetarea_pr2[j] = -1.;//added
		jetarea_pr3[j] = -1.;//added

		massdrop_pr[j] = -1.;
		massdrop_pr1[j] = -1.;// Added
		massdrop_pr2[j] = -1.;// added
		massdrop_pr3[j] = -1.; //added

		jetpt_tr_uncorr[j] = -1.;
		jetpt_tr1_uncorr[j] = -1.;// added 
		jetpt_tr2_uncorr[j] = -1.;// added
		jetpt_tr3_uncorr[j] = -1.;// added  

		jetpt_tr[j] = -1.;
		jetpt_tr1[j] = -1.; // added
		jetpt_tr2[j] = -1.; // added
		jetpt_tr3[j] = -1.; // added

		// jetpt_tr1_uncorr[j] = -1.;  ////
		// jetpt_tr1[j] = -1.; //added ////

		jeteta_tr[j] = -10.;
		jetphi_tr[j] = -10.;
		jete_tr[j] = -1.;

		jeteta_tr1[j] = -10.; //added
		jetphi_tr1[j] = -10.; //added
		jete_tr1[j] = -1.;  //added

		jeteta_tr2[j] = -10.; //added
		jetphi_tr2[j] = -10.; //added
		jete_tr2[j] = -1.;  //added

		jeteta_tr3[j] = -10.; //added
		jetphi_tr3[j] = -10.; //added
		jete_tr3[j] = -1.;  //added

		jetpt_ft_uncorr[j] = -1.;
		jetpt_ft[j] = -1.;
		jeteta_ft[j] = -10.;
		jetphi_ft[j] = -10.;
		jete_ft[j] = -1.;

		jetpt_pr_uncorr[j] = -1.;
		jetpt_pr[j] = -1.;
		jeteta_pr[j] = -10.;
		jetphi_pr[j] = -10.;
		jete_pr[j] = -1.;

		jetpt_pr1_uncorr[j] = -1.; //Added
		jetpt_pr1[j] = -1.;
		jeteta_pr1[j] = -10.;
		jetphi_pr1[j] = -10.;
		jete_pr1[j] = -1.;

		jetpt_pr2_uncorr[j] = -1.;// Added
		jetpt_pr2[j] = -1.;
		jeteta_pr2[j] = -10.;
		jetphi_pr2[j] = -10.;
		jete_pr2[j] = -1.;

		jetpt_pr3_uncorr[j] = -1.;// Added
		jetpt_pr3[j] = -1.;
		jeteta_pr3[j] = -10.;
		jetphi_pr3[j] = -10.;
		jete_pr3[j] = -1.;



		jetconstituents[j] = 0;
		jetcharge[j] = 0;

		prsubjet1_px[j] = 0.;
		prsubjet1_py[j] = 0.;
		prsubjet1_pz[j] = 0.;
		prsubjet1_e[j] = 0.;        
		prsubjet2_px[j] = 0.;
		prsubjet2_py[j] = 0.;
		prsubjet2_pz[j] = 0.;
		prsubjet2_e[j] = 0.;        

		for (int k = 0; k < 11; ++k){
			rcores[k][j] = -1.; ptcores[k][j] = -1.;
			planarflow[k][j] = -1.;
		}

		for (int k = 0; k < mQJetsN; ++k){
			qjetmass[k] = 0; qjetmassdrop[k] = 0;
		}
	}

	for(Int_t k=0;k<NUM_JETCLEANSING_MODE_MAX ;k++){
		jetmass_JetCleansing_DiffMode[k]= -1;
		jetpt_JetCleansing_DiffMode[k]  = -1;
		jeteta_JetCleansing_DiffMode[k]  = -1;
		jetphi_JetCleansing_DiffMode[k]  = -1;
		tau2tau1_JetCleansing_DiffMode[k]= -1;
		tau1_JetCleansing_DiffMode[k]= -1;
		tau2_JetCleansing_DiffMode[k]= -1;
		tau3_JetCleansing_DiffMode[k]= -1;
		tau4_JetCleansing_DiffMode[k]= -1;
	};
	nPV_ = 0.; // number of Primery vertex
	rhoVal_ = -99.; //rho of KT6PF from SW reco
	rhoVal_hand = -99.;
	rhoVal_hand2 = -99.;
	rhoVal_grid = -99.;

	// PID and charge
	PF_id_handle.clear();
	charge_handle_Gen.clear();


}

//////////////////////////////////////////////////////////////////

// ------------ method called to produce the data  ------------
void ewk::GroomedJetFiller::fill(
			//const edm::Event& iEvent,
			std::vector<fastjet::PseudoJet> FJparticles,
			bool doJetCleansing, std::vector<fastjet::PseudoJet> FJparticles_hardcharge, 
			std::vector<fastjet::PseudoJet> FJparticles_pileupcharge, 
			std::vector<fastjet::PseudoJet> FJparticles_fullneutral ) 
{

	// ------ start processing ------    
	if (FJparticles.size() < 1){ return;}
	//else{ std::cout << "FJparticles.size() = " << FJparticles.size() << std::endl; }

	Init();
	for ( Int_t k=0; k< Int_t(vect_groomtools.size()); k++){
		vect_groomtools[k].Init();
	}


	//nPV_ = getnPV(iEvent); // ------ get nPV: primary/secondary vertices------ 
	nPV_ = 0; // ------ get nPV: primary/secondary vertices------ 

	//for jet constituents PID and charge
	for(size_t i = 0; i < FJparticles.size(); ++ i) {
		fastjet::PseudoJet  P = FJparticles[i];
		//print_p4(P, "constituents", 1);
		PF_id_handle.push_back(P.user_info<PseudoJetUserInfo>().pdg_id());
		if(isGenJ)charge_handle_Gen.push_back(P.user_info<PseudoJetUserInfo>().charge());
	}

	// do re-clustering
	std::auto_ptr<fastjet::ClusterSequenceArea> thisClustering_area(new fastjet::ClusterSequenceArea(FJparticles, *mJetDef, *mAreaDefinition) );
	std::auto_ptr<fastjet::ClusterSequence> thisClustering_basic(new fastjet::ClusterSequence(FJparticles, *mJetDef));

	fastjet::Selector selected_eta = fastjet::SelectorAbsEtaMax(2.4);
	std::vector<fastjet::PseudoJet> out_jets       = sorted_by_pt( selected_eta(thisClustering_area->inclusive_jets(25.0)) );
	std::vector<fastjet::PseudoJet> out_jets_basic = sorted_by_pt( selected_eta(thisClustering_basic->inclusive_jets(25.0)) );

	fastjet::Subtractor* subtractor_medi=NULL; fastjet::Subtractor* subtractor_grid=NULL;

	//rhoVal_ = getrho(iEvent); // ------ get rho from SW reco --------    
	rhoVal_ = 1; // ------ get rho from SW reco --------    
	rhoVal_hand = getrho_Hand(FJparticles );


	rhoVal_hand2 = getrho_Hand2(FJparticles, &subtractor_medi ); // ------ get rho by hand2--------    
	rhoVal_grid = getrho_Grid(FJparticles, &subtractor_grid); // ------ get rho by grid--------    
	//std::cout<<"(kt6PF rho in SW) = "<<rhoVal_<<std::endl; std::cout << "(kt6PF rho by hand) = " << rhoVal_hand<<endl; std::cout<<"medi rho = "<<rhoVal_hand2<<endl; std::cout<<"grid rho = "<<rhoVal_grid<<endl;

	// define groomers
	fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
	fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)));
	fastjet::Filter trimmer2( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.03)));
	fastjet::Filter trimmer3( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.1), fastjet::SelectorPtFractionMin(0.03)));
	fastjet::Filter filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3)));
	fastjet::Pruner pruner(fastjet::cambridge_algorithm, 0.1, 0.5);
	fastjet::Pruner pruner1(fastjet::cambridge_algorithm, 0.05, 0.5);
	fastjet::Pruner pruner2(fastjet::cambridge_algorithm, 0.1, 0.75);
	fastjet::Pruner pruner3(fastjet::cambridge_algorithm, 0.05, 0.75);

	std::vector<fastjet::Transformer const *> transformers;
	transformers.push_back(&trimmer);
	transformers.push_back(&trimmer1);
	transformers.push_back(&trimmer2);
	transformers.push_back(&trimmer3);
	transformers.push_back(&filter);
	transformers.push_back(&pruner);
	transformers.push_back(&pruner1);
	transformers.push_back(&pruner2);
	transformers.push_back(&pruner3);

	// ----------------start loop on jets-----------------------
	number_jet_central = out_jets.size();
	for(int j = 0; j < number_jet_central && j<NUM_JET_MAX; j++) {
		//print_p4(out_jets.at(j), Form("%i jet",j));
		if (mSaveConstituents && j==0){
			if (out_jets_basic.at(j).constituents().size() >= 100) nconstituents0 = 100;
			else nconstituents0 = (int) out_jets_basic.at(j).constituents().size();
			std::vector<fastjet::PseudoJet> cur_constituents = sorted_by_pt(out_jets_basic.at(j).constituents());
			for (int aa = 0; aa < nconstituents0; aa++){        
				constituents0_eta[aa] = cur_constituents.at(aa).eta();
				constituents0_phi[aa] = cur_constituents.at(aa).phi();                
				constituents0_e[aa] = cur_constituents.at(aa).e();                                
			}
		}

		jetarea[j] = out_jets.at(j).area();
		//std::cout<<"jetarea="<<out_jets.at(j).area()<<std::endl; 
		jetconstituents[j] = out_jets_basic.at(j).constituents().size();

		jetmass_uncorr[j] = out_jets.at(j).m();
		jetpt_uncorr[j] = out_jets.at(j).pt();

		// SW JEC corretion
		TLorentzVector jet_corr = getCorrectedJet(out_jets.at(j), 100, 0);
		jetmass[j] = jet_corr.M();
		jetpt[j] = jet_corr.Pt();
		jeteta[j] = jet_corr.Eta();
		jetphi[j] = jet_corr.Phi();
		jete[j]   = jet_corr.Energy();

		// SW JEC L1 corretion
		TLorentzVector jet_corrL1 = getCorrectedJet(out_jets.at(j), 1, 0);
		jetmass_JECL1[j] = jet_corrL1.M();
		jetpt_JECL1[j] = jet_corrL1.Pt();
		jeteta_JECL1[j] = jet_corrL1.Eta();
		jetphi_JECL1[j] = jet_corrL1.Phi();
		jete_JECL1[j]   = jet_corrL1.Energy();
		//print_p4(out_jets.at(j),"  jet before    corr"); print_p4(jet_corr,      "--jet after JEC corr");

		//do rhoA correction
		fastjet::PseudoJet jet_rhoSW   = do_rhoA_correction(out_jets.at(j), rhoVal_, out_jets.at(j).area() );
		fastjet::PseudoJet jet_rhoHand = do_rhoA_correction(out_jets.at(j), rhoVal_hand, out_jets.at(j).area() );
		fastjet::PseudoJet jet_rhoHand2= do_rhoA_correction(out_jets.at(j), rhoVal_hand2, out_jets.at(j).area() );
		fastjet::PseudoJet jet_rhoGrid = do_rhoA_correction(out_jets.at(j), rhoVal_grid, out_jets.at(j).area() );

		jetpt_L1_rhoSW[j]   =jet_rhoSW.pt(); 
		jetpt_L1_rhoHand[j] =jet_rhoHand.pt(); 
		jetpt_L1_rhoHand2[j]=jet_rhoHand2.pt(); 
		jetpt_L1_rhoGrid[j] =jet_rhoGrid.pt(); 
		jetmass_rhoArea[j]= jet_rhoHand.m();
		jetmass_rhoGArea[j]= jet_rhoGrid.m();
		//std::cout<<"L1 p1={ "<<jetpt_L1_rhoSW[j]<<","<<jetpt_L1_rhoHand[j]<<","<<jetpt_L1_rhoHand2[j]<<","<<jetpt_L1_rhoGrid[j]<<"}"<<std::endl;
		//cout<<"jetmass_rhoA="<<jetmass_rhoArea[j]<<","<<jetmass_rhoGArea[j]<<endl; 

		//rho*area_4vector correction
		fastjet::PseudoJet jet_corr_medi = (*subtractor_medi)(out_jets.at(j));
		fastjet::PseudoJet jet_corr_grid = (*subtractor_grid)(out_jets.at(j));
		//print_p4(jet_corr_medi,      "--jet after medi rho4A corr"); print_p4(jet_corr_grid,      "--jet after grid rho4A corr");

		//print_p4(jet_corr_medi, "rho medi");
		jetpt_rho4A[j]= jet_corr_medi.pt();
		jetmass_rho4A[j]= jet_corr_medi.m();
		jetmass_rhoG4A[j]= jet_corr_grid.m();
		//cout<<"rho*A4: "<<jetpt_rho4A[j]<<" , "<<jetmass_rho4A[j]<<endl;

		//Generic Shape correction, arXiv 1211.2811
		Double_t rhom=0.0;
		do_GenericShapeSubtract_correction(out_jets.at(j), mBgeMedi.get(), jetpt_shapesubtraction[j], jetmass_shapesubtraction[j], tau2tau1_shapesubtraction[j], rhom);
		//cout<<"shapesubtraction: "<<jetpt_shapesubtraction[j]<<" , "<<jetmass_shapesubtraction[j]<<endl;

		//modified 4A for jet mass: p_sub = p_raw -rho*A4 - rhom*A4m
		fastjet::PseudoJet area_4vector=out_jets.at(j).area_4vector();
		fastjet::PseudoJet jet_corr_rhom4Am;
		jet_corr_rhom4Am.reset( jet_corr_medi.px(), jet_corr_medi.py(), jet_corr_medi.pz()-rhom*area_4vector.pz(), jet_corr_medi.E()-rhom*area_4vector.E());
		//print_p4(jet_corr_rhom4Am, "rhom4Am");
		jetmass_rhom4Am[j]= jet_corr_rhom4Am.m();

		//Trimming+Shapesubtraction
		fastjet::PseudoJet out_jets_trimming=trimmer( out_jets.at(j) );
		do_GenericShapeSubtract_correction(out_jets_trimming, mBgeMedi.get(), jetpt_trimmingshapesubtraction[j], jetmass_trimmingshapesubtraction[j], tau2tau1_trimmingshapesubtraction[j], rhom);
		//cout<<"trimmingshapesubtraction: "<<jetpt_trimmingshapesubtraction[j]<<" , "<<jetmass_trimmingshapesubtraction[j]<<endl;

		for ( Int_t k=0; k< Int_t(vect_groomtools.size());k++){
			vect_groomtools[k].Groom(out_jets.at(j), j); //cout<<"vect_groom "<<k<<endl;
		}


		// pruning, trimming, filtering  -------------
		int transctr = 0;
		//print_p4(out_jets.at(j), "raw_jet");
		for ( std::vector<fastjet::Transformer const *>::const_iterator 
					itransf = transformers.begin(), itransfEnd = transformers.end(); 
					itransf != itransfEnd; ++itransf ) {  

			fastjet::PseudoJet transformedJet = out_jets.at(j);
			transformedJet = (**itransf)(transformedJet);

			fastjet::PseudoJet transformedJet_basic = out_jets_basic.at(j);
			transformedJet_basic = (**itransf)(transformedJet_basic);

			if (transctr == 0){ // trimmed
				jetmass_tr_uncorr[j] = transformedJet.m();
				jetpt_tr_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_tr_corr = getCorrectedJet(transformedJet);
				if(jet_tr_corr.Energy()>0){
					jetmass_tr[j] = jet_tr_corr.M();
					jetpt_tr[j] = jet_tr_corr.Pt();
					jeteta_tr[j] = jet_tr_corr.Eta();
					jetphi_tr[j] = jet_tr_corr.Phi();
					jete_tr[j]   = jet_tr_corr.Energy();
					jetarea_tr[j] = transformedJet.area();
				}

			}
			else if (transctr == 1){        // trimmed1 Added
				jetmass_tr1_uncorr[j] = transformedJet.m();
				jetpt_tr1_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_tr1_corr = getCorrectedJet(transformedJet);
				if(jet_tr1_corr.Energy()>0){
					jetmass_tr1[j] = jet_tr1_corr.M();
					jetpt_tr1[j] = jet_tr1_corr.Pt();
					jeteta_tr1[j] = jet_tr1_corr.Eta();
					jetphi_tr1[j] = jet_tr1_corr.Phi();
					jete_tr1[j]   = jet_tr1_corr.Energy();
					jetarea_tr1[j] = transformedJet.area();
				}

			}

			else if(transctr ==2){        // Added
				jetmass_tr2_uncorr[j] = transformedJet.m();
				jetpt_tr2_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_tr2_corr = getCorrectedJet(transformedJet);
				if(jet_tr2_corr.Energy()>0){
					jetmass_tr2[j] = jet_tr2_corr.M();
					jetpt_tr2[j] = jet_tr2_corr.Pt();
					jeteta_tr2[j] = jet_tr2_corr.Eta();
					jetphi_tr2[j] = jet_tr2_corr.Phi();
					jete_tr2[j]   = jet_tr2_corr.Energy();
					jetarea_tr2[j] = transformedJet.area();
				}
			}

			else if(transctr==3){         // Added
				jetmass_tr3_uncorr[j] = transformedJet.m();
				jetpt_tr3_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_tr3_corr = getCorrectedJet(transformedJet);
				//print_p4(jet_tr3_corr,"jet_tr3_corr");
				if(jet_tr3_corr.Energy()>0){
					jetmass_tr3[j] = jet_tr3_corr.M();
					jetpt_tr3[j] = jet_tr3_corr.Pt();
					jeteta_tr3[j] = jet_tr3_corr.Eta();
					jetphi_tr3[j] = jet_tr3_corr.Phi();
					jete_tr3[j]   = jet_tr3_corr.Energy();
					jetarea_tr3[j] = transformedJet.area();
				}
			}

			else if (transctr == 4){ // filtered
				jetmass_ft_uncorr[j] = transformedJet.m();
				jetpt_ft_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_ft_corr = getCorrectedJet(transformedJet);
				if(jet_ft_corr.Energy()>0){
					jetmass_ft[j] = jet_ft_corr.M();
					jetpt_ft[j] = jet_ft_corr.Pt();
					jeteta_ft[j] = jet_ft_corr.Eta();
					jetphi_ft[j] = jet_ft_corr.Phi();
					jete_ft[j]   = jet_ft_corr.Energy();
					jetarea_ft[j] = transformedJet.area();                    
				}
			}
			else if (transctr == 5){ // pruned
				jetmass_pr_uncorr[j] = transformedJet.m();
				jetpt_pr_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_pr_corr = getCorrectedJet(transformedJet);
				if(jet_pr_corr.Energy()>0){
					jetmass_pr[j] = jet_pr_corr.M();
					jetpt_pr[j] = jet_pr_corr.Pt();
					jeteta_pr[j] = jet_pr_corr.Eta();
					jetphi_pr[j] = jet_pr_corr.Phi();
					jete_pr[j]   = jet_pr_corr.Energy();
					jetarea_pr[j] = transformedJet.area();
				}
			}        
			else if(transctr ==6){         // Added
				jetmass_pr1_uncorr[j] = transformedJet.m();
				jetpt_pr1_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_pr1_corr = getCorrectedJet(transformedJet);
				if(jet_pr1_corr.Energy()>0){
					jetmass_pr1[j] = jet_pr1_corr.M();
					jetpt_pr1[j] = jet_pr1_corr.Pt();
					jeteta_pr1[j] = jet_pr1_corr.Eta();
					jetphi_pr1[j] = jet_pr1_corr.Phi();
					jete_pr1[j]   = jet_pr1_corr.Energy();
					jetarea_pr1[j] = transformedJet.area();
				}
			}   
			else if(transctr ==7){         //Added
				jetmass_pr2_uncorr[j] = transformedJet.m();
				jetpt_pr2_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_pr2_corr = getCorrectedJet(transformedJet);
				if(jet_pr2_corr.Energy()>0){
					jetmass_pr2[j] = jet_pr2_corr.M();
					jetpt_pr2[j] = jet_pr2_corr.Pt();
					jeteta_pr2[j] = jet_pr2_corr.Eta();
					jetphi_pr2[j] = jet_pr2_corr.Phi();
					jete_pr2[j]   = jet_pr2_corr.Energy();
					jetarea_pr2[j] = transformedJet.area();
				}
			}  
			else if(transctr==8){         // Added
				jetmass_pr3_uncorr[j] = transformedJet.m();
				jetpt_pr3_uncorr[j] = transformedJet.pt();
				TLorentzVector jet_pr3_corr = getCorrectedJet(transformedJet);
				if(jet_pr3_corr.Energy()>0){
					jetmass_pr3[j] = jet_pr3_corr.M();
					jetpt_pr3[j] = jet_pr3_corr.Pt();
					jeteta_pr3[j] = jet_pr3_corr.Eta();
					jetphi_pr3[j] = jet_pr3_corr.Phi();
					jete_pr3[j]   = jet_pr3_corr.Energy();
					jetarea_pr3[j] = transformedJet.area();       

					jetarea_pr[j] = transformedJet.area();                    
					//decompose into requested number of subjets:
					if (transformedJet_basic.constituents().size() > 1){
						int nsubjetstokeep = 2;
						std::vector<fastjet::PseudoJet> subjets = transformedJet_basic.associated_cluster_sequence()->exclusive_subjets(transformedJet_basic,nsubjetstokeep);    

						//                    for (unsigned k = 0; k < subjets.size(); k++) {
						//                        std::cout << "subjet " << k << ": mass = " << subjets.at(k).m() << " and pt = " << subjets.at(k).pt() << std::endl;
						//                    }
						TLorentzVector sj1( subjets.at(0).px(),subjets.at(0).py(),subjets.at(0).pz(),subjets.at(0).e());
						TLorentzVector sj2( subjets.at(1).px(),subjets.at(1).py(),subjets.at(1).pz(),subjets.at(1).e());     
						prsubjet1_px[j] = subjets.at(0).px(); prsubjet1_py[j] = subjets.at(0).py(); prsubjet1_pz[j] = subjets.at(0).pz(); prsubjet1_e[j] = subjets.at(0).e();
						prsubjet2_px[j] = subjets.at(1).px(); prsubjet2_py[j] = subjets.at(1).py(); prsubjet2_pz[j] = subjets.at(1).pz(); prsubjet2_e[j] = subjets.at(1).e();                    
						TLorentzVector fullj = sj1 + sj2; 


						if (subjets.at(0).m() >= subjets.at(1).m()){
							massdrop_pr_uncorr[j] = subjets.at(0).m()/transformedJet.m();
							massdrop_pr[j] = (subjets.at(0).m()/jetmass_pr[j]);                        
						}
						else{
							massdrop_pr_uncorr[j] = subjets.at(1).m()/transformedJet.m();
							massdrop_pr[j] = (subjets.at(1).m()/jetmass_pr[j]);                                    
						}



						//Added
						if (subjets.at(0).m() >= subjets.at(1).m()){
							massdrop_pr1_uncorr[j] = subjets.at(0).m()/transformedJet.m();
							massdrop_pr1[j] = (subjets.at(0).m()/jetmass_pr1[j]);
						}
						else{
							massdrop_pr1_uncorr[j] = subjets.at(1).m()/transformedJet.m();
							massdrop_pr1[j] = (subjets.at(1).m()/jetmass_pr1[j]);
						}

						//Added
						if (subjets.at(0).m() >= subjets.at(1).m()){
							massdrop_pr2_uncorr[j] = subjets.at(0).m()/transformedJet.m();
							massdrop_pr2[j] = (subjets.at(0).m()/jetmass_pr2[j]);
						}
						else{
							massdrop_pr2_uncorr[j] = subjets.at(1).m()/transformedJet.m();
							massdrop_pr2[j] = (subjets.at(1).m()/jetmass_pr2[j]);
						}

						// Added
						if (subjets.at(0).m() >= subjets.at(1).m()){
							massdrop_pr3_uncorr[j] = subjets.at(0).m()/transformedJet.m();
							massdrop_pr3[j] = (subjets.at(0).m()/jetmass_pr3[j]);
						}
						else{
							massdrop_pr3_uncorr[j] = subjets.at(1).m()/transformedJet.m();
							massdrop_pr3[j] = (subjets.at(1).m()/jetmass_pr3[j]);
						}

					}

					// pruned tests
					if (mSaveConstituents && j==0){
						if (transformedJet_basic.constituents().size() >= 100) nconstituents0pr = 100;
						else nconstituents0pr = (int) transformedJet_basic.constituents().size();
						std::vector<fastjet::PseudoJet> cur_constituentspr = sorted_by_pt(transformedJet_basic.constituents());
						for (int aa = 0; aa < nconstituents0pr; aa++){        
							constituents0pr_eta[aa] = cur_constituentspr.at(aa).eta();
							constituents0pr_phi[aa] = cur_constituentspr.at(aa).phi();                
							constituents0pr_e[aa] = cur_constituentspr.at(aa).e();                                
						}
					}
				}
			}
			else{ std::cout << "error in number of transformers" << std::endl;}                    
			transctr++;
		}        

		//nsubjettiness
		get_nsubjettiness(out_jets.at(j), tau1[j], tau2[j], tau3[j], tau4[j], tau2tau1[j]);
		//cout<<"raw tau2tau1 = "<<tau2tau1[j]<<endl;

		// cores computation  -------------
		//std::cout<< "Beging the core computation" << endl;
		std::vector<fastjet::PseudoJet> constits = thisClustering_area->constituents(out_jets.at(j));
		for (int kk = 0; kk < 11; ++kk){
			double coreCtr = (double) kk;    
			if (coreCtr < mJetRadius*10.){
				float tmpm = 0, tmppt = 0;
				computeCore( constits, coreCtr/10., tmpm, tmppt );
				if (tmpm > 0){ rcores[kk][j] = tmpm/out_jets.at(j).m(); }
				if (tmppt > 0) ptcores[kk][j] = tmppt/out_jets.at(j).pt();
			}
		}
		//std::cout<< "Ending the core computation" << endl;

		//std::cout<< "Beging the planarflow computation" << endl;
		//planarflow computation
		for (int kk = 0; kk < 11; ++kk){
			double coreCtr = (double) (kk + 1);
			if (coreCtr < mJetRadius*10.){
				float tmppflow = 0;
				computePlanarflow(constits,coreCtr/10.,out_jets.at(j),mJetAlgo,tmppflow);
				planarflow[kk][j] = tmppflow;
			}
		}
		//std::cout<< "Ending the planarflow computation" << endl;

		// qjets computation  -------------
		if ((mDoQJets)&&(j == 0)){ // do qjets only for the hardest jet in the event!
			double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1);                
			QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity);
			fastjet::JetDefinition qjet_def(&qjet_plugin);
			vector<fastjet::PseudoJet> constits;
			unsigned int nqjetconstits = out_jets_basic.at(j).constituents().size();
			if (nqjetconstits < (unsigned int) mQJetsPreclustering) constits = out_jets_basic.at(j).constituents();
			else constits = out_jets_basic.at(j).associated_cluster_sequence()->exclusive_subjets_up_to(out_jets_basic.at(j),mQJetsPreclustering);
			for(unsigned int ii = 0 ; ii < (unsigned int) mQJetsN ; ii++){
				fastjet::ClusterSequence qjet_seq(constits, qjet_def);
				vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(20.0));
				//if(mJetAlgo == "AK" && fabs(mJetRadius-0.5)<0.001) inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(20.0));

				if (inclusive_jets2.size()>0) {
					qjetmass[ii] = inclusive_jets2[0].m();
					if (inclusive_jets2[0].constituents().size() > 1){
						vector<fastjet::PseudoJet> subjets_qjet = qjet_seq.exclusive_subjets(inclusive_jets2[0],2);
						if (subjets_qjet.at(0).m() >= subjets_qjet.at(1).m()){
							qjetmassdrop[ii] = (subjets_qjet.at(0).m()/inclusive_jets2[0].m());                        
						}
						else{
							qjetmassdrop[ii] = (subjets_qjet.at(1).m()/inclusive_jets2[0].m());                                    
						}
					}
					else{
						qjetmassdrop[ii] = 1.;
					}
				}else{
					qjetmassdrop[ii] = 1.;
				}

				//qjet_seq.delete_self_when_unused();
			}
		}
		// jet charge try (?) computation  -------------
		std::vector< float > pdgIds;
		for (unsigned ii = 0; ii < out_jets_basic.at(j).constituents().size(); ii++){
			for (unsigned jj = 0; jj < FJparticles.size(); jj++){
				//std::cout << ii << ", " << jj << ": " << FJparticles.at(jj).pt() << ", " << out_jets_basic.at(j).constituents().at(ii).pt() << std::endl;
				if (FJparticles.at(jj).pt() == out_jets_basic.at(j).constituents().at(ii).pt()){
					pdgIds.push_back(PF_id_handle.at(jj));
					//if(!isGenJ) {
					//  if(mJetAlgo == "AK" && fabs(mJetRadius-0.5)<0.001) {
					//  pdgIds.push_back(PF_id_handle_AK5.at(jj));
					//  }else{
					//  pdgIds.push_back(PF_id_handle->at(jj));
					//  }
					//  }else{
					//  pdgIds.push_back(PF_id_handle_Gen.at(jj));
					//  }
					break;
				}
			}
		}
		//jetcharge[j] = computeJetCharge(out_jets_basic.at(j).constituents(),pdgIds,out_jets_basic.at(j).e());
	}

	delete subtractor_medi;
	delete subtractor_grid;
	// jet JetCleansing
	if(doJetCleansing) { DoJetCleansing( *mJetDef, FJparticles, FJparticles_hardcharge, FJparticles_pileupcharge, FJparticles_fullneutral, out_jets.at(0)); }
}  

double ewk::GroomedJetFiller::getJEC(double curJetEta, double curJetPt, double curJetE, double curJetArea){
	//// Jet energy corrections, something like this...
	//jec_->setJetEta( curJetEta );
	//jec_->setJetPt ( curJetPt );
	//jec_->setJetE  ( curJetE );
	//jec_->setJetA  ( curJetArea );
	//jec_->setRho   ( rhoVal_ );
	//jec_->setNPV   ( nPV_ );
	//double corr = jec_->getCorrection();
	//return corr;
	return 1;
}
double ewk::GroomedJetFiller::getJECL1(double curJetEta, double curJetPt, double curJetE, double curJetArea){
	//// Jet energy corrections, something like this...
	//jecL1_->setJetEta( curJetEta );
	//jecL1_->setJetPt ( curJetPt );
	//jecL1_->setJetE  ( curJetE );
	//jecL1_->setJetA  ( curJetArea );
	//jecL1_->setRho   ( rhoVal_ );
	//jecL1_->setNPV   ( nPV_ );
	//double corr = jecL1_->getCorrection();
	//return corr;
	return 1;
}
double ewk::GroomedJetFiller::getJECL2(double curJetEta, double curJetPt, double curJetE, double curJetArea){
	//// Jet energy corrections, something like this...
	//jecL2_->setJetEta( curJetEta );
	//jecL2_->setJetPt ( curJetPt );
	//jecL2_->setJetE  ( curJetE );
	//jecL2_->setJetA  ( curJetArea );
	//jecL2_->setRho   ( rhoVal_ );
	//jecL2_->setNPV   ( nPV_ );
	//double corr = jecL2_->getCorrection();
	//return corr;
	return 1;
}
double ewk::GroomedJetFiller::getJECL3(double curJetEta, double curJetPt, double curJetE, double curJetArea){
	//// Jet energy corrections, something like this...
	//jecL3_->setJetEta( curJetEta );
	//jecL3_->setJetPt ( curJetPt );
	//jecL3_->setJetE  ( curJetE );
	//jecL3_->setJetA  ( curJetArea );
	//jecL3_->setRho   ( rhoVal_ );
	//jecL3_->setNPV   ( nPV_ );
	//double corr = jecL3_->getCorrection();
	//return corr;
	return 1;
}

TLorentzVector ewk::GroomedJetFiller::getCorrectedJet(fastjet::PseudoJet& jet, Int_t jeclevel, bool debug) {
	double jecVal = 1.0;

	if(applyJECToGroomedJets_ && !isGenJ){ 
		if (jeclevel>=10) jecVal = getJEC( jet.eta(), jet.pt(), jet.e(), jet.area() );
		if (jeclevel==1) jecVal = getJECL1( jet.eta(), jet.pt(), jet.e(), jet.area() );
		if (jeclevel==2) jecVal = getJECL2( jet.eta(), jet.pt(), jet.e(), jet.area() );
		if (jeclevel==3) jecVal = getJECL3( jet.eta(), jet.pt(), jet.e(), jet.area() );
	}   

	if(debug){
		std::cout<<"jecVal="<<jecVal<<std::endl;
		std::cout<<"jec LA="<<getJEC( jet.eta(), jet.pt(), jet.e(), jet.area() )<<std::endl;
		std::cout<<"jec L1="<<getJECL1( jet.eta(), jet.pt(), jet.e(), jet.area() )<<std::endl;
		std::cout<<"jec L2="<<getJECL2( jet.eta(), jet.pt(), jet.e(), jet.area() )<<std::endl;
		std::cout<<"jec L3="<<getJECL3( jet.eta(), jet.pt(), jet.e(), jet.area() )<<std::endl;
	}
	TLorentzVector jet_corr(jet.px() * jecVal, 
				jet.py() * jecVal, 
				jet.pz() * jecVal, 
				jet.e() * jecVal);
	return jet_corr;
}

fastjet::PseudoJet ewk::GroomedJetFiller::getScaledJet(fastjet::PseudoJet& jet, double scale) {
	fastjet::PseudoJet jet_scale(jet.px() * scale, jet.py() * scale, jet.pz() * scale, jet.e()  * scale);
	return jet_scale;
}

void ewk::GroomedJetFiller::computeCore( std::vector<fastjet::PseudoJet> constits, double Rval, float &m_core, float &pt_core ){

	fastjet::JetDefinition jetDef_rcore(fastjet::cambridge_algorithm, Rval);
	fastjet::ClusterSequence thisClustering(constits, jetDef_rcore);

	std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering.inclusive_jets(0.0));
	m_core = out_jets.at(0).m();
	pt_core = out_jets.at(0).pt();
	//thisClustering.delete_self_when_unused();
}

void ewk::GroomedJetFiller::computePlanarflow(std::vector<fastjet::PseudoJet> constits, double Rval, fastjet::PseudoJet jet,std::string mJetAlgo, float &planarflow){

	fastjet::JetDefinition jetDef_rplanarflow(fastjet::cambridge_algorithm,Rval);
	if (mJetAlgo == "AK") jetDef_rplanarflow.set_jet_algorithm( fastjet::antikt_algorithm );
	else if (mJetAlgo == "CA") jetDef_rplanarflow.set_jet_algorithm( fastjet::cambridge_algorithm );
	//else throw cms::Exception("GroomedJetFiller") << " unknown jet algorithm " << std::endl;
	else cout<<"GroomedJetFiller" << " unknown jet algorithm " << std::endl;
	fastjet::ClusterSequence thisClustering(constits, jetDef_rplanarflow);

	//reclustering jets
	std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering.inclusive_jets(0.0));

	//leading sub jet constits mass not equal Zero
	float mJ = jet.m();
	if(mJ != 0)
	{
		std::vector<fastjet::PseudoJet> subconstits = thisClustering.constituents(out_jets.at(0)); 

		TLorentzVector jetp4;
		//jetp4.SetPxPyPzE(out_jets.at(0).px(),out_jets.at(0).py(),out_jets.at(0).pz(),out_jets.at(0).e());
		jetp4.SetPxPyPzE(jet.px(),jet.py(),jet.pz(),jet.e());

		TVector3 zaxis = jetp4.Vect().Unit();
		TVector3 zbeam(0, 0, 1);

		//Transverse component (X, Y) relative to the jet(Z) axis
		TVector3 xaxis = (zaxis.Cross(zbeam)).Unit();
		TVector3 yaxis = (xaxis.Cross(zaxis)).Unit();

		double I[3][3];
		for (int i = 0; i < 3; i ++) for (int j = 0; j < 3; j ++) I[i][j] = 0;

		int matrixsize = subconstits.size();

		for(int k = 0; k < matrixsize; k++)
		{   
			TLorentzVector tmpjetk;
			tmpjetk.SetPxPyPzE(subconstits.at(k).px(),subconstits.at(k).py(),subconstits.at(k).pz(),subconstits.at(k).e());
			float tmp_px = tmpjetk.Vect().Dot(xaxis);
			float tmp_py = tmpjetk.Vect().Dot(yaxis);
			//Avoid Too Samll Energy
			if(subconstits.at(k).e() >= 0.001){

				I[1][1] += tmp_px * tmp_px / (mJ * subconstits.at(k).e());
				I[1][2] += tmp_px * tmp_py / (mJ * subconstits.at(k).e());
				I[2][1] += tmp_py * tmp_px / (mJ * subconstits.at(k).e());
				I[2][2] += tmp_py * tmp_py / (mJ * subconstits.at(k).e());

			}
		}

		//From arXiv 1012.2077
		planarflow = 4*(I[1][1]*I[2][2] - I[1][2]*I[2][1])/((I[1][1]+I[2][2])*(I[1][1]+I[2][2])); 
	}
	//thisClustering.delete_self_when_unused();
}

float ewk::GroomedJetFiller::computeJetCharge( std::vector<fastjet::PseudoJet> constits, std::vector<float> pdgIds, float Ejet ){

	float val = 0.;
	for (unsigned int i = 0; i < pdgIds.size(); i++){
		float qq ;
		if(isGenJ) {
			qq = charge_handle_Gen.at(i);
		}else{
			qq = getPdgIdCharge( pdgIds.at(i) );
		}
		val += qq*pow(constits.at(i).e(),mJetChargeKappa);
	}
	val /= Ejet;
	return val;

}

float ewk::GroomedJetFiller::getPdgIdCharge( float fid ){

	float qq = -99.;
	int id = (int) fid;
	if (std::find(neutrals.begin(), neutrals.end(), id) != neutrals.end()){
		qq = 0.;
	}
	else if (std::find(positives.begin(), positives.end(), id) != positives.end()){
		qq = 1.;
	}
	else if (std::find(negatives.begin(), negatives.end(), id) != negatives.end()){
		qq = -1.;
	}
	else{
		cout<<"id="<<id<<endl;
		//BREAK("unknown PDG id");
		//throw cms::Exception("GroomedJetFiller") << " unknown PDG id " << id << std::endl;
		cout<<"GroomedJetFiller" << " unknown PDG id " << id << std::endl;
	}
	return qq;
}

Double_t ewk::GroomedJetFiller::getnPV( ){//const edm::Event& iEvent ){
	// ------ get nPV: primary/secondary vertices------ 
	double nPVval = 0;
	//edm::Handle <edm::View<reco::Vertex> > recVtxs;
	//iEvent.getByLabel( mPrimaryVertex, recVtxs);
	//for(unsigned int ind=0;ind<recVtxs->size();ind++){
	//	if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>=4) 
	//				&& (fabs((*recVtxs)[ind].z())<=24.0) &&  
	//				((*recVtxs)[ind].position().Rho()<=2.0) ) {
	//		nPVval += 1;
	//	}
	//}
	return nPVval;
}


Double_t ewk::GroomedJetFiller::getrho(){// const edm::Event& iEvent ){ // ------ get rho --------    
	//edm::Handle<double> rho;
	//const edm::InputTag eventrho(JetsFor_rho, "rho");//kt6PFJetsPFlow
	//iEvent.getByLabel(eventrho,rho);
	//return *rho;
	return 0.5;
}

Double_t ewk::GroomedJetFiller::getrho_Hand(std::vector<fastjet::PseudoJet>  FJparticles ){ // ------ get rho --------    

	std::auto_ptr<double> rho_on_the_fly(new double(0.0));
	std::auto_ptr<double> sigma_on_the_fly(new double(0.0));
	double mean_area = 0;

	//ClusterSequencePtr fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceArea( FJparticles, fastjet::JetDefinition(fastjet::kt_algorithm,0.6), *mAreaDefinition ) );
	ClusterSequencePtr fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceArea( FJparticles, 
					fastjet::JetDefinition(fastjet::kt_algorithm,0.6), fastjet::VoronoiAreaSpec(0.9) ) );
	fastjet::ClusterSequenceAreaBase const* clusterSequenceWithArea = 
		dynamic_cast<fastjet::ClusterSequenceAreaBase const *> ( &*fjClusterSeq_ );

	double rhoEtaMax=4.4;
	RangeDefPtr fjRangeDef_ = RangeDefPtr( new fastjet::RangeDefinition(rhoEtaMax) );

	clusterSequenceWithArea->get_median_rho_and_sigma(*fjRangeDef_,false,*rho_on_the_fly,*sigma_on_the_fly,mean_area);
	if((*rho_on_the_fly < 0)|| (std::isnan(*rho_on_the_fly))) { //edm::LogError("BadRho") << "rho_on_the_fly value is " << *rho_on_the_fly << " area:" << mean_area << " and n_empty_jets: " << clusterSequenceWithArea->n_empty_jets(*fjRangeDef_) << " with range " << fjRangeDef_->description() <<". Setting rho_on_the_fly to rezo.";
		*rho_on_the_fly = 0;
	}

	//fjClusterSeq_->delete_self_when_unused();
	return *rho_on_the_fly;
}

Double_t ewk::GroomedJetFiller::getrho_Hand2(std::vector<fastjet::PseudoJet>  FJparticles, fastjet::Subtractor** subtractor )
{
	double rhoEtaMax=4.4;
	mBgeMedi.reset(new fastjet::JetMedianBackgroundEstimator(fastjet::SelectorAbsRapMax(rhoEtaMax), fastjet::JetDefinition(fastjet::kt_algorithm,0.6), fastjet::VoronoiAreaSpec(0.9) ) );
	mBgeMedi->set_particles(FJparticles);
	if(*subtractor) delete *subtractor;
	*subtractor= new fastjet::Subtractor(mBgeMedi.get());
	return mBgeMedi->rho();
}
Double_t ewk::GroomedJetFiller::getrho_Grid(std::vector<fastjet::PseudoJet>  FJparticles, fastjet::Subtractor** subtractor )
{
	double rhoEtaMax=4.4;
	mBgeGrid.reset(new fastjet::GridMedianBackgroundEstimator(rhoEtaMax, 0.55) );
	mBgeGrid->set_particles(FJparticles);
	if(*subtractor) delete *subtractor;
	*subtractor= new fastjet::Subtractor(mBgeGrid.get());
	return mBgeGrid->rho();
}

fastjet::PseudoJet ewk::GroomedJetFiller::do_rhoA_correction(fastjet::PseudoJet jet_origin, double rho, double area){
	double pt_new= jet_origin.pt() - rho*area;
	double factor_jec = pt_new / jet_origin.pt();
	fastjet::PseudoJet jet_new = getScaledJet(jet_origin, factor_jec);
	return jet_new;
}




void ewk::GroomedJetFiller::do_GenericShapeSubtract_correction(fastjet::PseudoJet jet_origin, fastjet::BackgroundEstimatorBase* bge_rho, float& jetpt_new, float& jetmass_new, float& tau2tau1_shapesubtraction, Double_t& rhom)
{
	Mass jetshape_mass; 
	Pt   jetshape_pt;

	double beta = mNsubjettinessKappa; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
	double R0 = mJetRadius; // Characteristic jet radius for normalization	      
	double Rcut = mJetRadius; // maximum R particles can be from axis to be included in jet	      
	NSubjettinessRatio tau21(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	//NSubjettinessRatio tau21(2);

	//GenericSubtractor gen_sub(mBgeMedi);
	GenericSubtractor gen_sub(bge_rho);
	GenericSubtractorInfo info;

	//std::cout << gen_sub.description() << std::endl; std::cout << setprecision(4);
	// uncomment this if you also want rho_m to be estimated (using the
	// same background estimator)
	gen_sub.use_common_bge_for_rho_and_rhom(true);

	// compute the subtracted shape, and retrieve additional information

	jetpt_new = gen_sub(jetshape_pt, jet_origin, info);
	/*cout<< "uncorr pt = " << jetshape_pt(jet_origin) << endl;
	  cout<< "jetshape_pt corr  = " << jetpt_new << endl;
	  cout << "  rho  = " << info.rho() << endl;
	  cout << "  rhom = " << info.rhom() << endl;
	  cout << "  1st derivative: " << info.first_derivative() << endl;
	  cout << "  2nd derivative: " << info.second_derivative() << endl;
	  cout << "  unsubtracted: " << info.unsubtracted() << endl;
	  cout << "  1st order: " << info.first_order_subtracted() << endl;
	  cout << "# step used: " << info.ghost_scale_used() << endl;
	  */


	jetmass_new = gen_sub(jetshape_mass, jet_origin, info);
	/*cout<< "uncorr mass = " << jetshape_mass(jet_origin) << endl;
	  cout<< "jetshape_mass corr  = " << jetmass_new << endl;
	  cout << "  rho  = " << info.rho() << endl;
	  cout << "  rhom = " << info.rhom() << endl;
	  cout << "  1st derivative: " << info.first_derivative() << endl;
	  cout << "  2nd derivative: " << info.second_derivative() << endl;
	  cout << "  unsubtracted: " << info.unsubtracted() << endl;
	  cout << "  1st order: " << info.first_order_subtracted() << endl;
	  cout << "# step used: " << info.ghost_scale_used() << endl;
	  */



	//cout<<"unsubtracted jets: "<<tau21(jet_origin)<<endl;
	//cout<<"  subtracted jets: "<<gen_sub(tau21, jet_origin)<<endl;
	tau2tau1_shapesubtraction = gen_sub(tau21, jet_origin); 
	rhom=info.rhom();
}

void  ewk::GroomedJetFiller::get_nsubjettiness(fastjet::PseudoJet jet_origin, float &tau1, float &tau2, float &tau3, float &tau4, float & tau2tau1){
	// Defining Nsubjettiness parameters
	double beta = mNsubjettinessKappa; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
	double R0 = mJetRadius; // Characteristic jet radius for normalization	      
	double Rcut = mJetRadius; // maximum R particles can be from axis to be included in jet	      

	fastjet::Nsubjettiness nSub1KT(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	fastjet::Nsubjettiness nSub2KT(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	fastjet::Nsubjettiness nSub3KT(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	fastjet::Nsubjettiness nSub4KT(4, Njettiness::onepass_kt_axes, beta, R0, Rcut);

	tau1 = nSub1KT(jet_origin);
	tau2 = nSub2KT(jet_origin);
	tau3 = nSub3KT(jet_origin);
	tau4 = nSub4KT(jet_origin);
	tau2tau1 = tau2/tau1;
}


//Jet Cleansing
JetCleanser ewk::GroomedJetFiller::makeJVFCleanser(fastjet::JetDefinition subjet_def, std::string projectmode, double fcut, int nsj )//projectmode: CMS or ATLAS
{
	JetCleanser::input_mode tmp_projectmode=JetCleanser::input_nc_separate;
	if(projectmode=="ATLAS") tmp_projectmode=JetCleanser::input_nc_together;
	JetCleanser tmpCleanser(subjet_def, JetCleanser::jvf_cleansing, tmp_projectmode);

	if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
	else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
	else if( fcut>0 && nsj>0 ){ tmpCleanser.SetGroomingParameters(fcut,nsj);}

	return tmpCleanser;
}
JetCleanser ewk::GroomedJetFiller::makeLinearCleanser(fastjet::JetDefinition subjet_def, double linear_para0,std::string projectmode, double fcut, int nsj)//projectmode: CMS or ATLAS
{
	JetCleanser::input_mode tmp_projectmode=JetCleanser::input_nc_separate;
	if(projectmode=="ATLAS") tmp_projectmode=JetCleanser::input_nc_together;
	JetCleanser tmpCleanser(subjet_def, JetCleanser::linear_cleansing, tmp_projectmode);
	tmpCleanser.SetLinearParameters(linear_para0);

	if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
	else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
	else if( fcut>0 && nsj>0 ){ tmpCleanser.SetGroomingParameters(fcut,nsj);}

	return tmpCleanser;
}
JetCleanser ewk::GroomedJetFiller::makeGausCleanser(fastjet::JetDefinition subjet_def, double gaus_para0, double gaus_para1, double gaus_para2, double gaus_para3, std::string projectmode, double fcut, int nsj )//projectmode: CMS or ATLAS
{
	JetCleanser::input_mode tmp_projectmode=JetCleanser::input_nc_separate;
	if(projectmode=="ATLAS") tmp_projectmode=JetCleanser::input_nc_together;
	JetCleanser tmpCleanser(subjet_def, JetCleanser::gaussian_cleansing, tmp_projectmode);
	tmpCleanser.SetGaussianParameters(gaus_para0, gaus_para1, gaus_para2, gaus_para3);

	if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
	else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
	else if( fcut>0 && nsj>0 ){ tmpCleanser.SetGroomingParameters(fcut,nsj);}

	return tmpCleanser;
}

void ewk::GroomedJetFiller::DoJetCleansing(fastjet::JetDefinition jetDef, std::vector<fastjet::PseudoJet> FJparticles, 
			std::vector<fastjet::PseudoJet> FJparticles_hardcharge,
			std::vector<fastjet::PseudoJet> FJparticles_pileupcharge,
			std::vector<fastjet::PseudoJet> FJparticles_fullneutral,
			fastjet::PseudoJet recoJet)
{
	//cout<<"DoJetCleansing:  "<<" FJparticles.size="<<FJparticles.size()<<" FJparticles_hardcharge.size="<<FJparticles_hardcharge.size()<<" FJparticles_pileupcharge.size="<<FJparticles_pileupcharge.size()<<" FJparticles_fullneutral.size="<<FJparticles_fullneutral.size()<<endl;
	// find jets
	vector< vector<fastjet::PseudoJet> > sets;
	sets.push_back( FJparticles );           // calorimeter cells
	sets.push_back( FJparticles_hardcharge );   // tracks from primary interaction
	sets.push_back( FJparticles_pileupcharge );       // tracks from pileup
	sets.push_back( FJparticles_fullneutral );   // neutral particles

	// collect jets
	vector< vector<fastjet::PseudoJet> > jet_sets = ClusterSets(jetDef, FJparticles, sets, 25.0);
	vector<fastjet::PseudoJet> jets_plain     = jet_sets[0];
	vector<fastjet::PseudoJet> jets_tracks_LV = jet_sets[1];
	vector<fastjet::PseudoJet> jets_tracks_PU = jet_sets[2];
	vector<fastjet::PseudoJet> jets_neutrals  = jet_sets[3];

	//cout<<"jets_tracks_LV.size()="<<jets_tracks_LV.size()<<endl;
	//cout<<"jets_tracks_PU.size()="<<jets_tracks_PU.size()<<endl;
	//cout<<"jets_neutrals.size()="<<jets_neutrals.size()<<endl;

	Int_t num_matching_with_reco=-1;//all reco jet eta<2.4
	for(Int_t i=0; i< Int_t(jets_plain.size()); i++){
		//cout<<i<<"st jet:"<<endl; print_p4( jets_plain[i], "jets_plain");
		if ( isMatching(jets_plain[i], recoJet, 0.1 ) ){
			num_matching_with_reco=i;
			break;
		}
	}
	//if( num_matching_with_reco==-1) throw cms::Exception("JetCleansing Failed") << " couldn't matching"<< std::endl;
	if( num_matching_with_reco==-1) cout<<"JetCleansing Failed" << " couldn't matching"<< std::endl;
	// Jet cleansing
	vector<JetCleanser> jetcleanser_vect;
	fastjet::JetDefinition subjet_def_kt03(fastjet::kt_algorithm, 0.3);
	fastjet::JetDefinition subjet_def_kt025(fastjet::kt_algorithm, 0.25);
	fastjet::JetDefinition subjet_def_kt02(fastjet::kt_algorithm, 0.2);
	fastjet::JetDefinition subjet_def_kt015(fastjet::kt_algorithm, 0.15);

	JetCleanser jetcleanser1=makeJVFCleanser(subjet_def_kt03, "CMS"); jetcleanser_vect.push_back(jetcleanser1);
	JetCleanser jetcleanser2=makeJVFCleanser(subjet_def_kt02, "CMS"); jetcleanser_vect.push_back(jetcleanser2);
	JetCleanser jetcleanser3=makeLinearCleanser(subjet_def_kt03,0.55, "CMS"); jetcleanser_vect.push_back(jetcleanser3);
	JetCleanser jetcleanser4=makeLinearCleanser(subjet_def_kt02,0.55, "CMS"); jetcleanser_vect.push_back(jetcleanser4);
	JetCleanser jetcleanser5=makeLinearCleanser(subjet_def_kt03,0.60, "CMS"); jetcleanser_vect.push_back(jetcleanser5);
	JetCleanser jetcleanser6=makeLinearCleanser(subjet_def_kt02,0.60, "CMS"); jetcleanser_vect.push_back(jetcleanser6);


	/*// define groomers
	  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)));
	  fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
	  fastjet::Filter trimmer2( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.03)));
	  fastjet::Filter trimmer3( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.1), fastjet::SelectorPtFractionMin(0.03)));

	  fastjet::Filter filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3)));

	  if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
	  else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
	  */

	// linear + trimming
	for(Int_t trim_par=0; trim_par<=8; trim_par++){
		JetCleanser jetcleanser_lineartrim1=makeLinearCleanser(subjet_def_kt03, 0.55, "CMS", 0.01+0.01*trim_par); jetcleanser_vect.push_back(jetcleanser_lineartrim1);
		JetCleanser jetcleanser_lineartrim2=makeLinearCleanser(subjet_def_kt02, 0.55, "CMS", 0.01+0.01*trim_par); jetcleanser_vect.push_back(jetcleanser_lineartrim2);
	}
	// linear + filtering
	for(Int_t filt_par=0; filt_par<=5; filt_par++){
		JetCleanser jetcleanser_linearfilt1=makeLinearCleanser(subjet_def_kt03, 0.55, "CMS", -1, 1+1*filt_par); jetcleanser_vect.push_back(jetcleanser_linearfilt1);
		JetCleanser jetcleanser_linearfilt2=makeLinearCleanser(subjet_def_kt02, 0.55, "CMS", -1, 1+1*filt_par); jetcleanser_vect.push_back(jetcleanser_linearfilt2);
	}

	/*
	// jvf
	JetCleanser jetcleanser01=makeJVFCleanser(subjet_def_kt03, "CMS"); jetcleanser_vect.push_back(jetcleanser01);
	JetCleanser jetcleanser02=makeJVFCleanser(subjet_def_kt025, "CMS"); jetcleanser_vect.push_back(jetcleanser02);
	JetCleanser jetcleanser03=makeJVFCleanser(subjet_def_kt02, "CMS"); jetcleanser_vect.push_back(jetcleanser03);
	JetCleanser jetcleanser04=makeJVFCleanser(subjet_def_kt015, "CMS"); jetcleanser_vect.push_back(jetcleanser04);
	// linear
	for(Int_t linear_par=0; linear_par<=30; linear_par++){
	JetCleanser jetcleanser1=makeLinearCleanser(subjet_def_kt03,0.4+0.01*linear_par, "CMS"); jetcleanser_vect.push_back(jetcleanser1);
	}
	for(Int_t linear_par=0; linear_par<=30; linear_par++){
	JetCleanser jetcleanser1=makeLinearCleanser(subjet_def_kt025,0.4+0.01*linear_par, "CMS"); jetcleanser_vect.push_back(jetcleanser1);
	}
	for(Int_t linear_par=0; linear_par<=30; linear_par++){
	JetCleanser jetcleanser1=makeLinearCleanser(subjet_def_kt02,0.4+0.01*linear_par, "CMS"); jetcleanser_vect.push_back(jetcleanser1);
	}
	for(Int_t linear_par=0; linear_par<=30; linear_par++){
	JetCleanser jetcleanser1=makeLinearCleanser(subjet_def_kt015,0.4+0.01*linear_par, "CMS"); jetcleanser_vect.push_back(jetcleanser1);
	} 
	// gaus
	JetCleanser jetcleanser2=makeGausCleanser(subjet_def_kt03, 0.67, 0.62, 0.20, 0.25, "CMS");
	jetcleanser_vect.push_back(jetcleanser2);
	*/
	if( int(jetcleanser_vect.size()) >= NUM_JETCLEANSING_MODE_MAX ){ std::cout<<"Error! Jet Cleansing Mode is too many!"<<endl; BREAK(); }
	//cout<<"jetcleanser_vect.size()="<<jetcleanser_vect.size()<<endl;
	for(Int_t j=0;j<int(jetcleanser_vect.size());j++){
		//fastjet::PseudoJet tmp_cleansed_jet = jetcleanser_vect[j]( jets_neutrals[0].constituents(), jets_tracks_LV[0].constituents(), jets_tracks_PU[0].constituents() );//CMS model
		fastjet::PseudoJet tmp_cleansed_jet = jetcleanser_vect[j]( jets_neutrals[num_matching_with_reco].constituents(), jets_tracks_LV[num_matching_with_reco].constituents(), jets_tracks_PU[num_matching_with_reco].constituents() );//CMS model
		jetmass_JetCleansing_DiffMode[j]= tmp_cleansed_jet.m();
		jetpt_JetCleansing_DiffMode[j]  = tmp_cleansed_jet.pt();
		jeteta_JetCleansing_DiffMode[j]  = tmp_cleansed_jet.eta();
		jetphi_JetCleansing_DiffMode[j]  = tmp_cleansed_jet.phi();

		float tmp1; float tmp2; float tmp3; float tmp4; float tmp5;
		get_nsubjettiness(tmp_cleansed_jet, tmp1, tmp2, tmp3, tmp4, tmp5 );
		tau1_JetCleansing_DiffMode[j]= tmp1;
		tau2_JetCleansing_DiffMode[j]= tmp2;
		tau3_JetCleansing_DiffMode[j]= tmp3;
		tau4_JetCleansing_DiffMode[j]= tmp4;
		tau2tau1_JetCleansing_DiffMode[j]= tmp5;


		//if( tmp_cleansed_jet.pt() <10 ){
		//cout<<"jet cleansing pt="<<tmp_cleansed_jet.pt()<<endl;
		//cout<<"jet cleansing mass="<<tmp_cleansed_jet.m()<<endl;
		//cout<<"jet cleansing tau1="<<tmp1<<endl;
		//cout<<"jet cleansing tau2="<<tmp2<<endl;
		//cout<<"jet cleansing tau3="<<tmp3<<endl;
		//cout<<"jet cleansing tau4="<<tmp4<<endl;
		//cout<<"jet cleansing tau2tau1="<<tmp5<<endl;
		//}
		//cout<<jetcleanser_vect[j].description()<<endl<<"jet mass="<< tmp_cleansed_jet.m()<<" jet pt="<< tmp_cleansed_jet.pt()<<endl; cout<<"================"<<endl;
		//cout<<j<<" mass="<< tmp_cleansed_jet.m()<<" pt="<< tmp_cleansed_jet.pt()<<endl; 
	};


}

