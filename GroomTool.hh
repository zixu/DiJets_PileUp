/*
 * =====================================================================================
 * 
 *       Filename:  GroomTool.hh
 * 
 *    Description:  Groom Tool: prunning, filtering, trimming
 * 
 *        Version:  1.0
 *        Created:  05/19/14 14:56:56 CDT
 *       Revision:  none
 *       Compiler:  gcc, root
 * 
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
 * 
 * =====================================================================================
 */

#ifndef  GROOMTOOL_INC
#define  GROOMTOOL_INC

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLorentzVector.h"


#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"

using namespace fastjet;
using namespace std;

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

#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;


class GroomTool{
	public:
		static const int NUM_JET_MAX = 4;
		fastjet::Transformer* groomer;
		std::string groomer_label;
		float jetmass_groomed[NUM_JET_MAX];
		float jetpt_groomed[NUM_JET_MAX];
		float jeteta_groomed[NUM_JET_MAX];
		float jetphi_groomed[NUM_JET_MAX];
		float jete_groomed[NUM_JET_MAX];
		float jetarea_groomed[NUM_JET_MAX];
		float tau2tau1[NUM_JET_MAX];

		GroomTool(){};
		GroomTool(string in_groom_label, fastjet::Transformer* in_groomer);
		~GroomTool(){};
		void Init();
		void SetBranchs(TTree *t1, string title);
		void SetBranch(TTree *t1, float* obs, string title, string obs_label);
		void Groom(fastjet::PseudoJet jet_raw, Int_t number_jet);
		fastjet::PseudoJet getJECJet(fastjet::PseudoJet jet_raw);

};





#endif   /* ----- #ifndef GROOMTOOL_INC  ----- */

