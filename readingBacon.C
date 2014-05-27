#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPFPart.hh"
#include <vector>
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "../../fjcontrib-1.009/include/fastjet/contrib/EnergyCorrelator.hh"

#include "tools.h"
#include "GroomedJetFiller.h"

using namespace fastjet;
using namespace fastjet::contrib;
using namespace baconhep;
using namespace ewk;

void readingBacon(){
    
    //TFile* fIn = new TFile("/eos/uscms/store/user/ntran/PUPPI/bacon/rsgww1000_62x_PU40BX50/ntuple_1_1_VQC.root");
    TFile* fIn = new TFile("ntuple_1_1_VQC.root");
	TFile fout("fout.root","recreate");
    TTree* tree = (TTree*) fIn->Get("Events");

	TTree* myTree = new TTree("Tree","DiJets");
	GroomedJetFiller* gf=new GroomedJetFiller("GroomedJetFiller", myTree, "CA10", "_PF");


    TClonesArray *fPFPart = new TClonesArray("baconhep::TPFPart");
    tree->SetBranchAddress("PFPart",       &fPFPart);

    //for(int i0 = 0; i0 <10 ; i0++)
    for(int i0 = 0; i0 < tree->GetEntriesFast(); i0++) 
	{ 
        
		tree->GetEntry(i0);

		std::vector<fastjet::PseudoJet> fjinputs_pfs; fjinputs_pfs.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_charge;  fjinputs_pfs_charge.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_neutral; fjinputs_pfs_neutral.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_noLep_noCHS; fjinputs_pfs_noLep_noCHS.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_PileUp; fjinputs_pfs_PileUp.clear();

		std::cout << "i0 = " << i0 << ", and N PF candidates = " << fPFPart->GetEntriesFast() << std::endl;

		for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//9,entries loop,fill the vector fjinputs_pfs with PF fjinputs_pfs
			//cout<<"p"<<i1<<endl;

			baconhep::TPFPart *pPartTmp = (baconhep::TPFPart*)((*fPFPart)[i1]);

			double Px = pPartTmp->pt*cos(pPartTmp->phi);
			double Py = pPartTmp->pt*sin(pPartTmp->phi);
			double theta = 2*atan(exp(-pPartTmp->eta)); //eta = -ln(tan(theta/2))
			double Pz = pPartTmp->pt/tan(theta);
			double E  = pPartTmp->e;
			double pdgId = pPartTmp->pfType;
			int charge = pPartTmp->q;
			int vtxId = pPartTmp->vtxId;
			fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
			tmp_psjet.set_user_info(new PseudoJetUserInfo(pdgId, charge) );
			//cout<<"tmp_"<<i1<<" pdg="<<pdgId<<" charge="<<charge<<" ("<<Px<<","<<Py<<","<<Pz<<","<<E<<")"<<endl;
			//print_p4(tmp_psjet,"tmp", 1 );
			fjinputs_pfs.push_back(tmp_psjet);
			fjinputs_pfs_noLep_noCHS.push_back(tmp_psjet);
			if(charge==0){
				fjinputs_pfs_neutral.push_back(tmp_psjet);
			}else{
				if(vtxId==0){
					fjinputs_pfs_charge.push_back(tmp_psjet);
				}else{
					fjinputs_pfs_PileUp.push_back(tmp_psjet );
				}
			}
/*
	}
	fjinputs_pfs_noLep_noCHS = fastjet::sorted_by_pt(fjinputs_pfs_noLep_noCHS);
	// Now, we get 3 collections:
	// fjinputs_pfs_PileUp: PF charged, PileUp  
	// fjinputs_pfs_charge: PF charged, NoPileUp
	// fjinputs_pfs_neutral:PF neutral
	// fjinputs_pfs: PF neutral + charge, after CHS
	// fjinputs_pfs_noLep_noCHS: PF neutral + charge, no CHS
*/




		}//9,entries loop ,fill the vector fjinputs_pfs with PFparticles


		bool doJetCleansing=1;
		gf->fill( fjinputs_pfs_noLep_noCHS, doJetCleansing, fjinputs_pfs_charge,fjinputs_pfs_PileUp,fjinputs_pfs_neutral);
		myTree->Fill();

	}

	fout.Write();
	myTree->Print();

}
