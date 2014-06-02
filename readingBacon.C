#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPFPart.hh"
#include "TGenParticle.hh"
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

	//TFile* fIn = new TFile("/eos/uscms/store/user/ntran/PUPPI/bacon/qcd300-470_62x_PU40BX50/ntuple_1_1_IR6.root");
	TFile* fIn = new TFile("ntuple_1_1_IR6.root");
	TFile fout("fout.root","recreate");
	TTree* tree = (TTree*) fIn->Get("Events");

	TTree* myTree = new TTree("DiJets","DiJets");
	GroomedJetFiller* gf_GEN=new GroomedJetFiller("genGroomedJetFiller", myTree, "CA10", "_GEN", 1);
	GroomedJetFiller* gf_PF=new GroomedJetFiller("GroomedJetFiller", myTree, "CA10", "_PF");
	GroomedJetFiller* gf_PFCHS=new GroomedJetFiller("GroomedJetFiller", myTree, "CA10", "_PFCHS");


	//RECO
	TClonesArray *fPFPart = new TClonesArray("baconhep::TPFPart");
	tree->SetBranchAddress("PFPart",       &fPFPart);
	//GEN
	TClonesArray *fGenParticle = new TClonesArray("baconhep::TGenParticle");
	tree->SetBranchAddress("GenParticle",   &fGenParticle);


	for(int i0 = 2; i0 <100; i0++)
	//for(int i0 = 0; i0 < tree->GetEntriesFast(); i0++) 
	{ 

		tree->GetEntry(i0);

		std::vector<fastjet::PseudoJet> fjinputs_pfs_noLep_CHS; fjinputs_pfs_noLep_CHS.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_charge;  fjinputs_pfs_charge.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_neutral; fjinputs_pfs_neutral.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_noLep_noCHS; fjinputs_pfs_noLep_noCHS.clear();
		std::vector<fastjet::PseudoJet> fjinputs_pfs_PileUp; fjinputs_pfs_PileUp.clear();
		// fjinputs_pfs_noLep_noCHS = fastjet::sorted_by_pt(fjinputs_pfs_noLep_noCHS);
		// Now, we get 3 collections:
		// fjinputs_pfs_PileUp: PF charged, PileUp  
		// fjinputs_pfs_charge: PF charged, NoPileUp
		// fjinputs_pfs_neutral:PF neutral
		// fjinputs_pfs_noLep_CHS: PF neutral + charge(vtxID<1), after CHS
		// fjinputs_pfs_noLep_noCHS: PF neutral + charge, no CHS


		std::cout << "i0 = " << i0 << ", and N PF candidates = " << fPFPart->GetEntriesFast()  <<", Gen_Particles = " << fGenParticle->GetEntriesFast() << std::endl;

		for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){

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
			fjinputs_pfs_noLep_noCHS.push_back(tmp_psjet);
			if(charge==0){
				fjinputs_pfs_neutral.push_back(tmp_psjet);
				fjinputs_pfs_noLep_CHS.push_back(tmp_psjet);
			}else{
				if(vtxId<1){//0 and -1 are for 1st primary vertex; >0 for pileup
					fjinputs_pfs_charge.push_back(tmp_psjet);
					fjinputs_pfs_noLep_CHS.push_back(tmp_psjet);
				}else{
					fjinputs_pfs_PileUp.push_back(tmp_psjet );
				}
			}
		}

		fjinputs_pfs_noLep_noCHS= fastjet::sorted_by_pt(fjinputs_pfs_noLep_noCHS);
		fjinputs_pfs_noLep_CHS  = fastjet::sorted_by_pt(fjinputs_pfs_noLep_CHS);
		fjinputs_pfs_charge     = fastjet::sorted_by_pt(fjinputs_pfs_charge);
		fjinputs_pfs_neutral    = fastjet::sorted_by_pt(fjinputs_pfs_neutral);
		fjinputs_pfs_PileUp     = fastjet::sorted_by_pt(fjinputs_pfs_PileUp);

		//Gen
		std::vector<fastjet::PseudoJet> fjinputs_gen; fjinputs_gen.clear();

		for( int i1 = 0; i1 < fGenParticle->GetEntriesFast(); i1++){

			baconhep::TGenParticle *pPartTmp = (baconhep::TGenParticle*)((*fGenParticle)[i1]);

			Int_t  status = pPartTmp->status;
			double pt = pPartTmp->pt;
			double y  = pPartTmp->y;
			double phi= pPartTmp->phi;
			double mass=pPartTmp->mass;

			double pdgId = pPartTmp->pdgId;
			int charge = status;

			fastjet::PseudoJet tmp_psjet;
			tmp_psjet.reset_PtYPhiM(pt, y, phi, mass);
			tmp_psjet.set_user_info(new PseudoJetUserInfo(pdgId, charge) );

			if (status ==1 || status ==2) {
				if( pt > 1e-5){
					fjinputs_gen.push_back(tmp_psjet);
					print_p4(tmp_psjet,"tmp", 1);
				}
			}else{
				continue;
			}
		}
		fjinputs_gen= fastjet::sorted_by_pt(fjinputs_gen);


		//GEN
		gf_GEN->fill( fjinputs_gen);
		//RECO
		bool doJetCleansing=1;
		gf_PF->fill( fjinputs_pfs_noLep_noCHS, doJetCleansing, fjinputs_pfs_charge,fjinputs_pfs_PileUp,fjinputs_pfs_neutral);
		gf_PFCHS->fill( fjinputs_pfs_noLep_CHS);
		myTree->Fill();

	}

	fout.Write();
	myTree->Print();

}
