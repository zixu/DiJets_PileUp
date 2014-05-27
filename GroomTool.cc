/*
 * =====================================================================================
 *
 *       Filename:  GroomTool.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/19/14 14:58:05 CDT
 *       Revision:  none
 *       Compiler:  gcc, root
 *
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
 *
 * =====================================================================================
 */
#include "GroomTool.hh"
GroomTool::GroomTool(string in_groom_label, fastjet::Transformer* in_groomer):groomer_label(in_groom_label) { 
	groomer = in_groomer;
	Init();
}

void GroomTool::Init(){
	for(Int_t i=0;i< NUM_JET_MAX; i++){
		jetmass_groomed[i]=-1;
		jetpt_groomed[i]=-1;
		jeteta_groomed[i]=-10;
		jetphi_groomed[i]=-10;
		jete_groomed[i]=-1;
		jetarea_groomed[i]=-1;
		tau2tau1[i]=-1;
	}
}

void GroomTool::Groom(fastjet::PseudoJet jet_raw, Int_t number_jet){
	//cout<<"Begin Groom"<<endl;
	if(number_jet< NUM_JET_MAX){
		//print_p4(jet_raw, "jet_raw");
		fastjet::PseudoJet jet_groomed= (*groomer)(jet_raw);
		//cout<<"groomer="<<groomer->description()<<endl;
		jetarea_groomed[number_jet]=jet_groomed.area();
		//print_p4(jet_groomed, "jet_groomed");

		//TLorentzVector jet_groomed_corr=getJECJet( jet_groomed );
		fastjet::PseudoJet jet_groomed_corr=getJECJet( jet_groomed );
		//print_p4(jet_groomed_corr, "jet_groomed_corr");
		jetmass_groomed[number_jet]=jet_groomed_corr.m();
		jetpt_groomed[number_jet]=jet_groomed_corr.pt();
		jeteta_groomed[number_jet]=jet_groomed_corr.eta();
		jetphi_groomed[number_jet]=jet_groomed_corr.phi();
		jete_groomed[number_jet]=jet_groomed_corr.e();

		//if(number_jet==1) cout<<groomer_label<<" #jet = "<<number_jet<<" area= "<<jetarea_groomed[number_jet]<<" mass= "<<jetmass_groomed[number_jet]<<endl; 
	}
	//cout<<"End Groom"<<endl;
}
fastjet::PseudoJet GroomTool::getJECJet(fastjet::PseudoJet jet_raw){
	//TLorentzVector jet_corr = getCorrectedJet(out_jets.at(j), 100, 0);
	if( jet_raw.e()<=0.){ fastjet::PseudoJet jtmp=jet_raw; return jtmp; }

	Double_t jecVal=1;
	fastjet::PseudoJet jet_corr( jet_raw.px()*jecVal, jet_raw.py()*jecVal, jet_raw.pz()*jecVal, jet_raw.e()*jecVal);
	return jet_corr;
}

void GroomTool::SetBranch(TTree *t1, float* obs, string title, string obs_label){
	t1->Branch( Form("%s_%s_%s", title.c_str(), obs_label.c_str(), groomer_label.c_str()), 
				obs, Form("%s_%s_%s[%i]/F", title.c_str(), obs_label.c_str(), groomer_label.c_str(), NUM_JET_MAX) );
}

void GroomTool::SetBranchs(TTree *t1, string title){
	SetBranch(t1, jetpt_groomed, title, "pt");
	SetBranch(t1, jetmass_groomed, title, "mass");
	SetBranch(t1, jeteta_groomed, title, "eta");
	SetBranch(t1, jetphi_groomed, title, "phi");
	SetBranch(t1, jete_groomed, title, "e");
	SetBranch(t1, jetarea_groomed, title, "area");
}


