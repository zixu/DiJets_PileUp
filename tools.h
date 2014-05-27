/*
 * =====================================================================================
 * 
 *       Filename:  tools.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  08/07/13 02:04:15 CDT
 *       Revision:  none
 *       Compiler:  gcc, root
 * 
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
 * 
 * =====================================================================================
 */

#ifndef  TOOLS_INC
#define  TOOLS_INC
#include "fastjet/PseudoJet.hh"
#include "TROOT.h"
#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"

void print_p4(fastjet::PseudoJet tmpJ, std::string tmpName="",bool extra_info=0);
void BREAK(std::string info="");
Bool_t isMatching( fastjet::PseudoJet j1, fastjet::PseudoJet j2, Double_t deltaR=0.3);
#endif   /* ----- #ifndef TOOLS_INC  ----- */

