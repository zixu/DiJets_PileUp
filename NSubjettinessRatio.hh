// $Id: ExampleShapes.hh 3001 2013-01-29 10:41:40Z soyez $
//
// Copyright (c) 2012-, Matteo Cacciari, Jihun Kim, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __NSubjettinessRatio__HH__
#define __NSubjettinessRatio__HH__

#include "fastjet/FunctionOfPseudoJet.hh"
#include "../../fjcontrib-1.009/GenericSubtractor/ShapeWithPartition.hh"
#include "../../fjcontrib-1.009/GenericSubtractor/ShapeWithComponents.hh"
#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h"

//FASTJET_BEGIN_NAMESPACE

//namespace contrib{
class NSubjettinessRatio : public contrib::ShapeWithComponents{
	public:
		NSubjettinessRatio(int N, Njettiness::AxesMode in_mode, double in_beta, double in_R0, double in_Rcut) : _N(N), _mode(in_mode), _beta(in_beta), _R0(in_R0), _Rcut(in_Rcut) { assert(_N>1); }
		virtual std::string description() const{ return "N-subjettiness ratio from components";}

		virtual unsigned int n_components() const { return 2;}

		virtual std::vector<double> components(const fastjet::PseudoJet &jet) const{
			std::vector<double> comp(n_components());
			comp[0] = fastjet::Nsubjettiness(_N  , _mode, _beta, _R0, _Rcut)(jet);
			comp[1] = fastjet::Nsubjettiness(_N-1, _mode, _beta, _R0, _Rcut)(jet);

			return comp;
		}

		virtual double result_from_components(const std::vector <double> &components) const{
			return components[0]/components[1];
		}

		virtual FunctionOfPseudoJet<double> * component_shape(unsigned int index) const{
			return new fastjet::Nsubjettiness(_N-index, _mode, _beta, _R0, _Rcut);
		}

	protected:
		const int _N;
		const Njettiness::AxesMode _mode;
		const double _beta;
		const double _R0;
		const double _Rcut;

};

//} // namespace contrib

//FASTJET_END_NAMESPACE

#endif // __NSubjettinessRatio__HH__
