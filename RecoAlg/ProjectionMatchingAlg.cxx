////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ProjectionMatchingAlg.h"

pma::ProjectionMatchingAlg::ProjectionMatchingAlg(const fhicl::ParameterSet& pset)
{
  this->reconfigure(pset); 
}

pma::ProjectionMatchingAlg::~ProjectionMatchingAlg(void)
{
}

void pma::ProjectionMatchingAlg::reconfigure(const fhicl::ParameterSet& p)
{
  //fBlurN               = p.get<int>   ("BlurN");
}

