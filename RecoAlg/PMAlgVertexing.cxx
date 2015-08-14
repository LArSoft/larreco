////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgVertexing
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), August 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "RecoAlg/PMAlgVertexing.h"

#include "RecoAlg/PMAlg/Utilities.h"

pma::PMAlgVertexing::PMAlgVertexing(const fhicl::ParameterSet& pset)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

pma::PMAlgVertexing::~PMAlgVertexing(void)
{
	cleanTracks();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::cleanTracks(void)
{
	for (auto t : fOutTracks) delete t;
	fOutTracks.clear();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::reconfigure(const fhicl::ParameterSet& pset)
{
	fInputVtxDist2D = pset.get< double >("InputVtxDist2D");
	fInputVtxDistY = pset.get< double >("InputVtxDistY");
}
// ------------------------------------------------------

size_t pma::PMAlgVertexing::run(
	const std::vector< pma::Track3D* >& trk_input,
	const std::vector< TVector3 >& vtx_input)
{
	cleanTracks();

	for (auto t : trk_input) fOutTracks.push_back(new pma::Track3D(*t));

	return 0;
}
// ------------------------------------------------------

std::vector< recob::Vertex > pma::PMAlgVertexing::getVertices(void) const
{
	std::vector< recob::Vertex > vtxs;

	return vtxs;
}
// ------------------------------------------------------

