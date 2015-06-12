////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ProjectionMatchingAlg.h"

#include "RecoAlg/PMAlg/Utilities.h"

pma::ProjectionMatchingAlg::ProjectionMatchingAlg(const fhicl::ParameterSet& pset)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

pma::ProjectionMatchingAlg::~ProjectionMatchingAlg(void)
{
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fOptimizationEps = p.get< double >("OptimizationEps");
	fFineTuningEps = p.get< double >("FineTuningEps");

	fTrkValidationDist2D = 1.0;
	fHitTestingDist2D = 0.4;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView) const
{
	double step = 0.3;
	double max_d = fTrkValidationDist2D;
	double d2, max_d2 = max_d * max_d;
	unsigned int nAll = 0, nPassed = 0;

	TVector3 p(trk.front()->Point3D());
	for (size_t i = 0; i < trk.Nodes().size() - 1; i++)
	{
		unsigned int tpc = trk.Nodes()[i]->TPC();
		unsigned int cryo = trk.Nodes()[i]->Cryo();

		TVector3 vNext(trk.Nodes()[i + 1]->Point3D());
		TVector3 vThis(trk.Nodes()[i]->Point3D());

		TVector3 dc(vNext); dc -= vThis;
		dc *= step / dc.Mag();

		double f = pma::GetSegmentProjVector(p, vThis, vNext);
		while (f < 1.0)
		{
			TVector2 p2d = pma::GetProjectionToPlane(p, testView, tpc, cryo);

			for (const auto& h : hits)
				if (h->WireID().Plane == testView)
			{
				d2 = pma::Dist2(p2d, pma::WireDriftToCm(h->WireID().Wire, h->PeakTime(), testView, tpc, cryo));
				if (d2 < max_d2) { nPassed++; break; }
			}
			nAll++;

			p += dc; f = pma::GetSegmentProjVector(p, vThis, vNext);
		}
		p = vNext;
	}
	if (nAll > 3) // validate actually only if there are some hits in testView
	{
		double v = nPassed / (double)nAll;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  trk fraction ok: " << v;
		return v;
	}
	else return 1.0;
}
// ------------------------------------------------------

size_t pma::ProjectionMatchingAlg::getSegCount(size_t trk_size)
{
	int nSegments = trk_size / (2 * (int)sqrt( trk_size ));

	if (nSegments > 1) return (size_t)nSegments;
	else return 1;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildTrack(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2,
	unsigned int testView) const
{
	pma::Track3D* trk = new pma::Track3D(); // track candidate
	trk->AddHits(hits_1);
	trk->AddHits(hits_2);

	mf::LogVerbatim("ProjectionMatchingAlg") << "track size: " << trk->size();
	std::vector< unsigned int > tpcs = trk->TPCs();
	for (size_t t = 0; t < tpcs.size(); ++t)
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "  tpc:" << tpcs[t];
	}
	mf::LogVerbatim("ProjectionMatchingAlg")
		<< "  #coll:" << trk->NHits(geo::kZ)
		<< "  #ind2:" << trk->NHits(geo::kV)
		<< "  #ind1:" << trk->NHits(geo::kU);

	size_t nSegments = getSegCount(trk->size());
	size_t nNodes = (size_t)( nSegments - 1 ); // n nodes to add

	mf::LogVerbatim("ProjectionMatchingAlg") << "  initialize trk";
	trk->Initialize();

	mf::LogVerbatim("ProjectionMatchingAlg") << "  optimize trk (" << nSegments << " seg)";
	if (nNodes) trk->Optimize(nNodes, fOptimizationEps);   // build nodes
	double g = trk->Optimize(0, fFineTuningEps);           // final tuning
	mf::LogVerbatim("ProjectionMatchingAlg") << "  done, g = " << g;
	mf::LogVerbatim("ProjectionMatchingAlg") << "  sort trk";

	trk->SortHits();
	return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::extendTrack(
	const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	bool add_nodes) const
{
	pma::Track3D* copy = new pma::Track3D(trk);
	copy->AddHits(hits);

	mf::LogVerbatim("ProjectionMatchingAlg") << "ext. track size: " << copy->size()
		<< "  #coll:" << copy->NHits(geo::kZ)
		<< "  #ind2:" << copy->NHits(geo::kV)
		<< "  #ind1:" << copy->NHits(geo::kU);

	if (add_nodes)
	{
		size_t nSegments = getSegCount(copy->size());
		int nNodes = nSegments - copy->Nodes().size() + 1; // n nodes to add
		if (nNodes < 0) nNodes = 0;

		if (nNodes)
		{
			mf::LogVerbatim("ProjectionMatchingAlg") << "  add " << nNodes << " nodes";
			copy->Optimize(nNodes, fOptimizationEps);
		}
	}
	double g = copy->Optimize(0, fFineTuningEps);
	mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt done, g = " << g;

	return copy;
}
// ------------------------------------------------------


