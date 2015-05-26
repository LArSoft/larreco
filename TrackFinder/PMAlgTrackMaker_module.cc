////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Makes 3D tracks with the Projection Matching Algorithm.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoAlg/PMAlg/PmaTrack3D.h"

#include <memory>

namespace trkf {

class PMAlgTrackMaker : public art::EDProducer {
public:
  explicit PMAlgTrackMaker(fhicl::ParameterSet const & p);

  PMAlgTrackMaker(PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker(PMAlgTrackMaker &&) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker &&) = delete;

  void reconfigure(fhicl::ParameterSet const& p);

  void produce(art::Event & e) override;

private:
  // initial functionality, only to make the PMA running - to be removed
  int maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindMany< recob::Hit >& fb,
	geo::View_t view, int* tpc);
  // --------------------------------------------------- - to be removed


  std::string fCluModuleLabel; // label for input cluster collection

};


PMAlgTrackMaker::PMAlgTrackMaker(fhicl::ParameterSet const & p)
{
	this->reconfigure(p);
	produces< std::vector<recob::Track> >();
}

void PMAlgTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
{
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
}

void PMAlgTrackMaker::produce(art::Event & evt)
{
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		art::FindMany< recob::Hit > fb(cluListHandle, evt, fCluModuleLabel);

		int max_clu_tpc = -1;
		int max_clu_idx = maxCluster(cluListHandle, fb, geo::kZ, &max_clu_tpc);
		if (max_clu_idx >= 0)
		{
			//recob::Cluster const& clu = (*cluListHandle)[max_clu_idx];

			std::vector< recob::Hit const* > v_coll, v_ind2, v_ind1;
			fb.get(max_clu_idx, v_coll);

			max_clu_idx = maxCluster(cluListHandle, fb, geo::kV, &max_clu_tpc);
			if (max_clu_idx >= 0) fb.get(max_clu_idx, v_ind2);

			max_clu_idx = maxCluster(cluListHandle, fb, geo::kU, &max_clu_tpc);
			if (max_clu_idx >= 0) fb.get(max_clu_idx, v_ind1);

			pma::Track3D trk;
			if (v_coll.size() > 10) { trk.AddHits(v_coll); std::cout << "add coll:" << v_coll.size() << std::endl; }
			if (v_ind2.size() > 10) { trk.AddHits(v_ind2); std::cout << "add ind2:" << v_ind2.size() << std::endl; }
			if (v_ind1.size() > 10) { trk.AddHits(v_ind1); std::cout << "add ind1:" << v_ind1.size() << std::endl; }

			std::cout << "*************** track size: " << trk.size() << std::endl;
			std::vector< int > tpcs = trk.TPCs();
			for (size_t t = 0; t < tpcs.size(); ++t)
			{
				std::cout << "  tpc:" << tpcs[t] << std::endl;
			}
			std::cout << "  #coll:" << trk.NHits(geo::kZ)
				<< " #ind2:" << trk.NHits(geo::kV)
				<< " #ind1:" << trk.NHits(geo::kU)
				<< std::endl;

			trk.Initialize();

		}
		else std::cout << "*************** small clusters" << std::endl;
	}
	else std::cout << "*************** no clusters" << std::endl;
}

// initial functionality, only to make the PMA running - to be removed
int PMAlgTrackMaker::maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindMany< recob::Hit >& fb,
	geo::View_t view, int* tpc)
{
	int idx = -1, tpc_result = -1;
	size_t min_clu_size = 10, s_max = 0, s;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		std::vector< recob::Hit const* > v;
		fb.get(i, v);

		if ((v.front()->WireID().Plane == view) &&
		    ((*tpc == -1) || ((int)(v.front()->WireID().TPC) == *tpc)))
		{
			s = v.size();
			if ((s > min_clu_size) && (s > s_max))
			{
				tpc_result = v.front()->WireID().TPC;
				s_max = s; idx = i;
			}
		}
	}
	if (*tpc == -1) *tpc = tpc_result;
	return idx;
}
// --------------------------------------------------- - to be removed

DEFINE_ART_MODULE(PMAlgTrackMaker)

} // namespace trkf
