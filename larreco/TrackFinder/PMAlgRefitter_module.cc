////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgRefitter
// Module Type: producer
// File:        PMAlgRefitter_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), Oct. 2016
//
// Testing ideas to use PMA to re-fit 3D trajectories. Started for ArgoNeuT - MINOS matching problem
// where non-uniformities in E field result with spuriously curved reconstructed trajectory.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "canvas/Utilities/InputTag.h"
//#include "lardata/Utilities/PtrMaker.h"

#include "larreco/RecoAlg/PMAlgTracking.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include <memory>

namespace trkf {

class PMAlgRefitter : public art::EDProducer {
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<pma::ProjectionMatchingAlg::Config> ProjectionMatchingAlg {
			Name("ProjectionMatchingAlg")
		};

		fhicl::Table<pma::PMAlgRefitter::Config> PMAlgRefitting {
			Name("PMAlgRefitting")
		};

		fhicl::Atom<art::InputTag> HitModuleLabel {
			Name("HitModuleLabel"),
			Comment("tag of unclustered hits, which were used to produce PFPs and tracks")
		};

		fhicl::Atom<art::InputTag> PfpModuleLabel {
			Name("PfpModuleLabel"),
			Comment("tag of the input PFParticles and associated clusters")
		};
    };
    using Parameters = art::EDProducer::Table<Config>;

	explicit PMAlgRefitter(Parameters const& config);

	PMAlgRefitter(PMAlgRefitter const &) = delete;
	PMAlgRefitter(PMAlgRefitter &&) = delete;
	PMAlgRefitter & operator = (PMAlgRefitter const &) = delete;
	PMAlgRefitter & operator = (PMAlgRefitter &&) = delete;

	void produce(art::Event & e) override;

private:
  // ******************** fcl parameters ***********************
  art::InputTag fHitModuleLabel; // tag for unclustered hit collection
  art::InputTag fPfpModuleLabel; // tag for PFParticle and cluster collections

  pma::ProjectionMatchingAlg::Config fPmaConfig;
  pma::PMAlgRefitter::Config fPmaRefitterConfig;

  // *********************** geometry **************************
  art::ServiceHandle< geo::Geometry > fGeom;
};

PMAlgRefitter::PMAlgRefitter(PMAlgRefitter::Parameters const& config) :
    fHitModuleLabel(config().HitModuleLabel()),
	fPfpModuleLabel(config().PfpModuleLabel()),

	fPmaConfig(config().ProjectionMatchingAlg()),
	fPmaRefitterConfig(config().PMAlgRefitting())
{
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();

	produces< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();

	produces< art::Assns<recob::PFParticle, recob::Track> >();
}
// ------------------------------------------------------

void PMAlgRefitter::produce(art::Event& evt)
{
	// ---------------- Create data products ------------------
	auto tracks = std::make_unique< std::vector<recob::Track> >();
	auto allsp = std::make_unique< std::vector<recob::SpacePoint> >();

	auto trk2hit_oldway = std::make_unique< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	auto trk2hit = std::make_unique< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	auto trk2sp = std::make_unique< art::Assns<recob::Track, recob::SpacePoint> >();
	auto sp2hit = std::make_unique< art::Assns<recob::SpacePoint, recob::Hit> >();

	auto pfp2trk = std::make_unique< art::Assns< recob::PFParticle, recob::Track> >();

	// ------------------- Collect inputs ---------------------
	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
	art::Handle< std::vector<recob::Track> > trkHandle;
	if (!(evt.getByLabel(fHitModuleLabel, allHitListHandle) &&
	      evt.getByLabel(fPfpModuleLabel, pfparticleHandle) &&
	      evt.getByLabel(fPfpModuleLabel, trkHandle)))
	{
		mf::LogError("PMAlgRefitter") << "Inputs not found in the event.";
		return;
	}

    std::vector< art::Ptr<recob::Hit> > allhitlist;
    art::fill_ptr_vector(allhitlist, allHitListHandle);

	art::FindManyP< recob::Track > tracksFromPfps(pfparticleHandle, evt, fPfpModuleLabel);
	art::FindManyP< recob::Vertex > vtxFromPfps(pfparticleHandle, evt, fPfpModuleLabel);
	art::FindManyP< recob::Hit > hitsFromTracks(trkHandle, evt, fPfpModuleLabel);

	// -------------- PMA Refitter for this event ---------------
	auto pmalgRefitter = pma::PMAlgRefitter(
	    allhitlist, *pfparticleHandle,
		tracksFromPfps, vtxFromPfps, hitsFromTracks,
		fPmaConfig, fPmaRefitterConfig);

	// ------------------ Do the job here: --------------------
	int retCode = pmalgRefitter.build();
	// --------------------------------------------------------
	switch (retCode)
	{
		case -2: mf::LogError("Summary") << "problem"; break;
		case  0: mf::LogVerbatim("Summary") << "no tracks done"; break;
		default:
			if (retCode < 0) mf::LogVerbatim("Summary") << "unknown result";
			else if (retCode == 1) mf::LogVerbatim("Summary") << retCode << " track ready";
			else mf::LogVerbatim("Summary") << retCode << " tracks ready";
			break;
	}

	// ---------- Translate output to data products: ----------
	auto const & result = pmalgRefitter.result();

		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6];
		for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

		// use the following to create PFParticle <--> Track associations;
		std::map< size_t, std::vector< art::Ptr<recob::Track> > > pfPartToTrackVecMap;

		//auto const make_trkptr = lar::PtrMaker<recob::Track>(evt, *this); // PtrMaker Step #1

		tracks->reserve(result.size());
		for (size_t trkIndex = 0; trkIndex < result.size(); ++trkIndex)
		{
			pma::Track3D* trk = result[trkIndex].Track();

			trk->SelectHits();  // just in case, set all to enabled
			unsigned int itpc = trk->FrontTPC(), icryo = trk->FrontCryo();
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kU)) trk->CompleteMissingWires(geo::kU);
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kV)) trk->CompleteMissingWires(geo::kV);
			if (fGeom->TPC(itpc, icryo).HasPlane(geo::kZ)) trk->CompleteMissingWires(geo::kZ);

			tracks->push_back(pma::convertFrom(*trk, trkIndex));

			//auto const trkPtr = make_trkptr(tracks->size() - 1); // PtrMaker Step #2

			size_t trkIdx = tracks->size() - 1; // stuff for assns:
			art::ProductID trkId = getProductID< std::vector<recob::Track> >(evt);
			art::Ptr<recob::Track> trkPtr(trkId, trkIdx, evt.productGetter(trkId));

			// which idx from start, except disabled, really....
			unsigned int hIdxs[trk->size()];
			for (size_t h = 0, cnt = 0; h < trk->size(); h++)
			{
				if ((*trk)[h]->IsEnabled()) hIdxs[h] = cnt++;
				else hIdxs[h] = 0;
			}

			art::PtrVector< recob::Hit > sp_hits;
			spStart = allsp->size();
			for (int h = trk->size() - 1; h >= 0; h--)
			{
				pma::Hit3D* h3d = (*trk)[h];
				if (!h3d->IsEnabled()) continue;

				recob::TrackHitMeta metadata(hIdxs[h], h3d->Dx());
				trk2hit->addSingle(trkPtr, h3d->Hit2DPtr(), metadata);
				trk2hit_oldway->addSingle(trkPtr, h3d->Hit2DPtr()); // ****** REMEMBER to remove when FindMany improved ******

				double hx = h3d->Point3D().X();
				double hy = h3d->Point3D().Y();
				double hz = h3d->Point3D().Z();

				if ((h == 0) || (sp_pos[0] != hx) || (sp_pos[1] != hy) || (sp_pos[2] != hz))
				{
					if (sp_hits.size()) // hits assigned to the previous sp
					{
						util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
						sp_hits.clear();
					}
					sp_pos[0] = hx; sp_pos[1] = hy; sp_pos[2] = hz;
					allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
				}
				sp_hits.push_back(h3d->Hit2DPtr());
			}

			if (sp_hits.size()) // hits assigned to the last sp
			{
				util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
			}
			spEnd = allsp->size();

			if (spEnd > spStart) util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);

			// if there is a PFParticle collection then recover PFParticle and add info to map
			if (result[trkIndex].Key() > -1)
			{
				size_t trackIdx = tracks->size() - 1;
				art::ProductID trackId = getProductID< std::vector<recob::Track> >(evt);
				art::Ptr<recob::Track> trackPtr(trackId, trackIdx, evt.productGetter(trackId));
				pfPartToTrackVecMap[result[trkIndex].Key()].push_back(trackPtr);
			}
		}

    for (size_t pidx = 0; pidx < pfparticleHandle->size(); ++pidx)
    {
        art::Ptr<recob::PFParticle> pfParticle(pfparticleHandle, pidx);
        if (pfParticle.isNull())
        {
            mf::LogError("PMAlgRefitter") << "Error in PFParticle table, index: " << pidx;
            continue;
		}

        const auto piter = pfPartToTrackVecMap.find(pidx);
        if (piter != pfPartToTrackVecMap.end()) // track was refitted
        {
			util::CreateAssn(*this, evt, pfParticle, (*piter).second, *pfp2trk);
        }
        else // otherwise use the old tracks
        {
            util::CreateAssn(*this, evt, pfParticle, tracksFromPfps.at(pidx), *pfp2trk);
        }
	}

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));

	evt.put(std::move(trk2hit_oldway)); // ****** REMEMBER to remove when FindMany improved ******
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(sp2hit));

	evt.put(std::move(pfp2trk));
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

DEFINE_ART_MODULE(PMAlgRefitter)

} // namespace trkf

