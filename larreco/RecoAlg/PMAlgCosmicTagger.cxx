////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgCosmicTagger
// Author:      ..., R.Sulej (..., robert.sulej@cern.ch) March 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgCosmicTagger.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

void pma::PMAlgCosmicTagger::tag(pma::TrkCandidateColl& tracks)
{
    mf::LogInfo("pma::PMAlgCosmicTagger") << "Passed " << tracks.size() << " tracks for tagging cosmics.";

    size_t n = 0;

    if (fTagOutOfDriftTracks) n += outOfDriftWindow(tracks);
    if (fTagTopDownTracks) n += topDownCrossing(tracks);
		if (fTagNonBeamT0Tracks) n += nonBeamT0Tag(tracks);

    mf::LogInfo("pma::PMAlgCosmicTagger") << "...tagged " << n << " cosmic-like tracks.";
}


size_t pma::PMAlgCosmicTagger::outOfDriftWindow(pma::TrkCandidateColl& tracks)
{
    mf::LogInfo("pma::PMAlgCosmicTagger") << "   - tag tracks out of 1 drift window;";
    size_t n = 0;

    auto const* geom = lar::providerFrom<geo::Geometry>();

    for (auto & t : tracks.tracks())
    {
		
        // If this track is already tagged then don't try again!
		    if(t.Track()->GetTag() == pma::Track3D::kCosmic) continue;

        pma::Track3D & trk = *(t.Track());

        double min, max, p;
        bool node_out_of_drift = false;
        for (size_t nidx = 0; nidx < trk.Nodes().size(); ++nidx)
        {
            auto const & node = *(trk.Nodes()[nidx]);
            auto const & tpcGeo = geom->TPC(node.TPC(), node.Cryo());
            // DetectDriftDirection returns a short int, but switch requires an int
            int driftDir = abs(tpcGeo.DetectDriftDirection());
            p = node.Point3D()[driftDir-1];
            switch (driftDir)
            {
                case 1: min = tpcGeo.MinX(); max = tpcGeo.MaxX(); break;
                case 2: min = tpcGeo.MinY(); max = tpcGeo.MaxY(); break;
                case 3: min = tpcGeo.MinZ(); max = tpcGeo.MaxZ(); break;
                default: throw cet::exception("PMAlgCosmicTagger") << "Drift direction unknown: " << driftDir << std::endl;
            }
            if ((p < min - fOutOfDriftMargin) || (p > max + fOutOfDriftMargin)) { node_out_of_drift = true; break; }
        }

        if (node_out_of_drift) { trk.SetTagFlag(pma::Track3D::kCosmic); ++n; }
    }

		mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " tracks out of 1 drift window.";

    return n;
}

size_t pma::PMAlgCosmicTagger::topDownCrossing(pma::TrkCandidateColl& tracks)
{
	mf::LogInfo("pma::PMAlgCosmicTagger") << "   - tag tracks crossing full Y range.";
	size_t n = 0;

	auto const* geom = lar::providerFrom<geo::Geometry>();

	// Need to find the minimum and maximum height values from the geometry.
	double minY = 1.e6;
	double maxY = -1.e6;

	// Since we can stack TPCs, we can't just use geom::TPCGeom::Height()
  for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
    geo::TPCGeo const& TPC = geom->TPC(tID);

    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    
    if( center[1] - tpcDim[1] < minY ) minY = center[1] - tpcDim[1];
    if( center[1] + tpcDim[1] > maxY ) maxY = center[1] + tpcDim[1];
  } // for all TPC

	double detHeight = maxY - minY;

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - y range of detector = " << minY << ", " << maxY;

	// Loop over the tracks
	for(auto & t : tracks.tracks()){

		// If this track is already tagged then don't try again!
		if(t.Track()->GetTag() == pma::Track3D::kCosmic) continue;

		// Get the first and last y-positions from the track.
		auto const & node0 = *(t.Track()->Nodes()[0]);
		auto const & node1 = *(t.Track()->Nodes()[t.Track()->Nodes().size()-1]);

		double trkHeight = fabs(node0.Point3D()[1]-node1.Point3D()[1]);
	
		if((detHeight - trkHeight) < fTopDownMargin){
			++n;
			t.Track()->SetTagFlag(pma::Track3D::kCosmic);
		}
	}

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " tracks crossing the full Y range.";
	return n;
}

// Leigh: Make use of the fact that our cathode and anode crossing tracks have a reconstructed T0.
// Check to see if this time is consistent with the beam
size_t pma::PMAlgCosmicTagger::nonBeamT0Tag(pma::TrkCandidateColl &tracks){

	size_t n = 0;

	auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

	// Search through all of the tracks
	for(auto & t : tracks.tracks()){
		
		// If this track is already tagged then don't try again!
		if(t.Track()->GetTag() == pma::Track3D::kCosmic) continue;

		// Non zero T0 means we reconstructed it
		if(t.Track()->GetT0() != 0.0){
		mf::LogInfo("pma::PMAlgCosmicTagger") << " - track with T0 = " << t.Track()->GetT0();

			if(fabs(t.Track()->GetT0() - detprop->TriggerOffset()) > fNonBeamT0Margin){
				++n;
				t.Track()->SetTagFlag(pma::Track3D::kCosmic);
			}

		}

	}

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " non-beam T0 tracks.";
	return n;

}


