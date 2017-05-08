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
		// Get the detector dimensions
		GetDimensions();

    mf::LogInfo("pma::PMAlgCosmicTagger") << "Passed " << tracks.size() << " tracks for tagging cosmics.";

    size_t n = 0;

    if (fTagOutOfDriftTracks) n += outOfDriftWindow(tracks);
    if (fTagFullHeightTracks) n += fullHeightCrossing(tracks);
    if (fTagFullWidthTracks) n += fullWidthCrossing(tracks);
    if (fTagFullLengthTracks) n += fullLengthCrossing(tracks);
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

size_t pma::PMAlgCosmicTagger::fullHeightCrossing(pma::TrkCandidateColl& tracks){

	// Just use the first tpc to get the height dimension (instead of assuming y).
	auto const* geom = lar::providerFrom<geo::Geometry>();
	TVector3 dir = geom->TPC(0,0).HeightDir();

	size_t n = fullCrossingTagger(tracks,ConvertDirToInt(dir));

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " tracks crossing the full detector height";
	return n;

}

size_t pma::PMAlgCosmicTagger::fullWidthCrossing(pma::TrkCandidateColl& tracks){

	// Just use the first tpc to get the width dimension (instead of assuming x).
	auto const* geom = lar::providerFrom<geo::Geometry>();
	TVector3 dir = geom->TPC(0,0).WidthDir();

	size_t n = fullCrossingTagger(tracks,ConvertDirToInt(dir));

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " tracks crossing the full detector width";
	return n;

}

size_t pma::PMAlgCosmicTagger::fullLengthCrossing(pma::TrkCandidateColl& tracks){

	// Just use the first tpc to get the length dimension (instead of assuming z).
	auto const* geom = lar::providerFrom<geo::Geometry>();
	TVector3 dir = geom->TPC(0,0).LengthDir();

	size_t n = fullCrossingTagger(tracks,ConvertDirToInt(dir));

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " tracks crossing the full detector length";
	return n;

}

size_t pma::PMAlgCosmicTagger::fullCrossingTagger(pma::TrkCandidateColl& tracks, int direction){

	if(direction == -1){
		mf::LogWarning("pma::PMAlgCosmicTagger") << " - Could not recognise direction, not attempting to perform fullCrossingTagger.";
	}

	size_t n = 0;

	double detDim = fDimensions[direction];

	// Loop over the tracks
	for(auto & t : tracks.tracks()){

		// If this track is already tagged then don't try again!
		if(t.Track()->GetTag() == pma::Track3D::kCosmic) continue;

		// Get the first and last y-positions from the track.
		auto const & node0 = *(t.Track()->Nodes()[0]);
		auto const & node1 = *(t.Track()->Nodes()[t.Track()->Nodes().size()-1]);

		// Get the length of the track in the requested direction
		double trkDim = fabs(node0.Point3D()[direction]-node1.Point3D()[direction]);
	
		if((detDim - trkDim) < fFullCrossingMargin){
			++n;
			t.Track()->SetTagFlag(pma::Track3D::kCosmic);
			mf::LogInfo("pma::PMAlgCosmicTagger") << " -- track tagged in direction " << direction << " with " << trkDim << " (c.f. " << detDim << ")";
		}
	}

	return n;
}

void pma::PMAlgCosmicTagger::GetDimensions(){

	// Need to find the minimum and maximum height values from the geometry.
	double minX = 1.e6;
	double maxX = -1.e6;
	double minY = 1.e6;
	double maxY = -1.e6;
	double minZ = 1.e6;
	double maxZ = -1.e6;
	
	auto const* geom = lar::providerFrom<geo::Geometry>();

	// Since we can stack TPCs, we can't just use geom::TPCGeom::Height()
  for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
    geo::TPCGeo const& TPC = geom->TPC(tID);

    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    
    if( center[0] - tpcDim[0] < minX ) minX = center[0] - tpcDim[0];
    if( center[0] + tpcDim[0] > maxX ) maxX = center[0] + tpcDim[0];
    if( center[1] - tpcDim[1] < minY ) minY = center[1] - tpcDim[1];
    if( center[1] + tpcDim[1] > maxY ) maxY = center[1] + tpcDim[1];
    if( center[2] - tpcDim[2] < minZ ) minZ = center[2] - tpcDim[2];
    if( center[2] + tpcDim[2] > maxZ ) maxZ = center[2] + tpcDim[2];
  } // for all TPC

	fDimensions.clear();
	fDimensions.push_back(maxX - minX);
	fDimensions.push_back(maxY - minY);
	fDimensions.push_back(maxZ - minZ);

}

int pma::PMAlgCosmicTagger::ConvertDirToInt(TVector3 &dir){

	if(dir.X() > 0.99) return 0; 
	if(dir.Y() > 0.99) return 1; 
	if(dir.Z() > 0.99) return 2; 

	else return -1;
}
