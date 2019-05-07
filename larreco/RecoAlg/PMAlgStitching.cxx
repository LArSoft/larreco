////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgStitching
// Author:      L.Whitehead (leigh.howard.whitehead@cern.ch) January 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgStitching.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TVector3.h"

// Constructor
pma::PMAlgStitching::PMAlgStitching(const pma::PMAlgStitching::Config &config)
{

  // Set parameters from the config.
  fStitchingThreshold = config.StitchingThreshold();
  fNodesFromEnd = config.NodesFromEnd();

  // Get CPA and APA positions.
  GetTPCXOffsets();

}

// Destructor
pma::PMAlgStitching::~PMAlgStitching(){


}

// CPA stitching wrapper
void pma::PMAlgStitching::StitchTracksCPA(pma::TrkCandidateColl &tracks){
  mf::LogInfo("pma::PMAlgStitching") << "Passed " << tracks.size() << " tracks for CPA stitching.";
  StitchTracks(tracks,true);
}

// APA stitching wrapper
void pma::PMAlgStitching::StitchTracksAPA(pma::TrkCandidateColl &tracks){
  mf::LogInfo("pma::PMAlgStitching") << "Passed " << tracks.size() << " tracks for APA stitching.";
  StitchTracks(tracks,false);
}

// Main function of the algorithm
// isCPA = true  : attempt to stitch tracks across the cathode.
//       = false : attempt to stitch tracks across the anode.
void pma::PMAlgStitching::StitchTracks(pma::TrkCandidateColl &tracks, bool isCPA){

  unsigned int minTrkLength = 2*fNodesFromEnd + 3;
  // Special case for fNodesFromEnd = 0
  if(minTrkLength < 6) minTrkLength = 6;

  // Loop over the track collection
  unsigned int t = 0;
  while(t < tracks.size()){

    pma::Track3D* t1 = tracks[t].Track();
    if(t1->Nodes().size() < minTrkLength) { ++t; continue; }

    // Look through the following tracks for one to stitch
    pma::Track3D* bestTrkMatch = 0x0;

    // Don't use the very end points of the tracks in case of scatter or distortion.
    TVector3 trk1Front = t1->Nodes()[(fNodesFromEnd)]->Point3D();
    TVector3 trk1Back = t1->Nodes()[t1->Nodes().size()-1-fNodesFromEnd]->Point3D();
    TVector3 trk1FrontDir = (trk1Front - t1->Nodes()[(fNodesFromEnd+1)]->Point3D()).Unit();
    TVector3 trk1BackDir = (trk1Back - t1->Nodes()[t1->Nodes().size()-1-(fNodesFromEnd+1)]->Point3D()).Unit();

    // For stitching, we need to consider both ends of the track.
    double offsetFront1 = GetTPCOffset(t1->FrontTPC(),t1->FrontCryo(),isCPA);
    double offsetBack1 = GetTPCOffset(t1->BackTPC(),t1->BackCryo(),isCPA);

    bool isBestFront1 = false;
    bool isBestFront2 = false;
    double xBestShift = 0;
    double frontShift1 = t1->Nodes()[0]->Point3D().X() - offsetFront1;
    double backShift1 = t1->Nodes()[t1->Nodes().size()-1]->Point3D().X() - offsetBack1;

    double bestMatchScore = 99999;

    for(unsigned int u = t+1; u < tracks.size(); ++u){

      pma::Track3D* t2 = tracks[u].Track();
      if(t2->Nodes().size() < minTrkLength) continue;

      // Don't use the very end points of the tracks in case of scatter or distortion.
      TVector3 trk2Front = t2->Nodes()[(fNodesFromEnd)]->Point3D();
      TVector3 trk2Back = t2->Nodes()[t2->Nodes().size()-1-fNodesFromEnd]->Point3D();
      TVector3 trk2FrontDir = (trk2Front - t2->Nodes()[(fNodesFromEnd+1)]->Point3D()).Unit();
      TVector3 trk2BackDir = (trk2Back - t2->Nodes()[t2->Nodes().size()-1-(fNodesFromEnd+1)]->Point3D()).Unit();

      // For stitching, we need to consider both ends of the track.
      double offsetFront2 = GetTPCOffset(t2->FrontTPC(),t2->FrontCryo(),isCPA);
      double offsetBack2 = GetTPCOffset(t2->BackTPC(),t2->BackCryo(),isCPA);

      // If the points to match are in the same TPC, then don't bother.
      // Remember we have 4 points to consider here.
      bool giveUp = false;
      geo::TPCID tpc1(0,0);
      geo::TPCID tpc2(0,0);
      // Front-to-front
      tpc1 = geo::TPCID(t1->FrontCryo(),t1->FrontTPC());
      tpc2 = geo::TPCID(t2->FrontCryo(),t2->FrontTPC());
      if(tpc1 == tpc2) giveUp = true;
      // Front-to-back
      tpc2 = geo::TPCID(t2->BackCryo(),t2->BackTPC());
      if(tpc1 == tpc2) giveUp = true;
      // Back-to-front
      tpc1 = geo::TPCID(t1->BackCryo(),t1->BackTPC());
      tpc2 = geo::TPCID(t2->FrontCryo(),t2->FrontTPC());
      if(tpc1 == tpc2) giveUp = true;
      // Back-to-back
      tpc2 = geo::TPCID(t2->BackCryo(),t2->BackTPC());
      if(tpc1 == tpc2) giveUp = true;

      // If the tracks have one end in the same TPC, give up.
      if(giveUp) continue;

      // Also check that these tpcs do meet at the stitching surface (not a problem for protoDUNE).
      double surfaceGap = 10.0;
      bool carryOn[4] = {true,true,true,true};
      if(fabs(offsetFront1 - offsetFront2) > surfaceGap) carryOn[0] = false;
      if(fabs(offsetFront1 - offsetBack2) > surfaceGap) carryOn[1] = false;
      if(fabs(offsetBack1 - offsetFront2) > surfaceGap) carryOn[2] = false;
      if(fabs(offsetBack1 - offsetBack2) > surfaceGap) carryOn[3] = false;

      // Loop over the four options
      for(int i = 0; i < 4; ++i){

        if(!carryOn[i]) continue;

        TVector3 t1Pos;
        TVector3 t2Pos;
        TVector3 t1Dir;
        TVector3 t2Dir;
        double xShift1;
        if(i < 2){
          t1Pos = trk1Front;
          t1Dir = trk1FrontDir;
          xShift1 = frontShift1;
        }
        else{
          t1Pos = trk1Back;
          t1Dir = trk1BackDir;
          xShift1 = backShift1;
        }
        if(i%2 == 0){
          t2Pos = trk2Front;
          t2Dir = trk2FrontDir;
        }
        else{
          t2Pos = trk2Back;
          t2Dir = trk2BackDir;
        }

        // Make sure the x directions point towards eachother (could be an issue for matching a short track)
        if(t1Dir.X() * t2Dir.X() > 0){
          continue;
        }
//        t1Pos.SetX(t1Pos.X() - xShift1);
//        t2Pos.SetX(t2Pos.X() + xShift1);

        double score = 0;
//        score = GetTrackPairDelta(t1Pos,t2Pos,t1Dir,t2Dir);
        score = GetOptimalStitchShift(t1Pos,t2Pos,t1Dir,t2Dir,xShift1);

        if(score < fStitchingThreshold && score < bestMatchScore){

          bestTrkMatch = t2;
          xBestShift = xShift1;
          bestMatchScore = score;
          if(i < 2){
            isBestFront1 = true;
          }
          else{
            isBestFront1 = false;
          }
          if(i % 2 == 0){
            isBestFront2 = true;
          }
          else{
            isBestFront2 = false;
          }
          mf::LogInfo("pma::PMAlgStitcher") << "Tracks " << t << " and " << u << " matching score = " << score << std::endl
            << " - " << t1Pos.X() << ", " << t1Pos.Y() << ", " << t1Pos.Z() << " :: " << t1Dir.X() << ", " << t1Dir.Y() << ", " << t1Dir.Z() << std::endl
            << " - " << t2Pos.X() << ", " << t2Pos.Y() << ", " << t2Pos.Z() << " :: " << t2Dir.X() << ", " << t2Dir.Y() << ", " << t2Dir.Z() << std::endl
            << " - " << t1->FrontCryo() << ", " << t1->FrontTPC() << " :: " << t1->BackCryo() << ", " << t1->BackTPC() << std::endl
            << " - " << t2->FrontCryo() << ", " << t2->FrontTPC() << " :: " << t2->BackCryo() << ", " << t2->BackTPC() << std::endl
            << " - " << isBestFront1 << " :: " << isBestFront2 << std::endl;
        } // End successful match if
      } // Loop over matching options
    } // Loop over track 2

    // If we found a match, do something about it.
    if(bestTrkMatch != 0x0){

      bool flip1 = false;
      bool flip2 = false;
      bool reverse = false;

      // Front-to-front match
      if(isBestFront1 && isBestFront2){
        flip1 = true;
      }
      // Front-to-back match
      else if(isBestFront1 && !isBestFront2){
        reverse = true;
      }
      // Back-to-back match
      else if(!isBestFront1 && !isBestFront2){
        flip2 = true;
      }
      // Back-to-front match (do nothing)

      int tid1 = tracks.getCandidateTreeId(t1);
      int tid2 = tracks.getCandidateTreeId(bestTrkMatch);

      bool canMerge = true;
      if ((tid1 < 0) || (tid2 < 0))
      {
        throw cet::exception("pma::PMAlgStitching") << "Track not found in the collection." << std::endl;
      }
      if (flip1){

        std::vector< pma::Track3D* > newTracks;
        if(t1->Flip(newTracks)){
          mf::LogInfo("pma::PMAlgStitching") << "Track 1 flipped.";
        }
        else{
          mf::LogInfo("pma::PMAlgStitching") << "Unable to flip Track 1.";
          canMerge = false;
        }
        for (const auto ts : newTracks){ // there may be a new track even if entire flip was not possible
          tracks.tracks().emplace_back(ts, -1, tid1);
        }
      }
      if (flip2){

        std::vector< pma::Track3D* > newTracks;
        if(bestTrkMatch->Flip(newTracks)){
          mf::LogInfo("pma::PMAlgStitching") << "Track 2 flipped.";
        }
        else{
          mf::LogInfo("pma::PMAlgStitching") << "Unable to flip Track 1.";
          canMerge = false;
        }
        for (const auto ts : newTracks){ // there may be a new track even if entire flip was not possible
          tracks.tracks().emplace_back(ts, -1, tid2);
        }
      }

      t1->GetRoot()->ApplyDriftShiftInTree(-xBestShift);
      bestTrkMatch->GetRoot()->ApplyDriftShiftInTree(+xBestShift);

      if (canMerge)
      {
        mf::LogInfo("pma::PMAlgStitching") << "Merging tracks...";
        int idx1 = tracks.getCandidateIndex(t1);
        int idx2 = tracks.getCandidateIndex(bestTrkMatch);
        if (reverse) // merge current track to another track, do not increase the outer loop index t (next after current track jumps in at t)
        {
          if (tracks.setTreeOriginAtFront(t1)) { tracks.merge((size_t)idx2, (size_t)idx1); }
          else { mf::LogWarning("pma::PMAlgStitching") << "   could not merge."; ++t; }
        }
        else // merge to the current track, do not increase the outer loop index t (maybe something else will match to the extended track)
        {
          if (tracks.setTreeOriginAtFront(bestTrkMatch)) { tracks.merge((size_t)idx1, (size_t)idx2); }
          else { mf::LogWarning("pma::PMAlgStitching") << "   could not merge."; ++t; }
        }
        mf::LogInfo("pma::PMAlgStitching") << "...done";
      }
      else { ++t; } // track matched, but not merged, go to the next track in the outer loop
    }
    else { ++t; } // no track matched, go to the next track in the outer loop
  }
}

// Perform the matching, allowing the shift to vary within +/- 5cm.
double pma::PMAlgStitching::GetOptimalStitchShift(TVector3 &pos1, TVector3 &pos2, TVector3 &dir1, TVector3 &dir2, double &shift){

  double stepSize = 0.1;
  double minShift = shift - (50. * stepSize);
  double maxShift = shift + (50. * stepSize);
  double bestShift = 99999;
  double bestScore = 99999;

  for(shift = minShift; shift <= maxShift; shift += stepSize){
    TVector3 newPos1 = pos1;
    TVector3 newPos2 = pos2;
    newPos1.SetX(pos1.X() - shift);
    newPos2.SetX(pos2.X() + shift);
    double thisScore = GetTrackPairDelta(newPos1,newPos2,dir1,dir2);
    if(thisScore < bestScore){
      bestShift = shift;
      bestScore = thisScore;
    }
  }

  shift = bestShift;
  return bestScore;
}

// Perform the extrapolation between the two vectors and return the distance between them.
double pma::PMAlgStitching::GetTrackPairDelta(TVector3 &pos1, TVector3 &pos2, TVector3 &dir1, TVector3 &dir2){

  double delta = -999.;

  // Calculate number of steps to reach the merge point in x.
  double steps1 = (pos2.X() - pos1.X()) / dir1.X();
  double steps2 = (pos1.X() - pos2.X()) / dir2.X();

  // Extrapolate each vector to the other's extrapolation point
  TVector3 trk1Merge = pos1 + steps1*dir1;
  TVector3 trk2Merge = pos2 + steps2*dir2;

  // Find the difference between each vector and the extrapolation
  delta = (trk1Merge-pos2).Mag() + (trk2Merge-pos1).Mag();

  return delta;

}

// Get the CPA and APA positions from the geometry
void pma::PMAlgStitching::GetTPCXOffsets(){

  // Grab hold of the geometry
  auto const* geom = lar::providerFrom<geo::Geometry>();

  // Loop over each TPC and store the require information
  for (geo::TPCID const& tID: geom->IterateTPCIDs()) {

    geo::TPCGeo const& aTPC = geom->TPC(tID);

    // Loop over the 3 possible readout planes to find the x position
    unsigned int plane = 0;
    bool hasPlane = false;
    for(;plane < 4; ++plane){
      hasPlane = aTPC.HasPlane(plane);
      if(hasPlane){
        break;
      }
    }

    if(!hasPlane){
      continue;
    }

    // Get the x position of the readout plane
    double xAnode = aTPC.PlaneLocation(plane)[0];
    fTPCXOffsetsAPA.insert(std::make_pair(tID,xAnode));

    // For the cathode, we have to try a little harder. Firstly, find the
    // min and max x values for the TPC.
    double origin[3] = {0.};
    double center[3] = {0.};
    aTPC.LocalToWorld(origin, center);
    double tpcDim[3] = {aTPC.HalfWidth(), aTPC.HalfHeight(), 0.5*aTPC.Length() };
    double xmin = center[0] - tpcDim[0];
    double xmax  = center[0] + tpcDim[0];
    double xCathode = 0.;
    // Now check which is further from the APA and use that as the cathode value.
    if(fabs(xmin - xAnode) > fabs(xmax-xAnode)){
      xCathode = xmin;
    }
    else{
      xCathode = xmax;
    }
    fTPCXOffsetsCPA.insert(std::make_pair(tID,xCathode));
  }
}

// Interface to get the CPA and APA positions from the maps in which they are stored.
double pma::PMAlgStitching::GetTPCOffset(unsigned int tpc, unsigned int cryo, bool isCPA){

  geo::TPCID thisTPCID(cryo,tpc);
  double offset = 0.0;
  if(isCPA){
    offset = fTPCXOffsetsCPA[thisTPCID];
  }
  else{
    offset = fTPCXOffsetsAPA[thisTPCID];
  }
  return offset;
}

