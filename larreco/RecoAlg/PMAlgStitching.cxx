////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgStitching
// Author:      L.Whitehead (leigh.howard.whitehead@cern.ch) January 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgStitching.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"

// Constructor
//pma::PMAlgStitching::PMAlgStitching(pma::TrkCandidateColl &inputTracks, const pma::PMAlgStitching::Config &config):fInputTracks(inputTracks)
pma::PMAlgStitching::PMAlgStitching(pma::TrkCandidateColl &inputTracks):fInputTracks(inputTracks)
{

  GetTPCXOffsets();

}

// Destructor
pma::PMAlgStitching::~PMAlgStitching(){


}

// Main function of the algorithm
void pma::PMAlgStitching::StitchTracks(){

  // Loop over the track collection
  for(unsigned int t = 0; t < fInputTracks.size(); ++t){

    pma::Track3D* t1 = fInputTracks[t].Track();

    // Look through the following tracks for one to stitch
    pma::Track3D* bestTrkMatch = 0x0;
//    double bestScore = 0;
//    bool reverse = false;
//    bool flip1 = false;
//    bool flip2 = false;

    // Don't use the very end points of the tracks in case of scatter or distortion.
    TVector3 trk1Front = t1->Nodes()[2]->Point3D();
    TVector3 trk1Back = t1->Nodes()[t1->Nodes().size()-3]->Point3D();
    TVector3 trk1FrontDir = (t1->Nodes()[1]->Point3D() - trk1Front).Unit();
    TVector3 trk1BackDir = (t1->Nodes()[t1->Nodes().size()-2]->Point3D() - trk1Back).Unit();

    // For stitching, we only need the end nearest the stitching surface.
    double offset1 = GetTPCOffset(t1->FrontTPC(),t1->FrontCryo(),true);
    TVector3 t1Pos, t1Dir;
    bool isFront1 = false;
    GetBestPosAndDir(trk1Front,trk1FrontDir,trk1Back,trk1BackDir,t1Pos,t1Dir,offset1,isFront1);
    // Shift the best position in x to move the end point to the stitching surface
    double xShift1 = 0;
    if(isFront1){
      xShift1 = t1->Nodes()[0]->Point3D().X() - offset1;
    }
    else{
      xShift1 = t1->Nodes()[t1->Nodes().size()-1]->Point3D().X() - offset1;
    }
    t1Pos.X() -= xShift1;
    
    for(unsigned int u = t+1; u < fInputTracks.size(); ++u){

      pma::Track3D* t2 = fInputTracks[u].Track();

      // Don't use the very end points of the tracks in case of scatter or distortion.
      TVector3 trk2Front = t2->Nodes()[2]->Point3D();
      TVector3 trk2Back = t2->Nodes()[t1->Nodes().size()-3]->Point3D();
      TVector3 trk2FrontDir = (t2->Nodes()[1]->Point3D() - trk2Front).Unit();
      TVector3 trk2BackDir = (t2->Nodes()[t1->Nodes().size()-2]->Point3D() - trk2Back).Unit();

      // For stitching, we only need the end nearest the stitching surface.
      double offset2 = GetTPCOffset(t2->FrontTPC(),t2->FrontCryo(),true);
      TVector3 t2Pos, t2Dir;
      bool isFront2 = false;
      GetBestPosAndDir(trk2Front,trk2FrontDir,trk2Back,trk2BackDir,t2Pos,t2Dir,offset2,isFront2);
      // Shift the best position in x to move the end point to the stitching surface
      double xShift2 = 0;
      if(isFront2){
        xShift2 = t2->Nodes()[0]->Point3D().X() - offset2;
      }
      else{
        xShift2 = t2->Nodes()[t2->Nodes().size()-1]->Point3D().X() - offset2;
      }
      t2Pos.X() -= xShift2;

      double score = 0;
      score = GetTrackPairDelta(t1Pos,t2Pos,t1Dir,t2Dir,offset1);
      std::cout << score << std::endl;
    }

    // If we found a match, do something about it.
    if(bestTrkMatch != 0x0){


    }

  }

}

void pma::PMAlgStitching::GetBestPosAndDir(TVector3 &pos1, TVector3 &dir1, TVector3 &pos2, TVector3 &dir2, TVector3 &bestPos, TVector3 &bestDir, double offset, bool& isFront){ 

  if(fabs(pos1.X()-offset) < fabs(pos2.X()-offset)){
    bestPos = pos1;
    bestDir = dir1;
    isFront = true;
  }
  else{
    bestPos = pos2;
    bestDir = dir2;
    isFront = false;
  }

}

double pma::PMAlgStitching::GetTrackPairDelta(TVector3 &pos1, TVector3 &pos2, TVector3 &dir1, TVector3 &dir2, double mergePointX){

  double delta = -999.;

  // Calculate number of steps to reach the merge point in x.
  double steps1 = (pos1.X() - mergePointX) / dir1.X();
  double steps2 = (pos2.X() - mergePointX) / dir2.X();

  // Vectors at the merge point:
  TVector3 trk1Merge = pos1 + steps1*dir1;
  TVector3 trk2Merge = pos2 + steps2*dir2;  

  delta = (trk1Merge-trk2Merge).Mag();

  // Also need a way to punish the directions disagreeing.
  double angle = (180.0/TMath::Pi())*dir1.Angle(dir2);

  if(angle > 10) delta += 100;

  return delta;

}

void pma::PMAlgStitching::GetTPCXOffsets(){

  // Grab hold of the geometry
  auto const* geom = lar::providerFrom<geo::Geometry>();

  // Loop over each TPC and store the require information
  for (geo::TPCID const& tID: geom->IterateTPCIDs()) {

    geo::TPCGeo const& aTPC = geom->TPC(tID);

    // Loop over the 3 possible readout planes to find the x position
    unsigned int plane = 0;
    while ((plane < 3) && aTPC.HasPlane(plane)) ++plane;
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

double pma::PMAlgStitching::GetTPCOffset(unsigned int tpc, unsigned int cryo, bool isCPA){

  geo::TPCID thisTPCID(tpc,cryo);
  double offset = 0.0;
  if(isCPA){
    offset = fTPCXOffsetsCPA[thisTPCID];
  }
  else{
    offset = fTPCXOffsetsAPA[thisTPCID];
  }
  return offset;

}

