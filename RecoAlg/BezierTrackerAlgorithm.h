#include "art/Persistency/Common/PtrVector.h"
#include "TVector3.h"

#ifndef BEZIERTRACKERALG_H
#define BEZIERTRACKERALG_H

//
// Name: BezierTrackerAlgorithm.h
//
// Purpose: Header file for module BezierTrackerAlgorithm.  This 
//  algorithm contains tools for producing and combining bezier
//  tracks made up of seed segments.
////
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SeedFinderAlgorithm.h"

namespace recob
{
  class Seed;
  class Track;
  class Hit;
  class Vertex;
}



namespace trkf {

  class BezierTrack;
 
  class BezierTrackerAlgorithm
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackerAlgorithm(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackerAlgorithm();


    
    void FilterOverlapTracks(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);

    void MakeDirectJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);
    
    void AddPtrVectors(art::PtrVector<recob::Hit>& Receiever, art::PtrVector<recob::Hit> const & ToAdd);
      
      
    
    trkf::SeedFinderAlgorithm * GetSeedFinderAlgorithm() { return fTheSeedFinder;}

    void  MakeVertexJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<recob::Vertex>& Vertices, std::vector<std::vector<int> >& Mapping);
    

    void MakeOverlapJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);


    std::vector<trkf::BezierTrack> MakeTracks(std::vector<std::vector<art::PtrVector<recob::Hit> > >& SortedHits, std::vector<art::PtrVector<recob::Hit> >& HitAssocs);
     
    void GetTracksForCombo(std::vector<recob::Seed>& Seeds, art::PtrVector<recob::Hit>& UHits, art::PtrVector<recob::Hit>& VHits, art::PtrVector<recob::Hit>& WHits);

    std::vector<std::vector< recob::Seed > > OrganizeSeedsIntoTracks(std::vector<recob::Seed >& AllSeeds, std::vector<art::PtrVector<recob::Hit> * >& AllHits, std::vector<art::PtrVector<recob::Hit> >& WhichHitsPerSeed, std::vector<std::vector< std::vector<int> >* >& OrgHits, std::vector<std::vector<std::vector<int> > >& WhichHitsPerTrack);

    void GetSeedDirProjected(recob::Seed const& TheSeed, std::vector<double>& WireCoord, std::vector<double>& TimeCoord);

    bool EvaluateOccupancy(recob::Seed& Seed1, recob::Seed& Seed2, double dThresh,  std::vector<art::PtrVector<recob::Hit>*>& AllHits,  std::vector<std::vector< std::vector<int> >* >& OrgHits,  std::vector<uint32_t>& LowChan, std::vector<uint32_t>& HighChan, std::vector<std::vector<int> >& HitStatus, std::vector<std::vector<int> >& TheseHits);

    void SortTracksByLength(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);
    
    void CalculateGeometricalElements();
    
    trkf::BezierTrack JoinTracks(trkf::BezierTrack& BT1, trkf::BezierTrack& BT2);


    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    
    

  private:

    // Fcl Attributes.

    

    
    double fOverlapCut;
    double fDirectJoinDistance;
    double fTrackJoinAngle;
    double fOccupancyThresh;
    double fTrackResolution;

    double fVertexImpactThreshold;
    double fVertexExtrapDistance;

    std::vector<double>   fPitches;
    std::vector<double>   fPlaneTimeOffsets;
    std::vector<TVector3> fPitchDir;
    std::vector<TVector3> fWireDir;
    std::vector<double>   fWireZeroOffset;
    TVector3              fXDir, fYDir, fZDir;


    void GetImpact(TVector3 t1pt, TVector3 t1dir, TVector3 t2pt, TVector3 t2dir, double& ImpactParam, double& Dist1, double& Dist2);

 


    SeedFinderAlgorithm * fTheSeedFinder;

  };
}

#endif
