#ifndef SEEDFINDERALG_H
#define SEEDFINDERALG_H

//
// Name: SeedFinderAlgorithm.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SpacePointAlg.h"
#include "TVector3.h"
#include "Geometry/Geometry.h"
#include "TTree.h"

namespace recob
{
  class SpacePoint;
  class Seed;
  class Hit;
}

namespace trkf {

  class SeedFinderAlgorithm
  {
  public:

    //--------------------------------------
    // Constructors, destructor, reconfigure
    //--------------------------------------

    SeedFinderAlgorithm(const fhicl::ParameterSet& pset);
   ~SeedFinderAlgorithm();

    void reconfigure(fhicl::ParameterSet const& pset);

 

    //----------------------
    // Seedfinding methods
    //----------------------


    std::vector<std::vector<recob::Seed> > GetSeedsFromSortedHits( std::vector< std::vector<art::PtrVector<recob::Hit> > >  const& SortedHits, 
								   std::vector<std::vector<art::PtrVector<recob::Hit> > >& HitsPerSeed, unsigned int StopAfter=0);
                                    // Return a vector of vectors of seeds, one vector for each supplied cluster 
                                    //   combination which has sufficient overlap. The second argument returns
                                    //   the hits sorted by combo and by seed
   
    

    std::vector<recob::Seed>    GetSeedsFromUnSortedHits(art::PtrVector<recob::Hit> const &, std::vector<art::PtrVector<recob::Hit> >&, unsigned int StopAfter=0);
                                    // Return a vector of seeds formed from an unstructured collection of hits    



    //----------------------
    // Alg passing 
    //----------------------


    SpacePointAlg *             GetSpacePointAlg() const { return fSptalg; }
                                    // Return the SpacePointAlg, as configured for the Seed Finding 
   
    


    
  private:

    //----------------------
    // Internal methods
    //----------------------


    std::vector<recob::Seed>    FindSeeds( art::PtrVector<recob::Hit> const& HitsFlat, std::vector<art::PtrVector<recob::Hit> >& CataloguedHits, unsigned int StopAfter);
                                    // Find a collection of seeds, based on the supplied set of hits.
                                    //  The second argument returns the hits catalogued by which
                                    //  seed they fell into (if any) 
   


    recob::Seed                 FindSeedAtEnd(std::vector<recob::SpacePoint> const&, std::vector<char>&, std::vector<int>&,
					      art::PtrVector<recob::Hit> const& HitsFlat, std::vector< std::vector< std::vector<int> > >& OrgHits);
                                    // Find one seed at high Z from the spacepoint collection given. Latter arguments are 
                                    //  for internal book keeping.



    size_t                      CountHits(std::vector<recob::SpacePoint> const& Points);
                                   // Counting the number of hits in each view which are associated with a set of SPs


  
    void                        GetCenterAndDirection(art::PtrVector<recob::Hit> const& HitsFlat, std::vector<int>& HitsToUse, TVector3& Center, TVector3& Direction, std::vector<double>& ViewRMS, std::vector<int>& HitsPerView);

  
    void                        ConsolidateSeed(recob::Seed& TheSeed, art::PtrVector<recob::Hit> const&, std::vector<char>& HitStatus,
						std::vector< std::vector< std::vector<int> > >& OrgHits, bool Extend);

    void                        GetHitDistAndProj( recob::Seed const& ASeed,  art::Ptr<recob::Hit> const& AHit, double& disp, double& s);
    

    void                        CalculateGeometricalElements();
                                   // Pre-calculating geometrical factors
  
                       
    // Fcl Attributes.

    SpacePointAlg       *  fSptalg;            
 
    double                fInitSeedLength;    
                                        
    int                   fMinPointsInSeed;   
                                        
    int                   fRefits;            
    
    std::vector<double>   fMaxViewRMS;

    float                 fHitResolution;

    float                 fOccupancyCut;
    
    double                fLengthCut;
    
    bool                  fExtendSeeds;


    std::vector<double>   fPitches;
    std::vector<TVector3> fPitchDir;
    std::vector<TVector3> fWireDir;
    std::vector<double>   fWireZeroOffset;
    TVector3              fXDir, fYDir, fZDir;
    size_t                fNChannels;
  };
  
}

#endif // SEEDFINDER_H
