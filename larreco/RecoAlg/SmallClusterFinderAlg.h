////////////////////////////////////////////////////////////////////////
//
// \file SmallClusterFinderAlg.h
//
// \author corey.adams@yale.edu
//
// This algorithm is designed to find small clusters that could correspond to gammas
// or low energy electrons.
//
/*	There are two parameters that matter from the fcl file:
                fNHitsInClust is the number of hits that should be in these small clusters
                        ^-- Gamma events seem to rarely have more than 4 hits in the cluster
                        ^-- SN events are unclear.  Should this even be used for SN?
                fRadiusSizePar is the distance (in cm) between the small clusters and any other hits.

        This algorithm sorts the hits by plane, and then looks at each hit individually.  If
        there is a hit within RadiusSizePar, it gets added to a local list.  All other hits
        are ignored.  Then, if the number of hits that got added to the local list is greater
        then NHitsInClust, the original hit is ignored.  If it's less, the original hit is
        presumed to be part of a very small (or single hit) cluster.  So its added to the list
        of hits in the small cluster.

        All of the small clusters are then split apart into groups in the way you would expect.
        Each cluster is assigned an ID number to distinguish it, and the hits that aren't
        identified as small clusters all end up in the "leftover" cluster.  The numbering scheme
        is ID = 100*iPlane + Cluster on that plane, and the leftover hits are the first (0th)
        cluster written out.

        -Corey
*/
//
//
///////////////////////////////////////////////////////////////////////

#ifndef SMALLCLUSTERFINDERALG_H
#define SMALLCLUSTERFINDERALG_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/fwd.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace util {
  class GeometryUtilities;
}

namespace cluster {

  class SmallClusterFinderAlg {
  public:
    explicit SmallClusterFinderAlg(fhicl::ParameterSet const& pset);

    void SelectLocalHitlist(util::GeometryUtilities const& gser,
                            std::vector<art::Ptr<recob::Hit>> hitlist,
                            std::vector<art::Ptr<recob::Hit>>& hitlistlocal,
                            double wire_start,
                            double time_start,
                            double radlimit) const;
    void SelectLocalHitlist(util::GeometryUtilities const& gser,
                            std::vector<art::Ptr<recob::Hit>> hitlist,
                            std::vector<art::Ptr<recob::Hit>>& hitlistlocal,
                            double wire_start,
                            double time_start,
                            double radlimit,
                            std::vector<int>&) const;

    // pass in a hit, all the others are filled in the function.
    int GetPlaneAndTPC(art::Ptr<recob::Hit> a,
                       unsigned int& plane,
                       unsigned int& cryostat,
                       unsigned int& time,
                       unsigned int& wire) const;

    //passes in a hitlist, returns a list of gammas and puts everything else in hitlistleftover
    std::vector<art::Ptr<recob::Hit>> CreateHighHitlist(
      util::GeometryUtilities const& gser,
      std::vector<art::Ptr<recob::Hit>> const& hitlist,
      std::vector<art::Ptr<recob::Hit>>& hitlistleftover) const;

    // this routine actually finds the clusters
    // It stores them in smallClustList
    // It puts the leftover hits in hitlistleftover
    void FindSmallClusters(util::GeometryUtilities const& gser,
                           detinfo::DetectorClocksData const& dataClocks,
                           detinfo::DetectorPropertiesData const& detProp,
                           std::vector<art::Ptr<recob::Hit>> allHits);

    // This function is solely to sort the leftover hits.  It's not done
    // automatically because it's not a short process if the number of hits is
    // big.  For a few hits its not that bad. not yet functional!  by this I
    // mean, does not yet exist
    void sortLeftoverHits();

    std::vector<std::vector<art::Ptr<recob::Hit>>> GetSmallClustersByPlane(unsigned int iPlane);
    //This function just returns the clusters by plane, and they are split already.

    std::vector<art::Ptr<recob::Hit>> GetLeftoversByPlane(unsigned int iPlane);
    //This function returns the leftover hits.

  private:
    //Special function to make sure the vectors are all the right size
    //can't do the hits, but can clear them and do the planes.
    //Gets called each time FindSmallClusters is called, as well as in the constructor.
    void ClearandResizeVectors();

    //**************All of the constants needed to do this process:
    double fWiretoCm, fTimetoCm, fWireTimetoCmCm;
    float fTimeTick;                   // time sample in us
    float fDriftVelocity;              // used to convert time to cm
    double fWirePitch;                 // wire pitch in cm
    std::vector<unsigned int> fNWires; // Number of wires on each plane
    unsigned int fNPlanes;             // number of planes
    double fNTimes;

    // The various art grabbed Larsoft stuff to get physical properties.
    art::ServiceHandle<geo::Geometry const> geom;

    // The vectors that actually get used to hold hits, sets of hits, etc.
    // All have the plane as the first index
    std::vector<std::vector<std::vector<art::Ptr<recob::Hit>>>> smallClustList;
    // This list is meant to contain the small clusters.
    // It should be indexed by smallClustList[plane][smallcluster][hit]
    std::vector<std::vector<art::Ptr<recob::Hit>>>
      hitlistrefined; // list, by plane, of hits that ought to be gammas
    std::vector<std::vector<art::Ptr<recob::Hit>>> hitlistbyplane; // list of all hits on each plane
    std::vector<std::vector<art::Ptr<recob::Hit>>>
      hitlistleftover; // list of all hits on each plane that are NOT gammas

    // input parameters, grabbed from the parameter set
    bool verbose;          // Make lots of printouts?
    double fRadiusSizePar; // Determines the max radius of the cluster, must be
                           // separated
    double fNHitsInClust;  // Forces cluster to have a max number of hits
                           // Remember, this is the *small* cluster finder
  };                       // class SmallAngleFinder
}
#endif // SMALLCLUSTERFINDERALG_H
