////////////////////////////////////////////////////////////////////////
///
/// \file   SpacePointAlg.h
///
/// \brief  Algorithm for generating space points from hits.
///
/// \author H. Greenlee 
///
/// This class calculates space points (recob::SpacePoint) from an 
/// unsorted collection of hits (recob::Hit).  The resulting space
/// points will contain one hit from from two or three views.
///
/// FCL parameters:
///
/// MaxDT - The maximum time difference (ticks) between any pair of hits.
/// MaxS  - The maximum 3-view wire separation parameter S (cm).
/// MinViews - Minimum number of views to make a space point (2 or 3).
/// EnableU - Use U view hits.
/// EnableV - Use V view hits.
/// EnableW - Use W view hits.
/// Filter - Filter space points flag.
/// Merge - Merge space points flag.
/// PreferColl - Collection view will be used for filtering and merging, and
///              space points will be sorted by collection wire.
///
/// The parameters fMaxDT and fMaxS are used to implement a notion of whether
/// the input hits are compatible with being a space point.  Parameter
/// MaxS is a cut on the 3-plane wire separation parameter S, which is
/// defined as follows:
///
/// S = sin(theta_vw)*u + sin(theta_wu)*v + sin(theta_uv)*w
///
/// where wire coordinates (u,v,w) are measured in cm with respect to a 
/// common origin.
///
/// The time offsets are subtracted from times embedded in hits before
/// comparing times in different planes, and before converting time
/// to distance.
///
/// If enabled, filtering eliminates multiple space points with similar
/// times on the same wire of the most popluated plane.
///
/// If enabled, merging combines multiple space points with similar 
/// times on the same wire of the most populated plane (potentially
/// producing space points with more hits than the number of planes).
///
/// There should eventually be a better way to specify time offsets.
///
////////////////////////////////////////////////////////////////////////

#ifndef SPACEPOINTALG_H
#define SPACEPOINTALG_H

#include <vector>
#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"

class TH1F;
namespace sim {
  class IDE;
}
namespace trkf{
  class KHitTrack;
}
namespace recob{
  class Hit;
  class SpacePoint;
}

namespace trkf {

  class SpacePointAlg {
  public:

    // Constructor.
    SpacePointAlg(const fhicl::ParameterSet& pset);

    // Destructor.
    ~SpacePointAlg();

    // Configuration Accessors.

    bool filter() const {return fFilter;}
    bool merge() const {return fMerge;}
    double maxDT() const {return fMaxDT;}
    double maxS() const {return fMaxS;}
    int minViews() const {return fMinViews;}
    bool enableU() const {return fEnableU;}
    bool enableV() const {return fEnableV;}
    bool enableW() const {return fEnableW;}

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    // Print constants obtained from geometry and properties services.
    void update() const;

    // Corrected time accessors.
    double correctedTime(const recob::Hit& hit) const;

    // Spatial separation of hits (zero if two or fewer).
    double separation(const art::PtrVector<recob::Hit>& hits) const;

    // Test whether the specified hits are compatible with a space point.
    // The last two arguments can be used to override the default cuts.
    bool compatible(const art::PtrVector<recob::Hit>& hits,
		    bool useMC = false) const;

    // Fill a single simple space point using the specified hits.
    // Hits are assumed to be compatible.
    void fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
			std::vector<recob::SpacePoint>& sptv,
			int sptid) const;

    /// Fill a collection of space points.
    void fillSpacePoints(std::vector<recob::SpacePoint>& spts,
			 std::multimap<double, KHitTrack> const& trackMap) const;

    // Fill a single complex space point using the specified hits.
    // Complex space points allow multiple hits in one plane.
    // Hits are assumed to be compatible.
    void fillComplexSpacePoint(const art::PtrVector<recob::Hit>& hits,
			       std::vector<recob::SpacePoint> &sptv,
			       int sptid) const;

    // Fill a vector of space points from an unsorted collection of hits.
    // Space points are generated for all compatible combinations of hits.
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts) const;

    // Fill a vector of space points compatible with mc truth information 
    void makeMCTruthSpacePoints(const art::PtrVector<recob::Hit>& hits,
				std::vector<recob::SpacePoint>& spts) const;

    // Get hits associated with a particular space point, based on most recent 
    // invocation of any make*SpacePoints method.
    const art::PtrVector<recob::Hit>& getAssociatedHits(const recob::SpacePoint& spt) const;

    // Clear space point to Hit associations.
    void clearHitMap() const {fSptHitMap.clear();}

    // Return number of space point to Hit associations.
    int numHitMap() const {return fSptHitMap.size();}

  private:

    // This is the real method for calculating space points (each of
    // the public make*SpacePoints methods comes here).
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts,
			 bool useMC) const;

    // Configuration paremeters.

    double fMaxDT;          ///< Maximum time difference between planes.
    double fMaxS;           ///< Maximum space separation between wires.
    int fMinViews;          ///< Mininum number of views per space point.
    bool fEnableU;          ///< Enable flag (U).
    bool fEnableV;          ///< Enable flag (V).
    bool fEnableW;          ///< Enable flag (W).
    bool fFilter;           ///< Filter flag.
    bool fMerge;            ///< Merge flag.
    bool fPreferColl;       ///< Sort by collection wire.

    // Temporary variables.

    struct HitMCInfo
    {
      std::vector<int> trackIDs;       ///< Parent trackIDs.
      std::vector<double> xyz;         ///< Location of ionization (all tracks).
      std::vector<const recob::Hit*> pchit;   ///< Pointer to nearest neighbor hit (indexed by plane).
      std::vector<double> dist2;              ///< Distance to nearest neighbor hit (indexed by plane).
    };
    mutable std::map<const recob::Hit*, HitMCInfo> fHitMCMap;
    mutable std::map<int, art::PtrVector<recob::Hit> > fSptHitMap;
  };
}

#endif
