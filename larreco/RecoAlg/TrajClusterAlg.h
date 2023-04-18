////////////////////////////////////////////////////////////////////////
//
//
// TrajClusterAlg
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALG_H
#define TRAJCLUSTERALG_H

// C/C++ standard libraries
#include <string>
#include <utility> // std::pair<>
#include <vector>

// framework libraries
namespace fhicl {
  class ParameterSet;
}

// LArSoft libraries
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larevt/CalibrationDBI/Interface/CalibrationDBIFwd.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
namespace detinfo {
  class DetectorClocksData;
}

// ROOT libraries
#include "TMVA/Reader.h"
class TTree;

namespace tca {

  class TrajClusterAlg {
  public:
    explicit TrajClusterAlg(fhicl::ParameterSet const& pset);

    bool SetInputHits(std::vector<recob::Hit> const& inputHits,
                      unsigned int run,
                      unsigned int event);
    void SetInputSpts(std::vector<recob::SpacePoint> const& sptHandle)
    {
      evt.sptHandle = &sptHandle;
    }
    void SetSourceHits(std::vector<recob::Hit> const& srcHits);
    void ExpectSlicedHits() { evt.expectSlicedHits = true; }
    void RunTrajClusterAlg(detinfo::DetectorClocksData const& clockData,
                           detinfo::DetectorPropertiesData const& detProp,
                           std::vector<unsigned int>& hitsInSlice,
                           int sliceID,
                           lariov::DBTimeStamp_t ts);
    bool CreateSlice(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     std::vector<unsigned int>& hitsInSlice,
                     int sliceID,
                     lariov::DBTimeStamp_t ts);
    void FinishEvent();

    void DefineShTree(TTree* t);

    unsigned short GetSlicesSize() const { return slices.size(); }
    TCSlice const& GetSlice(unsigned short sliceIndex) const { return slices[sliceIndex]; }
    void MergeTPHits(std::vector<unsigned int>& tpHits,
                     std::vector<recob::Hit>& newHitCol,
                     std::vector<unsigned int>& newHitAssns) const;

    std::vector<unsigned int> const& GetAlgModCount() const { return fAlgModCount; }
    std::vector<std::string> const& GetAlgBitNames() const { return AlgBitNames; }

    /// Deletes all the results
    void ClearResults()
    {
      slices.resize(0);
      evt.sptHits.resize(0);
      evt.wireHitRange.resize(0);
    }

  private:
    recob::Hit MergeTPHitsOnWire(std::vector<unsigned int>& tpHits) const;

    // SHOWER VARIABLE TREE
    TTree* showertree;

    calo::CalorimetryAlg fCaloAlg;
    TMVA::Reader fMVAReader;

    std::vector<unsigned int> fAlgModCount;

    void ReconstructAllTraj(detinfo::DetectorPropertiesData const& detProp,
                            TCSlice& slc,
                            CTP_t inCTP);
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj(TCSlice& slc, CTP_t inCTP);
    // Check allTraj -> inTraj associations
    void ChkInTraj(std::string someText, TCSlice& slc);

  }; // class TrajClusterAlg

} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
