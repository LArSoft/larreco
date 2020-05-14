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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT libraries
#include "TMVA/Reader.h"
class TTree;

namespace tca {

  class TrajClusterAlg {
  public:
    /// @{
    /// @name Data structures for the reconstruction results

    TrajClusterAlg(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset);

    bool SetInputHits(std::vector<recob::Hit> const& inputHits,
                      unsigned int run,
                      unsigned int event);
    void
    SetInputSpts(std::vector<recob::SpacePoint> const& sptHandle)
    {
      evt.sptHandle = &sptHandle;
    }
    void SetSourceHits(std::vector<recob::Hit> const& srcHits);
    void
    ExpectSlicedHits()
    {
      evt.expectSlicedHits = true;
    }
    void RunTrajClusterAlg(std::vector<unsigned int>& hitsInSlice, int sliceID);
    bool CreateSlice(std::vector<unsigned int>& hitsInSlice, int sliceID);
    void FinishEvent();

    void DefineShTree(TTree* t);

    unsigned short
    GetSlicesSize() const
    {
      return slices.size();
    }
    TCSlice const&
    GetSlice(unsigned short sliceIndex) const
    {
      return slices[sliceIndex];
    }
    void MergeTPHits(std::vector<unsigned int>& tpHits,
                     std::vector<recob::Hit>& newHitCol,
                     std::vector<unsigned int>& newHitAssns) const;

    std::vector<unsigned int> const&
    GetAlgModCount() const
    {
      return fAlgModCount;
    }
    std::vector<std::string> const&
    GetAlgBitNames() const
    {
      return AlgBitNames;
    }

    /// Deletes all the results
    void
    ClearResults()
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

    void ReconstructAllTraj(TCSlice& slc, CTP_t inCTP);
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj(TCSlice& slc, CTP_t inCTP);
    // Check allTraj -> inTraj associations
    void ChkInTraj(std::string someText, TCSlice& slc);

  }; // class TrajClusterAlg

} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
