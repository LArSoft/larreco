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

#include "larreco/RecoAlg/TCAlg/StepUtils.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/TCCR.h"

// C/C++ standard libraries
#include <array>
#include <vector>
#include <utility> // std::pair<>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larreco/Calorimetry/CalorimetryAlg.h"

//#include "TH1F.h"
//#include "TH2F.h"
//#include "TProfile.h"
#include "TTree.h"

namespace tca {

  class TrajClusterAlg {
    public:
    
 
    /// @{
    /// @name Data structures for the reconstruction results
    
    TrajClusterAlg(fhicl::ParameterSet const& pset);

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    void SetMCPHandle(std::vector<simb::MCParticle> const& mcpHandle) { evt.mcpHandle = &mcpHandle; }
    bool SetInputHits(std::vector<recob::Hit> const& inputHits);
    void RunTrajClusterAlg(std::vector<unsigned int>& hitsInSlice, int sliceID);
    bool CreateSlice(std::vector<unsigned int>& hitsInSlice);
    void FinishEvent();
    

    void DefineShTree(TTree* t);
    
//    void DefineCRTree(TTree* t);

    unsigned short GetSlicesSize() { return slices.size(); }
    TCSlice const& GetSlice(unsigned short sliceIndex) const {return slices[sliceIndex]; }
    void MergeTPHits(std::vector<unsigned int>& tpHits, std::vector<recob::Hit>& newHitCol, 
                     std::vector<unsigned int>& newHitAssns);
    
    std::vector<unsigned int> const& GetAlgModCount() const {return fAlgModCount; }
    std::vector<std::string> const& GetAlgBitNames() const {return AlgBitNames; }
    
    /// Deletes all the results
    void ClearResults() { slices.resize(0); evt.allHitsMCPIndex.resize(0); evt.allHitsRanges.resize(0);}
    TruthMatcher fTM;
    
    private:

    recob::Hit MergeTPHitsOnWire(std::vector<unsigned int>& tpHits);

    // SHOWER VARIABLE TREE
    TTree* showertree;

    // Cosmic Removal Variable Tree
//    TTree* crtree;
    
    calo::CalorimetryAlg fCaloAlg;
    TMVA::Reader fMVAReader;
        
    std::vector<unsigned int> fAlgModCount;

    void ReconstructAllTraj(TCSlice& slc, CTP_t inCTP);
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj(TCSlice& slc, CTP_t inCTP);
    // Check allTraj -> inTraj associations
    void ChkInTraj(std::string someText, TCSlice& slc);
    void FindMissedVxTjs(TCSlice& slc);

  }; // class TrajClusterAlg

} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
