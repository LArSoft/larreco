/**
 * @file   DumpHits_module.cc
 * @brief  Dumps on screen the content of the hits
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   Match 9th, 2015
 */

// C//C++ standard libraries
#include <string>

// support libraries
#include "fhiclcpp/ParameterSet.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// ... plus see below ...

namespace hit {
  
  /**
   * @brief Prints the content of all the hits on screen
   *
   * This analyser prints the content of all the hits into the
   * LogInfo/LogVerbatim stream.
   * 
   * Configuration parameters
   * =========================
   * 
   * - *HitModuleLabel* (string): label of the producer used to create the
   *   recob::Hit collection
   * - *OutputCategory* (string, default: "DumpHits"): the category
   *   used for the output (useful for filtering)
   *
   */
  class DumpHits: public art::EDAnalyzer {
      public:
    
    /// Default constructor
    explicit DumpHits(fhicl::ParameterSet const& pset); 
    
    /// Does the printing
    void analyze (const art::Event& evt);
    
      private:
    
    std::string fHitsModuleLabel; ///< name of module that produced the hits
    std::string fOutputCategory; ///< category for LogInfo output
    
  }; // class DumpHits
  
} // namespace hit


//------------------------------------------------------------------------------
//---  module implementation
//---
// C//C++ standard libraries
#include <vector>

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"

// LArSoft includes
#include "RecoBase/Hit.h"


namespace hit {
  
  //-------------------------------------------------
  DumpHits::DumpHits(fhicl::ParameterSet const& pset) 
    : EDAnalyzer         (pset) 
    , fHitsModuleLabel   (pset.get<std::string>("HitModuleLabel"))
    , fOutputCategory    (pset.get<std::string>("OutputCategory", "DumpHits"))
    {}
  
  
  //-------------------------------------------------
  void DumpHits::analyze(const art::Event& evt) {
    
    // fetch the data to be dumped on screen
    art::InputTag HitInputTag(fHitsModuleLabel);
    
    art::ValidHandle<std::vector<recob::Hit>> Hits
      = evt.getValidHandle<std::vector<recob::Hit>>(HitInputTag);
    
    mf::LogInfo(fOutputCategory)
      << "The event contains " << Hits->size() << " '"
      << HitInputTag.encode() << "' hits";
    
    unsigned int iHit = 0;
    for (const recob::Hit& hit: *Hits) {
      
      // print a header for the cluster
      mf::LogVerbatim(fOutputCategory)
        << "Hit #" << (iHit++) << ": " << hit;
    
    } // for clusters
    
  } // DumpHits::analyze()
  
  DEFINE_ART_MODULE(DumpHits)
  
} // namespace hit
