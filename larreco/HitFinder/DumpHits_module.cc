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
#include "art/Framework/Core/FindOne.h"
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
   * - *CheckWireAssociation* (string, default: false): if set, verifies
   *   that the associated wire are on the same channel as the hit
   * - *CheckRawDigitAssociation* (string, default: false): if set, verifies
   *   that the associated raw digits are on the same channel as the hit
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
    std::string fOutputCategory;  ///< category for LogInfo output
    bool bCheckRawDigits;         ///< check associations with raw digits
    bool bCheckWires;             ///< check associations with wires
    
  }; // class DumpHits
  
} // namespace hit


//------------------------------------------------------------------------------
//---  module implementation
//---
// C//C++ standard libraries
#include <vector>
#include <memory> // std::unique_ptr<>

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"

// LArSoft includes
#include "larcore/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Wire.h"
#include "lardata/RawData/RawDigit.h"


namespace hit {
  
  //-------------------------------------------------
  DumpHits::DumpHits(fhicl::ParameterSet const& pset) 
    : EDAnalyzer         (pset) 
    , fHitsModuleLabel   (pset.get<std::string>("HitModuleLabel"))
    , fOutputCategory    (pset.get<std::string>("OutputCategory", "DumpHits"))
    , bCheckRawDigits    (pset.get<bool>       ("CheckRawDigitAssociation", false))
    , bCheckWires        (pset.get<bool>       ("CheckWireAssociation", false))
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
    
    std::unique_ptr<art::FindOne<raw::RawDigit>> HitToRawDigit;
    if (bCheckRawDigits) {
      HitToRawDigit.reset
        (new art::FindOne<raw::RawDigit>(Hits, evt, HitInputTag));
      if (!HitToRawDigit->isValid()) {
        throw art::Exception(art::errors::ProductNotFound)
          << "DumpHits: can't find associations between raw digits and hits from '"
          << HitInputTag << "'";
      }
    } // if check raw digits
    
    std::unique_ptr<art::FindOne<recob::Wire>> HitToWire;
    if (bCheckWires) {
      HitToWire.reset(new art::FindOne<recob::Wire>(Hits, evt, HitInputTag));
      if (!HitToWire->isValid()) {
        throw art::Exception(art::errors::ProductNotFound)
          << "DumpHits: can't find associations between wires and hits from '"
          << HitInputTag << "'";
      }
    } // if check wires
    
    unsigned int iHit = 0;
    for (const recob::Hit& hit: *Hits) {
      
      // print a header for the cluster
      mf::LogVerbatim(fOutputCategory)
        << "Hit #" << iHit << ": " << hit;
      
      if (HitToRawDigit) {
        raw::ChannelID_t assChannelID = HitToRawDigit->at(iHit).ref().Channel();
        if (assChannelID != hit.Channel()) {
          throw art::Exception(art::errors::DataCorruption)
            << "Hit #" << iHit << " on channel " << hit.Channel()
            << " is associated with raw digit on channel " << assChannelID
            << "!!";
        } // mismatch
      } // raw digit check
      
      if (HitToWire) {
        raw::ChannelID_t assChannelID = HitToWire->at(iHit).ref().Channel();
        if (assChannelID != hit.Channel()) {
          throw art::Exception(art::errors::DataCorruption)
            << "Hit #" << iHit << " on channel " << hit.Channel()
            << " is associated with wire on channel " << assChannelID
            << "!!";
        } // mismatch
      } // wire check
      
      ++iHit;
    } // for hits
    
  } // DumpHits::analyze()
  
  DEFINE_ART_MODULE(DumpHits)
  
} // namespace hit
