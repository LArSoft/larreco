/**
 * @file   DumpClusters_module.cc
 * @brief  Dumps on screen the content of the clusters
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 13th, 2014
 */

// C//C++ standard libraries
#include <string>
#include <sstream>
#include <iomanip> // std::setw()
#include <algorithm> // std::min(), std::sort()

// support libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"

// LArSoft includes
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"


namespace cluster {

  /**
   * @brief Prints the content of all the clusters on screen
   *
   * This analyser prints the content of all the clusters into the
   * LogInfo/LogVerbatim stream.
   * 
   * <b>Configuration parameters</b>
   * 
   * - <b>ClusterModuleLabel</b> (string, required): label of the
   *   producer used to create the recob::Cluster collection to be dumped
   * - <b>HitModuleLabel</b> (string): label of the producer used to create the
   *   recob::Hit collection the clusters are based on (not used)
   * - <b>OutputCategory</b> (string, default: "DumpClusters"): the category
   *   used for the output (useful for filtering)
   * - <b>HitsPerLine</b> (integer, default: 40): the dump of hits
   *   will put this many of them for each line
   *
   */
  class DumpClusters : public art::EDAnalyzer {
      public:
    
    /// Default constructor
    explicit DumpClusters(fhicl::ParameterSet const& pset); 
    
    /// Does the printing
    void analyze (const art::Event& evt);

      private:

    std::string fClusterModuleLabel; ///< name of module that produced clusters
//    std::string fHitsModuleLabel; ///< name of module that produced the hits
    std::string fOutputCategory; ///< category for LogInfo output
    unsigned int fHitsPerLine; ///< hits per line in the output

  }; // class DumpWires

} // End cluster namespace.


//------------------------------------------------------------------------------
namespace {
  
  /// Returns the length of the string representation of the specified object
  template <typename T>
  size_t StringLength(const T& value) {
    std::ostringstream sstr;
    sstr << value;
    return sstr.str().length();
  } // StringLength()
  
} // local namespace

namespace cluster {

  //-------------------------------------------------
  DumpClusters::DumpClusters(fhicl::ParameterSet const& pset) 
    : EDAnalyzer         (pset) 
    , fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel"))
//    , fHitsModuleLabel   (pset.get<std::string>("HitModuleLabel", ""))
    , fOutputCategory    (pset.get<std::string>("OutputCategory", "DumpClusters"))
    , fHitsPerLine       (pset.get<unsigned int>("HitsPerLine", 20))
    {}


  //-------------------------------------------------
  void DumpClusters::analyze(const art::Event& evt) {
    
    // fetch the data to be dumped on screen
  //  art::InputTag HitInputTag(fHitsModuleLabel);
    art::InputTag ClusterInputTag(fClusterModuleLabel);
    
  //  art::ValidHandle<std::vector<recob::Hit>> Hits
  //    = evt.getValidHandle<std::vector<recob::Hit>>(HitInputTag);
    art::ValidHandle<std::vector<recob::Cluster>> Clusters
      = evt.getValidHandle<std::vector<recob::Cluster>>(ClusterInputTag);
    
    // get cluster-hit associations
    art::FindManyP<recob::Hit> HitAssn(Clusters, evt, ClusterInputTag);
    
    mf::LogInfo(fOutputCategory)
      << "The event contains " << Clusters->size() << " '"
      << ClusterInputTag.encode() << "' clusters";
    
    unsigned int iCluster = 0;
    std::vector<size_t> HitBuffer(fHitsPerLine), LastBuffer;
    for (const recob::Cluster& cluster: *Clusters) {
      std::vector< art::Ptr<recob::Hit> > ClusterHits = HitAssn.at(iCluster);
      
      // print a header for the cluster
      mf::LogVerbatim(fOutputCategory)
        << "Cluster #" << (iCluster++) << " from " << ClusterHits.size()
        << " hits: " << cluster;
      
      
      // print the hits of the cluster
      if ((fHitsPerLine > 0) && !ClusterHits.empty()) {
        std::vector<size_t> HitIndices;
        for (art::Ptr<recob::Hit> pHit: ClusterHits)
          HitIndices.push_back(pHit.key());
        std::sort(HitIndices.begin(), HitIndices.end());
        
        unsigned int Padding = ::StringLength(HitIndices.back());
        
        mf::LogVerbatim(fOutputCategory) << "  hit indices:";
        
        std::vector<size_t>::const_iterator iHit = HitIndices.begin(),
          hend = HitIndices.end();
        size_t RangeStart = *iHit, RangeStop = RangeStart;
        std::ostringstream output_line;
        size_t nItemsInLine = 0;
        while (++iHit != hend) {
          
          if (*iHit == RangeStop + 1) {
            ++RangeStop;
          }
          else {
            // the new item does not belong to the current range:
            // - print the current range
            if (RangeStart == RangeStop) {
              output_line << "  " << std::setw(Padding) << RangeStart;
              ++nItemsInLine;
            }
            else {
              char fill = (RangeStart + 1 == RangeStop)? ' ': '-';
              output_line << "  " << std::setw(Padding) << RangeStart
                << fill << fill
                << std::setw(Padding) << std::setfill(fill) << RangeStop
                << std::setfill(' ');
              nItemsInLine += 2;
            }
            // - start a new one
            RangeStart = RangeStop = *iHit;
          } // if ... else
          
          // if we have enough stuff in the buffer, let's print it
          if (nItemsInLine >= fHitsPerLine) {
            nItemsInLine = 0;
            mf::LogVerbatim(fOutputCategory) << " " << output_line.str();
            output_line.str("");
          }
          
        } // while
        
        mf::LogVerbatim line_out(fOutputCategory);
        line_out << " " << output_line.str();
        if (RangeStart == RangeStop)
          line_out << "  " << std::setw(Padding) << RangeStart;
        else {
          char fill = (RangeStart + 1 == RangeStop)? ' ': '-';
          line_out << "  " << std::setw(Padding) << RangeStart
            << fill << fill
            << std::setw(Padding) << std::setfill(fill) << RangeStop;
        }
      } // if dumping the hits
    
    } // for clusters
    
  } // DumpClusters::analyze()

  DEFINE_ART_MODULE(DumpClusters)

} // namespace cluster
