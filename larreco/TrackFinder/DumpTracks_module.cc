/**
 * @file   DumpTracks_module.cc
 * @brief  Dumps on screen the content of the tracks
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 12th, 2014
 */

// C//C++ standard libraries
#include <string>
#include <sstream>
#include <iomanip> // std::setw()
#include <algorithm> // std::max(), std::sort(), std::transform()
#include <iterator> // std::back_inserter()
#include <functional> // std::mem_fn()
#include <memory> // std::unique_ptr()

// support libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

// LArSoft includes
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/PFParticle.h"


namespace {
  
  /// Returns the length of the string representation of the specified object
  template <typename T>
  size_t StringLength(const T& value) {
    std::ostringstream sstr;
    sstr << value;
    return sstr.str().length();
  } // StringLength()
  
  /**
   * @brief Prints a table with compact integers
   * @param log stream to send the output to
   * @param Indices sorted container of the indices to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndentStr string to be printed at the beginning of each new line
   *
   */
  template <typename STREAM, typename CONT>
  void PrintCompactIndexTable(
    STREAM& log, const CONT& Indices, unsigned int IndicesPerLine,
    std::string IndentStr = "")
  {
    unsigned int Padding = StringLength(Indices.back());
    
    typename CONT::const_iterator iIndex = Indices.begin(),
      iend = Indices.end();
    size_t RangeStart = *iIndex, RangeStop = RangeStart;
    std::ostringstream output_line;
    size_t nItemsInLine = 0;
    while (++iIndex != iend) {
      
      if (*iIndex == RangeStop + 1) {
        ++RangeStop;
      }
      else {
        // the new item does not belong to the current range:
        // - print the current range
        if (nItemsInLine) output_line << "  ";
        if (RangeStart == RangeStop) {
          output_line << std::setw(Padding) << RangeStart;
          ++nItemsInLine;
        }
        else {
          char fill = (RangeStart + 1 == RangeStop)? ' ': '-';
          output_line << std::setw(Padding) << RangeStart
            << fill << fill
            << std::setw(Padding) << std::setfill(fill) << RangeStop
            << std::setfill(' ');
          nItemsInLine += 2;
        }
        // - start a new one
        RangeStart = RangeStop = *iIndex;
      } // if ... else
      
      // if we have enough stuff in the buffer, let's print it
      if (nItemsInLine >= IndicesPerLine) {
        nItemsInLine = 0;
        log << IndentStr << output_line.str() << "\n";
        output_line.str("");
      }
      
    } // while
    
    // print what we have accumulated so far of the last line
    log << IndentStr << output_line.str();
    // print the last range (or single element)
    if (nItemsInLine) log << "  ";
    if (RangeStart == RangeStop)
      log << std::setw(Padding) << RangeStart;
    else {
      char fill = (RangeStart + 1 == RangeStop)? ' ': '-';
      log << std::setw(Padding) << RangeStart
        << fill << fill
        << std::setw(Padding) << std::setfill(fill) << RangeStop
        << std::setfill(' ');
    }
    log << std::endl;
  } // PrintCompactIndexTable()
  
  
  
  /**
   * @brief Prints a table with keys from a container of objects
   * @param log stream to send the output to
   * @param Indices sorted container of the indices to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndexExtractor functor extracting the index from a element
   * @param IndentStr string to be printed at the beginning of each new line
   *
   * The key extractor must oterate on elements of Indices and returning
   * something that can be converted to a size_t.
   */
  template <typename STREAM, typename CONT, typename GETINDEX>
  void PrintCompactIndexTable(
    STREAM& log, const CONT& Objects, unsigned int IndicesPerLine,
    GETINDEX IndexExtractor, std::string IndentStr)
  {
    if ((IndicesPerLine == 0) || Objects.empty()) return;
    
    std::vector<size_t> Indices;
    Indices.reserve(Objects.size());
    std::transform(Objects.begin(), Objects.end(), std::back_inserter(Indices),
      IndexExtractor);
    std::sort(Indices.begin(), Indices.end());
    PrintCompactIndexTable(log, Indices, IndicesPerLine, IndentStr);
    
  } // PrintCompactIndexTable()
  
  /**
   * @brief Prints a table with indices from a container of objects
   * @param log stream to send the output to
   * @param Objects container of art::Ptr to objects to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndentStr string to be printed at the beginning of each new line
   *
   * The class in the container must have a key() member returning something
   * that can be converted to a size_t.
   * This is designed for containers like std::vector<art::Ptr<T>>.
   */
  template <typename STREAM, typename T>
  inline void PrintAssociatedIndexTable(
    STREAM& log, const std::vector<art::Ptr<T>>& Objects,
    unsigned int IndicesPerLine, std::string IndentStr = ""
  ) {
    PrintCompactIndexTable(
      log, Objects, IndicesPerLine, std::mem_fn(&art::Ptr<T>::key),
      IndentStr
      );
  } // PrintAssociatedIndexTable()
  
  
  /**
   * @brief Prints a table with indices from a container of objects
   * @param log stream to send the output to
   * @param Objects container of art::Ptr to objects to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndentStr string to be printed at the beginning of each new line
   *
   * The class in the container must have a key() member returning something
   * that can be converted to a size_t.
   * This is designed for containers like std::vector<art::Ptr<T>>.
   */
  template <typename STREAM, typename T>
  inline void PrintAssociatedIDTable(
    STREAM& log, const std::vector<art::Ptr<T>>& Objects,
    unsigned int IndicesPerLine, std::string IndentStr = ""
  ) {
    PrintCompactIndexTable(
      log, Objects, IndicesPerLine,
      [](const art::Ptr<T>& ptr) { return ptr->ID(); },
      IndentStr
      );
  } // PrintAssociatedIDTable()
  
} // local namespace


namespace track {

  /**
   * @brief Prints the content of all the tracks on screen
   *
   * This analyser prints the content of all the tracks into the
   * LogInfo/LogVerbatim stream.
   * 
   * <b>Configuration parameters</b>
   * 
   * - <b>TrackModuleLabel</b> (string, required): label of the
   *   producer used to create the recob::Track collection to be dumped
   * - <b>OutputCategory</b> (string, default: "DumpTracks"): the category
   *   used for the output (useful for filtering)
   * - <b>WayPoints</b> (unsigned integer, default: 10): approximate number
   *   of way points printed in the output
   * - <b>SpacePointAssociations</b> (boolean, default: true): prints the number
   *   of space points associated with the tracks
   * - <b>PrintSpacePoints</b> (boolean, default: false): print the index of all hits
   *   associated with the tracks
   * - <b>HitAssociations</b> (boolean, default: true): prints the number of
   *   hits associated with the tracks
   * - <b>PrintHits</b> (boolean, default: false): print the index of all hits
   *   associated with the tracks
   * - <b>ParticleAssociations</b> (boolean, default: true): prints the number
   *   of particle-flow particles associated with the tracks
   *
   */
  class DumpTracks : public art::EDAnalyzer {
      public:
    
    /// Default constructor
    explicit DumpTracks(fhicl::ParameterSet const& pset); 
    
    /// Does the printing
    void analyze (const art::Event& evt);

      private:

    std::string fTrackModuleLabel; ///< name of module that produced the tracks
    std::string fOutputCategory; ///< category for LogInfo output
    unsigned int fPrintWayPoints; ///< number of printed way points
    
    bool fPrintNHits; ///< prints the number of associated hits
    bool fPrintNSpacePoints; ///< prints the number of associated space points
    bool fPrintNParticles; ///< prints the number of associated PFParticles
    bool fPrintHits; ///< prints the index of associated hits
    bool fPrintSpacePoints; ///< prints the index of associated space points
    bool fPrintParticles; ///< prints the index of associated PFParticles

  }; // class DumpTracks

} // End track namespace.


//------------------------------------------------------------------------------
namespace {
  
  template <typename STREAM>
  bool PrintdQdXinView(
    STREAM& log,
    const recob::Track& track, geo::View_t view, std::string ViewName
  ) {
    auto nData = 0;
    try { nData = track.NumberdQdx(view); }
    catch (std::out_of_range) { return false; }
    log << " " << ViewName << ": " << nData;
    return true;
  } // PrintdQdXinView()
  
} // local namespace


namespace track {

  //-------------------------------------------------
  DumpTracks::DumpTracks(fhicl::ParameterSet const& pset) 
    : EDAnalyzer        (pset) 
    , fTrackModuleLabel (pset.get<std::string>("TrackModuleLabel"))
    , fOutputCategory   (pset.get<std::string>("OutputCategory", "DumpTracks"))
    , fPrintWayPoints   (pset.get<unsigned int>("WayPoints", 10))
    , fPrintNHits       (pset.get<bool>("HitAssociations", true))
    , fPrintNSpacePoints(pset.get<bool>("SpacePointAssociations", true))
    , fPrintNParticles  (pset.get<bool>("ParticleAssociations", true))
    , fPrintHits        (pset.get<bool>("PrintHits", false))
    , fPrintSpacePoints (pset.get<bool>("PrintSpacePoints", false))
    , fPrintParticles   (pset.get<bool>("PrintParticles", false))
    {}
  //-------------------------------------------------
  void DumpTracks::analyze(const art::Event& evt) {
    
    // fetch the data to be dumped on screen
    art::InputTag TrackInputTag(fTrackModuleLabel);
    art::ValidHandle<std::vector<recob::Track>> Tracks
      = evt.getValidHandle<std::vector<recob::Track>>(TrackInputTag);
    
    mf::LogInfo(fOutputCategory)
      << "The event contains " << Tracks->size() << " '"
      << TrackInputTag.encode() << "'tracks";
    
    std::unique_ptr<art::FindManyP<recob::Hit>> pHits(
      fPrintNHits?
      new art::FindManyP<recob::Hit>(Tracks, evt, TrackInputTag):
      nullptr
      );
    if (pHits && !pHits->isValid()) {
      throw art::Exception(art::errors::ProductNotFound)
        << "No hit associated with '" << TrackInputTag.encode()
        << "' tracks.\n";
    }
    
    std::unique_ptr<art::FindManyP<recob::SpacePoint>> pSpacePoints(
      fPrintNSpacePoints?
      new art::FindManyP<recob::SpacePoint>(Tracks, evt, TrackInputTag):
      nullptr
      );
    if (pSpacePoints && !pSpacePoints->isValid()) {
      throw art::Exception(art::errors::ProductNotFound)
        << "No space point associated with '" << TrackInputTag.encode()
        << "' tracks.\n";
    }
    
    std::unique_ptr<art::FindManyP<recob::PFParticle>> pPFParticles(
      fPrintNParticles?
      new art::FindManyP<recob::PFParticle>(Tracks, evt, TrackInputTag):
      nullptr
      );
    if (pPFParticles && !pPFParticles->isValid()) {
      throw art::Exception(art::errors::ProductNotFound)
        << "No particle-flow particle associated with '"
        << TrackInputTag.encode() << "' tracks.\n";
    }
    
    for (unsigned int iTrack = 0; iTrack < Tracks->size(); ++iTrack) {
      const recob::Track& track = Tracks->at(iTrack);
      
      // print a header for the track
      const unsigned int nPoints = track.NumberTrajectoryPoints();
      mf::LogVerbatim log(fOutputCategory);
      log
        << "Track #" << iTrack << " ID: " << track.ID()
          << std::fixed << std::setprecision(3)
          << " theta: " << track.Theta() << " rad, phi: " << track.Phi()
          << " rad, length: " << track.Length() << " cm"
        << "\n  start at: ( " << track.Vertex().X()
          << " ; " << track.Vertex().Y()
          << " ; " << track.Vertex().Z()
          << " ), direction: ( " << track.VertexDirection().X()
          << " ; " << track.VertexDirection().Y()
          << " ; " << track.VertexDirection().Z() << " )"
        << "\n  end at:   ( " << track.End().X()
          << " ; " << track.End().Y()
          << " ; " << track.End().Z()
          << " ), direction: ( " << track.EndDirection().X()
          << " ; " << track.EndDirection().Y()
          << " ; " << track.EndDirection().Z()
          << " )"
        << "\n  with "
          << nPoints << " trajectory points, "
          << track.NumberCovariance() << " covariance matrices";
      unsigned int nViews = 0;
      std::ostringstream sstr;
      if (PrintdQdXinView(sstr, track, geo::kU, "U")) ++nViews;
      if (PrintdQdXinView(sstr, track, geo::kV, "V")) ++nViews;
      if (PrintdQdXinView(sstr, track, geo::kZ, "Z")) ++nViews;
      if (PrintdQdXinView(sstr, track, geo::k3D, "3D")) ++nViews;
      if (nViews)
        log << ", dQ/dx in " << nViews << " views: " << sstr.str();
      else
        log << ", no dQ/dx";
      
      if (fPrintWayPoints > 0) {
        // print up to 10 (actually, 8 or 9) way points
        log << "\n  passes through:";
        unsigned int skip = std::max(nPoints / fPrintWayPoints, 1U);
        unsigned int iPoint = 0;
        while ((iPoint += skip) < nPoints) {
          const TVector3& point = track.LocationAtPoint(iPoint);
          log << "\n    [#" << iPoint << "] ("
            << point.X() << ", " << point.Y() << ", " << point.Z()
            << ")";
        } // while (iPoint)
      } // if print way points
      
      if (pHits || pSpacePoints || pPFParticles) {
        log << "\n  associated with:";
        if (pHits)
          log << " " << pHits->at(iTrack).size() << " hits;";
        if (pSpacePoints)
          log << " " << pSpacePoints->at(iTrack).size() << " space points;";
        if (pPFParticles)
          log << " " << pPFParticles->at(iTrack).size() << " PF particles;";
      } // if we have any association
      
      if (pHits && fPrintHits) {
        const auto& Hits = pHits->at(iTrack);
        log << "\n  hit indices (" << Hits.size() << "):\n";
        PrintAssociatedIndexTable(log, Hits, 10 /* 10 hits per line */, "    ");
      } // if print individual hits
      
      if (pSpacePoints && fPrintSpacePoints) {
        const auto& SpacePoints = pSpacePoints->at(iTrack);
        log << "\n  space point IDs (" << SpacePoints.size() << "):\n";
        PrintAssociatedIDTable
          (log, SpacePoints, 10 /* 10 hits per line */, "    ");
      } // if print individual space points
      
      if (pPFParticles && fPrintParticles) {
        const auto& PFParticles = pPFParticles->at(iTrack);
        log << "\n  particle indices (" << PFParticles.size() << "):\n";
        // currently a particle has no ID
        PrintAssociatedIndexTable
          (log, PFParticles, 10 /* 10 hits per line */, "    ");
      } // if print individual particles
    } // for tracks
  } // DumpTracks::analyze()

  DEFINE_ART_MODULE(DumpTracks)

} // namespace track
