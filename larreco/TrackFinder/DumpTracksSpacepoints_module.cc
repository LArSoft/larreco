/**
 * @file   DumpTracksSpacepoints_module.cc
 * @brief  Dumps on screen the content of the tracks
 * @author Saba Sehrish
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
#include <iostream>

// support libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"

// LArSoft includes
#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/PFParticle.h"

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
     * - <b>OutputCategory</b> (string, default: "DumpTracksSpacepoints"): the category
     *   used for the output (useful for filtering)
     *
     */
    class DumpTracksSpacepoints : public art::EDAnalyzer {
    public:
        
        /// Default constructor
        explicit DumpTracksSpacepoints(fhicl::ParameterSet const& pset);
        
        /// Does the printing
        void analyze (const art::Event& evt);
        
    private:
        
        std::string fTrackModuleLabel; ///< name of module that produced the tracks
        std::string fOutputCategory; ///< category for LogInfo output
        
    }; // class DumpTracksSpacepoints
    
} // End track namespace.


//------------------------------------------------------------------------------
namespace {
    
    template <typename STREAM>
    bool PrintdQdXinView(STREAM& log,
                         const recob::Track& track,
                         geo::View_t view,
                         std::string ViewName)
    {
        auto nData = 0;
        try { nData = track.NumberdQdx(view); }
        catch (std::out_of_range) { return false; }
        log << " " << ViewName << ": " << nData;
        return true;
    } // PrintdQdXinView()
    
    
    template <typename STREAM>
    void PrintTrack(STREAM& log,
                    const unsigned int iTrack,
                    const recob::Track& track)
    {
        // print a header for the track
        const unsigned int nPoints = track.NumberTrajectoryPoints();
        log
        << "Track # " << iTrack << " ,ID: " << track.ID()
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
        
        log << "\n Track " << track.ID() << " passes through:";
        for (unsigned int iPoint=0; iPoint <nPoints; ++iPoint) {
            const TVector3& point = track.LocationAtPoint(iPoint);
            log << "\n(" << point.X() << "," << point.Y() << "," << point.Z() << ")";
        }
        
    } //PrintTrack
    
    //Print Assns Hits
    
    template <typename STREAM, typename T>
    void PrintAssns( STREAM& log,
                    const unsigned int iTrack,
                    const std::vector<art::Ptr<T>>& Objects)
    {
        
        for (unsigned int itz = 0; itz < Objects.size(); ++itz) {
            log << Objects[itz] << "   ";
        }
    }
    
    template <typename STREAM, typename T>
    void PrintSpAssns( STREAM& log,
                      const unsigned int iTrack,
                      const std::vector<art::Ptr<T>>& Objects,
                      const art::ValidHandle<std::vector<recob::SpacePoint>>& SpacePts)
    {
        for (unsigned int ipz = 0; ipz < Objects.size(); ++ipz) {
            log << Objects[ipz].key() << "," << Objects[ipz].id() << "\n";
            
            const recob::SpacePoint& sp = SpacePts->at(Objects[ipz].key());
            log << sp.XYZ()[0] <<","  << sp.XYZ()[1] << "," << sp.XYZ()[2] << "\n ";
        }//Print space pt.
    }
} // local namespace


namespace track {
    
    //-------------------------------------------------
    DumpTracksSpacepoints::DumpTracksSpacepoints(fhicl::ParameterSet const& pset)
    : EDAnalyzer        (pset)
    , fTrackModuleLabel (pset.get<std::string>("TrackModuleLabel"))
    , fOutputCategory   (pset.get<std::string>("OutputCategory", "DumpTracksSpacepoints"))
    {}
    
    //-------------------------------------------------
    void DumpTracksSpacepoints::analyze(const art::Event& evt) {
        // fetch the data to be dumped on screen
        art::InputTag TrackInputTag(fTrackModuleLabel);
        art::ValidHandle<std::vector<recob::Track>> Tracks
        = evt.getValidHandle<std::vector<recob::Track>>(TrackInputTag);
        
        art::ValidHandle<std::vector<recob::SpacePoint>> SpacePts
        = evt.getValidHandle<std::vector<recob::SpacePoint>>(TrackInputTag);
        
        mf::LogVerbatim(fOutputCategory)
        << evt.id() << " contains " << Tracks->size() << " '"
        << TrackInputTag.encode() << "' tracks and " <<
        SpacePts->size()  << TrackInputTag.encode() << "' space points";
        
        //Hits associated with the tracks
        std::unique_ptr<art::FindManyP<recob::Hit>> pHits
        (new art::FindManyP<recob::Hit>(Tracks, evt, TrackInputTag));
        
        if (pHits && !pHits->isValid()) {
            throw art::Exception(art::errors::ProductNotFound)
            << "No hit associated with '" << TrackInputTag.encode()
            << "' tracks.\n";
        }
        
        //Space Points associated with the tracks
        std::unique_ptr<art::FindManyP<recob::SpacePoint>> pSpacePoints
        (new art::FindManyP<recob::SpacePoint>(Tracks, evt, TrackInputTag));
        
        if (pSpacePoints && !pSpacePoints->isValid()) {
            throw art::Exception(art::errors::ProductNotFound)
            << "No space point associated with '" << TrackInputTag.encode()
            << "' tracks.\n";
        }
        
        //Hits associated with space points
        std::unique_ptr<art::FindManyP<recob::Hit>> spHits
        (new art::FindManyP<recob::Hit>(SpacePts, evt, TrackInputTag));
        if (spHits && !spHits->isValid()) {
            throw art::Exception(art::errors::ProductNotFound)
            << "No hit associated with '" << TrackInputTag.encode()
            << "' space points.\n";
        }
        
        //PFParticles associated with the tracks
        std::unique_ptr<art::FindManyP<recob::PFParticle>> pPFParticles
        (new art::FindManyP<recob::PFParticle>(Tracks, evt, TrackInputTag));
        if (pPFParticles && !pPFParticles->isValid()) {
            throw art::Exception(art::errors::ProductNotFound)
            << "No particle-flow particle associated with '"
            << TrackInputTag.encode() << "' tracks.\n";
        }
        
        for (unsigned int iTrack = 0; iTrack < Tracks->size(); ++iTrack) {
            const recob::Track& track = Tracks->at(iTrack);
            mf::LogVerbatim log(fOutputCategory);
            
            PrintTrack(log, iTrack, track);
            
            log << "\n Track " << track.ID() << " is associated with:"
            << " " << pHits->at(iTrack).size() << " hits,"
            << " " << pSpacePoints->at(iTrack).size() << " space points,"
            << " " << pPFParticles->at(iTrack).size() << " PF particles.";
            
            //now print the associated hits
            log << "\n Hits associted with Track " << track.ID() << " are:\n";
            if(pHits) {
                const auto& Hits = pHits->at(iTrack);
                PrintAssns(log, iTrack, Hits);
            }
            
            //now print the associated space points and their associated hits
            log << "\n Space Points associted with Track " << track.ID() << " are:\n";
            if(pSpacePoints) {
                const auto& SpacePoints = pSpacePoints->at(iTrack);
                PrintSpAssns(log, iTrack, SpacePoints, SpacePts);
                for (unsigned int isph = 0; isph < SpacePoints.size(); ++isph) {
                    if(spHits) {
                        const auto& pHits = spHits->at(isph);
                        PrintAssns(log, isph, pHits);
                    }
                }
            }
            
            //now print the associated PF Particles
            log << "\n PF Particles associted with Track " << track.ID() << " are:\n";
            if(pPFParticles) {
                const auto& PFParticles = pPFParticles->at(iTrack);
                PrintAssns(log, iTrack, PFParticles);
            }
        } // for tracks
    } // DumpTracksSpacepoints::analyze()
    
    DEFINE_ART_MODULE(DumpTracksSpacepoints)
    
} // namespace track
