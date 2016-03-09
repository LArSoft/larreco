/**
 * @file   TrackUtils.h
 * @brief  Utility functions to extract information from recob::Track
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 8th, 2016
 * @see    TrackUtils.cxx
 * 
 * 
 * lar::TrackProjectedLength() and lar::TrackPitchInView() have been factored
 * out from recob::Track::ProjectedLength() and recob::Track::PitchInView()
 * respectively.
 */

#ifndef LARDATA_RECOALG_TRACKUTILS_H
#define LARDATA_RECOALG_TRACKUTILS_H 1

// LArSoft libraries
#include "lardata/RecoBase/Track.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h" // geo::View_t


namespace lar {
   
   /**
    * @brief Returns the length of the projection of a track on a view
    * @param track the track to be projected on a view
    * @param view the view to project the track on
    * @return length of the projection, in centimetres
    * 
    * 
    * @todo CAREFUL: using view to determine projected length does not work for DUNE.
    * Need to think more about this.
    * 
    */
   double TrackProjectedLength(recob::Track const& track, geo::View_t view);
   
   
   /**
    * @brief Provides projected wire pitch for the view
    * @param track the track to be projected on a view
    * @param view the view for track projection
    * @param trajectory_point at which point of the track to look for the pitch (default: 0)
    * @return wire pitch along the track direction at its specified point, in centimetres
    * 
    * The point of the trajectory with the specified index (by default, the
    * first point of the trajectory) is projected into the
    * specified view, and the effective distance between two wires along the
    * direction the track is pointing to at that point is computed and returned.
    */
   double TrackPitchInView
     (recob::Track const& track, geo::View_t view, size_t trajectory_point = 0);
   
   
} // namespace recob


#endif // LARDATA_RECOALG_TRACKUTILS_H
