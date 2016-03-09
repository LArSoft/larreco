/**
 * @file   TrackUtils.cxx
 * @brief  Utility functions to extract information from recob::Track - implementation
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 8th, 2016
 * @see    TrackUtils.h
 * 
 */

// out header
#include "larreco/RecoAlg/TrackUtils.h"

// LArSoft libraries
#include "larcore/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Track.h"

// framework libraries
#include "cetlib/exception.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <cmath>



//------------------------------------------------------------------------------
double lar::TrackProjectedLength(recob::Track const& track, geo::View_t view) {
                  
   if(view == geo::kUnknown) {
      throw cet::exception("TrackProjectedLength") << "cannot provide projected length for "
        << "unknown view\n";
   }
   double length = 0.;
   
   auto const* geom = lar::providerFrom<geo::Geometry>();
   double angleToVert = 0.;
   for(unsigned int i = 0; i < geom->Nplanes(); ++i){
      if(geom->Plane(i).View() == view){
         angleToVert = geom->Plane(i).Wire(0).ThetaZ(false) - 0.5*util::pi<>();
         break;
      }
   }

   // now loop over all points in the trajectory and add the contribution to the
   // to the desired view

   for(size_t p = 1; p < track.NumberTrajectoryPoints(); ++p){
      const TVector3& pos_cur = track.LocationAtPoint(p);
      const TVector3& pos_prev = track.LocationAtPoint(p - 1);
      double dist = std::sqrt( std::pow(pos_cur.x() - pos_prev.x(), 2) +
         std::pow(pos_cur.y() - pos_prev.y(), 2) +
         std::pow(pos_cur.z() - pos_prev.z(), 2) );
      
      // (sin(angleToVert),cos(angleToVert)) is the direction perpendicular to wire
      // fDir[p-1] is the direction between the two relevant points
      const TVector3& dir_prev = track.DirectionAtPoint(p - 1);
      double cosgamma = std::abs(std::sin(angleToVert)*dir_prev.Y() +
         std::cos(angleToVert)*dir_prev.Z() );
      
      /// @todo is this right, or should it be dist*cosgamma???
      length += dist/cosgamma;
   } // end loop over distances between trajectory points
   
   return length;
} // lar::TrackProjectedLength()



//------------------------------------------------------------------------------
double lar::TrackPitchInView
  (recob::Track const& track, geo::View_t view, size_t trajectory_point /* = 0 */)
{
   
   if(view == geo::kUnknown) {
      throw cet::exception("TrackPitchInView") << "Warning cannot obtain pitch for unknown view\n";
   }
   
   if(trajectory_point >= track.NumberTrajectoryPoints()) {
      cet::exception("TrackPitchInView") << "ERROR: Asking for trajectory point #" 
         << trajectory_point << " when trajectory vector size is of size "
         << track.NumberTrajectoryPoints() << ".\n";
   }
   
   auto const* geom = lar::providerFrom<geo::Geometry>();
   const TVector3& pos = track.LocationAtPoint(trajectory_point);
   const double Position[3] = { pos.X(), pos.Y(), pos.Z() };
   geo::TPCID tpcid = geom->FindTPCAtPosition ( Position );
   if (!tpcid.isValid) {
      cet::exception("TrackPitchInView") << "ERROR: invalid TPC " << std::string(tpcid)
         << " for trajectory point #" << trajectory_point  << ".\n";
   }
   double wirePitch   = geom->WirePitch(view, tpcid.TPC, tpcid.Cryostat);
   double angleToVert = geom->WireAngleToVertical(view, tpcid.TPC, tpcid.Cryostat) - 0.5*util::pi<>();
                                   
   const TVector3& dir = track.DirectionAtPoint(trajectory_point);
   //(sin(angleToVert),cos(angleToVert)) is the direction perpendicular to wire
   double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() +
     std::cos(angleToVert)*dir.Z());
   
   if(cosgamma < 1.e-5)
      throw cet::exception("Track") << "cosgamma is basically 0, that can't be right\n";
   return wirePitch/cosgamma;
                                                                                                  
} // lar::TrackPitchInView()
   
//------------------------------------------------------------------------------
