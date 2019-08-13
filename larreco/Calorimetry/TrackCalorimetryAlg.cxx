/*!
 * Title:   Track Calorimetry Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code in the Calorimetry_module
 *
 * Description: Algorithm that produces a calorimetry object given a track
 * Input:       recob::Track, Assn<recob::Spacepoint,recob::Track>, Assn<recob::Hit,recob::Track>
 * Output:      anab::Calorimetry, (and Assn<anab::Calorimetry,recob::Track>)
*/

#include "fhiclcpp/ParameterSet.h"
#include "TrackCalorimetryAlg.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()

calo::TrackCalorimetryAlg::TrackCalorimetryAlg(fhicl::ParameterSet const& p):
  caloAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  this->reconfigure(p);
}

void calo::TrackCalorimetryAlg::reconfigure(fhicl::ParameterSet const& p){
  caloAlg.reconfigure(p.get<fhicl::ParameterSet>("CalorimetryAlg"));
  fNHitsToDetermineStart = p.get<unsigned int>("NHitsToDetermineStart",3);
}

void calo::TrackCalorimetryAlg::ExtractCalorimetry(std::vector<recob::Track> const& trackVector,
						   std::vector<recob::Hit> const& hitVector,
						   std::vector< std::vector<size_t> > const& hit_indices_per_track,
						   std::vector<anab::Calorimetry>& caloVector,
						   std::vector<size_t>& assnTrackCaloVector,
						   Providers_t providers)
{
  auto const& geom = *(providers.get<geo::GeometryCore>());
//  auto const& larp = *(providers.get<detinfo::LArProperties>());
  auto const& detprop = *(providers.get<detinfo::DetectorProperties>());

  //loop over the track list
  for(size_t i_track=0; i_track<trackVector.size(); i_track++){

    recob::Track const& track = trackVector[i_track];
    std::vector<float> path_length_fraction_vec(CreatePathLengthFractionVector(track));

    //sort hits into each plane
    std::vector< std::vector<size_t> > hit_indices_per_plane(geom.Nplanes());
    for(auto const& i_hit : hit_indices_per_track[i_track])
      hit_indices_per_plane[hitVector[i_hit].WireID().Plane].push_back(i_hit);

    //loop over the planes
    for(size_t i_plane=0; i_plane<geom.Nplanes(); i_plane++){

      ClearInternalVectors();
      ReserveInternalVectors(hit_indices_per_plane[i_plane].size());

      //project down the track into wire/tick space for this plane
      std::vector< std::pair<geo::WireID,float> > traj_points_in_plane(track.NumberTrajectoryPoints());
      for(size_t i_trjpt=0; i_trjpt<track.NumberTrajectoryPoints(); i_trjpt++){
	double x_pos = track.LocationAtPoint(i_trjpt).X();
	float tick = detprop.ConvertXToTicks(x_pos,(int)i_plane,0,0);
	traj_points_in_plane[i_trjpt] = std::make_pair(geom.NearestWireID(track.LocationAtPoint(i_trjpt),i_plane),
						       tick);
      }

      HitPropertiesMultiset_t HitPropertiesMultiset;
      //now loop through hits
      for(auto const& i_hit : hit_indices_per_plane[i_plane])
	AnalyzeHit(hitVector[i_hit],
		   track,
		   traj_points_in_plane,
		   path_length_fraction_vec,
		   HitPropertiesMultiset,
		   geom);


      //PrintHitPropertiesMultiset(HitPropertiesMultiset);
      geo::PlaneID planeID(0,0,i_plane);
      MakeCalorimetryObject(HitPropertiesMultiset, track, i_track, caloVector, assnTrackCaloVector, planeID);

    }//end loop over planes

  }//end loop over tracks

}//end ExtractCalorimetry


class dist_projected{
public:
  dist_projected(recob::Hit const& h, geo::GeometryCore const& g):
    hit(h), geom(g){}
  bool operator() (std::pair<geo::WireID,float> i, std::pair<geo::WireID,float> j)
  {
    float dw_i = ((int)(i.first.Wire) - (int)(hit.WireID().Wire))*geom.WirePitch(i.first.Plane);
    float dw_j = ((int)(j.first.Wire) - (int)(hit.WireID().Wire))*geom.WirePitch(j.first.Plane);
    float dt_i = i.second - hit.PeakTime();
    float dt_j = j.second - hit.PeakTime();

    return (std::sqrt(dw_i*dw_i + dt_i*dt_i) < std::sqrt(dw_j*dw_j + dt_j*dt_j));
  }
private:
  recob::Hit const& hit;
  geo::GeometryCore const& geom;

};

std::vector<float> calo::TrackCalorimetryAlg::CreatePathLengthFractionVector(recob::Track const& track){

  std::vector<float> trk_path_length_frac_vec(track.NumberTrajectoryPoints());

  float cumulative_path_length=0;
  const float total_path_length = track.Length();
  for(size_t i_trj=1; i_trj<track.NumberTrajectoryPoints(); i_trj++){
    cumulative_path_length+=(track.LocationAtPoint(i_trj)-track.LocationAtPoint(i_trj-1)).R();
    trk_path_length_frac_vec[i_trj]=cumulative_path_length/total_path_length;
  }

  return trk_path_length_frac_vec;
}

void calo::TrackCalorimetryAlg::AnalyzeHit(recob::Hit const& hit,
					   recob::Track const& track,
					   std::vector< std::pair<geo::WireID,float> > const& traj_points_in_plane,
					   std::vector<float> const& path_length_fraction_vec,
					   HitPropertiesMultiset_t & HitPropertiesMultiset,
					   geo::GeometryCore const& geom){

  size_t traj_iter = std::distance(traj_points_in_plane.begin(),
				   std::min_element(traj_points_in_plane.begin(),
						    traj_points_in_plane.end(),
						    dist_projected(hit,geom)));
  float pitch = lar::util::TrackPitchInView(track, geom.View(hit.WireID().Plane),traj_iter);

  HitPropertiesMultiset.emplace(hit.Integral(),
				hit.Integral()/pitch,
				caloAlg.dEdx_AREA(hit,pitch),
				pitch,
				track.LocationAtPoint<TVector3>(traj_iter),
				path_length_fraction_vec[traj_iter]);
}

bool calo::TrackCalorimetryAlg::IsInvertedTrack(HitPropertiesMultiset_t const& hpm){

  if(hpm.size() <= fNHitsToDetermineStart) return false;

  float charge_start=0,charge_end=0;
  unsigned int counter=0;
  for(HitPropertiesMultiset_t::iterator it_hpm=hpm.begin();
      it_hpm!=hpm.end();
      it_hpm++)
    {
      charge_start += (*it_hpm).charge;
      counter++;
      if(counter==fNHitsToDetermineStart) break;
    }

  counter=0;
  for(HitPropertiesMultiset_t::reverse_iterator it_hpm=hpm.rbegin();
      it_hpm!=hpm.rend();
      it_hpm++)
    {
      charge_end += (*it_hpm).charge;
      counter++;
      if(counter==fNHitsToDetermineStart) break;
    }

  return (charge_start > charge_end);
}

void calo::TrackCalorimetryAlg::MakeCalorimetryObject(HitPropertiesMultiset_t const& hpm,
						      recob::Track const& track,
						      size_t const& i_track,
						      std::vector<anab::Calorimetry>& caloVector,
						      std::vector<size_t>& assnTrackCaloVector,
						      geo::PlaneID const& planeID){
  size_t n_hits = hpm.size();
  std::vector<float> dEdxVector,dQdxVector,resRangeVector,deadWireVector,pitchVector;
  std::vector<TVector3> XYZVector;

  dEdxVector.reserve(n_hits);
  dQdxVector.reserve(n_hits);
  resRangeVector.reserve(n_hits);
  deadWireVector.reserve(n_hits);
  pitchVector.reserve(n_hits);
  XYZVector.reserve(n_hits);

  float kinetic_energy=0,track_length=track.Length();
  if(IsInvertedTrack(hpm)){

    for(HitPropertiesMultiset_t::iterator it_hpm=hpm.begin();
	it_hpm!=hpm.end();
	it_hpm++)
      {
	dEdxVector.push_back((*it_hpm).dEdx);
	dQdxVector.push_back((*it_hpm).dQdx);
	resRangeVector.push_back((*it_hpm).path_fraction*track_length);
	pitchVector.push_back((*it_hpm).pitch);
	XYZVector.push_back((*it_hpm).xyz);
	kinetic_energy += dEdxVector.back()*pitchVector.back();
      }

  }
  else{

    for(HitPropertiesMultiset_t::reverse_iterator it_hpm=hpm.rbegin();
	it_hpm!=hpm.rend();
	it_hpm++)
      {
	dEdxVector.push_back((*it_hpm).dEdx);
	dQdxVector.push_back((*it_hpm).dQdx);
	resRangeVector.push_back((1-(*it_hpm).path_fraction)*track_length);
	pitchVector.push_back((*it_hpm).pitch);
	XYZVector.push_back((*it_hpm).xyz);
	kinetic_energy += dEdxVector.back()*pitchVector.back();
      }

  }

  caloVector.emplace_back(kinetic_energy,
			  dEdxVector,
			  dQdxVector,
			  resRangeVector,
			  deadWireVector,
			  track_length,
			  pitchVector,
			  recob::tracking::convertCollToPoint(XYZVector),
			  planeID);
  assnTrackCaloVector.emplace_back(i_track);
}

void calo::TrackCalorimetryAlg::PrintHitPropertiesMultiset(HitPropertiesMultiset_t const& hpm){

  for(auto const& hit : hpm)
    hit.Print();

  std::cout << "Inverted? " << IsInvertedTrack(hpm) << std::endl;

}
