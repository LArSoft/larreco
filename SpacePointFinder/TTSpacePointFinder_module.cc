#ifndef TTSPACEPOINTFINDER_H
#define TTSPACEPOINTFINDER_H
/*!
 * Title:   TTSpacePointFinder class
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Hit
 * Outputs: recob::SpacePoint
 *
 * Description:
 * This module, TimeTickSpacePointFinder (or TTSpacePointFinder for short) is 
 * designed to produce a spacepoint object based on hits from TTHitFinder.
 * There is intention to allow for a significant number of ghost spacepoints, 
 * with some downstream application dealing with the results.
 *
 * This code is totally microboone specific, btw.
 */
#include <string>
#include <math.h>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Core/EDProducer.h" 

// LArSoft Includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"

//boost includes
#include "boost/multi_array.hpp"

bool  HitTimeComparison(art::Ptr<recob::Hit> a, art::Ptr<recob::Hit> b) { return a->PeakTime() < b->PeakTime(); }

namespace sppt{

  class TTSpacePointFinder : public art::EDProducer {
    
  public:
    
    explicit TTSpacePointFinder(fhicl::ParameterSet const& pset); 
    virtual ~TTSpacePointFinder();
         
    void produce(art::Event& evt); 
    void beginJob();
    void beginRun(art::Run& run);
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                

  private:
        
    std::string    fHitModuleLabel;     /// Input hit module name
    std::string    fUHitsInstanceLabel; /// Input U hits instance name
    std::string    fVHitsInstanceLabel; /// Input V hits instance name
    std::string    fYHitsInstanceLabel; /// Input Y hits instance name
    float          fTimeDiffMax;        /// Maximum allowed time difference
    float          fYDiffMax;           /// Maximum allowed y-coordinate difference
    float          fZDiffMax;           /// Maximum allowed z-coordinate difference

    void  sortHitsByTime( std::vector< art::Ptr<recob::Hit> > &hits_handle);
   
    double TIME_OFFSET_U;
    double TIME_OFFSET_V;
    double TIME_OFFSET_Y;
    double TICKS_TO_X;

    boost::multi_array<double, 2> coordinates_UV_y;//(boost::extents[1][1]);
    boost::multi_array<double, 2> coordinates_UV_z;//(boost::extents[1][1]);
    boost::multi_array<double, 2> coordinates_UY_y;//(boost::extents[1][1]);
    boost::multi_array<double, 2> coordinates_UY_z;//(boost::extents[1][1]);

  protected:     
  
  }; // class TTSpacePointFinder  
  
  //-------------------------------------------------
  TTSpacePointFinder::TTSpacePointFinder(fhicl::ParameterSet const& pset) {
    this->reconfigure(pset);
    produces< std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::SpacePoint, recob::Hit>       >();
  }

  //-------------------------------------------------
  TTSpacePointFinder::~TTSpacePointFinder(){}

  //-------------------------------------------------
  void TTSpacePointFinder::reconfigure(fhicl::ParameterSet const& p) {
    fHitModuleLabel     = p.get< std::string >("HitModuleLabel","tthit");
    fUHitsInstanceLabel = p.get< std::string >("UHitsInstaceLabel","uhits");
    fVHitsInstanceLabel = p.get< std::string >("VHitsInstaceLabel","vhits");
    fYHitsInstanceLabel = p.get< std::string >("YHitsInstaceLabel","yhits");
    fTimeDiffMax        = p.get< float       >("TimeDiffMax");
    fZDiffMax           = p.get< float       >("ZDiffMax");
    fYDiffMax           = p.get< float       >("YDiffMax");

    //enforce a minimum time diff
    if(fTimeDiffMax<0){
      mf::LogError("TTSpacePointFinder") << "Time difference must be greater than zero.";
      return;
    }
    if(fZDiffMax<0){
      mf::LogError("TTSpacePointFinder") << "Z-coordinate difference must be greater than zero.";
      return;
    }
    if(fYDiffMax<0){
      mf::LogError("TTSpacePointFinder") << "Y-coordinate difference must be greater than zero.";
      return;
    }

  }

  //-------------------------------------------------
  void TTSpacePointFinder::beginJob(){}

  //-------------------------------------------------
  void TTSpacePointFinder::beginRun(art::Run &run){
        art::ServiceHandle<util::DetectorProperties> detprop;

    TIME_OFFSET_U = -1*detprop->GetXTicksOffset(geo::View_t::kU,0,0);
    TIME_OFFSET_V = -1*detprop->GetXTicksOffset(geo::View_t::kV,0,0);
    TIME_OFFSET_Y = -1*detprop->GetXTicksOffset(geo::View_t::kZ,0,0);

    TICKS_TO_X = detprop->GetXTicksCoefficient();
    
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int nwires_u = geom->Nwires(geo::View_t::kU);
    unsigned int nwires_v = geom->Nwires(geo::View_t::kV);
    unsigned int nwires_y = geom->Nwires(geo::View_t::kZ);
    
    coordinates_UV_y.resize(boost::extents[nwires_v][nwires_u]);
    coordinates_UV_z.resize(boost::extents[nwires_v][nwires_u]);
    coordinates_UY_y.resize(boost::extents[nwires_y][nwires_u]);
    coordinates_UY_z.resize(boost::extents[nwires_y][nwires_u]);
    for(unsigned int iu=0; iu<nwires_u; iu++){
      for(unsigned int iv=0; iv<nwires_v; iv++){
	geom->IntersectionPoint(iu,iv,
				geo::View_t::kU,geo::View_t::kV,
				0,0,
				coordinates_UV_y[iv][iu], //y
				coordinates_UV_z[iv][iu]); //z
      }
      for(unsigned int iy=0; iy<nwires_y; iy++){
	geom->IntersectionPoint(iu,iy,
				geo::View_t::kU,geo::View_t::kZ,
				0,0,
				coordinates_UY_y[iy][iu],
				coordinates_UY_z[iy][iu]);
      }
    }	    
  }

  //-------------------------------------------------
  void TTSpacePointFinder::endJob(){}

  //-------------------------------------------------
  void TTSpacePointFinder::produce(art::Event& evt)
  { 

    //initialize our spacepoint collection
    std::unique_ptr<std::vector<recob::SpacePoint> > spptCollection(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::SpacePoint,recob::Hit> > spptAssns(new art::Assns<recob::SpacePoint,recob::Hit>);

    // Read in the hits. Note, we will reorder hit vector, so we do in fact need a copy.
    art::Handle< std::vector<recob::Hit> > hitHandle_U;
    evt.getByLabel(fHitModuleLabel,fUHitsInstanceLabel,hitHandle_U);
    std::vector< art::Ptr<recob::Hit> > hitVec_U;
    art::fill_ptr_vector(hitVec_U,hitHandle_U);

    art::Handle< std::vector<recob::Hit> > hitHandle_V;
    evt.getByLabel(fHitModuleLabel,fVHitsInstanceLabel,hitHandle_V);
    std::vector< art::Ptr<recob::Hit> > hitVec_V;
    art::fill_ptr_vector(hitVec_V,hitHandle_V);

    art::Handle< std::vector<recob::Hit> > hitHandle_Y;
    evt.getByLabel(fHitModuleLabel,fYHitsInstanceLabel,hitHandle_Y);
    std::vector< art::Ptr<recob::Hit> > hitVec_Y;
    art::fill_ptr_vector(hitVec_Y,hitHandle_Y);

    mf::LogInfo("TTSpacePointFinderDetail") 
      << "Got handles to hits."; 
    
    //sort the hits by the time
    sortHitsByTime(hitVec_U);
    sortHitsByTime(hitVec_V);
    sortHitsByTime(hitVec_Y);

    mf::LogInfo("TTSpacePointFinderDetail") 
      << "Sorted " 
      << hitVec_U.size() << " u hits, "
      << hitVec_V.size() << " v hits, "
      << hitVec_Y.size() << " y hits.";

    //now, do the loop to search for like-timed hits across the three planes
    std::vector< art::Ptr<recob::Hit> >::iterator ihitu = hitVec_U.begin();
    std::vector< art::Ptr<recob::Hit> >::iterator ihitv = hitVec_V.begin();
    std::vector< art::Ptr<recob::Hit> >::iterator ihity = hitVec_Y.begin();
    std::vector< art::Ptr<recob::Hit> >::iterator ihitv_inner,ihity_inner;
    double time_hitu = (*ihitu)->PeakTime() + TIME_OFFSET_U;
    double time_hitv = (*ihitv)->PeakTime() + TIME_OFFSET_V;
    double time_hity = (*ihity)->PeakTime() + TIME_OFFSET_Y;
    double time_hitv_inner,time_hity_inner;
    while(ihitu != hitVec_U.end()){
      time_hitu = (*ihitu)->PeakTime() + TIME_OFFSET_U;

      mf::LogInfo("TTSpacePointFinderDetail") 
	<< "Hit times (u,v,y)=("
	<< time_hitu << ","
	<< time_hitv << ","
	<< time_hity << ")"; 

      //if time_hitu is too much bigger than time_hitv, need to advance hitv iterator
      while( (time_hitu-time_hitv)>fTimeDiffMax ){
	ihitv++;
	if(ihitv==hitVec_V.end()) break;
	time_hitv = (*ihitv)->PeakTime() + TIME_OFFSET_V;
      }
      if(ihitv==hitVec_V.end()) break;
      
      //same thing with time_hitu and time_hity
      while( (time_hitu-time_hity)>fTimeDiffMax ){
	ihity++;
	if(ihity==hitVec_Y.end()) break;
	time_hity = (*ihity)->PeakTime() + TIME_OFFSET_Y;
      }
      if(ihity==hitVec_Y.end()) break;

      //OK, now we know time_hitu <= time_hitv and time_hitu <= time_hity.
      //Next, check if time_hitu is near time_hitv and time_hit y. If not, 
      //we have to increment ihitu.
      if(std::abs(time_hitu-time_hitv)>fTimeDiffMax || std::abs(time_hitu-time_hity)>fTimeDiffMax){
	ihitu++;
	continue;
      }

      //OK! Note we KNOW that these three match in time: 
      //  -- time_hitu is within fTimeDiffMax of both time_hitv and time_hity; and
      //  -- time_hitu <= time_hitv AND time_hitu <=time_hity, so time_hitv and time_hity are near too
      
      mf::LogInfo("TTSpacePointFinderDetail") 
	<< "Matching hit times (u,v,y)=("
	<< time_hitu << ","
	<< time_hitv << ","
	<< time_hity << ")"; 

      //Next thing to do, we need to loop over all possible 3-hit matches for our given u-hit.
      //We need new iterators in v and y at this location, and will loop over those
      ihitv_inner = ihitv;
      time_hitv_inner = (*ihitv_inner)->PeakTime() + TIME_OFFSET_V;
      ihity_inner = ihity;
      time_hity_inner = (*ihity_inner)->PeakTime() + TIME_OFFSET_Y;

      while(std::abs(time_hitu-time_hitv_inner)<fTimeDiffMax && std::abs(time_hitu-time_hity_inner)<fTimeDiffMax){

	unsigned int uwire = (*ihitu)->WireID().Wire;
	unsigned int vwire = (*ihitv_inner)->WireID().Wire;
	unsigned int ywire = (*ihity_inner)->WireID().Wire;

	mf::LogInfo("TTSpacePointFinderDetail") 
	  << "(y,z) coordinate for uv/uy: ("
	  << coordinates_UV_y[vwire][uwire] << ","
	  << coordinates_UV_z[vwire][uwire] << ")/("
	  << coordinates_UY_y[ywire][uwire] << ","
	  << coordinates_UY_z[ywire][uwire] << ")";

	if(std::abs(coordinates_UV_y[vwire][uwire]-coordinates_UY_y[ywire][uwire])<fYDiffMax &&
	   std::abs(coordinates_UV_z[vwire][uwire]-coordinates_UY_z[ywire][uwire])<fZDiffMax){

	  double xyz[3];
	  double xyz_err[6];
	  //triangular error matrix:
	  // | 0.         |
	  // | 0.  0.     |
	  // | 0.  0.  0. |

	  //get average y and z, with errors
	  xyz[1] = (coordinates_UV_y[vwire][uwire] + coordinates_UY_y[ywire][uwire])*0.5;
	  xyz_err[2] = std::abs(coordinates_UV_y[vwire][uwire] - xyz[1]);
	  xyz[2] = (coordinates_UV_z[vwire][uwire] + coordinates_UY_z[ywire][uwire])*0.5;
	  xyz_err[5] = std::abs(coordinates_UV_z[vwire][uwire] - xyz[2]);

	  double t_val = (time_hitu + time_hitv_inner + time_hity_inner)/3.;
	  double t_err = 0.5*std::sqrt( (time_hitu-t_val)*(time_hitu-t_val) +
					(time_hitv_inner-t_val)*(time_hitv_inner-t_val) +
					(time_hity_inner-t_val)*(time_hity_inner-t_val));
	  xyz[0] = TICKS_TO_X * t_val;
	  xyz_err[0] = TICKS_TO_X * t_err;
	  
	  //make space point to put on event
	  recob::SpacePoint spt(xyz, xyz_err, 0., spptCollection->size());
	  spptCollection->push_back(spt);

	  //make association with hits
	  std::vector< art::Ptr<recob::Hit> > myhits = {*ihitu,*ihitv_inner,*ihity_inner};
	  util::CreateAssn(*this, evt, *spptCollection, myhits, *spptAssns);
	}

	//now increment the v or y hit, whichever is smalles (closest to u hit) in time
	if(time_hitv_inner <= time_hity_inner){
	  ihitv_inner++;
	  if(ihitv_inner==hitVec_V.end()) break;
	  time_hitv_inner = (*ihitv_inner)->PeakTime() + TIME_OFFSET_V;
	}
	else{
	  ihity_inner++;
	  if(ihity_inner==hitVec_Y.end()) break;
	  time_hity_inner = (*ihity_inner)->PeakTime() + TIME_OFFSET_Y;
	}

      }

      ihitu++;
    }
    
    //finally, put things on the event
    evt.put(std::move(spptCollection));
    evt.put(std::move(spptAssns));

  } // End of produce()  
    
  //-------------------------------------------------
  void TTSpacePointFinder::sortHitsByTime(std::vector< art::Ptr<recob::Hit> > &hitVec){
    std::sort(hitVec.begin(),hitVec.end(),HitTimeComparison);
  }
  

  DEFINE_ART_MODULE(TTSpacePointFinder)

} // end of hit namespace


#endif // TTSPACEPOINTFINDER_H
