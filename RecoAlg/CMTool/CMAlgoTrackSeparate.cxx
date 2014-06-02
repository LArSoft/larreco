#ifndef CMALGOTRACKSEPARATE_CXX
#define CMALGOTRACKSEPARATE_CXX

#include "CMAlgoTrackSeparate.h"

namespace cluster {

  //----------------------------------------
  CMAlgoTrackSeparate::CMAlgoTrackSeparate() : CBoolAlgoBase()
  //----------------------------------------
  {

    //this just sets default values    
    SetVerbose(true);
    SetDebug(true);
    
    //1e9 is huge; everything will be merged
    SetMinNumHits(30);
    SetMinAngleDiff(15.); //in degrees
    SetMaxOpeningAngle(12.0); //in deg (parameter in rad!!)
    SetMinLength(10.);
    SetMinPolyHitDensity(10.);
    SetMaxWidth(10.);

  }

  //--------------------------------------------------------
  bool CMAlgoTrackSeparate::Bool(const ClusterParamsAlg &cluster1,
				 const ClusterParamsAlg &cluster2)
  //--------------------------------------------------------
  {

    //two clusters are considered un-mergable if:
    //1) both have more than _MinNumHits
    //2) opening angle for both < _MAxOpeningAngle
    //3) diff. in direction of both < _MinAngleDiff

    int N_Hits1 = cluster1.GetParams().N_Hits;
    int N_Hits2 = cluster2.GetParams().N_Hits;
    util::PxPoint start_point1 = cluster1.GetParams().start_point;
    util::PxPoint start_point2 = cluster2.GetParams().start_point;
    double angle_2d1 = cluster1.GetParams().angle_2d;
    double angle_2d2 = cluster2.GetParams().angle_2d;
    double opening_angle1 = cluster1.GetParams().opening_angle;
    double opening_angle2 = cluster2.GetParams().opening_angle;
    Polygon2D PolyObject1 = cluster1.GetParams().PolyObject;
    Polygon2D PolyObject2 = cluster2.GetParams().PolyObject;
    double length1 = cluster1.GetParams().length;
    double length2 = cluster2.GetParams().length;
    double width1 = cluster1.GetParams().width;
    double width2 = cluster2.GetParams().width;

    //first filter out low hits clusters
    if ( (N_Hits1 > _MinNumHits) and
	 (N_Hits2 > _MinNumHits) ) {
      if (_debug) {
	std::cout << "Cluster1 Num Hits: " << N_Hits1 << std::endl;
	std::cout << "\t Start: (" << start_point1.w << " " << start_point1.t << " )" << std::endl;
	std::cout << "\t Opening ANgle " << opening_angle1*(360/(2*3.14)) << std::endl;
	std::cout << "\t Angle2D: " << angle_2d1 << std::endl;
	std::cout << "\t Length: " << length1 << std::endl;
	std::cout << "\t Width: " << width1 << std::endl;
	std::cout << "Cluster2 Num Hits: " << N_Hits2 << std::endl;
	std::cout << "\t Start: (" << start_point2.w  << " " << start_point2.t << " )" << std::endl;
	std::cout << "\t Opening ANgle " << opening_angle2*(360/(2*3.14)) << std::endl;
	std::cout << "\t Angle2D: " << angle_2d2 << std::endl;
	std::cout << "\t Length: " << length2 << std::endl;
	std::cout << "\t Width: " << width2 << std::endl;
      }
      if ( (N_Hits1 > _MinNumHits) and
	   (N_Hits2 > _MinNumHits) and
	   ( abs(angle_2d1 - angle_2d2) > _MinAngleDiff ) and
	   //( PolyObject1.Area()/N_Hits1 > _MinDensity ) and
	   //( PolyObject2.Area()/N_Hits2 > _MinDensity ) and
	   (opening_angle1 < _MaxOpeningAngle/(360/(2*3.14))) and
	   (opening_angle2 < _MaxOpeningAngle/(360/(2*3.14))) and
	   (width1 < _MaxWidth) and
	   (width2 < _MaxWidth) and
	   (length1 > _MinLength) and
	   (length2 > _MinLength) ){
	if (_verbose) { std::cout << "*****************************************Separate with TrackSeparate!" << std::endl; }
	return true;
      }
    }
    
    return false;

  }

  //-----------------------
  void CMAlgoTrackSeparate::Report()
  //-----------------------
  {
  }

}

#endif
