////////////////////////////////////////////////////////////////////////
//
// CornerFinderAlg class
//
// wketchum@fnal.gov
//
// CornerFinder is meant to use image-processing techniques (mainly Harris-Stephens
// corner-finding) to find "corners" using the information from calibrated wires.
//
//  Conversion_algorithm options:
//     standard --- basically a copy of the calibrated wires
//     skeleton --- a thinned copy of the calibrated wires
//     binary   --- ticks above threshold get assigned a value 10*threshold, everything else = threshold
//     function --- apply a function (like a double-Gaussian) to a neighborhood around each tick
//
//  Derivative options:
//     Sobel --- apply a Sobel mask (neighborhood of 1 or 2 supported)
//     local --- take slope from immediately neighboring bins (neighborhood of 1 supported)
//
//  Derivative_BlurFunc options:  none. You're stuck with a double gaussian.
//
//  CornerScore_algorithm options:
//     Noble  --- determinant / (trace + Noble_epsilon)
//     Harris --- determinant - (trace)^2 * Harris_kappa
////////////////////////////////////////////////////////////////////////


#include "larreco/RecoAlg/CornerFinderAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

// NOTE: In the .h file I assumed this would belong in the cluster class....if
// we decide otherwise we will need to search and replace for this


//-----------------------------------------------------------------------------
corner::CornerFinderAlg::CornerFinderAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

corner::CornerFinderAlg::~CornerFinderAlg(){
  this->CleanCornerFinderAlg();
}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::CleanCornerFinderAlg()
{

  //for (auto wd_histo : WireData_histos) delete wd_histo;
  //for (auto wd_histo : WireData_histos_ProjectionX) delete wd_histo;
  //for (auto wd_histo : WireData_histos_ProjectionY) delete wd_histo;
  //for (auto histo : fConversion_histos) delete histo;
  //for (auto histo : fDerivativeX_histos) delete histo;
  //for (auto histo : fDerivativeY_histos) delete histo;
  //for (auto histo : fCornerScore_histos) delete histo;
  //for (auto histo : fMaxSuppress_histos) delete histo;

  WireData_histos.clear();
  WireData_histos_ProjectionX.clear();
  WireData_histos_ProjectionY.clear();
  WireData_IDs.clear();

  //for (auto wd_histo : WireData_trimmed_histos) delete std::get<1>(wd_histo);
  WireData_trimmed_histos.clear();

}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::reconfigure(fhicl::ParameterSet const& p)
{
  // ### These are all the tuneable .fcl file parameters from the event ###
  fCalDataModuleLabel  			 = p.get< std::string 	 >("CalDataModuleLabel");
  fTrimming_threshold     		 = p.get< float    	 >("Trimming_threshold");
  fTrimming_totalThreshold                = p.get< double         >("Trimming_totalThreshold");
  fConversion_threshold     		 = p.get< float    	 >("Conversion_threshold");
  fConversion_bins_per_input_x  	 = p.get< int      	 >("Conversion_bins_per_input_x");
  fConversion_bins_per_input_y       	 = p.get< int      	 >("Conversion_bins_per_input_y");
  fConversion_algorithm                  = p.get< std::string    >("Conversion_algorithm");
  fConversion_func                       = p.get< std::string    >("Conversion_function");
  fConversion_func_neighborhood     	 = p.get< int		 >("Conversion_func_neighborhood");
  fDerivative_method        		 = p.get< std::string    >("Derivative_method");
  fDerivative_neighborhood     	         = p.get< int		 >("Derivative_neighborhood");
  fDerivative_BlurFunc        		 = p.get< std::string    >("Derivative_BlurFunc");
  fDerivative_BlurNeighborhood           = p.get< int		 >("Derivative_BlurNeighborhood");
  fCornerScore_neighborhood     	 = p.get< int		 >("CornerScore_neighborhood");
  fCornerScore_algorithm		 = p.get< std::string    >("CornerScore_algorithm");
  fCornerScore_Noble_epsilon		 = p.get< float          >("CornerScore_Noble_epsilon");
  fCornerScore_Harris_kappa		 = p.get< float          >("CornerScore_Harris_kappa");
  fMaxSuppress_neighborhood		 = p.get< int		 >("MaxSuppress_neighborhood");
  fMaxSuppress_threshold		 = p.get< int		 >("MaxSuppress_threshold");
  fIntegral_bin_threshold                = p.get< float          >("Integral_bin_threshold");
  fIntegral_fraction_threshold           = p.get< float          >("Integral_fraction_threshold");

  int neighborhoods[] = { fConversion_func_neighborhood,
			  fDerivative_neighborhood,
			  fDerivative_BlurNeighborhood,
			  fCornerScore_neighborhood,
			  fMaxSuppress_neighborhood };
  fTrimming_buffer = *std::max_element(neighborhoods,neighborhoods+5);

}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::InitializeGeometry(geo::Geometry const& my_geometry){

  CleanCornerFinderAlg();

  // set the sizes of the WireData_histos and WireData_IDs
  unsigned int nPlanes = my_geometry.Nplanes();
  WireData_histos.resize(nPlanes);
  WireData_histos_ProjectionX.resize(nPlanes);
  WireData_histos_ProjectionY.resize(nPlanes);
  //  fConversion_histos.resize(nPlanes);
  //  fDerivativeX_histos.resize(nPlanes);
  //  fDerivativeY_histos.resize(nPlanes);
  //  fCornerScore_histos.resize(nPlanes);
  //  fMaxSuppress_histos.resize(nPlanes);

  /* For now, we need something to associate each wire in the histogram with a wire_id.
     This is not a beautiful way of handling this, but for now it should work. */
  WireData_IDs.resize(nPlanes);
  for(unsigned int i_plane=0; i_plane < nPlanes; ++i_plane)
    WireData_IDs.at(i_plane).resize(my_geometry.Nwires(i_plane));

  WireData_trimmed_histos.resize(0);

}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::GrabWires( std::vector<recob::Wire> const& wireVec, geo::Geometry const& my_geometry){

  InitializeGeometry(my_geometry);

  const unsigned int nTimeTicks = wireVec.at(0).NSignal();

  // Initialize the histograms.
  // All of this should eventually be changed to not need to use histograms...
  for (unsigned int i_plane=0; i_plane < my_geometry.Nplanes(); i_plane++){

    std::stringstream ss_tmp_name,ss_tmp_title;
    ss_tmp_name << "h_WireData_" << i_plane;
    ss_tmp_title << fCalDataModuleLabel << " wire data for plane " << i_plane << ";Wire Number;Time Tick";

    if( (unsigned int)(WireData_histos.at(i_plane).GetNbinsX()) == (my_geometry.Nwires(i_plane)) ) {
      WireData_histos.at(i_plane).Reset();
      WireData_histos.at(i_plane).SetName(ss_tmp_name.str().c_str());
      WireData_histos.at(i_plane).SetTitle(ss_tmp_title.str().c_str());
    }
    else
      WireData_histos.at(i_plane) = TH2F(ss_tmp_name.str().c_str(),
					 ss_tmp_title.str().c_str(),
					 my_geometry.Nwires(i_plane),
					 0,
					 my_geometry.Nwires(i_plane),
					 nTimeTicks,
					 0,
					 nTimeTicks);
  }


  /* Now do the loop over the wires. */
  for (std::vector<recob::Wire>::const_iterator iwire = wireVec.begin(); iwire < wireVec.end(); iwire++) {

    std::vector<geo::WireID> possible_wireIDs = my_geometry.ChannelToWire(iwire->Channel());
    geo::WireID this_wireID;
    try { this_wireID = possible_wireIDs.at(0);}
    catch(cet::exception& excep) { mf::LogError("CornerFinderAlg") << "Bail out! No Possible Wires!\n"; }

    unsigned int i_plane = this_wireID.Plane;
    unsigned int i_wire = this_wireID.Wire;

    WireData_IDs.at(i_plane).at(i_wire) = this_wireID;

    for(unsigned int i_time = 0; i_time < nTimeTicks; i_time++){
      WireData_histos.at(i_plane).SetBinContent(i_wire,i_time,(iwire->Signal()).at(i_time));
    }//<---End time loop

  }//<-- End loop over wires


  for (unsigned int i_plane=0; i_plane < my_geometry.Nplanes(); i_plane++){
    WireData_histos_ProjectionX.at(i_plane) = *(WireData_histos.at(i_plane).ProjectionX());
    WireData_histos_ProjectionY.at(i_plane) = *(WireData_histos.at(i_plane).ProjectionY());
  }


}


//-----------------------------------------------------------------------------------
// This gives us a vecotr of EndPoint2D objects that correspond to possible corners
void corner::CornerFinderAlg::get_feature_points(std::vector<recob::EndPoint2D> & corner_vector,
						 geo::Geometry const& my_geometry){



  for(auto const& pid : my_geometry.IteratePlaneIDs()){
    attach_feature_points(WireData_histos.at(pid.Plane),
			  WireData_IDs.at(pid.Plane),
			  my_geometry.View(pid),
			  corner_vector);
  }

}

//-----------------------------------------------------------------------------------
// This gives us a vector of EndPoint2D objects that correspond to possible corners, but quickly!
void corner::CornerFinderAlg::get_feature_points_fast(std::vector<recob::EndPoint2D> & corner_vector,
						      geo::Geometry const& my_geometry){



  create_smaller_histos(my_geometry);

  for(unsigned int cstat = 0; cstat < my_geometry.Ncryostats(); ++cstat){
    for(unsigned int tpc = 0; tpc < my_geometry.Cryostat(cstat).NTPC(); ++tpc){
      for(size_t histos=0; histos!= WireData_trimmed_histos.size(); histos++){

	int plane = std::get<0>(WireData_trimmed_histos.at(histos));
	int startx = std::get<2>(WireData_trimmed_histos.at(histos));
	int starty = std::get<3>(WireData_trimmed_histos.at(histos));

	MF_LOG_DEBUG("CornerFinderAlg")
	  << "Doing histogram " << histos
	  << ", of plane " << plane
	  << " with start points " << startx << " " << starty;

	attach_feature_points(std::get<1>(WireData_trimmed_histos.at(histos)),
			      WireData_IDs.at(plane),my_geometry.Cryostat(cstat).TPC(tpc).Plane(plane).View(),corner_vector,startx,starty);

	MF_LOG_DEBUG("CornerFinderAlg") << "Total feature points now is " << corner_vector.size();
      }

      //remove_duplicates(corner_vector);

    }
  }

}

//-----------------------------------------------------------------------------------
// This gives us a vecotr of EndPoint2D objects that correspond to possible corners
// Uses line integral score as corner strength
void corner::CornerFinderAlg::get_feature_points_LineIntegralScore(std::vector<recob::EndPoint2D> & corner_vector,
								   geo::Geometry const& my_geometry){



  for(auto const& pid : my_geometry.IteratePlaneIDs()){
    attach_feature_points_LineIntegralScore(WireData_histos.at(pid.Plane),
					    WireData_IDs.at(pid.Plane),
					    my_geometry.View(pid),
					    corner_vector);
  }

}

//-----------------------------------------------------------------------------------
// void corner::CornerFinderAlg::remove_duplicates(std::vector<recob::EndPoint2D> & corner_vector){

//   int i_wire, j_wire;
//   float i_time, j_time;
//   for(size_t i=0; i != corner_vector.size(); i++){

//     i_wire = corner_vector.at(i).WireID().Wire;
//     i_time = corner_vector.at(i).DriftTime();

//     for(size_t j=i+1; j != corner_vector.size(); j++){

//       j_wire = corner_vector.at(j).WireID().Wire;
//       j_time = corner_vector.at(j).DriftTime();

//       if(std::abs(i_wire-j_wire) < 5 && std::abs(i_time - j_time) < 10){
//         corner_vector.erase(corner_vector.begin()+j);
//         j--;
//       }

//     }

//   }

// }

struct compare_to_value{

  compare_to_value(int b) {this->b = b;}
  bool operator() (int i, int j) {
    return std::abs(b-i)<std::abs(b-j);
  }

  int b;

};

struct compare_to_range{

  compare_to_range(int a, int b) {this->a = a; this->b = b;}
  bool operator() (int i, int j) {

    int mid = (b-a)/2 + a;
    if(i>=a && i<=b && j>=a && j<=b)
      return std::abs(mid-i)<std::abs(mid-j);

    else if(j>=a && j<=b && (i<a || i>b) )
      return false;

    else if(i>=a && i<=b && (j<a || j>b) )
      return true;

    else
      return true;
  }

  int a;
  int b;

};

//-----------------------------------------------------------------------------
// This looks for areas of the wires that are non-noise, to speed up evaluation
void corner::CornerFinderAlg::create_smaller_histos(geo::Geometry const& my_geometry){

  for(auto const& pid : my_geometry.IteratePlaneIDs() ){

    MF_LOG_DEBUG("CornerFinderAlg")
      << "Working plane " << pid.Plane << ".";

    int x_bins = WireData_histos_ProjectionX.at(pid.Plane).GetNbinsX();
    int y_bins = WireData_histos_ProjectionY.at(pid.Plane).GetNbinsX();

    std::vector<int> cut_points_x {0};
    std::vector<int> cut_points_y {0};

    for (int ix=1; ix<=x_bins; ix++){

      float this_value = WireData_histos_ProjectionX.at(pid.Plane).GetBinContent(ix);

      if(ix<fTrimming_buffer || ix>(x_bins-fTrimming_buffer)) continue;

      int jx=ix-fTrimming_buffer;
      while(this_value<fTrimming_threshold){
	if(jx==ix+fTrimming_buffer) break;
	this_value = WireData_histos_ProjectionX.at(pid.Plane).GetBinContent(jx);
	jx++;
      }
      if(this_value<fTrimming_threshold){
	cut_points_x.push_back(ix);
      }

    }

    for (int iy=1; iy<=y_bins; iy++){

      float this_value = WireData_histos_ProjectionY.at(pid.Plane).GetBinContent(iy);

      if(iy<fTrimming_buffer || iy>(y_bins-fTrimming_buffer)) continue;

      int jy=iy-fTrimming_buffer;
      while(this_value<fTrimming_threshold){
	if(jy==iy+fTrimming_buffer) break;
	this_value = WireData_histos_ProjectionY.at(pid.Plane).GetBinContent(jy);
	jy++;
      }
      if(this_value<fTrimming_threshold){
	cut_points_y.push_back(iy);
      }

    }

    MF_LOG_DEBUG("CornerFinderAlg")
      << "We have a total of " << cut_points_x.size() << " x cut points."
      << "\nWe have a total of " << cut_points_y.size() << " y cut points.";

    std::vector<int> x_low{1};
    std::vector<int> x_high{x_bins};
    std::vector<int> y_low{1};
    std::vector<int> y_high{y_bins};
    bool x_change = true;
    bool y_change = true;
    while(x_change || y_change){

      x_change = false;
      y_change = false;

      size_t current_size = x_low.size();

      for(size_t il=0; il<current_size; il++){

	int comp_value = (x_high.at(il) + x_low.at(il)) / 2;
	std::sort(cut_points_x.begin(),cut_points_x.end(),compare_to_value(comp_value));

	if(cut_points_x.at(0) <= x_low.at(il) || cut_points_x.at(0) >= x_high.at(il))
	  continue;

	double integral_low = WireData_histos.at(pid.Plane).Integral(x_low.at(il),cut_points_x.at(0),y_low.at(il),y_high.at(il));
	double integral_high = WireData_histos.at(pid.Plane).Integral(cut_points_x.at(0),x_high.at(il),y_low.at(il),y_high.at(il));
	if(integral_low > fTrimming_totalThreshold && integral_high > fTrimming_totalThreshold){
	  x_low.push_back(cut_points_x.at(0));
	  x_high.push_back(x_high.at(il));
	  y_low.push_back(y_low.at(il));
	  y_high.push_back(y_high.at(il));

	  x_high[il] = cut_points_x.at(0);
	  x_change = true;
	}
	else if(integral_low > fTrimming_totalThreshold && integral_high < fTrimming_totalThreshold){
	  x_high[il] = cut_points_x.at(0);
	  x_change = true;
	}
	else if(integral_low < fTrimming_totalThreshold && integral_high > fTrimming_totalThreshold){
	  x_low[il] = cut_points_x.at(0);
	  x_change = true;
	}
      }

      current_size = x_low.size();

      for(size_t il=0; il<current_size; il++){

	int comp_value = (y_high.at(il) - y_low.at(il)) / 2;
	std::sort(cut_points_y.begin(),cut_points_y.end(),compare_to_value(comp_value));

	if(cut_points_y.at(0) <= y_low.at(il) || cut_points_y.at(0) >= y_high.at(il))
	  continue;

	double integral_low = WireData_histos.at(pid.Plane).Integral(x_low.at(il),x_high.at(il),y_low.at(il),cut_points_y.at(0));
	double integral_high = WireData_histos.at(pid.Plane).Integral(x_low.at(il),x_high.at(il),cut_points_y.at(0),y_high.at(il));
	if(integral_low > fTrimming_totalThreshold && integral_high > fTrimming_totalThreshold){
	  y_low.push_back(cut_points_y.at(0));
	  y_high.push_back(y_high.at(il));
	  x_low.push_back(x_low.at(il));
	  x_high.push_back(x_high.at(il));

	  y_high[il] = cut_points_y.at(0);
	  y_change = true;
	}
	else if(integral_low > fTrimming_totalThreshold && integral_high < fTrimming_totalThreshold){
	  y_high[il] = cut_points_y.at(0);
	  y_change = true;
	}
	else if(integral_low < fTrimming_totalThreshold && integral_high > fTrimming_totalThreshold){
	  y_low[il] = cut_points_y.at(0);
	  y_change = true;
	}
      }

    }

    MF_LOG_DEBUG("CornerFinderAlg")
      << "First point in x is " << cut_points_x.at(0);

    std::sort(cut_points_x.begin(),cut_points_x.end(),compare_to_value(x_bins/2));

    MF_LOG_DEBUG("CornerFinderAlg")
      << "Now the first point in x is " << cut_points_x.at(0);

    MF_LOG_DEBUG("CornerFinderAlg")
       << "First point in y is " << cut_points_y.at(0);

    std::sort(cut_points_y.begin(),cut_points_y.end(),compare_to_value(y_bins/2));

    MF_LOG_DEBUG("CornerFinderAlg")
      << "Now the first point in y is " << cut_points_y.at(0);

    MF_LOG_DEBUG("CornerFinderAlg")
      << "\nIntegral on the SW side is "
      << WireData_histos.at(pid.Plane).Integral(1,cut_points_x.at(0),1,cut_points_y.at(0))
      << "\nIntegral on the SE side is "
      << WireData_histos.at(pid.Plane).Integral(cut_points_x.at(0),x_bins,1,cut_points_y.at(0))
      << "\nIntegral on the NW side is "
      << WireData_histos.at(pid.Plane).Integral(1,cut_points_x.at(0),cut_points_y.at(0),y_bins)
      << "\nIntegral on the NE side is "
      << WireData_histos.at(pid.Plane).Integral(cut_points_x.at(0),x_bins,cut_points_y.at(0),y_bins);


    for(size_t il=0; il<x_low.size(); il++){

      std::stringstream h_name;
      h_name << "h_" << pid.Cryostat << "_" << pid.TPC << "_" << pid.Plane << "_sub" << il;
      TH2F h_tmp( (h_name.str()).c_str(),"",
		  x_high.at(il)-x_low.at(il)+1,x_low.at(il),x_high.at(il),
		  y_high.at(il)-y_low.at(il)+1,y_low.at(il),y_high.at(il));

      for(int ix=1; ix<=(x_high.at(il)-x_low.at(il)+1); ix++){
	for(int iy=1; iy<=(y_high.at(il)-y_low.at(il)+1); iy++){
	  h_tmp.SetBinContent(ix,iy,WireData_histos.at(pid.Plane).GetBinContent(x_low.at(il)+(ix-1),y_low.at(il)+(iy-1)));
	}
      }

      WireData_trimmed_histos.push_back(std::make_tuple(pid.Plane,h_tmp,x_low.at(il)-1,y_low.at(il)-1));
    }

  }// end loop over PlaneIDs


}

//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void corner::CornerFinderAlg::attach_feature_points( TH2F const& h_wire_data,
						     std::vector<geo::WireID> wireIDs,
						     geo::View_t view,
						     std::vector<recob::EndPoint2D> & corner_vector,
						     int startx,
						     int starty){


  const int x_bins = h_wire_data.GetNbinsX();
  const float x_min = h_wire_data.GetXaxis()->GetBinLowEdge(1);
  const float x_max = h_wire_data.GetXaxis()->GetBinUpEdge(x_bins);

  const int y_bins = h_wire_data.GetNbinsY();
  const float y_min = h_wire_data.GetYaxis()->GetBinLowEdge(1);
  const float y_max = h_wire_data.GetYaxis()->GetBinUpEdge(y_bins);

  const int converted_y_bins = y_bins/fConversion_bins_per_input_y;
  const int converted_x_bins = x_bins/fConversion_bins_per_input_x;

  std::stringstream conversion_name;  conversion_name  << "h_conversion_"   << view << "_" << run_number << "_" << event_number;
  std::stringstream dx_name;          dx_name          << "h_derivative_x_" << view << "_" << run_number << "_" << event_number;
  std::stringstream dy_name;          dy_name          << "h_derivative_y_" << view << "_" << run_number << "_" << event_number;
  std::stringstream cornerScore_name; cornerScore_name << "h_cornerScore_"  << view << "_" << run_number << "_" << event_number;
  std::stringstream maxSuppress_name; maxSuppress_name << "h_maxSuppress_"  << view << "_" << run_number << "_" << event_number;

  TH2F conversion_histo(conversion_name.str().c_str(),"Image Conversion Histogram",
			 converted_x_bins,x_min,x_max,
			 converted_y_bins,y_min,y_max);

  TH2F derivativeX_histo(dx_name.str().c_str(),"Partial Derivatives (x)",
			  converted_x_bins,x_min,x_max,
			  converted_y_bins,y_min,y_max);

  TH2F derivativeY_histo(dy_name.str().c_str(),"Partial Derivatives (y)",
			  converted_x_bins,x_min,x_max,
			  converted_y_bins,y_min,y_max);

  TH2D cornerScore_histo(cornerScore_name.str().c_str(),"Corner Score",
			 converted_x_bins,x_min,x_max,
			 converted_y_bins,y_min,y_max);

  TH2D maxSuppress_histo(maxSuppress_name.str().c_str(),"Corner Points (Maximum Suppressed)",
			 converted_x_bins,x_min,x_max,
			 converted_y_bins,y_min,y_max);

  create_image_histo(h_wire_data,conversion_histo);
  create_derivative_histograms(conversion_histo,derivativeX_histo,derivativeY_histo);
  create_cornerScore_histogram(derivativeX_histo,derivativeY_histo,cornerScore_histo);
  perform_maximum_suppression(cornerScore_histo,corner_vector,wireIDs,view,maxSuppress_histo,startx,starty);
}


//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void corner::CornerFinderAlg::attach_feature_points_LineIntegralScore(TH2F const& h_wire_data,
								       std::vector<geo::WireID> wireIDs,
								       geo::View_t view,
								       std::vector<recob::EndPoint2D> & corner_vector){


  const int   x_bins = h_wire_data.GetNbinsX();
  const float x_min  = h_wire_data.GetXaxis()->GetBinLowEdge(1);
  const float x_max  = h_wire_data.GetXaxis()->GetBinUpEdge(x_bins);

  const int   y_bins = h_wire_data.GetNbinsY();
  const float y_min  = h_wire_data.GetYaxis()->GetBinLowEdge(1);
  const float y_max  = h_wire_data.GetYaxis()->GetBinUpEdge(y_bins);

  const int converted_y_bins = y_bins/fConversion_bins_per_input_y;
  const int converted_x_bins = x_bins/fConversion_bins_per_input_x;

  std::stringstream conversion_name;  conversion_name  << "h_conversion_"   << view << "_" << run_number << "_" << event_number;
  std::stringstream dx_name;          dx_name          << "h_derivative_x_" << view << "_" << run_number << "_" << event_number;
  std::stringstream dy_name;          dy_name          << "h_derivative_y_" << view << "_" << run_number << "_" << event_number;
  std::stringstream cornerScore_name; cornerScore_name << "h_cornerScore_"  << view << "_" << run_number << "_" << event_number;
  std::stringstream maxSuppress_name; maxSuppress_name << "h_maxSuppress_"  << view << "_" << run_number << "_" << event_number;

  TH2F h_conversion  ((conversion_name.str()).c_str(),
		      "Image Conversion Histogram",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2F h_derivative_x((dx_name.str()).c_str(),
		      "Partial Derivatives (x)",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2F h_derivative_y((dy_name.str()).c_str(),
		      "Partial Derivatives (y)",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2D h_cornerScore ((cornerScore_name.str()).c_str(),
		      "Feature Point Corner Score",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2D h_maxSuppress ((maxSuppress_name.str()).c_str(),
		      "Corner Points (Maximum Suppressed)",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);

  create_image_histo(h_wire_data,h_conversion);
  create_derivative_histograms(h_conversion,h_derivative_x,h_derivative_y);
  create_cornerScore_histogram(h_derivative_x,h_derivative_y,h_cornerScore);

  std::vector<recob::EndPoint2D> corner_vector_tmp;
  perform_maximum_suppression(h_cornerScore,corner_vector_tmp,wireIDs,view,h_maxSuppress);

  std::stringstream LI_name; LI_name << "h_lineIntegralScore_" << view << "_" << run_number << "_" << event_number;
  TH2F h_lineIntegralScore((LI_name.str()).c_str(),
			   "Line Integral Score",
			   x_bins,x_min,x_max,
			   y_bins,y_min,y_max);
  calculate_line_integral_score(h_wire_data,corner_vector_tmp,corner_vector,h_lineIntegralScore);

}


//-----------------------------------------------------------------------------
// Convert to pixel
void corner::CornerFinderAlg::create_image_histo(TH2F const& h_wire_data, TH2F & h_conversion) const {

  double temp_integral=0;

  const TF2 fConversion_TF2("fConversion_func",fConversion_func.c_str(),-20,20,-20,20);

  for(int ix=1; ix<=h_conversion.GetNbinsX(); ix++){
    for(int iy=1; iy<=h_conversion.GetNbinsY(); iy++){

      temp_integral = h_wire_data.Integral(ix,ix,iy,iy);

      if( temp_integral > fConversion_threshold){

	if(fConversion_algorithm.compare("binary")==0)
	  h_conversion.SetBinContent(ix,iy,10*fConversion_threshold);
	else if(fConversion_algorithm.compare("standard")==0)
	  h_conversion.SetBinContent(ix,iy,temp_integral);

	else if(fConversion_algorithm.compare("function")==0){

	  temp_integral = 0;
	  for(int jx=ix-fConversion_func_neighborhood; jx<=ix+fConversion_func_neighborhood; jx++){
	    for(int jy=iy-fConversion_func_neighborhood; jy<=iy+fConversion_func_neighborhood; jy++){
	      temp_integral += h_wire_data.GetBinContent(jx,jy)*fConversion_TF2.Eval(ix-jx,iy-jy);
	    }
	  }
	  h_conversion.SetBinContent(ix,iy,temp_integral);
	}

	else if(fConversion_algorithm.compare("skeleton")==0){

	  if( (temp_integral > h_wire_data.GetBinContent(ix-1,iy) && temp_integral > h_wire_data.GetBinContent(ix+1,iy))
	      || (temp_integral > h_wire_data.GetBinContent(ix,iy-1) && temp_integral > h_wire_data.GetBinContent(ix,iy+1)))
	    h_conversion.SetBinContent(ix,iy,temp_integral);
	  else
	    h_conversion.SetBinContent(ix,iy,fConversion_threshold);
	}
	else if(fConversion_algorithm.compare("sk_bin")==0){

	  if( (temp_integral > h_wire_data.GetBinContent(ix-1,iy) && temp_integral > h_wire_data.GetBinContent(ix+1,iy))
	      || (temp_integral > h_wire_data.GetBinContent(ix,iy-1) && temp_integral > h_wire_data.GetBinContent(ix,iy+1)))
	    h_conversion.SetBinContent(ix,iy,10*fConversion_threshold);
	  else
	    h_conversion.SetBinContent(ix,iy,fConversion_threshold);
	}
	else
	  h_conversion.SetBinContent(ix,iy,temp_integral);
      }

      else
	h_conversion.SetBinContent(ix,iy,fConversion_threshold);

    }
  }

}

//-----------------------------------------------------------------------------
// Derivative

void corner::CornerFinderAlg::create_derivative_histograms(TH2F const& h_conversion, TH2F & h_derivative_x, TH2F & h_derivative_y){

  const int x_bins = h_conversion.GetNbinsX();
  const int y_bins = h_conversion.GetNbinsY();

  for(int iy=1+fDerivative_neighborhood; iy<=(y_bins-fDerivative_neighborhood); iy++){
    for(int ix=1+fDerivative_neighborhood; ix<=(x_bins-fDerivative_neighborhood); ix++){

      if(fDerivative_method.compare("Sobel")==0){

	if(fDerivative_neighborhood==1){
	  h_derivative_x.SetBinContent(ix,iy,
					0.5*(h_conversion.GetBinContent(ix+1,iy)-h_conversion.GetBinContent(ix-1,iy))
					+ 0.25*(h_conversion.GetBinContent(ix+1,iy+1)-h_conversion.GetBinContent(ix-1,iy+1))
					+ 0.25*(h_conversion.GetBinContent(ix+1,iy-1)-h_conversion.GetBinContent(ix-1,iy-1)));
	  h_derivative_y.SetBinContent(ix,iy,
					0.5*(h_conversion.GetBinContent(ix,iy+1)-h_conversion.GetBinContent(ix,iy-1))
					+ 0.25*(h_conversion.GetBinContent(ix-1,iy+1)-h_conversion.GetBinContent(ix-1,iy-1))
					+ 0.25*(h_conversion.GetBinContent(ix+1,iy+1)-h_conversion.GetBinContent(ix+1,iy-1)));
	}
	else if(fDerivative_neighborhood==2){
	  h_derivative_x.SetBinContent(ix,iy,
					12*(h_conversion.GetBinContent(ix+1,iy)-h_conversion.GetBinContent(ix-1,iy))
					+ 8*(h_conversion.GetBinContent(ix+1,iy+1)-h_conversion.GetBinContent(ix-1,iy+1))
					+ 8*(h_conversion.GetBinContent(ix+1,iy-1)-h_conversion.GetBinContent(ix-1,iy-1))
					+ 2*(h_conversion.GetBinContent(ix+1,iy+2)-h_conversion.GetBinContent(ix-1,iy+2))
					+ 2*(h_conversion.GetBinContent(ix+1,iy-2)-h_conversion.GetBinContent(ix-1,iy-2))
					  + 6*(h_conversion.GetBinContent(ix+2,iy)-h_conversion.GetBinContent(ix-2,iy))
					+ 4*(h_conversion.GetBinContent(ix+2,iy+1)-h_conversion.GetBinContent(ix-2,iy+1))
					+ 4*(h_conversion.GetBinContent(ix+2,iy-1)-h_conversion.GetBinContent(ix-2,iy-1))
					+ 1*(h_conversion.GetBinContent(ix+2,iy+2)-h_conversion.GetBinContent(ix-2,iy+2))
					  + 1*(h_conversion.GetBinContent(ix+2,iy-2)-h_conversion.GetBinContent(ix-2,iy-2)));
	  h_derivative_y.SetBinContent(ix,iy,
					12*(h_conversion.GetBinContent(ix,iy+1)-h_conversion.GetBinContent(ix,iy-1))
					+ 8*(h_conversion.GetBinContent(ix-1,iy+1)-h_conversion.GetBinContent(ix-1,iy-1))
					+ 8*(h_conversion.GetBinContent(ix+1,iy+1)-h_conversion.GetBinContent(ix+1,iy-1))
					+ 2*(h_conversion.GetBinContent(ix-2,iy+1)-h_conversion.GetBinContent(ix-2,iy-1))
					+ 2*(h_conversion.GetBinContent(ix+2,iy+1)-h_conversion.GetBinContent(ix+2,iy-1))
					+ 6*(h_conversion.GetBinContent(ix,iy+2)-h_conversion.GetBinContent(ix,iy-2))
					+ 4*(h_conversion.GetBinContent(ix-1,iy+2)-h_conversion.GetBinContent(ix-1,iy-2))
					+ 4*(h_conversion.GetBinContent(ix+1,iy+2)-h_conversion.GetBinContent(ix+1,iy-2))
					+ 1*(h_conversion.GetBinContent(ix-2,iy+2)-h_conversion.GetBinContent(ix-2,iy-2))
					+ 1*(h_conversion.GetBinContent(ix+2,iy+2)-h_conversion.GetBinContent(ix+2,iy-2)));
	}
	else{
	  mf::LogError("CornerFinderAlg") << "Sobel derivative not supported for neighborhoods > 2.";
	  return;
	}

	} //end if Sobel

      else if(fDerivative_method.compare("local")==0){

	if(fDerivative_neighborhood==1){
	  h_derivative_x.SetBinContent(ix,iy,
					(h_conversion.GetBinContent(ix+1,iy)-h_conversion.GetBinContent(ix-1,iy)));
	  h_derivative_y.SetBinContent(ix,iy,
					(h_conversion.GetBinContent(ix,iy+1)-h_conversion.GetBinContent(ix,iy-1)));
	}
	else{
	  mf::LogError("CornerFinderAlg") << "Local derivative not yet supported for neighborhoods > 1.";
	  return;
	}
      } //end if local

      else{
	mf::LogError("CornerFinderAlg") << "Bad derivative algorithm! " << fDerivative_method;
	return;
      }

    }
  }


  //this is just a double Gaussian
  float func_blur[11][11];
  func_blur[0][0] = 0.000000;
  func_blur[0][1] = 0.000000;
  func_blur[0][2] = 0.000000;
  func_blur[0][3] = 0.000001;
  func_blur[0][4] = 0.000002;
  func_blur[0][5] = 0.000004;
  func_blur[0][6] = 0.000002;
  func_blur[0][7] = 0.000001;
  func_blur[0][8] = 0.000000;
  func_blur[0][9] = 0.000000;
  func_blur[0][10] = 0.000000;
  func_blur[1][0] = 0.000000;
  func_blur[1][1] = 0.000000;
  func_blur[1][2] = 0.000004;
  func_blur[1][3] = 0.000045;
  func_blur[1][4] = 0.000203;
  func_blur[1][5] = 0.000335;
  func_blur[1][6] = 0.000203;
  func_blur[1][7] = 0.000045;
  func_blur[1][8] = 0.000004;
  func_blur[1][9] = 0.000000;
  func_blur[1][10] = 0.000000;
  func_blur[2][0] = 0.000000;
  func_blur[2][1] = 0.000004;
  func_blur[2][2] = 0.000123;
  func_blur[2][3] = 0.001503;
  func_blur[2][4] = 0.006738;
  func_blur[2][5] = 0.011109;
  func_blur[2][6] = 0.006738;
  func_blur[2][7] = 0.001503;
  func_blur[2][8] = 0.000123;
  func_blur[2][9] = 0.000004;
  func_blur[2][10] = 0.000000;
  func_blur[3][0] = 0.000001;
  func_blur[3][1] = 0.000045;
  func_blur[3][2] = 0.001503;
  func_blur[3][3] = 0.018316;
  func_blur[3][4] = 0.082085;
  func_blur[3][5] = 0.135335;
  func_blur[3][6] = 0.082085;
  func_blur[3][7] = 0.018316;
  func_blur[3][8] = 0.001503;
  func_blur[3][9] = 0.000045;
  func_blur[3][10] = 0.000001;
  func_blur[4][0] = 0.000002;
  func_blur[4][1] = 0.000203;
  func_blur[4][2] = 0.006738;
  func_blur[4][3] = 0.082085;
  func_blur[4][4] = 0.367879;
  func_blur[4][5] = 0.606531;
  func_blur[4][6] = 0.367879;
  func_blur[4][7] = 0.082085;
  func_blur[4][8] = 0.006738;
  func_blur[4][9] = 0.000203;
  func_blur[4][10] = 0.000002;
  func_blur[5][0] = 0.000004;
  func_blur[5][1] = 0.000335;
  func_blur[5][2] = 0.011109;
  func_blur[5][3] = 0.135335;
  func_blur[5][4] = 0.606531;
  func_blur[5][5] = 1.000000;
  func_blur[5][6] = 0.606531;
  func_blur[5][7] = 0.135335;
  func_blur[5][8] = 0.011109;
  func_blur[5][9] = 0.000335;
  func_blur[5][10] = 0.000004;
  func_blur[6][0] = 0.000002;
  func_blur[6][1] = 0.000203;
  func_blur[6][2] = 0.006738;
  func_blur[6][3] = 0.082085;
  func_blur[6][4] = 0.367879;
  func_blur[6][5] = 0.606531;
  func_blur[6][6] = 0.367879;
  func_blur[6][7] = 0.082085;
  func_blur[6][8] = 0.006738;
  func_blur[6][9] = 0.000203;
  func_blur[6][10] = 0.000002;
  func_blur[7][0] = 0.000001;
  func_blur[7][1] = 0.000045;
  func_blur[7][2] = 0.001503;
  func_blur[7][3] = 0.018316;
  func_blur[7][4] = 0.082085;
  func_blur[7][5] = 0.135335;
  func_blur[7][6] = 0.082085;
  func_blur[7][7] = 0.018316;
  func_blur[7][8] = 0.001503;
  func_blur[7][9] = 0.000045;
  func_blur[7][10] = 0.000001;
  func_blur[8][0] = 0.000000;
  func_blur[8][1] = 0.000004;
  func_blur[8][2] = 0.000123;
  func_blur[8][3] = 0.001503;
  func_blur[8][4] = 0.006738;
  func_blur[8][5] = 0.011109;
  func_blur[8][6] = 0.006738;
  func_blur[8][7] = 0.001503;
  func_blur[8][8] = 0.000123;
  func_blur[8][9] = 0.000004;
  func_blur[8][10] = 0.000000;
  func_blur[9][0] = 0.000000;
  func_blur[9][1] = 0.000000;
  func_blur[9][2] = 0.000004;
  func_blur[9][3] = 0.000045;
  func_blur[9][4] = 0.000203;
  func_blur[9][5] = 0.000335;
  func_blur[9][6] = 0.000203;
  func_blur[9][7] = 0.000045;
  func_blur[9][8] = 0.000004;
  func_blur[9][9] = 0.000000;
  func_blur[9][10] = 0.000000;
  func_blur[10][0] = 0.000000;
  func_blur[10][1] = 0.000000;
  func_blur[10][2] = 0.000000;
  func_blur[10][3] = 0.000001;
  func_blur[10][4] = 0.000002;
  func_blur[10][5] = 0.000004;
  func_blur[10][6] = 0.000002;
  func_blur[10][7] = 0.000001;
  func_blur[10][8] = 0.000000;
  func_blur[10][9] = 0.000000;
  func_blur[10][10] = 0.000000;

  double temp_integral_x = 0;
  double temp_integral_y = 0;

  if(fDerivative_BlurNeighborhood>0){

    if(fDerivative_BlurNeighborhood>10){
      mf::LogWarning("CornerFinderAlg") << "WARNING...BlurNeighborhoods>10 not currently allowed. Shrinking to 10.";
      fDerivative_BlurNeighborhood=10;
    }

    TH2F *h_clone_derivative_x = (TH2F*)h_derivative_x.Clone("h_clone_derivative_x");
    TH2F *h_clone_derivative_y = (TH2F*)h_derivative_y.Clone("h_clone_derivative_y");

    temp_integral_x = 0;
    temp_integral_y = 0;

    for(int ix=1; ix<=h_derivative_x.GetNbinsX(); ix++){
      for(int iy=1; iy<=h_derivative_y.GetNbinsY(); iy++){

	temp_integral_x = 0;
	temp_integral_y = 0;

	for(int jx=ix-fDerivative_BlurNeighborhood; jx<=ix+fDerivative_BlurNeighborhood; jx++){
	  for(int jy=iy-fDerivative_BlurNeighborhood; jy<=iy+fDerivative_BlurNeighborhood; jy++){
	    temp_integral_x += h_clone_derivative_x->GetBinContent(jx,jy)*func_blur[(ix-jx)+5][(iy-jy)+5];
	    temp_integral_y += h_clone_derivative_y->GetBinContent(jx,jy)*func_blur[(ix-jx)+5][(iy-jy)+5];
	  }
	}
	h_derivative_x.SetBinContent(ix,iy,temp_integral_x);
	h_derivative_y.SetBinContent(ix,iy,temp_integral_y);

      }
    }

    delete h_clone_derivative_x;
    delete h_clone_derivative_y;


  } //end if blur


}


//-----------------------------------------------------------------------------
// Corner Score

void corner::CornerFinderAlg::create_cornerScore_histogram(TH2F const& h_derivative_x, TH2F const& h_derivative_y, TH2D & h_cornerScore){

  const int x_bins = h_derivative_x.GetNbinsX();
  const int y_bins = h_derivative_y.GetNbinsY();

  //the structure tensor elements
  double st_xx = 0., st_xy = 0., st_yy = 0.;

  for(int iy=1+fCornerScore_neighborhood; iy<=(y_bins-fCornerScore_neighborhood); iy++){
    for(int ix=1+fCornerScore_neighborhood; ix<=(x_bins-fCornerScore_neighborhood); ix++){

      if(ix==1+fCornerScore_neighborhood){
	st_xx=0.; st_xy=0.; st_yy=0.;

	for(int jx=ix-fCornerScore_neighborhood; jx<=ix+fCornerScore_neighborhood; jx++){
	  for(int jy=iy-fCornerScore_neighborhood; jy<=iy+fCornerScore_neighborhood; jy++){

	    st_xx += h_derivative_x.GetBinContent(jx,jy)*h_derivative_x.GetBinContent(jx,jy);
	    st_yy += h_derivative_y.GetBinContent(jx,jy)*h_derivative_y.GetBinContent(jx,jy);
	    st_xy += h_derivative_x.GetBinContent(jx,jy)*h_derivative_y.GetBinContent(jx,jy);

	  }
	}
      }

      // we do it this way to reduce computation time
      else{
	for(int jy=iy-fCornerScore_neighborhood; jy<=iy+fCornerScore_neighborhood; jy++){

	  st_xx -= h_derivative_x.GetBinContent(ix-fCornerScore_neighborhood-1,jy)*h_derivative_x.GetBinContent(ix-fCornerScore_neighborhood-1,jy);
	  st_xx += h_derivative_x.GetBinContent(ix+fCornerScore_neighborhood,jy)*h_derivative_x.GetBinContent(ix+fCornerScore_neighborhood,jy);

	  st_yy -= h_derivative_y.GetBinContent(ix-fCornerScore_neighborhood-1,jy)*h_derivative_y.GetBinContent(ix-fCornerScore_neighborhood-1,jy);
	  st_yy += h_derivative_y.GetBinContent(ix+fCornerScore_neighborhood,jy)*h_derivative_y.GetBinContent(ix+fCornerScore_neighborhood,jy);

	  st_xy -= h_derivative_x.GetBinContent(ix-fCornerScore_neighborhood-1,jy)*h_derivative_y.GetBinContent(ix-fCornerScore_neighborhood-1,jy);
	  st_xy += h_derivative_x.GetBinContent(ix+fCornerScore_neighborhood,jy)*h_derivative_y.GetBinContent(ix+fCornerScore_neighborhood,jy);
	}
      }

      if( fCornerScore_algorithm.compare("Noble")==0 ) {
	h_cornerScore.SetBinContent(ix,iy,
				     (st_xx*st_yy-st_xy*st_xy) / (st_xx+st_yy + fCornerScore_Noble_epsilon));
      }
	else if( fCornerScore_algorithm.compare("Harris")==0 ) {
	  h_cornerScore.SetBinContent(ix,iy,
				       (st_xx*st_yy-st_xy*st_xy) - ((st_xx+st_yy)*(st_xx+st_yy)*fCornerScore_Harris_kappa));
	}
	else{
	  mf::LogError("CornerFinderAlg") << "BAD CORNER ALGORITHM: " << fCornerScore_algorithm;
	  return;
	}

    } // end for loop over x bins
  } // end for loop over y bins

}


//-----------------------------------------------------------------------------
// Max Supress
size_t corner::CornerFinderAlg::perform_maximum_suppression(TH2D const& h_cornerScore,
							    std::vector<recob::EndPoint2D> & corner_vector,
							    std::vector<geo::WireID> wireIDs,
							    geo::View_t view,
							    TH2D & h_maxSuppress,
							    int startx,
                                                            int starty) const {

  const int x_bins = h_cornerScore.GetNbinsX();
  const int y_bins = h_cornerScore.GetNbinsY();

  double temp_max;
  bool temp_center_bin;

  for(int iy=1; iy<=y_bins; iy++){
    for(int ix=1; ix<=x_bins; ix++){

      if(h_cornerScore.GetBinContent(ix,iy) < fMaxSuppress_threshold)
	continue;

      temp_max = -1000;
      temp_center_bin = false;

      for(int jx=ix-fMaxSuppress_neighborhood; jx<=ix+fMaxSuppress_neighborhood; jx++){
	for(int jy=iy-fMaxSuppress_neighborhood; jy<=iy+fMaxSuppress_neighborhood; jy++){

	  if(h_cornerScore.GetBinContent(jx,jy) > temp_max){
	    temp_max = h_cornerScore.GetBinContent(jx,jy);
	    if(jx==ix && jy==iy) temp_center_bin=true;
	    else{ temp_center_bin=false; }
	  }

	}
      }

      if(temp_center_bin){

	float time_tick = 0.5 * (float)((2*(iy+starty)) * fConversion_bins_per_input_y);
	int wire_number = ( (2*(ix+startx))*fConversion_bins_per_input_x ) / 2;
	double totalQ = 0;
	int id = 0;
	recob::EndPoint2D corner(time_tick,
				 wireIDs[wire_number],
				 h_cornerScore.GetBinContent(ix,iy),
				 id,
				 view,
				 totalQ);
	corner_vector.push_back(corner);

	h_maxSuppress.SetBinContent(ix,iy,h_cornerScore.GetBinContent(ix,iy));
      }

    }
  }

  return corner_vector.size();

}


/* Silly little function for doing a line integral type thing. Needs improvement. */
float corner::CornerFinderAlg::line_integral(TH2F const& hist, int begin_x, float begin_y, int end_x, float end_y, float threshold) const{

  int x1 = hist.GetXaxis()->FindBin( begin_x );
  int y1 = hist.GetYaxis()->FindBin( begin_y );
  int x2 = hist.GetXaxis()->FindBin( end_x );
  int y2 = hist.GetYaxis()->FindBin( end_y );

  if(x1==x2 && abs(y1-y2)<1e-5)
    return 0;

  if(x2<x1){
    int tmp = x2;
    x2 = x1;
    x1 = tmp;

    int tmp_y = y2;
    y2 = y1;
    y1 = tmp_y;
  }

  float fraction = 0;
  int bin_counter = 0;

  if(x2!=x1){

    float slope = (y2-y1)/((float)(x2-x1));

    for(int ix=x1; ix<=x2; ix++){

      int y_min,y_max;

      if(slope>=0){
	y_min = y1 + slope*(ix-x1);
	y_max = y1 + slope*(ix+1-x1);
      }
      else {
	y_max = (y1+1) + slope*(ix-x1);
	y_min = (y1+1) + slope*(ix+1-x1);
      }

      for(int iy=y_min; iy<=y_max; iy++){
	bin_counter++;

	if( hist.GetBinContent(ix,iy) > threshold )
	  fraction += 1.;
      }

    }
  }
  else{

    int y_min,y_max;
    if(y1<y2){
      y_min=y1; y_max=y2;
    }
    else{
      y_min=y2; y_max=y1;
    }
    for(int iy=y_min; iy<=y_max; iy++){
	bin_counter++;
	if( hist.GetBinContent(x1,iy) > threshold)
	  fraction += 1.;
      }

  }

  return fraction/bin_counter;
}


//-----------------------------------------------------------------------------
// Do the silly little line integral score thing
size_t corner::CornerFinderAlg::calculate_line_integral_score( TH2F const& h_wire_data,
								std::vector<recob::EndPoint2D> const & corner_vector,
								std::vector<recob::EndPoint2D> & corner_lineIntegralScore_vector,
                                                                TH2F & h_lineIntegralScore) const {

  float score;

  for(auto const i_corner : corner_vector){

    score=0;

    for(auto const j_corner : corner_vector){


      if( line_integral(h_wire_data,
			i_corner.WireID().Wire,i_corner.DriftTime(),
			j_corner.WireID().Wire,j_corner.DriftTime(),
			fIntegral_bin_threshold) > fIntegral_fraction_threshold)
	{
	  score+=1.;
	}

    }

    recob::EndPoint2D corner(i_corner.DriftTime(),
			     i_corner.WireID(),
			     score,
			     i_corner.ID(),
			     i_corner.View(),
			     i_corner.Charge());

    corner_lineIntegralScore_vector.push_back(corner);


    h_lineIntegralScore.SetBinContent(h_wire_data.GetXaxis()->FindBin(i_corner.WireID().Wire),
				      h_wire_data.GetYaxis()->FindBin(i_corner.DriftTime()),
				      score);

  }

  return corner_lineIntegralScore_vector.size();
}



TH2F const& corner::CornerFinderAlg::GetWireDataHist(unsigned int i_plane) const {
  return WireData_histos.at(i_plane);
}
/*
TH2F* corner::CornerFinderAlg::GetConversionHist(unsigned int i_plane){

  if(i_plane >= fConversion_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fConversion_histos.at(i_plane);
}

TH2F* corner::CornerFinderAlg::GetDerivativeXHist(unsigned int i_plane){

  if(i_plane >= fDerivativeX_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fDerivativeX_histos.at(i_plane);
}

TH2F* corner::CornerFinderAlg::GetDerivativeYHist(unsigned int i_plane){

  if(i_plane >= fDerivativeY_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fDerivativeY_histos.at(i_plane);
}

TH2D* corner::CornerFinderAlg::GetCornerScoreHist(unsigned int i_plane){

  if(i_plane >= fCornerScore_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fCornerScore_histos.at(i_plane);
}

TH2D* corner::CornerFinderAlg::GetMaxSuppressHist(unsigned int i_plane){

  if(i_plane >= fMaxSuppress_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fMaxSuppress_histos.at(i_plane);
}
*/
