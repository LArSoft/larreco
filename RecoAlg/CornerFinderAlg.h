////////////////////////////////////////////////////////////////////////
/// \file  CornerFinderAlgAlg.h
/// \brief algorithm to find feature 2D points
///
/// \author  
////////////////////////////////////////////////////////////////////////

#ifndef CORNERFINDERALG_H
#define CORNERFINDERALG_H

#include "fhiclcpp/ParameterSet.h" 

#include "TMath.h"
#include "TH2.h"
#include "TF2.h"
#include "TH1D.h"
#include <vector>
#include <string>
#include "RecoBase/Wire.h"
#include "RecoBase/EndPoint2D.h"
#include "Geometry/Geometry.h"


namespace trkf {
  class BezierTrack;
}

namespace corner { //<---Not sure if this is the right namespace

   class CornerFinderAlg {
   
   public:
    
     explicit CornerFinderAlg(fhicl::ParameterSet const& pset); 
     virtual ~CornerFinderAlg();        
     
     void   reconfigure(fhicl::ParameterSet const& pset);
    
     
     void GrabWires( std::vector<recob::Wire> const& wireVec ,
		     geo::Geometry const&);                                      //this one creates the histograms we want to use
     
     
     void get_feature_points(std::vector<recob::EndPoint2D> &, 
			     geo::Geometry const&);                              //here we get feature points with corner score

     void get_feature_points_LineIntegralScore(std::vector<recob::EndPoint2D> &,
					       geo::Geometry const&);            //here we get feature points with LineIntegral score

     void get_feature_points_fast(std::vector<recob::EndPoint2D> &,
				  geo::Geometry const&);                         //here we get feature points with corner score
     
     float line_integral(TH2F const& hist, int x1, float y1, int x2, float y2, float threshold);				   
     
     std::vector<float> line_integrals(trkf::BezierTrack&, size_t Steps, float threshold);
     
     TH2F const& GetWireDataHist(unsigned int);
     TH2F const& GetConversionHist(unsigned int);
     TH2F const& GetDerivativeXHist(unsigned int);
     TH2F const& GetDerivativeYHist(unsigned int);
     TH2D const& GetCornerScoreHist(unsigned int);
     TH2D const& GetMaxSuppressHist(unsigned int);
     
    private:
     
     void CleanCornerFinderAlg();
     void InitializeGeometry(geo::Geometry const&);
     
     // Need to list the things we will take in from the .fcl file
     
     std::string  fCalDataModuleLabel;
     std::string    fConversion_algorithm;
     std::string    fConversion_func;
     float          fTrimming_threshold;
     int            fTrimming_buffer;
     double         fTrimming_totalThreshold;
     int            fConversion_func_neighborhood;
     float          fConversion_threshold;
     int            fConversion_bins_per_input_x;
     int            fConversion_bins_per_input_y;
     std::string    fDerivative_method;
     int            fDerivative_neighborhood;
     std::string    fDerivative_BlurFunc;
     int            fDerivative_BlurNeighborhood;
     int            fCornerScore_neighborhood;
     std::string    fCornerScore_algorithm;
     float          fCornerScore_Noble_epsilon;
     float          fCornerScore_Harris_kappa;
     int            fMaxSuppress_neighborhood;
     int            fMaxSuppress_threshold;
     float          fIntegral_bin_threshold;
     float          fIntegral_fraction_threshold;
     
     // Making a vector of histograms
     std::vector<TH2F> WireData_histos;
     std::vector<TH1D> WireData_histos_ProjectionX;
     std::vector<TH1D> WireData_histos_ProjectionY;
     std::vector< std::tuple<int,TH2F,int,int> > WireData_trimmed_histos;
     std::vector< std::vector<geo::WireID> > WireData_IDs;
     std::vector<TH2F> fConversion_histos;
     std::vector<TH2F> fDerivativeX_histos;
     std::vector<TH2F> fDerivativeY_histos;
     std::vector<TH2D> fCornerScore_histos;
     std::vector<TH2D> fMaxSuppress_histos;
     
     unsigned int event_number;
     unsigned int run_number;
     
     void create_image_histo(TH2F const& h_wire_data, TH2F & h_conversion);
     void create_derivative_histograms(TH2F const& h_conversion, TH2F & h_derivative_x, TH2F & h_derivative_y);
     void create_cornerScore_histogram(TH2F const& h_derivative_x, TH2F const& h_derivative_y, TH2D & h_cornerScore);
     size_t perform_maximum_suppression(TH2D const& h_cornerScore, 
					std::vector<recob::EndPoint2D> & corner_vector,
					std::vector<geo::WireID> wireIDs, 
					geo::View_t view,
					TH2D & h_maxSuppress,
					int startx=0,
					int starty=0);
     
     size_t calculate_line_integral_score( TH2F const& h_wire_data, 
					   std::vector<recob::EndPoint2D> const & corner_vector, 
					   std::vector<recob::EndPoint2D> & corner_lineIntegralScore_vector,
					   TH2F & h_lineIntegralScore);
     
     void attach_feature_points(TH2F const& h_wire_data, 
				std::vector<geo::WireID> wireIDs, 
				geo::View_t view,
				std::vector<recob::EndPoint2D>&,
				int startx=0,int starty=0);
     void attach_feature_points_LineIntegralScore(TH2F const& h_wire_data, 
						  std::vector<geo::WireID> wireIDs, 
						  geo::View_t view,
						  std::vector<recob::EndPoint2D>&);
     
     
     void create_smaller_histos(geo::Geometry const&);
     void remove_duplicates(std::vector<recob::EndPoint2D>&);
     
   };//<---End of class CornerFinderAlg
   
   
}

#endif //CORNERFINDERALG_H
