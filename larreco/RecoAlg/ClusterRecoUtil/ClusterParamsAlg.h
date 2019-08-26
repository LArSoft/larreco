////////////////////////////////////////////////////////////////////////
// ClusterParamsAlg.h
//
// ClusterParamsAlg class
//
// Andrzej Szelc (andrzej.szelc@yale.edu)
//
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTERPARAMSALG_H
#define CLUSTERPARAMSALG_H

//--- std/stl include ---//
#include <vector>
#include <string>

//--- LArSoft include ---//
#include "ClusterParams.h"
#include "lardata/Utilities/PxUtils.h"
namespace util { class GeometryUtilities; }

// ... more #include in tempimplementation below

namespace cluster {

  class ClusterParamsAlg {

  public:

    /// Default constructor
    ClusterParamsAlg();

    /// Alternative constructor with larutil::PxHit vector
    ClusterParamsAlg(const std::vector<util::PxHit>&);

    void Initialize();

    //void SetHits(const std::vector<larutil::PxHit*>&);

    void SetMinNHits(size_t nhit) { fMinNHits = nhit; }

    size_t MinNHits() const { return fMinNHits; }

   // int SetHits(const std::vector<const ::larlite::hit*> &);

    int SetHits(const std::vector<util::PxHit> &);

    void SetRefineDirectionQMin(double qmin){ fQMinRefDir = qmin; }

    void SetVerbose(bool yes=true){ verbose = yes;}

   // void SetArgoneutGeometry();

    template <typename Stream>
    void Report(Stream& stream) const;

    template <typename Stream>
    void TimeReport(Stream& stream) const;

    /**
     * This function returns a feature vector suitable for a neural net
     * This function uses the data from cluster_params but packages it
     * up in a different way, and so is inappropriate to include in
     * clusterParams.hh.   That's why it's here.
     * @param  data  takes a reference to a vector< float>
     */
    void  GetFANNVector(std::vector<float> & data);
    // std::vector<float> & GetFANNVector();

    /**
     * For debugging purposes, prints the result of GetFANNVector
     * in a nicely formatted form.
     * @return [description]
     */
    void  PrintFANNVector();


    /**
      Runs all the functions which calculate cluster params
      and stashes the results in the private ClusterParams
      struct.

      @param override_DoGetAverages       force re-execution of GetAverages()
      @param override_DoGetRoughAxis      force re-execution of GetRoughAxis()
      @param override_DoGetProfileInfo    force re-execution of GetProfileInfo()
      @param override_DoRefineStartPoints force re-execution of RefineStartPoints()
      @param override_DoGetFinalSlope     force re-execution of GetFinalSlope()
      @param override_DoEndCharge         force re-execution of GetEndCharges()
    */
    void FillParams(bool override_DoGetAverages      =false,
                    bool override_DoGetRoughAxis     =false,
                    bool override_DoGetProfileInfo   =false,
                    bool override_DoRefineStartPointsAndDirection=false,
            		    // bool override_DoRefineDirection  =false,
                    bool override_DoGetFinalSlope    =false,
                    bool override_DoTrackShowerSep   =false,
                    bool override_DoEndCharge = false);

    const cluster_params& GetParams() const
    { return fParams;}

    /**
       Calculates the following variables:
       mean_charge
       mean_x
       mean_y
       charge_wgt_x
       charge_wgt_y
       eigenvalue_principal
       eigenvalue_secondary
       multi_hit_wires
       N_Wires
       @param override force recalculation of variables
    */
    void GetAverages(bool override=false);


    /**
      Calculates the following variables:
      verticalness
      fRough2DSlope
      fRough2DIntercept
      @param override [description]
    */
    //void GetRoughAxis(bool override=false);
    void GetRoughAxis(bool override=false);


    /**
       Calculates the following variables:
       opening_angle
       opening_angle_highcharge
       closing_angle
       closing_angle_highcharge
       offaxis_hits
       @param override [description]
    */
    void GetProfileInfo(bool override=false);


    /**
       Calculates the following variables:
       length
       width
       @param override [description]
    */
    void RefineStartPoints(bool override=false);

    /**
       Calculates the following variables:
       hit_density_1D
       hit_density_2D
       angle_2d
       direction
       @param override [description]
    */
    void GetFinalSlope(bool override=false);

    /**
       Calculates the following variables:
       start_charge
       end_charge
       @param override_ force recompute the variables
       @see StartCharge(), EndCharge()
    */
    void GetEndCharges(bool override_ = false);

    void RefineDirection(bool override=false);

    void RefineStartPointAndDirection(bool override=false);

    void TrackShowerSeparation(bool override=false);

    void setNeuralNetPath(std::string s){fNeuralNetPath = s;}

    void FillPolygon();

    void GetOpeningAngle();

    const util::PxPoint& RoughStartPoint() {return fRoughBeginPoint;}
    const util::PxPoint& RoughEndPoint() {return fRoughEndPoint;}

    double RoughSlope() {return fRough2DSlope;}
    double RoughIntercept() {return fRough2DIntercept;}

    /**
     * @brief Returns the expected charge at the beginning of the cluster
     * @param nbins use at least this number of charge bins from charge profile
     * @param length space at the start of cluster where to collect charge, in cm
     * @brief the expected charge at the beginning of the cluster
     * @see EndCharge(), IntegrateFitCharge()
     *
     * ClusterParamsAlg extracts a binned charge profile, parametrized versus
     * the distance from the start of the cluster.
     * All the charge on the plane orthogonal to cluster axis is collapsed into
     * the point where that plane intersects the axis.
     * The resulting 1D distribution is then binned.
     *
     * This method returns the charge under the first length cm of the cluster.
     *
     * This method considers the first nbins of this charge distribution and
     * through a linear fit determines the expected charge at the first bin.
     * Then, it scales the result to reflect how much charge would be deposited
     * in a space of length centimetres, according to this linear fit.
     *
     * Note that length may be 0 (charge will be 0) or negative (sort of
     * extrapolation ahead of the cluster start).
     *
     * For more details, see IntegrateFitCharge().
     */
    double StartCharge(float length = 1., unsigned int nbins = 10);

    /**
     * @brief Returns the expected charge at the end of the cluster
     * @param nbins use at least this number of charge bins from charge profile
     * @param length space before the end of cluster where to collect charge, in cm
     * @brief the expected charge at the end of the cluster
     * @see StartCharge(), IntegrateFitCharge()
     *
     * This method returns the charge under the last length cm of the cluster.
     * See StartCharge() for a detailed explanation.
     * For even more details, see IntegrateFitCharge().
     */
    double EndCharge(float length = 1., unsigned int nbins = 10);

    /**
     * @brief Returns the number of multiple hits per wire
     * @return the number of multiple hits per wire
     *
     * This returns the fraction of wires that have more than one hit belonging
     * to this cluster.
     */
    float MultipleHitWires();

    /**
     * @brief Returns the number of multiple hits per wire
     * @return the number of multiple hits per wire
     *
     * This returns the number of wires with mmore than one hit belonging
     * to this cluster, divided by the cluster length in cm.
     */
    float MultipleHitDensity();


    void EnableFANN();

    void DisableFANN(){enableFANN = false;}

    size_t GetNHits() const {return fHitVector.size();}
    const std::vector<util::PxHit>& GetHitVector() const {return fHitVector;}
    int Plane() const {return fPlane;}
    void SetPlane(int p);

  protected:

    util::GeometryUtilities  *fGSer;

    /// Cut value for # hits: below this value clusters are not evaluated
    size_t fMinNHits;

    /**
       This vector holds the pointer to hits.
       This should be used for computation for speed.
    */
    std::vector<util::PxHit> fHitVector;

    // bool to control debug/verbose mode
    // defaults to off.
    bool verbose;

    //settable parameters:
    std::vector<double> fChargeCutoffThreshold;
    int fPlane;

    //this is required in RefineDirection
    double fQMinRefDir;

    std::vector< double > fChargeProfile;
    std::vector< double > fCoarseChargeProfile;

    std::vector< double > fChargeProfileNew;
   // double fMaxLinLength;
   // double fLinBins;

    int fCoarseNbins;
    int fProfileNbins;
    int fProfileMaximumBin;
    double fProfileIntegralForward;
    double fProfileIntegralBackward;
    double fProjectedLength;

    //extreme intercepts using the rough_2d_slope
   // double fInterHigh;
   // double fInterLow;
    double fBeginIntercept;
    double fEndIntercept;
    double fInterHigh_side;
    double fInterLow_side;

    // book keeping variables to validate completion of methods:
    bool fFinishedGetAverages;
    bool fFinishedGetRoughAxis;
    bool fFinishedGetProfileInfo;
    bool fFinishedRefineStartPoints;
    bool fFinishedRefineDirection;
    bool fFinishedGetFinalSlope;
    bool fFinishedRefineStartPointAndDirection;
    bool fFinishedTrackShowerSep;
    bool fFinishedGetEndCharges;

    double fRough2DSlope;        // slope
    double fRough2DIntercept;    // slope
    util::PxPoint fRoughBeginPoint;
    util::PxPoint fRoughEndPoint;
    bool enableFANN;

    /**
     * @brief Integrates the charge between two positions in the cluster axis
     * @param from_length position on the axis to start integration from, in cm
     * @param to_length position on the axis to end integration at, in cm
     * @param fit_first_bin first bin for the charge fit
     * @param fit_end_bin next-to-last bin for the charge fit
     * @return the charged fit integrated in the specified range, in ADC counts
     * @see StartCharge(), EndCharge()
     *
     * This function provides an almost-punctual charge at a position in the
     * axis. Since the effective punctual charge is 0 ADC counts by definition,
     * the charge can be integrated for some length.
     * The procedure is made of two steps:
     * 1. the charge profile is parametrized with a linear fit within the
     *    specified region
     * 2. an integration of that fit is performed along the segment specified.
     *
     * The region at point 1. is from `fit_first_bin` to `fit_end_bin`. These
     * are specified in bin units. The binning is the one of the charge profile.
     * It is suggested that a few bins are always kept, say 5 to 10, to reduce
     * statistical fluctuations but maintaining a decent hypothesis of linearity
     * along the range.
     * The linear fit weighs all the bins in the profile the same.
     *
     * The region at point to is from `from_length` to `to_length`, and it is
     * measured in cm along the cluster axis, starting at the start of the
     * cluster.
     */
    double IntegrateFitCharge(
      double from_length, double to_length,
      unsigned int fit_first_bin, unsigned int fit_end_bin
      );

    /// Returns the integral of f(x) = mx + q defined in [x1, x2]
    static double LinearIntegral(double m, double q, double x1, double x2);

    public:

    cluster::cluster_params fParams;

    std::string fNeuralNetPath;

    std::vector<std::string> fTimeRecord_ProcName;
    std::vector<double> fTimeRecord_ProcTime;

  }; //class ClusterParamsAlg

} //namespace cluster


//------------------------------------------------------------------------------
//--- template implementation
//---

#include <ostream> // std::endl

namespace cluster {

  template <typename Stream>
  void ClusterParamsAlg::TimeReport(Stream& stream) const {

    stream << "  <<ClusterParamsAlg::TimeReport>> starts..."<<std::endl;
    for(size_t i=0; i<fTimeRecord_ProcName.size(); ++i){

      stream << "    Function: "
        << fTimeRecord_ProcName[i].c_str()
        << " ... Time = "
        << fTimeRecord_ProcTime[i]
        << " [s]"
        << std::endl;

    }
    stream<< "  <<ClusterParamsAlg::TimeReport>> ends..."<<std::endl;
  }

} //namespace cluster

#endif
