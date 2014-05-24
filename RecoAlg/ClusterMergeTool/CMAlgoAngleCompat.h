/**
 * \file CMAlgoAngleCompat.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoAngleCompat
 *
 * @author davidkaleko
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOANGLECOMPAT_H
#define CMALGOANGLECOMPAT_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoAngleCompat
     User defined class CMAlgoAngleCompat ... these comments are used to generate
     doxygen documentation!
  */
  class CMAlgoAngleCompat : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoAngleCompat();
    
    /// Default destructor
    virtual ~CMAlgoAngleCompat(){};
        
    /// Overloaded (from CBoolAlgoBase) Bool function
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Overloaded function called at beginning of iterating over all pairs of clusters
    /// Currently using just to initialize some debugging histograms
    /// should rewrite it to not take an input argument...
    virtual void IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /// Method to set verbose mode
    void SetVerbose(bool on) { _verbose = on; }

    /// Method to set whether you allow angles to match with +/- 180 deg difference
    void SetAllow180Ambig(bool on) { _allow_180_ambig = on; }

    /// Method to set cut value in degrees for angle compatibility test
    void SetAngleCut(double angle) { _max_allowed_2D_angle_diff = angle; }

    /// Method to set angle cut value to be based on opening angle
    void SetUseOpeningAngle(bool on) { _use_opening_angle = on; }

    /// Set Minimum Number of Hits to consider Cluster
    void SetMinHits(int n) { _minHits = n; }

    TH1F* GetAngleDistHisto() const{ return angle_dist_histo; };

  protected:

    bool _verbose;    /// bool to suppress lots of output if you want

    ///bool to allow "backwards" clusters (swapped start/end points)
    ///to still match in angle, even though they are 180 degrees apart
    ///only valid for _use_opening_angle==false
    bool _allow_180_ambig; 

    /// hard shower-axis angle cutoff (only valid for _use_opening_angle==false)
    double _max_allowed_2D_angle_diff; //in degrees

    /// bool set to true if you want to use opening angle as the cutoff
    /// angle instead of whatever you set with SetAngleCut
    bool _use_opening_angle;

    int _minHits;        /// Min Number of hits for cluster to be considered

    /// Histogram used for debugging/cut value settings
    TH1F *angle_dist_histo;


  };
  
} // end namespace cluster

#endif
  /** @} */ // end of doxygen group 
  
