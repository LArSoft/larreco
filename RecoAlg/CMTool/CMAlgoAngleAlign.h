/**
 * \file CMAlgoAngleAlign.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMAlgoAngleCompat
 *
 * @author davidkaleko
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMALGOANGLEALIGN_HH
#define CMALGOANGLEALIGN_HH

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoAngleCompat
     User defined class CMAlgoAngleCompat ... these comments are used to generate
     doxygen documentation!
  */
  class CMAlgoAngleAlign : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoAngleAlign();
    
    /// Default destructor
    virtual ~CMAlgoAngleAlign(){};
        
    /// Overloaded (from CBoolAlgoBase) Bool function
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Method to set debug mode
    void SetDebug(bool on) { _debug = on; }

    /// Method to set whether you allow angles to match with +/- 180 deg difference
    void SetAllow180Ambig(bool on) { _allow_180_ambig = on; }

    /// Method to set cut value in degrees for angle compatibility test
    void SetAngleCut(double angle) { _MaxAngleSep = angle; }

    void SetMinNHits(size_t n) { _MinNHits = n; }

  protected:

    bool _debug;
    size_t _MinNHits;    /// minimum number of hits for cluster to be considered
    
    ///bool to allow "backwards" clusters (swapped start/end points)
    ///to still match in angle, even though they are 180 degrees apart
    ///only valid for _use_opening_angle==false
    bool _allow_180_ambig; 

    /// hard shower-axis angle cutoff (only valid for _use_opening_angle==false)
    double _MaxAngleSep;


  };
  
} // end namespace cluster

#endif
  /** @} */ // end of doxygen group 
  
