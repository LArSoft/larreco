/**
 * \file ShowerRecoAlgBase.h
 *
 * \ingroup ShowerReco3D
 * 
 * \brief Class def header for a class ShowerRecoAlgBase
 *
 * @author kazuhiro
 */

/** \addtogroup ShowerReco3D

    @{*/
#ifndef RECOTOOL_SHOWERRECOALGBASE_H
#define RECOTOOL_SHOWERRECOALGBASE_H

#include <iostream>
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include "ShowerRecoException.h"
#include <limits>
#include <climits>
namespace showerreco {

  struct ShowerCluster_t {

    ::util::PxPoint start_point;
    ::util::PxPoint end_point;
    double             angle_2d;
    unsigned short     plane_id;
    std::vector< ::util::PxHit> hit_vector;

    ShowerCluster_t() : hit_vector()
    {}

  };
  
  /**
     \class ShowerRecoAlgBase
     User defined class ShowerRecoAlgBase ... these comments are used to generate
     doxygen documentation!
  */
  class ShowerRecoAlgBase{
    
  public:
    
    /// Default constructor
    ShowerRecoAlgBase();
    
    /// Default destructor
    virtual ~ShowerRecoAlgBase(){}

    /// Function to reset algorithm, to be called @ beginning of each event
    virtual void Reset() = 0;

    /// Setter for a matched combination of clusters
    virtual void AppendInputClusters(const std::vector<cluster::ClusterParamsAlg>& cpan_v);

    /// Execute reconstruction
    std::vector<recob::Shower> Reconstruct();

    /// Verbosity switch
    virtual void Verbose(bool on=true) { fVerbosity=on; }

    /// Calorimetry algorithm setter
    void CaloAlgo(::calo::CalorimetryAlg* alg) { fCaloAlg = alg; }

  protected:

    /// Function to reorganize input cluster information
    virtual void ProcessInputClusters()
    { return; }

    /// Function to reconstruct one shower
    virtual ::recob::Shower RecoOneShower(const std::vector<showerreco::ShowerCluster_t>& clusters) = 0;
    
  protected:
    
    /// Verbosity flag
    bool fVerbosity;

    /// Calorimetry algorithm
    ::calo::CalorimetryAlg *fCaloAlg;

    /// Input clusters
    std::vector<std::vector<showerreco::ShowerCluster_t> > fInputClusters;
  };
}

#endif
/** @} */ // end of doxygen group 

