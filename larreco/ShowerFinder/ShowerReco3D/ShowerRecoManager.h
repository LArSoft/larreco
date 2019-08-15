/**
 * \file ShowerRecoManager.h
 *
 * \ingroup ShowerReco3D
 *
 * \brief Class def header for a class ShowerRecoManager
 *
 * @author kazuhiro
 */

/** \addtogroup ShowerReco3D

    @{*/
#ifndef LARLITE_SHOWERRECOMANAGER_H
#define LARLITE_SHOWERRECOMANAGER_H

#include <vector>

#include "lardata/Utilities/PxUtils.h"
namespace cmtool { class CMatchManager; }
namespace recob { class Shower; }

namespace showerreco {

  class ShowerRecoAlgBase;

  typedef std::vector<std::vector<unsigned int> > ClusterAss_t;
  typedef std::vector< ::util::PxHit> PxHitSet_t;

  /**
     \class ShowerRecoManager
     User defined class ShowerRecoManager ... these comments are used to generate
     doxygen documentation!
  */
  class ShowerRecoManager{

  public:

    /// Default constructor
    ShowerRecoManager();

    /// Default destructor
    ~ShowerRecoManager(){}

    void Algo (ShowerRecoAlgBase* alg) { fShowerAlgo = alg; }

    const ShowerRecoAlgBase* Algo() const { return fShowerAlgo; }

    void Reset();

    ClusterAss_t Reconstruct (const std::vector<std::vector<util::PxHit> >& clusters,
			      std::vector< ::recob::Shower>& showers);

    void Reconstruct (const std::vector<std::vector<util::PxHit> >& clusters,
		      const ClusterAss_t& ass,
		      std::vector< ::recob::Shower>& showers);

    ::cmtool::CMatchManager& MatchManager() { return *fMatchMgr; }

  private:

    /// Boolean flag to whether or not to run matching
    bool fMatch;

    /// Shower reconstruction algorithm
    ::showerreco::ShowerRecoAlgBase *fShowerAlgo;

    /// Cluster matching code
    ::cmtool::CMatchManager *fMatchMgr;

    void Process(const ClusterAss_t& ass,
		 std::vector< ::recob::Shower>& showers);

  };
}

#endif
/** @} */ // end of doxygen group
