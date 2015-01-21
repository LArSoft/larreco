////////////////////////////////////////////////////////////////////////
// \file ClusterMergeHelper.h
//
// \brief ClusterMergeHelper header file
//
// \author kazuhiro@nevis.columbia.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef CLUSTERMERGEHELPER_H
#define CLUSTERMERGEHELPER_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include "RecoAlg/CMTool/CMToolBase/CMergeManager.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/Geometry.h"

// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace cluster
{

  class ClusterMergeHelper {
    
  public:

    /// Default constructor with fhicl parameters
    ClusterMergeHelper(){};

    /// Default destructor
    virtual ~ClusterMergeHelper(){}

    /// A method to retrieve Manager
    ::cmtool::CMergeManager& GetManager() { return fMgr; }

    /// Utility method to set cluster input information to CMergeManager from LArSoft data product (vector of recob::Hit art::Ptr)
    void SetClusters(const std::vector<std::vector<art::Ptr<recob::Hit> > > &clusters);

    /// Utility method to set cluster input information to CMerteManager from art::Event and cluster data product label
    void SetClusters(const art::Event &evt, const std::string &cluster_module_label);

    /// Function to execute CMergeManager::Process()
    void Process();

    /// Utility method to retrieve merged clusters in terms of a vector of art::Ptr<recob::Hit>
    const std::vector<std::vector<art::Ptr<recob::Hit> > >& GetMergedClusterHits() const;

    /// Utility method to retrieve merged clusters in terms of a vector of CPAN
    const std::vector<cluster::ClusterParamsAlg>& GetMergedCPAN() const;

    /// Utility method to append result set to user's data product storage
    void AppendResult(art::EDProducer &ed,
		      art::Event      &ev,
		      std::vector<recob::Cluster>           &out_clusters,
		      art::Assns<recob::Cluster,recob::Hit> &assns) const;
      
  protected:

    /// Internal method to transfer input cluster information in the right format to CMergeManager
    void SetClusters(const std::vector<std::vector<util::PxHit> > &clusters)
    { 
      fMgr.Reset();
      fMgr.SetClusters(clusters); 
    }

  protected:

    /// CMergeManager instance
    ::cmtool::CMergeManager fMgr;

    /// GeometryUtilities 
    ::util::GeometryUtilities fGeoU;

    /// Input clusters in terms of a vector of art::Ptr<recob::Hit> collection
    std::vector<std::vector<art::Ptr<recob::Hit> > > fInputClusters;

    /// Output clusters in terms of a vector of art::Ptr<recob::Hit> collection
    std::vector<std::vector<art::Ptr<recob::Hit> > > fOutputClusters;

  }; // class ClusterMergeHelper
  
} //namespace cluster
#endif
