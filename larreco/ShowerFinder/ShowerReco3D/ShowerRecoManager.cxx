#include "ShowerRecoManager.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/CMTool/CMToolBase/CMTException.h"
#include "larreco/RecoAlg/CMTool/CMToolBase/CMatchBookKeeper.h"
#include "larreco/RecoAlg/CMTool/CMToolBase/CMatchManager.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include "larreco/ShowerFinder/ShowerReco3D/ShowerRecoAlgBase.h"

namespace showerreco {

  ShowerRecoManager::ShowerRecoManager() : fShowerAlgo(nullptr), fMatchMgr(nullptr)
  {
    fMatch = true;
    art::ServiceHandle<geo::Geometry const> geom;
    fMatchMgr = new ::cmtool::CMatchManager(geom->Nplanes());
  }

  void ShowerRecoManager::Reset()
  {
    if (fShowerAlgo) fShowerAlgo->Reset();
    fMatchMgr->Reset();
  }

  ClusterAss_t ShowerRecoManager::Reconstruct(geo::GeometryCore const& geom,
                                              detinfo::DetectorClocksData const& clockData,
                                              detinfo::DetectorPropertiesData const& detProp,
                                              const std::vector<std::vector<util::PxHit>>& clusters,
                                              std::vector<::recob::Shower>& showers)
  {
    util::GeometryUtilities const gser{geom, clockData, detProp};
    showers.clear();
    fMatchMgr->SetClusters(gser, clusters);

    ClusterAss_t res_ass;
    // Run matching & retrieve matched cluster indices
    try {
      fMatchMgr->Process(gser);
    }
    catch (::cmtool::CMTException& e) {
      e.what();
      return res_ass;
    }
    res_ass = fMatchMgr->GetBookKeeper().GetResult();

    Process(geom, clockData, detProp, res_ass, showers);

    return res_ass;
  }

  void ShowerRecoManager::Reconstruct(geo::GeometryCore const& geom,
                                      detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      const std::vector<std::vector<util::PxHit>>& clusters,
                                      const ClusterAss_t& ass,
                                      std::vector<::recob::Shower>& showers)
  {
    showers.clear();
    util::GeometryUtilities const gser{geom, clockData, detProp};
    fMatchMgr->SetClusters(gser, clusters);

    Process(geom, clockData, detProp, ass, showers);
  }

  void ShowerRecoManager::Process(geo::GeometryCore const& geom,
                                  detinfo::DetectorClocksData const& clockData,
                                  detinfo::DetectorPropertiesData const& detProp,
                                  const ClusterAss_t& ass,
                                  std::vector<::recob::Shower>& showers)
  {

    for (auto const& pair : ass) {
      std::vector<::cluster::ClusterParamsAlg> cpans;

      cpans.reserve(pair.size());

      for (auto const& index : pair)
        cpans.push_back(fMatchMgr->GetInputClusters()[index]);

      fShowerAlgo->AppendInputClusters(cpans);
    }

    // Run shower reco
    showers = fShowerAlgo->Reconstruct(geom, clockData, detProp);
  }

}
