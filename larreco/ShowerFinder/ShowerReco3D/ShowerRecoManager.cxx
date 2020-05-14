#include "ShowerRecoManager.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/Geometry.h"
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
    //auto geom = ::larutil::Geometry::GetME();
    art::ServiceHandle<geo::Geometry const> geom;
    fMatchMgr = new ::cmtool::CMatchManager(geom->Nplanes());
  }

  void
  ShowerRecoManager::Reset()
  {
    if (fShowerAlgo) fShowerAlgo->Reset();
    fMatchMgr->Reset();
  }

  ClusterAss_t
  ShowerRecoManager::Reconstruct(const std::vector<std::vector<util::PxHit>>& clusters,
                                 std::vector<::recob::Shower>& showers)
  {
    showers.clear();
    fMatchMgr->SetClusters(clusters);

    ClusterAss_t res_ass;
    // Run matching & retrieve matched cluster indices
    try {
      fMatchMgr->Process();
    }
    catch (::cmtool::CMTException& e) {
      e.what();
      return res_ass;
    }
    res_ass = fMatchMgr->GetBookKeeper().GetResult();

    Process(res_ass, showers);

    return res_ass;
  }

  void
  ShowerRecoManager::Reconstruct(const std::vector<std::vector<util::PxHit>>& clusters,
                                 const ClusterAss_t& ass,
                                 std::vector<::recob::Shower>& showers)
  {
    showers.clear();
    fMatchMgr->SetClusters(clusters);

    Process(ass, showers);
  }

  void
  ShowerRecoManager::Process(const ClusterAss_t& ass, std::vector<::recob::Shower>& showers)
  {

    for (auto const& pair : ass) {

      std::vector<::cluster::ClusterParamsAlg> cpans;

      cpans.reserve(pair.size());

      for (auto const& index : pair)

        cpans.push_back(fMatchMgr->GetInputClusters()[index]);

      fShowerAlgo->AppendInputClusters(cpans);
    }

    // Run shower reco
    showers = fShowerAlgo->Reconstruct();
  }

}
