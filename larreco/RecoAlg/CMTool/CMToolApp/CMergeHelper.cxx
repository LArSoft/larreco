#include "larreco/RecoAlg/CMTool/CMToolApp/CMergeHelper.h"

#include "larreco/RecoAlg/CMTool/CMToolBase/CMTException.h"

namespace cmtool {

  ::cmtool::CMergeManager& CMergeHelper::GetManager(size_t mgr_id)
  {
    if (_mgr_v.size() <= mgr_id) _mgr_v.resize(mgr_id + 1);
    return _mgr_v[mgr_id];
  }

  void CMergeHelper::SetAnaFile(TFile* fout)
  {
    for (auto& mgr : _mgr_v)
      mgr.SetAnaFile(fout);
  }

  void CMergeHelper::Process(util::GeometryUtilities const& gser,
                             const std::vector<std::vector<::util::PxHit>>& clusters)
  {
    _bk = ::cmtool::CMergeBookKeeper(clusters.size());

    for (size_t i = 0; i < _mgr_v.size(); ++i) {
      auto& mgr = _mgr_v[i];

      mgr.Reset();

      if (!i)
        mgr.SetClusters(gser, clusters);
      else
        mgr.SetClusters(_mgr_v[i - 1].GetClusters());

      mgr.Process(gser);

      auto const& new_bk = mgr.GetBookKeeper();

      if (!i)
        _bk = new_bk;
      else if (new_bk.GetResult().size() < new_bk.size())
        _bk.Combine(new_bk);
    }
  }

  const std::vector<::cluster::ClusterParamsAlg>& CMergeHelper::GetClusters() const
  {
    if (!(_mgr_v.size())) throw CMTException("No manager = no output clusters...");
    return _mgr_v.back().GetClusters();
  }
}
