//////////////////////////////////////////////////////////////////////
///
/// ClusterMatchTQ class
///
/// tjyang@fnal.gov
///
/// Algorithm for matching clusters between different views
/// based on time and charge information
///
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/ClusterMatchTQ.h"

#include "TH1D.h"

namespace {
  struct CluLen {
    int index;
    float length;
  };

  bool SortByLength(CluLen const& c1, CluLen const& c2) { return c1.length > c2.length; }

  bool SortByWire(art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2)
  {
    return h1->WireID().Wire < h2->WireID().Wire;
  }
}

namespace cluster {

  ClusterMatchTQ::ClusterMatchTQ(fhicl::ParameterSet const& pset)
  {
    fKSCut = pset.get<double>("KSCut");
    fEnableU = pset.get<bool>("EnableU");
    fEnableV = pset.get<bool>("EnableV");
    fEnableZ = pset.get<bool>("EnableZ");
  }

  //---------------------------------------------------------------------
  std::vector<std::vector<unsigned int>> ClusterMatchTQ::MatchedClusters(
    const detinfo::DetectorClocksData& clockdata,
    const detinfo::DetectorPropertiesData& detProp,
    const std::vector<art::Ptr<recob::Cluster>>& clusterlist,
    const art::FindManyP<recob::Hit>& fm) const
  {
    std::vector<std::vector<unsigned int>> matchedclusters;

    // get services
    art::ServiceHandle<geo::Geometry const> geom;
    int nplanes = geom->Nplanes();
    int nts = detProp.NumberTimeSamples();

    std::vector<std::vector<TH1D>> signals(nplanes);

    std::vector<std::vector<unsigned int>> Cls(nplanes);
    std::vector<std::vector<CluLen>> clulens(nplanes);

    for (size_t iclu = 0; iclu < clusterlist.size(); ++iclu) {

      float wire_pitch = geom->WirePitch(clusterlist[iclu]->Plane());

      float w0 = clusterlist[iclu]->StartWire();
      float w1 = clusterlist[iclu]->EndWire();
      float t0 = clusterlist[iclu]->StartTick();
      float t1 = clusterlist[iclu]->EndTick();

      CluLen clulen;
      clulen.index = iclu;

      auto const x0 = detProp.ConvertTicksToX(t0,
                                              clusterlist[iclu]->Plane().Plane,
                                              clusterlist[iclu]->Plane().TPC,
                                              clusterlist[iclu]->Plane().Cryostat);
      auto const x1 = detProp.ConvertTicksToX(t1,
                                              clusterlist[iclu]->Plane().Plane,
                                              clusterlist[iclu]->Plane().TPC,
                                              clusterlist[iclu]->Plane().Cryostat);
      clulen.length = std::hypot((w0 - w1) * wire_pitch, x0 - x1);

      switch (clusterlist[iclu]->View()) {
      case geo::kU:
        if (fEnableU) clulens[0].push_back(clulen);
        break;
      case geo::kV:
        if (fEnableV) clulens[1].push_back(clulen);
        break;
      case geo::kZ:
        if (fEnableZ) clulens[2].push_back(clulen);
        break;
      default: break;
      }
    }

    //sort clusters based on 2D length
    for (size_t i = 0; i < clulens.size(); ++i) {
      std::sort(clulens[i].begin(), clulens[i].end(), SortByLength);
      for (size_t j = 0; j < clulens[i].size(); ++j) {
        Cls[i].push_back(clulens[i][j].index);
      }
    }

    for (int i = 0; i < nplanes; ++i) {
      for (size_t ic = 0; ic < Cls[i].size(); ++ic) {
        TH1D sig(
          Form("sig_%d_%d", i, int(ic)), Form("sig_%d_%d", i, int(ic)), nts + 100, -100, nts);
        TH1D sigint(
          Form("sigint_%d_%d", i, int(ic)), Form("sigint_%d_%d", i, int(ic)), nts + 100, -100, nts);
        std::vector<art::Ptr<recob::Hit>> hitlist = fm.at(Cls[i][ic]);
        std::sort(hitlist.begin(), hitlist.end(), SortByWire);

        for (auto theHit = hitlist.begin(); theHit != hitlist.end(); theHit++) {

          double time = (*theHit)->PeakTime();
          time -= detProp.GetXTicksOffset(
            (*theHit)->WireID().Plane, (*theHit)->WireID().TPC, (*theHit)->WireID().Cryostat);

          double charge = (*theHit)->Integral();
          int bin = sig.FindBin(time);
          sig.SetBinContent(bin, sig.GetBinContent(bin) + charge);
          for (int j = bin; j <= sig.GetNbinsX(); ++j) {
            sigint.SetBinContent(j, sigint.GetBinContent(j) + charge);
          }
        }
        if (sigint.Integral()) sigint.Scale(1. / sigint.GetBinContent(sigint.GetNbinsX()));
        signals[i].push_back(sigint);
      }
    }

    //matching clusters between different views
    std::vector<int> matched(clusterlist.size());
    for (size_t i = 0; i < clusterlist.size(); ++i)
      matched[i] = 0;

    for (int i = 0; i < nplanes - 1; ++i) {
      for (int j = i + 1; j < nplanes; ++j) {
        for (size_t c1 = 0; c1 < Cls[i].size(); ++c1) {
          for (size_t c2 = 0; c2 < Cls[j].size(); ++c2) {

            // check if both are the same view
            if (clusterlist[Cls[i][c1]]->View() == clusterlist[Cls[j][c2]]->View()) continue;
            // check if both are in the same cryostat and tpc
            if (clusterlist[Cls[i][c1]]->Plane().Cryostat !=
                clusterlist[Cls[j][c2]]->Plane().Cryostat)
              continue;
            if (clusterlist[Cls[i][c1]]->Plane().TPC != clusterlist[Cls[j][c2]]->Plane().TPC)
              continue;
            // check if both are already in the matched list
            if (matched[Cls[i][c1]] == 1 && matched[Cls[j][c2]] == 1) continue;
            // KS test between two views in time
            double ks = 0;
            if (signals[i][c1].Integral() && signals[j][c2].Integral())
              ks = signals[i][c1].KolmogorovTest(&signals[j][c2]);
            else {
              mf::LogWarning("ClusterMatchTQ")
                << "One of the two clusters appears to be empty: " << clusterlist[Cls[i][c1]]->ID()
                << " " << clusterlist[Cls[j][c2]]->ID() << " " << i << " " << j << " " << c1 << " "
                << c2 << " " << signals[i][c1].Integral() << " " << signals[j][c2].Integral();
            }
            //hks->Fill(ks);
            int imatch = -1;   //track candidate index
            int iadd = -1;     //cluster index to be inserted
            if (ks > fKSCut) { //pass KS test
              // check both clusters with all matched clusters
              // if one is already matched,
              // check if need to add the other to the same track candidate
              for (size_t l = 0; l < matchedclusters.size(); ++l) {
                for (size_t m = 0; m < matchedclusters[l].size(); ++m) {
                  if (matchedclusters[l][m] == Cls[i][c1]) {
                    imatch = l; //track candidate
                    iadd = j;   //consider the other cluster
                  }
                  else if (matchedclusters[l][m] == Cls[j][c2]) {
                    imatch = l; //track candidate
                    iadd = i;   //consider the other cluster
                  }
                }
              }
              if (imatch >= 0) {
                if (iadd == i) {
                  bool matchview = false;
                  // check if one matched cluster has the same view
                  for (size_t ii = 0; ii < matchedclusters[imatch].size(); ++ii) {
                    if (clusterlist[matchedclusters[imatch][ii]]->View() ==
                        clusterlist[Cls[i][c1]]->View()) {
                      matchview = true;
                      //replace if the new cluster has more hits
                      if (fm.at(Cls[i][c1]).size() > fm.at(matchedclusters[imatch][ii]).size()) {
                        matched[matchedclusters[imatch][ii]] = 0;
                        matchedclusters[imatch][ii] = Cls[i][c1];
                        matched[Cls[i][c1]] = 1;
                      }
                    }
                  }
                  if (!matchview) { //not matched view found, just add
                    matchedclusters[imatch].push_back(Cls[i][c1]);
                    matched[Cls[i][c1]] = 1;
                  }
                }
                else {
                  bool matchview = false;
                  for (size_t jj = 0; jj < matchedclusters[imatch].size(); ++jj) {
                    if (clusterlist[matchedclusters[imatch][jj]]->View() ==
                        clusterlist[Cls[j][c2]]->View()) {
                      matchview = true;
                      //replace if it has more hits
                      if (fm.at(Cls[j][c2]).size() > fm.at(matchedclusters[imatch][jj]).size()) {
                        matched[matchedclusters[imatch][jj]] = 0;
                        matchedclusters[imatch][jj] = Cls[j][c2];
                        matched[Cls[j][c2]] = 1;
                      }
                    }
                  }
                  if (!matchview) {
                    matchedclusters[imatch].push_back(Cls[j][c2]);
                    matched[Cls[j][c2]] = 1;
                  }
                }
              }
              else {
                std::vector<unsigned int> tmp;
                tmp.push_back(Cls[i][c1]);
                tmp.push_back(Cls[j][c2]);
                matchedclusters.push_back(tmp);
                matched[Cls[i][c1]] = 1;
                matched[Cls[j][c2]] = 1;
              }
            } //pass KS test
          }   //c2
        }     //c1
      }       //j
    }         //i

    for (size_t i = 0; i < matchedclusters.size(); ++i) {
      if (matchedclusters[i].size())
        mf::LogVerbatim("ClusterMatchTQ") << "Cluster group " << i << ":";
      for (size_t j = 0; j < matchedclusters[i].size(); ++j) {
        mf::LogVerbatim("ClusterMatchTQ") << matchedclusters[i][j];
      }
    }

    return matchedclusters;
  }
} //namespace cluster
