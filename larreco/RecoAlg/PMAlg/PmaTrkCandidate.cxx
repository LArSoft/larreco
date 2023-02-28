/**
 *  @file   PmaTrkCandidate.cxx
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Track finding helper for the Projection Matching Algorithm
 *
 *          Candidate for 3D track. Used to test 2D cluster associations, validadion result, MSE value.
 *          See PmaTrack3D.h file for details.
 */

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "larreco/RecoAlg/PMAlg/PmaHit3D.h"
#include "larreco/RecoAlg/PMAlg/PmaNode3D.h"
#include "larreco/RecoAlg/PMAlg/PmaSegment3D.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TVector3.h"

pma::TrkCandidate::TrkCandidate()
  : fParent(-1), fTrack(0), fKey(-1), fTreeId(-1), fMse(0), fValidation(0), fGood(false)
{}
// ------------------------------------------------------

pma::TrkCandidate::TrkCandidate(pma::Track3D* trk, int key, int tid)
  : fParent(-1), fTrack(trk), fKey(key), fTreeId(tid), fMse(0), fValidation(0), fGood(false)
{}
// ------------------------------------------------------

void pma::TrkCandidate::SetTrack(pma::Track3D* trk)
{
  if (fTrack) delete fTrack;
  fTrack = trk;
}
// ------------------------------------------------------

void pma::TrkCandidate::DeleteTrack()
{
  if (fTrack) delete fTrack;
  fTrack = 0;
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

int pma::TrkCandidateColl::getCandidateIndex(pma::Track3D const* candidate) const
{
  for (size_t t = 0; t < fCandidates.size(); ++t)
    if (fCandidates[t].Track() == candidate) return t;
  return -1;
}

int pma::TrkCandidateColl::getCandidateTreeId(pma::Track3D const* candidate) const
{
  int id = getCandidateIndex(candidate);
  if (id >= 0)
    return fCandidates[id].TreeId();
  else
    return -1;
}

void pma::TrkCandidateColl::setParentDaughterConnections()
{
  fParents.clear();

  size_t t = 0;
  while (t < fCandidates.size()) {
    if (fCandidates[t].IsValid()) {
      fCandidates[t].SetParent(-1);
      fCandidates[t].Daughters().clear();
      t++;
    }
    else
      fCandidates.erase(fCandidates.begin() + t);
  }

  for (t = 0; t < fCandidates.size(); ++t) {
    if (!fCandidates[t].IsValid()) continue;

    pma::Track3D const* trk = fCandidates[t].Track();
    pma::Node3D const* firstNode = trk->Nodes().front();
    if (firstNode->Prev()) // parent is a reconstructed track
    {
      pma::Track3D const* parentTrk = static_cast<pma::Segment3D*>(firstNode->Prev())->Parent();
      fCandidates[t].SetParent(getCandidateIndex(parentTrk));
    }
    else if (fCandidates[t].Parent() < 0) // parent not reconstructed and not yet set, add empty
                                          // candidate as a "primary" particle
    {
      fParents.push_back(pma::TrkCandidate());
      fParents.back().SetTreeId(fCandidates[t].TreeId());
      size_t pri_idx = fCandidates.size() + fParents.size() - 1;

      for (size_t i = 0; i < firstNode->NextCount(); ++i) {
        pma::Track3D const* daughterTrk =
          static_cast<pma::Segment3D*>(firstNode->Next(i))->Parent();
        int idx = getCandidateIndex(daughterTrk);
        if (idx >= 0) {
          fCandidates[(size_t)idx].SetParent(pri_idx);
          fParents.back().Daughters().push_back((size_t)idx);
        }
      }
    }

    for (size_t n = 1; n < trk->Nodes().size(); ++n) {
      auto node = trk->Nodes()[n];
      for (size_t i = 0; i < node->NextCount(); ++i) {
        pma::Track3D const* daughterTrk = static_cast<pma::Segment3D*>(node->Next(i))->Parent();
        if (daughterTrk != trk) {
          int idx = getCandidateIndex(daughterTrk);
          if (idx >= 0) fCandidates[t].Daughters().push_back((size_t)idx);
        }
      }
    }
  }
}
// ------------------------------------------------------

void pma::TrkCandidateColl::setTreeId(int id, size_t trkIdx, bool isRoot)
{
  pma::Track3D* trk = fCandidates[trkIdx].Track();
  pma::Node3D* vtx = trk->Nodes().front();
  pma::Segment3D* segThis = 0;
  pma::Segment3D* seg = 0;

  if (!isRoot) {
    segThis = trk->NextSegment(vtx);
    if (segThis) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  while (vtx) {
    segThis = trk->NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) {
        int idx = getCandidateIndex(seg->Parent());

        if (idx >= 0)
          setTreeId(id, idx, false);
        else
          mf::LogError("pma::setTreeId") << "Branch of the tree not found in tracks collection.";
      }
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  fCandidates[trkIdx].SetTreeId(id);
}

int pma::TrkCandidateColl::setTreeIds()
{
  for (auto& t : fCandidates)
    t.SetTreeId(-1);

  int id = 0;
  for (auto& t : fCandidates) {
    if (!t.IsValid() || (t.TreeId() >= 0)) continue;

    // index of a valid (reconstructed) track without reconstructed parent track
    int rootTrkIdx = getCandidateIndex(t.Track()->GetRoot());

    if (rootTrkIdx >= 0)
      setTreeId(id, rootTrkIdx);
    else
      mf::LogError("pma::setTreeIds") << "Root of the tree not found in tracks collection.";

    id++;
  }

  return id;
}
// ------------------------------------------------------

void pma::TrkCandidateColl::flipTreesToCoordinate(detinfo::DetectorPropertiesData const& detProp,
                                                  size_t coordinate)
{
  std::map<int, std::vector<pma::Track3D*>> toFlip;
  std::map<int, double> minVal;

  setTreeIds();
  for (auto& t : fCandidates) {
    if (!t.IsValid()) continue;

    int tid = t.TreeId();
    if (minVal.find(tid) == minVal.end()) minVal[tid] = 1.0e12;

    TVector3 pFront(t.Track()->front()->Point3D());
    pFront.SetX(-pFront.X());
    pFront.SetY(-pFront.Y());
    TVector3 pBack(t.Track()->back()->Point3D());
    pBack.SetX(-pBack.X());
    pBack.SetY(-pBack.Y());

    bool pushed = false;
    if (pFront[coordinate] < minVal[tid]) {
      minVal[tid] = pFront[coordinate];
      toFlip[tid].push_back(t.Track());
      pushed = true;
    }
    if (pBack[coordinate] < minVal[tid]) {
      minVal[tid] = pBack[coordinate];
      if (!pushed) toFlip[tid].push_back(t.Track());
    }
  }

  for (auto& tEntry : toFlip)
    if (tEntry.first >= 0) {
      size_t attempts = 0;
      while (!tEntry.second.empty()) {
        pma::Track3D* trk = tEntry.second.back();
        tEntry.second.pop_back();

        TVector3 pFront(trk->front()->Point3D());
        pFront.SetX(-pFront.X());
        pFront.SetY(-pFront.Y());
        TVector3 pBack(trk->back()->Point3D());
        pBack.SetX(-pBack.X());
        pBack.SetY(-pBack.Y());

        if (pFront[coordinate] > pBack[coordinate]) {
          if (setTreeOriginAtBack(detProp, trk)) { break; }
          else {
            mf::LogWarning("pma::TrkCandidateColl") << "Flip to coordinate failed.";
          }
        }
        else {
          setTreeOriginAtFront(detProp, trk);
          break; // good orientation, go to the next tree
        }

        if (attempts++ > 2) break; // do not try all the tracks in the queue...
      }
    }
}
// ------------------------------------------------------

bool pma::TrkCandidateColl::setTreeOriginAtFront(detinfo::DetectorPropertiesData const& detProp,
                                                 pma::Track3D* trk)
{
  int trkIdx = getCandidateIndex(trk);
  int treeId = getCandidateTreeId(trk);
  if (trkIdx < 0) {
    throw cet::exception("pma::TrkCandidateColl")
      << "Track not found in the collection." << std::endl;
  }

  bool done = true;
  pma::Node3D* n = trk->Nodes().front();
  if (n->Prev()) {
    pma::Segment3D* seg = static_cast<pma::Segment3D*>(n->Prev());
    pma::Track3D* incoming = seg->Parent();
    std::vector<pma::Track3D*> newTracks;
    if (incoming->NextSegment(n)) // upfff, need to break the parent track
    {
      int idx = incoming->index_of(n);
      if (idx >= 0) {
        pma::Track3D* u = incoming->Split(detProp, idx, false); // u is in front of incoming
        if (u) {
          newTracks.push_back(u);
          done = u->Flip(detProp, newTracks);
        }
        else {
          done = false;
        }
      }
      else {
        throw cet::exception("pma::Track3D") << "Node not found." << std::endl;
      }
    }
    else {
      done = incoming->Flip(detProp, newTracks);
    } // just flip incoming

    for (const auto ts : newTracks) {
      fCandidates.emplace_back(ts, -1, treeId);
    }
  }
  return done;
}
// ------------------------------------------------------

bool pma::TrkCandidateColl::setTreeOriginAtBack(detinfo::DetectorPropertiesData const& detProp,
                                                pma::Track3D* trk)
{
  int trkIdx = getCandidateIndex(trk);
  int treeId = getCandidateTreeId(trk);
  if (trkIdx < 0) {
    throw cet::exception("pma::TrkCandidateColl")
      << "Track not found in the collection." << std::endl;
  }

  pma::Track3D* incoming = fCandidates[trkIdx].Track();
  std::vector<pma::Track3D*> newTracks;
  bool done = incoming->Flip(detProp, newTracks);
  for (const auto ts : newTracks) {
    fCandidates.emplace_back(ts, -1, treeId);
  }
  return done;
}
// ------------------------------------------------------

void pma::TrkCandidateColl::flipTreesByDQdx()
{
  std::map<int, std::vector<pma::Track3D*>> trkMap;

  setTreeIds();
  for (auto const& t : fCandidates) {
    if (t.IsValid()) trkMap[t.TreeId()].push_back(t.Track());
  }

  for (auto& tEntry : trkMap) {
    std::sort(tEntry.second.begin(), tEntry.second.end(), pma::bTrack3DLonger());
    for (size_t i = tEntry.second.size(); i > 0; --i) {
      pma::Track3D* trk = tEntry.second[i - 1];
      if (trk->CanFlip()) { trk->AutoFlip(pma::Track3D::kForward, 0.15); }
    }
  }
}
// ------------------------------------------------------

void pma::TrkCandidateColl::merge(size_t idx1, size_t idx2)
{
  fCandidates[idx1].Track()->ExtendWith(fCandidates[idx2].Track()); // deletes track at idx2

  for (auto c : fCandidates[idx2].Clusters()) {
    fCandidates[idx1].Clusters().push_back(c);
  }

  fCandidates.erase(fCandidates.begin() + idx2);

  setTreeId(fCandidates[idx1].TreeId(), idx1);
}

pma::Track3D* pma::TrkCandidateColl::getTreeCopy(pma::TrkCandidateColl& dst,
                                                 size_t trkIdx,
                                                 bool isRoot)
{
  pma::Track3D* trk = fCandidates[trkIdx].Track();
  pma::Node3D* vtx = trk->Nodes().front();
  pma::Segment3D* segThis = 0;
  pma::Segment3D* seg = 0;

  int key = fCandidates[trkIdx].Key();
  int tid = fCandidates[trkIdx].TreeId();

  pma::Track3D* trkCopy = new pma::Track3D(*trk);
  pma::Node3D* vtxCopy = trkCopy->Nodes().front();
  pma::Segment3D* segThisCopy = 0;

  dst.tracks().emplace_back(trkCopy, key, tid);

  if (!isRoot) {
    segThis = trk->NextSegment(vtx);
    if (segThis) vtx = static_cast<pma::Node3D*>(segThis->Next());

    segThisCopy = trkCopy->NextSegment(vtxCopy);
    if (segThisCopy) vtxCopy = static_cast<pma::Node3D*>(segThisCopy->Next());
  }

  while (vtx) {
    segThis = trk->NextSegment(vtx);
    segThisCopy = trkCopy->NextSegment(vtxCopy);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) {
        int idx = getCandidateIndex(seg->Parent());

        if (idx >= 0) {
          pma::Track3D* branchCopy = getTreeCopy(dst, idx, false);
          if (!branchCopy->AttachTo(vtxCopy, true)) // no flip
            mf::LogError("pma::getTreeCopy") << "Branch copy cannot be attached to the tree.";
        }
        else
          mf::LogError("pma::getTreeCopy") << "Branch of the tree not found in source collection.";
      }
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;

    if (segThisCopy) vtxCopy = static_cast<pma::Node3D*>(segThisCopy->Next());
  }

  return trkCopy;
}
