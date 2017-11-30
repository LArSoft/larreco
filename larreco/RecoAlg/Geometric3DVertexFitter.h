#ifndef GEOMETRIC3DVERTEXFITTER_H
#define GEOMETRIC3DVERTEXFITTER_H

////////////////////////////////////////////////////////////////////////
// Class:       Geometric3DVertexFitter
// File:        Geometric3DVertexFitter.h
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/VertexWrapper.h"
#include "lardataobj/RecoBase/TrackVertexMeta.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

namespace trkf {
  using SMatrixSym22 = recob::tracking::SMatrixSym22;
  using SVector2     = recob::tracking::SVector2;
  using SMatrixSym33 = recob::tracking::SMatrixSym33;
  using SVector3     = recob::tracking::SVector3;
  //
  class Geometric3DVertexFitter {
    //
  public:
    //
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> debugLevel {
	Name("debugLevel"),
	Comment("Debugging level: 0 for no printouts, 1 for minimal, 2 for full.")
      };
      fhicl::Atom<double> sipCut {
	Name("sipCut"),
	Comment("Cut on maximum impact parameter significance to use the track in the vertex fit.")
      };
    };
    //
    struct TracksFromVertexSorter {
      TracksFromVertexSorter( const recob::tracking::Point_t& vtxPos) : vtxPos_(vtxPos) {}
      bool operator()(std::reference_wrapper<const recob::Track> a, std::reference_wrapper<const recob::Track> b) const {
    	return ( (a.get().Trajectory().Start()-vtxPos_).Mag2()<(b.get().Trajectory().Start()-vtxPos_).Mag2() );
      }
      private:
      const recob::tracking::Point_t& vtxPos_;
    };
    //
    struct ParsCovsOnPlane {//fixme understand copies
    ParsCovsOnPlane(SVector2 p1, SVector2 p2, SMatrixSym22 c1, SMatrixSym22 c2, recob::tracking::Plane p)
    : par1(p1), par2(p2), cov1(c1), cov2(c2), plane(p) {}
      SVector2 par1, par2;
      SMatrixSym22 cov1, cov2;
      recob::tracking::Plane plane;
    };

    // Constructor
    Geometric3DVertexFitter(const fhicl::Table<Config>& o, const fhicl::Table<TrackStatePropagator::Config>& p)
      : debugLevel(o().debugLevel()), sipCut(o().sipCut())
      {
	prop = std::make_unique<TrackStatePropagator>(p);
      }

    VertexWrapper fitPFP(size_t iPF, const art::ValidHandle<std::vector<recob::PFParticle> >& inputPFParticle, 
			const std::unique_ptr<art::FindManyP<recob::Track> >& assocTracks) const;
    VertexWrapper fitTracks(const std::vector< art::Ptr<recob::Track> >& arttracks) const;
    VertexWrapper fitTracks(TrackRefVec& tracks) const;
    VertexWrapper fitTracksWithVtx(const std::vector< art::Ptr<recob::Track> >& tracks, const recob::tracking::Point_t& vtxPos) const;
    VertexWrapper fitTracksWithVtx(TrackRefVec& tracks, const recob::tracking::Point_t& vtxPos) const;
    VertexWrapper closestPointAlongTrack(const recob::Track& track, const recob::Track& other) const;
    VertexWrapper fitTwoTracks(const recob::Track& tk1, const recob::Track& tk2) const;
    //
    void addTrackToVertex(VertexWrapper& vtx, const recob::Track& tk) const;
    //
    std::vector<recob::TrackVertexMeta> computeMeta(const VertexWrapper& vtx);
    std::vector<recob::TrackVertexMeta> computeMeta(const VertexWrapper& vtx, const std::vector< art::Ptr<recob::Track> >& arttracks);
    std::vector<recob::TrackVertexMeta> computeMeta(const VertexWrapper& vtx, const TrackRefVec& trks);
    //
    double chi2 (const VertexWrapper& vtx, const recob::Track& tk) const;
    double ip   (const VertexWrapper& vtx, const recob::Track& tk) const;
    double ipErr(const VertexWrapper& vtx, const recob::Track& tk) const;
    double sip  (const VertexWrapper& vtx, const recob::Track& tk) const;
    double pDist(const VertexWrapper& vtx, const recob::Track& tk) const;
    //
    VertexWrapper unbiasedVertex(const VertexWrapper& vtx, const recob::Track& tk) const;
    double chi2Unbiased (const VertexWrapper& vtx, const recob::Track& tk) const;
    double ipUnbiased   (const VertexWrapper& vtx, const recob::Track& tk) const;
    double ipErrUnbiased(const VertexWrapper& vtx, const recob::Track& tk) const;
    double sipUnbiased  (const VertexWrapper& vtx, const recob::Track& tk) const;
    double pDistUnbiased(const VertexWrapper& vtx, const recob::Track& tk) const;
  private:
    std::unique_ptr<TrackStatePropagator> prop;
    int debugLevel;
    double sipCut;
    //
    double chi2 (const ParsCovsOnPlane& pcp) const;
    double ip   (const ParsCovsOnPlane& pcp) const;
    double ipErr(const ParsCovsOnPlane& pcp) const;
    double sip  (const ParsCovsOnPlane& pcp) const;
    ParsCovsOnPlane getParsCovsOnPlane(const trkf::VertexWrapper& vtx, const recob::Track& tk) const;
    std::pair<TrackState, double> weightedAverageState(ParsCovsOnPlane& pcop) const { return weightedAverageState(pcop.par1,pcop.par2,pcop.cov1,pcop.cov2,pcop.plane); };
    std::pair<TrackState, double> weightedAverageState(SVector2& par1, SVector2& par2, SMatrixSym22& cov1, SMatrixSym22& cov2, recob::tracking::Plane& target) const;
    //
  };
  //
}

#endif
