#ifndef GEOMETRIC3DVERTEXFITTER_H
#define GEOMETRIC3DVERTEXFITTER_H

#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/VertexAssnMeta.h"
#include "larreco/RecoAlg/VertexWrapper.h"

namespace detinfo {
  class DetectorPropertiesData;
}

namespace trkf {

  using SMatrixSym22 = recob::tracking::SMatrixSym22;
  using SVector2 = recob::tracking::SVector2;
  using SMatrixSym33 = recob::tracking::SMatrixSym33;
  using SVector3 = recob::tracking::SVector3;

  /**
   * @file  larreco/RecoAlg/Geometric3DVertexFitter.h
   * @class trkf::Geometric3DVertexFitter
   *
   * @brief 3D vertex fitter based on the geometric properties (start position, direction, covariance) of the tracks.
   *
   * This algorithm fits vertices with following procedure.
   * First, tracks are sorted based on their start positions and the number of hits.
   * A vertex is created from the first two tracks: it is defined as the weighted average of the points of closest approaches of the two tracks.
   * Then the other tracks are added, to the vertex: the updated vertex is defined as the weighted average
   * of the n-1 track vertex position and the point of closest approach of the n-th track.
   * Methods to obtain the (unbiased) propagation distance, impact parameter, impact parameter error, impact parameter significance, and chi2
   * of a track with respect to the vertex are provided.
   *
   * Inputs are: a set of tracks; interface is provided allowing these to be passed directly of through a PFParticle hierarchy.
   *
   * Outputs are: a VertexWrapper, containing the vertex and the reference to the tracks actually used in the fit;
   * also methods to produce recob::VertexAssnMeta are provided.
   *
   * For configuration options see Geometric3DVertexFitter#Config
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  class Geometric3DVertexFitter {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> debugLevel{
        Name("debugLevel"),
        Comment("Debugging level: 0 for no printouts, 1 for minimal, 2 for full.")};
      fhicl::Atom<double> sipCut{
        Name("sipCut"),
        Comment(
          "Cut on maximum impact parameter significance to use the track in the vertex fit.")};
    };

    struct TracksFromVertexSorter {
      TracksFromVertexSorter(const recob::tracking::Point_t& vtxPos) : vtxPos_(vtxPos) {}
      bool
      operator()(std::reference_wrapper<const recob::Track> a,
                 std::reference_wrapper<const recob::Track> b) const
      {
        return ((a.get().Trajectory().Start() - vtxPos_).Mag2() <
                (b.get().Trajectory().Start() - vtxPos_).Mag2());
      }

    private:
      const recob::tracking::Point_t& vtxPos_;
    };

    struct ParsCovsOnPlane {
      ParsCovsOnPlane(const SVector2& p1,
                      const SVector2& p2,
                      const SMatrixSym22& c1,
                      const SMatrixSym22& c2,
                      const recob::tracking::Plane& p)
        : par1(p1), par2(p2), cov1(c1), cov2(c2), plane(p)
      {}
      SVector2 par1, par2;
      SMatrixSym22 cov1, cov2;
      recob::tracking::Plane plane;
    };

    // Constructor
    Geometric3DVertexFitter(const fhicl::Table<Config>& o,
                            const fhicl::Table<TrackStatePropagator::Config>& p)
      : debugLevel(o().debugLevel()), sipCut(o().sipCut())
    {
      prop = std::make_unique<TrackStatePropagator>(p);
    }

    VertexWrapper fitPFP(detinfo::DetectorPropertiesData const& detProp,
                         size_t iPF,
                         const art::ValidHandle<std::vector<recob::PFParticle>>& inputPFParticle,
                         const std::unique_ptr<art::FindManyP<recob::Track>>& assocTracks) const;
    VertexWrapper fitTracks(detinfo::DetectorPropertiesData const& detProp,
                            const std::vector<art::Ptr<recob::Track>>& arttracks) const;
    VertexWrapper fitTracks(detinfo::DetectorPropertiesData const& detProp,
                            TrackRefVec& tracks) const;
    VertexWrapper fitTracksWithVtx(detinfo::DetectorPropertiesData const& detProp,
                                   const std::vector<art::Ptr<recob::Track>>& tracks,
                                   const recob::tracking::Point_t& vtxPos) const;
    VertexWrapper fitTracksWithVtx(detinfo::DetectorPropertiesData const& detProp,
                                   TrackRefVec& tracks,
                                   const recob::tracking::Point_t& vtxPos) const;
    VertexWrapper closestPointAlongTrack(detinfo::DetectorPropertiesData const& detProp,
                                         const recob::Track& track,
                                         const recob::Track& other) const;
    VertexWrapper fitTwoTracks(detinfo::DetectorPropertiesData const& detProp,
                               const recob::Track& tk1,
                               const recob::Track& tk2) const;

    void addTrackToVertex(detinfo::DetectorPropertiesData const& detProp,
                          VertexWrapper& vtx,
                          const recob::Track& tk) const;

    std::vector<recob::VertexAssnMeta> computeMeta(detinfo::DetectorPropertiesData const& detProp,
                                                   const VertexWrapper& vtx);
    std::vector<recob::VertexAssnMeta> computeMeta(
      detinfo::DetectorPropertiesData const& detProp,
      const VertexWrapper& vtx,
      const std::vector<art::Ptr<recob::Track>>& arttracks);
    std::vector<recob::VertexAssnMeta> computeMeta(detinfo::DetectorPropertiesData const& detProp,
                                                   const VertexWrapper& vtx,
                                                   const TrackRefVec& trks);

    double chi2(detinfo::DetectorPropertiesData const& detProp,
                const VertexWrapper& vtx,
                const recob::Track& tk) const;
    double ip(detinfo::DetectorPropertiesData const& detProp,
              const VertexWrapper& vtx,
              const recob::Track& tk) const;
    double ipErr(detinfo::DetectorPropertiesData const& detProp,
                 const VertexWrapper& vtx,
                 const recob::Track& tk) const;
    double sip(detinfo::DetectorPropertiesData const& detProp,
               const VertexWrapper& vtx,
               const recob::Track& tk) const;
    double pDist(const VertexWrapper& vtx, const recob::Track& tk) const;

    VertexWrapper unbiasedVertex(detinfo::DetectorPropertiesData const& detProp,
                                 const VertexWrapper& vtx,
                                 const recob::Track& tk) const;
    double chi2Unbiased(detinfo::DetectorPropertiesData const& detProp,
                        const VertexWrapper& vtx,
                        const recob::Track& tk) const;
    double ipUnbiased(detinfo::DetectorPropertiesData const& detProp,
                      const VertexWrapper& vtx,
                      const recob::Track& tk) const;
    double ipErrUnbiased(detinfo::DetectorPropertiesData const& detProp,
                         const VertexWrapper& vtx,
                         const recob::Track& tk) const;
    double sipUnbiased(detinfo::DetectorPropertiesData const& detProp,
                       const VertexWrapper& vtx,
                       const recob::Track& tk) const;
    double pDistUnbiased(detinfo::DetectorPropertiesData const& detProp,
                         const VertexWrapper& vtx,
                         const recob::Track& tk) const;

  private:
    std::unique_ptr<TrackStatePropagator> prop;
    int debugLevel;
    double sipCut;

    double chi2(const ParsCovsOnPlane& pcp) const;
    double ip(const ParsCovsOnPlane& pcp) const;
    double ipErr(const ParsCovsOnPlane& pcp) const;
    double sip(const ParsCovsOnPlane& pcp) const;
    ParsCovsOnPlane getParsCovsOnPlane(detinfo::DetectorPropertiesData const& detProp,
                                       const trkf::VertexWrapper& vtx,
                                       const recob::Track& tk) const;
    std::pair<TrackState, double>
    weightedAverageState(ParsCovsOnPlane& pcop) const
    {
      return weightedAverageState(pcop.par1, pcop.par2, pcop.cov1, pcop.cov2, pcop.plane);
    };
    std::pair<TrackState, double> weightedAverageState(SVector2& par1,
                                                       SVector2& par2,
                                                       SMatrixSym22& cov1,
                                                       SMatrixSym22& cov2,
                                                       recob::tracking::Plane& target) const;
  };

}

#endif
