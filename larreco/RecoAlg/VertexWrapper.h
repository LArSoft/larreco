#ifndef VERTEXWRAPPER_H
#define VERTEXWRAPPER_H

#include <functional>
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"

namespace trkf {
  
  // use reference_wrapper instead of pointers: we do not want ownership of the tracks
  typedef std::vector<std::reference_wrapper<const recob::Track> > TrackRefVec;

  class VertexWrapper {
    
  public:

    VertexWrapper() { vtx_ = recob::Vertex(); }
    VertexWrapper(const recob::Vertex& vtx) : vtx_(vtx) {}
    VertexWrapper(const recob::tracking::Point_t& pos, const recob::tracking::SMatrixSym33& cov, double chi2, int ndof) { vtx_ = recob::Vertex(pos,cov,chi2,ndof); }
    //
    const recob::Vertex& vertex() const { return vtx_; }
    bool isValid() const {return vtx_.isValid();}
    const recob::tracking::Point_t& position() const { return vtx_.position(); }
    const recob::tracking::SMatrixSym33& covariance() const { return vtx_.covariance(); }
    //
    void addTrack(const recob::Track& tk) { vtxtks_.push_back(tk); }
    void addTrackAndUpdateVertex(const recob::tracking::Point_t& pos, const recob::tracking::SMatrixSym33& cov,
    				 double chi2, int ndof, const recob::Track& tk) {
      vtx_ = recob::Vertex(pos, cov, vtx_.chi2()+chi2, vtx_.ndof()+ndof);
      addTrack(tk);
    }
    //
    size_t findTrack(const recob::Track& tk) const {
      for (size_t it = 0; it!=vtxtks_.size(); ++it) {
	if (&tk==&vtxtks_[it].get()) return it;
      }
      return vtxtks_.size();
    }
    //
    size_t tracksSize() const { return vtxtks_.size(); }
    const TrackRefVec& tracks() const { return vtxtks_; }
    TrackRefVec tracksWithoutElement(size_t element) const {
      TrackRefVec tks = vtxtks_;
      tks.erase(tks.begin()+element);
      return tks;
    }

  private:
    recob::Vertex vtx_;
    TrackRefVec vtxtks_;
  };
}

#endif
