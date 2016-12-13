#include "TrackKalmanFitter.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/RecoObjects/SurfWireX.h"
#include "lardata/RecoObjects/KHitWireX.h"
#include "lardata/RecoObjects/Surface.h"
#include "lardata/RecoObjects/KGTrack.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/SurfYZPlane.h"
#include "lardata/RecoObjects/PropYZPlane.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

trkf::KTrack trkf::TrackKalmanFitter::convertRecobTrackIntoKTrack(const TVector3& pos, const TVector3&  dir,  const double pval, const int pdgid) {
  //build a KTrack on the surface orthogonal to the dir at the vertex
  const std::shared_ptr<const trkf::SurfXYZPlane> vtxsurf(new trkf::SurfXYZPlane(pos.X(),  pos.Y(),  pos.Z(), dir.X(), dir.Y(), dir.Z()));
  const double vtxsintheta = std::sin(vtxsurf->theta());
  const double vtxcostheta = std::cos(vtxsurf->theta());
  const double vtxsinphi = std::sin(vtxsurf->phi());
  const double vtxcosphi = std::cos(vtxsurf->phi());
  const double vtxpu = dir.X()*vtxcostheta + dir.Y()*vtxsintheta*vtxsinphi - dir.Z()*vtxsintheta*vtxcosphi;
  const double vtxpv = dir.Y()*vtxcosphi + dir.Z()*vtxsinphi;
  const double vtxpw = dir.X()*vtxsintheta - dir.Y()*vtxcostheta*vtxsinphi + dir.Z()*vtxcostheta*vtxcosphi;
  trkf::TrackVector vtxvec(5);
  vtxvec(0) = (pos.X()-vtxsurf->x0())*vtxcostheta + (pos.Y()-vtxsurf->y0())*vtxsintheta*vtxsinphi - (pos.Z()-vtxsurf->z0())*vtxsintheta*vtxcosphi;
  vtxvec(1) = (pos.Y()-vtxsurf->y0())*vtxcosphi + (pos.Z()-vtxsurf->z0())*vtxsinphi;
  vtxvec(2) = vtxpu/vtxpw;
  vtxvec(3) = vtxpv/vtxpw;
  vtxvec(4) = 1./pval;
  return trkf::KTrack(vtxsurf,vtxvec,trkf::Surface::FORWARD, pdgid);
}

struct CompareHits : std::binary_function<unsigned int, unsigned int, bool>
{
  CompareHits(const recob::Track& track, const TVector3& dir)
    : m_track(track), m_dir(dir)
  {}
  bool operator()(unsigned int Lhs, unsigned int Rhs)const
  {
    if (fabs(m_dir.X())>fabs(m_dir.Y()) && fabs(m_dir.X())>fabs(m_dir.Z())) {
      if (m_dir.X()>0)
	return m_track.LocationAtPoint(Lhs).X() < m_track.LocationAtPoint(Rhs).X();
      else
	return m_track.LocationAtPoint(Lhs).X() > m_track.LocationAtPoint(Rhs).X();
    } else if (fabs(m_dir.Y())>fabs(m_dir.X()) && fabs(m_dir.Y())>fabs(m_dir.Z())) {
      if (m_dir.Y()>0)
	return m_track.LocationAtPoint(Lhs).Y() < m_track.LocationAtPoint(Rhs).Y();
      else
	return m_track.LocationAtPoint(Lhs).Y() > m_track.LocationAtPoint(Rhs).Y();
    } else {
      if (m_dir.Z()>0)
	return m_track.LocationAtPoint(Lhs).Z() < m_track.LocationAtPoint(Rhs).Z();
      else
	return m_track.LocationAtPoint(Lhs).Z() > m_track.LocationAtPoint(Rhs).Z();
    }
  }
  const recob::Track& m_track;
  const TVector3& m_dir;
};

bool trkf::TrackKalmanFitter::fitTrack(const recob::Track& track, const std::vector<art::Ptr<recob::Hit> >& hits, 
				       const double pval, const int pdgid, const bool flipDirection,
				       recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits) {

  // LOG_DEBUG("TrackKalmanFitter") << "TMP POS --- " << track.Vertex().X() << " " << track.Vertex().Y() << " " << track.Vertex().Z();
  // LOG_DEBUG("TrackKalmanFitter") << "TMP DIR --- " << track.VertexDirection().X() << " " << track.VertexDirection().Y() << " " << track.VertexDirection().Z();
  // LOG_DEBUG("TrackKalmanFitter") << "TMP MOM --- " << track.VertexMomentum();

  auto position = track.Vertex();
  auto direction = track.VertexDirection();

  if (flipDirection) {
    position = track.End();
    direction = -track.EndDirection();
  }

  const trkf::KTrack vtxTrack = convertRecobTrackIntoKTrack(position, direction, pval, pdgid);

  if (vtxTrack.isValid()==0) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  //now add an error to this track: KETrack
  trkf::TrackError vtxerr;
  vtxerr.resize(5, false);
  vtxerr.clear();
  //if covariance exists use it (scaled by 100 to remove bias from previous fits), otherwise "dummy" initialization
  if (track.NumberCovariance()) {
    auto covariance = direction.Z()<0 ? track.VertexCovariance() : track.EndCovariance();
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
	vtxerr(i,j) = 100.*covariance(i,j);
      }
    }
  } else {
    vtxTrack.getSurface()->getStartingError(vtxerr);
  }
  const trkf::KETrack tre(vtxTrack, vtxerr);

  //setup the KFitTrack we'll use throughout the fit, it's initial status in INVALID
  trkf::KFitTrack trf = trkf::KFitTrack(tre, 0., 0., trkf::KFitTrack::INVALID);

  // LOG_DEBUG("TrackKalmanFitter") << "INITIAL TRACK\n" << trf.Print(std::cout);

  //figure out if hit vector is sorted along or opposite to track direction
  trkf::KFitTrack trfVtxF = trkf::KFitTrack(tre, 0., 0.);
  trkf::KFitTrack trfVtxB = trkf::KFitTrack(tre, 0., 0.);
  std::shared_ptr<const trkf::SurfWireX> vtxsurfF(new trkf::SurfWireX(hits.front()->WireID()));
  std::shared_ptr<const trkf::SurfWireX> vtxsurfB(new trkf::SurfWireX(hits.back()->WireID()));
  const KHitWireX hitF(hits.front(),vtxsurfF);
  const KHitWireX hitB(hits.back(),vtxsurfB);
  boost::optional<double> distF = prop_->err_prop(trfVtxF,hitF.getMeasSurface(),trkf::Propagator::UNKNOWN,false);
  boost::optional<double> distB = prop_->err_prop(trfVtxB,hitB.getMeasSurface(),trkf::Propagator::UNKNOWN,false);
  double dF = (distF ? *distF : 0);
  double dB = (distB ? *distB : 0);
  std::vector<art::Ptr<recob::Hit> > hitsv;
  if (dF>dB) for (auto hit = hits.rbegin(); hit<hits.rend(); ++hit) hitsv.push_back(*hit);
  else for (auto hit = hits.begin(); hit<hits.end(); ++hit) hitsv.push_back(*hit);

  //optional re-sorting, based on trajectry points and initial direction (there must be a 1-1 correspondance between hits and traj points)
  if (sortHits_) {
    assert(hitsv.size() == track.NumberTrajectoryPoints());
    std::vector<unsigned int> pos(hitsv.size());
    for (unsigned int i = 0; i != pos.size(); ++i){
      pos[i] = i;
    }
    std::sort(pos.begin(), pos.end(), CompareHits(track,direction));
    std::vector<art::Ptr<recob::Hit> > tmp;//is there a way to avoid using this tmp vector?
    for (unsigned int i = 0; i != pos.size(); ++i){
      tmp.push_back(hitsv[pos[i]]);
    }
    tmp.swap(hitsv);
  }

  //setup the KGTrack we'll use for smoothing and to output the recob::Track
  trkf::KGTrack fittedTrack(0);

    for (auto ihit : hitsv) {

      std::shared_ptr<const trkf::SurfWireX> hsurf(new trkf::SurfWireX(ihit->WireID()));
      trkf::KHitWireX khit(ihit,hsurf);
      if (useRMS_) khit.setMeasError(khit.getMeasError()*khit.getHit()->RMS()*khit.getHit()->RMS()/(std::max(khit.getHit()->SigmaPeakTime()*khit.getHit()->SigmaPeakTime(),0.08333333333f)));//0.0833333333 is 1/12, see KHitWireX.cxx line 69

      if (skipNegProp_) {
	auto trftmp = trf;
	boost::optional<double> pdisttest = prop_->noise_prop(trftmp,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
	if (pdisttest.get_value_or(-1.)<0.) {
	  mf::LogWarning("TrackKalmanFitter") << "WARNING: negative propagation distance. Skip this hit...";
	  continue;
	}
      }

      //propagate to measurement surface
      boost::optional<double> pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
      if (!pdist) {
	//in case of zero distance the prediction surface is the one from the original track
	//this is not good, we need it to be on the hit measurement surface, do a dummy propagation to fix it
	boost::optional<double> dist = prop_->err_prop(trf,khit.getMeasSurface(),trkf::Propagator::UNKNOWN,false);
      }
      trf.setStat(trkf::KFitTrack::FORWARD_PREDICTED);
      bool okpred = khit.predict(trf, prop_);
      if (khit.getPredSurface()!=khit.getMeasSurface()) {
	mf::LogWarning("TrackKalmanFitter") << "WARNING: khit.getPredSurface()!=khit.getMeasSurface(). Skip this hit...";
	continue;
      }

      //now update the forward fitted track
      if (okpred) {
	khit.update(trf);
	trf.setStat(trkf::KFitTrack::FORWARD);
	//store this track for the backward fit+smooth
	const std::shared_ptr< const KHitBase > strp(new trkf::KHitWireX(khit));
	trkf::KHitTrack khitTrack(trf, strp);
	fittedTrack.addTrack(khitTrack);
      } else {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }

      if (trf.isValid()==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }

    } //for (auto ihit : hitsv)
      /*
    } //for (auto ihit : itergroup.getHits())
  }//for (auto itergroup : hitsToProcess)
      */

  // LOG_DEBUG("TrackKalmanFitter") << "AFTER FORWARD\n" << trf.Print(std::cout);

  //reinitialize trf for backward fit
  trf.setError(100.*trf.getError());
  trf.setStat(trkf::KFitTrack::BACKWARD_PREDICTED);
  trf.setPath(0.);
  trf.setChisq(0.);

  //backward loop over track states and hits in fittedTrack: use hits for backward fit and fwd track states for smoothing
  for (auto itertrack = fittedTrack.getTrackMap().rbegin(); itertrack != fittedTrack.getTrackMap().rend(); ++itertrack) {
    trkf::KHitTrack& fwdTrack = itertrack->second;
    trkf::KHitWireX khit(dynamic_cast<const trkf::KHitWireX&>(*fwdTrack.getHit().get()));//need a non const copy in case we want to modify the error
    if (useRMS_) khit.setMeasError(khit.getMeasError()*khit.getHit()->RMS()*khit.getHit()->RMS()/(std::max(khit.getHit()->SigmaPeakTime()*khit.getHit()->SigmaPeakTime(),0.08333333333f)));

    boost::optional<double> pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::BACKWARD,true);
    if (!pdist) {
      //in case of zero distance the prediction surface is the one from the original track
      //this is not good, we need it to be on the hit measurement surface, do a dummy propagation to fix it
      boost::optional<double> dist = prop_->err_prop(trf,khit.getMeasSurface(),trkf::Propagator::UNKNOWN,false);
    }
    bool okpred = khit.predict(trf, prop_);
    if (okpred) {
      trf.setStat(trkf::KFitTrack::BACKWARD_PREDICTED);
      //combine forward updated and backward predicted, add this to the output track
      fwdTrack.combineFit(trf);
      //now update the backward fitted track
      khit.update(trf);
      trf.setStat(trkf::KFitTrack::BACKWARD);
    } else {
      mf::LogWarning("TrackKalmanFitter") << "WARNING invalid predicted surface";
    }
    if (trf.isValid()==0) {
      mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
      return false;
    }
  }//for (auto itertrack = fittedTrack.getTrackMap().rbegin(); itertrack != fittedTrack.getTrackMap().rend(); ++itertrack)

  // LOG_DEBUG("TrackKalmanFitter") << "AFTER BACKWARD\n" << trf.Print(std::cout);

  if (fittedTrack.getTrackMap().size()==0) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }
  
  bool zeromom = false;
  for (auto itertrack = fittedTrack.getTrackMap().rbegin(); itertrack != fittedTrack.getTrackMap().rend(); ++itertrack) {
    trkf::KHitTrack& trh = itertrack->second;
    double mom[3];
    trh.getMomentum(mom);
    double p = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
    if (p == 0.) zeromom = true;
  }
    if (zeromom) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
    return false;
  }

  //fill return objects with smoothed track and its hits
  fittedTrack.fillTrack(outTrack,track.ID());
  std::vector<unsigned int> hittpindex;
  fittedTrack.fillHits(outHits, hittpindex);

  return true;

}
