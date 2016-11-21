#include "TrackKalmanFitter.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/RecoObjects/KHitWireX.h"
#include "lardata/RecoObjects/KHitContainerWireX.h"
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

bool trkf::TrackKalmanFitter::fitTrack(const recob::Track& track, const std::vector<art::Ptr<recob::Hit> >& hits, const double pval, const int pdgid,
				       recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits) {

  // LOG_DEBUG("TrackKalmanFitter") << "TMP POS --- " << track.Vertex().X() << " " << track.Vertex().Y() << " " << track.Vertex().Z();
  // LOG_DEBUG("TrackKalmanFitter") << "TMP DIR --- " << track.VertexDirection().X() << " " << track.VertexDirection().Y() << " " << track.VertexDirection().Z();
  // LOG_DEBUG("TrackKalmanFitter") << "TMP MOM --- " << track.VertexMomentum();

  auto position = track.Vertex();
  auto direction = track.VertexDirection();

  if (direction.Z()<0) {
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

  //fill the hit container and sort it according to the track direction
  trkf::KHitContainerWireX hitCont;
  //there must be a better way to convert a vector<Ptr> into a PtrVector
  art::PtrVector<recob::Hit> hitsv;
  for (auto hit : hits) {
    hitsv.push_back(hit);
  }
  hitCont.fill(hitsv, -1);
  hitCont.sort(vtxTrack,true,prop_,trkf::Propagator::FORWARD);
  auto sortedHits = hitCont.getSorted();

  //setup the KGTrack we'll ue for smoothing and to output the recob::Track
  trkf::KGTrack fittedTrack(hitCont.getPreferredPlane());

  //loop over groups and then on hits in groups
  for (auto itergroup : sortedHits) {
    for (auto ihit : itergroup.getHits()) {

      const trkf::KHitWireX* khitp = dynamic_cast<const trkf::KHitWireX*>(&*ihit);
      assert(khitp);
      trkf::KHitWireX khit(*khitp);//need a non const copy in case we want to modify the error
      if (useRMS_) khit.setMeasError(khit.getMeasError()*khit.getHit()->RMS()*khit.getHit()->RMS()/(khit.getHit()->SigmaPeakTime()*khit.getHit()->SigmaPeakTime()));
      // khit.setMeasError(khit.getMeasError()*25.);

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
	trkf::KHitTrack khitTrack(trf, *(itergroup.getHits().begin()));
	fittedTrack.addTrack(khitTrack);
      } else {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }

      if (trf.isValid()==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }

    } //for (auto ihit : itergroup.getHits())
  }//for (auto itergroup : sortedHits)

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
    if (useRMS_) khit.setMeasError(khit.getMeasError()*khit.getHit()->RMS()*khit.getHit()->RMS()/(khit.getHit()->SigmaPeakTime()*khit.getHit()->SigmaPeakTime()));

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
  fittedTrack.fillHits(outHits);

  return true;

}
