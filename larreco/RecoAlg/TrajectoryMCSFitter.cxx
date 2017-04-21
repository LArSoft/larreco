#include "TrajectoryMCSFitter.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

using namespace std;
using namespace trkf;
using namespace recob::tracking;

recob::MCSFitResult TrajectoryMCSFitter::fitMcs(const recob::TrackTrajectory& traj) const {
  //
  const double trajlen = traj.Length();
  const int nseg = trajlen/segLen_;
  const double thisSegLen = trajlen/double(nseg);
  std::cout << "track with length=" << trajlen << " broken in nseg=" << nseg << " of length=" << thisSegLen << " where segLen_=" << segLen_ << std::endl;
  //
  // Break the trajectory in segments of length approximately equal to segLen_
  //
  constexpr double lar_radl_inv = 1./14.0;
  //
  vector<size_t> breakpoints;
  vector<double> segradlengths;
  vector<double> cumseglens;
  cumseglens.push_back(0.);//first segment has zero cumulative length from previous segments
  double thislen = 0.;
  auto nextValid=traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  auto pos0 = traj.LocationAtPoint(nextValid);
  nextValid = traj.NextValidPoint(nextValid+1);
  //double totl = 0.;
  int npoints = 0;
  while (nextValid!=recob::TrackTrajectory::InvalidIndex) {
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += ( (pos1-pos0).R() );
    //const double dlen = (pos1-pos0).R();
    //totl += dlen;
    pos0=pos1;
    npoints++;
    if (thislen>=thisSegLen) {
      //cout << "thislen=" << thislen << " dlen=" << dlen << " npoints=" << npoints << " thisSegLen=" << thisSegLen << " totl=" << totl << endl;
      breakpoints.push_back(nextValid);
      if (npoints>=10) segradlengths.push_back(thislen*lar_radl_inv);
      else segradlengths.push_back(-999.);
      cumseglens.push_back(cumseglens.back()+thislen);
      thislen = 0.;
      npoints = 0;
    }
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //then add last segment
  if (thislen>0.) {
    breakpoints.push_back(nextValid);
    segradlengths.push_back(thislen*lar_radl_inv);
  }
  //
  // Fit segment directions, and get 3D angles between them
  //
  if (segradlengths.size()<2) return recob::MCSFitResult();
  vector<double> dtheta;
  Vector_t pcdir0;
  Vector_t pcdir1;
  for (unsigned int p = 0; p<segradlengths.size(); p++) {
    linearRegression(traj, breakpoints[p], breakpoints[p+1], pcdir1);
    if (p>0) {
      if (segradlengths[p]<-100. || segradlengths[p-1]<-100.) {
	dtheta.push_back(-999.);
      } else { 
	const double cosval = pcdir0.X()*pcdir1.X()+pcdir0.Y()*pcdir1.Y()+pcdir0.Z()*pcdir1.Z();
	assert(std::abs(cosval)<=1);
	//units are mrad
	double dt = 1000.*acos(cosval);//fixme, use expansion for small angles?
	dtheta.push_back(dt);
      }
    }
    pcdir0 = pcdir1;
  }
  //
  // Perform likelihood scan in forward and backward directions
  //
  constexpr double kcal2 = 0.002105*0.002105;//fixme input parameter?
  vector<double> dE2Fwd;
  vector<double> dE2Bwd;
  for (unsigned int i = 0; i<cumseglens.size()-1; i++) {
    dE2Fwd.push_back(kcal2*cumseglens[i]*cumseglens[i+1]);
    dE2Bwd.push_back(kcal2*(cumseglens.back()-cumseglens[i])*(cumseglens.back()-cumseglens[i+1]));
  }
  const ScanResult fwdResult = doLikelihoodScan(dtheta, segradlengths, dE2Fwd, true);
  const ScanResult bwdResult = doLikelihoodScan(dtheta, segradlengths, dE2Bwd, false);
  //
  return recob::MCSFitResult(pIdHyp_,
			     fwdResult.p,fwdResult.pUnc,fwdResult.logL,
			     bwdResult.p,bwdResult.pUnc,bwdResult.logL,
			     segradlengths,dtheta);
}

const TrajectoryMCSFitter::ScanResult TrajectoryMCSFitter::doLikelihoodScan(std::vector<double>& dtheta, std::vector<double>& seg_nradlengths, std::vector<double>& eLoss2, bool fwdFit) const {
  int    best_idx  = -1;
  double best_logL = std::numeric_limits<double>::max();
  double best_p    = -1.0;
  std::vector<double> vlogL;
  for (double p_test = pMin_; p_test <= pMax_; p_test+=pStep_) {
    double logL = mcsLikelihood(p_test, angResol_, dtheta, seg_nradlengths, eLoss2, fwdFit, true);
    if (logL < best_logL) {
      best_p    = p_test;
      best_logL = logL;
      best_idx  = vlogL.size();
    }
    vlogL.push_back(logL);
  }
  //
  //uncertainty from left side scan
  double lunc = -1.0;
  if (best_idx>0) {
    for (int j=best_idx-1;j>=0;j--) {
      double dLL = vlogL[j]-vlogL[best_idx];
      if ( dLL<0.5 ) {
	lunc = (best_idx-j)*pStep_;
      } else break;
    }
  }
  //uncertainty from right side scan
  double runc = -1.0;
  if (best_idx<int(vlogL.size()-1)) {  
    for (unsigned int j=best_idx+1;j<vlogL.size();j++) {
      double dLL = vlogL[j]-vlogL[best_idx];
      if ( dLL<0.5 ) {
	runc = (j-best_idx)*pStep_;
      } else break;
    }
  }
  return ScanResult(best_p, std::max(lunc,runc), best_logL);
}

void TrajectoryMCSFitter::linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
  //
  int npoints = 0;
  geo::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  while (nextValid!=lastPoint) {
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    nextValid = traj.NextValidPoint(nextValid+1);
    npoints++;
  }
  const auto avgpos = middlePointCalc.middlePoint();
  const double norm = 1./double(npoints);
  //
  assert(npoints>0);
  //
  TMatrixDSym m(3);
  nextValid = firstPoint;
  while (nextValid!=lastPoint) {
    const auto p = traj.LocationAtPoint(nextValid);
    const double xxw0 = p.X()-avgpos.X();
    const double yyw0 = p.Y()-avgpos.Y();
    const double zzw0 = p.Z()-avgpos.Z();
    m(0, 0) += xxw0*xxw0*norm;
    m(0, 1) += xxw0*yyw0*norm;
    m(0, 2) += xxw0*zzw0*norm;
    m(1, 0) += yyw0*xxw0*norm;
    m(1, 1) += yyw0*yyw0*norm;
    m(1, 2) += yyw0*zzw0*norm;
    m(2, 0) += zzw0*xxw0*norm;
    m(2, 1) += zzw0*yyw0*norm;
    m(2, 2) += zzw0*zzw0*norm;
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i=1; i<3; ++i) {
    if (eigenval(i)>maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  } 
  //
  pcdir = Vector_t(eigenvec(0, maxevalidx),eigenvec(1, maxevalidx),eigenvec(2, maxevalidx));
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir)<0.) pcdir*=-1.;
  //
}

double TrajectoryMCSFitter::mcsLikelihood(double p, double theta0x, std::vector<double>& dthetaij, std::vector<double>& seg_nradl, std::vector<double>& energyLoss2, bool fwd, bool momDepConst) const {
  //
  const int beg  = (fwd ? 0 : (dthetaij.size()-1));
  const int end  = (fwd ? dthetaij.size() : -1);
  const int incr = (fwd ? +1 : -1);
  //
  constexpr double fixedterm = 0.5 * std::log( 2.0 * M_PI );
  double result = std::abs(end-beg)*fixedterm;
  for (int i = beg; i != end; i+=incr ) {
    if (dthetaij[i]<0) {
      //cout << "skip segment with too few points" << endl;
      continue;
    }
    const double m2 = mass2();
    const double Etot = sqrt(p*p + m2);     //Initial
    const double Eij = Etot - sqrt(energyLoss2[i]);//energy at this segment
    const double Eij2 = Eij*Eij;
    if ( Eij2 < m2 ) {
      result = std::numeric_limits<double>::max();
      break;
    }
    const double pij = sqrt(Eij2 - m2);//momentum at this segment
    const double beta = sqrt( 1. - ((m2)/(pij*pij + m2)) );
    constexpr double tuned_HL_term1 = 11.0038;//11.17; fixme
    constexpr double HL_term2 = 0.038;
    const double tH0 = ( (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) / (pij*beta) ) * ( 1.0 + HL_term2 * std::log( seg_nradl[i] ) ) * sqrt( seg_nradl[i] );
    const double rms = sqrt( 2.0*( tH0 * tH0 + pow(theta0x, 2.0) ) );
    if (rms==0.0) {
      std::cout << " Error : RMS cannot be zero ! " << std::endl;
      return std::numeric_limits<double>::max();
    } 
    const double arg = dthetaij[i]/rms;
    result += ( std::log( rms ) + 0.5 * arg * arg );
  }
  return result;
}
