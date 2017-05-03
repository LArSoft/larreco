#include "TrajectoryMCSFitter.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

using namespace std;
using namespace trkf;
using namespace recob::tracking;

recob::MCSFitResult TrajectoryMCSFitter::fitMcs(const recob::TrackTrajectory& traj, bool momDepConst) const {
  //
  const double trajlen = traj.Length();
  const int nseg = std::max(minNSegs_,int(trajlen/segLen_));
  const double thisSegLen = trajlen/double(nseg);
  //std::cout << "track with length=" << trajlen << " broken in nseg=" << nseg << " of length=" << thisSegLen << " where segLen_=" << segLen_ << std::endl;
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
  int npoints = 0;
  while (nextValid!=recob::TrackTrajectory::InvalidIndex) {
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += ( (pos1-pos0).R() );
    pos0=pos1;
    npoints++;
    if (thislen>=thisSegLen) {
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
    cumseglens.push_back(cumseglens.back()+thislen);
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
  vector<double> cumLenFwd;
  vector<double> cumLenBwd;
  for (unsigned int i = 0; i<cumseglens.size()-2; i++) {
    cumLenFwd.push_back(cumseglens[i]);
    cumLenBwd.push_back(cumseglens.back()-cumseglens[i+2]);
  }
  const ScanResult fwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenFwd, true,  momDepConst);
  const ScanResult bwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenBwd, false, momDepConst);
  //
  return recob::MCSFitResult(pIdHyp_,
			     fwdResult.p,fwdResult.pUnc,fwdResult.logL,
			     bwdResult.p,bwdResult.pUnc,bwdResult.logL,
			     segradlengths,dtheta);
}

const TrajectoryMCSFitter::ScanResult TrajectoryMCSFitter::doLikelihoodScan(std::vector<double>& dtheta, std::vector<double>& seg_nradlengths, std::vector<double>& cumLen, bool fwdFit, bool momDepConst) const {
  int    best_idx  = -1;
  double best_logL = std::numeric_limits<double>::max();
  double best_p    = -1.0;
  std::vector<double> vlogL;
  for (double p_test = pMin_; p_test <= pMax_; p_test+=pStep_) {
    double logL = mcsLikelihood(p_test, angResol_, dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst);
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

double TrajectoryMCSFitter::mcsLikelihood(double p, double theta0x, std::vector<double>& dthetaij, std::vector<double>& seg_nradl, std::vector<double>& cumLen, bool fwd, bool momDepConst) const {
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
    const double m = mass();
    const double m2 = m*m;
    const double Etot = sqrt(p*p + m2);     //Initial
    //
    // test: MIP constant
    // constexpr double kcal2 = 0.002105*0.002105;
    // const double Eij = Etot - sqrt(kcal2*cumLen[i]*cumLen[i]);//energy at this segment
    //
    // test: Bethe-Bloch
    // const double elbb = energyLossBetheBloch(sqrt(m2),p);
    // const double Eij = Etot - sqrt(elbb*elbb*cumLen[i]*cumLen[i]);//energy at this segment
    //
    // MPV of Landau energy loss distribution
    const double elL = energyLossLandau(m,p,cumLen[i]);
    const double Eij = Etot - elL;//energy at this segment
    //
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

double TrajectoryMCSFitter::energyLossLandau(const double mass,const double p, const double x) const {
  //
  // eq. (33.11) in http://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (except density correction is ignored)
  //
  if (x<=0.) return 0.;
  constexpr double I = 188.E-6;
  constexpr double matConst = 1.4*18./40.;//density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  constexpr double j = 0.200;
  //
  double beta2 = p*p/(mass*mass+p*p);
  double gamma2 = 1./(1.0 - beta2);
  //
  double epsilon = 0.5*kappa*x*matConst/beta2;
  double deltaE = epsilon*( log(2.*me*beta2*gamma2/I) + log(epsilon/I) + j - beta2 );
  return deltaE*0.001;
}
//
//FIXME, test: to be removed eventually
double TrajectoryMCSFitter::energyLossBetheBloch(const double mass,const double p) const {
  // stolen, mostly, from GFMaterialEffects.
  const double charge(1.0);
  const double mEE(188.); // eV
  const double matZ(18.);
  const double matA(40.);
  const double matDensity(1.4);
  const double me(0.000511);
  //
  double beta = p/std::sqrt(mass*mass+p*p);
  double gammaSquare = 1./(1.0 - beta*beta);
  // 4pi.r_e^2.N.me = 0.307075, I think.
  double dedx = 0.307075*matDensity*matZ/matA/(beta*beta)*charge*charge;
  double massRatio = me/mass;
  // me=0.000511 here is in GeV. So mEE comes in here in eV.
  double argument = gammaSquare*beta*beta*me*1.E3*2./((1.E-6*mEE) * std::sqrt(1+2*std::sqrt(gammaSquare)*massRatio + massRatio*massRatio));
  //
  if (mass==0.0) return(0.0);
  if (argument <= exp(beta*beta))
    {
      dedx = 0.;
    }
  else{
    dedx *= (log(argument)-beta*beta); // Bethe-Bloch [MeV/cm]
    dedx *= 1.E-3;  // in GeV/cm, hence 1.e-3
    if (dedx<0.) dedx = 0.;
  }
  return dedx;
}
//
