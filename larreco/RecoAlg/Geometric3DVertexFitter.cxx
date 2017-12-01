#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

trkf::VertexWrapper trkf::Geometric3DVertexFitter::fitPFP(size_t iPF, const art::ValidHandle<std::vector<recob::PFParticle> >& inputPFParticle, 
							 const std::unique_ptr<art::FindManyP<recob::Track> >& assocTracks) const
{
  using namespace std;
  //
  art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPF);
  //
  if (debugLevel>1) std::cout << "pfp#" << iPF << " PdgCode=" << pfp->PdgCode() 
			      << " IsPrimary=" << pfp->IsPrimary()
			      << " NumDaughters=" << pfp->NumDaughters()
			      << std::endl;
  if (pfp->IsPrimary()==false || pfp->NumDaughters()<2) VertexWrapper();
  
  TrackRefVec tracks;
  
  auto& pfd = pfp->Daughters();
  for (auto ipfd : pfd) {
    vector< art::Ptr<recob::Track> > pftracks = assocTracks->at(ipfd);
    for (auto t : pftracks) {
      tracks.push_back(*t);
    }
  }
  if (tracks.size()<2) return VertexWrapper();
  //
  return fitTracks(tracks);
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::fitTracks(const std::vector< art::Ptr<recob::Track> >& arttracks) const
{
  TrackRefVec tracks;
  for (auto t : arttracks) {
    tracks.push_back(*t);
  }
  //
  return fitTracks(tracks);
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::fitTracks(TrackRefVec& tracks) const
{
  if (debugLevel>0) std::cout << "fitting vertex with ntracks=" << tracks.size() << std::endl;
  if (tracks.size()<2) return VertexWrapper();
  // sort tracks by number of hits
  std::sort(tracks.begin(), tracks.end(), [](std::reference_wrapper<const recob::Track> a, std::reference_wrapper<const recob::Track> b) {
      return a.get().CountValidPoints() > b.get().CountValidPoints();}
    );
  //find pair with closest start positions and put them upfront
  unsigned int tk0 = tracks.size();
  unsigned int tk1 = tracks.size();
  float mind = FLT_MAX;
  for (unsigned int i=0;i<tracks.size()-1;i++) {
    for (unsigned int j=i+1;j<tracks.size();j++) {
      float d = (tracks[i].get().Trajectory().Start()-tracks[j].get().Trajectory().Start()).Mag2();
      if (debugLevel>1) std::cout << "test i=" << i << " start=" << tracks[i].get().Trajectory().Start() << " j=" << j << " start=" << tracks[j].get().Trajectory().Start() << " d=" << d << " mind=" << mind << " tk0=" << tk0 << " tk1=" << tk1 << std::endl;
      if (d<mind) {
	mind=d;
	tk0 = i;
	tk1 = j;
      }
    }
  }
  if (tk0>1 || tk1>1) {
    if (tk0>1 && tk1>1) {
      std::swap(tracks[0],tracks[tk0]);
      std::swap(tracks[1],tracks[tk1]);
    }
    if (tk0==0) std::swap(tracks[1],tracks[tk1]);
    if (tk0==1) std::swap(tracks[0],tracks[tk1]);
    if (tk1==0) std::swap(tracks[1],tracks[tk0]);
    if (tk1==1) std::swap(tracks[0],tracks[tk0]);
  }
  //
  // find vertex between the first two tracks
  VertexWrapper vtx = fitTwoTracks(tracks[0], tracks[1]);
  if (vtx.isValid()==false || vtx.tracksSize()<2) return vtx;
  //
  // then add other tracks and update vertex measurement
  for (auto tk = tracks.begin()+2; tk<tracks.end(); ++tk) {
    auto sipv = sip(vtx, *tk);
    if (debugLevel>1) std::cout << "sip=" << sipv << std::endl;
    if (sipv>sipCut) continue;
    addTrackToVertex(vtx, *tk);
  }
  return vtx;
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::fitTracksWithVtx(const std::vector< art::Ptr<recob::Track> >& arttracks,
								    const recob::tracking::Point_t& vtxPos) const
{
  TrackRefVec tracks;
  for (auto t : arttracks) {
    tracks.push_back(*t);
  }
  //
  return fitTracksWithVtx(tracks,vtxPos);
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::fitTracksWithVtx(TrackRefVec& tracks,
								    const recob::tracking::Point_t& vtxPos) const
{
  if (debugLevel>0) std::cout << "fitting vertex with ntracks=" << tracks.size() << std::endl;
  if (tracks.size()<2) return VertexWrapper();
  // sort tracks by proximity to input vertex
  std::sort(tracks.begin(), tracks.end(), TracksFromVertexSorter(vtxPos) );
  //
  // find vertex between the first two tracks
  VertexWrapper vtx = fitTwoTracks(tracks[0], tracks[1]);
  if (vtx.isValid()==false || vtx.tracks().size()<2) return vtx;
  //
  // then add other tracks and update vertex measurement
  for (auto tk = tracks.begin()+2; tk<tracks.end(); ++tk) {
    auto sipv = sip(vtx, *tk);
    if (debugLevel>1) std::cout << "sip=" << sipv << std::endl;
    if (sipv>sipCut) continue;
    addTrackToVertex(vtx, *tk);
  }
  return vtx;
}

std::pair<trkf::TrackState, double> trkf::Geometric3DVertexFitter::weightedAverageState(SVector2& par1, SVector2& par2,
											SMatrixSym22& cov1, SMatrixSym22& cov2,
											recob::tracking::Plane& target) const
{
  SVector2 deltapar = par2 - par1;
  SMatrixSym22 covsum = (cov2+cov1);
  //
  if (debugLevel>1) {
    std::cout << "par1=" << par1 << std::endl;
    std::cout << "par2=" << par2 << std::endl;
    std::cout << "deltapar=" << deltapar << std::endl;
    //
    std::cout << "cov1=\n" << cov1 << std::endl;
    std::cout << "cov2=\n" << cov2 << std::endl;
    std::cout << "covsum=\n" << covsum << std::endl;
  }
  //
  if (debugLevel>1) {
    double det1;
    bool d1ok = cov1.Det2(det1);
    std::cout << "cov1 det=" << det1 << " ok=" << d1ok << std::endl;
    double det2;
    bool d2ok = cov2.Det2(det2);
    std::cout << "cov2 det=" << det2 << " ok=" << d2ok << std::endl;
    double detsum;
    bool dsok = covsum.Det2(detsum);
    std::cout << "covsum det=" << detsum << " ok=" << dsok << std::endl;
  }
  //
  bool invertok = covsum.Invert();
  if (!invertok) {
    SVector5 vtxpar5(0,0,0,0,0);
    SMatrixSym55 vtxcov55;
    TrackState vtxstate(vtxpar5, vtxcov55, target, true, 0);
    return std::make_pair(vtxstate, util::kBogusD);
  }
  auto k = cov1*covsum;

  if (debugLevel>1) {
    std::cout << "inverted covsum=\n" << covsum << std::endl;
    std::cout << "k=\n" << k << std::endl;
  }

  SVector2 vtxpar2 = par1 + k*deltapar;
  SMatrixSym22 vtxcov22 = cov1 - ROOT::Math::SimilarityT(cov1,covsum);

  if (debugLevel>1) {
    std::cout << "vtxpar2=" << vtxpar2 << std::endl;
    std::cout << "vtxcov22=\n" << vtxcov22 << std::endl;
  }
  
  auto chi2 = ROOT::Math::Similarity(deltapar,covsum);

  SVector5 vtxpar5(vtxpar2[0],vtxpar2[1],0,0,0);
  SMatrixSym55 vtxcov55;vtxcov55.Place_at(vtxcov22,0,0);
  TrackState vtxstate(vtxpar5, vtxcov55, target, true, 0);

  return std::make_pair(vtxstate, chi2);
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::closestPointAlongTrack(const recob::Track& track, const recob::Track& other) const
{
  // find the closest approach point along track
  const auto& start1 = track.Trajectory().Start();
  const auto& start2 = other.Trajectory().Start();

  const auto dir1 = track.Trajectory().StartDirection();
  const auto dir2 = other.Trajectory().StartDirection();

  if (debugLevel>0) {
    std::cout << "track1 with start=" << start1 << " dir=" << dir1
	      << " length=" << track.Length() << " points=" << track.CountValidPoints()
	      << std::endl;
    std::cout << "covariance=\n" << track.VertexCovarianceGlobal6D() << std::endl;
    std::cout << "track2 with start=" << start2 << " dir=" << dir2
              << " length=" << other.Length() << " points=" << other.CountValidPoints()
	      << std::endl;
    std::cout << "covariance=\n" << other.VertexCovarianceGlobal6D() << std::endl;
  }

  const auto dpos = start1-start2;
  const auto dotd1d2 = dir1.Dot(dir2);
  const auto dotdpd1 = dpos.Dot(dir1);
  const auto dotdpd2 = dpos.Dot(dir2);
  const auto dist2 = ( dotd1d2*dotdpd1 - dotdpd2 )/( dotd1d2*dotd1d2 - 1 );
  const auto dist1 = ( dotd1d2*dist2 - dotdpd1 );

  if (debugLevel>0) {
    std::cout << "track1 pca=" << start1+dist1*dir1 << " dist=" << dist1 << std::endl;
    std::cout << "track2 pca=" << start2+dist2*dir2 << " dist=" << dist2 << std::endl;
  }

  //by construction both point of closest approach on the two lines lie on this plane
  recob::tracking::Plane target(start1+dist1*dir1, dir1);

  recob::tracking::Plane plane1(start1, dir1);
  trkf::TrackState state1(track.VertexParametersLocal5D(), track.VertexCovarianceLocal5D(), plane1, true, track.ParticleId());
  bool propok1 = true;
  state1 = prop->propagateToPlane(propok1, state1, target, true, true, trkf::TrackStatePropagator::UNKNOWN);
  if (!propok1) {
    std::cout << "failed propagation, return track1 start pos=" << track.Start() << std::endl;
    VertexWrapper vtx;
    vtx.addTrackAndUpdateVertex(track.Start(), track.VertexCovarianceGlobal6D().Sub<SMatrixSym33>(0,0), 0, 0, track);
    return vtx;
  } else {
    VertexWrapper vtx;
    vtx.addTrackAndUpdateVertex(state1.position(), state1.covariance6D().Sub<SMatrixSym33>(0,0), 0, 0, track);
    return vtx;
  }
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::fitTwoTracks(const recob::Track& tk1, const recob::Track& tk2) const
{
  // find the closest approach points
  auto start1 = tk1.Trajectory().Start();
  auto start2 = tk2.Trajectory().Start();

  auto dir1 = tk1.Trajectory().StartDirection();
  auto dir2 = tk2.Trajectory().StartDirection();

  if (debugLevel>0) {
    std::cout << "track1 with start=" << start1 << " dir=" << dir1
	      << " length=" << tk1.Length() << " points=" << tk1.CountValidPoints()
	      << std::endl;
    std::cout << "covariance=\n" << tk1.VertexCovarianceGlobal6D() << std::endl;
    std::cout << "track2 with start=" << start2 << " dir=" << dir2
              << " length=" << tk2.Length() << " points=" << tk2.CountValidPoints()
  	      << std::endl;
    std::cout << "covariance=\n" << tk2.VertexCovarianceGlobal6D() << std::endl;
  }

  auto dpos = start1-start2;
  auto dotd1d2 = dir1.Dot(dir2);

  auto dotdpd1 = dpos.Dot(dir1);
  auto dotdpd2 = dpos.Dot(dir2);

  auto dist2 = ( dotd1d2*dotdpd1 - dotdpd2 )/( dotd1d2*dotd1d2 - 1 );
  auto dist1 = ( dotd1d2*dist2 - dotdpd1 );

  //by construction both point of closest approach on the two lines lie on this plane
  recob::tracking::Plane target(start1+dist1*dir1, dir1);

  recob::tracking::Plane plane1(start1, dir1);
  trkf::TrackState state1(tk1.VertexParametersLocal5D(), tk1.VertexCovarianceLocal5D(), plane1, true, tk1.ParticleId());
  bool propok1 = true;
  state1 = prop->propagateToPlane(propok1, state1, target, true, true, trkf::TrackStatePropagator::UNKNOWN);
  if (!propok1) {
    std::cout << "failed propagation, return track1 start pos=" << tk1.Start() << std::endl;
    VertexWrapper vtx;
    vtx.addTrackAndUpdateVertex(tk1.Start(), tk1.VertexCovarianceGlobal6D().Sub<SMatrixSym33>(0,0), 0, 0, tk1);
    return vtx;
  }

  recob::tracking::Plane plane2(start2, dir2);
  trkf::TrackState state2(tk2.VertexParametersLocal5D(), tk2.VertexCovarianceLocal5D(), plane2, true, tk2.ParticleId());
  bool propok2 = true;
  state2 = prop->propagateToPlane(propok2, state2, target, true, true, trkf::TrackStatePropagator::UNKNOWN);
  if (!propok2) {
    std::cout << "failed propagation, return track1 start pos=" << tk1.Start() << std::endl;
    VertexWrapper vtx;
    vtx.addTrackAndUpdateVertex(tk1.Start(), tk1.VertexCovarianceGlobal6D().Sub<SMatrixSym33>(0,0), 0, 0, tk1);
    return vtx;
  }

  if (debugLevel>0) {
    std::cout << "track1 pca=" << start1+dist1*dir1 << " dist=" << dist1 << std::endl;
    std::cout << "track2 pca=" << start2+dist2*dir2 << " dist=" << dist2 << std::endl;
  }

  //here we neglect that to find dist1 and dist2 one track depends on the other...
  SMatrixSym22 cov1 = state1.covariance().Sub<SMatrixSym22>(0,0);
  SMatrixSym22 cov2 = state2.covariance().Sub<SMatrixSym22>(0,0);

  if (debugLevel>1) {
    std::cout << "target pos=" << target.position() << " dir=" << target.direction() << std::endl;
    //test if orthogonal
    auto dcp = state1.position()-state2.position();
    std::cout << "dot dcp-dir1=" << dcp.Dot(tk1.Trajectory().StartDirection()) << std::endl;
    std::cout << "dot dcp-dir2=" << dcp.Dot(tk2.Trajectory().StartDirection()) << std::endl;
    //
    std::cout << "cov1=" << cov1 << std::endl;
    std::cout << "cov2=" << cov2 << std::endl;
  }

  SVector2 par1(state1.parameters()[0],state1.parameters()[1]);
  SVector2 par2(state2.parameters()[0],state2.parameters()[1]);

  std::pair<TrackState, double> was = weightedAverageState(par1, par2, cov1, cov2, target);
  if (was.second <= (util::kBogusD-1)) {
    std::cout << "failed inversion, return track1 start pos=" << tk1.Start() << std::endl;
    VertexWrapper vtx;
    vtx.addTrackAndUpdateVertex(tk1.Start(), tk1.VertexCovarianceGlobal6D().Sub<SMatrixSym33>(0,0), 0, 0, tk1);
    return vtx;
  }

  const int ndof = 1; //1=2*2-3; each measurement is 2D because it is defined on a plane
  trkf::VertexWrapper vtx(was.first.position(), was.first.covariance6D().Sub<SMatrixSym33>(0,0), was.second, ndof);
  vtx.addTrack(tk1);
  vtx.addTrack(tk2);

  if (debugLevel>0) {
    std::cout << "vtxpos=" << vtx.position() << std::endl;
    std::cout << "vtxcov=\n" << vtx.covariance() << std::endl;
    std::cout << "chi2=" << was.second << std::endl;
  }

  return vtx;
}

trkf::Geometric3DVertexFitter::ParsCovsOnPlane trkf::Geometric3DVertexFitter::getParsCovsOnPlane(const trkf::VertexWrapper& vtx, const recob::Track& tk) const {
  auto start = tk.Trajectory().Start();
  auto dir = tk.Trajectory().StartDirection();

  auto vtxpos = vtx.position();
  auto vtxcov = vtx.covariance();

  auto dist = dir.Dot(vtxpos-start);

  //by construction vtxpos also lies on this plane
  recob::tracking::Plane target(start+dist*dir, dir);

  recob::tracking::Plane plane(start, dir);
  trkf::TrackState state(tk.VertexParametersLocal5D(), tk.VertexCovarianceLocal5D(), plane, true, tk.ParticleId());
  bool propok = true;
  state = prop->propagateToPlane(propok, state, target, true, true, trkf::TrackStatePropagator::UNKNOWN);

  if (debugLevel>0) {
    std::cout << "input vtx=" << vtxpos << std::endl;
    std::cout << "input cov=\n" << vtxcov << std::endl;
    std::cout << "track pca=" << start+dist*dir << " dist=" << dist << std::endl;
    std::cout << "target pos=" << target.position() << " dir=" << target.direction() << std::endl;
  }

  //rotate the vertex position and covariance to the target plane
  SVector6 vtxpar6(vtxpos.X(),vtxpos.Y(),vtxpos.Z(),0,0,0);
  SMatrixSym66 vtxcov66;vtxcov66.Place_at(vtxcov,0,0);

  auto vtxpar5 = target.Global6DToLocal5DParameters(vtxpar6);
  SVector2 par1(vtxpar5[0],vtxpar5[1]);
  SVector2 par2(state.parameters()[0],state.parameters()[1]);

  //here we neglect that to find dist, the vertex is used...
  SMatrixSym22 cov1 = target.Global6DToLocal5DCovariance(vtxcov66,false,Vector_t()).Sub<SMatrixSym22>(0,0);
  SMatrixSym22 cov2 = state.covariance().Sub<SMatrixSym22>(0,0);

  if (debugLevel>1) {
    std::cout << "vtxpar5=" << vtxpar5 << std::endl;
    std::cout << "state.parameters()=" << state.parameters() << std::endl;
    std::cout << "vtx covariance=\n" << vtxcov66 << std::endl;
    std::cout << "trk covariance=\n" << state.covariance6D() << std::endl;
    std::cout << "cov1=\n" << cov1 << std::endl;
    std::cout << "cov2=\n" << cov2 << std::endl;
  }

  return ParsCovsOnPlane(par1,par2,cov1,cov2,target);
}

void trkf::Geometric3DVertexFitter::addTrackToVertex(trkf::VertexWrapper& vtx, const recob::Track& tk) const
{

  if (debugLevel>0) {
    std::cout << "adding track with start=" << tk.Start() << " dir=" << tk.StartDirection()
	      << " length=" << tk.Length() << " points=" << tk.CountValidPoints()
	      << std::endl;
    std::cout << "covariance=\n" << tk.VertexCovarianceGlobal6D() << std::endl;
  }

  ParsCovsOnPlane pcp = getParsCovsOnPlane(vtx, tk);
  std::pair<TrackState, double> was = weightedAverageState(pcp);
  if (was.second <= (util::kBogusD-1.)) {
    return;
  }

  const int ndof = 2;// Each measurement is 2D because it is defined on a plane
  vtx.addTrackAndUpdateVertex(was.first.position(), was.first.covariance6D().Sub<SMatrixSym33>(0,0), was.second, ndof, tk);

  if (debugLevel>0) {
    std::cout << "updvtxpos=" << vtx.position() << std::endl;
    std::cout << "updvtxcov=\n" << vtx.covariance() << std::endl;
    std::cout << "add chi2=" << was.second << std::endl;
  }

  return;  
}

double trkf::Geometric3DVertexFitter::chi2(const trkf::Geometric3DVertexFitter::ParsCovsOnPlane& pcp) const
{
  const SVector2 deltapar = pcp.par2 - pcp.par1;
  SMatrixSym22 covsum = (pcp.cov2+pcp.cov1);
  //
  bool invertok = covsum.Invert();
  if (!invertok) return -1.;
  //
  return ROOT::Math::Similarity(deltapar,covsum);
}

double trkf::Geometric3DVertexFitter::chi2(const VertexWrapper& vtx, const recob::Track& tk) const 
{
  return chi2(getParsCovsOnPlane(vtx, tk));
}

double trkf::Geometric3DVertexFitter::ip(const trkf::Geometric3DVertexFitter::ParsCovsOnPlane& pcp) const
{
  const SVector2 deltapar = pcp.par2 - pcp.par1;
  return std::sqrt( deltapar[0]*deltapar[0] + deltapar[1]*deltapar[1] );
}

double trkf::Geometric3DVertexFitter::ip(const VertexWrapper& vtx, const recob::Track& tk) const
{
  return ip(getParsCovsOnPlane(vtx, tk));
}

double trkf::Geometric3DVertexFitter::ipErr(const trkf::Geometric3DVertexFitter::ParsCovsOnPlane& pcp) const
{
  SVector2 deltapar = pcp.par2 - pcp.par1;
  deltapar/=std::sqrt( deltapar[0]*deltapar[0] + deltapar[1]*deltapar[1] );
  SMatrixSym22 covsum = (pcp.cov2+pcp.cov1);
  return std::sqrt(ROOT::Math::Similarity(deltapar,covsum));
}

double trkf::Geometric3DVertexFitter::ipErr(const VertexWrapper& vtx, const recob::Track& tk) const
{
  return ipErr(getParsCovsOnPlane(vtx, tk));
}

double trkf::Geometric3DVertexFitter::sip(const trkf::Geometric3DVertexFitter::ParsCovsOnPlane& pcp) const
{
  return ip(pcp)/ipErr(pcp);
}

double trkf::Geometric3DVertexFitter::sip(const VertexWrapper& vtx, const recob::Track& tk) const
{
  return sip(getParsCovsOnPlane(vtx, tk));
}

double trkf::Geometric3DVertexFitter::pDist(const VertexWrapper& vtx, const recob::Track& tk) const
{
  return tk.Trajectory().StartDirection().Dot(vtx.position()-tk.Trajectory().Start());
}

trkf::VertexWrapper trkf::Geometric3DVertexFitter::unbiasedVertex(const trkf::VertexWrapper& vtx, const recob::Track& tk) const 
{
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracksSize()) {
    return vtx;
  } else {
    auto tks = vtx.tracksWithoutElement(ittoerase);
    return fitTracks(tks);
  }
}

double trkf::Geometric3DVertexFitter::chi2Unbiased(const trkf::VertexWrapper& vtx, const recob::Track& tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracksSize()) {
    return chi2(vtx,tk);
  } else {
    auto tks = vtx.tracksWithoutElement(ittoerase);
    return chi2(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::ipUnbiased(const trkf::VertexWrapper& vtx, const recob::Track& tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracksSize()) {
    return ip(vtx,tk);
  } else {
    auto tks = vtx.tracksWithoutElement(ittoerase);
    return ip(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::ipErrUnbiased(const trkf::VertexWrapper& vtx, const recob::Track& tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracksSize()) {
    return ipErr(vtx,tk);
  } else {
    auto tks = vtx.tracksWithoutElement(ittoerase);
    return ipErr(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::sipUnbiased(const trkf::VertexWrapper& vtx, const recob::Track& tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracksSize()) {
    return sip(vtx,tk);
  } else {
    auto tks = vtx.tracksWithoutElement(ittoerase);
    return sip(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::pDistUnbiased(const trkf::VertexWrapper& vtx, const recob::Track& tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracksSize()) {
    return pDist(vtx,tk);
  } else {
    auto tks = vtx.tracksWithoutElement(ittoerase);
    return pDist(fitTracks(tks),tk);
  }
}

std::vector<recob::VertexAssnMeta> trkf::Geometric3DVertexFitter::computeMeta(const VertexWrapper& vtx)
{
  return computeMeta(vtx, vtx.tracks());
}

std::vector<recob::VertexAssnMeta> trkf::Geometric3DVertexFitter::computeMeta(const VertexWrapper& vtx, const std::vector< art::Ptr<recob::Track> >& arttracks)
{
  TrackRefVec tracks;
  for (auto t : arttracks) tracks.push_back(*t);
  return computeMeta(vtx, tracks);
}

std::vector<recob::VertexAssnMeta> trkf::Geometric3DVertexFitter::computeMeta(const VertexWrapper& vtx, const TrackRefVec& trks)
{
  std::vector<recob::VertexAssnMeta> result;
  for (auto tk : trks) {
    float d = util::kBogusF;
    float i = util::kBogusF;
    float e = util::kBogusF;
    float c = util::kBogusF;
    auto ittoerase = vtx.findTrack(tk);
    if (debugLevel>1) std::cout << "computeMeta for vertex with ntracks=" << vtx.tracksSize() << std::endl;
    auto ubvtx = unbiasedVertex(vtx,tk.get());
    if (debugLevel>1) std::cout << "got unbiased vertex with ntracks=" << ubvtx.tracksSize() << " isValid=" << ubvtx.isValid() << std::endl;
    if (ubvtx.isValid()) {
      d = pDist(ubvtx, tk.get());
      auto pcop = getParsCovsOnPlane(ubvtx, tk.get());
      i = ip(pcop);
      e = ipErr(pcop);
      c = chi2(pcop);
      if (debugLevel>1) std::cout << "unbiasedVertex d=" << d << " i=" << i << " e=" << e << " c=" << c << std::endl;
    } else if (vtx.tracksSize()==2 && ittoerase != vtx.tracksSize()) {
      auto tks = vtx.tracksWithoutElement(ittoerase);
      auto fakevtx = closestPointAlongTrack(tks[0],tk);
      d = pDist(fakevtx, tk.get());
      // these will be identical for the two tracks (modulo numerical instabilities in the matrix inversion for the chi2)
      auto pcop = getParsCovsOnPlane(fakevtx, tk.get());
      i = ip(pcop);
      e = ipErr(pcop);
      c = chi2(pcop);
      if (debugLevel>1) std::cout << "closestPointAlongTrack d=" << d << " i=" << i << " e=" << e << " c=" << c << std::endl;
    }
    if (ittoerase == vtx.tracksSize()) {
      result.push_back(recob::VertexAssnMeta(d,i,e,c,recob::VertexAssnMeta::NotUsedInFit));
    } else {
      result.push_back(recob::VertexAssnMeta(d,i,e,c,recob::VertexAssnMeta::IncludedInFit));
    }
  }
  return result;
}
