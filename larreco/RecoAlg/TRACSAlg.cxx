#include "larreco/RecoAlg/TRACSAlg.h"

shower::TRACSAlg::TRACSAlg(const fhicl::ParameterSet& pset):
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
  fUseCollectionOnly     = pset.get<bool>("UseCollectionOnly");
  fPFParticleModuleLabel = pset.get<art::InputTag> ("PFParticleModuleLabel");
  fHitModuleLabel        = pset.get<art::InputTag> ("HitModuleLabel");
}


//Order the shower hits with regards to their projected length onto
//the shower direction from the shower start position. This is done
//in the 2D coordinate system (wire direction, x)
void shower::TRACSAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> >& hits,
    TVector3 const& ShowerStartPosition,
    TVector3 const& ShowerDirection
    ) const {

  std::map<double, art::Ptr<recob::Hit> > OrderedHits;
  art::Ptr<recob::Hit> startHit = hits.front();

  //Get the wireID
  const geo::WireID startWireID = startHit->WireID();

  //Get the plane
  const geo::PlaneID planeid = startWireID.asPlaneID();

  //Get the pitch
  double pitch = fGeom->WirePitch(planeid);

  TVector2 Shower2DStartPosition = {
    fGeom->WireCoordinate(ShowerStartPosition, startHit->WireID().planeID())*pitch,
    ShowerStartPosition.X()
  };

  //Vector of the plane
  TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection();

  //get the shower 2D direction
  TVector2 Shower2DDirection = {
    ShowerDirection.Dot(PlaneDirection),
    ShowerDirection.X()
  };


  Shower2DDirection = Shower2DDirection.Unit();

  for(auto const& hit: hits){

    //Get the wireID
    const geo::WireID WireID = hit->WireID();

    if (WireID.asPlaneID() != startWireID.asPlaneID()) {
      break;
    }

    //Get the hit Vector.
    TVector2 hitcoord = HitCoordinates(hit);

    //Order the hits based on the projection
    TVector2 pos = hitcoord - Shower2DStartPosition;
    OrderedHits[pos*Shower2DDirection] = hit;
  }

  //Transform the shower.
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(OrderedHits.begin(), OrderedHits.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

  //Sometimes get the order wrong. Depends on direction compared to the plane Correct for it here:
  art::Ptr<recob::Hit> frontHit = showerHits.front();
  art::Ptr<recob::Hit> backHit  = showerHits.back();

  //Get the hit Vector.
  TVector2 fronthitcoord = HitCoordinates(frontHit);
  TVector2 frontpos = fronthitcoord - Shower2DStartPosition;


  //Get the hit Vector.
  TVector2 backhitcoord  = HitCoordinates(backHit);
  TVector2 backpos = backhitcoord - Shower2DStartPosition;

  double frontproj = frontpos*Shower2DDirection;
  double backproj  = backpos*Shower2DDirection;
  if (TMath::Abs(backproj) < TMath::Abs(frontproj)){
    std::reverse(showerHits.begin(),showerHits.end());
  }

  hits = showerHits;
  return;
}


//Orders the shower spacepoints with regards to there prejected length from
//the shower start position in the shower direction.
void shower::TRACSAlg::OrderShowerSpacePoints( std::vector<art::Ptr<recob::SpacePoint> >&
    showersps, TVector3 const& vertex,
    TVector3 const& direction) const {

  std::map<double,art::Ptr<recob::SpacePoint> > OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.
  for(auto const& sp: showersps){

    // Get the projection of the space point along the direction
    double len = SpacePointProjection(sp, vertex, direction);

    //Add to the list
    OrderedSpacePoints[len] = sp;
  }

  //Return an ordered list.
  showersps.clear();
  for(auto const& sp: OrderedSpacePoints){
    showersps.push_back(sp.second);
  }
  return;
}

TVector3 shower::TRACSAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> > const&
    showerspcs, art::FindManyP<recob::Hit> const& fmh) const {
  float totalCharge=0;
  TVector3 centre =  shower::TRACSAlg::ShowerCentre(showerspcs,fmh,totalCharge);
  return centre;

}

//Returns the vector to the shower centre and the total charge of the shower.
TVector3 shower::TRACSAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> > const& showersps,
    art::FindManyP<recob::Hit> const& fmh, float& totalCharge) const {

  TVector3 pos, chargePoint = TVector3(0,0,0);

  //Loop over the spacepoints and get the charge weighted center.
  for(auto const& sp: showersps){

    //Get the position of the spacepoint
    pos = SpacePointPosition(sp);

    //Get the associated hits
    std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());

    //Average the charge unless sepcified.
    float charge  = 0;
    float charge2 = 0;
    for(auto const& hit: hits){

      if(fUseCollectionOnly){
        if(hit->SignalType() == geo::kCollection){
          charge = hit->Integral();
          //Correct for the lifetime: Need to do other detproperites
          charge *= TMath::Exp((fDetProp->SamplingRate() * (hit)->PeakTime()) / (fDetProp->ElectronLifetime()*1e3));
          break;
        }
      }

      //Check if any of the points are not withing 2 sigma.
      if(!fUseCollectionOnly){

        //Correct for the lifetime FIX: Need  to do other detproperties somehow
        double Q = hit->Integral()*TMath::Exp( (fDetProp->SamplingRate() * (hit)->PeakTime()) / (fDetProp->ElectronLifetime()*1e3));

        charge  += Q;
        charge2 += Q*Q;
      }
    }

    if(!fUseCollectionOnly){
      //Calculate the unbiased standard deviation and mean.
      float mean = charge/((float) hits.size());

      float rms = 1;

      if(hits.size() > 1){
        rms  = TMath::Sqrt((charge2 - charge*charge)/((float)(hits.size()-1)));
      }

      charge = 0;
      int n = 0;
      for(auto const& hit: hits){
        double lifetimecorrection =  TMath::Exp((fDetProp->SamplingRate() * (hit)->PeakTime()) / (fDetProp->ElectronLifetime()*1e3));
        if(hit->Integral()*lifetimecorrection > (mean - 2*rms) && hit->Integral()*lifetimecorrection < (mean + 2*rms)){
          charge += hit->Integral()*lifetimecorrection;
          ++n;
        }
      }

      if(n==0){
        mf::LogWarning("TRACSAlg") <<
          "no points used to make the charge value. \n";
      }

      charge /= n;
    }

    chargePoint += charge * pos;
    totalCharge += charge;

    if(charge == 0){
      mf::LogWarning("TRACSAlg") <<
        "Averaged charge, within 2 sigma, for a spacepoint is zero, Maybe this not a good method. \n";
    }
  }

  double intotalcharge = 1/totalCharge;
  TVector3 centre = chargePoint *  intotalcharge;
  return centre;

}

//Return the spacepoint position in 3D cartesian coordinates.
TVector3 shower::TRACSAlg::SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const {

  const Double32_t* sp_xyz = sp->XYZ();
  TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
  return sp_postiion;
}


//Return the charge of the spacepoint in ADC.
double shower::TRACSAlg::SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
    art::FindManyP<recob::Hit> const& fmh) const {

  double Charge = 0;

  //Average over the charge even though there is only one
  std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
  for(auto const& hit: hits){
    Charge += hit->Integral();
  }

  Charge /= (float) hits.size();

  return Charge;
}

//Return the spacepoint time.
double shower::TRACSAlg::SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
    art::FindManyP<recob::Hit> const& fmh) const {

  double Time = 0;

  //Avergae over the hits
  std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
  for(auto const& hit: hits){
    Time += hit->PeakTime();
  }

  Time /= (float) hits.size();
  return Time;
}


//Return the cooordinates of the hit in cm in wire direction and x.
TVector2 shower::TRACSAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) const {

  //Get the pitch
  const geo::WireID  WireID = hit->WireID();
  const geo::PlaneID planeid = WireID.asPlaneID();
  double pitch = fGeom->WirePitch(planeid);

  return TVector2(WireID.Wire*pitch ,fDetProp->ConvertTicksToX(hit->PeakTime(),planeid));
}

double shower::TRACSAlg::SpacePointProjection(const art::Ptr<recob::SpacePoint>&sp,
    TVector3 const& vertex, TVector3 const& direction) const {

  // Get the position of the spacepoint
  TVector3 pos = shower::TRACSAlg::SpacePointPosition(sp) - vertex;

  // Get the the projected length
  double projLen = pos.Dot(direction);

  return projLen;
}

double shower::TRACSAlg::SpacePointPerpendiular(art::Ptr<recob::SpacePoint> const &sp,
    TVector3 const& vertex, TVector3 const& direction,
    double proj) const {

  // Get the position of the spacepoint
  TVector3 pos = shower::TRACSAlg::SpacePointPosition(sp) - vertex;

  // Take away the projection * distance to find the perpendicular vector

  pos = pos - proj * direction;

  // Get the the projected length
  double perpLen = pos.Mag();

  return perpLen;
}


void shower::TRACSAlg::DebugEVD(art::Ptr<recob::PFParticle> const& pfparticle,
    art::Event const& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder) const {

  std::cout<<"Making Debug Event Display"<<std::endl;

  //Function for drawing reco showers to check direction and initial track selection

  // Get run info to make unique canvas names
  int run    = Event.run();
  int subRun = Event.subRun();
  int event  = Event.event();
  int PFPID  = pfparticle->Self();

  // Create the canvas
  TString canvasName = Form("canvas_%i_%i_%i_%i",run,subRun,event,PFPID);
  TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);

  std::cout << "canvasName: " << canvasName << std::endl;

  // Initialise variables
  float x;
  float y;
  float z;

  // Get a bunch of associations (again)
  // N.B. this is a horribly inefficient way of doing things but as this is only
  // going to be used to debug I don't care, I would rather have generality in this case

  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
    throw cet::exception("Shower3DTrackFinderEMShower") << "Could not get the pandora pf particles\
      . Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
    return;
  }

  // Get the spacepoint - PFParticle assn
  art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
  if (!fmspp.isValid()){
    throw cet::exception("Shower3DTrackFinder") << "Trying to get the spacepoint and failed. Somet\
      hing is not configured correctly. Stopping ";
    return;
  }

  // Get the SpacePoints
  std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

  //We cannot progress with no spacepoints.
  if(spacePoints.size() == 0){
    //throw cet::exception("Shower3DTrackFinder") << "No Space Points. Stopping.";
    return;
  }

  if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
    mf::LogError("Shower3DTrackFinder") << "Start position not set, returning "<< std::endl;
    return;
  }
  if(!ShowerEleHolder.CheckElement("ShowerDirection")){
    mf::LogError("Shower3DTrackFinder") << "Direction not set, returning "<< std::endl;
    return;
  }
  if(!ShowerEleHolder.CheckElement("InitialTrackSpacePoints")){
    mf::LogError("Shower3DTrackFinder") << "TrackSpacePoints not set, returning "<< std::endl;
    return;
  }

  // Get info from shower property holder
  TVector3 showerStartPosition = {-999,-999,-999};
  TVector3 showerDirection = {-999,-999,-999};
  std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

  ShowerEleHolder.GetElement("ShowerStartPosition",showerStartPosition);
  ShowerEleHolder.GetElement("ShowerDirection", showerDirection);
  ShowerEleHolder.GetElement("InitialTrackSpacePoints",trackSpacePoints);

  // Create 3D point at vertex, chosed to be origin for ease of use of display
  double startXYZ[3] = {0,0,0};
  TPolyMarker3D* startPoly = new TPolyMarker3D(1,startXYZ);

  // Get the min and max projections along the direction to know how long to draw
  // the direction line
  double minProj=-1;
  double maxProj=1;

  //initialise counter point
  int point = 0;
  // Make 3D points for each spacepoint in the shower
  TPolyMarker3D* allPoly = new TPolyMarker3D(spacePoints.size());
  for (auto spacePoint : spacePoints){
    TVector3 pos = shower::TRACSAlg::SpacePointPosition(spacePoint) - showerStartPosition;

    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    allPoly->SetPoint(point,x,y,z);
    ++point;

    // Calculate the projection of (point-startpoint) along the direction
    double proj = shower::TRACSAlg::SpacePointProjection(spacePoint, showerStartPosition,
        showerDirection);
    if (proj>maxProj) {
      maxProj = proj;
    } else if (proj<minProj) {
      minProj = proj ;
    }

  } // loop over spacepoints

  // Create TPolyLine3D arrays
  double xDirPoints[3] = {minProj*showerDirection.X(), 0, maxProj*showerDirection.X()};
  double yDirPoints[3] = {minProj*showerDirection.Y(), 0, maxProj*showerDirection.Y()};
  double zDirPoints[3] = {minProj*showerDirection.Z(), 0, maxProj*showerDirection.Z()};

  TPolyLine3D* dirPoly = new TPolyLine3D(3,xDirPoints,yDirPoints,zDirPoints);

  point = 0; // re-initialise counter
  TPolyMarker3D* trackPoly = new TPolyMarker3D(trackSpacePoints.size());
  for (auto spacePoint : trackSpacePoints){
    TVector3 pos = shower::TRACSAlg::SpacePointPosition(spacePoint) - showerStartPosition;
    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    trackPoly->SetPoint(point,x,y,z);
    ++point;
  } // loop over track spacepoints

  //  we want to draw all of the PFParticles in the event
  //Get the PFParticles
  std::vector<art::Ptr<recob::PFParticle> > pfps;
  if (Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
    art::fill_ptr_vector(pfps, pfpHandle);
  }
  else {
    throw cet::exception("TRACS") << "pfps not loaded." << std::endl;
  }
  // initialse counters
  // Split into tracks and showers to make it clearer what pandora is doing
  int pfpTrackCounter = 0;
  int pfpShowerCounter = 0;

  // initial loop over pfps to find nuber of spacepoints for tracks and showers
  for(auto const& pfp: pfps){
    std::vector<art::Ptr<recob::SpacePoint> > sps = fmspp.at(pfp.key());
    if (pfp->PdgCode()==11) { pfpShowerCounter += sps.size();
    } else if (pfp->PdgCode()==13) { pfpTrackCounter += sps.size(); }
  }

  TPolyMarker3D* pfpPolyTrack = new TPolyMarker3D(pfpTrackCounter);
  TPolyMarker3D* pfpPolyShower = new TPolyMarker3D(pfpShowerCounter);

  // initialise counters
  int trackPoints  = 0;
  int showerPoints = 0;

  for(auto const& pfp: pfps){
    std::vector<art::Ptr<recob::SpacePoint> > sps = fmspp.at(pfp.key());
    int pdg = pfp->PdgCode(); // Track or shower
    for (auto sp : sps){
      TVector3 pos = shower::TRACSAlg::SpacePointPosition(sp) - showerStartPosition;
      x = pos.X();
      y = pos.Y();
      z = pos.Z();
      if (pdg==11){
        pfpPolyShower->SetPoint(showerPoints,x,y,z);
        ++showerPoints;
      } else if (pdg==13){
        pfpPolyTrack->SetPoint(trackPoints,x,y,z);
        ++trackPoints;
      }
    } // loop over sps

    pfpPolyShower->SetMarkerStyle(20);
    pfpPolyShower->SetMarkerColor(4);
    pfpPolyShower->Draw();
    pfpPolyTrack->SetMarkerStyle(20);
    pfpPolyTrack->SetMarkerColor(6);
    pfpPolyTrack->Draw();

  } // if (fDrawAllPFPs)

  // Draw all of the things
  allPoly->SetMarkerStyle(20);
  allPoly->Draw();
  trackPoly->SetMarkerStyle(20);
  trackPoly->SetMarkerColor(2);
  trackPoly->Draw();
  startPoly->SetMarkerStyle(21);
  startPoly->SetMarkerSize(2);
  startPoly->SetMarkerColor(3);
  startPoly->Draw();
  dirPoly->SetLineWidth(1);
  dirPoly->SetLineColor(6);
  dirPoly->Draw();

  // Save the canvas. Don't usually need this when using TFileService but this in the alg
  // not a module and didn't work without this so im going with it.
  canvas->Write();
}
