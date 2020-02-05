#include "larreco/RecoAlg/TRACSCheatingAlg.h"

shower::TRACSCheatingAlg::TRACSCheatingAlg(const fhicl::ParameterSet& pset):
  fTRACSAlg(pset.get<fhicl::ParameterSet>("TRACSAlg"))
{
  fPFParticleModuleLabel = pset.get<art::InputTag> ("PFParticleModuleLabel");
  fHitModuleLabel        = pset.get<art::InputTag> ("HitModuleLabel");
  fShowerStartPositionInputLabel = pset.get<std::string>("ShowerStartPositionInputFile");
  fShowerDirectionInputLabel     = pset.get<std::string>("ShowerDirectionInputFile");
  fInitialTrackSpacePointsInputLabel = pset.get<std::string>("InitialTrackSpacePointsInputLabel");
}

std::map<int,const simb::MCParticle*>  shower::TRACSCheatingAlg::GetTrueParticleMap() const {

  const sim::ParticleList& particles = particleInventory->ParticleList();

  std::map<int,const simb::MCParticle*> trueParticles;
  // Make a map of track id to particle
  for (sim::ParticleList::const_iterator particleIt = particles.begin();
      particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    trueParticles[particle->TrackId()] = particle;
    //std::cout<<"Particle ID: "<<particle->TrackId()<<" and PDG: "<<particle->PdgCode()<<std::endl;
  }
  return trueParticles;
}


std::map<int,std::vector<int> > shower::TRACSCheatingAlg::GetTrueChain(
    std::map<int,const simb::MCParticle*>& trueParticles) const {

  // Roll up showers if not already done:
  std::map<int,std::vector<int> > showerMothers;

  // Loop over daughters and find th`e mothers
  for (const auto &particleIt : trueParticles ){
    const simb::MCParticle *particle = particleIt.second;
    const simb::MCParticle *mother   = particle;
    // While the grand mother exists and is an electron or photon
    // Note the true mother will skip this loop and fill itself into the map
    while (mother->Mother()!=0 && trueParticles.find(mother->Mother())!= trueParticles.end() &&
        (TMath::Abs(trueParticles[mother->Mother()]->PdgCode()==11) ||
         TMath::Abs(trueParticles[mother->Mother()]->PdgCode()==22))){
      mother = trueParticles[mother->Mother()];
    }
    //std::cout<<"Mother Done"<<std::endl;
    // This if statement double checks mothers and discards non-shower primary particles
    showerMothers[mother->TrackId()].push_back(particle->TrackId());
  }
  return showerMothers;
}

void shower::TRACSCheatingAlg::CheatDebugEVD(const simb::MCParticle* trueParticle,
    art::Event const& Event, reco::shower::ShowerElementHolder& ShowerEleHolder,
    const art::Ptr<recob::PFParticle>& pfparticle) const {

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

  // Initialise variables
  float x;
  float y;
  float z;

  std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints;
  std::vector<art::Ptr<recob::SpacePoint> > otherSpacePoints;

  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(Event.getByLabel(fHitModuleLabel, hitHandle)){
    art::fill_ptr_vector(hits, hitHandle);
  }


  // Get the hits associated with the space points
  art::FindManyP<recob::SpacePoint> fmsph(hitHandle, Event, fPFParticleModuleLabel);
  if(!fmsph.isValid()){
    throw cet::exception("ShowerTrackFinderCheater") << "Spacepoint and hit association not valid. Stopping.";
    return;
  }

  std::map< art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit> > spacePointHitMap;
  //Get the hits from the true particle
  for (auto hit : hits){
    int trueParticleID = TMath::Abs(TrueParticleID(hit));
    std::vector<art::Ptr<recob::SpacePoint> > sps = fmsph.at(hit.key());
    if (sps.size() == 1){
      art::Ptr<recob::SpacePoint> sp = sps.front();
      if (trueParticleID == trueParticle->TrackId()){
        showerSpacePoints.push_back(sp);
      } else {
        otherSpacePoints.push_back(sp);
      }
    }
  }


  if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
    mf::LogError("Shower3DTrackFinder") << "Start position not set, returning "<< std::endl;
    return;
  }
  if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
    mf::LogError("Shower3DTrackFinder") << "Direction not set, returning "<< std::endl;
    return;
  }
  if(!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)){
    mf::LogError("Shower3DTrackFinder") << "TrackSpacePoints not set, returning "<< std::endl;
    return;
  }

  // Get info from shower property holder
  TVector3 showerStartPosition = {-999,-999,-999};
  TVector3 showerDirection = {-999,-999,-999};
  std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

  ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,showerStartPosition);
  ShowerEleHolder.GetElement(fShowerDirectionInputLabel, showerDirection);
  ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel,trackSpacePoints);

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
  TPolyMarker3D* showerPoly = new TPolyMarker3D(showerSpacePoints.size());
  for (auto spacePoint : showerSpacePoints){
    TVector3 pos = fTRACSAlg.SpacePointPosition(spacePoint) - showerStartPosition;

    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    showerPoly->SetPoint(point,x,y,z);
    ++point;

    // Calculate the projection of (point-startpoint) along the direction
    double proj = fTRACSAlg.SpacePointProjection(spacePoint, showerStartPosition,
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
    TVector3 pos = fTRACSAlg.SpacePointPosition(spacePoint) - showerStartPosition;
    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    trackPoly->SetPoint(point,x,y,z);
    ++point;
  } // loop over track spacepoints

  //  we want to draw all of the PFParticles in the event
  //Get the PFParticles

  TPolyMarker3D* otherPoly = new TPolyMarker3D(otherSpacePoints.size());

  // initialise counters
  point = 0; // re-initialise counter

  for(auto const& sp: otherSpacePoints){
    TVector3 pos = fTRACSAlg.SpacePointPosition(sp) - showerStartPosition;
    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    otherPoly->SetPoint(point,x,y,z);
    ++point;
  }

  otherPoly->SetMarkerStyle(20);
  otherPoly->SetMarkerColor(4);
  otherPoly->Draw();


  // Draw all of the things
  showerPoly->SetMarkerStyle(20);
  showerPoly->Draw();
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

int shower::TRACSCheatingAlg::TrueParticleID(const art::Ptr<recob::Hit>& hit) const {
  double particleEnergy = 0;
  int likelyTrackID = 0;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = trackIDs.at(idIt).trackID;
    }
  }
  return likelyTrackID;
}

std::pair<int,double> shower::TRACSCheatingAlg::TrueParticleIDFromTrueChain(std::map<int,std::vector<int>> const& ShowersMothers,
                                                                            std::vector<art::Ptr<recob::Hit> > const& hits, int planeid) const {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //Find the energy for each track ID.
  std::map<int,double> trackIDToEDepMap;
  std::map<int,double> trackIDTo3EDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;

    //Get the plane ID
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDTo3EDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
      if(PlaneID == planeid){trackIDToEDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;}
    }
  }

  //Find the energy for each showermother.
  std::map<int,double> MotherIDtoEMap;
  std::map<int,double> MotherIDto3EMap;
  for(std::map<int,std::vector<int> >::const_iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
    for(std::vector<int>::const_iterator daughter=(showermother->second).begin(); daughter!=(showermother->second).end(); ++daughter){
      MotherIDtoEMap[showermother->first]  +=  trackIDToEDepMap[*daughter];
      MotherIDto3EMap[showermother->first] +=  trackIDTo3EDepMap[*daughter];
    }
  }

  //Loop over the mothers to find the most like candidate by identifying which true shower deposited the most energy in the hits.
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = MotherIDto3EMap.begin(); mapIt != MotherIDto3EMap.end(); mapIt++){
    double energy_three = mapIt->second;
    double trackid = mapIt->first;
    if (energy_three > maxenergy){
      maxenergy = energy_three;
      objectTrack = trackid;
    }
  }

  //If the none of the shower mother deposited no energy then we cannot match this.
  if(maxenergy == 0){
    return std::make_pair(-99999,-99999);
  }

  return std::make_pair(objectTrack,MotherIDtoEMap[objectTrack]);
}
