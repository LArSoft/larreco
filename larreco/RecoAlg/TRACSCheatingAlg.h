#ifndef TRACSCheatingAlg_hxx
#define TRACSCheatingAlg_hxx

//Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/ShowerElementHolder.hh"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larreco/RecoAlg/TRACSAlg.h"

//C++ Includes
#include <iostream>
#include <vector>
#include <map>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TString.h"

namespace shower {
  class TRACSCheatingAlg;
}

class shower::TRACSCheatingAlg {
  public:
    TRACSCheatingAlg(const fhicl::ParameterSet& pset);

    std::map<int,const simb::MCParticle*> GetTrueParticleMap() const;
    std::map<int,std::vector<int> > GetTrueChain(std::map<int,const simb::MCParticle*>& trueParticles) const;
    void CheatDebugEVD(const simb::MCParticle* trueParticle, art::Event const& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder,
        const art::Ptr<recob::PFParticle>& pfparticle) const;

    int TrueParticleID(const art::Ptr<recob::Hit>& hit) const;

    std::pair<int,double> TrueParticleIDFromTrueChain(std::map<int,std::vector<int> > const& ShowersMothers,
                                                      std::vector<art::Ptr<recob::Hit> > const& hits, int planeid) const;

  private:

    shower::TRACSAlg fTRACSAlg;

    art::InputTag                                       fHitModuleLabel;
    art::InputTag                                       fPFParticleModuleLabel;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<art::TFileService>   tfs;


    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fInitialTrackSpacePointsInputLabel;

};
#endif
