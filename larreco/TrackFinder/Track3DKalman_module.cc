////////////////////////////////////////////////////////////////////////
//
// \file Track3DKalman.cxx
//
// \author echurch@fnal.gov
//
//  This algorithm is designed to reconstruct 3D tracks through
//  GENFIT's Kalman filter.
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <algorithm>
#include <string>
#include <vector>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nurandom
#include "nurandom/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

// ROOT includes
#include "TMatrixT.h"
#include "TTree.h"
#include "TVector3.h"

// GENFIT includes
#include "larreco/Genfit/GFAbsTrackRep.h"
#include "larreco/Genfit/GFConstField.h"
#include "larreco/Genfit/GFException.h"
#include "larreco/Genfit/GFFieldManager.h"
#include "larreco/Genfit/GFKalman.h"
#include "larreco/Genfit/GFTrack.h"
#include "larreco/Genfit/PointHit.h"
#include "larreco/Genfit/RKTrackRep.h"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace trkf {

  class Track3DKalman : public art::EDProducer {
  public:
    explicit Track3DKalman(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;
    void beginJob() override;
    void endJob() override;

    std::string fSpacePtsModuleLabel; // label for input collection
    std::string fGenieGenModuleLabel; // label for input MC single particle generator
    std::string fG4ModuleLabel;       // label for input MC single particle generator
    bool fGenfPRINT;

    TTree* tree;

    TMatrixT<Double_t>* stMCT;
    TMatrixT<Double_t>* covMCT;
    TMatrixT<Double_t>* stREC;
    TMatrixT<Double_t>* covREC;
    Float_t chi2;
    Float_t chi2ndf;

    Float_t* fpRECt3D;
    Float_t* fpRECL;
    Float_t* fpREC;
    Float_t* fpMCT;
    int nfail;
    int ndf;
    unsigned int evtt;
    unsigned int nTrks;
    unsigned int fptsNo;
    Float_t* fshx;
    Float_t* fshy;
    Float_t* fshz;
    unsigned int fDimSize; // if necessary will get this from pset in constructor.

    std::vector<double> fPosErr;
    std::vector<double> fMomErr;
    std::vector<double> fMomStart;
    genf::GFAbsTrackRep* repMC;
    genf::GFAbsTrackRep* rep;
    CLHEP::HepRandomEngine& fEngine;
  }; // class Track3DKalman

} // end namespace

namespace {
  bool sp_sort_3dz(const art::Ptr<recob::SpacePoint>& h1, const art::Ptr<recob::SpacePoint>& h2)
  {
    const double* xyz1 = h1->XYZ();
    const double* xyz2 = h2->XYZ();
    return xyz1[2] < xyz2[2];
  }
}

namespace trkf {

  //-------------------------------------------------
  Track3DKalman::Track3DKalman(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fSpacePtsModuleLabel{pset.get<std::string>("SpacePtsModuleLabel")}
    , fGenieGenModuleLabel{pset.get<std::string>("GenieGenModuleLabel")}
    , fG4ModuleLabel{pset.get<std::string>("G4ModuleLabel")}
    , fGenfPRINT{pset.get<bool>("GenfPRINT")}
    , fPosErr{pset.get<std::vector<double>>("PosErr3")}     // resolution. cm
    , fMomErr{pset.get<std::vector<double>>("MomErr3")}     // GeV
    , fMomStart{pset.get<std::vector<double>>("MomStart3")} // Will be unit norm'd.
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0),
                                                                                 pset,
                                                                                 "Seed"))
  {
    produces<std::vector<recob::Track>>();
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Track, recob::Cluster>>();
    produces<art::Assns<recob::Track, recob::SpacePoint>>();
    produces<art::Assns<recob::Track, recob::Hit>>();
  }

  //-------------------------------------------------
  void Track3DKalman::beginJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;

    stMCT = new TMatrixT<Double_t>(5, 1);
    covMCT = new TMatrixT<Double_t>(5, 5);
    stREC = new TMatrixT<Double_t>(5, 1);
    covREC = new TMatrixT<Double_t>(5, 5);

    fpMCT = new Float_t[4];
    fpREC = new Float_t[4];
    fpRECL = new Float_t[4];
    fpRECt3D = new Float_t[4];
    fDimSize = 400; // if necessary will get this from pset in constructor.

    fshx = new Float_t[fDimSize];
    fshy = new Float_t[fDimSize];
    fshz = new Float_t[fDimSize];

    tree = tfs->make<TTree>("GENFITttree", "GENFITttree");
    //tree->Branch("stMCT",&stMCT,"stMCT[5]/F"); // "TMatrixT<Double_t>"

    tree->Branch("stMCT", "TMatrixD", &stMCT, 64000, 0);
    //tree->Branch("covMCT",&covMCT,"covMCT[25]/F");
    tree->Branch("covMCT", "TMatrixD", &covMCT, 64000, 0);
    //tree->Branch("stREC",&stREC,"stREC[5]/F");
    tree->Branch("stREC", "TMatrixD", &stREC, 64000, 0);
    //tree->Branch("covREC",&covREC,"covREC[25]/F");
    tree->Branch("covREC", "TMatrixD", &covREC, 64000, 0);

    tree->Branch("chi2", &chi2, "chi2/F");
    tree->Branch("nfail", &nfail, "nfail/I");
    tree->Branch("ndf", &ndf, "ndf/I");
    tree->Branch("evtNo", &evtt, "evtNo/I");
    tree->Branch("chi2ndf", &chi2ndf, "chi2ndf/F");

    tree->Branch("trkNo", &nTrks, "trkNo/I");
    tree->Branch("ptsNo", &fptsNo, "ptsNo/I");
    tree->Branch("shx", fshx, "shx[ptsNo]/F");
    tree->Branch("shy", fshy, "shy[ptsNo]/F");
    tree->Branch("shz", fshz, "shz[ptsNo]/F");

    tree->Branch("pMCT", fpMCT, "pMCT[4]/F");
    tree->Branch("pRECKalF", fpREC, "pRECKalF[4]/F");
    tree->Branch("pRECKalL", fpRECL, "pRECKalL[4]/F");
    tree->Branch("pRECt3D", fpRECt3D, "pRECt3D[4]/F");
  }

  //-------------------------------------------------
  void Track3DKalman::endJob()
  {
    if (!rep) delete rep;
    if (!repMC) delete repMC;
  }

  //------------------------------------------------------------------------------------//
  void Track3DKalman::produce(art::Event& evt)
  {
    rep = 0;
    repMC = 0;

    // get services
    art::ServiceHandle<geo::Geometry const> geom;
    CLHEP::RandGaussQ gauss(fEngine);

    //////////////////////////////////////////////////////
    // Make a std::unique_ptr<> for the thing you want to put into the event
    // because that handles the memory management for you
    //////////////////////////////////////////////////////
    std::unique_ptr<std::vector<recob::Track>> tcol(new std::vector<recob::Track>);
    std::unique_ptr<std::vector<recob::SpacePoint>> spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint>> tspassn(
      new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::Cluster>> tcassn(
      new art::Assns<recob::Track, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>> thassn(
      new art::Assns<recob::Track, recob::Hit>);

    // define TPC parameters
    TString tpcName = geom->GetLArTPCVolumeName();

    // get input Hit object(s).
    art::Handle<std::vector<recob::Track>> trackListHandle;
    evt.getByLabel(fSpacePtsModuleLabel, trackListHandle);

    art::PtrVector<simb::MCTruth> mclist;

    if (!evt.isRealData()) {
      art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
      evt.getByLabel(fGenieGenModuleLabel, mctruthListHandle);

      for (unsigned int ii = 0; ii < mctruthListHandle->size(); ++ii) {
        art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle, ii);
        mclist.push_back(mctparticle);
      }
    }

    //create collection of spacepoints that will be used when creating the Track object
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints;
    art::PtrVector<recob::Track> trackIn;
    // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
    mf::LogInfo("Track3DKalman: ")
      << "There are " << trackListHandle->size()
      << " Track3Dreco/SpacePt tracks/groups (whichever) in this event.";

    art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fSpacePtsModuleLabel);
    art::FindManyP<recob::Cluster> fmc(trackListHandle, evt, fSpacePtsModuleLabel);
    art::FindManyP<recob::Hit> fmh(trackListHandle, evt, fSpacePtsModuleLabel);

    for (unsigned int ii = 0; ii < trackListHandle->size(); ++ii) {
      art::Ptr<recob::Track> track(trackListHandle, ii);
      trackIn.push_back(track);
    }

    TVector3 MCOrigin;
    TVector3 MCMomentum;
    // TVector3 posErr(.05,.05,.05); // resolution. 0.5mm
    // TVector3 momErr(.1,.1,0.2);   // GeV
    TVector3 posErr(fPosErr[0], fPosErr[1], fPosErr[2]); // resolution. 0.5mm
    TVector3 momErr(fMomErr[0], fMomErr[1], fMomErr[2]); // GeV

    // This is strictly for MC
    ///\todo: remove this statement, there is no place for checks on isRealData in reconstruction code
    if (!evt.isRealData()) {
      // Below breaks are stupid, I realize. But rather than keep all the MC
      // particles I just take the first primary, e.g., muon and only keep its
      // info in the branches of the Ttree. I could generalize later, ...
      for (unsigned int ii = 0; ii < mclist.size(); ++ii) {
        //art::Ptr<const simb::MCTruth> mc(mctruthListHandle,i);
        art::Ptr<simb::MCTruth> mc(mclist[ii]);
        for (int jj = 0; jj < mc->NParticles(); ++jj) {
          simb::MCParticle part(mc->GetParticle(jj));
          mf::LogInfo("Track3DKalman: ")
            << "FROM MC TRUTH, the particle's pdg code is: " << part.PdgCode()
            << " with energy = " << part.E() << ", with energy = " << part.E()
            << " and vtx and momentum in Global (not volTPC) coords are ";
          MCOrigin.SetXYZ(part.Vx(), part.Vy(), part.Vz()); // V for Vertex
          MCMomentum.SetXYZ(part.Px(), part.Py(), part.Pz());
          MCOrigin.Print();
          MCMomentum.Print();
          repMC = new genf::RKTrackRep(MCOrigin, MCMomentum, posErr, momErr, part.PdgCode());
          break;
        }
        break;
      }
      //for saving of MC truth
      stMCT->ResizeTo(repMC->getState());
      *stMCT = repMC->getState();
      covMCT->ResizeTo(repMC->getCov());
      *covMCT = repMC->getCov();
      mf::LogInfo("Track3DKalman: ") << " repMC, covMC are ... ";
      repMC->getState().Print();
      repMC->getCov().Print();

    } // !isRealData

    art::PtrVector<recob::Track>::const_iterator trackIter = trackIn.begin();

    nTrks = 0;
    while (trackIter != trackIn.end()) {
      spacepoints.clear();
      spacepoints = fmsp.at(nTrks);

      mf::LogInfo("Track3DKalman: ") << "found " << spacepoints.size() << " 3D spacepoint(s).";

      // Add the 3D track to the vector of the reconstructed tracks
      if (spacepoints.size() > 0) {
        // Insert the GENFIT/Kalman stuff here then make the tracks. Units are cm, GeV.
        const double resolution = 0.5; // dunno, 5 mm
        const int numIT = 3;           // 3->1, EC, 6-Jan-2011. Back, 7-Jan-2011.

        //TVector3 mom(0.0,0.0,2.0);
        TVector3 mom(fMomStart[0], fMomStart[1], fMomStart[2]);
        //mom.SetMag(1.);
        TVector3 momM(mom);
        momM.SetX(gauss.fire(momM.X(), momErr.X() /* *momM.X() */));
        momM.SetY(gauss.fire(momM.Y(), momErr.Y() /* *momM.Y() */));
        momM.SetZ(gauss.fire(momM.Z(), momErr.Z() /* *momM.Z() */));
        //std::cout << "Track3DKalman: sort spacepoints by z

        std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dz);

        //std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dx); // Reverse sort!

        genf::GFFieldManager::getInstance()->init(new genf::GFConstField(0., 0., 0.0));
        genf::GFDetPlane planeG((TVector3)(spacepoints[0]->XYZ()), momM);

        //      std::cout<<"Track3DKalman about to do GAbsTrackRep."<<std::endl;
        // Initialize with 1st spacepoint location and a guess at the momentum.
        rep = new genf::RKTrackRep( //posM-.5/momM.Mag()*momM,
          (TVector3)(spacepoints[0]->XYZ()),
          momM,
          posErr,
          momErr,
          13); // mu- hypothesis
        //      std::cout<<"Track3DKalman: about to do GFTrack. repDim is " << rep->getDim() <<std::endl;

        genf::GFTrack fitTrack(rep); //initialized with smeared rep
        // Gonna sort in x cuz I want to essentially transform here to volTPC coords.
        // volTPC coords, cuz that's what the Geant3/Geane stepper wants, as that's its understanding
        // from the Geant4 geometry, which it'll use. EC, 7-Jan-2011.
        int ihit = 0;
        fptsNo = 0;
        for (unsigned int point = 0; point < spacepoints.size(); ++point) {

          TVector3 spt3(spacepoints[point]->XYZ());
          if (point > 0) {
            double eps(0.1);
            TVector3 magNew(spt3[0], spt3[1], spt3[2]);
            TVector3 magLast(spacepoints.back()->XYZ()[0],
                             spacepoints.back()->XYZ()[1],
                             spacepoints.back()->XYZ()[2]);
            if (!(magNew.Mag() >= magLast.Mag() + eps || magNew.Mag() <= magLast.Mag() - eps))
              continue;
          }

          if (point % 20) // Jump out of loop except on every 20th pt.
          {
            //continue;
            // Icarus paper suggests we may wanna decimate our data in order to give
            // trackfitter a better idea of multiple-scattering. EC, 7-Jan-2011.
            //if (std::abs(spt3[0]-spacepoints.at(point-1).XYZ()[0]) < 2) continue;
          }
          if (fptsNo < fDimSize) {
            fshx[fptsNo] = spt3[0];
            fshy[fptsNo] = spt3[1];
            fshz[fptsNo] = spt3[2];
          }
          fptsNo++;

          MF_LOG_DEBUG("Track3DKalman: ")
            << "ihit xyz..." << spt3[0] << "," << spt3[1] << "," << spt3[2];
          fitTrack.addHit(new genf::PointHit(spt3, resolution),
                          1, //dummy detector id
                          ihit++);
        }

        //      std::cout<<"Track3DKalman about to do GFKalman."<<std::endl;
        genf::GFKalman k;
        //k.setBlowUpFactor(500); // Instead of 500 out of box. EC, 6-Jan-2011.
        //k.setInitialDirection(+1); // Instead of 1 out of box. EC, 6-Jan-2011.
        k.setNumIterations(numIT);
        //      std::cout<<"Track3DKalman back from setNumIterations."<<std::endl;
        try {
          //	std::cout<<"Track3DKalman about to processTrack."<<std::endl;
          k.processTrack(&fitTrack);
          //std::cout<<"Track3DKalman back from processTrack."<<std::endl;
        }
        catch (GFException& e) {
          mf::LogError("Track3DKalman: ") << "just caught a GFException." << std::endl;
          e.what();
          mf::LogError("Track3DKalman: ")
            << "Exceptions won't be further handled ->exit(1) " << __LINE__;

          //	exit(1);
        }

        if (rep->getStatusFlag() == 0) // 0 is successful completion
        {
          MF_LOG_DEBUG("Track3DKalman: ") << __FILE__ << " " << __LINE__;
          MF_LOG_DEBUG("Track3DKalman: ") << "Track3DKalman.cxx: Original plane:";

          if (fGenfPRINT) planeG.Print();
          MF_LOG_DEBUG("Track3DKalman: ") << "Current (fit) reference Plane:";
          if (fGenfPRINT) rep->getReferencePlane().Print();

          MF_LOG_DEBUG("Track3DKalman: ") << "Track3DKalman.cxx: Last reference Plane:";
          if (fGenfPRINT) rep->getLastPlane().Print();

          if (fGenfPRINT) {
            if (planeG != rep->getReferencePlane())
              MF_LOG_DEBUG("Track3DKalman: ") << "Track3DKalman: Original hit plane (not "
                                                 "surprisingly) not current reference Plane!"
                                              << std::endl;
          }

          stREC->ResizeTo(rep->getState());
          *stREC = rep->getState();
          covREC->ResizeTo(rep->getCov());
          *covREC = rep->getCov();
          if (fGenfPRINT) {
            MF_LOG_DEBUG("Track3DKalman: ") << " Final State and Cov:";
            stREC->Print();
            covREC->Print();
          }
          chi2 = rep->getChiSqu();
          ndf = rep->getNDF();
          nfail = fitTrack.getFailedHits();
          chi2ndf = chi2 / ndf;
          TVector3 dircoss = (*trackIter)->VertexDirection<TVector3>();

          for (int ii = 0; ii < 3; ++ii) {
            fpMCT[ii] = MCMomentum[ii] / MCMomentum.Mag();
            fpREC[ii] = rep->getReferencePlane().getNormal()[ii];
            fpRECL[ii] = rep->getLastPlane().getNormal()[ii];
            fpRECt3D[ii] = dircoss[ii];
          }
          fpMCT[3] = MCMomentum.Mag();
          fpREC[3] = -1.0 / (*stREC)[0][0];

          evtt = (unsigned int)evt.id().event();
          mf::LogInfo("Track3DKalman: ")
            << "Track3DKalman about to do tree->Fill(). Chi2/ndf is " << chi2 / ndf
            << ". All in volTPC coords .... pMCT[0-3] is " << fpMCT[0] << ", " << fpMCT[1] << ", "
            << fpMCT[2] << ", " << fpMCT[3] << ". pREC[0-3] is " << fpREC[0] << ", " << fpREC[1]
            << ", " << fpREC[2] << ", " << fpREC[3] << ".";

          tree->Fill();

          // Get the clusters associated to each track in induction and collection view
          std::vector<art::Ptr<recob::Cluster>> clusters = fmc.at(nTrks);
          std::vector<art::Ptr<recob::Hit>> hits = fmh.at(nTrks);

          // Use new Track constructor... EC, 21-Apr-2011.
          std::vector<TVector3> dircos(spacepoints.size());
          dircos[0] = TVector3(fpREC);
          dircos.back() = TVector3(fpRECL);

          // fill a vector of TVector3 with the space points used to make this track
          std::vector<TVector3> xyz(spacepoints.size());
          size_t spStart = spcol->size();
          for (size_t tv = 0; tv < spacepoints.size(); ++tv) {
            xyz[tv] = TVector3(spacepoints[tv]->XYZ());
            spcol->push_back(*(spacepoints[tv].get()));
          }
          size_t spEnd = spcol->size();

          tcol->push_back(
            recob::Track(recob::TrackTrajectory(recob::tracking::convertCollToPoint(xyz),
                                                recob::tracking::convertCollToVector(dircos),
                                                recob::Track::Flags_t(xyz.size()),
                                                false),
                         0,
                         -1.,
                         0,
                         recob::tracking::SMatrixSym55(),
                         recob::tracking::SMatrixSym55(),
                         tcol->size() - 1));

          // associate the track with its clusters and tracks
          util::CreateAssn(*this, evt, *tcol, clusters, *tcassn);
          util::CreateAssn(*this, evt, *tcol, hits, *thassn);

          // associate the track to the space points
          util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);

        } // getStatusFlag
      }   // spacepoints.size()>0

      //
      //std::cout<<"Track3DKalman found "<< tcol->size() <<" 3D track(s)"<<std::endl;
      if (trackIter != trackIn.end()) trackIter++;

      nTrks++;

    } // end loop over Track3Dreco/SpacePt tracks/groups (whichever) we brought into this event.

    evt.put(std::move(tcol));
    evt.put(std::move(spcol));
    evt.put(std::move(tcassn));
    evt.put(std::move(thassn));
    evt.put(std::move(tspassn));

  } // end method

} // end namespace trkf

DEFINE_ART_MODULE(trkf::Track3DKalman)
