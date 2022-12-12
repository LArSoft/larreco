///////////////////////////////////////////////////////////////////////////////////////////
//
// FeatureVertexFinder designed to analyze 2d & 3d verticies found in the TPC
//
// jaasaadi@syr.edu
//
// Note: This ana module will utilize MC truth information for verticies and is not (yet)
// intended as a unit test for data...though many of the same methods will likely be used
// for real data (when the time comes)
//////////////////////////////////////////////////////////////////////////////////////////

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larsim/MCCheater/BackTrackerService.h"

// C++
#include <string>

// Framework
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TH1F.h"
#include "TH2D.h"

namespace vertex {

  class FeatureVertexFinderAna : public art::EDAnalyzer {
  public:
    explicit FeatureVertexFinderAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void beginJob() override;

    std::string fLArG4ModuleLabel;
    std::string fGenieModuleLabel;
    std::string fVertexModuleLabel;
    std::string fEndPoint2dModuleLabel; //<---2d Vertex Module Label (EndPoint2d)

    // Outputting histograms for analysis

    TH1F* fRun;
    TH1F* fEvt;
    TH1F* fTruthVtxXPos;
    TH1F* fTruthVtxYPos;
    TH1F* fTruthVtxZPos;
    TH1F* fTruthWireNumberPlane0;
    TH1F* fTruthTimeTickPlane0;
    TH1F* fTruthWireNumberPlane1;
    TH1F* fTruthTimeTickPlane1;
    TH1F* fTruthWireNumberPlane2;
    TH1F* fTruthTimeTickPlane2;
    TH1F* fTruthWireInCmPlane0;
    TH1F* fTruthTimeInCmPlane0;
    TH1F* fTruthWireInCmPlane1;
    TH1F* fTruthTimeInCmPlane1;
    TH1F* fTruthWireInCmPlane2;
    TH1F* fTruthTimeInCmPlane2;

    TH1F* fTwoDNVtxPlane0;
    TH1F* fTwoDNVtxPlane1;
    TH1F* fTwoDNVtxPlane2;

    TH1F* fTwoDWireNumberPlane0;
    TH1F* fTwoDTimeTickPlane0;
    TH1F* fTwoDWireNumberPlane1;
    TH1F* fTwoDTimeTickPlane1;
    TH1F* fTwoDWireNumberPlane2;
    TH1F* fTwoDTimeTickPlane2;
    TH1F* fTwoDWireInCmPlane0;
    TH1F* fTwoDTimeInCmPlane0;
    TH1F* fTwoDWireInCmPlane1;
    TH1F* fTwoDTimeInCmPlane1;
    TH1F* fTwoDWireInCmPlane2;
    TH1F* fTwoDTimeInCmPlane2;

    TH1F* fTwoDStrengthPlane0;
    TH1F* fTwoDStrengthPlane1;
    TH1F* fTwoDStrengthPlane2;

    TH1F* fRecoCheck2dWireNumPlane0;
    TH1F* fRecoCheck2dTimeTickPlane0;
    TH1F* fRecoCheck2dWireInCmPlane0;
    TH1F* fRecoCheck2dTimeInCmPlane0;

    TH1F* fRecoCheck2dWireNumPlane1;
    TH1F* fRecoCheck2dTimeTickPlane1;
    TH1F* fRecoCheck2dWireInCmPlane1;
    TH1F* fRecoCheck2dTimeInCmPlane1;

    TH1F* fRecoCheck2dWireNumPlane2;
    TH1F* fRecoCheck2dTimeTickPlane2;
    TH1F* fRecoCheck2dWireInCmPlane2;
    TH1F* fRecoCheck2dTimeInCmPlane2;

    TH1F* fRecoVtxN3d;
    TH1F* fRecoVtxXPos;
    TH1F* fRecoVtxYPos;
    TH1F* fRecoVtxZPos;

    TH1F* fRecoCheck3dVtxX;
    TH1F* fRecoCheck3dVtxY;
    TH1F* fRecoCheck3dVtxZ;

    TH2D* fRecoCheck3dVtxXvsX;
    TH2D* fRecoCheck3dVtxYvsY;
    TH2D* fRecoCheck3dVtxZvsZ;

  }; // End class VertexFinderAna

  FeatureVertexFinderAna::FeatureVertexFinderAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    fLArG4ModuleLabel = pset.get<std::string>("LArGeantModuleLabel");
    fGenieModuleLabel = pset.get<std::string>("GenieModuleLabel");
    fVertexModuleLabel = pset.get<std::string>("VertexModuleLabel");
    fEndPoint2dModuleLabel = pset.get<std::string>("EndPoint2dModuleLabel");
  }

  void FeatureVertexFinderAna::beginJob()
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    // Outputting TH1F Histograms
    fRun = tfs->make<TH1F>("fRun", "Run Number", 1000, 0, 1000);
    fEvt = tfs->make<TH1F>("fEvt", "Event Number", 1000, 0, 1000);
    fTruthVtxXPos = tfs->make<TH1F>("fTruthVtxXPos", "Truth Vertex X Position", 400, 0, 200);
    fTruthVtxYPos = tfs->make<TH1F>("fTruthVtxYPos", "Truth Vertex Y Position", 400, -100, 100);
    fTruthVtxZPos = tfs->make<TH1F>("fTruthVtxZPos", "Truth Vertex Z Position", 2000, 0, 1000);
    fTruthWireNumberPlane0 =
      tfs->make<TH1F>("fTruthWireNumberPlane0", "Truth Wire Number Plane 0", 3000, 0, 3000);
    fTruthTimeTickPlane0 =
      tfs->make<TH1F>("fTruthTimeTickPlane0", "Truth Time Tick Plane 0", 3200, 0, 3200);
    fTruthWireNumberPlane1 =
      tfs->make<TH1F>("fTruthWireNumberPlane1", "Truth Wire Number Plane 1", 3000, 0, 3000);
    fTruthTimeTickPlane1 =
      tfs->make<TH1F>("fTruthTimeTickPlane1", "Truth Time Tick Plane 1", 3200, 0, 3200);
    fTruthWireNumberPlane2 =
      tfs->make<TH1F>("fTruthWireNumberPlane2", "Truth Wire Number Plane 2", 3000, 0, 3000);
    fTruthTimeTickPlane2 =
      tfs->make<TH1F>("fTruthTimeTickPlane2", "Truth Time Tick Plane 2", 3200, 0, 3200);
    fTruthWireInCmPlane0 =
      tfs->make<TH1F>("fTruthWireInCmPlane0", "Truth Wire In CM Plane 0", 2000, 0, 1000);
    fTruthTimeInCmPlane0 =
      tfs->make<TH1F>("fTruthTimeInCmPlane0", "Truth Time In Cm Plane 0", 2000, 0, 1000);
    fTruthWireInCmPlane1 =
      tfs->make<TH1F>("fTruthWireInCmPlane1", "Truth Wire In CM Plane 1", 2000, 0, 1000);
    fTruthTimeInCmPlane1 =
      tfs->make<TH1F>("fTruthTimeInCmPlane1", "Truth Time In Cm Plane 1", 2000, 0, 1000);
    fTruthWireInCmPlane2 =
      tfs->make<TH1F>("fTruthWireInCmPlane2", "Truth Wire In CM Plane 2", 2000, 0, 1000);
    fTruthTimeInCmPlane2 =
      tfs->make<TH1F>("fTruthTimeInCmPlane2", "Truth Time In Cm Plane 2", 2000, 0, 1000);

    fTwoDNVtxPlane0 =
      tfs->make<TH1F>("fTwoDNVtxPlane0", "TwoD Number of Verticies Found in Plane 0", 400, 0, 200);
    fTwoDNVtxPlane1 =
      tfs->make<TH1F>("fTwoDNVtxPlane1", "TwoD Number of Verticies Found in Plane 1", 400, 0, 200);
    fTwoDNVtxPlane2 =
      tfs->make<TH1F>("fTwoDNVtxPlane2", "TwoD Number of Verticies Found in Plane 2", 400, 0, 200);

    fTwoDWireNumberPlane0 =
      tfs->make<TH1F>("fTwoDWireNumberPlane0", "TwoD Wire Number Plane 0", 3000, 0, 3000);
    fTwoDTimeTickPlane0 =
      tfs->make<TH1F>("fTwoDTimeTickPlane0", "TwoD Time Tick Plane 0", 3200, 0, 3200);
    fTwoDWireNumberPlane1 =
      tfs->make<TH1F>("fTwoDWireNumberPlane1", "TwoD Wire Number Plane 1", 3000, 0, 3000);
    fTwoDTimeTickPlane1 =
      tfs->make<TH1F>("fTwoDTimeTickPlane1", "TwoD Time Tick Plane 1", 3200, 0, 3200);
    fTwoDWireNumberPlane2 =
      tfs->make<TH1F>("fTwoDWireNumberPlane2", "TwoD Wire Number Plane 2", 3000, 0, 3000);
    fTwoDTimeTickPlane2 =
      tfs->make<TH1F>("fTwoDTimeTickPlane2", "TwoD Time Tick Plane 2", 3200, 0, 3200);
    fTwoDWireInCmPlane0 =
      tfs->make<TH1F>("fTwoDWireInCmPlane0", "TwoD Wire In CM Plane 0", 2000, 0, 1000);
    fTwoDTimeInCmPlane0 =
      tfs->make<TH1F>("fTwoDTimeInCmPlane0", "TwoD Time In Cm Plane 0", 2000, 0, 1000);
    fTwoDWireInCmPlane1 =
      tfs->make<TH1F>("fTwoDWireInCmPlane1", "TwoD Wire In CM Plane 1", 2000, 0, 1000);
    fTwoDTimeInCmPlane1 =
      tfs->make<TH1F>("fTwoDTimeInCmPlane1", "TwoD Time In Cm Plane 1", 2000, 0, 1000);
    fTwoDWireInCmPlane2 =
      tfs->make<TH1F>("fTwoDWireInCmPlane2", "TwoD Wire In CM Plane 2", 2000, 0, 1000);
    fTwoDTimeInCmPlane2 =
      tfs->make<TH1F>("fTwoDTimeInCmPlane2", "TwoD Time In Cm Plane 2", 2000, 0, 1000);

    fTwoDStrengthPlane0 =
      tfs->make<TH1F>("fTwoDStrengthPlane0", "TwoD Strength Plane 0", 1000, 0, 500);
    fTwoDStrengthPlane1 =
      tfs->make<TH1F>("fTwoDStrengthPlane1", "TwoD Strength Plane 1", 1000, 0, 500);
    fTwoDStrengthPlane2 =
      tfs->make<TH1F>("fTwoDStrengthPlane2", "TwoD Strength Plane 2", 1000, 0, 500);

    fRecoCheck2dWireNumPlane0 = tfs->make<TH1F>(
      "fRecoCheck2dWireNumPlane0", "Reco Wire Number - True Wire Number Plane 0", 400, -200, 200);
    fRecoCheck2dTimeTickPlane0 = tfs->make<TH1F>(
      "fRecoCheck2dTimeTickPlane0", "Reco Time Tick - True Time Tick Plane 0", 1000, -500, 500);
    fRecoCheck2dWireInCmPlane0 = tfs->make<TH1F>(
      "fRecoCheck2dWireInCmPlane0", "Reco Wire in CM - True Wire in CM Plane 0", 200, -50, 50);
    fRecoCheck2dTimeInCmPlane0 = tfs->make<TH1F>(
      "fRecoCheck2dTimeInCmPlane0", "Reco Time in CM - True Time in CM Plane 0", 200, -50, 50);

    fRecoCheck2dWireNumPlane1 = tfs->make<TH1F>(
      "fRecoCheck2dWireNumPlane1", "Reco Wire Number - True Wire Number Plane 1", 400, -200, 200);
    fRecoCheck2dTimeTickPlane1 = tfs->make<TH1F>(
      "fRecoCheck2dTimeTickPlane1", "Reco Time Tick - True Time Tick Plane 1", 1000, -500, 500);
    fRecoCheck2dWireInCmPlane1 = tfs->make<TH1F>(
      "fRecoCheck2dWireInCmPlane1", "Reco Wire in CM - True Wire in CM Plane 1", 200, -50, 50);
    fRecoCheck2dTimeInCmPlane1 = tfs->make<TH1F>(
      "fRecoCheck2dTimeInCmPlane1", "Reco Time in CM - True Time in CM Plane 1", 200, -50, 50);

    fRecoCheck2dWireNumPlane2 = tfs->make<TH1F>(
      "fRecoCheck2dWireNumPlane2", "Reco Wire Number - True Wire Number Plane 2", 400, -200, 200);
    fRecoCheck2dTimeTickPlane2 = tfs->make<TH1F>(
      "fRecoCheck2dTimeTickPlane2", "Reco Time Tick - True Time Tick Plane 2", 1000, -500, 500);
    fRecoCheck2dWireInCmPlane2 = tfs->make<TH1F>(
      "fRecoCheck2dWireInCmPlane2", "Reco Wire in CM - True Wire in CM Plane 2", 200, -50, 50);
    fRecoCheck2dTimeInCmPlane2 = tfs->make<TH1F>(
      "fRecoCheck2dTimeInCmPlane2", "Reco Time in CM - True Time in CM Plane 2", 200, -50, 50);

    fRecoVtxN3d = tfs->make<TH1F>("fRecoVtxN3d", "Number of 3d-Reco Verticies", 400, 0, 200);
    fRecoVtxXPos = tfs->make<TH1F>("fRecoVtxXPos", "Reco Vertex X Position", 400, -10, 200);
    fRecoVtxYPos = tfs->make<TH1F>("fRecoVtxYPos", "Reco Vertex Y Position", 400, -100, 100);
    fRecoVtxZPos = tfs->make<TH1F>("fRecoVtxZPos", "Reco Vertex Z Position", 2000, -10, 1000);

    fRecoCheck3dVtxX =
      tfs->make<TH1F>("fRecoCheck3dVtxX", "Reco X Position - True X Postion", 800, -100, 100);
    fRecoCheck3dVtxY =
      tfs->make<TH1F>("fRecoCheck3dVtxY", "Reco Y Position - True Y Postion", 800, -100, 100);
    fRecoCheck3dVtxZ =
      tfs->make<TH1F>("fRecoCheck3dVtxZ", "Reco Z Position - True Z Postion", 800, -100, 100);

    fRecoCheck3dVtxXvsX = tfs->make<TH2D>(
      "fRecoCheck3dVtxXvsX", "(Reco X - True X)/True X vs True X", 400, -10, 200, 120, -20, 20);
    fRecoCheck3dVtxYvsY = tfs->make<TH2D>(
      "fRecoCheck3dVtxYvsY", "(Reco Y - True Y)/True Y vs True Y", 400, -100, 100, 120, -20, 20);
    fRecoCheck3dVtxZvsZ = tfs->make<TH2D>(
      "fRecoCheck3dVtxZvsZ", "(Reco Z - True Z)/True Z vs True Z", 1000, -10, 1000, 120, -20, 20);
  }

  void FeatureVertexFinderAna::analyze(const art::Event& evt)
  {
    // Filling Run/Event Histo
    fRun->Fill(evt.run());
    fEvt->Fill(evt.id().event());

    // Getting the Geant Information Directly
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle);

    // Getting MC Truth Info from simb
    auto const& mctruthList = *evt.getValidHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);

    // Getting information from BackTrackerService
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

    art::ServiceHandle<geo::Geometry const> geom;
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

    // Getting 2d Vertex information (vertex2dHandle)
    art::Handle<std::vector<recob::EndPoint2D>> vertex2dHandle;
    evt.getByLabel(fEndPoint2dModuleLabel, vertex2dHandle);

    // Getting the 3d Vertex (vertex3dListHandle)
    art::Handle<std::vector<recob::Vertex>> vertex3dListHandle;
    evt.getByLabel(fVertexModuleLabel, vertex3dListHandle);

    // Detector numbers that are useful (to be set later)
    std::vector<double> WirePitch_CurrentPlane(geom->Views().size(),
                                               0.); //<---Setting the Wire pitch for each plane
                                                    // Right now assuming only 3 planes
    // get the wire pitch for each view
    size_t vn = 0;
    for (auto v : geom->Views()) {
      WirePitch_CurrentPlane[vn] = geom->Plane({0, 0, v}).WirePitch();
      ++vn;
    }

    // Calculating the Timetick to CM conversion
    float TimeTick =
      sampling_rate(clock_data) / 1000.; //<---To get units of microsecond...not nanosec
    float DriftVelocity = det_prop.DriftVelocity(det_prop.Efield(), det_prop.Temperature());
    float TimetoCm = TimeTick * DriftVelocity;

    // Looping over MC information

    // Variables for truth vertex information
    geo::Point_t truth_vertex{}; //<---Truth x,y,z information

    uint32_t VtxWireNum[3] = {0}; //<---Wire number in each plane ( WireNum[plane#] )
    double VtxTimeTick[3] = {0.}; //<---Time tick in each plane   ( TimeTick[plane#] )

    double VtxWireNum_InCM[3] = {0.};  //<---Wire number in each plane in CM ( WireNum[plane#] )
    double VtxTimeTick_InCM[3] = {0.}; //<---Time tick in each plane in CM   ( TimeTick[plane#] )

    // Finding the MC truth vertex

    // FIXME: Can this be right?  We're looping over all of the
    //        MCTruth particles just so that we can use the last one
    //        (i.e. truth_vertex is set each time but only the last
    //        iteration is relevant)?
    for (std::size_t i = 0; i != mctruthList.size(); ++i) {
      auto const& neut = mctruthList[i].GetParticle(i);
      truth_vertex.SetX(neut.Vx());
      truth_vertex.SetY(neut.Vy());
      truth_vertex.SetZ(neut.Vz());
    } // end i loop

    // Filling Histograms
    fTruthVtxXPos->Fill(truth_vertex.X());
    fTruthVtxYPos->Fill(truth_vertex.Y());
    fTruthVtxZPos->Fill(truth_vertex.Z());

    // Looping over geo::PlaneIDs
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {
      auto const& pid = plane.ID();
      // Calculating the nearest wire the vertex corresponds to in each plane
      try {
        VtxWireNum[pid.Plane] = plane.NearestWireID(truth_vertex).Wire;
      }
      catch (...) {
        mf::LogWarning("FeatureVertexFinderAna") << "Can't find nearest wire";
        continue;
      }
      VtxTimeTick[pid.Plane] =
        det_prop.ConvertXToTicks(truth_vertex.X(), pid) + det_prop.GetXTicksOffset(pid);

      // Translating each of these in cm
      VtxWireNum_InCM[pid.Plane] = VtxWireNum[pid.Plane] * WirePitch_CurrentPlane[pid.Plane];
      VtxTimeTick_InCM[pid.Plane] = VtxTimeTick[pid.Plane] * TimetoCm;
    } //<---End pid loop

    // Filling Histograms
    fTruthWireNumberPlane0->Fill(VtxWireNum[0]);
    fTruthTimeTickPlane0->Fill(VtxTimeTick[0]);
    fTruthWireNumberPlane1->Fill(VtxWireNum[1]);
    fTruthTimeTickPlane1->Fill(VtxTimeTick[1]);
    fTruthWireNumberPlane2->Fill(VtxWireNum[2]);
    fTruthTimeTickPlane2->Fill(VtxTimeTick[2]);

    fTruthWireInCmPlane0->Fill(VtxWireNum_InCM[0]);
    fTruthTimeInCmPlane0->Fill(VtxTimeTick_InCM[0]);
    fTruthWireInCmPlane1->Fill(VtxWireNum_InCM[1]);
    fTruthTimeInCmPlane1->Fill(VtxTimeTick_InCM[1]);
    fTruthWireInCmPlane2->Fill(VtxWireNum_InCM[2]);
    fTruthTimeInCmPlane2->Fill(VtxTimeTick_InCM[2]);

    // Looping over EndPoint2d information

    art::PtrVector<recob::EndPoint2D> vert2d;

    // Variables for Vertex2d
    double Vertex2d_TimeTick[10000] = {
      0.}; //<---Vertex2d Time Tick for the current plane ( TimeTick[#2d] )
    double Vertex2d_Wire[10000] = {0.}; //<---Veretx2d Wire # ( Wire[#2d] )

    double Vertex2d_TimeTick_InCM[10000] = {0.}; //<---Vertex 2d Time tick in CM ( TimeTick[#2d] )
    double Vertex2d_Wire_InCM[10000] = {0.};     //<---Veretx2d Wire in CM ( Wire[#2d] )
    int n2dVtx = 0;

    int n2dVtxPlane0 = 0, n2dVtxPlane1 = 0, n2dVtxPlane2 = 0;

    bool vertexWstrengthplane0 = false,
         vertexWstrengthplane1 =
           false; //, vertexWstrengthplane2 = false;   //commented out, Wes, 12.4.13

    // Loop over the EndPoint2d List
    for (size_t ii = 0; ii < vertex2dHandle->size(); ++ii) {
      art::Ptr<recob::EndPoint2D> vertex(vertex2dHandle, ii);
      vert2d.push_back(vertex);
    }

    // If we have 2d vertex, loop over them
    if (vert2d.size() > 0) {

      // Looping over geo::PlaneIDs
      for (auto const& pid : geom->Iterate<geo::PlaneID>()) {
        for (size_t ww = 0; ww < vert2d.size(); ++ww) {
          // Only look at this 2d vertex if it is in the current plane
          if (vert2d[ww]->WireID().planeID() != pid) { continue; }

          Vertex2d_TimeTick[n2dVtx] = vert2d[ww]->DriftTime();
          Vertex2d_Wire[n2dVtx] = vert2d[ww]->WireID().Wire;

          // Translating each of these in cm
          Vertex2d_Wire_InCM[n2dVtx] = Vertex2d_Wire[n2dVtx] * WirePitch_CurrentPlane[pid.Plane];
          Vertex2d_TimeTick_InCM[n2dVtx] = Vertex2d_TimeTick[n2dVtx] * TimetoCm;

          // Checking how well we did in reconstructing the vertex (Reco - True)

          float RecoCheck_TimeTick = Vertex2d_TimeTick[n2dVtx] - VtxTimeTick[pid.Plane];
          float RecoCheck_WireNum = Vertex2d_Wire[n2dVtx] - VtxWireNum[pid.Plane];

          float RecoCheck_TimeInCm = Vertex2d_TimeTick_InCM[n2dVtx] - VtxTimeTick_InCM[pid.Plane];
          float RecoCheck_WireInCm = Vertex2d_Wire_InCM[n2dVtx] - VtxWireNum_InCM[pid.Plane];

          if (vert2d[ww]->Strength() > -1) {
            if (pid.Plane == 0) {
              vertexWstrengthplane0 = true;

              fTwoDWireNumberPlane0->Fill(Vertex2d_Wire[n2dVtx]);
              fTwoDTimeTickPlane0->Fill(Vertex2d_TimeTick[n2dVtx]);
              fTwoDWireInCmPlane0->Fill(Vertex2d_Wire_InCM[n2dVtx]);
              fTwoDTimeInCmPlane0->Fill(Vertex2d_TimeTick_InCM[n2dVtx]);

              fTwoDStrengthPlane0->Fill(vert2d[ww]->Strength());

              fRecoCheck2dWireNumPlane0->Fill(RecoCheck_WireNum);
              fRecoCheck2dTimeTickPlane0->Fill(RecoCheck_TimeTick);
              fRecoCheck2dWireInCmPlane0->Fill(RecoCheck_WireInCm);
              fRecoCheck2dTimeInCmPlane0->Fill(RecoCheck_TimeInCm);

              n2dVtxPlane0++;

            } //<---End Plane 0

            if (pid.Plane == 1) {
              vertexWstrengthplane1 = true;
              fTwoDWireNumberPlane1->Fill(Vertex2d_Wire[n2dVtx]);
              fTwoDTimeTickPlane1->Fill(Vertex2d_TimeTick[n2dVtx]);
              fTwoDWireInCmPlane1->Fill(Vertex2d_Wire_InCM[n2dVtx]);
              fTwoDTimeInCmPlane1->Fill(Vertex2d_TimeTick_InCM[n2dVtx]);

              fTwoDStrengthPlane1->Fill(vert2d[ww]->Strength());

              fRecoCheck2dWireNumPlane1->Fill(RecoCheck_WireNum);
              fRecoCheck2dTimeTickPlane1->Fill(RecoCheck_TimeTick);
              fRecoCheck2dWireInCmPlane1->Fill(RecoCheck_WireInCm);
              fRecoCheck2dTimeInCmPlane1->Fill(RecoCheck_TimeInCm);

              n2dVtxPlane1++;

            } //<---End Plane 1

            if (pid.Plane == 2) {
              fTwoDWireNumberPlane2->Fill(Vertex2d_Wire[n2dVtx]);
              fTwoDTimeTickPlane2->Fill(Vertex2d_TimeTick[n2dVtx]);
              fTwoDWireInCmPlane2->Fill(Vertex2d_Wire_InCM[n2dVtx]);
              fTwoDTimeInCmPlane2->Fill(Vertex2d_TimeTick_InCM[n2dVtx]);

              fTwoDStrengthPlane2->Fill(vert2d[ww]->Strength());

              fRecoCheck2dWireNumPlane2->Fill(RecoCheck_WireNum);
              fRecoCheck2dTimeTickPlane2->Fill(RecoCheck_TimeTick);
              fRecoCheck2dWireInCmPlane2->Fill(RecoCheck_WireInCm);
              fRecoCheck2dTimeInCmPlane2->Fill(RecoCheck_TimeInCm);

              n2dVtxPlane2++;

            } //<---End Plane 2
          }

          ++n2dVtx;
        } //<--end ww loop

      } //<---End plane ID

    } //<--End checking if we have 2d end points

    fTwoDNVtxPlane0->Fill(n2dVtxPlane0);
    fTwoDNVtxPlane1->Fill(n2dVtxPlane1);
    fTwoDNVtxPlane2->Fill(n2dVtxPlane2);

    // Looping over 3dVertex information
    art::PtrVector<recob::Vertex> Vertexlist;

    double xyz[3] = {0.};

    for (unsigned int ii = 0; ii < vertex3dListHandle->size(); ++ii) {
      art::Ptr<recob::Vertex> vertex3d(vertex3dListHandle, ii);
      Vertexlist.push_back(vertex3d); // class member

    } //<---End ii loop

    if (Vertexlist.size() > 0 && vertexWstrengthplane0 && vertexWstrengthplane1) {

      fRecoVtxN3d->Fill(Vertexlist.size());
      for (unsigned int ww = 0; ww < Vertexlist.size(); ww++) {
        Vertexlist[ww]->XYZ(xyz);

        // Checking how well we did in reconstructing the vertex (Reco - True)

        // Finding the Delta X, Y, Z between Reco vtx and truth
        double DeltaX = xyz[0] - truth_vertex.X();
        double DeltaY = xyz[1] - truth_vertex.Y();
        double DeltaZ = xyz[2] - truth_vertex.Z();

        double DeltaXoverTrueX = DeltaX / truth_vertex.X();
        double DeltaYoverTrueY = DeltaY / truth_vertex.X();
        double DeltaZoverTrueZ = DeltaZ / truth_vertex.X();

        fRecoVtxXPos->Fill(xyz[0]);
        fRecoVtxYPos->Fill(xyz[1]);
        fRecoVtxZPos->Fill(xyz[2]);

        fRecoCheck3dVtxX->Fill(DeltaX);
        fRecoCheck3dVtxY->Fill(DeltaY);
        fRecoCheck3dVtxZ->Fill(DeltaZ);

        fRecoCheck3dVtxXvsX->Fill(xyz[0], DeltaXoverTrueX);
        fRecoCheck3dVtxYvsY->Fill(xyz[1], DeltaYoverTrueY);
        fRecoCheck3dVtxZvsZ->Fill(xyz[2], DeltaZoverTrueZ);
      }
    }
  } //<---End access the event

  DEFINE_ART_MODULE(FeatureVertexFinderAna)

} // end of vertex namespace
