////////////////////////////////////////////////////////////////////////
//
// \file ShowerReco_module.cc
//
// biagio.rossi@lhep.unibe.ch   (FWMK : argoslope.resize(fNPlanes);neut specific)
// thomas.strauss@lhep.unibe.ch (ART  : general detector)
//
// andrzej.szelc@yale.edu (port to detector agnostic version)
// jasaadi@fnal.gov (Try to rewrite the code to make it more readable and to be able to handle
//                   multiple TPC's and cryostats and more than one cluster per plane.)
//
// This algorithm is designed to reconstruct showers
//
///////////////////////////////////////////////////////////////////////

// ### Generic C++ includes ###
#include <cmath> // std::tan() ...
#include <string>
#include <vector>

// ### Framework includes ###
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

// ### ROOT includes ###
#include "TMath.h"
#include "TTree.h"

// ### LArSoft includes ###
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/ArtDataHelper/ToElement.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "range/v3/view.hpp"

#include <cmath>

namespace shwf {

  class ShowerReco : public art::EDProducer {
  public:
    explicit ShowerReco(fhicl::ParameterSet const& pset);

  private:
    void beginJob();
    void beginRun(art::Run& run);
    void produce(art::Event& evt); /**Actual routine that reconstruct the shower*/

    void GetVertexAndAnglesFromCluster(art::Ptr<recob::Cluster> clust,
                                       unsigned int plane); // Get shower vertex and slopes.

    void LongTransEnergy(geo::GeometryCore const* geom,
                         detinfo::DetectorClocksData const& clockData,
                         detinfo::DetectorPropertiesData const& detProp,
                         unsigned int set,
                         std::vector<art::Ptr<recob::Hit>> hitlist); // Longitudinal

    void ClearandResizeVectors(unsigned int nPlanes);

    int fRun, fEvent, fSubRun;
    // input labels:
    float slope[3]; // in cm, cm
    float angle[3]; // in radians

    std::string fClusterModuleLabel;

    float ftimetick; // time sample in us

    double fMean_wire_pitch; // wire pitch in cm
    fhicl::ParameterSet fCaloPSet;

    std::vector<double> fRMS_2cm;
    std::vector<int> fNpoints_2cm;
    std::vector<double> fCorr_MeV_2cm;
    std::vector<double> fCorr_Charge_2cm;

    std::vector<int> fNpoints_corr_ADC_2cm;
    std::vector<int> fNpoints_corr_MeV_2cm;

    std::vector<double> fTotChargeADC;      //Total charge in ADC/cm for each plane
    std::vector<double> fTotChargeMeV;      //Total charge in MeV/cm for each plane
    std::vector<double> fTotChargeMeV_MIPs; //Total charge in MeV/cm for each plane

    std::vector<double> fChargeADC_2cm; //Initial charge in ADC/cm for each plane first 4cm
    std::vector<double> fChargeMeV_2cm; //initial charge in MeV/cm for each plane first 4cm

    std::vector<double> fChargeMeV_2cm_refined;
    std::vector<double> fChargeMeV_2cm_axsum;

    std::vector<std::vector<double>> fDistribChargeADC; //vector with the first De/Dx points ADC
    std::vector<std::vector<double>>
      fDistribChargeMeV; //vector with the first De/Dx points converted energy
    std::vector<std::vector<double>> fDistribHalfChargeMeV;
    std::vector<std::vector<double>>
      fDistribChargeposition; //vector with the first De/Dx points' positions

    std::vector<std::vector<double>> fSingleEvtAngle;    //vector with the first De/Dx points
    std::vector<std::vector<double>> fSingleEvtAngleVal; //vector with the first De/Dx points

    std::vector<unsigned int> fWire_vertex; // wire coordinate of vertex for each plane
    std::vector<double> fTime_vertex;       // time coordinate of vertex for each plane

    std::vector<double> fWire_vertexError; // wire coordinate of vertex for each plane
    std::vector<double> fTime_vertexError; // time coordinate of vertex for each plane

    std::vector<unsigned int> fWire_last; // wire coordinate of last point for each plane
    std::vector<double> fTime_last;       // time coordinate of last point for each plane

    // reconstructed 3D position of shower start
    std::vector<double> xyz_vertex_fit;

    std::vector<std::vector<double>>
      fNPitch; // double array, to use each plane for each set of angles

    // calorimetry variables
    float Kin_En;
    std::vector<float> vdEdx;
    std::vector<float> vresRange;
    std::vector<float> vdQdx;
    std::vector<float> deadwire; // residual range for dead wires
    float Trk_Length;
    float fTrkPitchC;
    float fdEdxlength;     // distance that gets used to determine e/gamma
                           // separation
    float fcalodEdxlength; // cutoff distance for hits saved to the calo object.
    bool fUseArea;

    double xphi, xtheta;   // new calculated angles.
    unsigned int fNPlanes; // number of planes
    unsigned int fNAngles;
    TTree* ftree_shwf;

    // conversion and useful constants
    double fWirePitch; // wire pitch in cm
    double fTimeTick;
    double fDriftVelocity;
    double fWireTimetoCmCm;

    std::vector<int> fNhitsperplane;
    std::vector<double> fTotADCperplane;
  }; // class ShowerReco

  //------------------------------------------------------------------------------
  ShowerReco::ShowerReco(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fCaloPSet = pset.get<fhicl::ParameterSet>("CaloAlg");

    fdEdxlength =
      pset.get<double>("dEdxlength"); // distance that gets used to determine e/gamma separation
    fcalodEdxlength =
      pset.get<double>("calodEdxlength"); // cutoff distance for hits saved to the calo object.
    fUseArea = pset.get<bool>("UseArea");

    produces<std::vector<recob::Shower>>();
    produces<art::Assns<recob::Shower, recob::Cluster>>();
    produces<art::Assns<recob::Shower, recob::Hit>>();
    produces<std::vector<anab::Calorimetry>>();
    produces<art::Assns<recob::Shower, anab::Calorimetry>>();
  }

  // ***************** //
  void ShowerReco::beginJob()
  {
    art::ServiceHandle<geo::Geometry const> geo;

    /// \todo the call to geo->Nplanes() assumes this is a single cryostat and
    /// single TPC detector; need to generalize to multiple cryostats and
    /// TPCs
    fNPlanes = geo->Nplanes();
    fMean_wire_pitch = geo->WirePitch(); // wire pitch in cm

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    ftimetick = sampling_rate(clockData) / 1000.;

    /**Get TFileService and define output Histograms*/
    art::ServiceHandle<art::TFileService const> tfs;

    ftree_shwf = tfs->make<TTree>("ShowerReco",
                                  "Results"); /**All-knowing tree with reconstruction information*/
    ftree_shwf->Branch("run", &fRun, "run/I");
    ftree_shwf->Branch("subrun", &fSubRun, "subrun/I");
    ftree_shwf->Branch("event", &fEvent, "event/I");
    ftree_shwf->Branch("nplanes", &fNPlanes, "nplanes/I");
    ftree_shwf->Branch("nangles", &fNAngles, "nangles/I");
    ftree_shwf->Branch("xtheta", &xtheta, "xtheta/D");
    ftree_shwf->Branch("xphi", &xphi, "xphi/D");
    ftree_shwf->Branch("ftotChargeADC", "std::vector<double>", &fTotChargeADC);
    ftree_shwf->Branch("ftotChargeMeV", "std::vector<double>", &fTotChargeMeV);
    ftree_shwf->Branch("fTotChargeMeV_MIPs", "std::vector<double>", &fTotChargeMeV_MIPs);
    ftree_shwf->Branch("NPitch", "std::vector< std::vector<double> >", &fNPitch);

    // this should be temporary - until the omega is sorted out.
    ftree_shwf->Branch("RMS_2cm", "std::vector<double>", &fRMS_2cm);
    ftree_shwf->Branch("Npoints_2cm", "std::vector<int>", &fNpoints_2cm);
    ftree_shwf->Branch("ChargeADC_2cm", "std::vector<double>", &fChargeADC_2cm);
    ftree_shwf->Branch("ChargeMeV_2cm", "std::vector<double>", &fChargeMeV_2cm);
    ftree_shwf->Branch("ChargeMeV_2cm_refined", "std::vector<double>", &fChargeMeV_2cm_refined);
    ftree_shwf->Branch("ChargeMeV_2cm_axsum", "std::vector<double>", &fChargeMeV_2cm_axsum);
    ftree_shwf->Branch("fNhitsperplane", "std::vector<int>", &fNhitsperplane);
    ftree_shwf->Branch("fTotADCperplane", "std::vector<double>", &fTotADCperplane);
    ftree_shwf->Branch(
      "ChargedistributionADC", "std::vector<std::vector<double>>", &fDistribChargeADC);
    ftree_shwf->Branch(
      "ChargedistributionMeV", "std::vector<std::vector<double>>", &fDistribChargeMeV);
    ftree_shwf->Branch(
      "DistribHalfChargeMeV", "std::vector<std::vector<double>>", &fDistribHalfChargeMeV);
    ftree_shwf->Branch(
      "ChargedistributionPosition", "std::vector<std::vector<double>>", &fDistribChargeposition);
    ftree_shwf->Branch("xyz_vertex_fit", "std::vector<double>", &xyz_vertex_fit);
  }

  void ShowerReco::beginRun(art::Run&)
  {
    auto const* geom = lar::providerFrom<geo::Geometry>();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();

    fWirePitch = geom->WirePitch(); // wire pitch in cm
    fTimeTick = sampling_rate(clockData) / 1000.;
    fDriftVelocity = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
    fWireTimetoCmCm = (fTimeTick * fDriftVelocity) / fWirePitch;
  }

  void ShowerReco::ShowerReco::ClearandResizeVectors(unsigned int /*nPlanes*/)
  {
    //calculate factorial for number of angles
    int fact = 1;
    for (unsigned int i = 1; i <= fNPlanes; ++i)
      fact *= i;

    fNAngles = fact / 2;

    fDistribChargeADC.clear();
    fDistribChargeMeV.clear();
    fDistribChargeposition.clear();

    fDistribChargeADC.resize(fNPlanes);
    fDistribChargeMeV.resize(fNPlanes);
    fDistribChargeposition.resize(fNPlanes);

    fNPitch.clear();
    fDistribChargeADC.clear();
    fDistribChargeMeV.clear();
    fDistribHalfChargeMeV.clear();
    fDistribChargeposition.clear();

    fNPitch.resize(fNAngles);
    fDistribChargeADC.resize(fNPlanes);
    fDistribChargeMeV.resize(fNPlanes);
    fDistribHalfChargeMeV.resize(fNPlanes);
    fDistribChargeposition.resize(fNPlanes);

    for (unsigned int ii = 0; ii < fNAngles; ii++) {
      fNPitch[ii].resize(fNPlanes, -1);
    }

    fNpoints_corr_ADC_2cm.clear();
    fNpoints_corr_MeV_2cm.clear();

    fNpoints_corr_ADC_2cm.resize(fNAngles, -1);
    fNpoints_corr_MeV_2cm.resize(fNAngles, -1);

    for (unsigned int ii = 0; ii < fNPlanes; ii++) {
      fDistribChargeADC[ii].resize(0); //vector with the first De/Dx points
      fDistribChargeMeV[ii].resize(0); //vector with the first De/Dx points
      fDistribHalfChargeMeV[ii].resize(0);
      fDistribChargeposition[ii].resize(0); //vector with the first De/Dx points' positions
    }

    fWire_vertex.clear();
    fTime_vertex.clear();
    fWire_vertexError.clear();
    fTime_vertexError.clear();
    fWire_last.clear();
    fTime_last.clear();
    fTotChargeADC.clear();
    fTotChargeMeV.clear();
    fTotChargeMeV_MIPs.clear();
    fRMS_2cm.clear();
    fNpoints_2cm.clear();

    fNhitsperplane.clear();
    fTotADCperplane.clear();

    fWire_vertex.resize(fNAngles, -1);
    fTime_vertex.resize(fNAngles, -1);
    fWire_vertexError.resize(fNPlanes, -1);
    fTime_vertexError.resize(fNPlanes, -1);
    fWire_last.resize(fNAngles, -1);
    fTime_last.resize(fNAngles, -1);
    fTotChargeADC.resize(fNAngles, 0);
    fTotChargeMeV.resize(fNAngles, 0);
    fTotChargeMeV_MIPs.resize(fNAngles, 0);
    fRMS_2cm.resize(fNAngles, 0);
    fNpoints_2cm.resize(fNAngles, 0);

    fCorr_MeV_2cm.clear();
    fCorr_Charge_2cm.clear();
    xyz_vertex_fit.clear();

    fChargeADC_2cm.clear();
    fChargeMeV_2cm.clear();
    fChargeMeV_2cm_refined.clear();
    fChargeMeV_2cm_refined.clear();
    fChargeMeV_2cm_axsum.clear();

    fCorr_MeV_2cm.resize(fNAngles, 0);
    fCorr_Charge_2cm.resize(fNAngles, 0);
    xyz_vertex_fit.resize(3);

    fChargeADC_2cm.resize(fNAngles,
                          0); //Initial charge in ADC/cm for each plane angle calculation 4cm
    fChargeMeV_2cm.resize(fNAngles,
                          0); //initial charge in MeV/cm for each angle calculation first 4cm
    fChargeMeV_2cm_refined.resize(0);
    fChargeMeV_2cm_refined.resize(fNAngles, 0);
    fChargeMeV_2cm_axsum.resize(fNAngles, 0);

    fNhitsperplane.resize(fNPlanes, -1);
    fTotADCperplane.resize(fNPlanes, -1);

    vdEdx.clear();
    vresRange.clear();
    vdQdx.clear();
  }

  // ***************** //
  void ShowerReco::produce(art::Event& evt)
  {
    auto const* geom = lar::providerFrom<geo::Geometry>();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    util::GeometryUtilities const gser{*geom, clockData, detProp};
    constexpr geo::TPCID tpcid{0, 0};
    fNPlanes = geom->Nplanes(tpcid);
    auto Shower3DVector = std::make_unique<std::vector<recob::Shower>>();
    auto cassn = std::make_unique<art::Assns<recob::Shower, recob::Cluster>>();
    auto hassn = std::make_unique<art::Assns<recob::Shower, recob::Hit>>();
    auto calorimetrycol = std::make_unique<std::vector<anab::Calorimetry>>();
    auto calassn = std::make_unique<art::Assns<anab::Calorimetry, recob::Shower>>();

    /**Get Clusters*/

    art::Handle<std::vector<recob::Cluster>> clusterListHandle;
    evt.getByLabel(fClusterModuleLabel, clusterListHandle);

    art::Handle<std::vector<art::PtrVector<recob::Cluster>>> clusterAssociationHandle;
    evt.getByLabel(fClusterModuleLabel, clusterAssociationHandle);

    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);

    fRun = evt.id().run();
    fSubRun = evt.id().subRun();
    fEvent = evt.id().event();

    // find all the hits associated to all the clusters (once and for all);
    // the index of the query matches the index of the cluster in the collection
    // (conveniently carried around in its art pointer)
    art::FindManyP<recob::Hit> ClusterHits(clusterListHandle, evt, fClusterModuleLabel);

    std::vector<art::PtrVector<recob::Cluster>>::const_iterator clusterSet =
      clusterAssociationHandle->begin();
    // loop over vector of vectors (each size of NPlanes) and reconstruct showers from each of those
    for (size_t iClustSet = 0; iClustSet < clusterAssociationHandle->size(); iClustSet++) {

      const art::PtrVector<recob::Cluster> CurrentClusters = (*(clusterSet++));

      // do some error checking - i.e. are the clusters themselves present.
      if (clusterListHandle->size() < 2 || CurrentClusters.size() < 2) {
        ftree_shwf->Fill();
        return;
      }

      ClearandResizeVectors(fNPlanes);

      std::vector<std::vector<art::Ptr<recob::Hit>>> hitlist_all;
      hitlist_all.resize(fNPlanes);

      for (size_t iClust = 0; iClust < CurrentClusters.size(); iClust++) {
        art::Ptr<recob::Cluster> const& pclust = CurrentClusters[iClust];

        // get all the hits for this cluster;
        // pclust is a art::Ptr to the original cluster collection stored in the event;
        // its key corresponds to its index in the collection
        // (and therefore to the index in the query)
        std::vector<art::Ptr<recob::Hit>> const& hitlist = ClusterHits.at(pclust.key());

        unsigned int p(0); //c=channel, p=plane, w=wire

        if (hitlist.size() == 0) continue;

        p = (*hitlist.begin())->WireID().Plane;
        // get vertex position and slope information to start with - ii is the
        // posistion of the correct cluster:
        GetVertexAndAnglesFromCluster(pclust, p);

        double ADCcharge = 0;
        //loop over cluster hits
        for (art::Ptr<recob::Hit> const& hit : hitlist) {
          p = hit->WireID().Plane;
          hitlist_all[p].push_back(hit);
          ADCcharge += hit->PeakAmplitude();
        }
        fNhitsperplane[p] = hitlist_all[p].size();
        fTotADCperplane[p] = ADCcharge;
      } // End loop on clusters.
      // Now I have the Hitlists and the relevent clusters parameters saved.

      // find best set:
      unsigned int bp1 = 0, bp2 = 0;
      double minerror1 = 99999999, minerror2 = 9999999;
      for (unsigned int ii = 0; ii < fNPlanes; ++ii) {
        double locerror =
          fWire_vertexError[ii] * fWire_vertexError[ii] +
          fTime_vertexError[ii] * fTime_vertexError[ii]; // time coordinate of vertex for each plane

        if (minerror1 >= locerror) // the >= sign is to favorize collection
        {
          minerror1 = locerror;
          bp1 = ii;
        }
      }
      for (unsigned int ij = 0; ij < fNPlanes; ++ij) {
        double locerror =
          fWire_vertexError[ij] * fWire_vertexError[ij] +
          fTime_vertexError[ij] * fTime_vertexError[ij]; // time coordinate of vertex for each plane

        if (minerror2 >= locerror && ij != bp1) {
          minerror2 = locerror;
          bp2 = ij;
        }
      }

      gser.Get3DaxisN(bp1, bp2, angle[bp1], angle[bp2], xphi, xtheta);

      ///////////////////////////////////////////////////////////
      std::vector<geo::Point_t> position;
      position.reserve(fNPlanes);
      // get starting positions for all planes -- FIXME: only position[0] is used.
      for (auto const& plane : geom->Iterate<geo::PlaneGeo>(tpcid)) {
        position.push_back(plane.GetBoxCenter());
      }

      // Assuming there is no problem ( and we found the best pair that comes
      // close in time ) we try to get the Y and Z coordinates for the start of
      // the shower.
      double fTimeTick = sampling_rate(clockData) / 1000.;
      double fDriftVelocity = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
      try {
        int chan1 = geom->PlaneWireToChannel({0, 0, bp1, fWire_vertex[bp1]});
        int chan2 = geom->PlaneWireToChannel({0, 0, bp2, fWire_vertex[bp2]});

        double y, z;
        geom->ChannelsIntersect(chan1, chan2, y, z);

        xyz_vertex_fit[1] = y;
        xyz_vertex_fit[2] = z;
        xyz_vertex_fit[0] =
          (fTime_vertex[bp1] - trigger_offset(clockData)) * fDriftVelocity * fTimeTick +
          position[0].X();
      }
      catch (cet::exception const& e) {
        mf::LogWarning("ShowerReco") << "caught exception \n" << e;
        xyz_vertex_fit[1] = 0;
        xyz_vertex_fit[2] = 0;
        xyz_vertex_fit[0] = 0;
      }

      // if collection is not best plane, project starting point from that
      if (bp1 != fNPlanes - 1 && bp2 != fNPlanes - 1) {
        geo::PlaneID const lastPlaneID{0, 0, fNPlanes - 1};
        auto pos = geom->Plane(lastPlaneID).GetBoxCenter();
        pos.SetY(xyz_vertex_fit[1]);
        pos.SetZ(xyz_vertex_fit[2]);
        auto const wirevertex = geom->NearestWireID(pos, lastPlaneID).Wire;

        double drifttick =
          (xyz_vertex_fit[0] / detProp.DriftVelocity(detProp.Efield(), detProp.Temperature())) *
          (1. / fTimeTick);
        fWire_vertex[fNPlanes - 1] = wirevertex; // wire coordinate of vertex for each plane
        fTime_vertex[fNPlanes - 1] =
          drifttick -
          (pos.X() / detProp.DriftVelocity(detProp.Efield(), detProp.Temperature())) *
            (1. / fTimeTick) +
          trigger_offset(clockData);
      }

      if (fabs(xphi) < 5.) {
        xtheta = gser.Get3DSpecialCaseTheta(
          bp1, bp2, fWire_last[bp1] - fWire_vertex[bp1], fWire_last[bp2] - fWire_vertex[bp2]);
      }

      // zero the arrays just to make sure
      for (unsigned int i = 0; i < fNAngles; ++i) {
        fTotChargeADC[i] = 0;
        fTotChargeMeV[i] = 0;
        fTotChargeMeV_MIPs[i] = 0;
        fNpoints_corr_ADC_2cm[i] = 0;
        fNpoints_corr_MeV_2cm[i] = 0;

        fRMS_2cm[i] = 0;
        fNpoints_2cm[i] = 0;

        fCorr_MeV_2cm[i] = 0;
        fCorr_Charge_2cm[i] = 0;

        fChargeADC_2cm[i] = 0; //Initial charge in ADC/cm for each plane angle calculation 4cm
        fChargeMeV_2cm[i] = 0; //initial charge in MeV/cm for each angle calculation first 4cm
      }

      // do loop and choose Collection. With ful calorimetry can do all.
      if (!(fabs(xphi) > 89 && fabs(xphi) < 91)) // do not calculate pitch for extreme angles
        LongTransEnergy(geom,
                        clockData,
                        detProp,
                        0,
                        hitlist_all[fNPlanes - 1]); // temporary only plane 2.

      //////create spacepoints, and direction cosines for Shower creation

      // make an art::PtrVector of the clusters
      art::PtrVector<recob::Cluster> prodvec;
      for (unsigned int i = 0; i < clusterListHandle->size(); ++i) {
        art::Ptr<recob::Cluster> prod(clusterListHandle, i);
        prodvec.push_back(prod);
      }

      //create a singleSpacePoint at vertex.
      std::vector<recob::SpacePoint> spcpts;

      //get direction cosines and set them for the shower
      // TBD determine which angle to use for the actual shower
      double fPhi = xphi;
      double fTheta = xtheta;

      TVector3 dcosVtx(std::cos(fPhi * TMath::Pi() / 180) * std::sin(fTheta * TMath::Pi() / 180),
                       std::cos(fTheta * TMath::Pi() / 180),
                       std::sin(fPhi * TMath::Pi() / 180) * std::sin(fTheta * TMath::Pi() / 180));
      /// \todo really need to determine the values of the arguments of the
      /// recob::Shower ctor
      // fill with bogus values for now
      TVector3 dcosVtxErr(util::kBogusD, util::kBogusD, util::kBogusD);
      recob::Shower singShower;
      singShower.set_direction(dcosVtx);
      singShower.set_direction_err(dcosVtxErr);

      Shower3DVector->push_back(singShower);
      // associate the shower with its clusters
      util::CreateAssn(evt, *Shower3DVector, prodvec, *cassn);

      // get the hits associated with each cluster and associate those with the shower
      for (size_t p = 0; p < prodvec.size(); ++p) {
        std::vector<art::Ptr<recob::Hit>> hits = fmh.at(p);
        util::CreateAssn(evt, *Shower3DVector, hits, *hassn);
      }

      geo::PlaneID planeID(0, 0, fNPlanes - 1);
      calorimetrycol->emplace_back(
        Kin_En, vdEdx, vdQdx, vresRange, deadwire, Trk_Length, fTrkPitchC, planeID);

      art::PtrVector<recob::Shower> ssvec;

      art::ProductID aid = evt.getProductID<std::vector<recob::Shower>>();
      art::Ptr<recob::Shower> aptr(aid, 0, evt.productGetter(aid));
      ssvec.push_back(aptr);

      util::CreateAssn(evt, *calorimetrycol, ssvec, *calassn);
      ftree_shwf->Fill();
    } // end loop on Vectors of "Associated clusters"

    evt.put(std::move(Shower3DVector));
    evt.put(std::move(cassn));
    evt.put(std::move(hassn));
    evt.put(std::move(calorimetrycol));
    evt.put(std::move(calassn));
  }

  //------------------------------------------------------------------------------
  void ShowerReco::LongTransEnergy(geo::GeometryCore const* geom,
                                   detinfo::DetectorClocksData const& clockData,
                                   detinfo::DetectorPropertiesData const& detProp,
                                   unsigned int set,
                                   std::vector<art::Ptr<recob::Hit>> hitlist)
  {
    // alogorithm for energy vs dx of the shower (roto-translation) COLLECTION
    // VIEW

    calo::CalorimetryAlg calalg(fCaloPSet);

    double totCnrg = 0,
           totCnrg_corr = 0; // tot enegry of the shower in collection

    double time;
    unsigned int wire = 0, plane = fNPlanes - 1;

    double mevav2cm = 0.;
    double sum = 0.;
    double npoints_calo = 0;

    int direction = -1;

    //override direction if phi (XZ angle) is less than 90 degrees
    if (fabs(xphi) < 90) direction = 1;

    //variables to check whether a hit is close to the shower axis.
    double ortdist, linedist;
    double wire_on_line, time_on_line;

    //get effective pitch using 3D angles
    util::GeometryUtilities const gser{*geom, clockData, detProp};
    double newpitch = gser.PitchInView(plane, xphi, xtheta);

    using lar::to_element;
    using ranges::views::transform;
    for (auto const& hit : hitlist | transform(to_element)) {
      time = hit.PeakTime();
      wire = hit.WireID().Wire;
      plane = hit.WireID().Plane;

      double dEdx_new;

      if (fUseArea) { dEdx_new = calalg.dEdx_AREA(clockData, detProp, hit, newpitch); }
      else // this will hopefully go away, once all of the calibration factors
           // are calculated.
      {
        dEdx_new = calalg.dEdx_AMP(clockData, detProp, hit, newpitch);
      }

      //calculate total energy.
      totCnrg_corr += dEdx_new;

      // calculate the wire,time coordinates of the hit projection on to the 2D shower axis
      gser.GetPointOnLine(slope[plane] / fWireTimetoCmCm,
                          fWire_vertex[plane],
                          fTime_vertex[plane],
                          wire,
                          time,
                          wire_on_line,
                          time_on_line);
      linedist =
        gser.Get2DDistance(wire_on_line, time_on_line, fWire_vertex[plane], fTime_vertex[plane]);
      ortdist = gser.Get2DDistance(wire_on_line, time_on_line, wire, time);

      //calculate the distance from the vertex using the effective pitch metric
      double wdist = (((double)wire - (double)fWire_vertex[plane]) * newpitch) *
                     direction; //wdist is always positive

      if ((wdist < fcalodEdxlength) && (wdist > 0.2)) {

        vdEdx.push_back(dEdx_new);
        vresRange.push_back(fabs(wdist));
        vdQdx.push_back(hit.PeakAmplitude() / newpitch);
        Trk_Length = wdist;
        fTrkPitchC = fNPitch[set][plane];
        Kin_En += dEdx_new * newpitch;
        npoints_calo++;
        sum += dEdx_new;

        if (wdist < fdEdxlength &&
            ((direction == 1 && wire > fWire_vertex[plane]) // take no hits before vertex
                                                            // (depending on direction)
             || (direction == -1 && wire < fWire_vertex[plane])) &&
            ortdist < 4.5 && linedist < fdEdxlength) {
          fChargeMeV_2cm[set] += dEdx_new;
          fNpoints_2cm[set]++;
        }

        // fill out for 4cm preshower

        fDistribChargeMeV[set].push_back(dEdx_new); // vector with the first De/Dx points
        fDistribChargeposition[set].push_back(
          wdist); // vector with the first De/Dx points' positions

      } // end inside range if statement

    } // end first loop on hits.

    auto const signalType =
      hitlist.empty() ? geo::kMysteryType : geom->SignalType(hitlist.front()->WireID());

    if (signalType == geo::kCollection) {
      fTotChargeADC[set] = totCnrg * newpitch;
      fTotChargeMeV[set] = totCnrg_corr * newpitch;
    }

    // calculate average dE/dx
    if (fNpoints_2cm[set] > 0) { mevav2cm = fChargeMeV_2cm[set] / fNpoints_2cm[set]; }

    // second loop to calculate RMS
    for (auto const& hit : hitlist | transform(to_element)) {
      time = hit.PeakTime();
      wire = hit.WireID().Wire;
      plane = hit.WireID().Plane;
      double dEdx = 0;

      if (fUseArea) { dEdx = calalg.dEdx_AREA(clockData, detProp, hit, newpitch); }
      else // this will hopefully go away, once all of the calibration factors
           // are calculated.
      {
        dEdx = calalg.dEdx_AMP(clockData, detProp, hit, newpitch);
      }

      gser.GetPointOnLine(slope[plane] / fWireTimetoCmCm,
                          fWire_vertex[plane],
                          fTime_vertex[plane],
                          wire,
                          time,
                          wire_on_line,
                          time_on_line);
      linedist =
        gser.Get2DDistance(wire_on_line, time_on_line, fWire_vertex[plane], fTime_vertex[plane]);
      ortdist = gser.Get2DDistance(wire_on_line, time_on_line, wire, time);

      double wdist = (((double)wire - (double)fWire_vertex[plane]) * newpitch) * direction;

      if ((wdist < fcalodEdxlength) && (wdist > 0.2)) {
        if (wdist < fdEdxlength &&
            ((direction == 1 && wire > fWire_vertex[plane]) ||
             (direction == -1 && wire < fWire_vertex[plane])) &&
            ortdist < 4.5 && linedist < fdEdxlength) {
          fRMS_2cm[set] += (dEdx - mevav2cm) * (dEdx - mevav2cm);
        }

      } // end if on correct hits.
    }   // end RMS_calculating loop.

    if (fNpoints_2cm[set] > 0) { fRMS_2cm[set] = TMath::Sqrt(fRMS_2cm[set] / fNpoints_2cm[set]); }

    /// third loop to get only points inside of 1RMS of value.

    for (auto const& hit : hitlist | transform(to_element)) {
      time = hit.PeakTime();
      wire = hit.WireID().Wire;
      plane = hit.WireID().Plane;

      double dEdx = 0;
      if (fUseArea) { dEdx = calalg.dEdx_AREA(clockData, detProp, hit, newpitch); }
      else // this will hopefully go away, once all of the calibration factors
           // are calculated.
      {
        dEdx = calalg.dEdx_AMP(clockData, detProp, hit, newpitch);
      }

      gser.GetPointOnLine(slope[plane] / fWireTimetoCmCm,
                          fWire_vertex[plane],
                          fTime_vertex[plane],
                          wire,
                          time,
                          wire_on_line,
                          time_on_line);
      linedist =
        gser.Get2DDistance(wire_on_line, time_on_line, fWire_vertex[plane], fTime_vertex[plane]);
      ortdist = gser.Get2DDistance(wire_on_line, time_on_line, wire, time);

      double wdist = (((double)wire - (double)fWire_vertex[plane]) * newpitch) * direction;

      if ((wdist < fcalodEdxlength) && (wdist > 0.2 &&
                                        ((direction == 1 && wire > fWire_vertex[plane]) ||
                                         (direction == -1 && wire < fWire_vertex[plane])) &&
                                        ortdist < 4.5 && linedist < fdEdxlength)) {
        if (wdist < fdEdxlength) {
          if (((dEdx > (mevav2cm - fRMS_2cm[set])) && (dEdx < (mevav2cm + fRMS_2cm[set]))) ||
              (newpitch > 0.3 * fdEdxlength)) {
            fCorr_MeV_2cm[set] += dEdx;
            fNpoints_corr_MeV_2cm[set]++;
          }

        } // end if on good hits
      }
    } // end of third loop on hits

    if (fNpoints_corr_MeV_2cm[set] > 0) {
      fCorr_MeV_2cm[set] /= fNpoints_corr_MeV_2cm[set];
      fChargeMeV_2cm_refined[set] = fCorr_MeV_2cm[set];
    }
  }

  //------------------------------------------------------------------------------------//
  void ShowerReco::GetVertexAndAnglesFromCluster(art::Ptr<recob::Cluster> clust, unsigned int plane)
  // Get shower vertex and slopes.
  {
    // convert to cm/cm units needed in the calculation
    angle[plane] = clust->StartAngle();
    slope[plane] = std::tan(clust->StartAngle());
    fWire_vertex[plane] = clust->StartWire();
    fTime_vertex[plane] = clust->StartTick();

    fWire_vertexError[plane] = clust->SigmaStartWire(); // wire coordinate of vertex for each plane
    fTime_vertexError[plane] = clust->SigmaStartTick(); // time coordinate of vertex for each plane

    fWire_last[plane] = clust->EndWire(); // wire coordinate of last point for each plane
    fTime_last[plane] = clust->EndTick();
  }

  DEFINE_ART_MODULE(ShowerReco)
}
