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
#include <vector>
#include <string>
#include <cmath> // std::tan() ...

// ### Framework includes ###
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// ### ROOT includes ###
#include "TMath.h"
#include "TTree.h"

// ### LArSoft includes ###

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"


namespace shwf {

  class ShowerReco : public art::EDProducer {
  public:
    explicit ShowerReco(fhicl::ParameterSet const& pset);

  private:

    void beginJob();
    void beginRun(art::Run& run);
    void produce(art::Event& evt);                       /**Actual routine that reconstruct the shower*/

     void   GetVertexAndAnglesFromCluster(art::Ptr< recob::Cluster > clust,unsigned int plane); // Get shower vertex and slopes.
  //  void GetVertexN(art::Event& evt);
    void   Get2DVariables(art::PtrVector < recob::Hit> hitlist);

    void   LongTransEnergy(unsigned int set, std::vector< art::Ptr<recob::Hit> > hitlist, bool isData=false); //Longtudinal



  void ClearandResizeVectors(unsigned int nPlanes);


  int fRun,fEvent,fSubRun;
 //  calo::CalorimetryAlg calorim;
//input labels:
  float slope[3];       // in cm, cm
  float angle[3];       // in radians


  std::string fClusterModuleLabel;

  float ftimetick; // time sample in us
  double xyz_vertex[3];


  double totCnrg,totCnrg_corr;
  double fMean_wire_pitch ;   // wire pitch in cm
  fhicl::ParameterSet fCaloPSet;

  std::vector<double> fRMS_2cm;
  std::vector<int> fNpoints_2cm;
  std::vector<double> fCorr_MeV_2cm;
  std::vector<double> fCorr_Charge_2cm;

  std::vector<int> fNpoints_corr_ADC_2cm;
  std::vector<int> fNpoints_corr_MeV_2cm;

  std::vector<double> fTotChargeADC;   //Total charge in ADC/cm for each plane
  std::vector<double> fTotChargeMeV;  //Total charge in MeV/cm for each plane
  std::vector<double> fTotChargeMeV_MIPs;  //Total charge in MeV/cm for each plane

  std::vector<double> fChargeADC_2cm;   //Initial charge in ADC/cm for each plane first 4cm
  std::vector<double> fChargeMeV_2cm;  //initial charge in MeV/cm for each plane first 4cm

  std::vector<double> fChargeMeV_2cm_refined;
  std::vector<double> fChargeMeV_2cm_axsum;

  std::vector<std::vector<double> > fDistribChargeADC;  //vector with the first De/Dx points ADC
  std::vector<std::vector<double> > fDistribChargeMeV;  //vector with the first De/Dx points converted energy
  std::vector<std::vector<double> > fDistribHalfChargeMeV;
  std::vector<std::vector<double> > fDistribChargeposition;  //vector with the first De/Dx points' positions

  std::vector<std::vector<double> > fSingleEvtAngle;  //vector with the first De/Dx points
  std::vector<std::vector<double> > fSingleEvtAngleVal;  //vector with the first De/Dx points


  std::vector<unsigned int> fWire_vertex;  // wire coordinate of vertex for each plane
  std::vector<double> fTime_vertex;  // time coordinate of vertex for each plane

  std::vector<double> fWire_vertexError;  // wire coordinate of vertex for each plane
  std::vector<double> fTime_vertexError;  // time coordinate of vertex for each plane

  std::vector<unsigned int> fWire_last;  // wire coordinate of last point for each plane
  std::vector<double> fTime_last;  // time coordinate of last point for each plane




  //reconstructed 3D position of shower start
  std::vector<double> xyz_vertex_fit;

  std::vector< std::vector<double> > fNPitch;   // double array, to use each plane for each set of angles

  //calorimetry variables
  float Kin_En;
  std::vector<float> vdEdx;
  std::vector<float> vresRange;
  std::vector<float> vdQdx;
  std::vector<float> deadwire; //residual range for dead wires
  float Trk_Length;
  float fTrkPitchC;
  float fdEdxlength;	   //distance that gets used to determine e/gamma separation
  float fcalodEdxlength;  // cutoff distance for hits saved to the calo object.
  bool fUseArea;

  double xphi,xtheta;   // new calculated angles.
  unsigned int fTPC;    //tpc type
  unsigned int fNPlanes; // number of planes
  unsigned int fNAngles;
  TTree* ftree_shwf;

 //conversion and useful constants
  double fWirePitch;    //wire pitch in cm
  double fTimeTick;
  double fDriftVelocity;
  double fWireTimetoCmCm;

  std::vector< int > fNhitsperplane;
  std::vector< double > fTotADCperplane;
  //services

  art::ServiceHandle<geo::Geometry const> geom;
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

//    //temporary
//     int mcpdg;
//     double mcenergy;
//     double mcphi;
//     double mctheta;
//     std::vector< unsigned int> mcwirevertex;  // wire coordinate of vertex for each plane
//     std::vector< double> mctimevertex;  // time coordinate of vertex for each plane




  }; // class ShowerReco



//------------------------------------------------------------------------------
ShowerReco::ShowerReco(fhicl::ParameterSet const& pset) :
    EDProducer{pset}
{
  fClusterModuleLabel = pset.get< std::string >("ClusterModuleLabel");
//  fVertexCLusterModuleLabel=pset.get<std::string > ("VertexClusterModuleLabel");
  fCaloPSet=pset.get< fhicl::ParameterSet >("CaloAlg");

  fdEdxlength= pset.get< double >("dEdxlength");   //distance that gets used to determine e/gamma separation
  fcalodEdxlength= pset.get< double >("calodEdxlength");  // cutoff distance for hits saved to the calo object.
  fUseArea= pset.get< bool >("UseArea");

  produces< std::vector<recob::Shower>                >();
  produces< art::Assns<recob::Shower, recob::Cluster> >();
  produces< art::Assns<recob::Shower, recob::Hit>     >();
  produces< std::vector<anab::Calorimetry>              >();
  produces< art::Assns<recob::Shower, anab::Calorimetry> >();
}

// struct SortByWire
//   {
//     bool operator() (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) const
//     {
//       return (h1->Wire()->RawDigit()->Channel() < h2->Wire()->RawDigit()->Channel());
//     }
//   };

// ***************** //
void ShowerReco::beginJob()
{



  /** Get Geometry*/
  art::ServiceHandle<geo::Geometry const> geo;

  detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  /// \todo the call to geo->Nplanes() assumes this is a single cryostat and single TPC detector
  /// \todo need to generalize to multiple cryostats and TPCs
  fNPlanes = geo->Nplanes();
  fMean_wire_pitch = geo->WirePitch(); //wire pitch in cm

  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService const> tfs;

 // auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  ftimetick=detprop->SamplingRate()/1000.;






  ftree_shwf =tfs->make<TTree>("ShowerReco","Results");/**All-knowing tree with reconstruction information*/

  // ftree_shwf->Branch("ftheta_Mean","std::vector<double>",&fTheta_Mean);
  // ftree_shwf->Branch("ftheta_RMS","std::vector<double>",&fTheta_RMS );

  ftree_shwf->Branch("run",&fRun,"run/I");
  ftree_shwf->Branch("subrun",&fSubRun,"subrun/I");
  ftree_shwf->Branch("event",&fEvent,"event/I");
  ftree_shwf->Branch("nplanes",&fNPlanes,"nplanes/I");
  ftree_shwf->Branch("nangles",&fNAngles,"nangles/I");

//   ftree_shwf->Branch("fthetaN","std::vector<double>",&fThetaN_ang);
//   ftree_shwf->Branch("fphiN","std::vector<double>",&fPhiN_ang);

  ftree_shwf->Branch("xtheta",&xtheta,"xtheta/D");
  ftree_shwf->Branch("xphi",&xphi,"xphi/D");

  ftree_shwf->Branch("ftotChargeADC","std::vector<double>",&fTotChargeADC);
  ftree_shwf->Branch("ftotChargeMeV","std::vector<double>",&fTotChargeMeV);
  ftree_shwf->Branch("fTotChargeMeV_MIPs","std::vector<double>",&fTotChargeMeV_MIPs);

  ftree_shwf->Branch("NPitch","std::vector< std::vector<double> >", &fNPitch);
  //  ftree_shwf->Branch("Pitch","std::vector<double>", &fPitch);

  // this should be temporary - until the omega is sorted out.
  ftree_shwf->Branch("RMS_2cm","std::vector<double>",&fRMS_2cm);
  ftree_shwf->Branch("Npoints_2cm","std::vector<int>",&fNpoints_2cm);
//    ftree_shwf->Branch("RMS_4cm","std::vector<double>",&fRMS_4cm);
//    ftree_shwf->Branch("Npoints_4cm","std::vector<int>",&fNpoints_4cm);

  ftree_shwf->Branch("ChargeADC_2cm","std::vector<double>",&fChargeADC_2cm);
  ftree_shwf->Branch("ChargeMeV_2cm","std::vector<double>",&fChargeMeV_2cm);
//    ftree_shwf->Branch("ChargeADC_4cm","std::vector<double>",&fChargeADC_4cm);
//    ftree_shwf->Branch("ChargeMeV_4cm","std::vector<double>",&fChargeMeV_4cm);

  ftree_shwf->Branch("ChargeMeV_2cm_refined","std::vector<double>",&fChargeMeV_2cm_refined);
//    ftree_shwf->Branch("ChargeMeV_4cm_refined","std::vector<double>",&fChargeMeV_4cm_refined);

  ftree_shwf->Branch("ChargeMeV_2cm_axsum","std::vector<double>",&fChargeMeV_2cm_axsum);
//    ftree_shwf->Branch("ChargeMeV_4cm_axsum","std::vector<double>",&fChargeMeV_4cm_axsum);


  ftree_shwf->Branch("fNhitsperplane","std::vector<int>",&fNhitsperplane);
  ftree_shwf->Branch("fTotADCperplane","std::vector<double>",&fTotADCperplane);



  ftree_shwf->Branch("ChargedistributionADC","std::vector<std::vector<double>>",&fDistribChargeADC);

  ftree_shwf->Branch("ChargedistributionMeV","std::vector<std::vector<double>>",&fDistribChargeMeV);
  ftree_shwf->Branch("DistribHalfChargeMeV","std::vector<std::vector<double>>",&fDistribHalfChargeMeV);
  ftree_shwf->Branch("ChargedistributionPosition","std::vector<std::vector<double>>",&fDistribChargeposition);
  ftree_shwf->Branch("xyz_vertex_fit","std::vector<double>", &xyz_vertex_fit);

}


void ShowerReco::beginRun(art::Run&)
{
 detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

 fWirePitch = geom->WirePitch(); //wire pitch in cm
 fTimeTick=detprop->SamplingRate()/1000.;
 fDriftVelocity=detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
 fWireTimetoCmCm=(fTimeTick*fDriftVelocity)/fWirePitch;

}




 void ShowerReco::ShowerReco::ClearandResizeVectors(unsigned int /*nPlanes*/) {
    //calculate factorial for number of angles
  int fact=1;
  for (unsigned int i = 1; i <= fNPlanes; ++i) fact *= i;

  fNAngles=fact/2;

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

  for(unsigned int ii=0;ii<fNAngles;ii++) {
     fNPitch[ii].resize(fNPlanes,-1);
   }

   fNpoints_corr_ADC_2cm.clear();
   fNpoints_corr_MeV_2cm.clear();

   fNpoints_corr_ADC_2cm.resize(fNAngles,-1);
   fNpoints_corr_MeV_2cm.resize(fNAngles,-1);
//   fNpoints_corr_ADC_4cm.resize(fNAngles,-1);
//   fNpoints_corr_MeV_4cm.resize(fNAngles,-1);



   for(unsigned int ii=0;ii<fNPlanes;ii++)
     {  fDistribChargeADC[ii].resize(0);  //vector with the first De/Dx points
        fDistribChargeMeV[ii].resize(0);  //vector with the first De/Dx points
        fDistribHalfChargeMeV[ii].resize(0);
        fDistribChargeposition[ii].resize(0);  //vector with the first De/Dx points' positions
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



    fWire_vertex.resize(fNAngles,-1);
    fTime_vertex.resize(fNAngles,-1);
    fWire_vertexError.resize(fNPlanes,-1);
    fTime_vertexError.resize(fNPlanes,-1);
    fWire_last.resize(fNAngles,-1);
    fTime_last.resize(fNAngles,-1);
    fTotChargeADC.resize(fNAngles,0);
    fTotChargeMeV.resize(fNAngles,0);
    fTotChargeMeV_MIPs.resize(fNAngles,0);
    fRMS_2cm.resize(fNAngles,0);
    fNpoints_2cm.resize(fNAngles,0);

    fCorr_MeV_2cm.clear();
    fCorr_Charge_2cm.clear();
    xyz_vertex_fit.clear();

    fChargeADC_2cm.clear();
    fChargeMeV_2cm.clear();
    fChargeMeV_2cm_refined.clear();
    fChargeMeV_2cm_refined.clear();
    fChargeMeV_2cm_axsum.clear();


    fCorr_MeV_2cm.resize(fNAngles,0);
    fCorr_Charge_2cm.resize(fNAngles,0);
    xyz_vertex_fit.resize(3);

    fChargeADC_2cm.resize(fNAngles,0);   //Initial charge in ADC/cm for each plane angle calculation 4cm
    fChargeMeV_2cm.resize(fNAngles,0);   //initial charge in MeV/cm for each angle calculation first 4cm
    fChargeMeV_2cm_refined.resize(0);
    fChargeMeV_2cm_refined.resize(fNAngles,0);;
    fChargeMeV_2cm_axsum.resize(fNAngles,0);;

    fNhitsperplane.resize(fNPlanes,-1);
    fTotADCperplane.resize(fNPlanes,-1);


    vdEdx.clear();
    vresRange.clear();
    vdQdx.clear();


 }





// ***************** //
void ShowerReco::produce(art::Event& evt)
{

  detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  util::GeometryUtilities gser;
  fNPlanes = geom->Nplanes();
  //fdriftvelocity = detprop->DriftVelocity(Efield_SI,Temperature);
   std::unique_ptr<std::vector<recob::Shower> > Shower3DVector(new std::vector<recob::Shower>);
   std::unique_ptr< art::Assns<recob::Shower, recob::Cluster> > cassn(new art::Assns<recob::Shower, recob::Cluster>);
   std::unique_ptr< art::Assns<recob::Shower, recob::Hit>     > hassn(new art::Assns<recob::Shower, recob::Hit>);
   std::unique_ptr< std::vector<anab::Calorimetry> > calorimetrycol(new std::vector<anab::Calorimetry>);
   std::unique_ptr< art::Assns< anab::Calorimetry,recob::Shower> > calassn(new art::Assns<anab::Calorimetry,recob::Shower>);


  /**Get Clusters*/




  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  art::Handle< std::vector<art::PtrVector < recob::Cluster> > > clusterAssociationHandle;
  evt.getByLabel(fClusterModuleLabel,clusterAssociationHandle);


  art::FindManyP<recob::Hit> fmh(clusterListHandle, evt,fClusterModuleLabel);

//   std::cout << " ---------- ShowerReco!!! -------------- " << std::endl;
  fRun = evt.id().run();
  fSubRun = evt.id().subRun();
  fEvent = evt.id().event();
//   unsigned int nCollections= clusterAssociationHandle->size();
//   std::vector < art::PtrVector<recob::Cluster> >::const_iterator clusterSet = clusterAssociationHandle->begin();
//    for(unsigned int iCol=0;iCol<nCollections;iCol++)
//   {
//    const art::PtrVector<recob::Cluster> pvcluster(*(clusterSet++));
//     auto it = pvcluster.begin();
//     int nClusts = pvcluster.size();
//     ClearandResizeVectors(nClusts);   // need to do this smarter (coutn at start and have an overall index?)
//     for(int iClust = 0; iClust < nClusts; ++iClust) {
//       const art::Ptr<recob::Cluster> pclust(*(it++));
//       auto pcoll { pclust };
//       art::FindManyP<recob::Hit> fs( pcoll, evt, fClusterModuleLabel);
//       std::vector< art::Ptr<recob::Hit> > hitlist = fs.at(0);
//       //std::cout << " hitlist size for coll " << iCol << " clust " << iClust << " " << hitlist.size() << std::endl;
//     }
//   }


  // find all the hits associated to all the clusters (once and for all);
  // the index of the query matches the index of the cluster in the collection
  // (conveniently carried around in its art pointer)
  art::FindManyP<recob::Hit> ClusterHits(clusterListHandle, evt, fClusterModuleLabel);

  std::vector < art::PtrVector<recob::Cluster> >::const_iterator clusterSet = clusterAssociationHandle->begin();
  // loop over vector of vectors (each size of NPlanes) and reconstruct showers from each of those
  for(size_t iClustSet = 0;iClustSet < clusterAssociationHandle->size(); iClustSet++){

    const art::PtrVector<recob::Cluster>  CurrentClusters=(*(clusterSet++));

    // do some error checking - i.e. are the clusters themselves present.
    if(clusterListHandle->size() < 2 || CurrentClusters.size() < 2) {
      //std::cout << "not enough clusters to reconstruct" << std::endl;
      //std::cout << "emergency filling tree @ run, evt," <<  fRun << " " << fEvent << std::endl;
      ftree_shwf->Fill();
      return;
    }
    //std::cout << " Cluster Set: " << iClustSet << " " << std::endl;

    ClearandResizeVectors( fNPlanes);

    std::vector< std::vector< art::Ptr<recob::Hit> > > hitlist_all;
    hitlist_all.resize(fNPlanes);

   // int nClusts = CurrentClusters.size();
    for(size_t iClust = 0; iClust < CurrentClusters.size(); iClust++){
      art::Ptr<recob::Cluster> const& pclust = CurrentClusters[iClust];
        //size_t ii=0;
	//std::cout << " clusterListHandle  " << clusterListHandle->size() << " fNPlanes " << fNPlanes << " "<< CurrentClusters.size() << std::endl;

      // get all the hits for this cluster;
      // pclust is a art::Ptr to the original cluster collection stored in the event;
      // its key corresponds to its index in the collection
      // (and therefore to the index in the query)
      std::vector< art::Ptr<recob::Hit> > const& hitlist = ClusterHits.at(pclust.key());
      //std::cout << " hitlist size for coll " << iClustSet << " clust " << iClust << " " << hitlist.size() << std::endl;

      unsigned int p(0); //c=channel, p=plane, w=wire

      //std::cout << " hitlist size " << hitlist.size() << std::endl;
      if(hitlist.size() == 0) continue;

      p=  (*hitlist.begin())->WireID().Plane;
      //art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
      //get vertex position and slope information to start with - ii is the posistion of the correct cluster:
      GetVertexAndAnglesFromCluster( pclust,p);


      double ADCcharge=0;
      //loop over cluster hits
      for(art::Ptr<recob::Hit> const& hit: hitlist){
        p=  hit->WireID().Plane;
        hitlist_all[p].push_back(hit);
        ADCcharge+= hit->PeakAmplitude();
      }
      fNhitsperplane[p]=hitlist_all[p].size();
      fTotADCperplane[p]=ADCcharge;
    //   unsigned int nCollections= clusterAssociationHandle->size();
//   std::vector < art::PtrVector<recob::Cluster> >::const_iterator clusterSet = clusterAssociationHandle->begin();
//    for(unsigned int iCol=0;iCol<nCollections;iCol++)
//   {
//    const art::PtrVector<recob::Cluster> pvcluster(*(clusterSet++));
//     auto it = pvcluster.begin();
//     int nClusts = pvcluster.size();
//     ClearandResizeVectors(nClusts);   // need to do this smarter (coutn at start and have an overall index?)
//     for(int iClust = 0; iClust < nClusts; ++iClust) {
//       const art::Ptr<recob::Cluster> pclust(*(it++));
//       auto pcoll { pclust };
//       art::FindManyP<recob::Hit> fs( pcoll, evt, fClusterModuleLabel);
//       std::vector< art::Ptr<recob::Hit> > hitlist = fs.at(0);
//       //std::cout << " hitlist size for coll " << iCol << " clust " << iClust << " " << hitlist.size() << std::endl;
//     }
//   }





   // Find the right cluster from the standard cluster list -> to get the hitlist associations.
//    for( ii = 0; ii < clusterListHandle->size(); ii++){
// 	art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
// 	if((*cl).ID() == (*CurrentClusters[iClust]).ID() )  //find the right cluster out of the list of associated clusters
// 	break;
//        }
//
    //get vertex position and slope information to start with - ii is the posistion of the correct cluster:
   // std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(ii);
  //  std::sort(hitlist.begin(), hitlist.end(), SortByWire());


  } // End loop on clusters.
  // Now I have the Hitlists and the relevent clusters parameters saved.


//   for(unsigned int i = 0; i < fNPlanes; ++i){
//     std::sort(hitlist_all[i].begin(), hitlist_all[i].end(),SortByWire());
//
//   }

  //find best set:
  unsigned int bp1 = 0,bp2 = 0;
  double minerror1=99999999,minerror2=9999999;
  for(unsigned int ii = 0; ii < fNPlanes; ++ii)
  {
    double locerror=fWire_vertexError[ii]*fWire_vertexError[ii]+ fTime_vertexError[ii]*fTime_vertexError[ii];  // time coordinate of vertex for each plane

    if(minerror1 >= locerror )  // the >= sign is to favorize collection
    {
     minerror1=locerror;
     bp1=ii;
    }
  }
  for(unsigned int ij = 0; ij < fNPlanes; ++ij)
  {
    double locerror=fWire_vertexError[ij]*fWire_vertexError[ij]+ fTime_vertexError[ij]*fTime_vertexError[ij];  // time coordinate of vertex for each plane

    if(minerror2 >= locerror && ij!=bp1 )
    {
     minerror2=locerror;
     bp2=ij;
    }
  }
  //bp1=0;
  //bp1=2;
  //std::cout << " best planes " << bp1 << " " << bp2 << std::endl;
for(unsigned int ij = 0; ij < fNPlanes; ++ij)
  {
   //std::cout << " wire distances: " << ij << " " << fabs((double)fWire_vertex[ij]-(double)fWire_last[ij]);

  }


 //std::cout << " angles in:  " << bp1 << " " << bp2 << " "   << slope[bp1]*TMath::Pi()/180. << " " << slope[bp2]*TMath::Pi()/180. << " " << slope[bp1] << " " << slope[bp2] << std::endl;
// gser.Get3DaxisN(bp1,bp2,slope[bp1]*TMath::Pi()/180.,slope[bp2]*TMath::Pi()/180.,xphi,xtheta);
 gser.Get3DaxisN(bp1,bp2,angle[bp1],angle[bp2],xphi,xtheta);

  ///////////////////////////////////////////////////////////
 const double origin[3] = {0.};
 std::vector <std::vector <  double > > position;
 // auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
 double fTimeTick=detprop->SamplingRate()/1000.;
 double fDriftVelocity=detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
 // get starting positions for all planes
 for(unsigned int xx=0;xx<fNPlanes;xx++){
      double pos1[3];
      geom->Plane(xx).LocalToWorld(origin, pos1);
      std::vector <double > pos2;
      pos2.push_back(pos1[0]);
      pos2.push_back(pos1[1]);
      pos2.push_back(pos1[2]);
      position.push_back(pos2);
    }
  // Assuming there is no problem ( and we found the best pair that comes close in time )
  // we try to get the Y and Z coordinates for the start of the shower.
    try{
	int chan1=geom->PlaneWireToChannel(bp1,fWire_vertex[bp1], 0);
	int chan2=geom->PlaneWireToChannel(bp2,fWire_vertex[bp2], 0);

	double y,z;
//	bool wires_cross = geom->ChannelsIntersect(chan1,chan2,y,z);
	geom->ChannelsIntersect(chan1,chan2,y,z);
 //       geom->ChannelsIntersect(chan1,chan2,y,z);

	xyz_vertex_fit[1]=y;
	xyz_vertex_fit[2]=z;
	xyz_vertex_fit[0]=(fTime_vertex[bp1]-detprop->TriggerOffset()) *fDriftVelocity*fTimeTick+position[0][0];


	//std::cout << ":::::: found x,y,z vertex " << wires_cross << " " << xyz_vertex_fit[0] << " " << y << " " << z << " " << wires_cross << std::endl;
    }
    catch(cet::exception const& e) {
      mf::LogWarning("ShowerReco") << "caught exception \n" << e;
      xyz_vertex_fit[1]=0;
      xyz_vertex_fit[2]=0;
      xyz_vertex_fit[0]=0;
     }


   // if collection is not best plane, project starting point from that
      if(bp1!=fNPlanes-1 && bp2!=fNPlanes-1)
      {
	double pos[3];
	unsigned int  wirevertex;

	geom->Plane(fNPlanes-1).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
	//std::cout << "plane X positionp " << 2 << " " << pos[0] << std::endl;

	pos[1]=xyz_vertex_fit[1];
	pos[2]=xyz_vertex_fit[2];
	wirevertex = geom->NearestWire(pos,fNPlanes-1);
	//geo->ChannelToWire(channel2,cs,t,p,wirevertex);   //\fixme!
        //wirevertex=  (*a)->WireID().Wire;




     double drifttick=(xyz_vertex_fit[0]/detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()))*(1./fTimeTick);
     fWire_vertex[fNPlanes-1]= wirevertex;  // wire coordinate of vertex for each plane
     fTime_vertex[fNPlanes-1] = drifttick-(pos[0]/detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()))*(1./fTimeTick)+detprop->TriggerOffset();


      }


   //  std::cout << "^^^^^^cross-check xphi and xtheta: " << xphi << " " << xtheta << std::endl;


      if(fabs(xphi) < 5. )
      { xtheta= gser.Get3DSpecialCaseTheta(bp1,bp2,fWire_last[bp1]-fWire_vertex[bp1], fWire_last[bp2]-fWire_vertex[bp2]);

	//std::cout << "xphi, xtheta1:,alt" << xphi << " " << xtheta  <<std::endl;
      }
    //}

// zero the arrays just to make sure
  for(unsigned int i = 0; i < fNAngles; ++i){
    fTotChargeADC[i]=0;
    fTotChargeMeV[i]=0;
    fTotChargeMeV_MIPs[i]=0;
    fNpoints_corr_ADC_2cm[i]=0;
    fNpoints_corr_MeV_2cm[i]=0;

    fRMS_2cm[i]=0;
    fNpoints_2cm[i]=0;

    fCorr_MeV_2cm[i]=0;
    fCorr_Charge_2cm[i]=0;

    fChargeADC_2cm[i]=0;  //Initial charge in ADC/cm for each plane angle calculation 4cm
    fChargeMeV_2cm[i]=0;  //initial charge in MeV/cm for each angle calculation first 4cm

  }


  // do loop and choose Collection. With ful calorimetry can do all.
  //for(unsigned int set=0;set<fNAngles;set++)
   if(!(fabs(xphi) >89 && fabs(xphi)<91)) // do not calculate pitch for extreme angles
    LongTransEnergy(0,hitlist_all[fNPlanes-1]); //temporary only plane 2.


   //////create spacepoints, and direction cosines for Shower creation

   //std::vector< recob::SpacePoint > 	spacepoints = std::vector<recob::SpacePoint>()



  // make an art::PtrVector of the clusters
  art::PtrVector<recob::Cluster> prodvec;
  for(unsigned int i = 0; i < clusterListHandle->size(); ++i){
    art::Ptr<recob::Cluster> prod(clusterListHandle, i);
    prodvec.push_back(prod);
  }

  //create a singleSpacePoint at vertex.
  std::vector< recob::SpacePoint > spcpts;

  //get direction cosines and set them for the shower
  // TBD determine which angle to use for the actual shower
  double fPhi=xphi;
  double fTheta=xtheta;

  TVector3 dcosVtx(TMath::Cos(fPhi*TMath::Pi()/180)*TMath::Sin(fTheta*TMath::Pi()/180),
		   TMath::Cos(fTheta*TMath::Pi()/180),
		   TMath::Sin(fPhi*TMath::Pi()/180)*TMath::Sin(fTheta*TMath::Pi()/180));
  /// \todo really need to determine the values of the arguments of the recob::Shower ctor
  // fill with bogus values for now
  TVector3 dcosVtxErr(util::kBogusD, util::kBogusD, util::kBogusD);
  //double maxTransWidth[2] = { util::kBogusD };
  //double distMaxWidth = util::kBogusD;
  //recob::Shower  singShower(dcosVtx, dcosVtxErr, maxTransWidth, distMaxWidth, 1);
  recob::Shower singShower;
  singShower.set_direction(dcosVtx);
  singShower.set_direction_err(dcosVtxErr);

  Shower3DVector->push_back(singShower);
  // associate the shower with its clusters
  util::CreateAssn(*this, evt, *Shower3DVector, prodvec, *cassn);

   // get the hits associated with each cluster and associate those with the shower
   for(size_t p = 0; p < prodvec.size(); ++p){
      std::vector< art::Ptr<recob::Hit> > hits = fmh.at(p);
      util::CreateAssn(*this, evt, *Shower3DVector, hits, *hassn);
   }


   geo::PlaneID planeID(0,0,fNPlanes-1);
   calorimetrycol->push_back(anab::Calorimetry(Kin_En,
					       vdEdx,
					       vdQdx,
					       vresRange,
					       deadwire,
					       Trk_Length,
					       fTrkPitchC,
					       planeID));

  art::PtrVector < recob::Shower >  ssvec;

    //for(unsigned int ip=0;ip<1;ip++)  {
        art::ProductID aid = evt.getProductID< std::vector < recob::Shower > >();
	art::Ptr< recob::Shower > aptr(aid, 0, evt.productGetter(aid));
	ssvec.push_back(aptr);
      //}


  //util::CreateAssn(*this, evt, *Shower3DVector, calorimetrycol, *calassn);
   util::CreateAssn(*this, evt, *calorimetrycol,ssvec,*calassn);
  /**Fill the output tree with all information */
  //std::cout << " filling tree @ run, evt," <<  fRun << " " << fEvent << std::endl;
  ftree_shwf->Fill();

  //for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane)
  //fh_theta[iplane]->Write(Form("fh_theta_%d_%d",iplane,evt.id().event()));
  // This needs work, clearly.
  //for(int p=0;p<2;p++)Shower3DVector->push_back(shower);

  } // end loop on Vectors of "Associated clusters"

  evt.put(std::move(Shower3DVector));
  evt.put(std::move(cassn));
  evt.put(std::move(hassn));
  evt.put(std::move(calorimetrycol));
  evt.put(std::move(calassn));



}



//------------------------------------------------------------------------------
void ShowerReco::LongTransEnergy(unsigned int set, std::vector < art::Ptr<recob::Hit> > hitlist, bool /*isData*/)
{
  // alogorithm for energy vs dx of the shower (roto-translation) COLLECTION VIEW
  // double  wire_cm, time_cm;
  // int loop_nrg = 0;

  calo::CalorimetryAlg calalg(fCaloPSet);
  util::GeometryUtilities gser;



  double totCnrg = 0,totCnrg_corr =0;//, totNewCnrg=0 ; // tot enegry of the shower in collection
//   art::ServiceHandle<geo::Geometry const> geom;
//   auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  double time;
  unsigned int wire=0,plane=fNPlanes-1;


  double mevav2cm=0.;
  double sum=0.;
  double npoints_calo=0;

  int direction=-1;

  //override direction if phi (XZ angle) is less than 90 degrees
  if(fabs(xphi)<90)
       direction=1;

  //variables to check whether a hit is close to the shower axis.
  double ortdist,linedist;
  double wire_on_line,time_on_line;

  //get effective pitch using 3D angles
  double newpitch=gser.PitchInView(plane,xphi,xtheta);



  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;

    wire=  theHit->WireID().Wire;
    plane=  theHit->WireID().Plane;

    double dEdx_new;
   // double dEdx_MIP;

    if(fUseArea)
    { dEdx_new = calalg.dEdx_AREA((*hitIter), newpitch );
     // dEdx_MIP = calalg.dEdx_AREA_forceMIP((*hitIter), newpitch );
    }
    else  //this will hopefully go away, once all of the calibration factors are calculated.
    {
      dEdx_new = calalg.dEdx_AMP((*hitIter), newpitch );
     // dEdx_MIP = calalg.dEdx_AMP_forceMIP((*hitIter), newpitch );
    }

    //calculate total energy.
    totCnrg_corr += dEdx_new;
   //totNewCnrg+=dEdx_MIP;

    // calculate the wire,time coordinates of the hit projection on to the 2D shower axis
    gser.GetPointOnLine(slope[plane]/fWireTimetoCmCm,fWire_vertex[plane],fTime_vertex[plane],wire,time,wire_on_line,time_on_line);
    linedist=gser.Get2DDistance(wire_on_line,time_on_line,fWire_vertex[plane],fTime_vertex[plane]);
    ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);



   //calculate the distance from the vertex using the effective pitch metric
    double wdist=(((double)wire-(double)fWire_vertex[plane])*newpitch)*direction;  //wdist is always positive


//   if( (fabs(wdist)<fcalodEdxlength)&&(fabs(wdist)>0.2)){
  if( (wdist<fcalodEdxlength)&&(wdist>0.2)){

    vdEdx.push_back(dEdx_new);
    vresRange.push_back(fabs(wdist));
    vdQdx.push_back((*hitIter)->PeakAmplitude()/newpitch);
    Trk_Length=wdist;
    fTrkPitchC=fNPitch[set][plane];
    Kin_En+=dEdx_new*newpitch;
    npoints_calo++;
    sum+=dEdx_new;



//std::cout << " CALORIMETRY:" << " Pitch " <<newpitch << " dist: " << wdist <<  " dE/dx: " << dEdx_new << "MeV/cm "  << " average: " <<  sum/npoints_calo  << "hit: wire, time " << wire << " " << time << " line,ort " << linedist << " " << ortdist<< " direction " << direction << std::endl;


     if(wdist<fdEdxlength
      && ((direction==1 && wire>fWire_vertex[plane])     //take no hits before vertex (depending on direction)
      || (direction==-1 && wire<fWire_vertex[plane])  )
	     && ortdist<4.5 && linedist < fdEdxlength ){
	      fChargeMeV_2cm[set]+= dEdx_new ;
	      fNpoints_2cm[set]++;
	     // std::cout << " CALORIMETRY:" << " Pitch " <<newpitch << " dist: " << wdist <<  " dE/dx: " << dEdx_new << "MeV/cm "  << " average: " <<  sum/npoints_calo  << "hit: wire, time " << wire << " " << time << " line,ort " << linedist << " " << ortdist<< " direction " << direction << std::endl;
	     }

        // fill out for 4cm preshower

        //fDistribChargeADC[set].push_back(ch_adc);  //vector with the first De/Dx points
	fDistribChargeMeV[set].push_back(dEdx_new);  //vector with the first De/Dx points
	//fDistribHalfChargeMeV[set].push_back(Bcorr_half);
        fDistribChargeposition[set].push_back(wdist);  //vector with the first De/Dx points' positions

     }//end inside range if statement

  }// end first loop on hits.

  auto const signalType
    = hitlist.empty()? geo::kMysteryType: geom->SignalType(hitlist.front()->WireID());

  if(signalType == geo::kCollection)
	{
      fTotChargeADC[set]=totCnrg*newpitch;
      fTotChargeMeV[set]=totCnrg_corr*newpitch;
      //fTotChargeMeV_MIPs[set]=totNewCnrg*newpitch;
	}


  //calculate average dE/dx
  if(fNpoints_2cm[set]>0)	{
    mevav2cm=fChargeMeV_2cm[set]/fNpoints_2cm[set];
   }
   //double RMS_ADC_2cm=0.;



  //second loop to calculate RMS
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end(); hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;
    wire=  theHit->WireID().Wire;
    plane=  theHit->WireID().Plane;
    double dEdx=0;

   if(fUseArea)
    { dEdx = calalg.dEdx_AREA((*hitIter), newpitch );
     }
    else  //this will hopefully go away, once all of the calibration factors are calculated.
    {
      dEdx = calalg.dEdx_AMP((*hitIter), newpitch );
     }



    gser.GetPointOnLine(slope[plane]/fWireTimetoCmCm,fWire_vertex[plane],fTime_vertex[plane],wire,time,wire_on_line,time_on_line);
    linedist=gser.Get2DDistance(wire_on_line,time_on_line,fWire_vertex[plane],fTime_vertex[plane]);
    ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);



    double wdist=(((double)wire-(double)fWire_vertex[plane])*newpitch)*direction;

 //    //std::cout << dEdx << " MeV, outside of if;; wd " << wdist << " ld " << linedist << " od " << ortdist << std::endl;

    if( (wdist<fcalodEdxlength)&&(wdist>0.2)){
      if(wdist<fdEdxlength
	  && ((direction==1 && wire>fWire_vertex[plane]) ||
	  (direction==-1 && wire<fWire_vertex[plane])  )
	     && ortdist<4.5 && linedist < fdEdxlength)
	{
//	  //std::cout << dEdx << " MeV " << std::endl;
	fRMS_2cm[set]+= (dEdx-mevav2cm)*(dEdx-mevav2cm);
	}


      }  //end if on correct hits.
  }  //end RMS_calculating loop.

  if(fNpoints_2cm[set]>0)
    {
    fRMS_2cm[set]=TMath::Sqrt(fRMS_2cm[set]/fNpoints_2cm[set]);
     }

  //std::cout << " average dE/dx: " << mevav2cm << " RMS::  " << fRMS_2cm[set] << " " << fNpoints_2cm[set] <<  std::endl;

  /// third loop to get only points inside of 1RMS of value.

  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end(); hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;
    wire=  theHit->WireID().Wire;
    plane=  theHit->WireID().Plane;

    double dEdx=0;
     if(fUseArea)
    { dEdx = calalg.dEdx_AREA((*hitIter), newpitch );
     }
    else  //this will hopefully go away, once all of the calibration factors are calculated.
    {
      dEdx = calalg.dEdx_AMP((*hitIter), newpitch );
     }

    gser.GetPointOnLine(slope[plane]/fWireTimetoCmCm,fWire_vertex[plane],fTime_vertex[plane],wire,time,wire_on_line,time_on_line);
    linedist=gser.Get2DDistance(wire_on_line,time_on_line,fWire_vertex[plane],fTime_vertex[plane]);
    ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);

    double wdist=(((double)wire-(double)fWire_vertex[plane])*newpitch)*direction;


    if( (wdist < fcalodEdxlength) && (wdist > 0.2
        && ((direction==1 && wire>fWire_vertex[plane]) ||
        (direction==-1 && wire<fWire_vertex[plane])  )
	     && ortdist<4.5 && linedist < fdEdxlength ))
      {
      if(wdist < fdEdxlength)
	  {
	  if( ((dEdx > (mevav2cm-fRMS_2cm[set]) )
	    && (dEdx < (mevav2cm+fRMS_2cm[set]) ))
		|| (newpitch > 0.3*fdEdxlength ) ) {
	    fCorr_MeV_2cm[set]+= dEdx;
	  fNpoints_corr_MeV_2cm[set]++;
	  }

      } // end if on good hits


    }
  } //end of third loop on hits

  if(fNpoints_corr_MeV_2cm[set]>0){
	//std::cout << " ++ NPoints 2cm, ADC and MeV "
	//			  << fNpoints_corr_MeV_2cm[set] << " "
	//			  << fNpoints_corr_ADC_2cm[set]
	//			  << " corrected average De/Dx, charge, MeV:  "
	//			  << fCorr_Charge_2cm[set]/fNpoints_corr_ADC_2cm[set]
	//			  << " " << fCorr_MeV_2cm[set]/fNpoints_corr_MeV_2cm[set] << std::endl;
	//fCorr_Charge_2cm[set]/=fNpoints_corr_ADC_2cm[set];
        fCorr_MeV_2cm[set]/=fNpoints_corr_MeV_2cm[set];
	fChargeMeV_2cm_refined[set]=fCorr_MeV_2cm[set];
  }


////std::cout << " total ENERGY, birks: " << fTotChargeMeV[set] << " MeV " << " assumeMIPs:  " << fTotChargeMeV_MIPs[set] << "MeV " <<  std::endl;
// std::cout << " total ENERGY, birks: " << fTotChargeMeV[set] << " MeV "  << " |average:  " << fChargeMeV_2cm_refined[set] <<   std::endl;
}



//------------------------------------------------------------------------------------//
void   ShowerReco::GetVertexAndAnglesFromCluster(art::Ptr< recob::Cluster > clust,unsigned int plane)
// Get shower vertex and slopes.
{
  //convert to cm/cm units needed in the calculation
//  slope[plane]=clust->dTdW();//*(ftimetick*fdriftvelocity)/fMean_wire_pitch;
  angle[plane]=clust->StartAngle();
  slope[plane]=std::tan(clust->StartAngle());
  fWire_vertex[plane]=clust->StartWire();
  fTime_vertex[plane]=clust->StartTick();

  fWire_vertexError[plane]=clust->SigmaStartWire();  // wire coordinate of vertex for each plane
  fTime_vertexError[plane]=clust->SigmaStartTick();  // time coordinate of vertex for each plane

  fWire_last[plane]=clust->EndWire();  // wire coordinate of last point for each plane
  fTime_last[plane]=clust->EndTick();

  ////////// insert detector offset

 // std::cout << "======= setting slope for view: " << plane
//			    << " " << slope[plane] << " " << fWire_vertex[plane]
//			    << " " << fTime_vertex[plane] << " " << std::endl;
	//		    <<  fWire_vertex[plane]+50<< " "
	//		    << fTime_vertex[plane] + slope[plane]*(fWire_vertex[plane]+50)<< std::endl;

}


// //------------------------------------------------------------------------------
// //to be moved to an _ana module
// void ShowerReco::GetVertexN(art::Event& evt){
//
//
//
//   double fTimeTick=detprop->SamplingRate()/1000.;
//
//   art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
//   evt.getByLabel("generator",mctruthListHandle);
//
//   art::PtrVector<simb::MCTruth> mclist;
//   for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
//     {
//       art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
//       mclist.push_back(mctparticle);
//     }
//
//
//   mf::LogVerbatim("ShowerAngleClusterAna")  << "%%%%%%% mc size size,  "<<mclist.size() <<    std::endl;
//
//
//   art::Ptr<simb::MCTruth> mc(mclist[0]);
//   simb::MCParticle neut(mc->GetParticle(0));
//
//   mcpdg=neut.PdgCode();
//   mcenergy=neut.E();
//   //std::cout << " ----------- mcenergy:: " << mcenergy << std::endl;
//   if (neut.P()){
//     double lep_dcosx_truth = neut.Px()/neut.P();
//     double lep_dcosy_truth = neut.Py()/neut.P();
//     double lep_dcosz_truth = neut.Pz()/neut.P();
//
//     mf::LogVerbatim("ShowerAngleClusterAna")  << "-----  cx,cy,cz " << lep_dcosx_truth << " " << lep_dcosy_truth << " " << lep_dcosz_truth << std::endl;
//
//
//     mcphi=  (lep_dcosx_truth == 0.0 && lep_dcosz_truth == 0.0) ? 0.0 : TMath::ATan2(lep_dcosx_truth,lep_dcosz_truth);
//     mctheta= (lep_dcosx_truth == 0.0 && lep_dcosy_truth == 0.0 && lep_dcosz_truth == 0.0) ? 0.0 : TMath::Pi()*0.5-TMath::ATan2(std::sqrt(lep_dcosx_truth*lep_dcosx_truth + lep_dcosz_truth*lep_dcosz_truth),lep_dcosy_truth);
//
//
//     mcphi=180*mcphi/TMath::Pi();
//     mctheta= 180*mctheta/TMath::Pi();
//     mf::LogVerbatim("ShowerAngleClusterAna")  << "-----  phi, theta " <<  mcphi << " " << mctheta << std::endl;
//
//   }
//
//   int npart=0;
//   //  while(&& npart < mc->NParticles() )
//   //     {
//   mf::LogVerbatim("ShowerAngleClusterAna")  << "%%%%%%%####### is PDG: "<< npart <<" " << neut.PdgCode() << std::endl;
//   //	neut=mc->GetParticle(npart++);
//
//   //   }
//
//   mf::LogVerbatim("ShowerAngleClusterAna")  << "%%%%%%%####### after loop is PDG: "<< npart <<" " << neut.PdgCode() << std::endl;
//   //if((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1){
//
//
//   xyz_vertex[0] =neut.Vx();
//   xyz_vertex[1] =neut.Vy();
//   xyz_vertex[2] =neut.Vz();
//
//   mf::LogVerbatim("ShowerAngleClusterAna") <<"neut.Vx()= "<<neut.Vx()<<" ,y= "<<neut.Vy()<<" ,z= "<<neut.Vz()<<std::endl;
//   //if(((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1))
//   //    break;
//
//
//   double drifttick=(xyz_vertex[0]/detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()))*(1./fTimeTick);
//
//   const double origin[3] = {0.};
//   for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
//     {
//       double pos[3];
//       geom->Plane(iplane).LocalToWorld(origin, pos);
//       //planex[p] = pos[0];
//       mf::LogVerbatim("ShowerAngleClusterAna")  << "plane X positionp " << iplane << " " << pos[0] << std::endl;
//
//       pos[1]=xyz_vertex[1];
//       pos[2]=xyz_vertex[2];
//       ///\todo: have to change to use cryostat and TPC in NearestWire too
//       unsigned int wirevertex = geom->NearestWire(pos,iplane);
//
//
//       mcwirevertex[iplane]=wirevertex;  // wire coordinate of vertex for each plane
//       mctimevertex[iplane]=drifttick-(pos[0]/detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()))*(1./fTimeTick)+detprop->TriggerOffset();  // time coordinate of vertex for each plane
//
//       //fWireVertex[p]=wirevertex;
//       //fTimeVertex[p]=drifttick-(pos[0]/detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()))*(1./fTimeTick)+60;
//       mf::LogVerbatim("ShowerAngleClusterAna") << "wirevertex= "<< wirevertex
// 					    << " timevertex " << mctimevertex[iplane]
// 					    << " correction "
// 					    << (pos[0]/detprop->DriftVelocity(detprop->Efield(),
// 									   detprop->Temperature()))*(1./fTimeTick)
// 					    << " " << pos[0];
//
//
//     }
//
//
//
//   return (void)0;
// }







  DEFINE_ART_MODULE(ShowerReco)


}
