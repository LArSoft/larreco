//
// Name: FeatureTracker.h
//
// Purpose:  This module takes features found in 2D and uses them
//            to produce seeds for 3D tracking.
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "TVector3.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Wire.h"

#include "larreco/RecoAlg/CornerFinderAlg.h"
#include "larreco/RecoAlg/SpacePointAlg.h"

namespace trkf {

  class FeatureTracker : public art::EDProducer {
  public:
    explicit FeatureTracker(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    // Fcl Attributes.

    std::string fTruthModuleLabel;
    std::string fHitModuleLabel;
    std::string fCalDataModuleLabel;

    void GetProjectedEnds(TVector3 xyz, std::vector<double>& uvw, std::vector<double>& t, int tpc=0, int cryo=0);

    std::map<int, std::vector<double> > ExtractEndPointTimes(std::vector<recob::EndPoint2D> EndPoints);

    std::vector<recob::SpacePoint> Get3DFeaturePoints(std::vector<recob::EndPoint2D> EndPoints, art::PtrVector<recob::Hit> Hits);

    std::vector<recob::Seed> GetValidLines(std::vector<recob::SpacePoint>& sps,
					   bool ApplyFilter = true);

    void RemoveDuplicatePaths(std::vector<recob::Seed>& Seeds,
			      std::vector<int>& ConnectionPoint1,
			      std::vector<int>& ConnectionPoint2);

    recob::Seed ExtendSeed(recob::Seed TheSeed);


    bool CheckSeedLineInt(recob::Seed& TheSeed);

    trkf::SpacePointAlg       fSP;
    corner::CornerFinderAlg  fCorner;


    double  fLineIntThreshold;
    double  fLineIntFraction;

    std::map<int, std::vector<double> > fEndPointTimes;
    art::ServiceHandle<geo::Geometry const> fGeometryHandle;
  };
}

#include "art/Framework/Core/ModuleMacros.h"

namespace trkf {
  DEFINE_ART_MODULE(FeatureTracker)
}

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"



namespace trkf {

  FeatureTracker::FeatureTracker(const fhicl::ParameterSet& pset):
    EDProducer{pset},
    fSP(pset.get<fhicl::ParameterSet>("SpacepointPset")),
    fCorner(pset.get<fhicl::ParameterSet>("CornerPset"))
  {
    fHitModuleLabel    = pset.get<std::string>("HitModuleLabel");
    fLineIntFraction   = pset.get<double>("LineIntFraction");
    fLineIntThreshold  = pset.get<double>("LineIntThreshold");
    fhicl::ParameterSet CornerPset = pset.get<fhicl::ParameterSet>("CornerPset");
    fCalDataModuleLabel = CornerPset.get<std::string>("CalDataModuleLabel");

    produces< std::vector<recob::Seed> >();
  }

  void FeatureTracker::produce(art::Event& evt)
  {

    // Extract hits PtrVector from event
    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(fHitModuleLabel, hith);

    art::PtrVector<recob::Hit> hitvec;
    for(unsigned int i = 0; i < hith->size(); ++i){
      art::Ptr<recob::Hit> prod(hith, i);
      hitvec.push_back(prod);
    }

    //We need to grab out the wires.
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel(fCalDataModuleLabel,wireHandle);
    std::vector<recob::Wire> const& wireVec(*wireHandle);

    //First, have it process the wires.
    fCorner.GrabWires(wireVec,*fGeometryHandle);

    std::vector<recob::EndPoint2D> EndPoints;
    fCorner.get_feature_points(EndPoints,*fGeometryHandle);

    fEndPointTimes = ExtractEndPointTimes(EndPoints);

    std::vector<recob::SpacePoint> sps = Get3DFeaturePoints(EndPoints, hitvec);

    std::vector<recob::Seed> SeedsToStore = GetValidLines( sps, true );

    std::unique_ptr< std::vector<recob::Seed > > seeds ( new std::vector<recob::Seed>);

    for(size_t i=0; i!=SeedsToStore.size(); ++i)
      seeds->push_back(SeedsToStore.at(i));

    evt.put(std::move(seeds));
  }




  //---------------------------------------------------------------------

  std::map<int, std::vector<double> > FeatureTracker::ExtractEndPointTimes(std::vector<recob::EndPoint2D> EndPoints)
  {
    std::map<int, std::vector<double> > EndPointTimesInPlane;
    for(size_t i=0; i!=EndPoints.size(); ++i)
      {
	EndPointTimesInPlane[EndPoints.at(i).View()].push_back(EndPoints.at(i).DriftTime());
      }

    for(std::map<int, std::vector<double> >::iterator itEpTime = EndPointTimesInPlane.begin();
	itEpTime != EndPointTimesInPlane.end(); ++itEpTime)
      {
	std::sort(itEpTime->second.begin(), itEpTime->second.end());
      }
    return EndPointTimesInPlane;
  }



  //---------------------------------------------------------------------

  std::vector<recob::Seed> FeatureTracker::GetValidLines(std::vector<recob::SpacePoint>& spts, bool ApplyFilter)
  {
    std::vector<recob::Seed> SeedsToReturn;

    std::vector<int> ConnectionPoint1;
    std::vector<int> ConnectionPoint2;
    std::map<int, std::vector<int> > SeedConnections;

    for(size_t i=0; i!=spts.size(); ++i)
      {
        for(size_t j=0; j!=i; ++j)
          {


            TVector3 xyz_i;
            TVector3 xyz_j;

	    std::vector<double> t_i, t_j;

	    std::vector<double> uvw_i;
	    std::vector<double> uvw_j;

            for(size_t d=0; d!=3; ++d)
              {
                xyz_i[d] = spts.at(i).XYZ()[d];
                xyz_j[d] = spts.at(j).XYZ()[d];
              }


            GetProjectedEnds(xyz_i, uvw_i, t_i, 0, 0);
            GetProjectedEnds(xyz_j, uvw_j, t_j, 0, 0);

	    bool ThisLineGood=true;

	    for(size_t p=0; p!=uvw_i.size(); ++p)
	      {
		TH2F const& RawHist = fCorner.GetWireDataHist(p);

		double lineint =
		  fCorner.line_integral(RawHist,
					uvw_i.at(p), t_i.at(p),
					uvw_j.at(p), t_j.at(p),
					fLineIntThreshold);

		if(lineint < fLineIntFraction)
		  {

		    ThisLineGood=false;
		  }
	      }
	    if(ThisLineGood)
	      {
		double Err[3];
		double Pos[3];
		double Dir[3];

		for(size_t d=0; d!=3; ++d)
		  {
		    Pos[d] = 0.5*(xyz_i[d] + xyz_j[d]);
		    Dir[d] = 0.5*(xyz_i[d] - xyz_j[d]);
		    Err[d] = 0;
		  }

		ConnectionPoint1.push_back(i);
		ConnectionPoint2.push_back(j);

		SeedsToReturn.push_back(recob::Seed(Pos,Dir,Err,Err));
	      }

	  }

      }

    if(ApplyFilter)
      {
	RemoveDuplicatePaths(SeedsToReturn, ConnectionPoint1, ConnectionPoint2);
	mf::LogInfo("FeatureTracker") <<"Seeds after filter " << SeedsToReturn.size()<<" seeds"<<std::endl;
      }

    return SeedsToReturn;
  }

  //--------------------------------------------------

  void FeatureTracker::RemoveDuplicatePaths(std::vector<recob::Seed>& Seeds,
				     std::vector<int>& ConnectionPoint1,
				     std::vector<int>& ConnectionPoint2)
  {

    std::map<int, bool> MarkedForRemoval;

    std::map<int, std::vector<int> > SeedsSharingPoint;
    for(size_t i=0; i!=Seeds.size(); ++i)
      {
	SeedsSharingPoint[ConnectionPoint1[i] ].push_back(i);
	SeedsSharingPoint[ConnectionPoint2[i] ].push_back(i);
      }

    for(size_t s=0; s!=Seeds.size(); ++s)
      {

	int StartPt = ConnectionPoint1.at(s);
	int EndPt   = ConnectionPoint2.at(s);
	int MidPt   = -1;

	for(size_t SeedsWithThisStart =0; SeedsWithThisStart!=SeedsSharingPoint[StartPt].size(); SeedsWithThisStart++)
	  {
	    int i = SeedsSharingPoint[StartPt].at(SeedsWithThisStart);
	    if(ConnectionPoint1.at(i)==StartPt)
	      MidPt = ConnectionPoint2.at(i);
	    else if(ConnectionPoint2.at(i)==StartPt)
	      MidPt = ConnectionPoint1.at(i);

	    for(size_t SeedsWithThisMid =0; SeedsWithThisMid!=SeedsSharingPoint[MidPt].size(); SeedsWithThisMid++)
	      {
		int j =  SeedsSharingPoint[MidPt].at(SeedsWithThisMid);
		if((ConnectionPoint1.at(j)==EndPt)||(ConnectionPoint2.at(j)==EndPt))
		  {

		    double Lengthi = Seeds.at(i).GetLength();
		    double Lengthj = Seeds.at(j).GetLength();
		    double Lengths = Seeds.at(s).GetLength();

		    if((Lengths>Lengthi)&&(Lengths>Lengthj))
		      {
			MarkedForRemoval[i]=true;
			MarkedForRemoval[j]=true;
		      }

		    if((Lengthi>Lengths)&&(Lengthi>Lengthj))
		      {
			MarkedForRemoval[s]=true;
			MarkedForRemoval[j]=true;
		      }
		    if((Lengthj>Lengthi)&&(Lengthj>Lengths))
		      {
			MarkedForRemoval[s]=true;
			MarkedForRemoval[i]=true;
		      }
		  }
	      }
	  }
      }
    for(std::map<int,bool>::reverse_iterator itrem = MarkedForRemoval.rbegin();
	itrem != MarkedForRemoval.rend(); ++itrem)
      {
	if(itrem->second==true)
	  {
	    Seeds.erase(             Seeds.begin()            + itrem->first);
	    ConnectionPoint1.erase(  ConnectionPoint1.begin() + itrem->first);
	    ConnectionPoint2.erase(  ConnectionPoint2.begin() + itrem->first);
	  }
      }

  }





  //---------------------------------------------------------------------

  void FeatureTracker::GetProjectedEnds(TVector3 xyz, std::vector<double>& uvw, std::vector<double>&t, int tpc, int cryo)
  {

    const detinfo::DetectorProperties* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<geo::Geometry const>              geo;

    int NPlanes=geo->Cryostat(cryo).TPC(tpc).Nplanes();

    uvw.resize(NPlanes);
    t.resize(NPlanes);

    for(int plane = 0; plane != NPlanes; plane++)
      {
	uvw[plane] = geo->NearestWire(xyz, plane, tpc, cryo);
	t[plane] = det->ConvertXToTicks(xyz[0], plane, tpc, cryo);
      }

  }

  //----------------------------------------------------------------------

  std::vector<recob::SpacePoint> FeatureTracker::Get3DFeaturePoints(std::vector<recob::EndPoint2D> EndPoints, art::PtrVector<recob::Hit> Hits)
    {

      art::PtrVector<recob::Hit> HitsForEndPoints;

      // Loop through the hits looking for the ones which match corners
      for(std::vector<recob::EndPoint2D>::const_iterator itEP=EndPoints.begin(); itEP!=EndPoints.end(); ++itEP)
	{
	  int EPMatchCount=0;

	  for(art::PtrVector<recob::Hit>::const_iterator itHit = Hits.begin(); itHit != Hits.end(); ++itHit)
	    {


	      if(
		 (itEP->View() == (*itHit)->View()) &&
		 (itEP->WireID().Wire == (*itHit)->WireID().Wire))
		{
		  HitsForEndPoints.push_back(*itHit);
		  EPMatchCount++;

		}
	    }

	}
      std::vector<recob::SpacePoint> spts;
      fSP.makeSpacePoints(HitsForEndPoints, spts);

      // This is temporary for debugging purposes
      const detinfo::DetectorProperties* det = lar::providerFrom<detinfo::DetectorPropertiesService>();

      for(size_t i=0; i!=spts.size(); ++i)
	{
	  for(size_t p=0; p!=3; ++p)
	    {
	      double Closest =100000;
	      double spt_t = det->ConvertXToTicks(spts.at(i).XYZ()[0], p, 0, 0);
	      for(size_t epTime=0; epTime!=fEndPointTimes[p].size(); ++epTime)
		{
		  if( fabs(fEndPointTimes[p].at(epTime) - spt_t)<Closest)
		    {
		      Closest =fabs(fEndPointTimes[p].at(epTime)-spt_t);
		    }
		}
	    }

	}
      return spts;
    }




  //----------------------------------------------------------------------

  recob::Seed FeatureTracker::ExtendSeed(recob::Seed TheSeed)
  {

    double ExtendFrac=1.1;

    // Extend Forward
    bool SuccessfulExtend = true;
    while(SuccessfulExtend)
      {
	recob::Seed NewSeed = TheSeed;
	double Dir[3], Err[3];
	double Pt[3];
	NewSeed.GetDirection( Dir, Err );
	NewSeed.GetPoint(     Pt,  Err );
	for(size_t i=0; i!=3; ++i)
	  {
	    Dir[i] *= ExtendFrac;
	    Pt[i]  += Dir[i] * ( ExtendFrac - 1. ) * 0.5;
	    NewSeed.SetPoint(     Pt,  Err );
	    NewSeed.SetDirection( Dir, Err );
	  }
	SuccessfulExtend  = CheckSeedLineInt(NewSeed);
	if(SuccessfulExtend) TheSeed = NewSeed;
      }

    // Extend backward
    SuccessfulExtend = true;
    while(SuccessfulExtend)
      {
	recob::Seed NewSeed = TheSeed;
	double Dir[3], Err[3];
	double Pt[3];
	NewSeed.GetDirection( Dir, Err );
	NewSeed.GetPoint(     Pt,  Err );
	for(size_t i=0; i!=3; ++i)
	  {
	    Dir[i] *= ExtendFrac;
	    Pt[i]  -= Dir[i] * ( ExtendFrac - 1. ) * 0.5;
	    NewSeed.SetPoint(     Pt,  Err );
	    NewSeed.SetDirection( Dir, Err );
	  }
	SuccessfulExtend  = CheckSeedLineInt(NewSeed);
	if(SuccessfulExtend) TheSeed = NewSeed;
      }
    return TheSeed;
  }

  //---------------------------------------------------

  bool FeatureTracker::CheckSeedLineInt(recob::Seed& TheSeed)
  {

    double Pt[3],Dir[3],Err[3];
    TheSeed.GetPoint(Pt,Err);
    TheSeed.GetDirection(Dir,Err);

    TVector3 endi, endj;

    for(size_t i=0; i!=3; ++i)
      {
	endi[i] = Pt[i]+Dir[i];
	endj[i] = Pt[i]-Dir[i];
      }

    std::vector<double> t_i, t_j;

    std::vector<double> uvw_i;
    std::vector<double> uvw_j;

    GetProjectedEnds(endi, uvw_i, t_i, 0, 0);
    GetProjectedEnds(endj, uvw_j, t_j, 0, 0);

    for(size_t p=0; p!=uvw_i.size(); ++p)
      {
	TH2F const& RawHist = fCorner.GetWireDataHist(p);

	double lineint =
	  fCorner.line_integral(RawHist,
				uvw_i.at(p), t_i.at(p),
				uvw_j.at(p), t_j.at(p),
				fLineIntThreshold);
	if(lineint < fLineIntFraction)
	  {
	    return false;
	  }
      }
    return true;
  }

}
