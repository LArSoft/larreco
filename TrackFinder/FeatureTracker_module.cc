#include "art/Persistency/Common/PtrVector.h"

#ifndef FEATURETRACKER_H
#define FEATURETRACKER_H

//
// Name: FeatureTracker.h
//
// Purpose:  This module takes features found in 2D and uses them
//            to produce seeds for 3D tracking.
// 
//
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SeedFinderAlgorithm.h"
#include "TVector3.h"
#include "RecoBase/EndPoint2D.h"

#include "RecoAlg/CornerFinderAlg.h"
#include "RecoAlg/SpacePointAlg.h"
#include "RecoObjects/BezierTrack.h"

namespace recob
{
  class Seed;
  class Hit;
}


namespace trkf {

  class BezierTrack;
  class BezierTrackerAlgorithm;
 
  class FeatureTracker : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit FeatureTracker(fhicl::ParameterSet const& pset);
    virtual ~FeatureTracker();

    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    
    

  private:

    // Fcl Attributes.

    std::string fTruthModuleLabel;
    std::string fHitModuleLabel;

    void GetProjectedEnds(TVector3 xyz, std::vector<double>& uvw, std::vector<double>& t, int tpc=0, int cryo=0);
    
    std::map<int, std::vector<double> > ExtractEndPointTimes(std::vector<recob::EndPoint2D> EndPoints);

    std::vector<recob::SpacePoint> Get3DFeaturePoints(std::vector<recob::EndPoint2D> EndPoints, art::PtrVector<recob::Hit> Hits);
    
    std::vector<recob::Seed> GetValidLines(std::vector<recob::SpacePoint>& sps,
					   bool ApplyFilter = true);
    
    void RemoveDuplicatePaths(std::vector<recob::Seed>& Seeds,
			      std::vector<int>& ConnectionPoint1,
			      std::vector<int>& ConnectionPoint2);

    recob::Seed ExtendSeed(recob::Seed TheSeed);


    std::map<int, std::map<int, double> >  GetConnectionMap(std::vector<recob::Seed>& Seeds, double ADCThresh, double FracThresh);
    
    std::vector<trkf::BezierTrack> GenerateBezierTracks(std::map<int,std::map<int,double> > , std::vector<recob::Seed>);

    bool CheckSeedLineInt(recob::Seed& TheSeed);

    trkf::SpacePointAlg       fSP;
    cluster::CornerFinderAlg  fCorner;

    
    double  fLineIntThreshold;
    double  fLineIntFraction;

    std::map<int, std::vector<double> > fEndPointTimes;
  };
}

#endif 

#include "art/Framework/Core/ModuleMacros.h" 


namespace trkf {
  DEFINE_ART_MODULE(FeatureTracker)
}

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"



namespace trkf {

  FeatureTracker::FeatureTracker(const fhicl::ParameterSet& pset):
    fSP(pset.get<fhicl::ParameterSet>("SpacepointPset")),
    fCorner(pset.get<fhicl::ParameterSet>("CornerPset"))

  {
    reconfigure(pset);
    produces< std::vector<recob::Seed> >();
    
  }

  FeatureTracker::~FeatureTracker()
  {
  }

  void FeatureTracker::reconfigure(fhicl::ParameterSet const& pset)
  {
    fHitModuleLabel    = pset.get<std::string>("HitModuleLabel");
    fLineIntFraction   = pset.get<double>("LineIntFraction");
    fLineIntThreshold  = pset.get<double>("LineIntThreshold");
    
  }

  void FeatureTracker::beginJob()
  {}


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
    
    
    fCorner.TakeInRaw(evt);
   
    std::vector<recob::EndPoint2D> EndPoints;
    fCorner.get_feature_points(EndPoints);
    
    fEndPointTimes = ExtractEndPointTimes(EndPoints);
    
    std::vector<recob::SpacePoint> sps = Get3DFeaturePoints(EndPoints, hitvec);    
       
    std::vector<recob::Seed> SeedsToStore = GetValidLines( sps, true );    
    
    std::map<int, std::map<int, double> > ConnMap = GetConnectionMap(SeedsToStore, 3, 0.90);
    
    /*    for(size_t i=0; i!=SeedsToStore.size(); ++i)
      {
	SeedsToStore[i] = ExtendSeed(SeedsToStore.at(i));
      }
    

    ConnMap = GetConnectionMap(SeedsToStore, 3, 0.90);
    */   
 
    std::vector<trkf::BezierTrack> BTracks = GenerateBezierTracks(ConnMap, SeedsToStore);

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
		TH2F * RawHist = fCorner.GetWireDataHist(p);
		
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
    
    art::ServiceHandle<util::DetectorProperties>   det;
    art::ServiceHandle<geo::Geometry>              geo;

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
    art::ServiceHandle<util::DetectorProperties>   det;
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
	TH2F * RawHist = fCorner.GetWireDataHist(p);
	
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
  

  //----------------------------------------------------------------------
  
  std::map<int, std::map<int, double> > FeatureTracker::GetConnectionMap(std::vector<recob::Seed>& Seeds, double ADCThresh, double FracThresh)
  {

    art::ServiceHandle<geo::Geometry> geom;
    std::vector<TVector3> WireDirs;
    
    double EndThresh = 0.5;

    std::map<int,bool> RedundantSeeds;

    std::map<int, std::map<int, double> > ConnectionMap;
   
    for(size_t i=0; i!=Seeds.size(); ++i)
      {
	for(size_t j=0; j!=i; ++j)
	  {
	    std::vector<recob::Seed > SeedsVec;
	    SeedsVec.push_back(Seeds.at(i));
	    SeedsVec.push_back(Seeds.at(j));
	    
	    trkf::BezierTrack TestTrack(SeedsVec);
	    
	    std::vector<float> LineScore = fCorner.line_integrals(TestTrack, 100, ADCThresh);
	    
	    bool HasConnection=true;
	    
	    for(size_t view =0; view!=LineScore.size(); ++view)
	      {
		if(LineScore.at(view)<FracThresh)
		  {
		    HasConnection=false;
		  }	      
	      }
	    double Seed1Pt[3], Seed2Pt[3], Seed1Dir[3], Seed2Dir[3], Err[3];
	    
	    Seeds.at(i).GetPoint(Seed1Pt,Err);
	    Seeds.at(j).GetPoint(Seed2Pt,Err);
	    Seeds.at(i).GetDirection(Seed1Dir,Err);
	    Seeds.at(j).GetDirection(Seed2Dir,Err);
	    
	    TVector3 Seed1End1, Seed1End2, Seed2End1, Seed2End2;

	    for(size_t index=0; index!=3; ++index)
	      {
		Seed1End1[index]=Seed1Pt[index]+Seed1Dir[index];
		Seed1End2[index]=Seed1Pt[index]-Seed1Dir[index];
		Seed2End1[index]=Seed2Pt[index]+Seed2Dir[index];
		Seed2End2[index]=Seed2Pt[index]-Seed2Dir[index];		
	      }
	    
	    if( ( (Seed1End1-Seed2End1).Mag() < EndThresh)
		||( (Seed1End2-Seed2End1).Mag() < EndThresh)
		||( (Seed1End1-Seed2End2).Mag() < EndThresh)
		||( (Seed1End2-Seed2End2).Mag() < EndThresh))
	      HasConnection=true;
	    
	    if(HasConnection)
	      {
		
		double Pt[3]; double Dir[3]; double Err[3];
		TVector3 End1; TVector3 End2;
		
		ConnectionMap[i][j] = ConnectionMap[j][i] = TestTrack.GetLength();  
		for(size_t isd=0; isd!=Seeds.size(); ++isd)
		  {
		    if((isd!=i)&&(isd!=j)&&(RedundantSeeds[i]!=true)&&(RedundantSeeds[j]!=true)&&(RedundantSeeds[isd]!=true))
		      {
			Seeds.at(isd).GetPoint(Pt,Err);
			Seeds.at(isd).GetDirection(Dir,Err);
			TVector3 End1;
			TVector3 End2;
			for(size_t index=0; index!=3; ++index)
			  {
			    End1[index]=Pt[index]+Dir[index];
			    End2[index]=Pt[index]-Dir[index];
			  }
			double s1, s2, d1,d2;
			TestTrack.GetClosestApproach(End1,s1,d1);
			TestTrack.GetClosestApproach(End2,s2,d2);
			
			//			std::cout<<"in 3D :  " << d1 << ", " <<d2<<std::endl;
			if((d1<1.)&&(d2<1.))
			  {
			    std::cout<<" meets 3D throw condition"<<std::endl;
			  }
			
			bool NoFails=true;
			for(size_t p=0; p!=3; ++p)
			  {
			    int t=0, c=0;
			    uint32_t wire1 = geom->NearestWire(End1, p, t, c);
			    uint32_t wire2 = geom->NearestWire(End2, p, t, c);
			    double dp1, sp1, dp2, sp2;
			    TestTrack.GetClosestApproach(wire1, p, t, c, End1[0], sp1,dp1);
			    TestTrack.GetClosestApproach(wire2, p, t, c, End2[0], sp2,dp2);
			    //  std::cout<<p<<": [ "<< dp1 <<", "<<dp2<<" ]     " ;
			    if((dp1>1.0)||(dp2>1.0)) NoFails = false;
			  }
				
			if(NoFails) 
			  {
			    std::cout<<"Propose throwing out seed " << isd<<std::endl;
			    RedundantSeeds[isd]=true;
			  }
		      }
		  }
	      }
	  }
	
	
      }
    
    // Now we need to throw out all the seeds we marked as redundant
    std::map<int, std::map<int, double> > FilteredMap;
   
    std::vector<int> OldNew;
    for(size_t i=0; i!=Seeds.size(); ++i)
      {
	OldNew.push_back(i);
      }

    for(int i=Seeds.size()-1; i!=-1; --i)
      {
	if(RedundantSeeds[i])
	  {
	    Seeds.erase(Seeds.begin()+i);
	    OldNew.erase(OldNew.begin()+i);
	  }
      }
    
    for(size_t i=0; i!=OldNew.size(); ++i)
      {
	for(size_t j=0; j!=OldNew.size(); ++j)
	  {
	    FilteredMap[i][j] = ConnectionMap[OldNew[i]][OldNew[j]];
	  }
      }
    
    ConnectionMap = FilteredMap;
    

    // Deal with loops by throwing out the longest paths.
    for(size_t i=0; i!=Seeds.size(); ++i)
      {
	for(size_t j=0; j!=i; ++j)
	  {
	    if(ConnectionMap[i][j]>0)
	      {
		for(size_t k=0; k!=Seeds.size(); ++k)
		  {
		    if((ConnectionMap[i][k]>0)&&(ConnectionMap[k][j]>0))
		      {
			if( ( ConnectionMap[i][k] > ConnectionMap[i][j]) && ( ConnectionMap[i][k] > ConnectionMap[j][k]))
			  ConnectionMap[i][j] = ConnectionMap[j][i] = ConnectionMap[k][j] = ConnectionMap[j][k] = 0; 
			else if( ( ConnectionMap[i][j] > ConnectionMap[i][k]) && ( ConnectionMap[i][j] > ConnectionMap[j][k]))
			  ConnectionMap[j][k] = ConnectionMap[k][j] = ConnectionMap[i][k] = ConnectionMap[k][i] = 0;
			else
			  ConnectionMap[i][j] = ConnectionMap[j][i] =  ConnectionMap[i][k] = ConnectionMap[k][i] = 0;
			
		      }
		  }
	      }
	  }
      }

    for(size_t i=0; i!=Seeds.size(); ++i)
      {
	int Count=0;
	for(size_t j=0; j!=Seeds.size(); ++j)
	  {
	    if(ConnectionMap[i][j]>0)
	      Count++;
	  }
	std::cout<<" After delooping, seed " << i << " is connected " << Count << " times"<<std::endl;

      }
    return FilteredMap;
    
  }

  //----------------------------------------------------------------------
  
  std::vector<trkf::BezierTrack> FeatureTracker::GenerateBezierTracks(std::map<int,std::map<int,double> > ConnMap, std::vector<recob::Seed> Seeds)
  {
    std::vector<trkf::BezierTrack> ReturnVec;
    
    std::vector<std::vector<int> > CollectedSeeds;
    std::map<int, bool>  AlreadyCounted;
    
    
  
    bool StillFinding=true;

    
    for(size_t baseseed=0; baseseed!=Seeds.size(); ++baseseed)
      {
	if(!AlreadyCounted[baseseed]) 
	  {
	    std::vector<int>     SeedsThisTrack;
	    SeedsThisTrack.clear();
	 
	    SeedsThisTrack.push_back(baseseed);	    
	    while(StillFinding)
	      {
		StillFinding=false;
		for(size_t i=0; i!=SeedsThisTrack.size(); ++i)
		  {
		    for(size_t j=0; j!=Seeds.size(); ++j)
		      {
			if((!AlreadyCounted[j])&&(ConnMap[SeedsThisTrack[i]][j]>0))
			  {
			    SeedsThisTrack.push_back(j);
			    AlreadyCounted[j]=true;
			    StillFinding=true;	    
			  }
		      }
		  }
	      }
	    CollectedSeeds.push_back(SeedsThisTrack);
	  }
      }
    std::cout<<"Found " << CollectedSeeds.size()<< " sensible collections.  Sizes:"<<std::endl;
    for(size_t i=0; i!=CollectedSeeds.size();  ++i)
      std::cout<<"  " << CollectedSeeds.at(i).size()<<std::endl;
    
    return std::vector<trkf::BezierTrack>();
  }



  //----------------------------------------------------------------------
  void FeatureTracker::endJob()
  {
    
  }
}

