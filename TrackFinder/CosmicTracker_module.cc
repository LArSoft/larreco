////////////////////////////////////////////////////////////////////////
//
//  CosmicTracker
//
//  Tracker to reconstruct cosmic ray muons and neutrino interactions
// 
//  tjyang@fnal.gov
// 
//  This algorithm is based on orignal idea in Track3Dreco
//
//  ** Modified by Muhammad Elnimr to check multiple TPCs for the 35t prototype.
//   
//  mmelnimr@as.ua.edu
//  April 2014
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/ClusterMatchTQ.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1D.h"
#include "TVirtualFitter.h"

//2-D weighted fit of time vs wire
std::vector<double> vwire;
std::vector<double> vtime;
std::vector<double> vph;

void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  //minimisation function computing the sum of squares of residuals
  f = 0;
  for (size_t i = 0; i<vwire.size(); ++i){
    double x = vwire[i];
    double y = par[0]+par[1]*x+par[2]*x*x;
    //using ph^2 to suppress low ph hits
    f += vph[i]*vph[i]*(vtime[i]-y)*(vtime[i]-y);
  }
}

bool SortByWire (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) { 
  return h1->WireID().Wire < h2->WireID().Wire;
}

bool sp_sort_x0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] < xyz2[0];
}

bool sp_sort_x1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] > xyz2[0];
}

bool sp_sort_y0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[1] < xyz2[1];
}

bool sp_sort_y1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[1] > xyz2[1];
}

bool sp_sort_z0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] < xyz2[2];
}

bool sp_sort_z1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] > xyz2[2];
}

namespace trkf {

  class CosmicTracker : public art::EDProducer {
    
  public:
    
    explicit CosmicTracker(fhicl::ParameterSet const& pset);
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:

    cluster::ClusterMatchTQ  fClusterMatch;

    //double          fKScut;              ///< tolerance for cluster matching based on KS test.

    double          ftmatch;             ///< tolerance for time matching (in time samples) 
    
    double          fsmatch;             ///< tolerance for distance matching (in cm)

    double          ftoler1;             ///< tolerance for using hits in fit

    double          ftoler2;             ///< tolerance for using hits in fit

    std::string     fClusterModuleLabel; ///< label for input cluster collection
    
    bool             fdebug;

    int             fisohitcut;
    
    std::string     fSortDir;            ///< sort space points 

    bool            fCleanUpHits;        ///< flag to remove outlier hits

    bool            fDirSPS;             ///< calculate direction cosine for each space point
  
  }; // class CosmicTracker

}

namespace trkf {

  //-------------------------------------------------
  CosmicTracker::CosmicTracker(fhicl::ParameterSet const& pset) :
    fClusterMatch(pset.get< fhicl::ParameterSet >("ClusterMatch"))
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Track>                        >();
    produces< std::vector<recob::SpacePoint>                   >();
    produces< art::Assns<recob::Track,      recob::Cluster>    >();
    produces< art::Assns<recob::Track,      recob::SpacePoint> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit>        >();
    produces< art::Assns<recob::Track,      recob::Hit>        >();

  }

  //-------------------------------------------------
  void CosmicTracker::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
    ftmatch                 = pset.get< double >("TMatch");
    fsmatch                 = pset.get< double >("SMatch");
    ftoler1                 = pset.get< double >("Toler1");
    ftoler2                 = pset.get< double >("Toler2");
    fdebug                  = pset.get< bool   >("Debug");
    fisohitcut              = pset.get< int    >("IsoHitCut");
    fSortDir                = pset.get< std::string >("SortDirection","+z");
    fCleanUpHits            = pset.get< bool   >("CleanUpHits");
    fDirSPS                 = pset.get< bool   >("DirSPS");
  }

  //-------------------------------------------------
  void CosmicTracker::beginJob()
  {
  }

  //-------------------------------------------------
  void CosmicTracker::endJob()
  {
  }

  //------------------------------------------------------------------------------------//
  void CosmicTracker::produce(art::Event& evt){
  
    // get services
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    std::unique_ptr<std::vector<recob::Track>      >              tcol (new std::vector<recob::Track>);           
    std::unique_ptr<std::vector<recob::SpacePoint> >                 spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::Cluster> >    tcassn(new art::Assns<recob::Track, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit> >        thassn(new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> >   shassn(new art::Assns<recob::SpacePoint, recob::Hit>);

    double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
    //double presamplings = detprop->TriggerOffset(); // presamplings in ticks  
    //double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
    double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
    double Efield_drift = larprop->Efield(0);  // Electric Field in the drift region in kV/cm
    double Temperature = larprop->Temperature();  // LAr Temperature in K

    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
    double timepitch = driftvelocity*timetick;                         //time sample (cm) 


    // get input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    std::vector<art::Ptr<recob::Cluster> > clusterlist;
    if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);

    art::FindManyP<recob::Hit> fm(clusterListHandle, evt, fClusterModuleLabel);

    fClusterMatch.ClusterMatch(clusterlist,fm);
    std::vector<std::vector<unsigned int> > &matchedclusters = fClusterMatch.matchedclusters;

    /////////////////////////////////////////////////////
    /////// 2D Track Matching and 3D Track Reconstruction
    /////////////////////////////////////////////////////
    
    ///Prepare fitter
    double arglist[10];
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    if (fCleanUpHits){
      fitter->SetFCN(myfcn);
      arglist[0] = -1;
      fitter->ExecuteCommand("SET PRIN",arglist,1);
    }
    //fit each cluster in 2D using pol2, iterate once to remove outliers
    for (size_t itrk = 0; itrk<matchedclusters.size(); ++itrk){//loop over tracks

      //all the clusters associated with the current track
      art::PtrVector<recob::Cluster> clustersPerTrack;
      for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){
        art::Ptr <recob::Cluster> cluster(clusterListHandle,matchedclusters[itrk][iclu]);
        clustersPerTrack.push_back(cluster);
      }

      //save time/hit information along track trajectory
      std::vector<std::map<int,double> > vtimemap;
      std::vector<std::map<int,art::Ptr<recob::Hit> > > vhitmap;

      for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){//loop over clusters

        vwire.clear();
        vtime.clear();
        vph.clear();
        //fit hits time vs wire with pol2
        std::vector< art::Ptr<recob::Hit> > hits = fm.at(matchedclusters[itrk][iclu]);
        std::sort(hits.begin(), hits.end(), SortByWire);
        if (fCleanUpHits){
          double dtdw = 0;
          if (clusterlist[matchedclusters[itrk][iclu]]->StartTick()-
              clusterlist[matchedclusters[itrk][iclu]]->EndTick()){
            dtdw = (clusterlist[matchedclusters[itrk][iclu]]->EndWire()-
                    clusterlist[matchedclusters[itrk][iclu]]->StartWire())/
              (clusterlist[matchedclusters[itrk][iclu]]->EndTick()-
               clusterlist[matchedclusters[itrk][iclu]]->StartTick());
          }
          fitter->SetParameter(0,"p0",clusterlist[matchedclusters[itrk][iclu]]->StartTick()-dtdw-detprop->GetXTicksOffset(clusterlist[matchedclusters[itrk][iclu]]->Plane().Plane,clusterlist[matchedclusters[itrk][iclu]]->Plane().TPC,clusterlist[matchedclusters[itrk][iclu]]->Plane().Cryostat),0.1,0,0);
          fitter->SetParameter(1,"p1",std::tan(clusterlist[matchedclusters[itrk][iclu]]->StartAngle()),0.1,0,0);
          fitter->SetParameter(2,"p2",0,0.1,0,0);
        }
        for (size_t ihit = 0; ihit<hits.size(); ++ihit){//loop over hits
          geo::WireID hitWireID = hits[ihit]->WireID();
          unsigned int w = hitWireID.Wire;
          vwire.push_back(w);
          double time = hits[ihit]->PeakTime();
          time -= detprop->GetXTicksOffset(hits[ihit]->WireID().Plane,
                                           hits[ihit]->WireID().TPC,
                                           hits[ihit]->WireID().Cryostat);
          vtime.push_back(time);
          vph.push_back(hits[ihit]->Integral());
        }

        if (fCleanUpHits){
          arglist[0] = 0;
          if (vwire.size()>2) fitter->ExecuteCommand("MIGRAD", arglist, 0);
          else{
            fitter->SetParameter(0,"p0",vtime[0],0.1,0,0);
            fitter->SetParameter(1,"p1",0,0.1,0,0);
            fitter->SetParameter(2,"p2",0,0.1,0,0);
          }
          //remove outliers
          for (auto iw = vwire.begin(), it = vtime.begin(), iph = vph.begin(); iw!=vwire.end(); ){
            double y = fitter->GetParameter(0)+
              fitter->GetParameter(1)*(*iw)+
              fitter->GetParameter(2)*(*iw)*(*iw);
            if (std::abs(*it-y)>ftoler1){
              iw = vwire.erase(iw);
              it = vtime.erase(it);
              iph = vph.erase(iph);
            }
            else{
              ++iw;
              ++it;
              ++iph;
            }
          }
          
          //refit
          if (vwire.size()>2) fitter->ExecuteCommand("MIGRAD", arglist, 0);
        }
        
        std::map<int,double> timemap;
        std::map<int,double> phmap;
        std::map<int,art::Ptr<recob::Hit> > hitmap;

        //find hit on each wire along the fitted line
        for (size_t ihit = 0; ihit<hits.size(); ++ihit){//loop over hits
          geo::WireID hitWireID = hits[ihit]->WireID();
          unsigned int w = hitWireID.Wire;
          vwire.push_back(w);
          double time = hits[ihit]->PeakTime();
          time -= detprop->GetXTicksOffset(hits[ihit]->WireID().Plane,
                                           hits[ihit]->WireID().TPC,
                                           hits[ihit]->WireID().Cryostat);
          double ph = hits[ihit]->Integral();
          if (fCleanUpHits){
            if (ph>(phmap[w])){
              double y = fitter->GetParameter(0)+
                fitter->GetParameter(1)*w+
                fitter->GetParameter(2)*w*w;
              if (std::abs(time-y)<ftoler2){
                phmap[w] = ph;
                timemap[w] = time;
                hitmap[w] = hits[ihit];
              }
            }
          }
          else{
            phmap[w] = ph;
            timemap[w] = time;
            hitmap[w] = hits[ihit];
            //std::cout<<w<<" "<<time<<" "<<ph<<" "<<hits[ihit]->WireID().Plane<<std::endl;
          }          
        }//ihit
        vtimemap.push_back(timemap);
        vhitmap.push_back(hitmap);
      }//iclu
      
      
      // Remove isolated hits
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
        auto ihit = vhitmap[iclu].begin();
        for (auto itime = vtimemap[iclu].begin(); itime!=vtimemap[iclu].end();){
          int diffw0 = 0;
          int diffw1 = 0;
          if (itime!=vtimemap[iclu].begin()){
            auto itime0 = std::prev(itime,1);
            diffw0 = std::abs(itime->first - itime0->first);
          }
          else diffw0 = 99999;
          auto itime1 = std::next(itime,1);
          if (itime1!=vtimemap[iclu].end()){
            diffw1 = abs(itime->first - itime1->first);
          }
          else diffw1 = 99999;
          if (diffw0>fisohitcut&&diffw1>fisohitcut){
            vtimemap[iclu].erase(itime++);
            vhitmap[iclu].erase(ihit++);
          }
          else{
            ++itime;
            ++ihit;
          }
        }
      }

      // Find two clusters with the most numbers of hits, and time ranges
      int iclu1 = -1;
      int iclu2 = -1;
      int iclu3 = -1;
      unsigned maxnumhits0 = 0;
      unsigned maxnumhits1 = 0;
    
      std::vector<double> tmin(vtimemap.size());
      std::vector<double> tmax(vtimemap.size());
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
        tmin[iclu] = 1e9;
        tmax[iclu] = -1e9;
      }
    
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
        for (auto itime = vtimemap[iclu].begin(); itime!=vtimemap[iclu].end(); ++itime){
          if (itime->second>tmax[iclu]){
            tmax[iclu] = itime->second;
          }
          if (itime->second<tmin[iclu]){
            tmin[iclu] = itime->second;
          }
        }
        if (vtimemap[iclu].size()>maxnumhits0){
          if (iclu1!=-1){
            iclu2 = iclu1;
            maxnumhits1 = maxnumhits0;
          }
          iclu1 = iclu;
          maxnumhits0 = vtimemap[iclu].size();
        }
        else if (vtimemap[iclu].size()>maxnumhits1){
          iclu2 = iclu;
          maxnumhits1 = vtimemap[iclu].size();
        }
      }
    
      std::swap(iclu1,iclu2); //now iclu1 has fewer hits than iclu2

      for (int iclu = 0; iclu<(int)vtimemap.size(); ++iclu){
        if (iclu!=iclu1&&iclu!=iclu2) iclu3 = iclu;
      }
    
      if (iclu1!=-1&&iclu2!=-1){
        //select hits in a common time range
        auto ihit = vhitmap[iclu1].begin();
        auto itime = vtimemap[iclu1].begin();
        while (itime!=vtimemap[iclu1].end()){
          if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
              itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
            vtimemap[iclu1].erase(itime++);
            vhitmap[iclu1].erase(ihit++);
          }
          else{
            ++itime;
            ++ihit;
          }
        }

        ihit = vhitmap[iclu2].begin();
        itime = vtimemap[iclu2].begin();
        while (itime!=vtimemap[iclu2].end()){
          if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
              itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
            vtimemap[iclu2].erase(itime++);
            vhitmap[iclu2].erase(ihit++);
          }
          else{
            ++itime;
            ++ihit;
          }
        }
      
        //if one cluster is empty, replace it with iclu3
        if (!vtimemap[iclu1].size()){
          if (iclu3!=-1){
            std::swap(iclu3,iclu1);
          }
        }
        if (!vtimemap[iclu2].size()){
          if (iclu3!=-1){
            std::swap(iclu3,iclu2);
            std::swap(iclu1,iclu2);
          }
        }
        if ((!vtimemap[iclu1].size())||(!vtimemap[iclu2].size())) continue;

        size_t spStart = spcol->size();
        std::vector<recob::SpacePoint> spacepoints;
        //TVector3 startpointVec,endpointVec, DirCos;

        bool rev = false;
        auto times1 = vtimemap[iclu1].begin();
        auto timee1 = vtimemap[iclu1].end();
        --timee1;
        auto times2 = vtimemap[iclu2].begin();
        auto timee2 = vtimemap[iclu2].end();
        --timee2;

        double ts1 = times1->second;
        double te1 = timee1->second;
        double ts2 = times2->second;
        double te2 = timee2->second;
            
        //find out if we need to flip ends
        if (std::abs(ts1-ts2)+std::abs(te1-te2)>std::abs(ts1-te2)+std::abs(te1-ts2)){
          rev = true;
        }
//        std::cout<<times1->first<<" "<<times1->second<<std::endl;
//        std::cout<<timee1->first<<" "<<timee1->second<<std::endl;
//        std::cout<<times2->first<<" "<<times2->second<<std::endl;
//        std::cout<<timee2->first<<" "<<timee2->second<<std::endl;
        std::vector<double> vtracklength;
      
        for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
        
          double tracklength = 0;
          if (vtimemap[iclu].size()==1){
            tracklength = wire_pitch;
          }
          else{
            double t0 = 0., w0 = 0.;
            for (auto iw = vtimemap[iclu].begin(); iw!=vtimemap[iclu].end(); ++iw){
              if (iw==vtimemap[iclu].begin()){
                w0 = iw->first;
                t0 = iw->second;
              }
              else{
                tracklength += std::sqrt(std::pow((iw->first-w0)*wire_pitch,2)+std::pow((iw->second-t0)*timepitch,2));
                w0 = iw->first;
                t0 = iw->second;             
              }
            }
          }
          vtracklength.push_back(tracklength);
        }
      
        std::map<int,int> maxhitsMatch;

	auto besthit2 = vhitmap[iclu2].end();

        auto ihit1 = vhitmap[iclu1].begin();
        for (auto itime1 = vtimemap[iclu1].begin(); 
             itime1 != vtimemap[iclu1].end(); 
             ++itime1, ++ihit1){//loop over min-hits
          art::PtrVector<recob::Hit> sp_hits;
          sp_hits.push_back(ihit1->second);
          double hitcoord[3];
          double length1 = 0;
          if (vtimemap[iclu1].size()==1){
            length1 = wire_pitch;
          }
          else{
            for (auto iw1 = vtimemap[iclu1].begin(); iw1!=itime1; ++iw1){
              auto iw2 = iw1;
              ++iw2;
              length1 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)
                                   +std::pow((iw1->second-iw2->second)*timepitch,2));
            }
          }
          double difference = 1e10; //distance between two matched hits
          auto matchedtime = vtimemap[iclu2].end();
          auto matchedhit  = vhitmap[iclu2].end();
        
          auto ihit2 = vhitmap[iclu2].begin();
          for (auto itime2 = vtimemap[iclu2].begin(); 
               itime2!=vtimemap[iclu2].end(); 
               ++itime2, ++ihit2){//loop over max-hits
            if (maxhitsMatch[itime2->first]) continue;
	    if (besthit2!=vhitmap[iclu2].end()){
	      if (std::abs(besthit2->first-ihit2->first)>100) continue;
	    }
            bool timematch = std::abs(itime1->second-itime2->second)<ftmatch;
            if (timematch){
	      double length2 = 0;
	      if (vtimemap[iclu2].size()==1){
		length2 = wire_pitch;
	      }
	      else{
		for (auto iw1 = vtimemap[iclu2].begin(); iw1!=itime2; ++iw1){
		  auto iw2 = iw1;
		  ++iw2;
		  length2 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)+std::pow((iw1->second-iw2->second)*timepitch,2));
		}
	      }
	      if (rev) length2 = vtracklength[iclu2] - length2;
	      length2 = vtracklength[iclu1]/vtracklength[iclu2]*length2;

	      if(std::abs(length2-length1)<difference){
		difference = std::abs(length2-length1);
		matchedtime = itime2;
		matchedhit = ihit2;
	      }
            }
          }//loop over hits2
          if (difference<fsmatch){
	    besthit2 = matchedhit;
	    //add back XTicksOffset which will be removed in ConvertTicksToX
            hitcoord[0] = detprop->ConvertTicksToX(matchedtime->second+
						   detprop->GetXTicksOffset(clusterlist[matchedclusters[itrk][iclu1]]->Plane().Plane,
									    clusterlist[matchedclusters[itrk][iclu1]]->Plane().TPC,
									    clusterlist[matchedclusters[itrk][iclu1]]->Plane().Cryostat),
						   clusterlist[matchedclusters[itrk][iclu1]]->Plane().Plane,
						   clusterlist[matchedclusters[itrk][iclu1]]->Plane().TPC,
						   clusterlist[matchedclusters[itrk][iclu1]]->Plane().Cryostat);
            hitcoord[1] = -1e10;
            hitcoord[2] = -1e10;

            //WireID is the exact segment of the wire where the hit is on (1 out of 3 for the 35t)
            geo::WireID c1=(ihit1->second)->WireID();
            geo::WireID c2=(matchedhit->second)->WireID();
            geo::WireIDIntersection tmpWIDI;            
	    bool sameTpcOrNot=geom->WireIDsIntersect(c1,c2, tmpWIDI);
                    
	    if(sameTpcOrNot){
	      hitcoord[1]=tmpWIDI.y;
	      hitcoord[2]=tmpWIDI.z;
	    }

            if (hitcoord[1]>-1e9&&hitcoord[2]>-1e9){
              maxhitsMatch[matchedtime->first] = 1;
              sp_hits.push_back(matchedhit->second);
            }
          }
	  else{
	    besthit2 = vhitmap[iclu2].end();
	  }
          if (sp_hits.size()>1){
            double err[6] = {util::kBogusD};
            recob::SpacePoint mysp(hitcoord, 
                                   err, 
                                   util::kBogusD, 
                                   spStart + spacepoints.size());//3d point at end of track
            spacepoints.push_back(mysp);
            spcol->push_back(mysp);        
            util::CreateAssn(*this, evt, *spcol, sp_hits, *shassn);
          }
        }//loop over hits1
      
        size_t spEnd = spcol->size();

        if (fSortDir == "+x"){
          std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_x0);
          std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_x0);
        }
        if (fSortDir == "-x"){
          std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_x1);
          std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_x1);
        }
        if (fSortDir == "+y"){
          std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_y0);
          std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_y0);
        }
        if (fSortDir == "-y"){
          std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_y1);
          std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_y1);
        }
        if (fSortDir == "+z"){
          std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_z0);
          std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_z0);
        }
        if (fSortDir == "-z"){
          std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_z1);
          std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_z1);
        }

        // Add the 3D track to the vector of the reconstructed tracks
        if(spacepoints.size()>0){

          // make a vector of the trajectory points along the track
          std::vector<TVector3> xyz(spacepoints.size());
          for(size_t s = 0; s < spacepoints.size(); ++s){
            xyz[s] = TVector3(spacepoints[s].XYZ());
          }        
          //Calculate track direction cosines 
          TVector3 startpointVec,endpointVec, DirCos;
          startpointVec = xyz[0];
          endpointVec = xyz.back();
          DirCos = endpointVec - startpointVec;
          //SetMag casues a crash if the magnitude of the vector is zero
          try
            {
              DirCos.SetMag(1.0);//normalize vector
            }
          catch(...){std::cout<<"The Spacepoint is infinitely small"<<std::endl;
            continue;
          }
          //std::cout<<DirCos.x()<<" "<<DirCos.y()<<" "<<DirCos.z()<<std::endl;
          std::vector<TVector3> dircos(spacepoints.size(), DirCos);          
          if (fDirSPS){//calculat direction for each spacepoint
            for (int s = 0; s < int(xyz.size()); ++s){
              int np = 0;
              std::vector<double> vx;
              std::vector<double> vy;
              std::vector<double> vz;
              std::vector<double> vs;
              vx.push_back(xyz[s].x());
              vy.push_back(xyz[s].y());
              vz.push_back(xyz[s].z());
              vs.push_back(0);
              ++np;
              for (int ip = 1; ip<int(xyz.size()); ++ip){
                if (s-ip>=0){
                  vx.push_back(xyz[s-ip].x());
                  vy.push_back(xyz[s-ip].y());
                  vz.push_back(xyz[s-ip].z());
                  double dis = 0;
                  for (int j = s-ip; j<s; ++j){
                    dis += -sqrt(pow(xyz[j].x()-xyz[j+1].x(),2)+
                                 pow(xyz[j].y()-xyz[j+1].y(),2)+
                                 pow(xyz[j].z()-xyz[j+1].z(),2));
                  }
                  vs.push_back(dis);
                  ++np;
                  if (np==5) break;
                }
                if (s+ip<int(xyz.size())){
                  vx.push_back(xyz[s+ip].x());
                  vy.push_back(xyz[s+ip].y());
                  vz.push_back(xyz[s+ip].z());
                  double dis = 0;
                  for (int j = s; j<s+ip; ++j){
                    dis += sqrt(pow(xyz[j].x()-xyz[j+1].x(),2)+
                                pow(xyz[j].y()-xyz[j+1].y(),2)+
                                pow(xyz[j].z()-xyz[j+1].z(),2));
                  }
                  vs.push_back(dis);
                  ++np;
                  if (np==5) break;
                }
              }
              double kx = 0, ky = 0, kz = 0;
              if (np>=2){//at least two points
                TGraph *xs = new TGraph(np,&vs[0],&vx[0]);
                //for (int i = 0; i<np; i++) std::cout<<i<<" "<<vs[i]<<" "<<vx[i]<<" "<<vy[i]<<" "<<vz[i]<<std::endl;
                try{
                  if (np>2){
                    xs->Fit("pol2","Q");
                  }
                  else{
                    xs->Fit("pol1","Q");
                  }
                  TF1 *pol = 0;
                  if (np>2) pol = (TF1*) xs->GetFunction("pol2");
                  else pol = (TF1*) xs->GetFunction("pol1");
                  kx = pol->GetParameter(1);
                  //std::cout<<xyz3d[0]<<" "<<kx<<std::endl;
                }
                catch(...){
                  mf::LogWarning("CosmicTracker") <<"Fitter failed";
                }
                delete xs;
                TGraph *ys = new TGraph(np,&vs[0],&vy[0]);
                try{
                  if (np>2){
                    ys->Fit("pol2","Q");
                  }
                  else{
                    ys->Fit("pol1","Q");
                  }
                  TF1 *pol = 0;
                  if (np>2) pol = (TF1*) ys->GetFunction("pol2");
                  else pol = (TF1*) ys->GetFunction("pol1");
                  ky = pol->GetParameter(1);
                  //std::cout<<xyz3d[1]<<" "<<ky<<std::endl;
                }
                catch(...){
                  mf::LogWarning("CosmicTracker") <<"Fitter failed";
                }
                delete ys;
                TGraph *zs = new TGraph(np,&vs[0],&vz[0]);
                try{
                  if (np>2){
                    zs->Fit("pol2","Q");
                  }
                  else{
                    zs->Fit("pol1","Q");
                  }
                  TF1 *pol = 0;
                  if (np>2) pol = (TF1*) zs->GetFunction("pol2");
                  else pol = (TF1*) zs->GetFunction("pol1");
                  kz = pol->GetParameter(1);
                  //std::cout<<xyz3d[2]<<" "<<kz<<std::endl;
                }
                catch(...){
                  mf::LogWarning("CosmicTracker") <<"Fitter failed";
                }
                delete zs;
                if (kx||ky||kz){
                  double tot = sqrt(kx*kx+ky*ky+kz*kz);
                  kx /= tot;
                  ky /= tot;
                  kz /= tot;
                  dircos[s].SetXYZ(kx,ky,kz);
                  //std::cout<<s<<" "<<kx<<" "<<ky<<" "<<kz<<std::endl;
                }
              }//np>=2
            }//loop over space points
          }
          std::vector< std::vector<double> > dQdx;
          std::vector<double> mom(2, util::kBogusD);
          tcol->push_back(recob::Track(xyz, dircos, dQdx, mom, tcol->size()));
        
          // make associations between the track and space points
          util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);
        
          // now the track and clusters
          util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *tcassn);
        
          // and the hits and track
          if (!fdebug){
            const std::vector<unsigned int>& TrackClusters = matchedclusters[itrk];
            for (size_t cpt = 0; cpt < TrackClusters.size(); ++cpt) {
              util::CreateAssn
                (*this, evt, *tcol, fm.at(TrackClusters[cpt]), *thassn);
            } // for clusters of the track
          }
          else{
            std::vector<art::Ptr<recob::Hit> > trkhits;
            for(auto ihit = vhitmap[iclu1].begin(); ihit!=vhitmap[iclu1].end();++ihit){
              trkhits.push_back(ihit->second);
            }
            for(auto ihit = vhitmap[iclu2].begin(); ihit!=vhitmap[iclu2].end();++ihit){
              trkhits.push_back(ihit->second);
            }
            util::CreateAssn(*this, evt, *tcol, trkhits, *thassn);
          }
        }
      }//if iclu1&&iclu2
    }//itrk                 

    mf::LogVerbatim("Summary") << std::setfill('-') 
                               << std::setw(175) 
                               << "-" 
                               << std::setfill(' ');
    mf::LogVerbatim("Summary") << "CosmicTracker Summary:";
    for(unsigned int i = 0; i<tcol->size(); ++i) mf::LogVerbatim("Summary") << tcol->at(i) ;
    mf::LogVerbatim("Summary") << "CosmicTracker Summary End:";
    
    evt.put(std::move(tcol));
    evt.put(std::move(spcol));
    evt.put(std::move(tspassn));
    evt.put(std::move(tcassn));
    evt.put(std::move(thassn));
    evt.put(std::move(shassn));
    
    return;
  }
  
  DEFINE_ART_MODULE(CosmicTracker)

} // namespace
