////////////////////////////////////////////////////////////////////////
// $Id: AggregateVertex_module.cc Exp $
//
// AggregateVertex module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef AGGREGATEVERTEX_H
#define AGGREGATEVERTEX_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "RecoBase/Cluster.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "Utilities/AssociationUtil.h"

#include <vector>
#include <string>

namespace vertex {


  class AggregateVertex : public art::EDProducer
  {

  public:

    explicit AggregateVertex(fhicl::ParameterSet const& pset);
    virtual ~AggregateVertex();

    void produce(art::Event& evt); 
    void beginJob(); 

    std::unique_ptr< std::vector<recob::Vertex> >  MatchV2T(art::Event& evt,
							  art::Assns<recob::Vertex, recob::Track>& vtassn,
							  art::Assns<recob::Vertex, recob::Shower>& vsassn,
							  art::Assns<recob::Vertex, recob::Hit>& vhassn);

  private:

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fEndPointModuleLabel;

    art::PtrVector<recob::EndPoint2D> feplist;
    art::PtrVector<recob::Track>      ftracklist;
    art::PtrVector<recob::EndPoint2D> feplistStrong;

  }; // class AggregateVertex

}  // Namespace vertex

namespace vertex {

  //-----------------------------------------------
  AggregateVertex::AggregateVertex(fhicl::ParameterSet const& pset) 
    : fDBScanModuleLabel  (pset.get< std::string >("DBScanModuleLabel"  ))
    , fHoughModuleLabel   (pset.get< std::string >("HoughModuleLabel"   ))
    , fTrack3DModuleLabel (pset.get< std::string >("Track3DModuleLabel" ))
    , fEndPointModuleLabel(pset.get< std::string >("EndPointModuleLabel"))
  {
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Vertex, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
  }
  
  //-----------------------------------------------
  AggregateVertex::~AggregateVertex()
  {
  }

  //-----------------------------------------------
  void AggregateVertex::beginJob()
  {
  }

  //-----------------------------------------------
  void AggregateVertex::produce(art::Event& evt) 
  {

    art::Handle< std::vector<recob::EndPoint2D> > epListHandle;
    evt.getByLabel(fEndPointModuleLabel,epListHandle);
    for(size_t ii = 0; ii < epListHandle->size(); ++ii){
      art::Ptr<recob::EndPoint2D> ep(epListHandle, ii);
      feplist.push_back(ep); // class member
    }

    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrack3DModuleLabel,trackListHandle);
    for(size_t ii = 0; ii < trackListHandle->size(); ++ii){
      art::Ptr<recob::Track> track(trackListHandle, ii);
      ftracklist.push_back(track); // class member
    }


    // Only match strong vertices to tracks.
    art::PtrVector<recob::EndPoint2D>::const_iterator epIter = feplist.begin();

    while (epIter != feplist.end()){            
      art::Ptr <recob::EndPoint2D> ep = (*epIter);  
      if (ep->ID() < 3) feplistStrong.push_back(ep); // -- cuz there's one less
      epIter++;
    }

    // We will suck up the hits out of each vertex and out of the tracks
    // and see if there's an overlap of a sufficient (>0) number of hits. If so,
    // call the track a match to the vtx. Then, .... stick the track pointer(s)
    // into the AggVertex object.  EC, 23-July-2010.
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track>  > vtassn(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> > vsassn(new art::Assns<recob::Vertex, recob::Shower>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit>    > vhassn(new art::Assns<recob::Vertex, recob::Hit>);

    std::unique_ptr< std::vector<recob::Vertex> > vcol (MatchV2T(evt, *vtassn, *vsassn, *vhassn));

    evt.put(std::move(vcol));
    evt.put(std::move(vtassn));
    evt.put(std::move(vsassn));
    evt.put(std::move(vhassn));
  }

  //-------------------------------------------------------------------------------
  std::unique_ptr< std::vector<recob::Vertex> >  AggregateVertex::MatchV2T(art::Event& evt,
									 art::Assns<recob::Vertex, recob::Track>& vtassn,
									 art::Assns<recob::Vertex, recob::Shower>& vsassn,
									 art::Assns<recob::Vertex, recob::Hit>& vhassn)
  {
    mf::LogInfo("AggregateVertex") << "AggregateEvent::MatchV2T(): (strong) vertexlistStrong"
				   << " and tracklist lengths are " 
				   << feplistStrong.size() << " and " <<  ftracklist.size() << ".";

    art::FindManyP<recob::Hit> fmhst(feplistStrong, evt, fEndPointModuleLabel);

    // Bail if there are no tracks or vertices.
    //  if (!((int)vertexlistStrong.size()) || !((int)tracklist.size())) return NULL;
    if (feplistStrong.isNull() || ftracklist.isNull()) {
      return std::unique_ptr< std::vector<recob::Vertex> > (new std::vector<recob::Vertex>);
    }

    // Loop on the vertices, and all the hits in each
    std::unique_ptr< std::vector<recob::Vertex> > verts(new std::vector<recob::Vertex>);

    art::FindManyP<recob::Hit> fmht(ftracklist, evt, fTrack3DModuleLabel);

    for(size_t epctr = 0; epctr < feplistStrong.size(); ++epctr){            
      art::PtrVector<recob::Track>  tlistAssoc; // Will fill with matching tracks.
      art::PtrVector<recob::Shower> slistAssoc; // Will fill with matching showers.
      
      std::vector< art::Ptr<recob::Hit> > hitvertexlistStrong = fmhst.at(epctr);
      
      // Should be just one hit per vtx, as per Josh, but we loop anyway.
      art::PtrVector<recob::Hit>::const_iterator hvIter = hitvertexlistStrong.begin();
      while (hvIter != hitvertexlistStrong.end())  {
	art::Ptr<recob::Hit> hitv = (*hvIter);

	// Now loop on the track hits
	for(size_t t = 0; t < ftracklist.size(); ++t){
	  
	  std::vector< art::Ptr<recob::Hit> > hittlist = fmht.at(t);

	  std::vector< art::Ptr<recob::Hit> >::const_iterator htIter = hittlist.begin();
	  while (htIter != hittlist.end())  {
	    
	    art::Ptr<recob::Hit> hitt = (*htIter);
	    
	    if (hitt==hitv){
	      //std::cout << "AggregateEvent::MatchV2T(): BooYiggity! Vtx/Trk hit match! hitt, hitv= " 
	      //<< hitt << " " <<hitv << std::endl;
	      //std::cout << "AggregateEvent::MatchV2T(): vtx, trk " << vtx<< " " <<trk << std::endl;
	      tlistAssoc.push_back(ftracklist[t]);
	      htIter = hittlist.end()-1; // jump to end of track hitlist, since we've satisfied match
	    }			
	    htIter++;
	  } // end hits associated to track
	} // end track loop
	hvIter++;
       }

      // Now if matching tracks were found for this vertex then create the recob::Vertex object.
      if (tlistAssoc.size()>0){
	/// \todo Really need to also determine the xyz position of the found vertex
	/// \todo Also should determine the ID for this vertex
	double xyz[3] = {-999., -999., -999.};
	verts->push_back(recob::Vertex(xyz, verts->size()));

	// associate the tracks to the vertex
	util::CreateAssn(*this, evt, *verts, tlistAssoc, vtassn);

	art::FindManyP<recob::Hit> fmhtl(tlistAssoc, evt, fTrack3DModuleLabel);

	//associate the track hits to the vertex
	for(size_t t = 0; t < tlistAssoc.size(); ++t){
	  std::vector< art::Ptr<recob::Hit> > hits = fmhtl.at(t);
	  util::CreateAssn(*this, evt, *verts, hits, vhassn);
	}
      }// end if there are tracks to be associated
    
      if (slistAssoc.size()>0){
	/// \todo Really need to also determine the xyz position of the found vertex
	/// \todo Also should determine the ID for this vertex
	double xyz[3] = {-999., -999., -999.};
	verts->push_back(recob::Vertex(xyz, verts->size()));

	// associate the showers to the vertex
	util::CreateAssn(*this, evt, *verts, slistAssoc, vsassn);

	art::FindManyP<recob::Hit> fmhsl(slistAssoc, evt, fTrack3DModuleLabel);

	//associate the shower hits to the vertex
	///\todo get a shower module label in line 196
	for(size_t t = 0; t < slistAssoc.size(); ++t){
	  std::vector< art::Ptr<recob::Hit> > hits = fmhsl.at(t);
	  util::CreateAssn(*this, evt, *verts, hits, vhassn);
	}
      }// end if there are showers to be associated

      //HnhinVtxes->Fill(hitvertexlistStrong.size(),1);
	      
    } // end loop over end points 
    
    return verts;
      
  } // end AggregateVertex::MatchV2T


}// end namespace

namespace vertex{

  DEFINE_ART_MODULE(AggregateVertex)

}
#endif // AGGREGATEVERTEX_H

