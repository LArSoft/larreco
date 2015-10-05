////////////////////////////////////////////////////////////////////////
// Class:       EMShower3D
// Module Type: producer
// File:        EMShower3D_module.cc
// Authors: dorota.stefan@cern.ch robert.sulej@cern.ch	
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"

#include "MCCheater/BackTracker.h"

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include <memory>

#include "DirOfGamma/DirOfGamma.h"

// ROOT includes
#include "TLorentzVector.h"
#include "TMathBase.h"

#include <iostream>
#include <fstream>

struct IniSeg
{
	size_t idcl1;
	size_t idcl2;
	size_t idcl3;
	size_t view1;
	size_t view2;
	size_t view3;
	pma::Track3D* track;
	std::vector< art::Ptr< recob::Hit > > hits1; 
	std::vector< art::Ptr< recob::Hit > > hits2; 
	std::vector< art::Ptr< recob::Hit > > hits3; 
};

namespace ems
{
	class EMShower3D;
}

class ems::EMShower3D : public art::EDProducer {
public:
  explicit EMShower3D(fhicl::ParameterSet const & p);
  
  EMShower3D(EMShower3D const &) = delete;
  EMShower3D(EMShower3D &&) = delete;
  EMShower3D & operator = (EMShower3D const &) = delete;
  EMShower3D & operator = (EMShower3D &&) = delete;

  void beginJob() override;

  void produce(art::Event & e) override;

  void reconfigure(fhicl::ParameterSet const& p);

	
private:
  	recob::Track ConvertFrom(pma::Track3D& src);
	recob::Track ConvertFrom2(pma::Track3D& src);
	recob::Cluster ConvertFrom(const std::vector< art::Ptr<recob::Hit> > & src);

	std::vector< ems::DirOfGamma* > CollectShower2D(art::Event const & e);
	void Link(art::Event const & e, std::vector< ems::DirOfGamma* > input);

	void Reoptimize();

	void Make3DSeg(art::Event const & e, std::vector< ems::DirOfGamma* > pair);

	// void Make3DPMA(art::Event const & e, unsigned int cryo, unsigned int tpc);

	// void SelectTrks(art::Event const & e);

	void MarkUsedCls(size_t idcl);

	void MarkUsedTrks(size_t idtrk);

	bool Validate(art::Event const & e, const pma::Track3D& src, size_t plane);
	bool Validate(std::vector< ems::DirOfGamma* > input, size_t id1, size_t id2, size_t c1, size_t c2, size_t plane3);

	void FilterOutSmallParts(
		double r2d, 
		const std::vector< art::Ptr<recob::Hit> >& hits_in,
					std::vector< art::Ptr<recob::Hit> >& hits_out);

	bool GetCloseHits(
		double r2d, 
		const std::vector< art::Ptr<recob::Hit> >& hits_in, 
		std::vector<size_t>& used,
		std::vector< art::Ptr<recob::Hit> >& hits_out);

	bool Has(const std::vector<size_t>& v, size_t idx);

	size_t LinkCandidates(art::Event const & e, std::vector< ems::DirOfGamma* > input, size_t id);
	
	std::vector< IniSeg > fInisegs;
	std::vector< IniSeg > fSeltracks;
	std::vector< IniSeg > fPMA3D;

	std::vector< std::vector< art::Ptr<recob::Hit> > > fClusters;

	std::vector< size_t > fClustersNotUsed;
	std::vector< size_t > fTracksNotUsed;

  	unsigned int fTrkIndex; unsigned int fClIndex;
	unsigned int fIniIndex;

	std::string fCluModuleLabel;
	std::string fTrk3DModuleLabel;

	pma::ProjectionMatchingAlg fProjectionMatchingAlg;

	art::Handle< std::vector< recob::Cluster > > fCluListHandle;

	// ofstream file;
};


ems::EMShower3D::EMShower3D(fhicl::ParameterSet const & p)
	: fProjectionMatchingAlg(p.get< fhicl::ParameterSet >("ProjectionMatchingAlg"))
{
	reconfigure(p);
  
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::Vertex> >();
	produces< std::vector<recob::Cluster> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< art::Assns<recob::Track, recob::Hit> >();
	produces< art::Assns<recob::Track, recob::Vertex> >();
	produces< art::Assns<recob::Cluster, recob::Hit> >();
	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
  	produces< art::Assns<recob::Track, recob::Cluster> >();

	// file.open("file.dat");
}

void ems::EMShower3D::beginJob()
{	
}

void ems::EMShower3D::reconfigure(fhicl::ParameterSet const & p)
{
	fCluModuleLabel = p.get< std::string >("ClustersModuleLabel");
  	fProjectionMatchingAlg.reconfigure(p.get< fhicl::ParameterSet >("ProjectionMatchingAlg"));
	fTrk3DModuleLabel = p.get< std::string >("Trk3DModuleLabel");
	
  return;
}

recob::Cluster ems::EMShower3D::ConvertFrom(const std::vector< art::Ptr<recob::Hit> > & src)
{
	
	return recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, src.size(), 0.0F, 0.0F, fClIndex,  src[0]->View(), src[0]->WireID().planeID());
}

recob::Track ems::EMShower3D::ConvertFrom(pma::Track3D& src)
{
	art::ServiceHandle<geo::Geometry> geom;
	size_t cryo = src.FrontCryo();
	const geo::CryostatGeo& cryostat = geom->Cryostat(cryo);

	std::vector< std::vector< double > > vdqdx;
	double dedxZ = fProjectionMatchingAlg.selectInitialHits(src, geo::kZ); //loop over planes
	double dedxV = fProjectionMatchingAlg.selectInitialHits(src, geo::kV);
	double dedxU = fProjectionMatchingAlg.selectInitialHits(src, geo::kU);

	std::vector< double > dqdx;
	dqdx.push_back(dedxZ); vdqdx.push_back(dqdx);
	dqdx.clear(); dqdx.push_back(dedxV); vdqdx.push_back(dqdx);
	dqdx.clear(); dqdx.push_back(dedxU); vdqdx.push_back(dqdx);

	std::vector< TVector3 > xyz, dircos; 

	for (size_t i = 0; i < src.size(); i++)
	{
		xyz.push_back(src[i]->Point3D());

		if (i < src.size() - 1)
		{
			TVector3 dc(src[i + 1]->Point3D());
			dc -= src[i]->Point3D();
			dc *= 1.0 / dc.Mag();
			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());
	}

	return recob::Track(xyz, dircos, vdqdx, std::vector< double >(2, util::kBogusD), fTrkIndex);
}

recob::Track ems::EMShower3D::ConvertFrom2(pma::Track3D& src)
{
	art::ServiceHandle<geo::Geometry> geom;
	size_t cryo = src.FrontCryo();
	const geo::CryostatGeo& cryostat = geom->Cryostat(cryo);

	std::vector< std::vector< double > > vdqdx;
	double dedxZ = fProjectionMatchingAlg.selectInitialHits(src, geo::kZ);
	double dedxV = fProjectionMatchingAlg.selectInitialHits(src, geo::kV);
	double dedxU = fProjectionMatchingAlg.selectInitialHits(src, geo::kU);

	std::vector< double > dqdx;
	dqdx.push_back(dedxZ); vdqdx.push_back(dqdx);
	dqdx.clear(); dqdx.push_back(dedxV); vdqdx.push_back(dqdx);
	dqdx.clear(); dqdx.push_back(dedxU); vdqdx.push_back(dqdx);

	std::vector< TVector3 > xyz, dircos;

	for (size_t i = 0; i < src.size(); i++)
	{
		xyz.push_back(src[i]->Point3D());

		if (i < src.size() - 1)
		{
			TVector3 dc(src[i + 1]->Point3D());
			dc -= src[i]->Point3D();
			dc *= 1.0 / dc.Mag();
			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());
	}

	return recob::Track(xyz, dircos, vdqdx, std::vector< double >(2, util::kBogusD), fIniIndex);
}

void ems::EMShower3D::produce(art::Event & e)
{
	art::ServiceHandle<geo::Geometry> geom;
	fSeltracks.clear();
	fInisegs.clear();
	fClusters.clear();
	fPMA3D.clear();
	fClustersNotUsed.clear();
	fTracksNotUsed.clear();

	std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
	std::unique_ptr< std::vector< recob::Vertex > > vertices(new std::vector< recob::Vertex >);
	std::unique_ptr< std::vector< recob::Cluster > > clusters(new std::vector< recob::Cluster >);
	std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);

	std::unique_ptr< art::Assns< recob::Track, recob::Hit > > trk2hit(new art::Assns< recob::Track, recob::Hit >);
	std::unique_ptr< art::Assns< recob::Track, recob::Vertex > > trk2vtx(new art::Assns< recob::Track, recob::Vertex >);
	std::unique_ptr< art::Assns< recob::Cluster, recob::Hit > > cl2hit(new art::Assns< recob::Cluster, recob::Hit >);
	std::unique_ptr< art::Assns< recob::Track, recob::Cluster > > trk2cl(new art::Assns< recob::Track, recob::Cluster >);
	std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
	std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);


  if (e.getByLabel(fCluModuleLabel, fCluListHandle))
	{
		fClustersNotUsed.clear();
		fTracksNotUsed.clear();
		fInisegs.clear();
		art::FindManyP< recob::Hit > fb(fCluListHandle, e, fCluModuleLabel);
	
		for (size_t id = 0; id < fCluListHandle->size(); id++)
		{
				std::vector< art::Ptr<recob::Hit> > hitlist;	
				hitlist = fb.at(id);

				if (hitlist.size() > 5) 
						fClustersNotUsed.push_back(id);
		}

		art::Handle< std::vector< recob::Track > > trkList;		
  		if ((e.getByLabel(fTrk3DModuleLabel, trkList)) && trkList->size())
				for (size_t id = 0; id < trkList->size(); id++)
					fTracksNotUsed.push_back(id);
	
		std::vector< ems::DirOfGamma* > showernviews = CollectShower2D(e); 	
		
		Link(e, showernviews);	

		while (fInisegs.size())
		{	
			fSeltracks.push_back(fInisegs[0]);
			fInisegs.erase(fInisegs.begin() + 0);
		}

		Reoptimize();

		/*SelectTrks(e);

		for (unsigned int c = 0; c < geom->Ncryostats(); c++)
		{
			const geo::CryostatGeo& cryo = geom->Cryostat(c);
			for (unsigned int t = 0; t < cryo.NTPC(); t++)
				Make3DPMA(e, c, t);
		}*/
	
		// conversion from pma track to recob::track

		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6], vtx_pos[3];
		for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

		fTrkIndex = 0;
		
		for (auto const trk : fSeltracks) 
		{
			tracks->push_back(ConvertFrom(*(trk.track)));

			vtx_pos[0] = trk.track->front()->Point3D().X();
			vtx_pos[1] = trk.track->front()->Point3D().Y();
			vtx_pos[2] = trk.track->front()->Point3D().Z();
			vertices->push_back(recob::Vertex(vtx_pos, fTrkIndex));

			fTrkIndex++;

			std::vector< art::Ptr< recob::Cluster > > cl2d; 
			cl2d.push_back( art::Ptr< recob::Cluster >(fCluListHandle, trk.idcl1) );
			cl2d.push_back( art::Ptr< recob::Cluster >(fCluListHandle, trk.idcl2) );

			std::vector< art::Ptr< recob::Hit > > hits2d; 
			art::PtrVector< recob::Hit > sp_hits;

			spStart = allsp->size();
			for (int h = trk.track->size() - 1; h >= 0; h--)
			{
				pma::Hit3D* h3d = (*trk.track)[h];
				hits2d.push_back(h3d->Hit2DPtr());

				if ((h == 0) ||
				      (sp_pos[0] != h3d->Point3D().X()) ||
				      (sp_pos[1] != h3d->Point3D().Y()) ||
				      (sp_pos[2] != h3d->Point3D().Z()))
				{
					if (sp_hits.size()) // hits assigned to the previous sp
					{
						util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
						sp_hits.clear();
					}
					sp_pos[0] = h3d->Point3D().X();
					sp_pos[1] = h3d->Point3D().Y();
					sp_pos[2] = h3d->Point3D().Z();
					allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
				}
				sp_hits.push_back(h3d->Hit2DPtr());
			}
			if (sp_hits.size()) // hits assigned to the last sp
			{
				util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
			}
			spEnd = allsp->size();

			if (vertices->size())
			{
				size_t vtx_idx = (size_t)(vertices->size() - 1);
				util::CreateAssn(*this, e, *tracks, *vertices, *trk2vtx, vtx_idx, vtx_idx + 1);
			}

			if (cl2d.size())
			{
				util::CreateAssn(*this, e, *tracks, cl2d, *trk2cl); 
			}

			if (hits2d.size())
			{
				util::CreateAssn(*this, e, *tracks, *allsp, *trk2sp, spStart, spEnd);
				util::CreateAssn(*this, e, *tracks, hits2d, *trk2hit); 
			}
		}

		fIniIndex = fTrkIndex + 1;
		for (auto const trk : fPMA3D) 
		{
			tracks->push_back(ConvertFrom2(*(trk.track)));

			fIniIndex++;

			std::vector< art::Ptr< recob::Cluster > > cl2d; 
			cl2d.push_back( art::Ptr< recob::Cluster >(fCluListHandle, trk.idcl1) );
			cl2d.push_back( art::Ptr< recob::Cluster >(fCluListHandle, trk.idcl2) );

			std::vector< art::Ptr< recob::Hit > > hits2d; 
			art::PtrVector< recob::Hit > sp_hits;

			spStart = allsp->size();
			for (int h = trk.track->size() - 1; h >= 0; h--)
			{
				pma::Hit3D* h3d = (*trk.track)[h];
				hits2d.push_back(h3d->Hit2DPtr());

				if ((h == 0) ||
				      (sp_pos[0] != h3d->Point3D().X()) ||
				      (sp_pos[1] != h3d->Point3D().Y()) ||
				      (sp_pos[2] != h3d->Point3D().Z()))
				{
					if (sp_hits.size()) // hits assigned to the previous sp
					{
						util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
						sp_hits.clear();
					}
					sp_pos[0] = h3d->Point3D().X();
					sp_pos[1] = h3d->Point3D().Y();
					sp_pos[2] = h3d->Point3D().Z();
					allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
				}
				sp_hits.push_back(h3d->Hit2DPtr());
			}
			if (sp_hits.size()) // hits assigned to the last sp
			{
				util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
			}
			spEnd = allsp->size();


			if (cl2d.size())
			{
				util::CreateAssn(*this, e, *tracks, cl2d, *trk2cl); 
			}

			if (hits2d.size())
			{
				util::CreateAssn(*this, e, *tracks, *allsp, *trk2sp, spStart, spEnd);
				util::CreateAssn(*this, e, *tracks, hits2d, *trk2hit); 
			}
		}

		// create cluster from hits, which were an input to find initial part of the cascade.
		fClIndex = 0;
		for (auto const& cl : fClusters) 
			if (cl.size())
			{
				clusters->push_back(ConvertFrom(cl)); 
				fClIndex++;

				util::CreateAssn(*this, e, *clusters, cl, *cl2hit); 
			}


		for (unsigned int i = 0; i < showernviews.size(); i++) 
				delete showernviews[i];
		
		for (unsigned int i = 0; i < fSeltracks.size(); i++)
				delete fSeltracks[i].track; 

		for (unsigned int i = 0; i < fInisegs.size(); i++)
				delete fInisegs[i].track; 

		for (unsigned int i = 0; i < fPMA3D.size(); i++)
				delete fPMA3D[i].track;  

}

		e.put(std::move(tracks));
		e.put(std::move(vertices));
		e.put(std::move(clusters));
		e.put(std::move(allsp));

		e.put(std::move(trk2hit));
		e.put(std::move(trk2vtx));
		e.put(std::move(cl2hit));
		e.put(std::move(trk2cl));
		e.put(std::move(trk2sp));
		e.put(std::move(sp2hit));
	
}

void ems::EMShower3D::Reoptimize()
{
	if (!fSeltracks.size()) return;
	const float min_dist = 3.0F;
	size_t ta = 0; 
	while (ta < (fSeltracks.size() - 1))
	{	
		size_t tb = ta + 1;
		bool found = false;
		while (tb < fSeltracks.size())
		{
			if (ta == tb) {tb++; continue;}
	
			TVector3 p1 = fSeltracks[ta].track->front()->Point3D();
			TVector3 p2 = fSeltracks[tb].track->front()->Point3D();
			float dist = std::sqrt(pma::Dist2(p1, p2));

			if (dist < min_dist)
				if ((fSeltracks[ta].idcl1 == fSeltracks[tb].idcl1) || (fSeltracks[ta].idcl1 == fSeltracks[tb].idcl2) ||
				    (fSeltracks[ta].idcl2 == fSeltracks[tb].idcl1) || (fSeltracks[ta].idcl2 == fSeltracks[tb].idcl2))
					{
						found = true;
						size_t view3 = fSeltracks[ta].view1; size_t idcl3 = fSeltracks[ta].idcl1;
						std::vector< art::Ptr< recob::Hit > > hits3 = fSeltracks[ta].hits1;
						std::vector< art::Ptr< recob::Hit > > hits = fSeltracks[ta].hits1;
						for (size_t h = 0; h < fSeltracks[ta].hits2.size(); ++h)
							hits.push_back(fSeltracks[ta].hits2[h]);

						if ((fSeltracks[tb].view1 != fSeltracks[ta].view1) && (fSeltracks[tb].view1 != fSeltracks[ta].view2))
						{
							view3 = fSeltracks[tb].view1;
							for (size_t h = 0; h < fSeltracks[tb].hits1.size(); ++h)	
								hits.push_back(fSeltracks[tb].hits1[h]); 
						}	
						if ((fSeltracks[tb].view2 != fSeltracks[ta].view1) && (fSeltracks[tb].view2 != fSeltracks[ta].view2))
						{
							view3 = fSeltracks[tb].view2;
							for (size_t h = 0; h < fSeltracks[tb].hits2.size(); ++h)
								hits.push_back(fSeltracks[tb].hits2[h]);
						}

						
						if ((view3 == fSeltracks[ta].view1) || (view3 == fSeltracks[ta].view2))
						{
							delete fSeltracks[ta].track;
							fSeltracks.erase(fSeltracks.begin() + ta);
						}
						else
						{
							pma::Track3D* track = fProjectionMatchingAlg.buildSegment(hits);

							if (pma::Dist2(track->back()->Point3D(), fSeltracks[ta].track->front()->Point3D()) < 
							    pma::Dist2(track->front()->Point3D(), fSeltracks[ta].track->front()->Point3D()))
									track->Flip();

							IniSeg initrack;
							initrack.idcl1 = fSeltracks[ta].idcl1; initrack.idcl3 = idcl3;
							initrack.view1 = fSeltracks[ta].view1; initrack.view3 = view3;
							initrack.hits1 = fSeltracks[ta].hits1; initrack.hits3 = hits3; 
							initrack.idcl2 = fSeltracks[ta].idcl2;
							initrack.view2 = fSeltracks[ta].view2;
							initrack.hits2 = fSeltracks[ta].hits2;
							initrack.track = track;	
		
							
							delete fSeltracks[tb].track; delete fSeltracks[ta].track;
							fSeltracks.erase(fSeltracks.begin() + tb); fSeltracks.erase(fSeltracks.begin() + ta);
							fSeltracks.push_back(initrack);
							
						}
					} 

			if (found) break;
			tb++;
		}

		if (!found) ta++;
	}
}

std::vector< ems::DirOfGamma* > ems::EMShower3D::CollectShower2D(art::Event const & e) 
{
	std::vector< ems::DirOfGamma* > input;

 	if (e.getByLabel(fCluModuleLabel, fCluListHandle))
	{
		art::FindManyP< recob::Hit > fb(fCluListHandle, e, fCluModuleLabel);
		for (unsigned int c = 0; c < fCluListHandle->size(); c++)
		{
			std::vector< art::Ptr<recob::Hit> > hitlist;	
			hitlist = fb.at(c);

			if (hitlist.size() > 5)
			{
				std::vector< art::Ptr<recob::Hit> > hits_out;
				FilterOutSmallParts(2.0, hitlist, hits_out);
				
				if (hits_out.size() > 5)
				{
					fClusters.push_back(hits_out);
					
					ems::DirOfGamma * sh = new ems::DirOfGamma(hits_out, 14, c); 
					
					if (sh->GetHits2D().size())
						input.push_back(sh);				
				}
			}
		}
	}

	/*if (input.size())
		for (size_t i = 0; i < input.size(); ++i)
		{
			for (size_t j = 0; j < input[i]->GetHits2D().size(); ++j)
				file << e.id().event() << " " << input[i]->GetHits2D().size() << " " << input[i]->GetFirstHit()->WireID().Plane << " " << input[i]->GetHits2D()[j]->GetPointCm().X() << " " << input[i]->GetHits2D()[j]->GetPointCm().Y() << " " << input[i]->GetFirstPoint().X() << " " << input[i]->GetFirstPoint().Y() << std::endl;

			for (size_t j = 0; j < input[i]->GetCandidates().size(); ++j)
				file << e.id().event() << " " << input[i]->GetHits2D().size() << " " << input[i]->GetFirstHit()->WireID().Plane << " " << input[i]->GetCandidates().size() << " " << input[i]->GetCandidates()[j].GetPosition().X() << " " << input[i]->GetCandidates()[j].GetPosition().Y() << std::endl;			
		}*/

	return input;
}

void ems::EMShower3D::Link(art::Event const & e, std::vector< ems::DirOfGamma* > input)
{
	art::ServiceHandle<util::DetectorProperties> detprop;
	art::ServiceHandle<geo::Geometry> geom;

	std::vector< std::vector< size_t > > saveids;
	std::vector< size_t > saveidsnotusedcls;
	size_t i = 0;
	
	while (i < input.size())	
	{
		if (!input[i]->GetCandidates().size()){i++; continue;}

		double mindist = 1.0; // cm 	
		std::vector< ems::DirOfGamma* > pairs; 

		size_t startview = input[i]->GetFirstHit()->WireID().Plane;
		size_t tpc = input[i]->GetFirstHit()->WireID().TPC;
		size_t cryo = input[i]->GetFirstHit()->WireID().Cryostat;

		float t1 = detprop->ConvertTicksToX(input[i]->GetFirstHit()->PeakTime(), startview, tpc, cryo);
		
		unsigned int idsave = 0; 
		for (unsigned int j = 0; j < input.size(); j++)
		{
			if (!input[j]->GetCandidates().size()) continue;

			size_t secondview = input[j]->GetFirstHit()->WireID().Plane;
			size_t tpc_j = input[j]->GetFirstHit()->WireID().TPC;
			size_t cryo_j = input[j]->GetFirstHit()->WireID().Cryostat;
		
			if ((i != j) && (secondview != startview) && (tpc == tpc_j) && (cryo == cryo_j))
			{				
				float t2 = detprop->ConvertTicksToX(input[j]->GetFirstHit()->PeakTime(), secondview, tpc_j, cryo_j);
				float dist = fabs(t2 - t1);
			
				if (dist < mindist)
				{
					mindist = dist;
					pairs.clear();
					pairs.push_back(input[i]); pairs.push_back(input[j]); 
					idsave = j;
				}
			}
		
		}

		bool exist = false;
		for (unsigned int v = 0; v < saveids.size(); v++)
				if ((saveids[v][0] == i) || (saveids[v][0] == idsave))
					if ((saveids[v][1] == i) || (saveids[v][1] == idsave)) 
						exist = true;


		if (pairs.size())
		{
			if (!exist) Make3DSeg(e, pairs);

			std::vector< size_t > ids;
			ids.push_back(i); ids.push_back(idsave);
			saveids.push_back(ids);
		}
		else
		{
			saveidsnotusedcls.push_back(i);
		}

		i++;
	}

	i = 0;
	while(i < saveidsnotusedcls.size())
	{
		LinkCandidates(e, input, i);
		i++;
	}
}

size_t ems::EMShower3D::LinkCandidates(art::Event const & e, std::vector< ems::DirOfGamma* > input, size_t id)
{
	art::ServiceHandle<geo::Geometry> geom;

	size_t index = id; bool found = false;

	if (input[id]->GetCandidates().size() < 2) { return index; }

	double mindist = 3.0; // cm 	
	std::vector< ems::DirOfGamma* > pairs; 

	size_t idcsave = 0; size_t idcjsave = 0;
	size_t c = 0; size_t idsave = 0; 
	while (c < input[id]->GetCandidates().size())	
	{

		size_t startview = input[id]->GetCandidates()[c].GetPlane();
		size_t tpc = input[id]->GetCandidates()[c].GetTPC();   
		size_t cryo = input[id]->GetCandidates()[c].GetCryo();

		float t1 = input[id]->GetCandidates()[c].GetPosition().Y(); // y --> drift in 2D space.
 
		// loop over 2D showers
		for (size_t j = 0; j < input.size(); ++j)
		{
			if (!input[j]->GetCandidates().size()) continue;
			if (j == id) continue;
			
			// loop over candidates
			for (size_t cj = 0; cj < input[j]->GetCandidates().size(); ++cj)
			{
				size_t secondview = input[j]->GetCandidates()[cj].GetPlane();
				size_t tpc_j = input[j]->GetCandidates()[cj].GetTPC();
				size_t cryo_j = input[j]->GetCandidates()[cj].GetCryo();

				size_t thirdview = startview;
				
				const geo::CryostatGeo& cryostat = geom->Cryostat(cryo);
				for (size_t p = 0; p < cryostat.MaxPlanes(); p++)
					if ((p == startview) || (p == secondview)) {continue;}
					else {thirdview = p; break;}
				
		
				if ((startview != secondview) && (tpc == tpc_j) && (cryo == cryo_j))// && Validate(input, id, cj, thirdview))
				{	
					float t2 = input[j]->GetCandidates()[cj].GetPosition().Y();
					float dist = fabs(t2 - t1);
				
					if ((dist < mindist) && Validate(input, id, j, c, cj, thirdview))
					{ 
						mindist = dist;
						pairs.clear();
						pairs.push_back(input[id]); pairs.push_back(input[j]);
						idsave = j; index = j;
						idcsave = c; idcjsave = cj;
						found = true;
					}
				}
			}
		}

		c++;
	}

	if (found && pairs.size())
	{
		input[id]->SetIdCandidate(idcsave);	
		input[idsave]->SetIdCandidate(idcjsave);
		Make3DSeg(e, pairs);
	}

	return index;
}

void ems::EMShower3D::Make3DSeg(art::Event const & e, std::vector< ems::DirOfGamma* > pair)
{	
	if (pair.size() < 2) return;

	// to build a track correctly 2d hits must belong to the same tpc	
	size_t tpc1 = pair[0]->GetFirstHit()->WireID().TPC;  
	size_t tpc2 = pair[1]->GetFirstHit()->WireID().TPC; 

	std::vector< art::Ptr< recob::Hit > > vec1 = pair[0]->GetIniHits();
	std::vector< art::Ptr< recob::Hit > > vec2 = pair[1]->GetIniHits();

	if ((vec1.size() < 3) && (vec2.size() < 3)) return;

	std::vector< art::Ptr<recob::Hit> > hitscl1uniquetpc;
	std::vector< art::Ptr<recob::Hit> > hitscl2uniquetpc;

	if (tpc1 == tpc2)
		for (size_t i = 0; i < vec1.size(); ++i)
			for (size_t j = 0; j < vec2.size(); ++j)
					if ((vec1[i]->WireID().TPC ==  vec2[j]->WireID().TPC) && (tpc1 == vec2[j]->WireID().TPC))
					{
						hitscl1uniquetpc.push_back(vec1[i]);
						hitscl2uniquetpc.push_back(vec2[j]);
					}

	if ((hitscl1uniquetpc.size() > 2) && (hitscl2uniquetpc.size() > 2))
	{
		pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(hitscl1uniquetpc, hitscl2uniquetpc);

		//pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(vec1, vec2);

		// turn the track that front is at vertex - easier to handle associations.
		if ((trk->back()->Hit2DPtr() == pair[0]->GetFirstHit()) 
			|| (trk->back()->Hit2DPtr() == pair[1]->GetFirstHit())) trk->Flip();	

	
		IniSeg initrack;
		initrack.idcl1 = pair[0]->GetIdCl();
		initrack.view1 = pair[0]->GetFirstHit()->WireID().Plane;
		initrack.hits1 = hitscl1uniquetpc;
		initrack.idcl2 = pair[1]->GetIdCl();
		initrack.view2 = pair[1]->GetFirstHit()->WireID().Plane;
		initrack.hits2 = hitscl2uniquetpc;
		initrack.track = trk;	
		
		fInisegs.push_back(initrack);

		
	}
}

/*void ems::EMShower3D::Make3DPMA(art::Event const & e, unsigned int cryo, unsigned int tpc)
{
	art::Handle< std::vector< recob::Track > > trkListHandle;	
	std::map< size_t, std::vector< size_t > > trkshs;

	if ((e.getByLabel(fTrk3DModuleLabel, trkListHandle)) && trkListHandle->size() &&
			 e.getByLabel(fCluModuleLabel, fCluListHandle))
	{
		art::FindManyP< recob::Hit > fbt(trkListHandle, e, fTrk3DModuleLabel);
		art::FindManyP< recob::Hit > fbc(fCluListHandle, e, fCluModuleLabel);		

		// chose unused cluster		
		for (size_t id = 0; id < fClustersNotUsed.size(); id++)
		{
			std::vector< art::Ptr<recob::Hit> > hitscl = fbc.at(fClustersNotUsed[id]);
			if (!hitscl.size()) continue;
			
				for (unsigned int t = 0; t < fTracksNotUsed.size(); t++)
				{
					std::vector< art::Ptr<recob::Hit> > hitstrk = fbt.at(fTracksNotUsed[t]);
					bool commonhits = false;
					for (unsigned int ht = 0; ht < hitstrk.size(); ht++)
					{
						if ((hitstrk[ht]->WireID().Cryostat == cryo) &&
								(hitstrk[ht]->WireID().TPC == tpc))
						{
							
							for (unsigned int c = 0; c < hitscl.size(); c++)
								if ((hitscl[c]->PeakTime() == hitstrk[ht]->PeakTime())
									 && (hitscl[c]->WireID().Wire == hitstrk[ht]->WireID().Wire) &&
											(hitscl[c]->WireID().Plane == hitstrk[ht]->WireID().Plane) &&
											(hitscl[c]->WireID().TPC == hitstrk[ht]->WireID().TPC) &&
											(hitscl[c]->WireID().Cryostat == hitstrk[ht]->WireID().Cryostat)) 
								{ 				
									commonhits = true; break; 
								}
						}

						if (commonhits)	
						{	
							trkshs[fClustersNotUsed[id]].push_back(fTracksNotUsed[t]);
							break;
						}		
					}
						
				}					
		}
		
		size_t ntracks = 1; 
		while (ntracks < 4)
		{
			auto it = trkshs.begin();
			while (it != trkshs.end())
			{	
				//size_t startsize = trkshs.size();
				
				if ((it->second.size() == ntracks) && fClustersNotUsed.size())
				{
					std::vector< art::Ptr<recob::Hit> > hitstrk = fbt.at(it->second[0]);
					if (!hitstrk.size()) continue;

					size_t plane1 = 0; size_t clid1 = it->first;
					size_t plane2 = plane1; size_t clid2 = clid1;
					size_t plane3 = plane1;	size_t clid3 = clid1;

					bool pl1 = false; bool pl2 = false; bool pl3 = false;

					for (size_t h = 0; h < hitstrk.size(); h++)
					{
						bool common = false;				
						std::vector< art::Ptr<recob::Hit> > hitscl = fbc.at(it->first);
						for (size_t hc = 0; hc < hitscl.size(); hc++)
							if ((hitscl[hc]->PeakTime() == hitstrk[h]->PeakTime()) &&
								(hitscl[hc]->WireID().Wire == hitstrk[h]->WireID().Wire) &&
								(hitscl[hc]->WireID().Plane == hitstrk[h]->WireID().Plane) &&
								(hitscl[hc]->WireID().TPC == hitstrk[h]->WireID().TPC) &&
								(hitscl[hc]->WireID().Cryostat == hitstrk[h]->WireID().Cryostat))
							{common = true; break;}

							if (common) {
							pl1 = true;
							clid1 = it->first; plane1 = hitstrk[h]->WireID().Plane; break;}		
					}
			
					for (size_t h = 0; h < hitstrk.size(); h++)
						if (pl1 && (hitstrk[h]->WireID().Plane != plane1)) 
						{
							bool common = false;
							for (size_t id = 0; id < fClustersNotUsed.size(); id++)
							{
								std::vector< art::Ptr<recob::Hit> > hitscl = fbc.at(fClustersNotUsed[id]);
								for (size_t hc = 0; hc < hitscl.size(); hc++)
									if ((hitscl[hc]->PeakTime() == hitstrk[h]->PeakTime()) &&
										(hitscl[hc]->WireID().Wire == hitstrk[h]->WireID().Wire) &&
										(hitscl[hc]->WireID().Plane == hitstrk[h]->WireID().Plane) &&
										(hitscl[hc]->WireID().TPC == hitstrk[h]->WireID().TPC) &&
										(hitscl[hc]->WireID().Cryostat == hitstrk[h]->WireID().Cryostat))
									{common = true; break;}

								if (common) {
								pl2 = true;
								clid2 = fClustersNotUsed[id]; plane2 = hitstrk[h]->WireID().Plane; break;}
							}
						
							if (common) break;
						}

					for (size_t h = 0; h < hitstrk.size(); h++)
						if (pl1 && pl2 && (hitstrk[h]->WireID().Plane != plane1) &&
								(hitstrk[h]->WireID().Plane != plane2))
						{
							bool common = false;
							for (size_t id = 0; id < fClustersNotUsed.size(); id++)
							{
								std::vector< art::Ptr<recob::Hit> > hitscl = fbc.at(fClustersNotUsed[id]);
								for (size_t hc = 0; hc < hitscl.size(); hc++)
									if ((hitscl[hc]->PeakTime() == hitstrk[h]->PeakTime()) &&
										(hitscl[hc]->WireID().Wire == hitstrk[h]->WireID().Wire) &&
										(hitscl[hc]->WireID().Plane == hitstrk[h]->WireID().Plane) &&
										(hitscl[hc]->WireID().TPC == hitstrk[h]->WireID().TPC) &&
										(hitscl[hc]->WireID().Cryostat == hitstrk[h]->WireID().Cryostat))
									{common = true; break;}

								if (common) {
									pl3 = true;
									clid3 = fClustersNotUsed[id]; plane3 = hitstrk[h]->WireID().Plane; break;}
							}
						
							if (common) break;
						}
					
						if (pl1 && pl2 && pl3)
						{
							if (trkshs[clid2].size() > trkshs[clid3].size()) 
							{
								clid2 = clid3;
								size_t temp = plane2;
								plane2 = plane3;
								plane3 = temp;
							}	

							std::vector< art::Ptr<recob::Hit> > hitscl1 = fbc.at(clid1);
							std::vector< art::Ptr<recob::Hit> > hitscl2 = fbc.at(clid2);

							std::vector< art::Ptr<recob::Hit> > hitscl1uniquetpc;
							std::vector< art::Ptr<recob::Hit> > hitscl2uniquetpc;
		
							for (size_t i = 0; i < hitscl1.size(); ++i)
								for (size_t j = 0; j < hitscl2.size(); ++j)
									if ((hitscl1[i]->WireID().TPC ==  hitscl2[j]->WireID().TPC) && (hitscl1[i]->WireID().TPC == tpc))
									{
										hitscl1uniquetpc.push_back(hitscl1[i]);
										hitscl2uniquetpc.push_back(hitscl2[j]);
									}

							if ((hitscl1uniquetpc.size() > 2) && (hitscl2uniquetpc.size() > 2))
							{						
							
								pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(hitscl1uniquetpc, hitscl2uniquetpc);

								if (Validate(e, *trk, plane3))
								{
									IniSeg initrack;
									initrack.idcl1 = clid1;
									initrack.view1 = plane1;
									initrack.hits1 = hitscl1uniquetpc;
									initrack.idcl2 = clid2;
									initrack.view2 = plane2;
									initrack.hits2 = hitscl2uniquetpc;
									initrack.track = trk;	
									fPMA3D.push_back(initrack);

									size_t cnt = 0; bool erase = false;
									while (cnt < fClustersNotUsed.size())
									{
										if (fClustersNotUsed[cnt] == clid1) 
										{
											fClustersNotUsed.erase(fClustersNotUsed.begin() + cnt);
											trkshs.erase(clid1); erase = true;	
										}
										else if (fClustersNotUsed[cnt] == clid2) 
										{
											fClustersNotUsed.erase(fClustersNotUsed.begin() + cnt);
											trkshs.erase(clid2); erase = true;
										}
										else cnt++;
									}

									if (erase) {it = trkshs.begin(); continue;}
								}
								else delete trk;
							}
						}
						else if (pl1 && pl2)
						{
							std::vector< art::Ptr<recob::Hit> > hitscl1 = fbc.at(clid1);
							std::vector< art::Ptr<recob::Hit> > hitscl2 = fbc.at(clid2);
							
							std::vector< art::Ptr<recob::Hit> > hitscl1uniquetpc;
							std::vector< art::Ptr<recob::Hit> > hitscl2uniquetpc;
		
							for (size_t i = 0; i < hitscl1.size(); ++i)
								for (size_t j = 0; j < hitscl2.size(); ++j)
									if ((hitscl1[i]->WireID().TPC ==  hitscl2[j]->WireID().TPC) && (hitscl1[i]->WireID().TPC == tpc))
									{
										hitscl1uniquetpc.push_back(hitscl1[i]);
										hitscl2uniquetpc.push_back(hitscl2[j]);
									}

							if ((hitscl1uniquetpc.size() > 2) && (hitscl2uniquetpc.size() > 2))
							{
								pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(hitscl1uniquetpc, hitscl2uniquetpc);
							
								art::ServiceHandle<geo::Geometry> geom;
								const geo::CryostatGeo& cryostat = geom->Cryostat(cryo);

								for (size_t p = 0; p < cryostat.MaxPlanes(); p++)
								{
									if ((plane1 == p) || (plane2 == p)) continue;
									else plane3 = p;
								}

								if (Validate(e, *trk, plane3))
								{
									IniSeg initrack;
									initrack.idcl1 = clid1;
									initrack.view1 = plane1;
									initrack.hits1 = hitscl1uniquetpc;
									initrack.idcl2 = clid2;
									initrack.view2 = plane2;
									initrack.hits2 = hitscl2uniquetpc;
									initrack.track = trk;	
									fPMA3D.push_back(initrack);
							
									size_t cnt = 0; bool erase = false;
									while (cnt < fClustersNotUsed.size())
									{
										if (fClustersNotUsed[cnt] == clid1) 
										{
											fClustersNotUsed.erase(fClustersNotUsed.begin() + cnt);
											trkshs.erase(clid1); erase = true;	
										}
										else if (fClustersNotUsed[cnt] == clid2) 
										{
											fClustersNotUsed.erase(fClustersNotUsed.begin() + cnt);
											trkshs.erase(clid2); erase = true;
										}
										else cnt++;
									}	
							

									if (erase) {it = trkshs.begin(); continue;}
								}
								else delete trk;
							}
						}
				}
				it++;
			}
			ntracks++;
		}
	}
}*/

bool ems::EMShower3D::Validate(std::vector< ems::DirOfGamma* > input, size_t id1, size_t id2, size_t c1, size_t c2, size_t plane3)
{
	bool result = false;
	if (id1 == id2) return false;

	std::vector< art::Ptr< recob::Hit > > vec1 = input[id1]->GetCandidates()[c1].GetIniHits();
	std::vector< art::Ptr< recob::Hit > > vec2 = input[id2]->GetCandidates()[c2].GetIniHits();

	if ((vec1.size() < 3) || (vec2.size() < 3)) return false;

	std::vector< art::Ptr<recob::Hit> > hitscl1uniquetpc;
	std::vector< art::Ptr<recob::Hit> > hitscl2uniquetpc;

	size_t tpc = vec1[0]->WireID().TPC;
	for (size_t i = 0; i < vec1.size(); ++i)
		for (size_t j = 0; j < vec2.size(); ++j)
			if ((vec1[i]->WireID().TPC ==  tpc) && (vec2[j]->WireID().TPC == tpc))
			{
				hitscl1uniquetpc.push_back(vec1[i]);
				hitscl2uniquetpc.push_back(vec2[j]);
			}
	
	if ((hitscl1uniquetpc.size() < 3) || (hitscl2uniquetpc.size() < 3)) return false;

	pma::Track3D* track = fProjectionMatchingAlg.buildSegment(hitscl1uniquetpc, hitscl2uniquetpc);
	for (size_t i = 0; i < input.size(); ++i)
	{	
		std::vector< Hit2D* > hits2dcl = input[i]->GetHits2D();
		for (size_t h = 0; h < hits2dcl.size(); ++h)
		{
			TVector2 pfront = pma::GetProjectionToPlane(track->front()->Point3D(), plane3, track->FrontTPC(), track->FrontCryo());
			TVector2 pback  = pma::GetProjectionToPlane(track->back()->Point3D(), plane3, track->BackTPC(), track->BackCryo());
			if ( (pma::Dist2(hits2dcl[h]->GetPointCm(), pfront) < 1.0F) && 
				(pma::Dist2(hits2dcl[h]->GetPointCm(), pback) < 1.0F) )
			{result = true; break;}
		}			
	}
	delete track;
	return result;
}



bool ems::EMShower3D::Validate(art::Event const & e, const pma::Track3D& src, size_t plane)
{
	bool result = false;

	art::FindManyP< recob::Hit > fbc(fCluListHandle, e, fCluModuleLabel);		
	std::vector< art::Ptr<recob::Hit> > hitscl;
	for (size_t id = 0; id < fClustersNotUsed.size(); id++)
	{
		std::vector< art::Ptr<recob::Hit> > hits = fbc.at(fClustersNotUsed[id]);
		for (size_t i = 0; i < hits.size(); i++) hitscl.push_back(hits[i]);
	}
	

	if (fProjectionMatchingAlg.validate(src, hitscl, plane) > 0.2) result = true;
		
	return result;
}

/*void ems::EMShower3D::SelectTrks(art::Event const & e)
{
	const float mindist2 = 1.0F; // 1 cm

	std::vector< size_t > trkids;
	art::Handle< std::vector< recob::Track > > trkListHandle;		
  if ((e.getByLabel(fTrk3DModuleLabel, trkListHandle)) && trkListHandle->size())
	{	
			for (unsigned int t = 0; t < trkListHandle->size(); t++)
			{
				const recob::Track & trk = (*trkListHandle)[t];

				double cosinemax = 0.8; bool found = false;
				unsigned int ids = 0;
				pma::Track3D* segsave;

				for (unsigned int s = 0; s < fInisegs.size(); s++)
				{
					pma::Track3D* seg = fInisegs[s].track;
					TVector3 segfront = seg->front()->Point3D();
					
					for (unsigned int v = 0; v < trk.NumberTrajectoryPoints(); v++)
					{
						float dist2 = pma::Dist2(segfront, trk.LocationAtPoint(v));		
					
						if (dist2 < mindist2)
						{
							TVector3 pos3d; TVector3 dir3d;
							trk.TrajectoryAtPoint(v, pos3d, dir3d);
							double inilength = (seg->back()->Point3D() - seg->front()->Point3D()).Mag();
							TVector3 dirseg = (seg->back()->Point3D() - seg->front()->Point3D()) * (1 / inilength);
							double cosine = dirseg * dir3d;
							if (fabs(cosine) > cosinemax)
							{
								cosinemax = fabs(cosine);
								found = true;
								ids = s;
								segsave = seg;
							}	
						}
					}
				}

				
				if (found)
				{
					double inilength = (segsave->back()->Point3D() - segsave->front()->Point3D()).Mag();
					if (inilength == 0) { fInisegs.erase(fInisegs.begin() + ids); continue;}

					fSeltracks.push_back(fInisegs[ids]);
					fInisegs.erase(fInisegs.begin() + ids);
					trkids.push_back(t);
				}
			}
	}

	for (size_t i = 0; i < trkids.size(); i++)
		MarkUsedTrks(trkids[i]);
	
	for (size_t i = 0; i < fSeltracks.size(); i++)
	{
		MarkUsedCls(fSeltracks[i].idcl1);
		MarkUsedCls(fSeltracks[i].idcl2);

		size_t plane1 = fSeltracks[i].view1;
		size_t plane2 = fSeltracks[i].view2;
		size_t plane3 = plane1;

		art::ServiceHandle<geo::Geometry> geom;

		for (unsigned int c = 0; c < geom->Ncryostats(); c++)
		{
			const geo::CryostatGeo& cryostat = geom->Cryostat(c);
			for (size_t p = 0; p < cryostat.MaxPlanes(); p++)
			{
				if ((plane1 == p) || (plane2 == p)) continue;
				else plane3 = p;
			}
		}

		// projection of track to elimate cluster from the 3rd view

		art::FindManyP< recob::Hit > fbc(fCluListHandle, e, fCluModuleLabel);
		bool has = false;  
		if (e.getByLabel(fCluModuleLabel, fCluListHandle))
			for (size_t c = 0; c < fCluListHandle->size(); c++)
			{
				std::vector< art::Ptr<recob::Hit> > hitlist;	
				hitlist = fbc.at(c);
				if (hitlist.size())
					for (size_t h = 0; h < hitlist.size(); h++)
					{
						if ((hitlist[h]->WireID().Plane == plane3) && (hitlist[h]->WireID().TPC == fSeltracks[i].track->FrontTPC()))
						{
							TVector2 projv = pma::GetProjectionToPlane(fSeltracks[i].track->front()->Point3D(), plane3, 
																							 					 fSeltracks[i].track->FrontTPC(), fSeltracks[i].track->FrontCryo());
							TVector2 wdcm = pma::WireDriftToCm (hitlist[h]->WireID().Wire, hitlist[h]->PeakTime(), plane3, hitlist[h]->WireID().TPC, hitlist[h]->WireID().Cryostat);
							if (pma::Dist2(wdcm, projv) < mindist2)
							{has = true; break;}
						}
						else if ((hitlist[h]->WireID().Plane == plane3) && (hitlist[h]->WireID().TPC == fSeltracks[i].track->BackTPC()))
						{
							TVector2 projv = pma::GetProjectionToPlane(fSeltracks[i].track->front()->Point3D(), plane3, 
																							 					 fSeltracks[i].track->BackTPC(), fSeltracks[i].track->BackCryo());
							TVector2 wdcm = pma::WireDriftToCm (hitlist[h]->WireID().Wire, hitlist[h]->PeakTime(), plane3, hitlist[h]->WireID().TPC, hitlist[h]->WireID().Cryostat);
							if (pma::Dist2(wdcm, projv) < mindist2)
							{has = true; break;}
						}		
					}

				if (has)
				{
					MarkUsedCls(c);
					break;
				}
				
			}
		

	}
}*/

void ems::EMShower3D::MarkUsedTrks(size_t idtrk)
{
	unsigned int cnt = 0;
	while (cnt < fTracksNotUsed.size())
	{
		bool has = false;

		if (idtrk == fTracksNotUsed[cnt])
		{
			fTracksNotUsed.erase(fTracksNotUsed.begin() + cnt);

			has = true;
		}

		if (has) continue;
		else cnt++;
	}
}

void ems::EMShower3D::MarkUsedCls(size_t idcl)
{
	unsigned int cnt = 0;
	while (cnt < fClustersNotUsed.size())
	{
		bool has = false;
		
		if (idcl == fClustersNotUsed[cnt]) 
		{
			fClustersNotUsed.erase(fClustersNotUsed.begin() + cnt);
				
			has = true;
		}

		if (has) continue;
		else cnt++;
	}
}

bool ems::EMShower3D::Has(const std::vector<size_t>& v, size_t idx)
{
    for (auto c : v) if (c == idx) return true;
    return false;
}

bool ems::EMShower3D::GetCloseHits(
		double r2d, 
		const std::vector< art::Ptr<recob::Hit> >& hits_in, 
		std::vector<size_t>& used,
		std::vector< art::Ptr<recob::Hit> >& hits_out)
{
	
	hits_out.clear();

	const double gapMargin = 5.0; // can be changed to f(id_tpc1, id_tpc2)
	size_t idx = 0;

	while ((idx < hits_in.size()) && Has(used, idx)) idx++;

	if (idx < hits_in.size())
	{		
		hits_out.push_back(hits_in[idx]);
		used.push_back(idx);
		

		double r2d2 = r2d*r2d;
		double gapMargin2 = sqrt(2 * gapMargin*gapMargin);
		gapMargin2 = (gapMargin2 + r2d)*(gapMargin2 + r2d);

		bool collect = true;
		while (collect)
		{
			collect = false;
			for (size_t i = 0; i < hits_in.size(); i++)
				if (!Has(used, i))
				{
					art::Ptr<recob::Hit> hi = hits_in[i];
					TVector2 hi_cm = pma::WireDriftToCm(hi->WireID().Wire, hi->PeakTime(), hi->WireID().Plane, hi->WireID().TPC, hi->WireID().Cryostat);

					bool accept = false;
					//for (auto const& ho : hits_out)
					for (size_t idx_o = 0; idx_o < hits_out.size(); idx_o++)					
					{
						art::Ptr<recob::Hit> ho = hits_out[idx_o];

						double d2 = pma::Dist2(
							hi_cm, pma::WireDriftToCm(ho->WireID().Wire, ho->PeakTime(), ho->WireID().Plane, ho->WireID().TPC, ho->WireID().Cryostat));
						
						if (hi->WireID().TPC == ho->WireID().TPC)
						{
							if (d2 < r2d2) { accept = true; break; }
						}
						else
						{
							if (d2 < gapMargin2) { accept = true; break; }
						}
					}
					if (accept)
					{
						collect = true;
						hits_out.push_back(hi);
						used.push_back(i);
					}
				}
		}
		return true;
	}
	else return false;
}


void ems::EMShower3D::FilterOutSmallParts(
		double r2d,
		const std::vector< art::Ptr<recob::Hit> >& hits_in,
		std::vector< art::Ptr<recob::Hit> >& hits_out)
{
	size_t min_size = hits_in.size() / 5;
	if (min_size < 3) min_size = 3;

	std::vector<size_t> used;
	std::vector< art::Ptr<recob::Hit> > close_hits;
	
	while (GetCloseHits(r2d, hits_in, used, close_hits))
	{
		if (close_hits.size() > min_size)
			for (auto h : close_hits) hits_out.push_back(h);
	}
}
DEFINE_ART_MODULE(ems::EMShower3D)
