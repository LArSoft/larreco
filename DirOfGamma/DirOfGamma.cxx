#include "DirOfGamma.h"

#include "RecoAlg/PMAlg/PmaHit3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include "Utilities/DetectorProperties.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

//class Hit2D

ems::Hit2D::Hit2D(art::Ptr< recob::Hit > src) :
fHit(src)
{
	art::ServiceHandle<geo::Geometry> geom;
	art::ServiceHandle<util::DetectorProperties> detprop;

	unsigned int plane = src->WireID().Plane;
	unsigned int tpc   = src->WireID().TPC;
	unsigned int cryo  = src->WireID().Cryostat;
	
	double wireCentre[3];
	geom->WireIDToWireGeo(src->WireID()).GetCenter(wireCentre);
	double x = detprop->ConvertTicksToX(src->PeakTime(), plane, tpc, cryo);
	double globalWire;
	
	if (tpc % 2 == 0) 
	{	
		globalWire = geom->WireCoordinate(wireCentre[1], wireCentre[2], plane, 0, cryo);
		fPoint.Set(globalWire, x);
	}
	else 
	{
		globalWire = geom->WireCoordinate(wireCentre[1], wireCentre[2], plane, 1, cryo);
		fPoint.Set(globalWire, x);
	}
	fCharge = src->SummedADC();
}

ems::Bin2D::Bin2D(const TVector2 & center) : 
fCenter2D(center), 
fTotCharge(0.0), 
fSize(0)
{
}

void ems::Bin2D::Add(Hit2D* hit)
{
	fHits2D.push_back(hit);
	fTotCharge += hit->GetCharge();
	fSize = fHits2D.size();
	SortLess();
}

void ems::Bin2D::Sort()
{
	return std::sort(fHits2D.begin(), fHits2D.end(), bDistCentMore2D(fCenter2D));
}

void ems::Bin2D::SortLess()
{
	return std::sort(fHits2D.begin(), fHits2D.end(), bDistCentLess2D(fCenter2D)); 
}

std::vector< art::Ptr< recob::Hit > > ems::Bin2D::GetIniHits(const double radius, const unsigned int nhits) const
{

	std::vector< art::Ptr< recob::Hit > > vec;
	for (unsigned int i = 0; i < fHits2D.size(); i++)
	{
		if (pma::Dist2(fHits2D[i]->GetPointCm(), fCenter2D) < radius*radius)
		{
			 vec.push_back(fHits2D[i]->GetHitPtr());
			if (vec.size() == nhits) break;
		}
	}

	return vec;
}

ems::EndPoint::EndPoint(const Hit2D & center, const std::vector< Hit2D* > & hits, unsigned int nbins) : 
fCenter2D(center), 
fPoints2D(hits),
fNbins(nbins)
{
	art::ServiceHandle<geo::Geometry> geom;
	art::ServiceHandle<util::DetectorProperties> detprop;
	
	for (unsigned int i = 0; i < fNbins; i++)
	{
		fBins.push_back(Bin2D(center.GetPointCm()));
	}

	FillBins();
	ComputeMaxCharge();
	ComputeMeanCharge();

	fPlane = center.GetHitPtr()->WireID().Plane;
	fTpc   = center.GetHitPtr()->WireID().TPC;
	fCryo  = center.GetHitPtr()->WireID().Cryostat;
}

void ems::EndPoint::FillBins() 
{
	TVector2 vstart(0, 1);

	unsigned int saveid = 0; bool exist = false;
	for (unsigned int i = 0; i < fPoints2D.size(); i++)
	{
		if (fPoints2D[i]->GetHitPtr().key() != fCenter2D.GetHitPtr().key())
		{
			TVector2 pos(fPoints2D[i]->GetPointCm());
			TVector2 centre(fCenter2D.GetPointCm());
			TVector2 vecp = pos - centre;
			float sinsign = (vstart.X() * vecp.Y()) - (vstart.Y() * vecp.X());
			float cosine  = (vstart * vecp) / vecp.Mod();
			float theta = 180.0F * (std::acos(cosine)) / TMath::Pi();

			unsigned int id = 0; double bin = double(360.0) / double(fNbins);

			if (sinsign >= 0) id = int (theta / bin); 
			else if (sinsign < 0) id = int (theta / bin) + (fNbins / 2); 
			if (id > (fNbins - 1)) id = (fNbins - 1);		

				fBins[id].Add(fPoints2D[i]); 
				fBins[(id + 1) % fNbins].Add(fPoints2D[i]);
		}
		else {saveid = i; exist = true;}
	}

	if (exist) 
		for (unsigned int id = 0; id < fNbins; id++) fBins[id].Add(fPoints2D[saveid]);	
}

void ems::EndPoint::ComputeMaxCharge()
{
	fMaxCharge = 0.0;
	unsigned int saveid = 0; 
	for (unsigned int i = 0; i < fNbins; i++)
		if (fBins[i].Size() && (fMaxCharge < fBins[i].GetTotCharge()))
		{
			fMaxCharge = fBins[i].GetTotCharge();
			saveid = i; 
		}

	fMaxChargeIdBin = saveid; 
}

void ems::EndPoint::ComputeMeanCharge()
{
	fMeanCharge = 0.0;
	if (fNbins == 0) return;

	unsigned int iprev, inext;

	if (fMaxChargeIdBin > 0) iprev = fMaxChargeIdBin - 1;
	else iprev = fNbins - 1;

	inext = (fMaxChargeIdBin + 1) % fNbins;

	double sumcharge = 0.0;
	for (unsigned int i = 0; i < fNbins; i++)
		if ((i!=fMaxChargeIdBin) && (i!=iprev) && (i!=inext)) sumcharge += fBins[i].GetTotCharge();

	fMeanCharge = sumcharge / double (fNbins);
}

double ems::EndPoint::GetAsymmetry() const
{
	if ((fMaxCharge + fMeanCharge) == 0) return 0.0;
	else return ((fMaxCharge - fMeanCharge) / (fMaxCharge + fMeanCharge));
}

ems::DirOfGamma::DirOfGamma(const std::vector< art::Ptr< recob::Hit > > & src, unsigned int nbins, unsigned int idcl) :
fNbins(nbins),
fIdCl(idcl),
fCandidateID(0),
fIsCandidateIDset(false)
{
	fHits = src;

	for (unsigned int i = 0; i < src.size(); i++)
	{
		Hit2D* hit = new Hit2D(src[i]);
		fPoints2D.push_back(hit);
	}

	ComputeBaryCenter();

	for (unsigned int i = 0; i < fNbins; i++)
		fBins.push_back(Bin2D(fBaryCenter));

	FillBins();
	ComputeMaxDist();
	if (FindCandidates())
	{
		ComputeMaxCharge();
		FindInitialPart();
	}
}

void ems::DirOfGamma::ComputeBaryCenter()
{
	double nomx = 0.0; double nomy = 0.0;
	double denom = 0.0;
	for (unsigned int i = 0; i < fPoints2D.size(); i++)
	{
		nomx  += fPoints2D[i]->GetPointCm().X() * fPoints2D[i]->GetCharge();
		nomy  += fPoints2D[i]->GetPointCm().Y() * fPoints2D[i]->GetCharge();
		denom += fPoints2D[i]->GetCharge();
	}

	double bx = nomx / denom; double by = nomy / denom;
	fBaryCenter.Set(bx, by);
}

void ems::DirOfGamma::FillBins()
{
	TVector2 vstart(0, 1);	

	for (unsigned int i = 0; i < fPoints2D.size(); i++)	
	{
		TVector2 pos(fPoints2D[i]->GetPointCm().X(), fPoints2D[i]->GetPointCm().Y());
		TVector2 vecp = pos - fBaryCenter;
		float sinsign = (vstart.X() * vecp.Y()) - (vstart.Y() * vecp.X());
		float cosine  = (vstart * vecp) / (vstart.Mod() * vecp.Mod());
		float theta = 180.0F * (std::acos(cosine)) / TMath::Pi();

		unsigned int id = 0; double bin = double(360.0) / double(fNbins);

		if (sinsign >= 0) id = int (theta / bin); 
		else if (sinsign < 0) id = int (theta / bin) + (fNbins/2); 
		if (id > (fNbins - 1)) id = (fNbins - 1);		

		fBins[id].Add(fPoints2D[i]); 
	}

	for (unsigned int id = 0; id < fBins.size(); id++) fBins[id].Sort();
	
}

void ems::DirOfGamma::ComputeMaxDist()
{
	double maxdist2 = 0.0;
	
	for (unsigned int id = 0; id < fBins.size(); id++)
	{
		
		if (!fBins[id].Size()) continue;
		
		Hit2D* candidate = fBins[id].GetHits2D().front();
		if (candidate)
		{
			double dist2 = pma::Dist2(candidate->GetPointCm(), fBaryCenter);
			if (dist2 > maxdist2) 
			{
				maxdist2 = dist2;
			}
		}
	}
	
	fNormDist = std::sqrt(maxdist2);
}

bool ems::DirOfGamma::FindCandidates()
{
	float rad = 0.5F * fNormDist; unsigned int nbins = fNbins * 4; 
	for (unsigned int id = 0; id < fNbins; id++)
	{
		
		if (!fBins[id].Size()) continue;

		std::vector< Hit2D* > points;
		Hit2D* candidate2D = fBins[id].GetHits2D().front(); 	

		for (unsigned int i = 0; i < fPoints2D.size(); i++)
		{
			double distnorm = std::sqrt(pma::Dist2(candidate2D->GetPointCm(), fBaryCenter)) / fNormDist;
			double dist2 = pma::Dist2(candidate2D->GetPointCm(), fPoints2D[i]->GetPointCm());

			if ((distnorm > 0.5) && (dist2 < rad*rad)) points.push_back(fPoints2D[i]);
		}

		
		if (fBins[id].Size() > 1)
		{
			EndPoint ep(*candidate2D, points, nbins); 
			fCandidates.push_back(ep);
		}
	}
	if (fCandidates.size()) return true;
	else return false;
}

void ems::DirOfGamma::ComputeMaxCharge()
{
	fNormCharge = 0.0;
	for (unsigned int i = 0; i < fCandidates.size(); i++)
	{
		if (fCandidates[i].GetMaxCharge() > fNormCharge) 
		{
			fNormCharge = fCandidates[i].GetMaxCharge();
		}
	}
}

void ems::DirOfGamma::FindInitialPart()
{
	double max_asymmetry = 0.0; 
	unsigned int saveid = 0; bool found = false;

	double maxdist2 = 0.0; double maxcharge = 0.0;
	unsigned int idmaxdist = 0; unsigned int idmaxcharge = 0;

	for (unsigned int i = 0; i < fCandidates.size(); i++)
	{
		double dist2 = pma::Dist2(fCandidates[i].GetPosition(), fBaryCenter);
		double charge = fCandidates[i].GetMaxCharge();
		if (dist2 > maxdist2) { maxdist2 = dist2; idmaxdist = i;}
		if (charge > maxcharge) { maxcharge = charge; idmaxcharge = i;}
	}

	maxdist2 = 0.0; unsigned int idmaxdistb = 0;
	for (size_t i = 0; i < fCandidates.size(); i++)
	{
		if ((i == idmaxdist) || (i == idmaxcharge)) continue;
		
		double dist2 = pma::Dist2(fCandidates[i].GetPosition(), fCandidates[idmaxdist].GetPosition());
		if (dist2 > maxdist2) { maxdist2 = dist2; idmaxdistb = i;}		
	}	

	if (fCandidates.size() > 2)
	{
		for (unsigned int i = 0; i < fCandidates.size(); i++)
		{
			double asymmetry = fCandidates[i].GetAsymmetry();

			if ((i == idmaxdist) || (i == idmaxcharge) || (i == idmaxdistb))
			{
				if (asymmetry > max_asymmetry)
				{			
					max_asymmetry = asymmetry;
					saveid = i; found = true;
				}
			}	
		}
	}
	else
	{
		for (unsigned int i = 0; i < fCandidates.size(); i++)
		{
			double asymmetry = fCandidates[i].GetAsymmetry();

			if ((i == idmaxdist) || (i == idmaxdistb))
			{
				if (asymmetry > max_asymmetry)
				{			
					max_asymmetry = asymmetry;
					saveid = i; found = true;
				}
			}	
		}
	}

	if (!found) mf::LogError("DirOfGamma") << fCandidates.size() << "DirOfGamma - Find Initial Part problem.";
	
	fStartHit = fCandidates[saveid].GetHit(); 
	fStartPoint = fCandidates[saveid].GetPosition();
	fIniHits  = fCandidates[saveid].MaxChargeBin().GetIniHits();
	fCandidateID = saveid;
	
}







