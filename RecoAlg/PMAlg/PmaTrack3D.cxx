/**
 *  @file   PmaTrack3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaTrack3D.h"

pma::Track3D::Track3D(void) :
	fMaxHitsPerSeg(70),
	fPenaltyFactor(2.0F),
	fMaxSegStopFactor(3.0F),

	fSegStopValue(2), fMinSegStop(2), fMaxSegStop(2),
	fSegStopFactor(0.2F),
	fPenaltyValue(0.1F),
	fEndSegWeight(0.05F),
	fHitsRadius(1.0F)
{
}

pma::Track3D::~Track3D(void)
{
	for (size_t i = 0; i < fHits.size(); i++) delete fHits[i];
	for (size_t i = 0; i < fAssignedPoints.size(); i++) delete fAssignedPoints[i];

	for (size_t i = 0; i < fSegments.size(); i++) delete fSegments[i];
	for (size_t i = 0; i < fNodes.size(); i++)
		if (!fNodes[i]->NextCount() && !fNodes[i]->Prev()) delete fNodes[i];
}

bool pma::Track3D::push_back(const recob::Hit& hit)
{
	pma::Hit3D* h3d = new pma::Hit3D(hit);
	for (size_t i = 0; i < fHits.size(); i++)
	{
		pma::Hit3D* hit_i = fHits[i];
		if ((h3d->PeakTime() == hit_i->PeakTime()) &&
		    (h3d->Wire() == hit_i->Wire()) &&
		    (h3d->View2D() == hit_i->View2D()) &&
		    (h3d->TPC() == hit_i->TPC()))
		{
			delete h3d; return false;
		}
	}
	fHits.push_back(h3d);
	return true;
}

void pma::Track3D::AddHits(const std::vector< recob::Hit const* >& hits)
{
	for (size_t i = 0; i < hits.size(); i++) push_back(*(hits[i]));
}

unsigned int pma::Track3D::NHits(unsigned int view) const
{
	unsigned int n = 0;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = static_cast< pma::Hit3D* >((*this)[i]);
		if (hit && (hit->View2D() == view)) n++;
	}
	return n;
}

unsigned int pma::Track3D::NEnabledHits(unsigned int view) const
{
	unsigned int n = 0;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = static_cast< pma::Hit3D* >((*this)[i]);
		if (hit && hit->IsEnabled() &&
		    ((view == geo::kUnknown) || (view == hit->View2D()))) n++;
	}
	return n;
}

pma::Segment3D* pma::Track3D::NextSegment(pma::Node3D* vtx) const
{
	pma::Segment3D* seg = 0;
	unsigned int nCount = vtx->NextCount();
	unsigned int k = 0;
	while (k < nCount)
	{
		seg = static_cast< pma::Segment3D* >(vtx->Next(k));
		if (seg && (seg->Parent() == this)) return seg;
		k++;
	}
	return 0;
}

pma::Segment3D* pma::Track3D::PrevSegment(pma::Node3D* vtx) const
{
	if (vtx->Prev())
	{
		pma::Segment3D* seg = static_cast< pma::Segment3D* >(vtx->Prev());
		if (seg->Parent() == this) return seg;
	}
	return 0;
}

bool pma::Track3D::HasRefPoint(TVector3* p) const
{
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
		if (fAssignedPoints[i] == p) return true;
	return false;
}

double pma::Track3D::GetObjFunction(bool suppressPenalty) const
{
	double sum = 0.0;
	float p = fPenaltyValue;
	if (suppressPenalty) p *= 0.33F;
	for (size_t i = 0; i < fNodes.size(); i++)
	{
		sum += fNodes[i]->GetObjFunction(p, fEndSegWeight);
	}
	return sum / fNodes.size();
}

void pma::Track3D::RebuildSegments(void)
{
	for (size_t i = 0; i < fSegments.size(); i++) delete fSegments[i];
	fSegments.clear();

	for (size_t i = 1; i < fNodes.size(); i++)
	{
		pma::Segment3D* s = new pma::Segment3D(this, fNodes[i - 1], fNodes[i]);
		fSegments.push_back(s);
	}
}

pma::Element3D* pma::Track3D::GetNearestElement(
	const TVector2& p2d, unsigned int view) const
{
	pma::Element3D* pe_min = fNodes.front();
	double dist, min_dist = pe_min->GetDistance2To(p2d, view);
	for (size_t i = 1; i < fNodes.size(); i++)
	{
		dist = fNodes[i]->GetDistance2To(p2d, view);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fNodes[i];
		}
	}
	for (size_t i = 0; i < fSegments.size(); i++)
	{
		dist = fSegments[i]->GetDistance2To(p2d, view);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fSegments[i];
		}
	}
	return pe_min;
}

pma::Element3D* pma::Track3D::GetNearestElement(const TVector3& p3d) const
{
	pma::Element3D* pe_min = fNodes.front();
	double dist, min_dist = pe_min->GetDistance2To(p3d);
	for (size_t i = 1; i < fNodes.size(); i++)
	{
		dist = fNodes[i]->GetDistance2To(p3d);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fNodes[i];
		}
	}
	for (size_t i = 0; i < fSegments.size(); i++)
	{
		dist = fSegments[i]->GetDistance2To(p3d);
		if (dist < min_dist)
		{
			min_dist = dist; pe_min = fSegments[i];
		}
	}
	return pe_min;
}

void pma::Track3D::SortHits(void)
{
	std::vector< pma::Hit3D* > hits_tmp;
	pma::Node3D* vtx = fNodes.front();
	pma::Segment3D* seg = NextSegment(vtx);
	while (vtx)
	{
		vtx->SortHits();
		for (size_t i = 0; i < vtx->NHits(); i++)
		{
			pma::Hit3D* h3d = &(vtx->Hit(i));
			for (size_t j = 0; j < size(); j++)
				if ((*this)[j] == h3d)
				{
					hits_tmp.push_back(h3d);
					break;
				}
		}

		if (seg)
		{
			seg->SortHits();
			for (size_t i = 0; i < seg->NHits(); i++)
			{
				pma::Hit3D* h3d = &(seg->Hit(i));
				for (size_t j = 0; j < size(); j++)
					if ((*this)[j] == h3d)
					{
						hits_tmp.push_back(h3d);
						break;
					}
			}

			vtx = static_cast< pma::Node3D* >(seg->Next());
			seg = NextSegment(vtx);
		}
		else break;
	}

	if (size() == hits_tmp.size())
	{
		for (size_t i = 0; i < size(); i++) (*this)[i] = hits_tmp[i];
	}
	else std::cout << "Hit sorting problem." << std::endl;
}

void pma::Track3D::MakeProjection(void)
{
	for (size_t i = 0; i < fNodes.size(); i++) fNodes[i]->ClearAssigned(this);
	for (size_t i = 0; i < fSegments.size(); i++) fSegments[i]->ClearAssigned(this);

	pma::Element3D* pe = 0;

	for (size_t i = 0; i < size(); i++) // assign hits to nodes/segments
	{
		pma::Hit3D* hit = static_cast< pma::Hit3D* >((*this)[i]);
		pe = GetNearestElement(hit->Point2D(), hit->View2D());
		pe->AddHit(hit);
	}

	for (size_t i = 0; i < fAssignedPoints.size(); i++) // assign ref points to nodes/segments
	{
		pe = GetNearestElement(*(fAssignedPoints[i]));
		pe->AddPoint(fAssignedPoints[i]);
	}

	// move hits to segments if not branching
	if (!(fNodes.front()->Prev()) && (fNodes.front()->NextCount() == 1))
	{
		for (size_t i = 0; i < fNodes.front()->NHits(); i++)
			fSegments.front()->AddHit(&(fNodes.front()->Hit(i)));
		fNodes.front()->ClearAssigned();
	}

	if (fNodes.back()->NextCount() == 0)
	{
		for (size_t i = 0; i < fNodes.back()->NHits(); i++)
			fSegments.back()->AddHit(&(fNodes.back()->Hit(i)));
		fNodes.back()->ClearAssigned();
	}

	for (unsigned int i = 0; i < fNodes.size(); i++) fNodes[i]->UpdateHitParams();
	for (unsigned int i = 0; i < fSegments.size(); i++) fSegments[i]->UpdateHitParams();
}

void pma::Track3D::UpdateProjection(void)
{
	for (size_t i = 0; i < fNodes.size(); i++) fNodes[i]->UpdateProjection();
	for (size_t i = 0; i < fSegments.size(); i++) fSegments[i]->UpdateProjection();
}

double pma::Track3D::AverageDist2(void) const
{
	double sum = 0.0;
	unsigned int count = 0;
	pma::Node3D* vtx = fNodes.front();
	pma::Segment3D* seg = NextSegment(vtx);
	while (vtx)
	{
		sum += vtx->SumDist2();
		count += vtx->NEnabledHits();
		if (seg)
		{
			sum += seg->SumDist2();
			count += seg->NEnabledHits();
			vtx = static_cast< pma::Node3D* >(seg->Next());
			seg = NextSegment(vtx);
		}
		else break;
	}
	return sum / count;
}

void pma::Track3D::UpdateParams(void)
{
	size_t n = size();
	if (!n) n = 1;

	float nCubeRoot = pow((double)n, 1.0/3.0);
	float avgDist2Root = sqrt(AverageDist2());	

	fPenaltyValue = fPenaltyFactor * pow((double)fSegments.size(), 1.8) * avgDist2Root / (fHitsRadius * nCubeRoot);

	fSegStopValue = (int)(fSegStopFactor * nCubeRoot * fHitsRadius / avgDist2Root);
	if (fSegStopValue < fMinSegStop) fSegStopValue = fMinSegStop;
}

bool pma::Track3D::SwapVertices(size_t v0, size_t v1)
{
	if (v0 == v1) return false;

	if (v0 > v1)
	{
		size_t vx = v0;
		v0 = v1; v1 = vx;
	}

	pma::Node3D* vtmp;
	if (v1 - v0 == 1)
	{
		pma::Segment3D* midSeg = NextSegment(fNodes[v0]);
		pma::Segment3D* prevSeg = PrevSegment(fNodes[v0]);
		pma::Segment3D* nextSeg = NextSegment(fNodes[v1]);

		fNodes[v1]->RemoveNext(nextSeg);
		midSeg->Disconnect();

		vtmp = fNodes[v0];
		fNodes[v0] = fNodes[v1];
		fNodes[v1] = vtmp;

		if (prevSeg) prevSeg->AddNext(fNodes[v0]);
		fNodes[v0]->AddNext(midSeg);
		midSeg->AddNext(fNodes[v1]);
		if (nextSeg) fNodes[v1]->AddNext(nextSeg);

		return false;
	}
	else
	{
		vtmp = fNodes[v0];
		fNodes[v0] = fNodes[v1];
		fNodes[v1] = vtmp;
		return true;
	}
}

bool pma::Track3D::CheckEndSegment(pma::Track3D::ETrackEnd endCode)
{
	unsigned int v1, v2;
	switch (endCode)
	{
		case pma::Track3D::kBegin:
			if (fSegments.front()->IsFrozen()) return false;
			if (fNodes.front()->NextCount() > 1) return false;
			v1 = 0; v2 = 1; break;
		case pma::Track3D::kEnd:
			if (fSegments.back()->IsFrozen()) return false;
			if (fNodes.back()->NextCount() > 1) return false;
			v1 = fNodes.size() - 1;
			v2 = fNodes.size() - 2;
			break;
		default: return false;
	}

	double g1, g0 = GetObjFunction();

	if (SwapVertices(v1, v2)) RebuildSegments();
	MakeProjection();
	g1 = GetObjFunction();

	if (g1 >= g0)
	{
		if (SwapVertices(v1, v2)) RebuildSegments();
		MakeProjection();
		return false;
	}
	else return true;
}

void pma::Track3D::UpdateHitsRadius(void)
{
	std::vector< pma::Hit3D* > hitsColl, hitsInd1, hitsInd2;
	for (size_t i = 0; i < size(); i++)
	{
		pma::Hit3D* hit = static_cast< pma::Hit3D* >((*this)[i]);
		switch (hit->View2D())
		{
			case geo::kZ: hitsColl.push_back(hit); break;
			case geo::kV: hitsInd2.push_back(hit); break;
			case geo::kU: hitsInd1.push_back(hit); break;
		}
	}
	fHitsRadius = pma::GetHitsRadius2D(hitsColl, true);
	float r = pma::GetHitsRadius2D(hitsInd2, true);
	if (r > fHitsRadius) fHitsRadius = r;
	r = pma::GetHitsRadius2D(hitsInd1, true);
	if (r > fHitsRadius) fHitsRadius = r;
}

