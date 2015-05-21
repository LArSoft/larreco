/**
 *  @file   PmaElement3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Base for 3D segments and track nodes. See PmaTrack3D.h file for details.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaElement3D.h"
#include "RecoAlg/PMAlg/SortedObjects.h"

// Impact factors on the objective function:  U     V     Z
float pma::Element3D::fOptFactors[3] =     { 0.2F, 0.8F, 1.0F };

pma::Element3D::Element3D(void) :
	fFrozen(false), fHitsRadius(0)
{
	for (unsigned int i = 0; i < 3; i++)
	{
		fNHits[i] = 0;
		fNThisHits[i] = 0;
		fSumHitsQ[i] = 0.0;
	}
}

size_t pma::Element3D::NEnabledHits(unsigned int view) const
{
	size_t n = 0;
	for (size_t i = 0; i < fAssignedHits.size(); i++)
		if (fAssignedHits[i]->IsEnabled() &&
		    ((view == geo::kUnknown) || (view == fAssignedHits[i]->View2D()))) n++;
	return n;
}

void pma::Element3D::ClearAssigned(pma::Track3D* trk)
{
	fAssignedPoints.clear();
	fAssignedHits.clear();
	fHitsRadius = 0.0;
}

void pma::Element3D::UpdateHitParams(void)
{
	std::vector< pma::Hit3D* > hitsColl, hitsInd1, hitsInd2;
	for (size_t i = 0; i < fAssignedHits.size(); i++)
	{
		switch (fAssignedHits[i]->View2D())
		{
			case geo::kZ: hitsColl.push_back(fAssignedHits[i]); break;
			case geo::kV: hitsInd2.push_back(fAssignedHits[i]); break;
			case geo::kU: hitsInd1.push_back(fAssignedHits[i]); break;
		}
	}
	fNThisHits[0] = hitsInd1.size();
	fNThisHits[1] = hitsInd2.size();
	fNThisHits[2] = hitsColl.size();

	pma::SortedObjectBase* chain = dynamic_cast< pma::SortedObjectBase* >(this);
	pma::Element3D* el = NULL;

	if (chain)
	{
		for (size_t b = 0; b < chain->NextCount(); b++)
		{
			el = dynamic_cast< pma::Element3D* >(chain->Next(b));
			if (el)
				for (size_t i = 0; i < el->fAssignedHits.size(); i++)
				{
					switch (el->fAssignedHits[i]->View2D())
					{
						case geo::kZ: hitsColl.push_back(el->fAssignedHits[i]); break;
						case geo::kV: hitsInd2.push_back(el->fAssignedHits[i]); break;
						case geo::kU: hitsInd1.push_back(el->fAssignedHits[i]); break;
					}
				}
		}
	}

	if (chain) el = dynamic_cast< pma::Element3D* >(chain->Prev());
	else el = NULL;

	if (el)
	{
		for (size_t i = 0; i < el->fAssignedHits.size(); i++)
		{
			switch (el->fAssignedHits[i]->View2D())
			{
				case geo::kZ: hitsColl.push_back(el->fAssignedHits[i]); break;
				case geo::kV: hitsInd2.push_back(el->fAssignedHits[i]); break;
				case geo::kU: hitsInd1.push_back(el->fAssignedHits[i]); break;
			}
		}
	}

	fHitsRadius = GetHitsRadius2D(hitsColl);
	double r = GetHitsRadius2D(hitsInd2);
	if (r > fHitsRadius) fHitsRadius = r;
	r = GetHitsRadius2D(hitsInd1);
	if (r > fHitsRadius) fHitsRadius = r;

	float amp, sigmaMax = 0.0F;
	fSumHitsQ[0] = 0.0; fNHits[0] = hitsInd1.size();
	for (size_t i = 0; i < hitsInd1.size(); i++)
	{
		amp = hitsInd1[i]->GetAmplitude();
		if (amp > sigmaMax) sigmaMax = amp;
		fSumHitsQ[0] += amp;
	}
	for (size_t i = 0; i < hitsInd1.size(); i++)
	{
		if (sigmaMax > 0.0F)
		{
			amp = hitsInd1[i]->GetAmplitude();
			if (amp > 0.0F)
				hitsInd1[i]->SetSigmaFactor((float)sqrt(amp / sigmaMax));
			else hitsInd1[i]->SetSigmaFactor(0.01F);
		}
		else hitsInd1[i]->SetSigmaFactor(1.0F);
	}

	sigmaMax = 0.0F;
	fSumHitsQ[1] = 0.0; fNHits[1] = hitsInd2.size();
	for (size_t i = 0; i < hitsInd2.size(); i++)
	{
		amp = hitsInd2[i]->GetAmplitude();
		if (amp > sigmaMax) sigmaMax = amp;
		fSumHitsQ[1] += amp;
	}
	for (size_t i = 0; i < hitsInd2.size(); i++)
	{
		if (sigmaMax > 0.0F)
		{
			amp = hitsInd2[i]->GetAmplitude();
			if (amp > 0.0F)
				hitsInd2[i]->SetSigmaFactor((float)sqrt(amp / sigmaMax));
			else hitsInd2[i]->SetSigmaFactor(0.01F);
		}
		else hitsInd2[i]->SetSigmaFactor(1.0F);
	}

	sigmaMax = 0.0F;
	fSumHitsQ[2] = 0.0; fNHits[2] = hitsColl.size();
	for (size_t i = 0; i < hitsColl.size(); i++)
	{
		amp = hitsColl[i]->SummedADC();
		if (amp > sigmaMax) sigmaMax = amp;
		fSumHitsQ[2] += amp;
	}
	for (size_t i = 0; i < hitsColl.size(); i++)
	{
		if (sigmaMax > 0.0F)
		{
			amp = hitsColl[i]->SummedADC();
			if (amp > 0.0F)
				hitsColl[i]->SetSigmaFactor((float)sqrt(amp / sigmaMax));
			else hitsColl[i]->SetSigmaFactor(0.01F);
		}
		else hitsColl[i]->SetSigmaFactor(1.0F);
	}
}

void pma::Element3D::UpdateProjection(void)
{
	for (size_t i = 0; i < fAssignedHits.size(); i++)
	{
		SetProjection(*(fAssignedHits[i]));
	}
}

void pma::Element3D::SortHits(void)
{
	std::sort(fAssignedHits.begin(), fAssignedHits.end(), pma::bTrajectory3DOrderLess());
}

double pma::Element3D::SumDist2(void) const
{
	double hit_sum = 0.0F;
	for (size_t i = 0; i < fAssignedHits.size(); i++)
	{
		if (fAssignedHits[i]->IsEnabled())
		{
			unsigned int hitView = fAssignedHits[i]->View2D();
/*			if (fAssignedHits[i]->IsSpanEnabled())
			{
				double s = 0.0F;

				double d1 = sqrt(GetDistance2To(fAssignedHits[i]->SpanPoint0(), hitView));
				double d2 = sqrt(GetDistance2To(fAssignedHits[i]->SpanPoint2(), hitView));
				double d3 = sqrt(GetDistance2To(fAssignedHits[i]->SpanPoint1(), hitView));
				double d, dx, dt;

				float driftSpan = fAssignedHits[i]->SpanDrift();

				if (driftSpan <= 0.0F) LOG_ERROR("Hit drift-length is <= zero.\n");

				dx = (d3 - d1) / 10.0;
				dt = (d2 - d1) / driftSpan;

				for (int x = 0; x < 10; x++)
				{
					d = d1 + x * dx;
					for (int t = 0; t < driftSpan; t++, d += dt)
						s += d * d;
				}

				hit_sum += s * fAssignedHits[i]->GetSigmaFactor() * OptFactor(hitView);
			}
			else */
			hit_sum += OptFactor(hitView) *                             // alpha_i
				fAssignedHits[i]->GetSigmaFactor() *                    // hit_amp / hit_max_amp
				GetDistance2To(fAssignedHits[i]->Point2D(), hitView);   // hit_to_fit_dist^2
		}
	}

	double ref_sum = 0.0F;
	if (fAssignedPoints.size())
	{
		double factor = 0.2 * fAssignedHits.size() / fAssignedPoints.size();
		for (size_t i = 0; i < fAssignedPoints.size(); i++)
		{
			ref_sum += GetDistance2To(*(fAssignedPoints[i]));
		}
		if (fAssignedHits.size()) ref_sum *= factor;
	}

	return hit_sum + ref_sum;
}

double pma::Element3D::HitsRadius3D(unsigned int view) const
{
	TVector3 mean3D(0, 0, 0);
	size_t nHits = 0;
	for (size_t i = 0; i < fAssignedHits.size(); i++)
		if (fAssignedHits[i]->View2D() == view)
			{ mean3D += fAssignedHits[i]->Point3D(); nHits++; }
	if (!nHits) return 0.0;
	mean3D *= (1.0 / nHits);

	double r2, maxR2 = 0.0;
	for (size_t i = 0; i < fAssignedHits.size(); i++)
		if (fAssignedHits[i]->View2D() == view)
		{
			r2 = pma::Dist2(fAssignedHits[i]->Point3D(), mean3D);
			if (r2 > maxR2) maxR2 = r2;
		}
	return sqrt(maxR2);
}

