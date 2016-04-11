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

#include "larreco/RecoAlg/PMAlg/PmaElement3D.h"
#include "larreco/RecoAlg/PMAlg/SortedObjects.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// Impact factors on the objective function:  U     V     Z
float pma::Element3D::fOptFactors[3] =     { 0.2F, 0.8F, 1.0F };

pma::Element3D::Element3D() :
	fTPC(-1), fCryo(-1),
	fFrozen(false),
	fHitsRadius(0)
{
	fNThisHitsEnabledAll = 0;
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

void pma::Element3D::SortHits(void)
{
	std::sort(fAssignedHits.begin(), fAssignedHits.end(), pma::bTrajectory3DOrderLess());
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
	for (size_t i = 0; i < 3; ++i) fNThisHitsEnabledAll = 0;
	for (auto h : fAssignedHits)
	{
		if (h->IsEnabled()) fNThisHitsEnabledAll++;
		switch (h->View2D())
		{
			case geo::kZ: hitsColl.push_back(h); break;
			case geo::kV: hitsInd2.push_back(h); break;
			case geo::kU: hitsInd1.push_back(h); break;
		}
	}
	fNThisHits[0] = hitsInd1.size();
	fNThisHits[1] = hitsInd2.size();
	fNThisHits[2] = hitsColl.size();

	pma::SortedObjectBase const * chain = dynamic_cast< pma::SortedObjectBase* >(this);
	pma::Element3D* el = 0;
	for (size_t b = 0; b < chain->NextCount(); b++)
	{
		el = dynamic_cast< pma::Element3D* >(chain->Next(b));
		if (el)
			for (auto h : el->fAssignedHits)
			{
				switch (h->View2D())
				{
					case geo::kZ: hitsColl.push_back(h); break;
					case geo::kV: hitsInd2.push_back(h); break;
					case geo::kU: hitsInd1.push_back(h); break;
				}
			}
	}
	el = dynamic_cast< pma::Element3D* >(chain->Prev());
	if (el)
	{
		for (auto h : el->fAssignedHits)
		{
			switch (h->View2D())
			{
				case geo::kZ: hitsColl.push_back(h); break;
				case geo::kV: hitsInd2.push_back(h); break;
				case geo::kU: hitsInd1.push_back(h); break;
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

double pma::Element3D::SumDist2(void) const
{
	if (fTPC < 0)
	{
		if (!fAssignedHits.empty()) mf::LogWarning("pma::Element3D") << "Hits assigned to TPC-crossing element.";
		return 0.0F;
	}

	double hit_sum = 0.0F;
	for (auto h : fAssignedHits)
	{
		if (h->IsEnabled())
		{
			unsigned int hitView = h->View2D();
/*			if (h->IsSpanEnabled())
			{
				double s = 0.0F;

				double d1 = sqrt(GetDistance2To(h->SpanPoint0(), hitView));
				double d2 = sqrt(GetDistance2To(h->SpanPoint2(), hitView));
				double d3 = sqrt(GetDistance2To(h->SpanPoint1(), hitView));
				double d, dx, dt;

				float driftSpan = h->SpanDrift();

				if (driftSpan <= 0.0F) mf::LogError("pma::Element3D") << "Hit drift-length is below zero.";

				dx = (d3 - d1) / 10.0;
				dt = (d2 - d1) / driftSpan;

				for (int x = 0; x < 10; x++)
				{
					d = d1 + x * dx;
					for (int t = 0; t < driftSpan; t++, d += dt)
						s += d * d;
				}

				hit_sum += s * h->GetSigmaFactor() * OptFactor(hitView);
			}
			else */
			hit_sum += OptFactor(hitView) *              // alpha_i
				h->GetSigmaFactor() *                    // hit_amp / hit_max_amp
				GetDistance2To(h->Point2D(), hitView);   // hit_to_fit_dist^2
		}
	}

	if (fAssignedPoints.size())
	{
		double d, ref_sum = 0.0F;
		for (auto p : fAssignedPoints)
		{
			d = sqrt( GetDistance2To(*p) ) - 0.5; // guide by ref points up to ~ 3D resolution
			if (d > 0.0) ref_sum += d * d;
		}
		if (fAssignedHits.size())
		{
			ref_sum *= 0.2 * fAssignedHits.size() / fAssignedPoints.size();
		}
		hit_sum += ref_sum;
	}

	return hit_sum;
}

double pma::Element3D::SumDist2(unsigned int view) const
{
	if (fTPC < 0)
	{
		if (!fAssignedHits.empty()) mf::LogWarning("pma::Element3D") << "Hits assigned to TPC-crossing element.";
		return 0.0F;
	}

	double hit_sum = 0.0F;
	for (auto h : fAssignedHits)
	{
		if (h->IsEnabled())
		{
			unsigned int hitView = h->View2D();
			if ((view == geo::kUnknown) || (view == hitView))
			{
				hit_sum += OptFactor(hitView) *              // alpha_i
					h->GetSigmaFactor() *                    // hit_amp / hit_max_amp
					GetDistance2To(h->Point2D(), hitView);   // hit_to_fit_dist^2
			}
		}
	}
	return hit_sum;
}

double pma::Element3D::HitsRadius3D(unsigned int view) const
{
	if (fTPC < 0)
	{
		if (!fAssignedHits.empty()) mf::LogWarning("pma::Element3D") << "Hits assigned to TPC-crossing element.";
		return 0.0F;
	}

	TVector3 mean3D(0, 0, 0);
	size_t nHits = 0;
	for (auto h : fAssignedHits)
		if (h->View2D() == view)
			{ mean3D += h->Point3D(); nHits++; }
	if (!nHits) return 0.0;
	mean3D *= (1.0 / nHits);

	double r2, maxR2 = 0.0;
	for (auto h : fAssignedHits)
		if (h->View2D() == view)
		{
			r2 = pma::Dist2(h->Point3D(), mean3D);
			if (r2 > maxR2) maxR2 = r2;
		}
	return sqrt(maxR2);
}

bool pma::Element3D::SelectRndHits(size_t nmax_per_view)
{
	if (!nmax_per_view) { return SelectAllHits(); }

	size_t nhits[3];
	for (size_t i = 0; i < 3; ++i) nhits[i] = NHits(i);

	int m[3], count[3];
	bool state[3];
	for (size_t i = 0; i < 3; ++i)
	{
		if (nhits[i] >= 2 * nmax_per_view)
		{
			m[i] = nhits[i] / nmax_per_view;
			state[i] = true;
		}
		else if (nhits[i] > nmax_per_view)
		{
			m[i] = nhits[i] / (nhits[i] - nmax_per_view);
			state[i] = false;
		}
		else { m[i] = 0; state[i] = false; }

		count[i] = 0;
	}

	bool b, changed = false;
	for (auto h : fAssignedHits)
	{
		b = h->IsEnabled();

		size_t view = h->View2D();
		if (m[view])
		{
			if (count[view] % m[view] == 0) h->SetEnabled(state[view]);
			else h->SetEnabled(!(state[view]));

			++count[view];
		}
		else h->SetEnabled(true);

		changed |= (b != h->IsEnabled());
	}
	return changed;
}

bool pma::Element3D::SelectAllHits(void)
{
	bool changed = false;
	for (auto h : fAssignedHits)
	{
		changed |= !(h->IsEnabled());
		h->SetEnabled(true);
	}
	return changed;
}

