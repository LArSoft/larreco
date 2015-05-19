/**
 *  @file   PmaElement3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Base for 3D segments and track nodes.
 */

#include "RecoAlg/PMAlg/PmaElement3D.h"

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
/*
	AF::HitSelection hitsColl, hitsInd1, hitsInd2;
	for (unsigned int i = 0; i < fAssignedHits.size(); i++)
	{
		switch (fAssignedHits[i]->View2D())
		{
			case T600::kViewColl: hitsColl.push_back(fAssignedHits[i]); break;
			case T600::kViewInd2: hitsInd2.push_back(fAssignedHits[i]); break;
			case T600::kViewInd1: hitsInd1.push_back(fAssignedHits[i]); break;
		}
	}
	fNThisHits[0] = hitsInd1.size();
	fNThisHits[1] = hitsInd2.size();
	fNThisHits[2] = hitsColl.size();

	SortedObjectBase* chain = dynamic_cast< SortedObjectBase* >(this);
	AF::PlElement3D* el = NULL;

	if (chain)
	{
		for (unsigned int b = 0; b < chain->NextCount(); b++)
		{
			el = dynamic_cast< AF::PlElement3D* >(chain->Next(b));
			if (el)
				for (unsigned int i = 0; i < el->fAssignedHits.size(); i++)
				{
					switch (el->fAssignedHits[i]->View2D())
					{
						case T600::kViewColl: hitsColl.push_back(el->fAssignedHits[i]); break;
						case T600::kViewInd2: hitsInd2.push_back(el->fAssignedHits[i]); break;
						case T600::kViewInd1: hitsInd1.push_back(el->fAssignedHits[i]); break;
					}
				}
		}
	}

	if (chain) el = dynamic_cast< AF::PlElement3D* >(chain->Prev());
	else el = NULL;
	
	if (el)
	{
		for (unsigned int i = 0; i < el->fAssignedHits.size(); i++)
		{
			switch (el->fAssignedHits[i]->View2D())
			{
				case T600::kViewColl: hitsColl.push_back(el->fAssignedHits[i]); break;
				case T600::kViewInd2: hitsInd2.push_back(el->fAssignedHits[i]); break;
				case T600::kViewInd1: hitsInd1.push_back(el->fAssignedHits[i]); break;
			}
		}
	}

	fHitsRadius = hitsColl.Radius();
	float r = hitsInd2.Radius();
	if (r > fHitsRadius) fHitsRadius = r;
	r = hitsInd1.Radius();
	if (r > fHitsRadius) fHitsRadius = r;

	float amp, sigmaMax = 0.0F;
	fSumHitsQ[0] = 0.0; fNHits[0] = hitsInd1.size();
	for (unsigned int i = 0; i < hitsInd1.size(); i++)
	{
		amp = ((PlHit3D*)hitsInd1[i])->GetAmplitude();
		if (amp > sigmaMax) sigmaMax = amp;
		fSumHitsQ[0] += amp;
	}
	for (unsigned int i = 0; i < hitsInd1.size(); i++)
	{
		if (sigmaMax > 0.0F)
		{
			amp = ((PlHit3D*)hitsInd1[i])->GetAmplitude();
			if (amp > 0.0F)
				((PlHit3D*)hitsInd1[i])->SetSigmaFactor((float)sqrt(amp / sigmaMax));
			else ((PlHit3D*)hitsInd1[i])->SetSigmaFactor(0.01F);
		}
		else ((PlHit3D*)hitsInd1[i])->SetSigmaFactor(1.0F);
	}

	sigmaMax = 0.0F;
	fSumHitsQ[1] = 0.0; fNHits[1] = hitsInd2.size();
	for (unsigned int i = 0; i < hitsInd2.size(); i++)
	{
		amp = ((PlHit3D*)hitsInd2[i])->GetAmplitude();
		if (amp > sigmaMax) sigmaMax = amp;
		fSumHitsQ[1] += amp;
	}
	for (unsigned int i = 0; i < hitsInd2.size(); i++)
	{
		if (sigmaMax > 0.0F)
		{
			amp = ((PlHit3D*)hitsInd2[i])->GetAmplitude();
			if (amp > 0.0F)
				((PlHit3D*)hitsInd2[i])->SetSigmaFactor((float)sqrt(amp / sigmaMax));
			else ((PlHit3D*)hitsInd2[i])->SetSigmaFactor(0.01F);
		}
		else ((PlHit3D*)hitsInd2[i])->SetSigmaFactor(1.0F);
	}

	sigmaMax = 0.0F;
	for (unsigned int i = 0; i < hitsColl.size(); i++)
	{
		amp = ((PlHit3D*)hitsColl[i])->GetEnergyDeposit();
		if (amp > sigmaMax) sigmaMax = amp;
	}
	for (unsigned int i = 0; i < hitsColl.size(); i++)
	{
		if (sigmaMax > 0.0F)
		{
			amp = ((PlHit3D*)hitsColl[i])->GetEnergyDeposit();
			if (amp > 0.0F)
				((PlHit3D*)hitsColl[i])->SetSigmaFactor((float)sqrt(amp / sigmaMax));
			else ((PlHit3D*)hitsColl[i])->SetSigmaFactor(0.01F);
		}
		else ((PlHit3D*)hitsColl[i])->SetSigmaFactor(1.0F);
	}
	fSumHitsQ[2] = hitsColl.GetEnergySum();
	fNHits[2] = hitsColl.size();

	hitsColl.release_all();
	hitsInd2.release_all();
	hitsInd1.release_all();
*/
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

