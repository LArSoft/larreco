/**
 *  @file   SimpleClustering.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Collect hits "touching" each other (next wire or consecutive ticks).
 */

#include "SimpleClustering.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoAlg/PMAlg/Utilities.h"

tss::Cluster2D::Cluster2D(const std::vector< const tss::Hit2D* > & hits) :
	fStartIdx(0), fEndIdx(0)
{
	fHits.reserve(hits.size());
	for (size_t h = 0; h < hits.size(); ++h) fHits.push_back(hits[h]);
}
// ------------------------------------------------------

const tss::Hit2D* tss::Cluster2D::release(size_t idx)
{
	const tss::Hit2D* hit = 0;
	if (idx < fHits.size())
	{
		hit = fHits[idx];
		fHits.erase(fHits.begin() + idx);
	}
	return hit;
}
// ------------------------------------------------------

const tss::Hit2D* tss::Cluster2D::closest(const TVector2 & p2d, size_t & idx) const
{
	idx = 0;
	if (!fHits.size()) return 0;

	const tss::Hit2D* hout = fHits.front();
	double d, dmin = pma::Dist2(hout->Point2D(), p2d);
	for (size_t h = 1; h < fHits.size(); ++h)
	{
		d = pma::Dist2(fHits[h]->Point2D(), p2d);
		if (d < dmin) { dmin = d; hout = fHits[h]; idx = h; }
	}
	return hout;
}
// ------------------------------------------------------

const tss::Hit2D* tss::Cluster2D::outermost(size_t & idx) const
{
	idx = 0;
	if (!fHits.size()) return 0;

	TVector2 mean(0., 0.);
	for (size_t h = 0; h < fHits.size(); ++h)
	{
		mean += fHits[h]->Point2D();
	}
	mean *= 1.0 / fHits.size();

	const tss::Hit2D* hout = fHits.front();
	double d, dmax = pma::Dist2(hout->Point2D(), mean);
	for (size_t h = 1; h < fHits.size(); ++h)
	{
		d = pma::Dist2(fHits[h]->Point2D(), mean);
		if (d > dmax) { dmax = d; hout = fHits[h]; idx = h; }
	}
	return hout;
}
// ------------------------------------------------------

bool tss::Cluster2D::Has(const tss::Hit2D* hit) const
{
	for (size_t i = 0; i < fHits.size(); ++i)
		if (fHits[i] == hit) return true;
	return false;
}
// ------------------------------------------------------

double tss::Cluster2D::dist2(const TVector2 & p2d) const
{
	if (fHits.size())
	{
		double d2, min_d2 = pma::Dist2(fHits.front()->Point2D(), p2d);
		for (size_t i = 1; i < fHits.size(); ++i)
		{
			d2 = pma::Dist2(fHits[i]->Point2D(), p2d);
			if (d2 < min_d2) { min_d2 = d2; }
		}
		return min_d2;
	}
	else return 0.0;
}
// ------------------------------------------------------

double tss::Cluster2D::dist2(const TVector2 & p2d, size_t & hIdx) const
{
	hIdx = 0;
	if (fHits.size())
	{
		double d2, min_d2 = pma::Dist2(fHits.front()->Point2D(), p2d);
		for (size_t i = 1; i < fHits.size(); ++i)
		{
			d2 = pma::Dist2(fHits[i]->Point2D(), p2d);
			if (d2 < min_d2) { min_d2 = d2; hIdx = i; }
		}
		return min_d2;
	}
	else return 0.0;
}
// ------------------------------------------------------

double tss::Cluster2D::dist2(const tss::Cluster2D & clu) const
{
	if (fHits.size())
	{
		double d2, min_d2 = clu.dist2(fHits.front()->Point2D());
		for (size_t i = 1; i < fHits.size(); ++i)
		{
			d2 = clu.dist2(fHits[i]->Point2D());
			if (d2 < min_d2) { min_d2 = d2; }
		}
		return min_d2;
	}
	else return 0.0;
}
// ------------------------------------------------------

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


bool tss::SimpleClustering::hitsTouching(const tss::Hit2D & h1, const tss::Hit2D & h2) const
{
	if ((h1.Wire() == h2.Wire()) &&
	    (h1.PeakTime() == h2.PeakTime())) return false;

	bool touches = false;
	if ((h1.Wire() >= h2.Wire() - 1) &&
	    (h1.Wire() <= h2.Wire() + 1))
	{
		if (((h2.StartTick() <= h1.StartTick()) && (h1.StartTick() <= h2.EndTick() + 1)) ||
			((h2.StartTick() <= h1.EndTick() + 1) && (h1.EndTick() <= h2.EndTick())) ||
			((h2.StartTick() >= h1.StartTick()) && (h1.EndTick() >= h2.EndTick())))
		{
			touches = true;
		}
	}
	return touches;
}
// ------------------------------------------------------

bool tss::SimpleClustering::hitsTouching(const tss::Cluster2D & c1, const tss::Hit2D & h2) const
{
	for (size_t i = 0; i < c1.size(); i++)
	{
		if (hitsTouching(c1[i], h2)) return true;	
	}
	return false;
}
// ------------------------------------------------------

bool tss::SimpleClustering::hitsTouching(const tss::Cluster2D & c1, const tss::Cluster2D & c2) const
{
	for (unsigned int i = 0; i < c1.size(); i++)
	{
		if (hitsTouching(c2, c1[i])) return true;
	}
	return false;
}
// ------------------------------------------------------

void tss::SimpleClustering::merge(std::vector< tss::Cluster2D > & clusters) const
{
	bool merged = true;
	while (merged)
	{
		merged = false;

		size_t i = 0;
		while (i < clusters.size() - 1)
		{
			size_t j = i + 1;
			while (j < clusters.size())
			{
				if (hitsTouching(clusters[i], clusters[j]))
				{
					clusters[i].hits().reserve(clusters[i].size() + clusters[i].size());
					for (size_t h = 0; h < clusters[j].size(); ++h)
						clusters[i].hits().push_back(clusters[j].hits()[h]);
					clusters.erase(clusters.begin() + j);
					merged = true;
				}
				else ++j;
			}
			++i;
		}
	}
}
// ------------------------------------------------------

std::vector< tss::Cluster2D > tss::SimpleClustering::run(const std::vector< tss::Hit2D > & inp)
{
	std::vector< tss::Cluster2D > result;
	for (size_t h = 0; h < inp.size(); ++h)
	{
		bool found = false;
		for (size_t r = 0; r < result.size(); ++r)
			if (hitsTouching(result[r], inp[h]))
		{
			result[r].hits().push_back(&(inp[h])); found = true; break;
		}
		if (!found)
		{
			result.push_back(tss::Cluster2D());
			result.back().hits().push_back(&(inp[h]));
		}
	}
	merge(result);

	return result;
}
// ------------------------------------------------------

std::vector< tss::Cluster2D > tss::SimpleClustering::run(const tss::Cluster2D & inp)
{
	std::vector< tss::Cluster2D > result;
	for (size_t h = 0; h < inp.size(); ++h)
	{
		bool found = false;
		for (size_t r = 0; r < result.size(); ++r)
			if (hitsTouching(result[r], inp[h]))
		{
			result[r].hits().push_back(inp.hits()[h]); found = true; break;
		}
		if (!found)
		{
			result.push_back(tss::Cluster2D());
			result.back().hits().push_back(inp.hits()[h]);
		}
	}
	merge(result);

	return result;
}
// ------------------------------------------------------

