/**
 *  @file   Segmentation2D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Split into linear clusters.
 */

#include "Segmentation2D.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoAlg/PMAlg/Utilities.h"

std::vector< tss::Cluster2D > tss::Segmentation2D::run(tss::Cluster2D & inp) const
{
	std::vector< tss::Cluster2D > result;
	while (inp.size() > 1)
	{
		size_t idx;
		const tss::Hit2D* hFirst = inp.outermost(idx);
		if (!hFirst) break;

		std::vector< TVector2 > centers;
		centers.emplace_back(hFirst->Point2D());

		while (centers.size())
		{
			run(inp, result, centers);
		}
	}
	return result;
}
// ------------------------------------------------------

void tss::Segmentation2D::run(
	tss::Cluster2D & inp,
	std::vector< tss::Cluster2D > & result,
	std::vector< TVector2 > & centers) const
{
	if (!centers.size()) return;

	TVector2 center(centers.front());
	centers.erase(centers.begin());

	size_t idx;
	const tss::Hit2D* hFirst = inp.closest(center, idx);

	const double dmax2 = 0.5 * 0.5; // does not look like startpoint selected before
	if (!hFirst || (pma::Dist2(hFirst->Point2D(), center) > dmax2)) return;
	center = hFirst->Point2D();

	tss::Cluster2D ring = selectRing(inp, center);
	if (ring.size())
	{
		std::vector< tss::Cluster2D > seeds = fSimpleClustering.run(ring);
		while (seeds.size())
		{
			size_t seedIdx = 0, hitIdx, h;
			double d2, min_d2 = seeds.front().dist2(center, hitIdx);
			for (size_t i = 1; i < seeds.size(); i++)
			{
				d2 = seeds[i].dist2(center, h);
				if (d2 < min_d2) { min_d2 = d2; seedIdx = i; hitIdx = h; }
			}

			tss::Cluster2D segment = buildSegment(inp, center, seeds[seedIdx][hitIdx].Point2D());
			if (segment.size())
			{
				result.emplace_back(segment);
				if (segment.size() > 1)
				{
					const tss::Hit2D* hEnd = segment.end();
					if (hEnd) centers.emplace_back(hEnd->Point2D());
				}
			}

			seeds.erase(seeds.begin() + seedIdx);
		}
	}
	else inp.release(hFirst);
}
// ------------------------------------------------------

tss::Cluster2D tss::Segmentation2D::buildSegment(tss::Cluster2D & inp, TVector2 center, TVector2 end) const
{
	const double max_d2 = fMaxLineDist * fMaxLineDist;
	TVector2 segDir = end - center;

	double dc, min_dc = 1.0e9;
	size_t firstIdx = 0;

	tss::Cluster2D candidates;
	for (auto h : inp.hits())
	{
		TVector2 proj = pma::GetProjectionToSegment(h->Point2D(), center, end);
		if (pma::Dist2(h->Point2D(), proj) < max_d2)
		{
			TVector2 hDir = h->Point2D() - center;
			dc = hDir.Mod();
			if ((hDir * segDir >= 0.0) || (dc < 0.1))
			{
				candidates.push_back(h);
				if (dc < min_dc)
				{
					min_dc = dc; firstIdx = candidates.size() - 1;
				}
			}
		}
	}
	if (candidates.size() > 1)
	{
		const tss::Hit2D* hFirst = candidates.hits()[firstIdx];
		candidates.hits()[firstIdx] = candidates.hits()[0];
		candidates.hits()[0] = hFirst;
		candidates.sort();
	}

	tss::Cluster2D segment;
	if (candidates.size())
	{
		segment.push_back(candidates.start());
		if (!inp.release(segment.start()))
		{
			mf::LogError("Segmentation2D") << "Hit not found in the input cluster.";
		}

		size_t i = 1;
		while ((i < candidates.size()) &&
		       fSimpleClustering.hitsTouching(segment, candidates[i]))
		{
			segment.push_back(candidates.hits()[i++]);
			if (!inp.release(segment.end()))
			{
				mf::LogError("Segmentation2D") << "Hit not found in the input cluster.";
			}
		}
	}

	return segment;
}
// ------------------------------------------------------

tss::Cluster2D tss::Segmentation2D::selectRing(const tss::Cluster2D & inp, TVector2 center) const
{
	double d2_min = fRadiusMin * fRadiusMin;
	double d2_max = fRadiusMax * fRadiusMax;

	tss::Cluster2D ring;
	for (size_t h = 0; h < inp.size(); h++)
	{
		double d2 = pma::Dist2(center, inp[h].Point2D());
		if ((d2 >= d2_min) && (d2 <= d2_max))
			ring.push_back(inp.hits()[h]);
	}
	return ring;
}
// ------------------------------------------------------

