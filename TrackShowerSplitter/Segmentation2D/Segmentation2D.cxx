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

	tagDenseEnds(result);
	mergeDenseParts(result);

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
	else
	{
		result.emplace_back(tss::Cluster2D());
		result.back().push_back(hFirst);
		inp.release(hFirst);
	}
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

void tss::Segmentation2D::tagDenseEnds(std::vector< tss::Cluster2D > & group) const
{
	const double rad2 = fDenseVtxRadius * fDenseVtxRadius;

	for (size_t i = 0; i < group.size(); i++)
	{
		bool denseStart = false, denseEnd = false;
		TVector2 start0(group[i].start()->Point2D()), end0(group[i].end()->Point2D());

		for (size_t j = 0; j < group.size(); j++)
		{
			if (i == j) continue;

			TVector2 start1(group[j].start()->Point2D()), end1(group[j].end()->Point2D());

			if (!group[j].isDenseStart())
			{
				if (pma::Dist2(start1, start0) < rad2)
				{
					group[j].tagDenseStart(true);
					denseStart = true;
				}
				if (pma::Dist2(start1, end0) < rad2)
				{
					group[j].tagDenseStart(true);
					denseEnd = true;
				}
			}

			if (!group[j].isDenseEnd())
			{
				if (pma::Dist2(end1, start0) < rad2)
				{
					group[j].tagDenseEnd(true);
					denseStart = true;
				}
				if (pma::Dist2(end1, end0) < rad2)
				{
					group[j].tagDenseEnd(true);
					denseEnd = true;
				}
			}
		}
					
		if (denseStart) group[i].tagDenseStart(true);
		if (denseEnd) group[i].tagDenseEnd(true);
	}
}
// ------------------------------------------------------

void tss::Segmentation2D::mergeDenseParts(std::vector< tss::Cluster2D > & group) const
{
	const double rad2 = fDenseVtxRadius * fDenseVtxRadius;

	bool merged = true;
	while (merged)
	{
		merged = false;

		size_t maxS = fDenseMinN, maxE = fDenseMinN;
		std::vector< size_t > toMergeS, toMergeE;
		int idxMaxS = -1, idxMaxE = -1;
		for (size_t i = 0; i < group.size(); i++)
		{
			if (group[i].isEM()) continue;

			if (group[i].isDenseStart())
			{
				size_t ns = 0;
				std::vector< size_t > toMerge;
				TVector2 start0(group[i].start()->Point2D());
				for (size_t j = 0; j < group.size(); j++)
				{
					if (group[j].isEM()) continue;

					if (i == j)
					{
						if ((group[j].size() > 1) &&
						    (group[j].length2() < rad2)) ns++;
					}
					else
					{
						bool tagged = false;
						if (pma::Dist2(start0, group[j].start()->Point2D()) < rad2)
						{
							ns++; toMerge.push_back(j); tagged = true;
						}
						if ((group[j].size() > 1) && (pma::Dist2(start0, group[j].end()->Point2D()) < rad2))
						{
							ns++; if (!tagged) toMerge.push_back(j);
						}
					}
				}
				if (ns > maxS) { maxS = ns; idxMaxS = i; toMergeS = toMerge; }
			}
			if ((group[i].size() > 1) && group[i].isDenseEnd())
			{
				size_t ne = 0;
				std::vector< size_t > toMerge;
				TVector2 end0(group[i].end()->Point2D());
				for (size_t j = 0; j < group.size(); j++)
				{
					if (group[j].isEM()) continue;

					if (i == j)
					{
						if ((group[j].size() > 1) &&
						    (group[j].length2() < rad2)) ne++;
					}
					else
					{
						bool tagged = false;
						if (pma::Dist2(end0, group[j].start()->Point2D()) < rad2)
						{
							ne++; toMerge.push_back(j); tagged = true;
						}
						if ((group[j].size() > 1) && (pma::Dist2(end0, group[j].end()->Point2D()) < rad2))
						{
							ne++; if (!tagged) toMerge.push_back(j);
						}
					}
				}
				if (ne > maxE) { maxE = ne; idxMaxE = i; toMergeE = toMerge; }
			}
		}

		int idx = idxMaxS;
		std::vector< size_t > toMergeIdxs = toMergeS;
		if (idxMaxE > idx) { idx = idxMaxE; toMergeIdxs = toMergeE; }
		if (idx > -1)
		{
			for (size_t i = 0; i < toMergeIdxs.size(); i++) // *** no merging, only tag instead ***
			{
				group[toMergeIdxs[i]].tagEM(true);
			}
			group[idx].tagEM(true);

			merged = true;
		}
	}
}
// ------------------------------------------------------

void tss::Segmentation2D::splitHits(
		const std::vector< tss::Cluster2D > & inp,
		std::vector< const tss::Hit2D* > & trackHits,
		std::vector< const tss::Hit2D* > & emHits) const
{
	trackHits.clear();
	emHits.clear();

	for (const auto & cx : inp)
	{
		if (!cx.size()) continue;

		if (cx.isEM())
		{
			for (auto h : cx.hits()) emHits.push_back(h);
		}
		else
		{
			for (auto h : cx.hits()) trackHits.push_back(h);
		}
	}
}
// ------------------------------------------------------

void tss::Segmentation2D::splitHitsNaive(
		const std::vector< tss::Cluster2D > & inp,
		std::vector< const tss::Hit2D* > & trackHits,
		std::vector< const tss::Hit2D* > & emHits) const
{
	const double rad2 = fDenseVtxRadius * fDenseVtxRadius;

	trackHits.clear();
	emHits.clear();

	for (const auto & cx : inp)
	{
		if (!cx.size()) continue;

		for (const auto hx : cx.hits())
		{
			size_t n = 0;
			for (const auto & cy : inp)
			{
				if (!cy.size()) continue;

				for (const auto hy : cy.hits())
				{
					if (hx->Hit2DPtr() == hy->Hit2DPtr()) continue;

					if (pma::Dist2(hx->Point2D(), hy->Point2D()) < rad2) n++;
				}
			}

			if (n > fDenseMinH)
			{
				emHits.push_back(hx);
			}
			else
			{
				trackHits.push_back(hx);
			}
		}
	}
}
// ------------------------------------------------------


