/**
 *  @file   SimpleClustering.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Collect hits "touching" each other (next wire or consecutive ticks).
 */

#include "SimpleClustering.h"

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

bool tss::SimpleClustering::hitsTouching(const std::vector< const tss::Hit2D* > & c1, const tss::Hit2D & h2) const
{
	for (size_t i = 0; i < c1.size(); i++)
	{
		if (hitsTouching(*(c1[i]), h2)) return true;	
	}
	return false;
}
// ------------------------------------------------------

bool tss::SimpleClustering::hitsTouching(const std::vector< const tss::Hit2D * > & c1, const std::vector< const tss::Hit2D* > & c2) const
{
	for (unsigned int i = 0; i < c1.size(); i++)
	{
		if (hitsTouching(c2, *(c1[i]))) return true;
	}
	return false;
}
// ------------------------------------------------------

void tss::SimpleClustering::merge(std::vector< std::vector< const tss::Hit2D* > > & clusters) const
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
					clusters[i].reserve(clusters[i].size() + clusters[i].size());
					for (size_t h = 0; h < clusters[j].size(); ++h)
						clusters[i].push_back(clusters[j][h]);
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

std::vector< std::vector< const tss::Hit2D* > > tss::SimpleClustering::run(const std::vector< tss::Hit2D > & inp)
{
	std::vector< std::vector< const tss::Hit2D* > > result;
	for (size_t h = 0; h < inp.size(); ++h)
	{
		bool found = false;
		for (size_t r = 0; r < result.size(); ++r)
			if (hitsTouching(result[r], inp[h]))
		{
			result[r].push_back(&(inp[h])); found = true; break;
		}
		if (!found)
		{
			result.push_back(std::vector< const tss::Hit2D* >());
			result.back().push_back(&(inp[h]));
		}
	}
	merge(result);

	return result;
}
// ------------------------------------------------------

std::vector< std::vector< const tss::Hit2D* > > tss::SimpleClustering::run(const std::vector< const tss::Hit2D* > & inp)
{
	std::vector< std::vector< const tss::Hit2D* > > result;
	for (size_t h = 0; h < inp.size(); ++h)
	{
		bool found = false;
		for (size_t r = 0; r < result.size(); ++r)
			if (hitsTouching(result[r], *(inp[h])))
		{
			result[r].push_back(inp[h]); found = true; break;
		}
		if (!found)
		{
			result.push_back(std::vector< const tss::Hit2D* >());
			result.back().push_back(inp[h]);
		}
	}
	merge(result);

	return result;
}
// ------------------------------------------------------

