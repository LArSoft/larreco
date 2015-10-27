/**
 *  @file   SimpleClustering.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Trivial, collect hits "touching" each other (next wire or consecutive ticks).
 */

#ifndef SimpleClustering_h
#define SimpleCustering_h

#include "TssHit2D.h"

namespace tss
{
	class Cluster2D;
	class SimpleClustering;
}

class tss::Cluster2D
{
public:

	Cluster2D(void) : fStartIdx(0), fEndIdx(0) {}
	Cluster2D(const std::vector< const tss::Hit2D* > & hits);

	size_t size(void) const { return fHits.size(); }

	const tss::Hit2D & operator [] (size_t index) const { return *(fHits[index]); }

	const std::vector< const tss::Hit2D* > & hits(void) const { return fHits; }
	std::vector< const tss::Hit2D* > & hits(void) { return fHits; }

	bool Has(const tss::Hit2D* hit) const;

	double dist2(const TVector2 & p2d) const;
	double dist2(const TVector2 & p2d, size_t & hIdx) const;
	double dist2(const tss::Cluster2D & clu) const;

	const tss::Hit2D* release(size_t idx);
	void push_back(const tss::Hit2D* hit) { fHits.push_back(hit); }
	void take_from(tss::Cluster2D & clu, size_t idx)
	{
		const tss::Hit2D* hit = clu.release(idx);
		if (hit) push_back(hit);
	}

	const tss::Hit2D* start(void) const
	{
		if (fStartIdx < fHits.size()) return fHits[fStartIdx];
		else return 0;
	}
	void set_start(size_t idx) { fStartIdx = idx; }

	const tss::Hit2D* end(void) const
	{
		if (fEndIdx < fHits.size()) return fHits[fEndIdx];
		else return 0;
	}

	const tss::Hit2D* closest(const TVector2 & p2d, size_t & idx) const;
	const tss::Hit2D* outermost(size_t & idx) const;

private:

	std::vector< const tss::Hit2D* > fHits;

	size_t fStartIdx, fEndIdx;

};

class tss::SimpleClustering
{
public:

	std::vector< tss::Cluster2D > run(const std::vector< tss::Hit2D > & inp);
	std::vector< tss::Cluster2D > run(const tss::Cluster2D & inp);

private:

	bool hitsTouching(const tss::Hit2D & h1, const tss::Hit2D & h2) const;
	bool hitsTouching(const tss::Cluster2D & c1, const tss::Hit2D & h2) const;
	bool hitsTouching(const tss::Cluster2D & c1, const tss::Cluster2D & c2) const;

	void merge(std::vector< tss::Cluster2D > & clusters) const;
};

#endif

