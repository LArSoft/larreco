/**
 *  @file   SimpleClustering.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Trivial, collect hits "touching" each other (next wire or consecutive ticks),
 *          plus Cluster2D class to hold data needed by tss algorithms.
 */

#ifndef SimpleClustering_h
#define SimpleCustering_h

#include "TssHit2D.h"
#include "RecoAlg/PMAlg/Utilities.h"

namespace tss
{
	struct bDistToPointLess;
	class Cluster2D;

	class SimpleClustering;
}

struct tss::bDistToPointLess :
	public std::binary_function<const tss::Hit2D*, const tss::Hit2D*, bool>
{
	bDistToPointLess(const TVector2& point) : p0(point) { }
	
	bool operator() (const tss::Hit2D* h1, const tss::Hit2D* h2)
	{
		if (h1 && h2) return pma::Dist2(h1->Point2D(), p0) < pma::Dist2(h2->Point2D(), p0);
		else return false;
	}
	
	private: TVector2 p0;
};

class tss::Cluster2D
{
public:

	Cluster2D(void) { }
	Cluster2D(const std::vector< const tss::Hit2D* > & hits);

	size_t size(void) const { return fHits.size(); }

	const tss::Hit2D & operator [] (size_t index) const { return *(fHits[index]); }

	const std::vector< const tss::Hit2D* > & hits(void) const { return fHits; }
	std::vector< const tss::Hit2D* > & hits(void) { return fHits; }

	bool has(const tss::Hit2D* hit) const;

	double dist2(const TVector2 & p2d) const;
	double dist2(const TVector2 & p2d, size_t & hIdx) const;
	double dist2(const tss::Cluster2D & clu) const;

	const tss::Hit2D* release_at(size_t idx);
	bool release(const tss::Hit2D* hit);

	void push_back(const tss::Hit2D* hit) { fHits.push_back(hit); }
	void take_from(tss::Cluster2D & clu, size_t idx)
	{
		const tss::Hit2D* hit = clu.release_at(idx);
		if (hit) push_back(hit);
	}

	const tss::Hit2D* start(void) const
	{
		if (fHits.size()) return fHits.front();
		else return 0;
	}
	const tss::Hit2D* end(void) const
	{
		if (fHits.size()) return fHits.back();
		else return 0;
	}
	void sort(void)
	{
		if (fHits.size() > 2)
			std::sort(fHits.begin() + 1, fHits.end(),
				tss::bDistToPointLess(fHits.front()->Point2D()));
	}

	const tss::Hit2D* closest(const TVector2 & p2d, size_t & idx) const;
	const tss::Hit2D* outermost(size_t & idx) const;

private:

	std::vector< const tss::Hit2D* > fHits;

};

class tss::SimpleClustering
{
public:

	std::vector< tss::Cluster2D > run(const std::vector< tss::Hit2D > & inp) const;
	std::vector< tss::Cluster2D > run(const tss::Cluster2D & inp) const;

	bool hitsTouching(const tss::Hit2D & h1, const tss::Hit2D & h2) const;
	bool hitsTouching(const tss::Cluster2D & c1, const tss::Hit2D & h2) const;
	bool hitsTouching(const tss::Cluster2D & c1, const tss::Cluster2D & c2) const;

private:

	void merge(std::vector< tss::Cluster2D > & clusters) const;
};

#endif

