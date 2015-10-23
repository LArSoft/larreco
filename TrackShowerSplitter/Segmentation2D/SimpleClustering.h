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

	Cluster2D(void) {}
	Cluster2D(const std::vector< const tss::Hit2D* > & hits);

	const std::vector< const tss::Hit2D* > & hits(void) const { return fHits; }
	std::vector< const tss::Hit2D* > & hits(void) { return fHits; }

	const tss::Hit2D* outermost(void) const;

private:

	std::vector< const tss::Hit2D* > fHits;

};

class tss::SimpleClustering
{
public:

	std::vector< tss::Cluster2D > run(const std::vector< tss::Hit2D > & inp);
	std::vector< tss::Cluster2D > run(const std::vector< const tss::Hit2D* > & inp);

private:

	bool hitsTouching(const tss::Hit2D & h1, const tss::Hit2D & h2) const;
	bool hitsTouching(const tss::Cluster2D & c1, const tss::Hit2D & h2) const;
	bool hitsTouching(const tss::Cluster2D & c1, const tss::Cluster2D & c2) const;

	void merge(std::vector< tss::Cluster2D > & clusters) const;
};

#endif

