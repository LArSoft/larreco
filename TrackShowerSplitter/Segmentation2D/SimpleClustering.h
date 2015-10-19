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
	class SimpleClustering;
}

class tss::SimpleClustering
{

public:

	std::vector< std::vector< const tss::Hit2D* > > run(const std::vector< tss::Hit2D > & inp);
	std::vector< std::vector< const tss::Hit2D* > > run(const std::vector< const tss::Hit2D* > & inp);

private:

	bool hitsTouching(const tss::Hit2D & h1, const tss::Hit2D & h2) const;
	bool hitsTouching(const std::vector< const tss::Hit2D* > & c1, const tss::Hit2D & h2) const;
	bool hitsTouching(const std::vector< const tss::Hit2D* > & c1, const std::vector< const tss::Hit2D* > & c2) const;

	void merge(std::vector< std::vector< const tss::Hit2D* > > & clusters) const;
};

#endif

