/**
 *  @file   Segmentation2D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Split into linear clusters.
 */

#ifndef Segmentation2D_h
#define Segmentation2D_h

#include "SimpleClustering.h"

namespace tss
{
	class Segmentation2D;
}

class tss::Segmentation2D
{
public:
	Segmentation2D(void) :
		fRadiusMin(1.6), fRadiusMax(3.2),
		fMaxLineDist(0.2)
	{ }

	std::vector< tss::Cluster2D > run(tss::Cluster2D & inp) const;

	void splitHits(
		const std::vector< tss::Cluster2D > & inp,
		std::vector< const tss::Hit2D* > & trackHits,
		std::vector< const tss::Hit2D* > & emHits) const;

private:

	void run(
		tss::Cluster2D & inp,
		std::vector< tss::Cluster2D > & result,
		std::vector< TVector2 > & centers) const;

	tss::Cluster2D buildSegment(tss::Cluster2D & inp, TVector2 center, TVector2 end) const;
	tss::Cluster2D selectRing(const tss::Cluster2D & inp, TVector2 center) const;

	tss::SimpleClustering fSimpleClustering;

	double fRadiusMin, fRadiusMax;
	double fMaxLineDist;
};

#endif

