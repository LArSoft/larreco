/**
 *  @file   Segmentation2D.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Split into linear clusters.
 */

#ifndef Segmentation2D_h
#define Segmentation2D_h

#include "fhiclcpp/fwd.h"

#include "SimpleClustering.h"

namespace tss
{
	class Segmentation2D;
}

class tss::Segmentation2D
{
public:

	Segmentation2D(const fhicl::ParameterSet& p) :
		fRadiusMin(0.5), fRadiusMax(1.0),
		fMaxLineDist(0.2),
		fDenseVtxRadius(1.0), fDenseHitRadius(5.0),
		fDenseMinN(5), fDenseMinH(100)
	{ reconfigure(p); }

	void reconfigure(const fhicl::ParameterSet& p);

	std::vector< tss::Cluster2D > run(tss::Cluster2D & inp) const;

	void splitHits(
		const std::vector< tss::Cluster2D > & inp,
		std::vector< const tss::Hit2D* > & trackHits,
		std::vector< const tss::Hit2D* > & emHits) const;

	void splitHitsNaive(
		const tss::Cluster2D & inp,
		std::vector< const tss::Hit2D* > & trackHits,
		std::vector< const tss::Hit2D* > & emHits) const;

	void splitHitsNaive(
		const std::vector< tss::Cluster2D > & inp,
		std::vector< const tss::Hit2D* > & trackHits,
		std::vector< const tss::Hit2D* > & emHits) const;

	int mergeClusters(
		std::vector< tss::Cluster2D > & group,
		const std::vector< size_t > & idxs) const;

private:

	void run(
		tss::Cluster2D & inp,
		std::vector< tss::Cluster2D > & result,
		std::vector< TVector2 > & centers) const;

	tss::Cluster2D buildSegment(tss::Cluster2D & inp, TVector2 center, TVector2 end) const;
	tss::Cluster2D selectRing(const tss::Cluster2D & inp, TVector2 center) const;

	void tagDenseEnds(std::vector< tss::Cluster2D > & group) const;
	void mergeDenseParts(std::vector< tss::Cluster2D > & group) const;

	bool Cl2InsideCl1(tss::Cluster2D& cl1, tss::Cluster2D& cl2) const;

	tss::SimpleClustering fSimpleClustering;

	double fRadiusMin, fRadiusMax;
	double fMaxLineDist;

	double fDenseVtxRadius, fDenseHitRadius;
	size_t fDenseMinN, fDenseMinH;
};

#endif
