/**
 *  @file   PmaTrack3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          Based on the "Precise 3D track reco..." AHEP (2013) 260820, with all the tricks that we
 *          developed later and with the work for the full-event topology optimization that is still
 *          under construction.
 *
 *          Current status: basic functionality for single segment optimization.
 */

#ifndef PmaTrack3D_h
#define PmaTrack3D_h

#include "RecoAlg/PMAlg/PmaHit3D.h"
#include "RecoAlg/PMAlg/PmaNode3D.h"
#include "RecoAlg/PMAlg/PmaSegment3D.h"

namespace pma
{
	class Track3D;
}

class pma::Track3D
{
public:
	enum ETrackEnd { kBegin = -1, kEnd = 1 };

	Track3D(void);
	~Track3D(void);

	bool push_back(const recob::Hit & hit);

	pma::Hit3D* & operator [] (size_t index) { return fHits[index]; }
	pma::Hit3D* const & operator [] (size_t index) const { return fHits[index]; }

	size_t size() const { return fHits.size(); }

	void AddHits(const std::vector< recob::Hit const* >& hits);
	unsigned int NHits(unsigned int view) const;
	unsigned int NEnabledHits(unsigned int view = geo::kUnknown) const;

	bool HasRefPoint(TVector3* p) const;

	double GetObjFunction(bool suppressPenalty = false) const;

	pma::Segment3D* NextSegment(pma::Node3D* vtx) const;
	pma::Segment3D* PrevSegment(pma::Node3D* vtx) const;

	std::vector< pma::Node3D* > const & Nodes(void) const { return fNodes; }
	pma::Node3D* FirstElement(void) const { return fNodes.front(); }
	pma::Node3D* LastElement(void) const { return fNodes.back(); }

	void MakeProjection(void);
	void UpdateProjection(void);
	void SortHits(void);

	float GetEndSegWeight(void) { return fEndSegWeight; }
	void SetEndSegWeight(float value) { fEndSegWeight = value; }

	float GetPenalty(void) { return fPenaltyFactor; }
	void SetPenalty(float value) { fPenaltyFactor = value; }

	unsigned int GetMaxHitsPerSeg(void) { return fMaxHitsPerSeg; }
	void SetMaxHitsPerSeg(unsigned int value) { fMaxHitsPerSeg = value; }

private:
	void UpdateHitsRadius(void);
	double AverageDist2(void) const;

	void RebuildSegments(void);
	bool SwapVertices(size_t v0, size_t v1);
	void UpdateParams(void);

	bool CheckEndSegment(pma::Track3D::ETrackEnd endCode);

	std::vector< pma::Hit3D* > fHits;
	std::vector< TVector3* > fAssignedPoints;

	pma::Element3D* GetNearestElement(const TVector2& p2d, unsigned int view) const;
	pma::Element3D* GetNearestElement(const TVector3& p3d) const;
	std::vector< pma::Node3D* > fNodes;
	std::vector< pma::Segment3D* > fSegments;

	unsigned int fMaxHitsPerSeg;
	float fPenaltyFactor;
	float fMaxSegStopFactor;

	unsigned int fSegStopValue;
	unsigned int fMinSegStop;
	unsigned int fMaxSegStop;

	float fSegStopFactor;
	float fPenaltyValue;
	float fEndSegWeight;
	float fHitsRadius;
};

#endif

