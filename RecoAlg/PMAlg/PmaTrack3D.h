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
 *          Progress:
 *             May-June 2015:  basic functionality for single (not-branching) track 3D optimization and dQ/dx.
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
	enum EDirection { kForward = -1, kBackward = 1 };

	Track3D(void);
	Track3D(const Track3D& src);
	~Track3D(void);

	void Initialize(float initEndSegW = 0.05F);

	bool push_back(art::Ptr< recob::Hit > hit);
	bool erase(art::Ptr< recob::Hit > hit);

	pma::Hit3D* & operator [] (size_t index) { return fHits[index]; }
	pma::Hit3D* const & operator [] (size_t index) const { return fHits[index]; }
	pma::Hit3D* & front() { return fHits.front(); }
	pma::Hit3D* const & front() const { return fHits.front(); }
	pma::Hit3D* & back() { return fHits.back(); }
	pma::Hit3D* const & back() const { return fHits.back(); }
	size_t size() const { return fHits.size(); }

	double Length(size_t step = 1) const { return Length(0, size() - 1, step); }
	double Length(size_t start, size_t stop, size_t step = 1) const;

	double Dist2(const TVector2& p2d, unsigned int view) const;
	double Dist2(const TVector3& p3d) const;

	/// Add hits; does not update hit->node/seg assignments nor hit projection to track,
	/// so MakeProjection() and SortHits() should be called as needed.
	void AddHits(const std::vector< art::Ptr<recob::Hit> >& hits);

	/// Remove hits; removes also hit->node/seg assignments.
	void RemoveHits(const std::vector< art::Ptr<recob::Hit> >& hits);

	unsigned int NHits(unsigned int view) const;
	unsigned int NEnabledHits(unsigned int view = geo::kUnknown) const;

	std::vector< unsigned int > TPCs(void) const;
	std::vector< unsigned int > Cryos(void) const;

	unsigned int FrontTPC(void) const { return fNodes.front()->TPC(); }
	unsigned int FrontCryo(void) const { return fNodes.front()->Cryo(); }

	unsigned int BackTPC(void) const { return fNodes.back()->TPC(); }
	unsigned int BackCryo(void) const { return fNodes.back()->Cryo(); }

	/// Invert the order of hits and vertices in the track.
	void Flip(void);

	/// MSE of 2D hits.
	double TestHitsMse(const std::vector< art::Ptr<recob::Hit> >& hits,
		bool normalized = true) const; // normalize to the number of hits

	/// Count close 2D hits.
	unsigned int TestHits(const std::vector< art::Ptr<recob::Hit> >& hits,
		double dist = 0.4) const; // max acceptable distance [cm]

	int NextHit(int index, unsigned int view = geo::kZ, bool inclDisabled = false) const;
	int PrevHit(int index, unsigned int view = geo::kZ, bool inclDisabled = false) const;

	/// Length of the track part associated with index'th hit. Calculated as a half distance to
	/// the preceding hit plus half distance to the subsequent hit. In case of the first (last)
	/// hit - missing part is estimated as 1/4 of the distance to the next (previous) hit.
	/// NOTE: only hits from a given view are considered; other hits are accounted for
	/// segment lengths but overall dx is calculated between hits in given view.
	double HitDxByView(size_t index, unsigned int view) const;

	/// Sequence of <hit_index, (wire, drift, X, Y, Z, dE, dx, range)> values for the track,
	/// hits tagged as outliers are skipped by default.
	/** Results are pushed into the dedx vector given in the function arguments:

	    hit (segment middle if many hits) 2D projection in view:
	      dedx[n][0] = wire;
	      dedx[n][1] = drift;

	    hit (segment middle if many hits) 3D position [cm]:
	      dedx[n][2] = X;
	      dedx[n][3] = Y;
	      dedx[n][4] = Z;

	      dedx[n][5] = dE [now ADC], energy assigned to the segment;

	      dedx[n][6] = dx [cm], length of the segment.

	      dedx[n][7] = range, total length to the track endpoint;

	    Parameters:
	      dedx  - vector to store results (empty at the begining);
	      view  - view (U, V or Z) from which dedx is created;
	      skip  - number of hits to skip at the begining (first hit has poorly estimated segment
	              length so it can be convenient to set skip=1 and handle first hit charge manually);
	      inclDisabled - if true then artificial hits added with CompleteMissingWires() are used,
	                     otherwise only true hits found in ADC are used.

	    Return value: sum of ADC's of hits skipped at the begining. */
	double GetRawdEdxSequence(std::map< size_t, std::vector<double> >& dedx, unsigned int view = geo::kZ,
		unsigned int skip = 0, bool inclDisabled = false) const;

	void AddRefPoint(const TVector3& p) { fAssignedPoints.push_back(new TVector3(p)); }
	bool HasRefPoint(TVector3* p) const;

	/// MSE of hits weighted with hit amplidudes and wire plane coefficients.
	double GetMse(void) const;

	/// Objective function optimized in track reconstruction.
	double GetObjFunction(float penaltyFactor = 1.0F) const;

	/// Main optimization method.
	double Optimize(int nNodes = -1, double eps = 0.01, bool selAllHits = true);

	/// Cut out tails with no hits assigned.
	void CleanupTails(void);

	/// Move the first/last Node3D to the first/last hit in the track;
	/// returns true if all OK, false if empty segments found.
	bool ShiftEndsToHits(void);

	pma::Segment3D* NextSegment(pma::Node3D* vtx) const;
	pma::Segment3D* PrevSegment(pma::Node3D* vtx) const;

	std::vector< pma::Node3D* > const & Nodes(void) const { return fNodes; }
	pma::Node3D* FirstElement(void) const { return fNodes.front(); }
	pma::Node3D* LastElement(void) const { return fNodes.back(); }

	void AddNode(TVector3 const & p3d, unsigned int tpc, unsigned int cryo);
	bool AddNode(void);

	void MakeProjection(void);
	void UpdateProjection(void);
	void SortHits(void);

	unsigned int DisableSingleViewEnds(void);
	void SelectHits(float fraction = 1.0F);

	float GetEndSegWeight(void) { return fEndSegWeight; }
	void SetEndSegWeight(float value) { fEndSegWeight = value; }

	float GetPenalty(void) { return fPenaltyFactor; }
	void SetPenalty(float value) { fPenaltyFactor = value; }

	unsigned int GetMaxHitsPerSeg(void) { return fMaxHitsPerSeg; }
	void SetMaxHitsPerSeg(unsigned int value) { fMaxHitsPerSeg = value; }

private:
	void ClearNodes(void);

	void InternalFlip(std::vector< pma::Track3D* >& toSort);

	void UpdateHitsRadius(void);
	double AverageDist2(void) const;

	bool PCEndpoints(TVector2 & start, TVector2 & stop,
		unsigned int view, double wpitch, double dpitch) const;
	bool InitFromHits(int tpc, int cryo, float initEndSegW = 0.05F);
	bool InitFromRefPoints(int tpc, int cryo);
	void InitFromMiddle(int tpc, int cryo);

	/// Distance to the nearest subsequent (dir = Track3D::kForward) or preceeding (dir = Track3D::kBackward)
	/// hit in given view. In case of last (first) hit in this view the half-distance in opposite direction is
	/// returned. Parameter secondDir is only for internal protection - please leave the default value.
	double HitDxByView(size_t index, unsigned int view, Track3D::EDirection dir, bool secondDir = false) const;

	void RebuildSegments(void);
	bool SwapVertices(size_t v0, size_t v1);
	void UpdateParams(void);

	bool CheckEndSegment(pma::Track3D::ETrackEnd endCode);

	int index_of(const pma::Hit3D* hit) const;
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

