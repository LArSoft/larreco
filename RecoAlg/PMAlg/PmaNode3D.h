/**
 *  @file   PmaNode3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track node. See PmaTrack3D.h file for details.
 */

#ifndef PmaNode3D_h
#define PmaNode3D_h

#include "RecoAlg/PMAlg/PmaElement3D.h"
#include "RecoAlg/PMAlg/SortedObjects.h"

#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"

#include "TVectorT.h"
#include "TMatrixT.h"

namespace pma
{
	class Node3D;
}

class pma::Node3D : public pma::Element3D, public pma::SortedBranchBase
{
public:
	Node3D(void);
	Node3D(const TVector3& p3d, unsigned int tpc, unsigned int cryo);
	virtual ~Node3D(void) {}

	unsigned int TPC(void) const { return fTPC; }
	unsigned int Cryo(void) const { return fCryo; }

	TVector3 const & Point3D(void) const { return fPoint3D; }

	void SetPoint3D(const TVector3& p3d);

	TVector2 const & Projection2D(unsigned int view) const { return fProj2D[view]; }

	/// Distance [cm] from the 3D point to the point 3D.
	virtual double GetDistance2To(const TVector3& p3d) const;

	/// Distance [cm] from the 2D point to the object's 2D projection in one of wire views.
	virtual double GetDistance2To(const TVector2& p2d, unsigned int view) const;

	/// Set hit 3D position and its 2D projection to the vertex.
	virtual void SetProjection(pma::Hit3D& h) const;

	/// Squared sum of half-lengths of connected 3D segments
	/// (used in the vertex position optimization).
	virtual double Length2(void) const;

	/// Cosine of 3D angle between connected segments.
	double SegmentCos(void) const;
	/// Cosine of 2D angle (in plane parallel to wire planes) between connected segments.
	/// Should be changed / generalized for horizontal wire planes (e.g. 2-phase LAr).
	double SegmentCosWirePlane(void) const;
	/// Cosine of 2D angle (in horizontal plane, parallel to drift) between connected segments.
	/// Should be changed / generalized for horizontal wire planes (e.g. 2-phase LAr).
	double SegmentCosTransverse(void) const;

	/// Objective function minimized during oprimization.
	double GetObjFunction(float penaltyValue, float endSegWeight) const;

	/// Optimize vertex 3D position with given penalty on connected
	/// segments angle and weight assigned to the outermost segments.
	/// Only MSE is used in case of branching nodes.
	void Optimize(float penaltyValue, float endSegWeight);

	virtual void ClearAssigned(pma::Track3D* trk = 0);

private:
	void LimitPoint3D(float margin = -3.0F); // default: let the node go out by 3cm
	void UpdateProj2D(void);

	double EndPtCos2Transverse(void) const;
	double PiInWirePlane(void) const;
	double PenaltyInWirePlane(void) const;

	double Pi(float endSegWeight) const;
	double Penalty(float endSegWeight) const;
	double Mse(void) const;

	double MakeGradient(float penaltyValue, float endSegWeight);
	double StepWithGradient(float alfa, float tol, float penalty, float weight);

	art::ServiceHandle<geo::Geometry> fGeom;
	art::ServiceHandle<util::DetectorProperties> fDetProp;

	unsigned int fTPC, fCryo;

	double fMinX, fMaxX, fMinY, fMaxY, fMinZ, fMaxZ; // TPC boundaries to limit the node position (+margin)
	double fWirePitch[3], fDriftPitch;               // TPC params to scale do [cm] domain

	TVector3 fPoint3D;       // node position in 3D space in [cm]
	TVector2 fProj2D[3];     // node projections to 2D views, scaled to [cm], updated on each change of 3D position

	TVector3 fGradient, fGDirX, fGDirY, fGDirZ;
	TVectorT<double> fGradV;
	TMatrixT<double> fGradM;

	static bool fGradFixed[3];
};

#endif

