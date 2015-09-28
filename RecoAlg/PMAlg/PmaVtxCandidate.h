/**
 *  @file   PmaVtxCandidate.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Vertex finding helper for the Projection Matching Algorithm
 *
 *          Candidate for 3D vertex. Used to test intersections and join tracks in vertices.
 *          See PmaTrack3D.h file for details.
 */

#ifndef VtxCandidate_h
#define VtxCandidate_h

#include "RecoAlg/PMAlg/PmaTrack3D.h"

namespace pma
{
	class VtxCandidate;
}

class pma::VtxCandidate
{
public:
	static const double kMaxDistToTrack;
	static const double kMinDistToNode;

	VtxCandidate(double segMinLength = 0.5) :
		tracksJoined(false),
		fSegMinLength(segMinLength),
		fMse(0.0), fMse2D(0.0),
		fCenter(0., 0., 0.),
		fErr(0., 0., 0.)
	{}

	bool Has(pma::Track3D* trk) const;

	bool Has(const VtxCandidate& other) const;

	bool IsAttached(pma::Track3D* trk) const;

	bool Add(pma::Track3D* trk);

	double ComputeMse2D(void);

	double Test(const VtxCandidate& other) const;

	double MaxAngle(double minLength = 0.0) const;

	size_t Size(double minLength = 0.0) const;

	bool MergeWith(const VtxCandidate& other);

	double Compute(void);

	bool JoinTracks(
		std::vector< pma::Track3D* >& tracks,
		std::vector< pma::Track3D* >& src);

	const TVector3& Center(void) const { return fCenter; }
	double Mse(void) const { return fMse; }
	double Mse2D(void) const { return fMse2D; }

	std::pair< pma::Track3D*, size_t > Track(size_t i) const { return fAssigned[i]; }

private:
	bool tracksJoined;
	double fSegMinLength, fMse, fMse2D;
	std::vector< std::pair< pma::Track3D*, size_t > > fAssigned;
	TVector3 fCenter, fErr;
};

#endif

