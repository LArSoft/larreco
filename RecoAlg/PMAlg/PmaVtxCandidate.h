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
	static const double kMaxDist;

	VtxCandidate(double segMinLength = 2.0) :
		fSegMinLength(segMinLength),
		fMse(0.0), fMse2D(0.0),
		fCenter(0., 0., 0.),
		fErr(0., 0., 0.),
		fTPC(-1), fCryo(-1)
	{}

	bool Has(pma::Track3D* trk) const;

	bool Has(const VtxCandidate& other) const;

	bool IsAttached(pma::Track3D* trk) const;

	bool Add(pma::Track3D* trk);

	double ComputeMse2D(void);

	double Test(const VtxCandidate& other) const;

	double MaxAngle(void) const;

	bool MergeWith(const VtxCandidate& other);

	double Compute(void);

	void JoinTracks(std::vector< pma::Track3D* >& tracks);

	const TVector3& Center(void) const { return fCenter; }
	double Weight(size_t i) const { return fWeights[i]; }
	size_t Size(void) const { return fAssigned.size(); }
	double Mse(void) const { return fMse; }
	double Mse2D(void) const { return fMse2D; }

	std::pair< pma::Track3D*, size_t > Track(size_t i) const { return fAssigned[i]; }

private:
	double fSegMinLength, fMse, fMse2D;
	std::vector< double > fWeights;
	std::vector< std::pair< pma::Track3D*, size_t > > fAssigned;
	TVector3 fCenter, fErr;

	int fTPC, fCryo;
};

#endif

