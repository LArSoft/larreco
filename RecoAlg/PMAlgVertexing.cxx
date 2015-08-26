////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgVertexing
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), August 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "RecoAlg/PMAlgVertexing.h"

#include "RecoAlg/PMAlg/PmaVtxCandidate.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::PMAlgVertexing::PMAlgVertexing(const fhicl::ParameterSet& pset)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

pma::PMAlgVertexing::~PMAlgVertexing(void)
{
	cleanTracks();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::reconfigure(const fhicl::ParameterSet& pset)
{
	fMinTrackLength = pset.get< double >("MinTrackLength");

	//fInputVtxDist2D = pset.get< double >("InputVtxDist2D");
	//fInputVtxDistY = pset.get< double >("InputVtxDistY");
}
// ------------------------------------------------------

void pma::PMAlgVertexing::cleanTracks(void)
{
	for (auto t : fOutTracks) delete t;
	fOutTracks.clear();

	for (auto t : fShortTracks) delete t;
	fShortTracks.clear();

	for (auto t : fEmTracks) delete t;
	fEmTracks.clear();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::collectTracks(
	std::vector< pma::Track3D* >& result)
{
	for (auto t : result) delete t;
	result.clear();

	for (auto t : fOutTracks) result.push_back(t);
	fOutTracks.clear();

	for (auto t : fShortTracks) result.push_back(t);
	fShortTracks.clear();

	for (auto t : fEmTracks) result.push_back(t);
	fEmTracks.clear();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::sortTracks(
	const std::vector< pma::Track3D* >& trk_input)
{
	cleanTracks();

	for (auto t : trk_input)
	{
		double l = t->Length();
		if ((t->GetTag() != pma::Track3D::kEmLike) || (l > 2 * fMinTrackLength))
		{
			fOutTracks.push_back(new pma::Track3D(*t));

			//if (l > fMinTrackLength) fOutTracks.push_back(new pma::Track3D(*t));
			//else fShortTracks.push_back(new pma::Track3D(*t));
		}
		else
		{
			fEmTracks.push_back(new pma::Track3D(*t));
		}
	}
}
// ------------------------------------------------------

bool pma::PMAlgVertexing::findOneVtx(void)
{
	if (fOutTracks.size() < 2) return false;

	std::vector< pma::VtxCandidate > candidates;
	for (size_t t = 0; t < fOutTracks.size() - 1; t++)
	{
		if (fOutTracks[t]->Length() < 2 * fMinTrackLength) continue;

		size_t u = t + 1;
		while (u < fOutTracks.size())
		{
			pma::VtxCandidate candidate;
			if (!candidate.Add(fOutTracks[t])) break; // no segments with length > thr

			// **************************** try Mse2D / or only Mse ************************************
			if (candidate.Add(fOutTracks[u++]) && (sqrt(candidate.Mse()) < 1.0))
			//if (candidate.Add(fOutTracks[u++]) && (sqrt(candidate.Mse()) < 2.0) && (candidate.Mse2D() < 1.0))
				candidates.push_back(candidate);
		}
	}

	bool merged = true;
	while (merged && (candidates.size() > 1))
	{
		merged = false;
		double d_thr = 1.0; // 1.0 = max weighted dist. threshold
		double d, dmin = d_thr;

		size_t k_best, l_best, k = 0;
		while (k < candidates.size() - 1)
		{
			size_t l = k + 1;
			while (l < candidates.size())
			{
				if (candidates[l].Has(candidates[k]))
				{
					candidates[k] = candidates[l];
					candidates.erase(candidates.begin() + l);
					mf::LogVerbatim("pma::PMAlgVertexing") << "removed (k)";
				}
				else if (candidates[k].Has(candidates[l]))
				{
					candidates.erase(candidates.begin() + l);
					mf::LogVerbatim("pma::PMAlgVertexing") << "removed (l)";
				}
				else
				{
					d = candidates[k].Test(candidates[l]);
					if (d < dmin) { dmin = d; k_best = k; l_best = l; }
					l++;
				}
			}
			k++;
		}
		if ((dmin < d_thr) && candidates[k_best].MergeWith(candidates[l_best]))
		{
			mf::LogVerbatim("pma::PMAlgVertexing") << "merged";
			candidates.erase(candidates.begin() + l_best);
			merged = true;
		}
	}

	int nmax = 0, c_best = -1;
	double amax = 0.0;

	mf::LogVerbatim("pma::PMAlgVertexing") << "*** Vtx candidates: " << candidates.size();
	for (size_t v = 0; v < candidates.size(); v++)
	{
		//vtxSel->push_back(new Vertex(candidates[v].Center()));

		if (((int)candidates[v].Size() > nmax) ||
			(((int)candidates[v].Size() == nmax) && (candidates[v].MaxAngle() > amax)))
		{
			nmax = candidates[v].Size();
			amax = candidates[v].MaxAngle();
			c_best = v;
		}

		mf::LogVerbatim("pma::PMAlgVertexing")
			<< "center x:" << candidates[v].Center().X()
			<< " y:" << candidates[v].Center().Y()
			<< " z:" << candidates[v].Center().Z();

		for (size_t i = 0; i < candidates[v].Size(); i++)
			mf::LogVerbatim("pma::PMAlgVertexing")
				<< "     trk:" << i << " "
				<< candidates[v].Track(i).first->size()
				<< "(" << candidates[v].Weight(i) << ")";

		mf::LogVerbatim("pma::PMAlgVertexing")
			<< " dist 3D:" << sqrt(candidates[v].Mse())
			<< " 2D:" << sqrt(candidates[v].Mse2D())
			<< " " << candidates[v].MaxAngle();
	}

	if (c_best >= 0)
	{
		candidates[c_best].JoinTracks(fOutTracks);
		fOutVertices.push_back(candidates[c_best].Center());
		return true;
	}
	else return false;
}
// ------------------------------------------------------

size_t pma::PMAlgVertexing::run(
	std::vector< pma::Track3D* >& trk_input)
{
	if (trk_input.size() < 2)
	{
		mf::LogWarning("pma::PMAlgVertexing") << "no source tracks!";
		return 0;
	}

	sortTracks(trk_input); // copy input and split by tag/size

	size_t nvtx = 0;
	while (findOneVtx()) nvtx++;
	mf::LogVerbatim("pma::PMAlgVertexing") << "  " << nvtx << " vertices.";

	collectTracks(trk_input);

	return nvtx;
}
// ------------------------------------------------------

size_t pma::PMAlgVertexing::run(
	std::vector< pma::Track3D* >& trk_input,
	const std::vector< TVector3 >& vtx_input)
{
	sortTracks(trk_input); // copy input and split by tag/size

	// ....

	//collectTracks(trk_input); // return output in place of (deleted) input

	return 0;
}
// ------------------------------------------------------

