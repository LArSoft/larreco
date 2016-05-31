////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTracking
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), June 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgTracking.h"

#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

recob::Track pma::convertFrom(const pma::Track3D& src, unsigned int tidx)
{
	std::vector< TVector3 > xyz, dircos;
	xyz.reserve(src.size()); dircos.reserve(src.size());

	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = src.FrontCryo();
	unsigned int tpc = src.FrontTPC();

	art::ServiceHandle<geo::Geometry> geom;
	std::map< unsigned int, pma::dedx_map > src_dQdx;
	if (geom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		src_dQdx[geo::kU] = pma::dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kU], geo::kU);
	}
	if (geom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		src_dQdx[geo::kV] = pma::dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kV], geo::kV);
	}
	if (geom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		src_dQdx[geo::kZ] = pma::dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kZ], geo::kZ);
	}

	TVector3 p3d;
	double xshift = src.GetXShift();
	bool has_shift = (xshift != 0.0);
	for (size_t i = 0; i < src.size(); i++)
		if (src[i]->IsEnabled())
	{
		p3d = src[i]->Point3D();
		if (has_shift) p3d.SetX(p3d.X() + xshift);
		xyz.push_back(p3d);

		if (i < src.size() - 1)
		{
			size_t j = i + 1;
			double mag = 0.0;
			TVector3 dc(0., 0., 0.);

			while ((mag == 0.0) && (j < src.size()))
			{
				dc = src[j]->Point3D();
				dc -= src[i]->Point3D();
				mag = dc.Mag();
				j++;
			}

			if (mag > 0.0) dc *= 1.0 / mag;
			else if (!dircos.empty()) dc = dircos.back();

			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());

		double dQ = 0., dx = 0.;
		dst_dQdx[geo::kU].push_back(0.);
		dst_dQdx[geo::kV].push_back(0.);
		dst_dQdx[geo::kZ].push_back(0.);

		double dQdx;
		for (auto const& m : src_dQdx)
		{
			auto it = m.second.find(i);
			if (it != m.second.end())
			{
				dQ = it->second[5];
				dx = it->second[6];
				if (dx > 0.) dQdx = dQ/dx;
				else dQdx = 0.;

				size_t backIdx = dst_dQdx[m.first].size() - 1;
				dst_dQdx[m.first][backIdx] = dQdx;

				break;
			}
		}
	}

	// 0 is track-like (long and/or very straight, well matching 2D hits);
	// 0x10000 is EM shower-like trajectory *** to be replaced with a better tag ***
	unsigned int pidTag = 0;
	if (src.GetTag() == pma::Track3D::kEmLike) pidTag = 0x10000;

	if (xyz.size() != dircos.size())
	{
		mf::LogError("PMAlgTrackMaker") << "pma::Track3D to recob::Track conversion problem.";
	}

	return recob::Track(xyz, dircos, dst_dQdx, std::vector< double >(2, util::kBogusD), tidx + pidTag);
}
// ------------------------------------------------------

