/**
 *  @file   PmaVtxCandidate.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Vertex finding helper for the Projection Matching Algorithm
 *
 *          Candidate for 3D vertex. Used to test intersections and join tracks in vertices.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaVtxCandidate.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include "Geometry/Geometry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

const double pma::VtxCandidate::kMaxDistToTrack = 4.0; // max. dist. track to center to create vtx
const double pma::VtxCandidate::kMinDistToNode = 3.0;  // min. dist. to node needed to split segment

bool pma::VtxCandidate::Has(pma::Track3D* trk) const
{
	for (size_t i = 0; i < fAssigned.size(); i++)
		if (trk == fAssigned[i].first) return true;
	return false;
}

bool pma::VtxCandidate::Has(const pma::VtxCandidate& other) const
{
	for (size_t t = 0; t < other.Size(); t++)
		if (!Has(other.fAssigned[t].first)) return false;
	return true;
}

bool pma::VtxCandidate::IsAttached(pma::Track3D* trk) const
{
	pma::Track3D const * rootTrk = trk->GetRoot();
	for (size_t i = 0; i < fAssigned.size(); i++)
	{
		pma::Track3D const * rootAssn = fAssigned[i].first->GetRoot();
		if (rootTrk->IsAttachedTo(rootAssn)) return true;
	}
	return false;
}

size_t pma::VtxCandidate::Size(double minLength) const
{
	size_t n = 0;
	for (auto const & c : fAssigned)
		if (c.first->Length() > minLength) n++;
	return n;
}

bool pma::VtxCandidate::Add(pma::Track3D* trk)
{
	if (IsAttached(trk)) return false;

	fAssigned.push_back(std::pair< pma::Track3D*, size_t >(trk, 0));

	double d, d_best;
	double mse, min_mse = kMaxDistToTrack * kMaxDistToTrack;
	if (fAssigned.size() > 2)
	{
		size_t n_best = 0;
		d_best = kMaxDistToTrack;
		for (size_t n = 0; n < trk->Nodes().size() - 1; n++)
		{
			pma::Segment3D* seg = trk->NextSegment(trk->Nodes()[n]);
			if (seg->Length() < fSegMinLength) continue;

			fAssigned.back().second = n;

			mse = Compute();
			if (mse < min_mse)
			{
				d = sqrt( seg->GetDistance2To(fCenter) );
				if (d < d_best)
				{
					min_mse = mse; n_best = n; d_best = d;
				}
			}
		}

		if (d_best < kMaxDistToTrack)
		{
			fAssigned.back().second = n_best;
			fMse = Compute();
			fMse2D = ComputeMse2D();
			return true;
		}
		else
		{
			fAssigned.pop_back();
			fMse = Compute();
			fMse2D = ComputeMse2D();
			return false;
		}
	}
	else if (fAssigned.size() == 2)
	{
		pma::Track3D* p0 = fAssigned.front().first;

		size_t n_best = 0, m_best = 0;
		d_best = kMaxDistToTrack;

		double lm, ln, l_best = 0;
		for (size_t m = 0; m < p0->Nodes().size() - 1; m++)
		{
			pma::Segment3D* seg_m = p0->NextSegment(p0->Nodes()[m]);
			lm = seg_m->Length();
			if (lm < fSegMinLength) continue;

			fAssigned.front().second = m;

			for (size_t n = 0; n < trk->Nodes().size() - 1; n++)
			{
				pma::Segment3D* seg_n = trk->NextSegment(trk->Nodes()[n]);
				ln = seg_n->Length();
				if (ln < fSegMinLength) continue;

				fAssigned.back().second = n;

				mse = Compute(); // std::cout << mse << std::endl;

				d = sqrt(ComputeMse2D());

				if (d < d_best)
				{
					double d_dist = (d_best - d) / d_best;
					if (lm + ln > 0.8 * d_dist * l_best) // take closer if not much shorter
					{
						min_mse = mse;
						n_best = n; m_best = m;
						d_best = d; l_best = lm + ln;
					}
				}
			}
		}

		if (d_best < kMaxDistToTrack)
		{
			fAssigned.front().second = m_best;
			fAssigned.back().second = n_best;
			fMse = Compute();

			fMse2D = ComputeMse2D();

			return true;
		}
		else
		{
			fAssigned.pop_back();
			fCenter.SetXYZ(0., 0., 0.);
			fMse = 0; fMse2D = 0;
			return false;
		}
	}
	else
	{
		for (size_t n = 0; n < trk->Nodes().size() - 1; n++)
		{
			pma::Segment3D* seg = trk->NextSegment(trk->Nodes()[n]);
			if (seg->Length() >= fSegMinLength)
			{
				return true;
			}
		}
		fAssigned.pop_back();
		fCenter.SetXYZ(0., 0., 0.);
		fMse = 0; fMse2D = 0;
		return false;
	}
}

double pma::VtxCandidate::ComputeMse2D(void)
{
	art::ServiceHandle< geo::Geometry > geom;

	double mse = 0.0;
	TVector2 center2d;
	for (size_t i = 0; i < fAssigned.size(); i++)
	{
		pma::Track3D* trk = fAssigned[i].first;
		pma::Segment3D* seg = trk->NextSegment(trk->Nodes()[fAssigned[i].second]);

		int tpc = trk->Nodes()[fAssigned[i].second]->TPC();
		int cryo = trk->Nodes()[fAssigned[i].second]->Cryo();

		size_t k = 0;
		double m = 0.0;
		if (geom->TPC(tpc, cryo).HasPlane(geo::kU))
		{
			center2d = GetProjectionToPlane(fCenter, geo::kU, tpc, cryo);
			m += seg->GetDistance2To(center2d, geo::kU); k++;
		}
		if (geom->TPC(tpc, cryo).HasPlane(geo::kV))
		{
			center2d = GetProjectionToPlane(fCenter, geo::kV, tpc, cryo);
			m += seg->GetDistance2To(center2d, geo::kV); k++;
		}
		if (geom->TPC(tpc, cryo).HasPlane(geo::kZ))
		{
			center2d = GetProjectionToPlane(fCenter, geo::kZ, tpc, cryo);
			m += seg->GetDistance2To(center2d, geo::kZ); k++;
		}
		mse += m / (double)k;
	}

	return mse / fAssigned.size();
}

double pma::VtxCandidate::Test(const VtxCandidate& other) const
{
	double dx = fCenter[0] - other.fCenter[0];
	double dy = fCenter[1] - other.fCenter[1];
	double dz = fCenter[2] - other.fCenter[2];
	double dw = fErr[0] * other.fErr[0] * dx * dx
		+ fErr[1] * other.fErr[1] * dy * dy
		+ fErr[2] * other.fErr[2] * dz * dz;
	return sqrt(dw);
}

double pma::VtxCandidate::MaxAngle(double minLength) const
{
	TVector3 dir_i;
	size_t max_l_idx = 0;
	double l, max_l = 0.0;
	for (size_t i = 0; i < fAssigned.size() - 1; i++)
	{
		l = fAssigned[i].first->Length();
		if (l > max_l)
		{
			max_l = l; max_l_idx = i;
			pma::Track3D* trk_i = fAssigned[i].first;
			pma::Node3D* vtx_i0 = trk_i->Nodes()[fAssigned[i].second];
			pma::Node3D* vtx_i1 = trk_i->Nodes()[fAssigned[i].second + 1];
			dir_i = vtx_i1->Point3D() - vtx_i0->Point3D();
			dir_i *= 1.0 / dir_i.Mag();
		}
	}

	double a, min = 1.0;
	for (size_t j = 0; j < fAssigned.size(); j++)
		if ((j != max_l_idx) && (fAssigned[j].first->Length() > minLength))
	{
		pma::Track3D* trk_j = fAssigned[j].first;
		pma::Node3D* vtx_j0 = trk_j->Nodes()[fAssigned[j].second];
		pma::Node3D* vtx_j1 = trk_j->Nodes()[fAssigned[j].second + 1];
		TVector3 dir_j = vtx_j1->Point3D() - vtx_j0->Point3D();
		dir_j *= 1.0 / dir_j.Mag();
		a = fabs(dir_i * dir_j);
		if (a < min) min = a;
	}

	return 180.0 * acos(min) / TMath::Pi();
}

bool pma::VtxCandidate::MergeWith(const pma::VtxCandidate& other)
{
	double d = sqrt( pma::Dist2(fCenter, other.fCenter) );
	if (d > 10.0)
	{
		mf::LogVerbatim("pma::VtxCandidate") << "too far..";
		return false;
	}

	double dw = Test(other);

	size_t ntrk = 0;
	for (size_t t = 0; t < other.Size(); t++)
	{
		if (IsAttached(other.fAssigned[t].first))
		{
			mf::LogVerbatim("pma::VtxCandidate") << "already attached..";
			return false;
		}
		if (!Has(other.fAssigned[t].first))
		{
			fAssigned.push_back(other.fAssigned[t]);
			ntrk++;
		}
	}
	if (ntrk)
	{
		double mse0 = fMse, mse1 = other.fMse;
		mf::LogVerbatim("pma::VtxCandidate")
			<< "try: " << d << " mse0:" << sqrt(mse0) << " mse1:" << sqrt(mse1);

		double mse = Compute();
		mf::LogVerbatim("pma::VtxCandidate")
			<< "out: " << Size() << " mse:" << sqrt(mse) << " dw:" << dw;

		if (mse < 1.0) // kMaxDistToTrack * kMaxDistToTrack)
		{
			fMse = mse;
			fMse2D = ComputeMse2D();
			return true;
		}
		else
		{
			mf::LogVerbatim("pma::VtxCandidate") << "high mse..";
			while (ntrk--) { fAssigned.pop_back(); }
			fMse = Compute();
			fMse2D = ComputeMse2D();
			return false;
		}
	}
	else
	{
		mf::LogVerbatim("pma::VtxCandidate") << "no tracks..";
		return false;
	}
}

double pma::VtxCandidate::Compute(void)
{
	std::vector< pma::Segment3D* > segments;
	std::vector< std::pair<TVector3, TVector3> > lines;
	std::vector< double > weights;
	for (size_t v = 0; v < fAssigned.size(); v++)
	{
		pma::Track3D* trk = fAssigned[v].first;
		int vIdx = fAssigned[v].second;

		pma::Node3D* vtx1 = trk->Nodes()[vIdx];
		pma::Segment3D* seg = trk->NextSegment(vtx1);
		double segLength = seg->Length();
		if (segLength >= fSegMinLength)
		{
			pma::Node3D* vtx2 = static_cast< pma::Node3D* >(seg->Next(0));

			std::pair<TVector3, TVector3> endpoints(vtx1->Point3D(), vtx2->Point3D());
			double dy = endpoints.first.Y() - endpoints.second.Y();
			double fy_norm = asin(fabs(dy) / segLength) / (0.5 * TMath::Pi());
			double w = 1.0 - pow(fy_norm - 1.0, 12);
			if (w < 0.3) w = 0.3;

			lines.push_back(endpoints);
			segments.push_back(seg);
			weights.push_back(w);
		}
	}

	fCenter.SetXYZ(0., 0., 0.); fErr.SetXYZ(0., 0., 0.);

	TVector3 result;
	double resultMse = pma::SolveLeastSquares3D(lines, result);
	if (resultMse < 0.0)
	{
		mf::LogWarning("pma::VtxCandidate") << "Cannot compute crossing point.";
		return 1.0E+6;
	}

	TVector3 pproj;
	//double dx, dy, dz
	double wsum = 0.0;
	for (size_t s = 0; s < segments.size(); s++)
	{
		pma::Node3D* vprev = static_cast< pma::Node3D* >(segments[s]->Prev());
		pma::Node3D* vnext = static_cast< pma::Node3D* >(segments[s]->Next(0));

		pproj = pma::GetProjectionToSegment(result, vprev->Point3D(), vnext->Point3D());

		//dx = weights[s] * (result.X() - pproj.X());
		//dy = result.Y() - pproj.Y();
		//dz = result.Z() - pproj.Z();

		fErr[0] += weights[s] * weights[s];
		fErr[1] += 1.0;
		fErr[2] += 1.0;

		fCenter[0] += weights[s] * pproj.X();
		fCenter[1] += pproj.Y();
		fCenter[2] += pproj.Z();
		wsum += weights[s];
	}
	fCenter[0] /= wsum;
	fCenter[1] /= segments.size();
	fCenter[2] /= segments.size();

	fErr *= 1.0 / segments.size();
	fErr[0] = sqrt(fErr[0]);
	fErr[1] = sqrt(fErr[1]);
	fErr[2] = sqrt(fErr[2]);

	return resultMse;
}

bool pma::VtxCandidate::JoinTracks(
	std::vector< pma::Track3D* >& tracks,
	std::vector< pma::Track3D* >& src)
{
	if (tracksJoined)
	{
		mf::LogError("pma::VtxCandidate") << "Tracks already attached to the vertex.";
		return false;
	}
	tracksJoined = true;

	pma::Node3D* vtxCenter = 0;

	mf::LogVerbatim("pma::VtxCandidate") << "JoinTracks at:"
		<< " vx:" << fCenter.X() << " vy:" << fCenter.Y() << " vz:" << fCenter.Z();

	for (size_t i = 0; i < fAssigned.size(); i++)
		for (size_t t = 0; t < src.size(); t++)
			if (fAssigned[i].first == src[t])
		{
			tracks.push_back(src[t]);
			src.erase(src.begin() + t);
			break;
		}

	for (size_t i = 0; i < fAssigned.size(); i++)
	{
		mf::LogVerbatim("pma::VtxCandidate") << "----------> track #" << i;

		pma::Track3D* trk = fAssigned[i].first;
		size_t idx = fAssigned[i].second;

		mf::LogVerbatim("pma::VtxCandidate") << "  track size:" << trk->size()
			<< " (nodes:" << trk->Nodes().size() << ")";

		TVector3 p0(trk->Nodes()[idx]->Point3D());
		TVector3 p1(trk->Nodes()[idx + 1]->Point3D());

		int tpc0 = trk->Nodes()[idx]->TPC();
		int tpc1 = trk->Nodes()[idx + 1]->TPC();

		int cryo0 = trk->Nodes()[idx]->Cryo();
		int cryo1 = trk->Nodes()[idx + 1]->Cryo();

		double d0 = sqrt( pma::Dist2(p0, fCenter) );
		double d1 = sqrt( pma::Dist2(p1, fCenter) );
		double ds = sqrt( pma::Dist2(p0, p1) );
		double f = pma::GetSegmentProjVector(fCenter, p0, p1);
		TVector3 proj = pma::GetProjectionToSegment(fCenter, p0, p1);

		//mf::LogVerbatim("pma::VtxCandidate") << "  seg len:" << ds << " d0:" << d0 << " d1:" << d1;
		//mf::LogVerbatim("pma::VtxCandidate") << "  idx:" << idx << " f:" << f << " df0:" << f*ds << " df1:" << (1.0-f)*ds;

		if ((idx == 0) && (f * ds <= kMinDistToNode))
		{
			if (i == 0)
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  new at front";
				vtxCenter = trk->Nodes().front();
				vtxCenter->SetPoint3D(fCenter);
			}
			else
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  front to center";
				trk->AttachTo(vtxCenter);
			}
		}
		else if ((idx + 2 == trk->Nodes().size()) && ((1.0 - f) * ds <= kMinDistToNode))
		{
			mf::LogVerbatim("pma::VtxCandidate") << "  flip trk";
			trk->Flip();
			if (i == 0)
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  ...to make new center";
				vtxCenter = trk->Nodes().front();
				vtxCenter->SetPoint3D(fCenter);
			}
			else
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  ...to attach back to center";
				trk->AttachTo(vtxCenter);
			}
		}
		else
		{
			mf::LogVerbatim("pma::VtxCandidate") << "  split track";

			pma::Track3D* trk0 = new pma::Track3D();
			for (size_t j = 0; j <= idx; j++)
				trk0->AddNode(trk->ExtractNodeCopy(j));

			if ((f >= 0.0F) && (f <= 1.0) &&
			    (f * ds > kMinDistToNode) && ((1.0 - f) * ds > kMinDistToNode))
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  add center inside segment";

				int tpc, cryo;
				if (f < 0.5) { tpc = tpc0; cryo = cryo0; }
				else { tpc = tpc1; cryo = cryo1; }

				trk->InsertNode(fCenter, ++idx, tpc, cryo);
				trk0->AddNode(fCenter, tpc, cryo);
			}
			else
			{
				if (d1 < d0)
				{
					mf::LogVerbatim("pma::VtxCandidate") << "  add center at end of segment";

					++idx;
					trk0->AddNode(trk->ExtractNodeCopy(idx));
				}
				else
				{
					mf::LogVerbatim("pma::VtxCandidate") << "  center at start of segment - no action";
				}
			}

			mf::LogVerbatim("pma::VtxCandidate") << "  flip trk0";
			trk0->Flip();
			if (i == 0)
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  center at trk0 front";
				vtxCenter = trk0->Nodes().front();
			}

			mf::LogVerbatim("pma::VtxCandidate") << "  remove vtxs from trk";
			for (size_t j = 0; j < idx; j++) trk->RemoveNode(0);

			mf::LogVerbatim("pma::VtxCandidate") << "  reassign hits";
			size_t j = 0;
			while (j < trk->size())
			{
				pma::Hit3D* h3d = (*trk)[j];
				double dist2D_old = trk->Dist2(h3d->Point2D(), h3d->View2D());
				double dist2D_new = trk0->Dist2(h3d->Point2D(), h3d->View2D());

				if (dist2D_new < dist2D_old) trk0->push_back(trk->release_at(j));
				else j++;
			}
			mf::LogVerbatim("pma::VtxCandidate") << "  trk size:" << trk->size() << " (nodes:" << trk->Nodes().size() << ")";
			mf::LogVerbatim("pma::VtxCandidate") << "  trk0 size:" << trk0->size() << " (nodes:" << trk0->Nodes().size() << ")";

			if (trk0->HasTwoViews())
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  attach trk0";
				trk0->AttachTo(vtxCenter);
			}
			else
			{
				while (trk0->size()) trk->push_back(trk0->release_at(0));
			}

			if (trk->HasTwoViews())
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  attach trk";
				trk->AttachTo(vtxCenter);
			}
			else
			{
				while (trk->size()) trk0->push_back(trk->release_at(0));
			}

			if (trk0->size())
			{
				trk0->MakeProjection();
				tracks.push_back(trk0);
			}
			else delete trk0;

			if (trk->size()) trk->MakeProjection();
			else
			{
				for (size_t t = 0; t < tracks.size(); ++t)
					if (tracks[t] == trk)
					{
						tracks.erase(tracks.begin() + t);
						delete trk;
						break;
					}
			}

			mf::LogVerbatim("pma::VtxCandidate") << "  done";
		}
	}

	if (vtxCenter && vtxCenter->NextCount())
	{
		fAssigned.clear();
		fCenter = vtxCenter->Point3D();
		fMse = 0.0; fMse2D = 0.0;

		pma::Segment3D* rootSeg = static_cast< pma::Segment3D* >(vtxCenter->Next(0));
		rootSeg->Parent()->GetRoot()->TuneFullTree();
		return true;
	}
	else
	{
		mf::LogError("pma::VtxCandidate") << "No common vertex";
		return false;
	}
}

