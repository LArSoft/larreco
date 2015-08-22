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

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

const double pma::VtxCandidate::kMaxDist = 4.0;  // maximum distance to center

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
	for (size_t i = 0; i < fAssigned.size(); i++)
		if ((trk == fAssigned[i].first) ||
		    trk->IsAttachedTo(fAssigned[i].first)) return true;
	return false;
}

bool pma::VtxCandidate::Add(pma::Track3D* trk)
{
	if (IsAttached(trk)) return false;

	fAssigned.push_back(std::pair< pma::Track3D*, size_t >(trk, 0));
	fWeights.push_back(1.0);

	double d, d_best;
	double mse, min_mse = kMaxDist * kMaxDist;
	if (fAssigned.size() > 2)
	{
		size_t n_best;
		d_best = kMaxDist;
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

		if (d_best < kMaxDist)
		{
			fAssigned.back().second = n_best;
			fMse = Compute();
			fMse2D = ComputeMse2D();
			return true;
		}
		else
		{
			fAssigned.pop_back();
			fWeights.pop_back();
			fMse = Compute();
			fMse2D = ComputeMse2D();
			return false;
		}
	}
	else if (fAssigned.size() == 2)
	{
		pma::Track3D* p0 = fAssigned.front().first;

		size_t n_best, m_best;
		d_best = kMaxDist;

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

		if (d_best < kMaxDist)
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
			fWeights.pop_back();
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
				fTPC = trk->FrontTPC();
				fCryo = trk->FrontCryo();
				return true;
			}
		}
		fAssigned.pop_back();
		fWeights.pop_back();
		fCenter.SetXYZ(0., 0., 0.);
		fMse = 0; fMse2D = 0;
		fTPC = -1; fCryo = -1;
		return false;
	}
}

double pma::VtxCandidate::ComputeMse2D(void)
{
	double mse = 0.0;
	TVector2 center2d_U, center2d_V, center2d_Z;
	for (size_t i = 0; i < fAssigned.size(); i++)
	{
		pma::Track3D* trk = fAssigned[i].first;
		pma::Segment3D* seg = trk->NextSegment(trk->Nodes()[fAssigned[i].second]);

		center2d_U = GetProjectionToPlane(fCenter, geo::kU, fTPC, fCryo);
		center2d_V = GetProjectionToPlane(fCenter, geo::kV, fTPC, fCryo);
		center2d_Z = GetProjectionToPlane(fCenter, geo::kZ, fTPC, fCryo);

		mse += seg->GetDistance2To(center2d_U, geo::kU);
		mse += seg->GetDistance2To(center2d_V, geo::kV);
		mse += seg->GetDistance2To(center2d_Z, geo::kZ);
	}
	return (mse / fAssigned.size()) / 3.0;
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

double pma::VtxCandidate::MaxAngle(void) const
{
	double a, min = 1.0;
	for (size_t i = 0; i < fAssigned.size() - 1; i++)
	{
		pma::Track3D* trk_i = fAssigned[i].first;
		pma::Node3D* vtx_i0 = trk_i->Nodes()[fAssigned[i].second];
		pma::Node3D* vtx_i1 = trk_i->Nodes()[fAssigned[i].second + 1];
		TVector3 dir_i = vtx_i1->Point3D() - vtx_i0->Point3D();
		dir_i *= 1.0 / dir_i.Mag();
		for (size_t j = i + 1; j < fAssigned.size(); j++)
		{
			pma::Track3D* trk_j = fAssigned[j].first;
			pma::Node3D* vtx_j0 = trk_j->Nodes()[fAssigned[j].second];
			pma::Node3D* vtx_j1 = trk_j->Nodes()[fAssigned[j].second + 1];
			TVector3 dir_j = vtx_j1->Point3D() - vtx_j0->Point3D();
			dir_j *= 1.0 / dir_j.Mag();
			a = fabs(dir_i * dir_j);
			if (a < min) min = a;
		}
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
		if (!Has(other.fAssigned[t].first))
		{
			fAssigned.push_back(other.fAssigned[t]);
			fWeights.push_back(other.fWeights[t]);
			ntrk++;
		}
	if (ntrk)
	{
		double mse0 = fMse, mse1 = other.fMse;
		mf::LogVerbatim("pma::VtxCandidate")
			<< "try: " << d << " mse0:" << sqrt(mse0) << " mse1:" << sqrt(mse1);

		double mse = Compute();
		mf::LogVerbatim("pma::VtxCandidate")
			<< "out: " << Size() << " mse:" << sqrt(mse) << " dw:" << dw;

		if (mse < 1.0) // kMaxDist * kMaxDist)
		{
			fMse = mse;
			fMse2D = ComputeMse2D();
			return true;
		}
		else
		{
			mf::LogVerbatim("pma::VtxCandidate") << "high mse..";
			while (ntrk--) { fAssigned.pop_back(); fWeights.pop_back(); }
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
	for (size_t v = 0; v < fAssigned.size(); v++)
	{
		pma::Track3D* trk = fAssigned[v].first;
		int vIdx = fAssigned[v].second;

		pma::Node3D* vtx = trk->Nodes()[vIdx];

		pma::Segment3D* seg = trk->NextSegment(vtx);
		segments.push_back(seg);

		pma::Node3D* vtx2 = static_cast< pma::Node3D* >(seg->Next(0));
		double dy = vtx->Point3D().Y() - vtx2->Point3D().Y();

		double fy_norm = asin(fabs(dy) / seg->Length()) / (0.5 * TMath::Pi());
		fWeights[v] = 1.0 - pow(fy_norm - 1.0, 12);
		if (fWeights[v] < 0.3) fWeights[v] = 0.3;
	}
	TVector3 result;
	double resultMse = pma::SolveLeastSquares3D(lines, result);

	fCenter.SetXYZ(0., 0., 0.); fErr.SetXYZ(0., 0., 0.);

	TVector3 pproj;
	double dx, dy, dz, wsum = 0.0;
	for (size_t s = 0; s < segments.size(); s++)
	{
		pma::Node3D* vprev = static_cast< pma::Node3D* >(segments[s]->Prev());
		pma::Node3D* vnext = static_cast< pma::Node3D* >(segments[s]->Next(0));

		pproj = pma::GetProjectionToSegment(result, vprev->Point3D(), vnext->Point3D()));

		dx = fWeights[s] * (result.X() - pproj.X());
		dy = result.Y() - pproj.Y();
		dz = result.Z() - pproj.Z();

		fErr[0] += fWeights[s] * fWeights[s];
		fErr[1] += 1.0;
		fErr[2] += 1.0;

		fCenter[0] += fWeights[s] * pproj.X();
		fCenter[1] += pproj.Y();
		fCenter[2] += pproj.Z();
		wsum += fWeights[s];
	}
	fCenter[0] /= wsum;
	fCenter[1] /= segments.size();
	fCenter[2] /= segments.size();

	fErr.Scale(1.0 / segments.size());
	fErr[0] = sqrt(fErr[0]);
	fErr[1] = sqrt(fErr[1]);
	fErr[2] = sqrt(fErr[2]);

	//std::cout << fAssigned[0].first->size() << " "
	//	<< fAssigned[1].first->size() << " "
	//	<< sqrt(resultMse) << std::endl;
	//std::cout
	//	<< " vx:" << fCenter.X()
	//	<< " vy:" << fCenter.Y()
	//	<< " vx:" << fCenter.Z() << std::endl;
	return resultMse;
}

void pma::VtxCandidate::JoinTracks(std::vector< pma::Track3D* >& tracks)
{
	pma::Node3D* vtxCenter = 0;

	mf::LogVerbatim("pma::VtxCandidate") << "*** JoinTracks at:"
		<< " vx:" << fCenter.X()
		<< " vy:" << fCenter.Y()
		<< " vx:" << fCenter.Z();
	for (size_t i = 0; i < fAssigned.size(); i++)
	{
		mf::LogVerbatim("pma::VtxCandidate") << "*** JoinTracks: track #" << i;

		pma::Track3D* trk = fAssigned[i].first;
		size_t idx = fAssigned[i].second;

		mf::LogVerbatim("pma::VtxCandidate") << "  track size:" << trk->size();

		if ((int)trk->FrontTPC() != fTPC) mf::LogError("pma::VtxCandidate") << "*** WRONG TPC ***";

		TVector3 p0(trk->Nodes()[idx]->Point3D());
		TVector3 p1(trk->Nodes()[idx + 1]->Point3D());

		double d0 = sqrt( pma::Dist2(p0, fCenter) );
		double d1 = sqrt( pma::Dist2(p1, fCenter) );
		double ds = sqrt( pma::Dist2(p0, p1) );

		double f = pma::GetSegmentProjVector(fCenter, p0, p1);
		TVector3 proj = pma::GetProjectionToSegment(fCenter, p0, p1);
		double d_proj = sqrt( pma::Dist2(fCenter, proj) );

		//p0.ls(); p1.ls(); proj.ls();
		mf::LogVerbatim("pma::VtxCandidate") << "  idx:" << idx << " f:" << f;
		if ((idx == 0) && (f * ds <= 1.0))
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
		else if ((idx + 2 == trk->Nodes().size()) && ((1.0 - f) * ds <= 1.0))
		{
			mf::LogVerbatim("pma::VtxCandidate") << "  flip";
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
				trk0->AddNode(
					trk->Nodes()[j]->Point3D(),
					trk->Nodes()[j]->TPC(),
					trk->Nodes()[j]->Cryo());

			if ((f >= 0.0F) && (f <= 1.0) && (f * ds > 1.0) && ((1.0 - f) * ds > 1.0))
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  add center inside segment";
				// --------------------- !!!
				//trk->InsertNode(fCenter, ++idx); // <----------------------------------------- ???
				// --------------------- !!!
				trk0->AddNode(fCenter, fTPC, fCryo);
			}
			else
			{
				if (d1 < d0)
				{
					mf::LogVerbatim("pma::VtxCandidate") << "  add center at end of segment";

					++idx;
					trk0->AddNode(
						trk->Nodes()[idx]->Point3D(),
						trk->Nodes()[idx]->TPC(),
						trk->Nodes()[idx]->Cryo());
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
			// --------------------- !!!
			//for (size_t j = 0; j < idx; j++) trk->RemoveVertex(0); // <--------------------- ???
			// --------------------- !!!

			mf::LogVerbatim("pma::VtxCandidate") << "  reassign hits";
			size_t j = 0;
			while (j < trk->size())
			{
				pma::Hit3D* h3d = (*trk)[j];
				double dist2D_old = trk->Dist2(h3d->Point2D(), h3d->View2D()); //pl->GetDistance2To(h3d->Point());
				double dist2D_new = trk0->Dist2(h3d->Point2D(), h3d->View2D()); //pl0->GetDistance2To(h3d->Point());

				// --------------------- !!!
				if (dist2D_new < dist2D_old) trk0->push_back(NULL); //trk0->push_back(trk->release_at(j)); // <--------------------- !!!!
				else j++;
				// --------------------- !!!
			}
			mf::LogVerbatim("pma::VtxCandidate") << "  trk size:" << trk->size() << " #vtxs:" << trk->Vertices().size();
			mf::LogVerbatim("pma::VtxCandidate") << "  trk0 size:" << trk0->size() << " #vtxs:" << trk0->Vertices().size();

			if (trk0->size() > 1)
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  attach trk0";
				trk0->AttachTo(vtxCenter);
			}
			else
			{
				// --------------------- !!!
				//while (trk0->size()) trk->push_back(trk0->release_at(0));  // <--------------------- !!!!
				// --------------------- !!!
			}

			if (trk->size() > 1)
			{
				mf::LogVerbatim("pma::VtxCandidate") << "  attach trk";
				trk->AttachTo(vtxCenter);
			}
			else
			{
				// --------------------- !!!
				while (trk->size()) trk0->push_back(NULL); // trk0->push_back(trk->release_at(0)); // <--------------------- !!!!
				// --------------------- !!!
			}

			if (trk0->size())
			{
				trk0->MakeProjection();
				tracks->push_back(trk0);
			}
			else delete trk0;

			if (trk->size()) trk->MakeProjection();
			else
			{
				// --------------------- !!!
				//tracks->release(trk); // <--------------------- !!!!
				// --------------------- !!!
				delete trk;
			}

			mf::LogVerbatim("pma::VtxCandidate") << "  done";
		}
	}
	// --------------------- !!!
	//fAssigned.front().first->TuneFullTree(); // <--------------------- !!!!
	// --------------------- !!!
	fCenter = fAssigned.front().first->Vertices().front()->Point();
	fMse = 0.0; fMse2D = 0.0;
}

