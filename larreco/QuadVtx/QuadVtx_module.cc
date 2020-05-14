// Chris Backhouse - c.backhouse@ucl.ac.uk - Oct 2019

#include "larreco/QuadVtx/HeatMap.h"

// C/C++ standard libraries
#include <iostream>
#include <random>
#include <string>

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "TGraph.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace quad {

  // ---------------------------------------------------------------------------
  struct Pt2D {
    Pt2D(double _x, double _z, int _view, double _energy)
      : x(_x), z(_z), view(_view), energy(_energy)
    {}

    bool
    operator<(const Pt2D& p) const
    {
      return z < p.z;
    }

    double x, z;
    int view;
    double energy;
  };

  // ---------------------------------------------------------------------------
  struct Line2D {
    Line2D(const Pt2D& a, const Pt2D& b)
      : m((b.x - a.x) / (b.z - a.z))
      , c(b.x - m * b.z)
      , // w(a.energy * b.energy),
      minz(std::min(a.z, b.z))
      , maxz(std::max(a.z, b.z))
    {
      assert(a.z != b.z); // no vertical lines
    }

    bool
    operator<(const Line2D& l) const
    {
      return m < l.m;
    }

    float m, c, /*w,*/ minz, maxz;
  };

  // ---------------------------------------------------------------------------
  class QuadVtx : public art::EDProducer {
  public:
    explicit QuadVtx(const fhicl::ParameterSet& pset);

    void beginJob() override;
    void produce(art::Event& evt) override;

  private:
    bool FindVtx(const std::vector<recob::Hit>& hits, TVector3& vtx, int evt) const;

    std::string fHitLabel;

    bool fSavePlots;

    const detinfo::DetectorProperties* detprop;
    const geo::GeometryCore* geom;
  };

  DEFINE_ART_MODULE(QuadVtx)

  // ---------------------------------------------------------------------------
  QuadVtx::QuadVtx(const fhicl::ParameterSet& pset)
    : EDProducer(pset)
    , fHitLabel(pset.get<std::string>("HitLabel"))
    , fSavePlots(pset.get<bool>("SavePlots"))
  {
    produces<std::vector<recob::Vertex>>();
  }

  // ---------------------------------------------------------------------------
  void
  QuadVtx::beginJob()
  {
    detprop = art::ServiceHandle<const detinfo::DetectorPropertiesService>()->provider();
    geom = art::ServiceHandle<const geo::Geometry>()->provider();
  }

  // ---------------------------------------------------------------------------
  // x = m*z+c. z1 and z2 are the two intercepts in case of returning true
  bool
  IntersectsCircle(float m, float c, float z0, float x0, float R, float& z1, float& z2)
  {
    // Change to the frame where (z0, x0) = (0, 0)
    c += m * z0 - x0;

    // z^2 + (m*z+c)^2 = R^2
    const float A = 1 + cet::square(m);
    const float B = 2 * m * c;
    const float C = cet::square(c) - cet::square(R);

    double desc = cet::square(B) - 4 * A * C;

    if (desc < 0) return false;

    desc = sqrt(desc);

    z1 = (-B - desc) / (2 * A);
    z2 = (-B + desc) / (2 * A);

    // Back to the original frame
    z1 += z0;
    z2 += z0;

    return true;
  }

  // ---------------------------------------------------------------------------
  void
  LinesFromPoints(const std::vector<Pt2D>& pts,
                  std::vector<Line2D>& lines,
                  float z0 = 0,
                  float x0 = 0,
                  float R = -1)
  {
    constexpr size_t kMaxLines = 10 * 1000 * 1000; // This is 150MB of lines...

    const size_t product = (pts.size() * (pts.size() - 1)) / 2;
    const int stride = product / kMaxLines + 1;

    lines.reserve(std::min(product, kMaxLines));

    for (int offset = 0; offset < stride; ++offset) {
      for (unsigned int i = 0; i < pts.size(); ++i) {
        for (unsigned int j = i + offset + 1; j < pts.size(); j += stride) {
          const Line2D l(pts[i], pts[j]);

          if (isinf(l.m) || isnan(l.m) || isinf(l.c) || isnan(l.c)) continue;

          if (R > 0) {
            float z1, z2;
            if (!IntersectsCircle(l.m, l.c, z0, x0, 2.5, z1, z2)) continue;
            if (l.minz < z1 && l.minz < z2 && l.maxz > z1 && l.maxz > z2) continue;
          }

          lines.push_back(std::move(l));
          if (lines.size() == kMaxLines) goto end; // break out of 3 loops
        }
      }
    }

  end:

    lines.shrink_to_fit();

    mf::LogInfo() << "Made " << lines.size() << " lines using stride " << stride
                  << " to fit under cap of " << kMaxLines << std::endl;

    // Lines are required to be sorted by gradient for a later optimization
    std::sort(lines.begin(), lines.end());
  }

  // ---------------------------------------------------------------------------
  inline bool
  CloseAngles(float ma, float mb)
  {
    const float cosCrit = cos(10 * M_PI / 180.);
    const float dot = 1 + ma * mb; // (1, ma)*(1, mb)
    return cet::square(dot) > (1 + cet::square(ma)) * (1 + cet::square(mb)) * cet::square(cosCrit);
  }

  // ---------------------------------------------------------------------------
  void
  MapFromLines(const std::vector<Line2D>& lines, HeatMap& hm)
  {
    // This maximum is driven by runtime
    constexpr size_t kMaxPts = 10 * 1000 * 1000;

    unsigned int j0 = 0;
    unsigned int jmax = 0;

    long npts = 0;
    for (unsigned int i = 0; i + 1 < lines.size(); ++i) {
      const Line2D a = lines[i];

      j0 = std::max(j0, i + 1);
      while (j0 < lines.size() && CloseAngles(a.m, lines[j0].m))
        ++j0;
      jmax = std::max(jmax, j0);
      while (jmax < lines.size() && !CloseAngles(a.m, lines[jmax].m))
        ++jmax;

      npts += jmax - j0;
    }

    const size_t product = (lines.size() * (lines.size() - 1)) / 2;
    const int stride = npts / kMaxPts + 1;

    mf::LogInfo() << "Combining lines to points with stride " << stride << std::endl;

    mf::LogInfo() << npts << " cf " << product << " ie " << double(npts) / product << std::endl;

    j0 = 0;
    jmax = 0;

    for (unsigned int i = 0; i + 1 < lines.size(); ++i) {
      const Line2D a = lines[i];

      j0 = std::max(j0, i + 1);
      while (j0 < lines.size() && CloseAngles(a.m, lines[j0].m))
        ++j0;
      jmax = std::max(jmax, j0);
      while (jmax < lines.size() && !CloseAngles(a.m, lines[jmax].m))
        ++jmax;

      for (unsigned int j = j0; j < jmax; j += stride) {
        const Line2D& b = lines[j];

        // x = mA * z + cA = mB * z + cB
        const float z = (b.c - a.c) / (a.m - b.m);
        const float x = a.m * z + a.c;

        // No solutions within a line
        if ((z < a.minz || z > a.maxz) && (z < b.minz || z > b.maxz)) {
          const int iz = hm.ZToBin(z);
          const int ix = hm.XToBin(x);
          if (iz >= 0 && iz < hm.Nz && ix >= 0 && ix < hm.Nx) { hm.map[iz * hm.Nx + ix] += stride; }
        }
      } // end for i
    }
  }

  // ---------------------------------------------------------------------------
  // Assumes that all three maps have the same vertical stride
  TVector3
  FindPeak3D(const std::vector<HeatMap>& hs, const std::vector<TVector3>& dirs) noexcept
  {
    assert(hs.size() == 3);
    assert(dirs.size() == 3);

    const int Nx = hs[0].Nx;

    TMatrixD M(2, 2);
    M(0, 0) = dirs[0].Y();
    M(0, 1) = dirs[0].Z();
    M(1, 0) = dirs[1].Y();
    M(1, 1) = dirs[1].Z();

    // Singular, and stupid setup of exceptions means we can't test any other way
    if (M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1) == 0) return TVector3(0, 0, 0);

    M.Invert();

    float bestscore = -1;
    TVector3 bestr;

    // Accumulate some statistics up front that will enable us to optimize
    std::vector<float> colMax[3];
    for (int view = 0; view < 3; ++view) {
      colMax[view].resize(hs[view].Nz);
      for (int iz = 0; iz < hs[view].Nz; ++iz) {
        colMax[view][iz] = *std::max_element(&hs[view].map[Nx * iz], &hs[view].map[Nx * (iz + 1)]);
      }
    }

    for (int iz = 0; iz < hs[0].Nz; ++iz) {
      const float z = hs[0].ZBinCenter(iz);
      const float bonus = 1; // works badly... exp((hs[0].maxz-z)/1000.);

      for (int iu = 0; iu < hs[1].Nz; ++iu) {
        const float u = hs[1].ZBinCenter(iu);
        // r.Dot(d0) = z && r.Dot(d1) = u
        TVectorD p(2);
        p(0) = z;
        p(1) = u;
        const TVectorD r = M * p;
        const float v = r[0] * dirs[2].Y() + r[1] * dirs[2].Z();
        const int iv = hs[2].ZToBin(v);
        if (iv < 0 || iv >= hs[2].Nz) continue;
        const double y = r(0);

        // Even if the maxes were all at the same x we couldn't beat the record
        if (colMax[0][iz] + colMax[1][iu] + colMax[2][iv] < bestscore) continue;

        // Attempt to micro-optimize the dx loop below
        const float* __restrict__ h0 = &hs[0].map[Nx * iz];
        const float* __restrict__ h1 = &hs[1].map[Nx * iu];
        const float* __restrict__ h2 = &hs[2].map[Nx * iv];

        int bestix = -1;
        for (int ix = 1; ix < Nx - 1; ++ix) {
          const float score = bonus * (h0[ix] + h1[ix] + h2[ix]);

          if (score > bestscore) {
            bestscore = score;
            bestix = ix;
          }
        } // end for dx

        if (bestix != -1) { bestr = TVector3(hs[0].XBinCenter(bestix), y, z); }
      } // end for u
    }   // end for z

    return bestr;
  }

  // ---------------------------------------------------------------------------
  void
  GetPts2D(const std::vector<recob::Hit>& hits,
           std::vector<std::vector<Pt2D>>& pts,
           std::vector<TVector3>& dirs,
           const geo::GeometryCore* geom,
           const detinfo::DetectorProperties* detprop)
  {
    pts.resize(3); // 3 views

    TVector3 dirZ(0, 0, 1);
    TVector3 dirU, dirV;

    for (const recob::Hit& hit : hits) {
      const geo::WireID wire = hit.WireID();

      const double xpos = detprop->ConvertTicksToX(hit.PeakTime(), wire);

      const TVector3 r0 = geom->WireEndPoints(wire).start();
      const TVector3 r1 = geom->WireEndPoints(wire).end();

      const double energy = hit.Integral();

      if (geom->View(hit.Channel()) == geo::kZ) {
        pts[0].emplace_back(xpos, r0.z(), 0, energy);
        continue;
      }

      // Compute the direction perpendicular to the wires
      TVector3 perp = (r1 - r0).Unit();
      perp = TVector3(0, -perp.z(), perp.y());
      // We want to ultimately have a positive z component in "perp"
      if (perp.z() < 0) perp *= -1;

      // TODO check we never get a 4th view a-la the bug we had in the 3D version

      // The "U" direction is the first one we see
      if (dirU.Mag2() == 0) { dirU = perp; }
      else if (dirV.Mag2() == 0 && fabs(dirU.Dot(perp)) < 0.99) {
        // If we still need a "V" and this direction differs from "U"
        dirV = perp;
      }

      // Hits belong to whichever view their perpendicular vector aligns with
      if (fabs(dirU.Dot(perp)) > 0.99) { pts[1].emplace_back(xpos, r0.Dot(dirU), 1, energy); }
      else {
        pts[2].emplace_back(xpos, r0.Dot(dirV), 2, energy);
      }
    } // end for hits

    dirs = {dirZ, dirU, dirV};

    std::default_random_engine gen;

    // In case we need to sub-sample they should be shuffled
    for (int view = 0; view < 3; ++view) {
      std::shuffle(pts[view].begin(), pts[view].end(), gen);
    }
  }

  // ---------------------------------------------------------------------------
  bool
  QuadVtx::FindVtx(const std::vector<recob::Hit>& hits, TVector3& vtx, int evt) const
  {
    if (hits.empty()) return false;

    std::vector<std::vector<Pt2D>> pts;
    std::vector<TVector3> dirs;

    GetPts2D(hits, pts, dirs, geom, detprop);

    double minx = +1e9;
    double maxx = -1e9;
    double minz[3] = {+1e9, +1e9, +1e9};
    double maxz[3] = {-1e9, -1e9, -1e9};
    for (int view = 0; view < 3; ++view) {
      for (const Pt2D& p : pts[view]) {
        minx = std::min(minx, p.x);
        maxx = std::max(maxx, p.x);
        minz[view] = std::min(minz[view], p.z);
        maxz[view] = std::max(maxz[view], p.z);
      }
    }

    // Add some padding
    for (int view = 0; view < 3; ++view) {
      minz[view] -= 100;
      maxz[view] += 100;
    }
    minx -= 20;
    maxx += 20;

    // Don't allow the vertex further downstream in z (view 0) than 25% of the
    // hits.
    std::vector<float> zs;
    zs.reserve(pts[0].size());
    for (const Pt2D& p : pts[0])
      zs.push_back(p.z);
    auto mid = zs.begin() + zs.size() / 4;
    if (mid != zs.end()) {
      std::nth_element(zs.begin(), mid, zs.end());
      maxz[0] = *mid;
    }

    std::vector<HeatMap> hms;
    hms.reserve(3);
    for (int view = 0; view < 3; ++view) {
      if (pts[view].empty()) return false;

      std::vector<Line2D> lines;
      LinesFromPoints(pts[view], lines);

      if (lines.empty()) return false;

      // Approximately cm bins
      hms.emplace_back(maxz[view] - minz[view], minz[view], maxz[view], maxx - minx, minx, maxx);
      MapFromLines(lines, hms.back());
    } // end for view

    vtx = FindPeak3D(hms, dirs);

    std::vector<HeatMap> hms_zoom;
    hms_zoom.reserve(3);
    for (int view = 0; view < 3; ++view) {
      const double x0 = vtx.X();
      const double z0 = vtx.Dot(dirs[view]);

      std::vector<Line2D> lines;
      LinesFromPoints(pts[view], lines, z0, x0, 2.5);

      if (lines.empty()) return false; // How does this happen??

      // mm granularity
      hms_zoom.emplace_back(50, z0 - 2.5, z0 + 2.5, 50, x0 - 2.5, x0 + 2.5);

      MapFromLines(lines, hms_zoom.back());
    }

    vtx = FindPeak3D(hms_zoom, dirs);

    if (fSavePlots) {
      art::TFileDirectory evt_dir =
        art::ServiceHandle<art::TFileService>()->mkdir(TString::Format("evt%d", evt).Data());

      for (int view = 0; view < 3; ++view) {
        art::TFileDirectory view_dir = evt_dir.mkdir(TString::Format("view%d", view).Data());

        TGraph* gpts = view_dir.makeAndRegister<TGraph>("hits", "");
        for (const Pt2D& p : pts[view])
          gpts->SetPoint(gpts->GetN(), p.z, p.x);

        view_dir.makeAndRegister<TH2F>("hmap", "", *hms[view].AsTH2());

        view_dir.makeAndRegister<TH2F>("hmap_zoom", "", *hms_zoom[view].AsTH2());

        const double x = vtx.X();
        const double z = vtx.Dot(dirs[view]);
        view_dir.makeAndRegister<TGraph>("vtx3d", "", 1, &z, &x);
      } // end for view
    }   // end if saving plots

    return true;
  }

  // ---------------------------------------------------------------------------
  void
  QuadVtx::produce(art::Event& evt)
  {
    auto vtxcol = std::make_unique<std::vector<recob::Vertex>>();

    art::Handle<std::vector<recob::Hit>> hits;
    evt.getByLabel(fHitLabel, hits);

    TVector3 vtx;
    if (FindVtx(*hits, vtx, evt.event())) {
      vtxcol->emplace_back(
        recob::Vertex::Point_t(vtx.X(), vtx.Y(), vtx.Z()), recob::Vertex::SMatrixSym33(), 0, 0);
    }

    evt.put(std::move(vtxcol));
  }

} // end namespace quad
