#ifndef LARRECO_QUADVTX_HEATMAP_H
#define LARRECO_QUADVTX_HEATMAP_H

#include <memory>
#include <vector>

class TH2F;

namespace quad {
  class HeatMap {
  public:
    HeatMap(int _Nz, double _minz, double _maxz, int _Nx, double _minx, double _maxx);

    std::unique_ptr<TH2F> AsTH2() const;

    double ZBinCenter(int iz) const { return minz + (iz + .5) * (maxz - minz) / Nz; }
    double XBinCenter(int ix) const { return minx + (ix + .5) * (maxx - minx) / Nx; }
    int ZToBin(double z) const { return fast_floor((z - minz) / (maxz - minz) * Nz); }
    int XToBin(double x) const { return fast_floor((x - minx) / (maxx - minx) * Nx); }

    const double minz, minx, maxz, maxx;
    const int Nx, Nz;

    std::vector<float> map;

  private:
    // The rounding functions in std:: are surprisingly slow
    inline int fast_floor(double x) const { return int(x + 100000) - 100000; }
  };
}

#endif
