#ifndef LARRECO_QUADVTX_HEATMAP_H
#define LARRECO_QUADVTX_HEATMAP_H

#include <vector>
#include <memory>

class TH2F;

namespace quad
{
  class HeatMap
  {
  public:
    HeatMap(int _Nz, float _minz, float _maxz,
            int _Nx, float _minx, float _maxx);

    std::unique_ptr<TH2F> AsTH2() const;

    inline float ZBinCenter(int iz) const {return minz + (iz+.5)*(maxz-minz)/Nz;}
    inline float XBinCenter(int ix) const {return minx + (ix+.5)*(maxx-minx)/Nx;}
    inline int ZToBin(float z) const {return fast_floor((z-minz)/(maxz-minz)*Nz);}
    inline int XToBin(float x) const {return fast_floor((x-minx)/(maxx-minx)*Nx);}

    float minz, minx, maxz, maxx;
    int Nx, Nz;

    std::vector<float> map;

  protected:
    // The rounding functions in std:: are surprisingly slow
    inline int fast_floor(float x) const {return int(x+100000)-100000;}
  };
}

#endif
