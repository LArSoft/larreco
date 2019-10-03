#include "larreco/QuadVtx/HeatMap.h"

#include "TH2.h"

namespace quad
{
  // -------------------------------------------------------------------------
  HeatMap::HeatMap(int _Nz, float _minz, float _maxz,
                   int _Nx, float _minx, float _maxx)
    : minz(_minz), minx(_minx), maxz(_maxz), maxx(_maxx), Nx(_Nx), Nz(_Nz),
      map(Nz*Nx, 0)
  {
  }

  // -------------------------------------------------------------------------
  std::unique_ptr<TH2F> HeatMap::AsTH2() const
  {
    auto h = std::make_unique<TH2F>("", "", Nz, minz, maxz, Nx, minx, maxx);
    for(unsigned int i = 0; i < map.size(); ++i){
      h->SetBinContent((i/Nx)+1, (i%Nx)+1, map[i]);
    }

    double maxz = 0;
    for(int i = 0; i < (h->GetNbinsX()+2)*(h->GetNbinsY()+2); ++i){
      maxz = std::max(maxz, h->GetBinContent(i));
    }

    // Slim out bins with low population (reduces output file size)
    for(int i = 0; i < (h->GetNbinsX()+2)*(h->GetNbinsY()+2); ++i){
      if(h->GetBinContent(i) < .05*maxz) h->SetBinContent(i, 0);
    }

    h->Sumw2(false); // Drop all the errors

    return h;
  }
}
