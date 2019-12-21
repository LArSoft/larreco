// Christopher Backhouse - bckhouse@fnal.gov

#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larsim/MCCheater/BackTrackerService.h"

#include "TGraph.h"
#include "TString.h"
#include "TVirtualPad.h"

namespace reco3d
{

class PlotSpacePoints: public art::EDAnalyzer
{
public:
  explicit PlotSpacePoints(const fhicl::ParameterSet& pset);

private:
  void analyze(const art::Event& evt);

  void Plot(const std::vector<recob::SpacePoint>& pts,
            const std::string& suffix) const;

  void Plot3D(const std::vector<recob::SpacePoint>& pts,
              const std::string& suffix) const;

  std::vector<recob::SpacePoint> TrueSpacePoints(art::Handle<std::vector<recob::Hit>> hits) const;

  art::InputTag fSpacePointTag;

  std::string fHitLabel;

  std::string fSuffix;

  bool fPlots;
  bool fPlots3D;
  bool fPlotsTrue;
};

DEFINE_ART_MODULE(PlotSpacePoints)

// ---------------------------------------------------------------------------
PlotSpacePoints::PlotSpacePoints(const fhicl::ParameterSet& pset)
  : EDAnalyzer(pset),
    fSpacePointTag(art::InputTag(pset.get<std::string>("SpacePointLabel"),
                                 pset.get<std::string>("SpacePointInstanceLabel"))),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fSuffix(pset.get<std::string>("Suffix")),
    fPlots(pset.get<bool>("Plots")),
    fPlots3D(pset.get<bool>("Plots3D")),
    fPlotsTrue(pset.get<bool>("PlotsTrue"))
{
  if(!fSuffix.empty()) fSuffix = "_"+fSuffix;
}

// ---------------------------------------------------------------------------
void PlotSpacePoints::Plot(const std::vector<recob::SpacePoint>& pts,
                           const std::string& suffix) const
{
  TGraph gZX;
  TGraph gYX;
  TGraph gZY;

  gZX.SetTitle(";z;x");
  gYX.SetTitle(";y;x");
  gZY.SetTitle(";z;y");

  for(const recob::SpacePoint& pt: pts){
    const double* xyz = pt.XYZ();
    const double x = xyz[0];
    const double y = xyz[1];
    const double z = xyz[2];
    gZX.SetPoint(gZX.GetN(), z, x);
    gYX.SetPoint(gYX.GetN(), y, x);
    gZY.SetPoint(gZY.GetN(), z, y);
  }

  if(gZX.GetN() == 0) gZX.SetPoint(0, 0, 0);
  if(gYX.GetN() == 0) gYX.SetPoint(0, 0, 0);
  if(gZY.GetN() == 0) gZY.SetPoint(0, 0, 0);

  gZX.Draw("ap");
  gPad->Print(("plots/evd"+suffix+".png").c_str());

  gYX.Draw("ap");
  gPad->Print(("plots/evd_ortho"+suffix+".png").c_str());
  gZY.Draw("ap");
  gPad->Print(("plots/evd_zy"+suffix+".png").c_str());
}

// ---------------------------------------------------------------------------
std::vector<recob::SpacePoint> PlotSpacePoints::
TrueSpacePoints(art::Handle<std::vector<recob::Hit>> hits) const
{
  std::vector<recob::SpacePoint> pts_true;

  const double err[6] = {0,};

  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  for(unsigned int i = 0; i < hits->size(); ++i){
    try{
      const std::vector<double> xyz = bt_serv->HitToXYZ(art::Ptr<recob::Hit>(hits, i));
      pts_true.emplace_back(&xyz[0], err, 0);
    }
    catch(...){} // some hits have no electrons?
  }

  return pts_true;
}

// ---------------------------------------------------------------------------
void PlotSpacePoints::Plot3D(const std::vector<recob::SpacePoint>& pts,
                             const std::string& suffix) const
{
  int frame = 0;
  for(int phase = 0; phase < 4; ++phase){
    const int Nang = 20;
    for(int iang = 0; iang < Nang; ++iang){
      const double ang = M_PI/2*iang/double(Nang);

      TGraph g;

      for(const recob::SpacePoint& p: pts){
        const double* xyz = p.XYZ();

        double x{};
        double y{};
        if(phase == 0){
          x = cos(ang)*xyz[1]+sin(ang)*xyz[2];
          y = xyz[0];
        }
        if(phase == 1){
          x = xyz[2];
          y = cos(ang)*xyz[0]+sin(ang)*xyz[1];
        }
        if(phase == 2){
          x = cos(ang)*xyz[2]-sin(ang)*xyz[0];
          y = xyz[1];
        }
        if(phase == 3){
          x = -cos(ang)*xyz[0] + sin(ang)*xyz[1];
          y = cos(ang)*xyz[1] + sin(ang)*xyz[0];
        }

        //        const double phi = phase/3.*M_PI/2 + ang/3;
        const double phi = 0;
        g.SetPoint(g.GetN(), cos(phi)*x+sin(phi)*y, cos(phi)*y-sin(phi)*x);
      }

      std::string fname = TString::Format(("anim/evd3d"+suffix+"_%03d.png").c_str(), frame++).Data();
      g.SetTitle(fname.c_str());
      if(g.GetN()) g.Draw("ap");
      gPad->Print(fname.c_str());
    }
  }
}

void PlotSpacePoints::analyze(const art::Event& evt)
{
  if(fPlots){
    art::Handle<std::vector<recob::SpacePoint>> pts;
    evt.getByLabel(fSpacePointTag, pts);

    const std::string suffix = TString::Format("%s_%d", fSuffix.c_str(), evt.event()).Data();

    if(!pts->empty()){
      Plot(*pts, suffix);

      if(fPlots3D) Plot3D(*pts, suffix);
    }
  }

  if(fPlotsTrue){
    art::Handle<std::vector<recob::Hit>> hits;
    evt.getByLabel(fHitLabel, hits);
    const std::vector<recob::SpacePoint> pts = TrueSpacePoints(hits);

    const std::string suffix = TString::Format("%s_true_%d", fSuffix.c_str(), evt.event()).Data();

    Plot(pts, suffix);

    if(fPlots3D) Plot3D(pts, suffix);
  }
}

} // namespace
