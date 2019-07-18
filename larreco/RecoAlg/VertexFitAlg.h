//////////////////////////////////////////////////////////////////////
///
/// VertexFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 3D vertex given a set of track hits
///
////////////////////////////////////////////////////////////////////////
#ifndef VERTEXFITALG_H
#define VERTEXFITALG_H

#include <vector>

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/VertexFitMinuitStruct.h"

// ROOT includes
#include "RtypesCore.h"
class TVector3;

namespace trkf {

  class VertexFitAlg {
    public:

    void VertexFit(std::vector<std::vector<geo::WireID>> const& hitWID,
                      std::vector<std::vector<double>> const& hitX,
                      std::vector<std::vector<double>> const& hitXErr,
                      TVector3& VtxPos, TVector3& VtxPosErr,
                      std::vector<TVector3>& TrkDir, std::vector<TVector3>& TrkDirErr,
                      float& ChiDOF) const;

    // Variables for minuit.
    static VertexFitMinuitStruct fVtxFitMinStr;

    static void fcnVtxPos(Int_t &, Double_t *, Double_t &fval, double *par, Int_t flag);

    private:

    art::ServiceHandle<geo::Geometry const> geom;


  }; // class VertexFitAlg

} // namespace trkf

#endif // ifndef VERTEXFITALG_H
