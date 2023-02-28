////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CCHITFINDERALG_H
#define CCHITFINDERALG_H

// C/C++ standard libraries
#include <memory>  // std::unique_ptr<>
#include <ostream> // std::endl
#include <vector>

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
namespace fhicl {
  class ParameterSet;
}

#include "canvas/Persistency/Provenance/Timestamp.h"
// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/RecoAlg/GausFitCache.h"

namespace hit {

  /**
   * @brief Hit finder algorithm designed to work with Cluster Crawler
   *
   * This algorithm used to store hits in a proprietary `CCHit` data structure.
   * It has now been changed to use `recob::Hit` class directly.
   * It is possible to translate the former into the latter, with one exception,
   * as follows:
   *
   *     // this is the original CCHit definition
   *     struct CCHit {
   *       float Charge;            // recob::Hit::Integral()
   *       float ChargeErr;         // recob::Hit::SigmaIntegral()
   *       float Amplitude;         // recob::Hit::PeakAmplitude()
   *       float AmplitudeErr;      // recob::Hit::SigmaPeakAmplitude()
   *       float Time;              // recob::Hit::PeakTime()
   *       float TimeErr;           // recob::Hit::SigmaPeakTime()
   *       float RMS;               // recob::Hit::RMS()
   *       float RMSErr;            // dropped
   *       float ChiDOF;            // recob::Hit::GoodnessOfFit()
   *       int   DOF;               // recob::Hit::DegreesOfFreedom()
   *       float ADCSum;            // recob::Hit::SummedADC()
   *       unsigned short WireNum;  // recob::Hit::WireID().Wire
   *       unsigned short numHits;  // recob::Hit::Multiplicity()
   *       unsigned int LoHitID;    // see below
   *       float LoTime;            // recob::Hit::StartTick()
   *       float HiTime;            // recob::Hit::EndTick()
   *       short InClus;            // dropped; see below
   *       geo::WireID WirID;       // recob::Hit::WireID()
   *       recob::Wire const* Wire; // dropped; see below
   *     };
   *
   * The uncertainty on RMS has been dropped for good.
   *
   * The `LoHitID` member used to mean the index of the first hit in the "hit
   * train" (that is the set of hits extracted from the same region of
   * interest). That is a concept that is not portable. If your hit list is
   * still the original one as produced by this algorithm, or if at least the
   * hits from the same train are stored sorted and contiguously, for a hit with
   * index `iHit`, the equivalent value of `LoHitID` is
   * `iHit - hit.LocalIndex()`.
   *
   * There is no pointer to the wire any more in `recob::Hit`. The wire can be
   * obtained through associations, that are typically produced by the art
   * module that runs CCHitFinderAlg (e.g. `CCHitFinder`). The channel ID is
   * also directly available as `recob::Hit::Channel()`.
   */

  class CCHitFinderAlg {

  public:
    std::vector<recob::Hit> allhits;

    CCHitFinderAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitFinderAlg() = default;

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitFinder(std::vector<recob::Wire> const& Wires, art::Timestamp t);

    /// Returns (and loses) the collection of reconstructed hits
    std::vector<recob::Hit>&& YieldHits() { return std::move(allhits); }

    /// Print the fit statistics
    template <typename Stream>
    void PrintStats(Stream& out) const;

  private:
    std::vector<float> fMinPeak;
    std::vector<float> fMinRMS;
    unsigned short fMaxBumps;    // make a crude hit if > MaxBumps are found in the RAT
    unsigned short fMaxXtraHits; // max num of hits in Region Above Threshold
    float fChiSplit;             ///<Estimated noise error on the Signal
                                 //  float ChgNorm;     // Area norm for the wire we are working on

    std::vector<float> fChiNorms;
    std::vector<float> fTimeOffsets;
    std::vector<float> fChgNorms;

    raw::ChannelID_t theChannel;
    unsigned short theWireNum;
    unsigned short thePlane;

    float chinorm;
    //  float timeoff;
    static constexpr float Sqrt2Pi = 2.5066;
    static constexpr float SqrtPi = 1.7725;

    bool fUseChannelFilter;

    //    bool prt;

    art::ServiceHandle<geo::Geometry const> geom;

    // fit n Gaussians possibly with bounds setting (parmin, parmax)
    void FitNG(unsigned short nGaus, unsigned short npt, float* ticks, float* signl);
    // parameters, errors, lower limit, upper limits for FitNG
    std::vector<double> par;
    std::vector<double> parerr;
    std::vector<double> parmin;
    std::vector<double> parmax;
    float chidof;
    int dof;
    std::vector<unsigned short> bumps;

    /// exchange data about the originating wire
    class HitChannelInfo_t {
    public:
      recob::Wire const* wire;
      geo::WireID wireID;
      geo::SigType_t sigType;

      HitChannelInfo_t(recob::Wire const* w, geo::WireID wid, geo::Geometry const& geom);
    }; // HitChannelInfo_t

    // make a cruddy hit if fitting fails
    void MakeCrudeHit(unsigned short npt, float* ticks, float* signl);
    // store the hits
    void StoreHits(unsigned short TStart, unsigned short npt, HitChannelInfo_t info, float adcsum);

    // study hit finding and fitting
    bool fStudyHits;
    std::vector<short> fUWireRange, fUTickRange;
    std::vector<short> fVWireRange, fVTickRange;
    std::vector<short> fWWireRange, fWTickRange;
    void StudyHits(unsigned short flag,
                   unsigned short npt = 0,
                   float* ticks = 0,
                   float* signl = 0,
                   unsigned short tstart = 0);
    std::vector<int> bumpCnt;
    std::vector<int> RATCnt;
    std::vector<float> bumpChi;
    std::vector<float> bumpRMS;
    std::vector<int> hitCnt;
    std::vector<float> hitRMS;
    // use to determine the slope of protons
    std::vector<float> loWire;
    std::vector<float> loTime;
    std::vector<float> hiWire;
    std::vector<float> hiTime;
    bool SelRAT; // set true if a Region Above Threshold should be studied

    bool fUseFastFit; ///< whether to attempt using a fast fit on single gauss.

    std::unique_ptr<GausFitCache> FitCache; ///< a set of functions ready to be used

    typedef struct {
      unsigned int FastFits;                   ///< count of single-Gaussian fast fits
      std::vector<unsigned int> MultiGausFits; ///< multi-Gaussian stats

      void Reset(unsigned int nGaus);

      void AddMultiGaus(unsigned int nGaus);

      void AddFast() { ++FastFits; }

    } FitStats_t;

    FitStats_t FinalFitStats; ///< counts of the good fits
    FitStats_t TriedFitStats; ///< counts of the tried fits

    /**
     * @brief Performs a "fast" fit
     * @param npt number of points to be fitted
     * @param ticks tick coordinates
     * @param signl signal amplitude
     * @param params an array where the fit parameters will be stored
     * @param paramerrors an array where the fit parameter errors will be stored
     * @param chidof a variable where to store chi^2 over degrees of freedom
     * @return whether the fit was successful or not
     *
     * Note that the fit will bail out and rteurn false if any of the input
     * signal amplitudes is zero or negative.
     *
     * Also note that currently the chi^2 is not the one from comparing the
     * Gaussian to the signal, but from comparing a fitted parabola to the
     * logarithm of the signal.
     */
    static bool FastGaussianFit(unsigned short npt,
                                float const* ticks,
                                float const* signl,
                                std::array<double, 3>& params,
                                std::array<double, 3>& paramerrors,
                                float& chidof);

    static constexpr unsigned int MaxGaussians = 20;

  }; // class CCHitFinderAlg

} // namespace hit

//==============================================================================
//===  Template implementation
//===
template <typename Stream>
void hit::CCHitFinderAlg::PrintStats(Stream& out) const
{

  out << "CCHitFinderAlg fit statistics:";
  if (fUseFastFit) {
    out << "\n  fast 1-Gaussian fits: " << FinalFitStats.FastFits << " succeeded ("
        << TriedFitStats.FastFits << " tried)";
  }
  else
    out << "\n  fast 1-Gaussian fits: disabled";

  for (unsigned int nGaus = 1; nGaus < MaxGaussians; ++nGaus) {
    if (TriedFitStats.MultiGausFits[nGaus - 1] == 0) continue;
    out << "\n  " << nGaus << "-Gaussian fits: " << FinalFitStats.MultiGausFits[nGaus - 1]
        << " accepted (" << TriedFitStats.MultiGausFits[nGaus - 1] << " tried)";
  } // for nGaus
  if (TriedFitStats.MultiGausFits.back() > 0) {
    out << "\n  " << FinalFitStats.MultiGausFits.size()
        << "-Gaussian fits or higher: " << FinalFitStats.MultiGausFits.back() << " accepted ("
        << TriedFitStats.MultiGausFits.back() << " tried)";
  }
  out << std::endl;

} // CCHitFinderAlg::FitStats_t::Print()

/////////////////////////////////////////

#endif // ifndef CCHITFINDERALG_H
