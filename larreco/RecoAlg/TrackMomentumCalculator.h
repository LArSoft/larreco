/// \file  TrackMomentumCalculator.h
//  \author  sowjanyag@phys.ksu.edu

#ifndef TrackMomentumCalculator_H
#define TrackMomentumCalculator_H

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Track.h"

#include "TGraph.h"
#include "TVector3.h"

#include <optional>
#include <tuple>
#include <vector>

class TPolyLine3D;

namespace trkf {

  class TrackMomentumCalculator {
  public:
    /**
    * @brief  Constructor
    *
    * Parameters are relevant for Multiple Coulomb Scattering
    *
    * @param  minLength minimum length in cm of tracks (length here is based on the amout of segments)
    * @param  maxLength maximum length in cm of tracks (length here is based on the amout of segments)
    * @param  stes_size size in cm of each segment to compute scattering
    */
    TrackMomentumCalculator(double minLength = 100.0,
                            double maxLength = 1350.0,
                            double steps_size = 10.);

    double GetTrackMomentum(double trkrange, int pdg) const;
    /**
    * @brief  Calculate muon momentum (MeV) using multiple coulomb scattering. Chi2 minimization of the Highland formula
    *
    * @param  trk the muon track
    * @param  checkValidPoints rather take into account only valid points or not
    * @param  maxMomentum_MeV maximum momentum in MeV for the minimization
    *
    * TODO: Add better description of the steps done.
    *
    * @return momentum in MeV
    */
    double GetMomentumMultiScatterChi2(art::Ptr<recob::Track> const& trk,
                                       const bool checkValidPoints = false,
                                       const int maxMomentum_MeV = 7500);
    /**
    * @brief  Calculate muon momentum (MeV) using multiple coulomb scattering by log likelihood
    *
    * @param  trk the muon track
    * @param  checkValidPoints rather take into account only valid points or not
    * @param  maxMomentum_MeV maximum momentum in MeV for the minimization
    * @param  MomentumStep_MeV energy steps for minimization
    * @param  max_resolution maximum angular resolution for fit. Setting to zero will cause the fit only over momentum and fixed resolution of 2 mrad
    *
    * TODO: Add better description of the steps done
    *
    * @return momentum in MeV
    */
    double GetMomentumMultiScatterLLHD(art::Ptr<recob::Track> const& trk,
                                       const bool checkValidPoints = false,
                                       const int maxMomentum_MeV = 7500,
                                       const int MomentumStep_MeV = 10,
                                       const int max_resolution = 0);
    double GetMuMultiScatterLLHD3(art::Ptr<recob::Track> const& trk, bool dir);
    TVector3 GetMultiScatterStartingPoint(art::Ptr<recob::Track> const& trk);

  private:
    bool plotRecoTracks_(std::vector<float> const& xxx,
                         std::vector<float> const& yyy,
                         std::vector<float> const& zzz);

    /**
    * @brief Computes the vector with most scattering inside a segment with size steps_size
    * @param segx, segy, segz segments points
    * @param segnx, segny, segnz vector components to be filled
    * @param vector used to control points to be used at segments
    *
    */
    void compute_max_fluctuation_vector(const std::vector<float> segx,
                                        const std::vector<float> segy,
                                        const std::vector<float> segz,
                                        std::vector<float>& segnx,
                                        std::vector<float>& segny,
                                        std::vector<float>& segnz,
                                        std::vector<float>& vx,
                                        std::vector<float>& vy,
                                        std::vector<float>& vz);
    /**
    * \struct Segments
    * @brief Struct to store segments.
    * x, y and z are the 3D points of the segment
    * nx, ny, nz forms the vector that point to the direction of most scattering
    * L is the length of the segment, using steps_size
    *
    */
    struct Segments {
      std::vector<float> x, nx;
      std::vector<float> y, ny;
      std::vector<float> z, nz;
      std::vector<float> L;
    };

    /**
    * @brief Split tracks into segments to calculate the scattered angle later. Check DOI 10.1088/1748-0221/12/10/P10010
    *
    * @param xxx 3D reconstructed points x-axis
    * @param yyy 3D reconstructed points y-axiy
    * @param zzz 3D reconstructed points z-axiz
    * @seg_size Segments size defined in class constructor
    *
    * TODO: Add better description of steps
    *
    * @return Segments
    */
    std::optional<Segments> getSegTracks_(std::vector<float> const& xxx,
                                          std::vector<float> const& yyy,
                                          std::vector<float> const& zzz,
                                          double seg_size);

    /**
    * @brief Gets the scattered angle RMS for a all segments
    *
    * @param segments segments computed
    * @param thick is the steps_size
    *
    * TODO: Add better description of steps
    *
    * @return tuple with mean value, rms and error of rms
    */
    std::tuple<double, double, double> getDeltaThetaRMS_(Segments const& segments,
                                                         double thick) const;

    /**
    * @brief Gets the scatterd angle for all the segments
    *
    * @param ei, ej will be filled with energy lost at point i and j
    * @param th delta theta to be filled
    * @param ind selection of scattering plane
    * @param segments
    * @param thick is the steps_size
    *
    * @return sucess or failure
    */
    int getDeltaThetaij_(std::vector<float>& ei,
                         std::vector<float>& ej,
                         std::vector<float>& th,
                         std::vector<float>& ind,
                         Segments const& segments,
                         double thick) const;

    /**
    * @brief chi square minizer using Minuit2, it will minize (xx-Q)/s
    *
    * @param xx
    * @param Q
    * @param s
    *
    * @return Momentum in MeV
    */
    double my_g(double xx, double Q, double s) const;

    /**
    * @brief  Minimizer of log likelihood for scattered angle
    *
    * @param  dEi energy at step i
    * @param  dEj energy at step j
    * @param  dthij scattered angle between points i and j
    * @param  ind selection of scattering plane
    * @param  x0 momentum to be fitted
    * @param  x1 resolution to be fitted
    *
    * TODO: Add better description of steps
    *
    * @return momentum in MeV
    */
    double my_mcs_llhd(std::vector<float> const& dEi,
                       std::vector<float> const& dEj,
                       std::vector<float> const& dthij,
                       std::vector<float> const& ind,
                       double x0,
                       double x1) const;

    float seg_stop{-1.};
    int n_seg{};

    float x_seg[50000];
    float y_seg[50000];
    float z_seg[50000];

    /**
    * @brief Gets angle between two vy and vz
    *
    * @param vz
    * @param vy
    *
    * @return angle in mrad
    */
    double find_angle(double vz, double vy) const;

    int n_steps{6};
    std::vector<float> steps;

    double minLength;
    double maxLength;
    double steps_size;
    double rad_length{14.0};

    // The following are objects that are created but not drawn or
    // saved.  This class should consider accepting a "debug"
    // parameter where if it is specified, then the graphs will be
    // created; otherwise, their creation is unnecessary and impedes
    // efficiency.
    //
    // N.B. TPolyLine3D objects are owned by ROOT, and we thus refer
    // to them by pointer.  It is important that 'delete' is not
    // called on the TPolyLine3D pointers during destruction of a
    // TrackMomentumCalculator object.
    TPolyLine3D* gr_reco_xyz{nullptr};
    TGraph gr_reco_xy{};
    TGraph gr_reco_yz{};
    TGraph gr_reco_xz{};

    TPolyLine3D* gr_seg_xyz{nullptr};
    TGraph gr_seg_xy{};
    TGraph gr_seg_yz{};
    TGraph gr_seg_xz{};
  };

} // namespace trkf

#endif // TrackMomentumCalculator_H
