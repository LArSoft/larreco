/** ****************************************************************************
 * @file   LazyClusterParamsAlg.h
 * @brief  Algorithm class inheriting pre-computed results
 * @author petrillo@fnal.gov
 * @date   February 3rd, 2015
 * @see    LazyClusterParamsAlg.cxx
 *
 * ****************************************************************************/

#ifndef LAZYCLUSTERPARAMSALG_H
#define LAZYCLUSTERPARAMSALG_H

// C/C++ standard library
#include <cstddef>
#include <vector>

// LArSoft libraries
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"
namespace recob {
  class Hit;
}

/// Cluster reconstruction namespace
namespace cluster {

  class cluster_params;

  /**
   * @brief Algorithm class inheriting cluster parameters
   * @see ClusterParamsAlg
   *
   * This class wraps ClusterParamsAlg class, designed in the context of shower
   * reconstruction, to expose a standard ClusterParamsBaseAlg interface.
   * All the information is supposed to have been computed already.
   *
   * In addition to the standard interface, GetParams() is also available.
   */
  class LazyClusterParamsAlg : public ClusterParamsAlgBase {
  public:
    using Measure_t = ClusterParamsAlgBase::Measure_t;

    /// Constructor: references to the parameters (no copy is performed!)
    LazyClusterParamsAlg(cluster_params const& new_params) : params(new_params) {}

    /// Restores the class to post-configuration, pre-initialization state; dummy
    void Clear() override {}

    /**
     * @brief Sets the list of input hits
     * @param hits list of hits
     * @throw undefined in case of error, this method can throw (anything)
     *
     * The parameters have already been computed.
     * This function is dummy.
     */
    void SetHits(util::GeometryUtilities const& gser,
                 std::vector<recob::Hit const*> const& hits) override
    {}

    //@{
    /**
     * @brief Computes the charge on the first and last wire of the track
     * @return the charge in ADC counts, with uncertainty
     *
     * The implementation in ClusterParamsAlg provides an estimation of the
     * charge collected in the first or last 1 cm of the cluster, using a linear
     * fit on the deposited charge to reduce fluctuations.
     */
    Measure_t StartCharge(util::GeometryUtilities const& gser) override;
    Measure_t EndCharge(util::GeometryUtilities const& gser) override;
    //@}

    //@{
    /**
     * @brief Computes the angle of the cluster
     * @return angle of the cluster in the wire x time space, in radians
     *
     * Uses the coordinates from the hits, weighted by charge (Hit::Integral())
     * to compute a slope in the homogenized wire x time space.
     * The homogenized space has both wires and ticks converted into a distance
     * (by using detector parameters: wire pitch and drift velocity).
     *
     * The angle is in the @f$ [ -\pi, \pi ] @f$ range, with 0 corresponding to
     * a cluster parallel to the wire plane and @f$ \pi @f$ to a cluster
     * orthogonal to the wire plane, going farther from it.
     *
     * @note Both the methods return the same value.
     */
    Measure_t StartAngle() override;
    Measure_t EndAngle() override;
    //@}

    //@{
    /**
     * @brief Computes the opening angle at the start or end of the cluster
     * @return angle at the start of the cluster, in radians
     *
     * This algorithm returns an opening angle after weighting the hits by
     * their charge (as defined bu Hit::Integral());
     */
    Measure_t StartOpeningAngle() override;
    Measure_t EndOpeningAngle() override;
    //@}

    /// @name Cluster charge
    /// @{
    /**
     * @brief Computes the total charge of the cluster from Hit::Integral()
     * @return total charge of the cluster, in ADC count units
     * @see IntegralStdDev(), SummedADC()
     *
     * ClusterParamsAlg computes the sum from all hits.
     */
    Measure_t Integral() override;

    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see Integral()
     *
     * ClusterParamsAlg computes the standard deviation of the sample of charges
     * from all hits.
     * Hit charge is obtained by recob::Hit::Integral().
     */
    Measure_t IntegralStdDev() override;

    /**
     * @brief Computes the total charge of the cluster from Hit::SummedADC()
     * @return total charge of the cluster, in ADC count units
     * @see SummedADCStdDev(), Integral()
     *
     * ClusterParamsAlg computes the sum from all hits.
     */
    Measure_t SummedADC() override;

    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see SummedADC()
     *
     * ClusterParamsAlg computes the standard deviation of the sample of charges
     * from all hits.
     * Hit charge is obtained by recob::Hit::SummedADC().
     */
    Measure_t SummedADCStdDev() override;

    /// @}

    /// Returns the number of hits in the cluster
    size_t NHits() override;

    /**
     * @brief Fraction of wires in the cluster with more than one hit
     * @return fraction of wires with more than one hit, or 0 if no wires
     *
     * Returns a quantity defined as NMultiHitWires / NWires,
     * where NWires is the number of wires hosting at least one hit of this
     * cluster, and NMultiHitWires is the number of wires which have more
     * than just one hit.
     */
    float MultipleHitDensity() override;

    /**
     * @brief Computes the width of the cluster
     * @return width of the cluster
     *
     * @todo provide a description of the algorithm by words
     */
    float Width(util::GeometryUtilities const&) override;

    /// Returns the original precomputed parameters
    cluster_params const& GetParams() const { return params; }

  protected:
    cluster_params const& params; ///< the parameters, already computed

  }; //class LazyClusterParamsAlg

} // namespace cluster

#endif // LAZYCLUSTERPARAMSALG_H
