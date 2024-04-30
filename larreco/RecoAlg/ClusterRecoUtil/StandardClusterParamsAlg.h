/** ****************************************************************************
 * @file   StandardClusterParamsAlg.h
 * @brief  Interface to class computing cluster parameters
 * @author petrillo@fnal.gov
 * @date   January 21, 2015
 * @see    StandardClusterParamsAlg.cxx
 *
 * Changes:
 * 20150121 Gianluca Petrillo (petrillo@fnal.gov)
 *   structure shaped on ClusterParamsAlg class by Andrzej Szelc
 *
 * ****************************************************************************/

#ifndef STANDARDCLUSTERPARAMSALG_H
#define STANDARDCLUSTERPARAMSALG_H

// C/C++ standard library
#include <vector>

// LArSoft libraries
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"

/// Cluster reconstruction namespace
namespace cluster {

  /**
   * @brief Algorithm collection class computing cluster parameters
   * @see ClusterParamsAlg
   *
   * This class wraps ClusterParamsAlg class, designed in the context of shower
   * reconstruction, to expose a standard ClusterParamsBaseAlg interface.
   */
  class StandardClusterParamsAlg : public ClusterParamsAlgBase {
  public:
    using Measure_t = ClusterParamsAlgBase::Measure_t;

    /// Constructor
    StandardClusterParamsAlg();

    /// Destructor
    virtual ~StandardClusterParamsAlg() = default;

    /// Restores the class to post-configuration, pre-initialization state
    void Clear() override;

    /**
     * @brief Sets the list of input hits
     * @param hits list of hits
     * @throw undefined in case of error, this method can throw (anything)
     *
     * Hits are translated into our own internal format.
     * The original hits are not used afterward, and their distruction will not
     * affect this object.
     * This method calls Clear() at the beginning (although the protocol does
     * not requires it).
     */
    void SetHitsFromPointers(util::GeometryUtilities const& gser,
                             std::vector<recob::Hit const*> const& hits) override;

    /**
     * @brief Sets the list of input hits
     * @param hits list of hits (hits will not be modified)
     * @throw undefined in case of error, this method can throw (anything)
     * @see ClusterParamsAlgBase::SetHits(std::vector<recob::Hit> const&)
     */
    void SetHits(util::GeometryUtilities const& gser, std::vector<recob::Hit> const& hits) override
    {
      ClusterParamsAlgBase::SetHits(gser, hits);
    }

    /// Set the verbosity level
    void SetVerbose(int level = 1) override;

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
    float Width(util::GeometryUtilities const& gser) override;

    /// Returns the number of input hits
    size_t NInputHits() const;

  protected:
    ClusterParamsAlg algo; ///< the actual algorithm class

  }; //class StandardClusterParamsAlg

} // namespace cluster

#endif // STANDARDCLUSTERPARAMSALG_H
