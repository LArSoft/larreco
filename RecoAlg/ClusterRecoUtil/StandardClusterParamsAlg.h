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
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"



/// Cluster reconstruction namespace
namespace cluster {
  
  /**
   * @brief Algorithm collection class computing cluster parameters
   * @see ClusterParamsAlg
   *
   * This class wraps ClusterParamsAlg class, designed in the context of shower
   * reconstruction, to expose a standard ClusterParamsBaseAlg interface.
   */
  class StandardClusterParamsAlg: public ClusterParamsAlgBase {
      public:
    using Measure_t = ClusterParamsAlgBase::Measure_t;
    
    /// Constructor
    StandardClusterParamsAlg();
    
    /// Destructor
    virtual ~StandardClusterParamsAlg() = default;
    
    
    /// Restores the class to post-configuration, pre-initialization state
    virtual void Clear() override;
    
    
    /**
     * @brief Sets the list of input hits
     * @param hits list of hits
     * @throw undefined in case of error, this method can throw (anything)
     * 
     * Hits are translated into our own internal format.
     * The original hits are not used afterward, and can disappear.
     * This method calls Clear() at the beginning.
     */
    virtual void SetHits(std::vector<recob::Hit const*> const& hits) override;
    
    
    /// Set the verbosity level
    virtual void SetVerbose(int level = 1) override;
    
    
    //@{
    /**
     * @brief Computes the charge on the first and last wire of the track
     * @return the charge in ADC counts, with uncertainty
     * 
     * The implementation in ClusterParamsAlg provides an estimation of the
     * charge collected in the first or last 1 cm of the cluster, using a linear
     * fit on the deposited charge to reduce fluctuations.
     */
    virtual Measure_t StartCharge() override;
    virtual Measure_t EndCharge() override;
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
    virtual Measure_t StartAngle() override;
    virtual Measure_t EndAngle() override;
    //@}
    
    //@{
    /**
     * @brief Computes the opening angle at the start or end of the cluster
     * @return angle at the start of the cluster, in radians
     * 
     * This algorithm returns an opening angle after weighting the hits by
     * their charge (as defined bu Hit::Integral());
     */
    virtual Measure_t StartOpeningAngle() override;
    virtual Measure_t EndOpeningAngle() override;
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
    virtual Measure_t Integral() override;
    
    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see Integral()
     * 
     * ClusterParamsAlg computes the standard deviation of the sample of charges
     * from all hits.
     * Hit charge is obtained by recob::Hit::Integral().
     */
    virtual Measure_t IntegralStdDev() override;
    
    /**
     * @brief Computes the total charge of the cluster from Hit::SummedADC()
     * @return total charge of the cluster, in ADC count units
     * @see SummedADCStdDev(), Integral()
     * 
     * ClusterParamsAlg computes the sum from all hits.
     */
    virtual Measure_t SummedADC() override;
    
    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see SummedADC()
     * 
     * ClusterParamsAlg computes the standard deviation of the sample of charges
     * from all hits.
     * Hit charge is obtained by recob::Hit::SummedADC().
     */
    virtual Measure_t SummedADCStdDev() override;
    
    /// @}
    
    
    /// Returns the number of hits in the cluster
    virtual size_t NHits() override;
    
    
    /**
     * @brief Computes the number of wires divided by the number of cluster hits
     * @return number of wires divided by the number of cluster hits
     * @see NHits()
     * 
     * @todo Ask for definition
     */
    virtual float NWiresOverNHits() override;
    
    /**
     * @brief Computes the width of the cluster
     * @return width of the cluster
     * 
     * @todo provide a description of the algorithm by words
     */
    virtual float Width() override;
    
    
      protected:
    ClusterParamsAlg algo; ///< the actual algorithm class
    
  }; //class StandardClusterParamsAlg
  
} // namespace cluster

#endif // STANDARDCLUSTERPARAMSALG_H
