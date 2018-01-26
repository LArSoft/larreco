/**
 *  @file   PrincipalComponentsAlg.h
 * 
 *  @brief  This header file defines the interface to a principal components analysis designed to
 *          be used within the 3D clustering
 * 
 */
#ifndef PrincipalComponentsAlg_h
#define PrincipalComponentsAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

// Eigen
#include <Eigen/Dense>

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>


//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{

/**
 *  @brief  Cluster3D class
 */
class PrincipalComponentsAlg
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    PrincipalComponentsAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~PrincipalComponentsAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    void reconfigure(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief Run the Principal Components Analysis
     */
    void PCAAnalysis(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, float doca3DScl = 3.)                    const;
    
    void PCAAnalysis_3D(const reco::HitPairListPtr& hitPairList, reco::PrincipalComponents& pca, bool skeletonOnly = false)               const;
    
    void PCAAnalysis_2D(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, bool updateAvePos = false)             const;
    
    void PCAAnalysis_calc3DDocas(const reco::HitPairListPtr& hitPairVector, const reco::PrincipalComponents& pca)                         const;
    
    void PCAAnalysis_calc2DDocas(const reco::Hit2DListPtr& hit2DVector, const reco::PrincipalComponents& pca)                             const;
    
    int  PCAAnalysis_reject2DOutliers(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, float aveHitDoca)       const;
    
    int  PCAAnalysis_reject3DOutliers(const reco::HitPairListPtr& hitPairVector, const reco::PrincipalComponents& pca, float aveHitDoca) const;
    
    

private:
    /**
     *  @brief This is used to get the poca, doca and arclen along cluster axis to 2D hit
     */
    void getHit2DPocaToAxis(const Eigen::Vector3f&    axisPos,
                            const Eigen::Vector3f&    axisDir,
                            const reco::ClusterHit2D* hit2D,
                            Eigen::Vector3f&          poca,
                            float&                    arcLenAxis,
                            float&                    arcLenWire,
                            float&                    doca);
    
    float                                 m_parallel;  ///< means lines are parallel
    
    geo::Geometry*                        m_geometry;  // pointer to the Geometry service
    const detinfo::DetectorProperties*    m_detector;  // Pointer to the detector properties
};

} // namespace lar_cluster3d
#endif
