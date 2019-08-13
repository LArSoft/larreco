/**
 *  @file   Cluster3D_module.cc
 *
 *  @brief  Producer module to create 3D clusters from input hits
 *
 */

#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"

// Framework Includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"

// std includes
#include <functional>
#include <iostream>
#include <numeric>

// Eigen includes
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/Jacobi"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

PrincipalComponentsAlg::PrincipalComponentsAlg(fhicl::ParameterSet const &pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PrincipalComponentsAlg::~PrincipalComponentsAlg()
{
}

void PrincipalComponentsAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_parallel = pset.get<float>("ParallelLines", 0.00001);
    m_geometry = art::ServiceHandle<geo::Geometry const>{}.get();
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

void PrincipalComponentsAlg::getHit2DPocaToAxis(const Eigen::Vector3f&    axisPos,
                                                const Eigen::Vector3f&    axisDir,
                                                const reco::ClusterHit2D* hit2D,
                                                Eigen::Vector3f&          poca,
                                                float&                    arcLenAxis,
                                                float&                    arcLenWire,
                                                float&                    doca)
{
    // Step one is to set up to determine the point of closest approach of this 2D hit to
    // the cluster's current axis.
    // Get this wire's geometry object
    const geo::WireID&  hitID     = hit2D->WireID();
    const geo::WireGeo& wire_geom = m_geometry->WireIDToWireGeo(hitID);

    // From this, get the parameters of the line for the wire
    double          wirePos[3] = {0.,0.,0.};
    Eigen::Vector3f wireDir(wire_geom.Direction().X(),wire_geom.Direction().Y(),wire_geom.Direction().Z());

    wire_geom.GetCenter(wirePos);

    // Correct the wire position in x to set to correspond to the drift time
    float hitPeak(hit2D->getHit()->PeakTime());

    wirePos[0] = m_detector->ConvertTicksToX(hitPeak, hitID.Plane, hitID.TPC, hitID.Cryostat);

    // Get a vector from the wire position to our cluster's current average position
    Eigen::Vector3f wVec(axisPos(0)-wirePos[0], axisPos(1)-wirePos[1], axisPos(2)-wirePos[2]);

    // Get the products we need to compute the arc lengths to the distance of closest approach
    float a(axisDir.dot(axisDir));
    float b(axisDir.dot(wireDir));
    float c(wireDir.dot(wireDir));
    float d(axisDir.dot(wVec));
    float e(wireDir.dot(wVec));

    float den(a*c - b*b);

    // Parallel lines is a special case
    if (fabs(den) > m_parallel)
    {
        arcLenAxis = (b*e - c*d) / den;
        arcLenWire = (a*e - b*d) / den;
    }
    else
    {
        mf::LogDebug("Cluster3D") << "** Parallel lines, need a solution here" << std::endl;
        arcLenAxis = 0.;
        arcLenWire = 0.;
    }

    // Now get the hit position we'll use for the pca analysis
    poca =      Eigen::Vector3f(wirePos[0] + arcLenWire * wireDir(0),
                                wirePos[1] + arcLenWire * wireDir(1),
                                wirePos[2] + arcLenWire * wireDir(2));
    Eigen::Vector3f axisPocaPos(axisPos[0] + arcLenAxis  * axisDir(0),
                                axisPos[1] + arcLenAxis  * axisDir(1),
                                axisPos[2] + arcLenAxis  * axisDir(2));

    float deltaX(poca(0) - axisPocaPos(0));
    float deltaY(poca(1) - axisPocaPos(1));
    float deltaZ(poca(2) - axisPocaPos(2));
    float doca2(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);

    doca = sqrt(doca2);

    return;
}


struct Sort3DHitsByDocaToAxis
{
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        return left->getDocaToAxis() < right->getDocaToAxis();
    }

};

struct Sort3DHitsByArcLen3D
{
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        return left->getArclenToPoca() < right->getArclenToPoca();
    }

};

struct Sort3DHitsByAbsArcLen3D
{
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        return fabs(left->getArclenToPoca()) < fabs(right->getArclenToPoca());
    }

};

void PrincipalComponentsAlg::PCAAnalysis(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, float doca3DScl) const
{
    // This is the controlling outside function for running
    // a Principal Components Analysis on the hits in our
    // input cluster.
    // There is a bootstrap process to be followed
    // 1) Get the initial results from all 3D hits associated
    //    with the cluster
    // 2) Refine the axis by using only the 2D hits on wires

    // Get the first axis from all 3D hits
    PCAAnalysis_3D(hitPairVector, pca);

    // Make sure first pass was good
    if (!pca.getSvdOK()) return;

    // First attempt to refine it using only 2D information
    reco::PrincipalComponents pcaLoop = pca;

    PCAAnalysis_2D(hitPairVector, pcaLoop);

    // If valid result then go to next steps
    if (pcaLoop.getSvdOK())
    {
        // Let's check the angle between the original and the updated axis
        float cosAngle = pcaLoop.getEigenVectors().row(0) * pca.getEigenVectors().row(0).transpose();

        // Set the scale factor for the outlier rejection
        float sclFctr(3.);

        // If we had a significant change then let's do some outlier rejection, etc.
        if (cosAngle < 1.0)   // pretty much everyone takes a turn
        {
            int   maxIterations(3);
            float maxRange = 3.*sqrt(pcaLoop.getEigenValues()[1]);
            float aveDoca  = pcaLoop.getAveHitDoca();                 // was 0.2

            maxRange = sclFctr * 0.5*(maxRange+aveDoca); // was std::max(maxRange, aveDoca);

            int numRejHits = PCAAnalysis_reject2DOutliers(hitPairVector, pcaLoop, maxRange);
            int totalRejects(numRejHits);
            int maxRejects(0.4*hitPairVector.size());

            // Try looping to see if we improve things
            while(maxIterations-- && numRejHits > 0 && totalRejects < maxRejects)
            {
                // Run the PCA
                PCAAnalysis_2D(hitPairVector, pcaLoop, true);

                maxRange = sclFctr * 0.5*(3.*sqrt(pcaLoop.getEigenValues()[1])+pcaLoop.getAveHitDoca());

                numRejHits = PCAAnalysis_reject2DOutliers(hitPairVector, pcaLoop, maxRange);
            }

        }

        // Ok at this stage copy the latest results back into the cluster
        pca = pcaLoop;

        // Now we make one last pass through the 3D hits to reject outliers there
        PCAAnalysis_reject3DOutliers(hitPairVector, pca, doca3DScl * pca.getAveHitDoca());
    }
    else pca = pcaLoop;

    return;
}

void PrincipalComponentsAlg::PCAAnalysis_3D(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, bool skeletonOnly) const
{
    // We want to run a PCA on the input TkrVecPoints...
    // The steps are:
    // 1) do a mean normalization of the input vec points
    // 2) compute the covariance matrix
    // 3) run the SVD
    // 4) extract the eigen vectors and values
    // see what happens

    // Run through the HitPairList and get the mean position of all the hits
    Eigen::Vector3d meanPos(Eigen::Vector3d::Zero());
    double          meanWeightSum(0.);
    int             numPairsInt(0);

//    const float minimumDeltaPeakSig(0.00001);
    double minimumDeltaPeakSig(0.00001);

    // Want to use the hit "chi square" to weight the hits but we need to put a lower limit on its value
    // to prevent a few hits being over counted.
    // This is a bit experimental until we can evaluate the cost (time to calculate) vs the benefit
    // (better fits)..
    std::vector<double> hitChiSquareVec;

    hitChiSquareVec.resize(hitPairVector.size());

    std::transform(hitPairVector.begin(),hitPairVector.end(),hitChiSquareVec.begin(),[](const auto& hit){return hit->getHitChiSquare();});
    std::sort(hitChiSquareVec.begin(),hitChiSquareVec.end());

    size_t numToKeep = 0.8 * hitChiSquareVec.size();

    hitChiSquareVec.resize(numToKeep);

    double aveValue = std::accumulate(hitChiSquareVec.begin(),hitChiSquareVec.end(),double(0.)) / double(hitChiSquareVec.size());
    double rms      = std::sqrt(std::inner_product(hitChiSquareVec.begin(),hitChiSquareVec.end(), hitChiSquareVec.begin(), 0.,std::plus<>(),[aveValue](const auto& left,const auto& right){return (left - aveValue) * (right - aveValue);}) / double(hitChiSquareVec.size()));

    minimumDeltaPeakSig = std::max(minimumDeltaPeakSig, aveValue - rms);

//    std::cout << "===>> Calculating PCA, ave chiSquare: " << aveValue << ", rms: " << rms << ", cut: " << minimumDeltaPeakSig << std::endl;

    for (const auto& hit : hitPairVector)
    {
        if (skeletonOnly && !((hit->getStatusBits() & reco::ClusterHit3D::SKELETONHIT) == reco::ClusterHit3D::SKELETONHIT)) continue;

        // Weight the hit by the peak time difference significance
        double weight = std::max(minimumDeltaPeakSig, double(hit->getHitChiSquare())); //hit->getDeltaPeakTime()); ///hit->getSigmaPeakTime());

        meanPos(0) += hit->getPosition()[0] * weight;
        meanPos(1) += hit->getPosition()[1] * weight;
        meanPos(2) += hit->getPosition()[2] * weight;
        numPairsInt++;

        meanWeightSum += weight;
    }

    meanPos /= meanWeightSum;

    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double xizi(0.0);
    double yi2(0.0);
    double yizi(0.0);
    double zi2(0.);
    double weightSum(0.);

    // Back through the hits to build the matrix
    for (const auto& hit : hitPairVector)
    {
        if (skeletonOnly && !((hit->getStatusBits() & reco::ClusterHit3D::SKELETONHIT) == reco::ClusterHit3D::SKELETONHIT)) continue;

        double weight = 1. / std::max(minimumDeltaPeakSig, double(hit->getHitChiSquare())); //hit->getDeltaPeakTime()); ///hit->getSigmaPeakTime());
        double x      = (hit->getPosition()[0] - meanPos(0)) * weight;
        double y      = (hit->getPosition()[1] - meanPos(1)) * weight;
        double z      = (hit->getPosition()[2] - meanPos(2)) * weight;

        weightSum += weight*weight;

        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z * z;
    }

    // Using Eigen package
    Eigen::Matrix3d sig;

    sig <<  xi2, xiyi, xizi,
           xiyi,  yi2, yizi,
           xizi, yizi,  zi2;

    sig *= 1./weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenMat(sig);

    if (eigenMat.info() == Eigen::ComputationInfo::Success)
    {
        // Now copy output
        // The returned eigen values and vectors will be returned in an xyz system where x is the smallest spread,
        // y is the next smallest and z is the largest. Adopt that convention going forward
        reco::PrincipalComponents::EigenValues  recobEigenVals = eigenMat.eigenvalues().cast<float>();
        reco::PrincipalComponents::EigenVectors recobEigenVecs = eigenMat.eigenvectors().transpose().cast<float>();

        // Check for a special case (which may have gone away with switch back to doubles for computation?)
        if (std::isnan(recobEigenVals[0]))
        {
            std::cout << "==> Third eigenvalue returns a nan" << std::endl;

            recobEigenVals[0] = 0.;

            // Assume the third axis is also kaput?
            recobEigenVecs.row(0) = recobEigenVecs.row(1).cross(recobEigenVecs.row(2));
        }

        // Store away
        pca = reco::PrincipalComponents(true, numPairsInt, recobEigenVals, recobEigenVecs, meanPos.cast<float>());
    }
    else
    {
        mf::LogDebug("Cluster3D") << "PCA decompose failure, numPairs = " << numPairsInt << std::endl;
        pca = reco::PrincipalComponents();
    }

    return;
}

void PrincipalComponentsAlg::PCAAnalysis_2D(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, bool updateAvePos) const
{
    // Once an axis has been found our goal is to refine it by using only the 2D hits
    // We'll get 3D information for each of these by using the axis as a reference and use
    // the point of closest approach as the 3D position

    // Define elements of our covariance matrix
    float xi2(0.);
    float xiyi(0.);
    float xizi(0.);
    float yi2(0.);
    float yizi(0.);
    float zi2(0.);

    float           aveHitDoca(0.);
    Eigen::Vector3f avePosUpdate(0.,0.,0.);
    int             nHits(0);

    // Recover existing line parameters for current cluster
    const reco::PrincipalComponents& inputPca = pca;
    Eigen::Vector3f                  avePosition(inputPca.getAvePosition());
    Eigen::Vector3f                  axisDirVec(inputPca.getEigenVectors().row(0));

    // We float loop here so we can use this method for both the first time through
    // and a second time through where we re-calculate the mean position
    // So, we need to keep track of the poca which we do with a float vector
    std::vector<Eigen::Vector3f> hitPosVec;

    // Outer loop over 3D hits
    for (const auto& hit3D : hitPairVector)
    {
        // Inner loop over 2D hits
        for (const auto& hit : hit3D->getHits())
        {
            // Step one is to set up to determine the point of closest approach of this 2D hit to
            // the cluster's current axis.
            // Get this wire's geometry object
            const geo::WireID&  hitID     = hit->WireID();
            const geo::WireGeo& wire_geom = m_geometry->WireIDToWireGeo(hitID);

            // From this, get the parameters of the line for the wire
            double wirePosArr[3] = {0.,0.,0.};
            wire_geom.GetCenter(wirePosArr);

            Eigen::Vector3f wireCenter(wirePosArr[0], wirePosArr[1], wirePosArr[2]);
            Eigen::Vector3f wireDirVec(wire_geom.Direction().X(),wire_geom.Direction().Y(),wire_geom.Direction().Z());

            // Correct the wire position in x to set to correspond to the drift time
            float hitPeak(hit->getHit()->PeakTime());

            Eigen::Vector3f wirePos(m_detector->ConvertTicksToX(hitPeak, hitID.Plane, hitID.TPC, hitID.Cryostat), wireCenter[1], wireCenter[2]);

            // Compute the wire plane normal for this view
            Eigen::Vector3f xAxis(1.,0.,0.);
            Eigen::Vector3f planeNormal = xAxis.cross(wireDirVec);   // This gives a normal vector in +z for a Y wire

            float docaInPlane(wirePos[0] - avePosition[0]);
            float arcLenToPlane(0.);
            float cosAxisToPlaneNormal = axisDirVec.dot(planeNormal);

            Eigen::Vector3f axisPlaneIntersection = wirePos;
            Eigen::Vector3f hitPosTVec            = wirePos;

            if (fabs(cosAxisToPlaneNormal) > 0.)
            {
                Eigen::Vector3f deltaPos = wirePos - avePosition;

                arcLenToPlane         = deltaPos.dot(planeNormal) / cosAxisToPlaneNormal;
                axisPlaneIntersection = avePosition + arcLenToPlane * axisDirVec;
                docaInPlane           = wirePos[0] - axisPlaneIntersection[0];

                Eigen::Vector3f axisToInter  = axisPlaneIntersection - wirePos;
                float           arcLenToDoca = axisToInter.dot(wireDirVec);

                hitPosTVec += arcLenToDoca * wireDirVec;
            }

            // Get a vector from the wire position to our cluster's current average position
            Eigen::Vector3f wVec = avePosition - wirePos;

            // Get the products we need to compute the arc lengths to the distance of closest approach
            float a(axisDirVec.dot(axisDirVec));
            float b(axisDirVec.dot(wireDirVec));
            float c(wireDirVec.dot(wireDirVec));
            float d(axisDirVec.dot(wVec));
            float e(wireDirVec.dot(wVec));

            float den(a*c - b*b);
            float arcLen1(0.);
            float arcLen2(0.);

            // Parallel lines is a special case
            if (fabs(den) > m_parallel)
            {
                arcLen1 = (b*e - c*d) / den;
                arcLen2 = (a*e - b*d) / den;
            }
            else
            {
                mf::LogDebug("Cluster3D") << "** Parallel lines, need a solution here" << std::endl;
                break;
            }

            // Now get the hit position we'll use for the pca analysis
            //float hitPos[]  = {wirePos[0]     + arcLen2 * wireDirVec[0],
            //                    wirePos[1]     + arcLen2 * wireDirVec[1],
            //                    wirePos[2]     + arcLen2 * wireDirVec[2]};
            //float axisPos[] = {avePosition[0] + arcLen1 * axisDirVec[0],
            //                    avePosition[1] + arcLen1 * axisDirVec[1],
            //                    avePosition[2] + arcLen1 * axisDirVec[2]};
            Eigen::Vector3f hitPos  = wirePos + arcLen2 * wireDirVec;
            Eigen::Vector3f axisPos = avePosition + arcLen1 * axisDirVec;
            float           deltaX  = hitPos(0) - axisPos(0);
            float           deltaY  = hitPos(1) - axisPos(1);
            float           deltaZ  = hitPos(2) - axisPos(2);
            float           doca2   = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            float           doca    = sqrt(doca2);

            docaInPlane = doca;

            aveHitDoca += fabs(docaInPlane);

            //Eigen::Vector3f deltaPos  = hitPos - hitPosTVec;
            //float   deltaDoca = doca - docaInPlane;

            //if (fabs(deltaPos[0]) > 1. || fabs(deltaPos[1]) > 1. || fabs(deltaPos[2]) > 1. || fabs(deltaDoca) > 2.)
            //{
            //    std::cout << "**************************************************************************************" << std::endl;
            //    std::cout << "Diff in doca: " << deltaPos[0] << "," << deltaPos[1] << "," << deltaPos[2] << ", deltaDoca: " << deltaDoca << std::endl;
            //    std::cout << "-- HitPosTVec: " << hitPosTVec[0] << "," << hitPosTVec[1] << "," << hitPosTVec[2] << std::endl;
            //    std::cout << "-- hitPos:     " << hitPos[0] << "," << hitPos[1] << "," << hitPos[2] << std::endl;
            //    std::cout << "-- WirePos:    " << wirePos[0] << "," << wirePos[1] << "," << wirePos[2] << std::endl;
            //}

            // Set the hit's doca and arclen
            hit->setDocaToAxis(fabs(docaInPlane));
            hit->setArcLenToPoca(arcLenToPlane);

            // If this point is considered an outlier then we skip
            // the accumulation for average position and covariance
            if (hit->getStatusBits() & 0x80)
            {
                continue;
            }

            //hitPosTVec = hitPos;

            avePosUpdate += hitPosTVec;

            hitPosVec.push_back(hitPosTVec);

            nHits++;
        }
    }

    // Get updated average position
    avePosUpdate /= float(nHits);

    // Get the average hit doca
    aveHitDoca /= float(nHits);

    if (updateAvePos)
    {
        avePosition = avePosUpdate;
    }

    // Now loop through the hits and build out the covariance matrix
    for(auto& hitPos : hitPosVec)
    {
        // And increment the values in the covariance matrix
        float x = hitPos[0] - avePosition[0];
        float y = hitPos[1] - avePosition[1];
        float z = hitPos[2] - avePosition[2];

        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z * z;
    }

    // Accumulation done, now do the actual work

    // Using Eigen package
    Eigen::Matrix3f sig;

    sig <<   xi2, xiyi, xizi,
            xiyi,  yi2, yizi,
            xizi, yizi,  zi2;

    sig *= 1./(nHits - 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

    if (eigenMat.info() == Eigen::ComputationInfo::Success)
    {
       // Now copy output
        // The returned eigen values and vectors will be returned in an xyz system where x is the smallest spread,
        // y is the next smallest and z is the largest. Adopt that convention going forward
        reco::PrincipalComponents::EigenValues  recobEigenVals = eigenMat.eigenvalues().cast<float>();
        reco::PrincipalComponents::EigenVectors recobEigenVecs = eigenMat.eigenvectors().transpose().cast<float>();

        // Store away
        pca = reco::PrincipalComponents(true, nHits, recobEigenVals, recobEigenVecs, avePosition, aveHitDoca);
    }
    else
    {
        mf::LogDebug("Cluster3D") << "PCA decompose failure, numPairs = " << nHits << std::endl;
        pca = reco::PrincipalComponents();
    }

    return;
}

void PrincipalComponentsAlg::PCAAnalysis_calc3DDocas(const reco::HitPairListPtr&      hitPairVector,
                                                     const reco::PrincipalComponents& pca) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.

    // We'll need the current PCA axis to determine doca and arclen
    Eigen::Vector3f avePosition(pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);
    Eigen::Vector3f axisDirVec(pca.getEigenVectors().row(2));

    // We want to keep track of the average
    float aveDoca3D(0.);

    // Outer loop over views
    for (const auto* clusterHit3D : hitPairVector)
    {
        // Always reset the existing status bit
        clusterHit3D->clearStatusBits(0x80);

        // Now we need to calculate the doca and poca...
        // Start by getting this hits position
        Eigen::Vector3f clusPos(clusterHit3D->getPosition()[0],clusterHit3D->getPosition()[1],clusterHit3D->getPosition()[2]);

        // Form a TVector from this to the cluster average position
        Eigen::Vector3f clusToHitVec = clusPos - avePosition;

        // With this we can get the arclength to the doca point
        float arclenToPoca = clusToHitVec.dot(axisDirVec);

        // Get the coordinates along the axis for this point
        Eigen::Vector3f docaPos = avePosition + arclenToPoca * axisDirVec;

        // Now get doca and poca
        Eigen::Vector3f docaPosToClusPos = clusPos - docaPos;
        float           docaToAxis       = docaPosToClusPos.norm();

        aveDoca3D += docaToAxis;

        // Ok, set the values in the hit
        clusterHit3D->setDocaToAxis(docaToAxis);
        clusterHit3D->setArclenToPoca(arclenToPoca);
    }

    // Compute the average and store
    aveDoca3D /= float(hitPairVector.size());

    pca.setAveHitDoca(aveDoca3D);

    return;
}

void PrincipalComponentsAlg::PCAAnalysis_calc2DDocas(const reco::Hit2DListPtr&        hit2DListPtr,
                                                     const reco::PrincipalComponents& pca         ) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.

    // We'll need the current PCA axis to determine doca and arclen
    Eigen::Vector3f avePosition(pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);
    Eigen::Vector3f axisDirVec(pca.getEigenVectors().row(2));

    // Recover the principle eigen value for range constraints
    float maxArcLen = 4.*sqrt(pca.getEigenValues()[0]);

    // We want to keep track of the average
    float aveHitDoca(0.);

    // Outer loop over views
    for (const auto* hit : hit2DListPtr)
    {
        // Step one is to set up to determine the point of closest approach of this 2D hit to
        // the cluster's current axis. We do that by finding the point of intersection of the
        // cluster's axis with a plane defined by the wire the hit is associated with.
        // Get this wire's geometry object
        const geo::WireID&  hitID     = hit->WireID();
        const geo::WireGeo& wire_geom = m_geometry->WireIDToWireGeo(hitID);

        // From this, get the parameters of the line for the wire
        double wirePosArr[3] = {0.,0.,0.};
        wire_geom.GetCenter(wirePosArr);

        Eigen::Vector3f wireCenter(wirePosArr[0], wirePosArr[1], wirePosArr[2]);
        Eigen::Vector3f wireDirVec(wire_geom.Direction().X(),wire_geom.Direction().Y(),wire_geom.Direction().Z());

        // Correct the wire position in x to set to correspond to the drift time
        Eigen::Vector3f wirePos(hit->getXPosition(), wireCenter[1], wireCenter[2]);

        // Compute the wire plane normal for this view
        Eigen::Vector3f xAxis(1.,0.,0.);
        Eigen::Vector3f planeNormal = xAxis.cross(wireDirVec);   // This gives a normal vector in +z for a Y wire

        float arcLenToPlane(0.);
        float docaInPlane(wirePos[0] - avePosition[0]);
        float cosAxisToPlaneNormal = axisDirVec.dot(planeNormal);

        Eigen::Vector3f axisPlaneIntersection = wirePos;

        // If current cluster axis is not parallel to wire plane then find intersection point
        if (fabs(cosAxisToPlaneNormal) > 0.)
        {
            Eigen::Vector3f deltaPos = wirePos - avePosition;

            arcLenToPlane         = std::min(float(deltaPos.dot(planeNormal) / cosAxisToPlaneNormal), maxArcLen);
            axisPlaneIntersection = avePosition + arcLenToPlane * axisDirVec;

            Eigen::Vector3f axisToInter  = axisPlaneIntersection - wirePos;
            float           arcLenToDoca = axisToInter.dot(wireDirVec);

            // If the arc length along the wire to the poca is outside the TPC then reset
            if (fabs(arcLenToDoca) > wire_geom.HalfL()) arcLenToDoca = wire_geom.HalfL();

            // If we were successful in getting to the wire plane then the doca is simply the
            // difference in x coordinates... but we hvae to worry about the special cases so
            // we calculate a 3D doca based on arclengths above...
            Eigen::Vector3f docaVec = axisPlaneIntersection - (wirePos + arcLenToDoca * wireDirVec);
            docaInPlane = docaVec.norm();
        }

        aveHitDoca += fabs(docaInPlane);

        // Set the hit's doca and arclen
        hit->setDocaToAxis(fabs(docaInPlane));
        hit->setArcLenToPoca(arcLenToPlane);
    }

    // Compute the average and store
    aveHitDoca /= float(hit2DListPtr.size());

    pca.setAveHitDoca(aveHitDoca);

    return;
}

int PrincipalComponentsAlg::PCAAnalysis_reject2DOutliers(const reco::HitPairListPtr& hitPairVector,
                                                         reco::PrincipalComponents&  pca,
                                                         float                       maxDocaAllowed) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.

    // First get the average doca scaled by some appropriate factor
    int    numRejHits(0);

    // Outer loop over views
    for (const auto& hit3D : hitPairVector)
    {
        // Inner loop over hits in this view
        for (const auto& hit : hit3D->getHits())
        {
            // Always reset the existing status bit
            hit->clearStatusBits(0x80);

            if (hit->getDocaToAxis() > maxDocaAllowed)
            {
                hit->setStatusBit(0x80);
                numRejHits++;
            }
        }
    }

    return numRejHits;
}

int PrincipalComponentsAlg::PCAAnalysis_reject3DOutliers(const reco::HitPairListPtr&      hitPairVector,
                                                         const reco::PrincipalComponents& pca,
                                                         float                            maxDocaAllowed) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.

    // First get the average doca scaled by some appropriate factor
    int    numRejHits(0);

    // We'll need the current PCA axis to determine doca and arclen
    Eigen::Vector3f avePosition(pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);
    Eigen::Vector3f axisDirVec(pca.getEigenVectors().row(2));

    // Outer loop over views
    for (const auto* clusterHit3D : hitPairVector)
    {
        // Always reset the existing status bit
        clusterHit3D->clearStatusBits(0x80);

        // Now we need to calculate the doca and poca...
        // Start by getting this hits position
        Eigen::Vector3f clusPos(clusterHit3D->getPosition()[0],clusterHit3D->getPosition()[1],clusterHit3D->getPosition()[2]);

        // Form a TVector from this to the cluster average position
        Eigen::Vector3f clusToHitVec = clusPos - avePosition;

        // With this we can get the arclength to the doca point
        float arclenToPoca = clusToHitVec.dot(axisDirVec);

        // Get the coordinates along the axis for this point
        Eigen::Vector3f docaPos = avePosition + arclenToPoca * axisDirVec;

        // Now get doca and poca
        Eigen::Vector3f docaPosToClusPos = clusPos - docaPos;
        float   docaToAxis       = docaPosToClusPos.norm();

        // Ok, set the values in the hit
        clusterHit3D->setDocaToAxis(docaToAxis);
        clusterHit3D->setArclenToPoca(arclenToPoca);

        // Check to see if this is a keeper
        if (clusterHit3D->getDocaToAxis() > maxDocaAllowed)
        {
            clusterHit3D->setStatusBit(0x80);
            numRejHits++;
        }
    }

    return numRejHits;
}



} // namespace lar_cluster3d
