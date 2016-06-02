/** ****************************************************************************
 * @file   LazyClusterParamsAlg.cxx
 * @brief  Implementation of interface to class computing cluster parameters
 * @author petrillo@fnal.gov
 * @date   January 22, 2015
 * @see    LazyClusterParamsAlg.h
 * 
 * ****************************************************************************/

// C/C++ standard library
#include <vector>

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::DegreesToRadians()
#include "lardata/Utilities/PxHitConverter.h"
#include "larreco/RecoAlg/ClusterRecoUtil/LazyClusterParamsAlg.h"


//==============================================================================
//===  cluster::LazyClusterParamsAlg
//===  

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::StartCharge()
{
  return { (float) params.start_charge };
} // LazyClusterParamsAlg::StartCharge()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::EndCharge()
{
  return { (float) params.end_charge };
} // LazyClusterParamsAlg::EndCharge()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::StartAngle()
{
  return { (float) util::DegreesToRadians(params.cluster_angle_2d) };
} // LazyClusterParamsAlg::StartAngle()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::EndAngle()
{
  return StartAngle();
} // LazyClusterParamsAlg::EndAngle()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::StartOpeningAngle()
{
  return { (float) params.opening_angle_charge_wgt };
} // LazyClusterParamsAlg::StartOpeningAngle()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::EndOpeningAngle()
{
  return { (float) params.closing_angle_charge_wgt };
} // LazyClusterParamsAlg::EndOpeningAngle()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::Integral()
{
  return { (float) params.sum_charge };
} // LazyClusterParamsAlg::Integral()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::IntegralStdDev()
{
  return { (float) params.rms_charge };
} // LazyClusterParamsAlg::IntegralStdDev()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::SummedADC()
{
  const double sumADC = params.sum_ADC;
  return { (float) sumADC, (float) std::sqrt(sumADC) };
} // LazyClusterParamsAlg::SummedADC()


//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::SummedADCStdDev()
{
  return { (float) params.rms_ADC };
} // LazyClusterParamsAlg::SummedADCStdDev()


//------------------------------------------------------------------------------
size_t cluster::LazyClusterParamsAlg::NHits() {
  return (size_t) params.N_Hits;
} // LazyClusterParamsAlg::NHits()


//------------------------------------------------------------------------------
float cluster::LazyClusterParamsAlg::MultipleHitDensity() {
  return params.N_Wires? params.multi_hit_wires / params.N_Wires: 0.;
} // LazyClusterParamsAlg::MultipleHitDensity()
    

//------------------------------------------------------------------------------
float cluster::LazyClusterParamsAlg::Width() {
  return params.width;
} // LazyClusterParamsAlg::Width()


//------------------------------------------------------------------------------
