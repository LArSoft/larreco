/** ****************************************************************************
 * @file   LazyClusterParamsAlg.cxx
 * @brief  Implementation of interface to class computing cluster parameters
 * @author petrillo@fnal.gov
 * @date   January 22, 2015
 * @see    LazyClusterParamsAlg.h
 *
 * ****************************************************************************/

// C/C++ standard library
#include <cmath>

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::DegreesToRadians()
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParams.h"
#include "larreco/RecoAlg/ClusterRecoUtil/LazyClusterParamsAlg.h"

//==============================================================================
//===  cluster::LazyClusterParamsAlg

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::StartCharge(util::GeometryUtilities const&)
{
  return {(float)params.start_charge};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::EndCharge(util::GeometryUtilities const&)
{
  return {(float)params.end_charge};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::StartAngle()
{
  return {(float)util::DegreesToRadians(params.cluster_angle_2d)};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::EndAngle()
{
  return StartAngle();
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::StartOpeningAngle()
{
  return {(float)params.opening_angle_charge_wgt};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::EndOpeningAngle()
{
  return {(float)params.closing_angle_charge_wgt};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::Integral()
{
  return {(float)params.sum_charge};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::IntegralStdDev()
{
  return {(float)params.rms_charge};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::SummedADC()
{
  const double sumADC = params.sum_ADC;
  return {(float)sumADC, (float)std::sqrt(sumADC)};
}

//------------------------------------------------------------------------------
cluster::LazyClusterParamsAlg::Measure_t
cluster::LazyClusterParamsAlg::SummedADCStdDev()
{
  return {(float)params.rms_ADC};
}

//------------------------------------------------------------------------------
size_t
cluster::LazyClusterParamsAlg::NHits()
{
  return (size_t)params.N_Hits;
}

//------------------------------------------------------------------------------
float
cluster::LazyClusterParamsAlg::MultipleHitDensity()
{
  return params.N_Wires ? params.multi_hit_wires / params.N_Wires : 0.;
}

//------------------------------------------------------------------------------
float
cluster::LazyClusterParamsAlg::Width(util::GeometryUtilities const&)
{
  return params.width;
}
