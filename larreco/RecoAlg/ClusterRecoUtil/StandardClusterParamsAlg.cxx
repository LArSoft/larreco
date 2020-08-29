/** ****************************************************************************
 * @file   StandardClusterParamsAlg.cxx
 * @brief  Implementation of interface to class computing cluster parameters
 * @author petrillo@fnal.gov
 * @date   January 22, 2015
 * @see    StandardClusterParamsAlg.h
 *
 * ****************************************************************************/

// C/C++ standard library

// LArSoft libraries
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::DegreesToRadians()
#include "lardata/Utilities/PxHitConverter.h"

//==============================================================================
//===  cluster::StandardClusterParamsAlg

cluster::StandardClusterParamsAlg::StandardClusterParamsAlg()
{
  SetVerbose(0);
}

//------------------------------------------------------------------------------
void
cluster::StandardClusterParamsAlg::SetVerbose(int level /* = 1 */)
{
  ClusterParamsAlgBase::SetVerbose(level);
  algo.SetVerbose(level > 0);
}

//------------------------------------------------------------------------------
void
cluster::StandardClusterParamsAlg::Clear()
{
  algo.Initialize();
}

//------------------------------------------------------------------------------
void
cluster::StandardClusterParamsAlg::SetHits(util::GeometryUtilities const& gser,
                                           std::vector<recob::Hit const*> const& hits)
{
  Clear();
  util::PxHitConverter pxhitconverter{gser};
  algo.SetHits(pxhitconverter.ToPxHitVector(hits));
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::StartCharge(util::GeometryUtilities const& gser)
{
  if (NInputHits() == 0) return {0.F};
  return {(float)algo.StartCharge(gser)};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::EndCharge(util::GeometryUtilities const& gser)
{
  if (NInputHits() == 0) return {0.F};
  return {(float)algo.EndCharge(gser)};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::StartAngle()
{
  if (NInputHits() < 2) return {0.F};

  algo.GetRoughAxis();
  return {(float)util::DegreesToRadians(algo.GetParams().cluster_angle_2d)};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::EndAngle()
{
  return StartAngle(); // Ummm...this doesn't look right. FIXME
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::StartOpeningAngle()
{
  if (NInputHits() < 3) return {0.F};

  algo.RefineDirection();
  return {(float)algo.GetParams().opening_angle_charge_wgt};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::EndOpeningAngle()
{
  if (NInputHits() < 3) return {0.F};

  algo.RefineDirection();
  return {(float)algo.GetParams().closing_angle_charge_wgt};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::Integral()
{
  if (NInputHits() == 0) return {0.F};

  algo.GetAverages();
  return {(float)algo.GetParams().sum_charge};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::IntegralStdDev()
{
  if (NInputHits() < 2) return {0.F};

  algo.GetAverages();
  return {(float)algo.GetParams().rms_charge};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::SummedADC()
{
  if (NInputHits() == 0) return {0.F};

  // compute all the averages
  algo.GetAverages();
  double sumADC = algo.GetParams().sum_ADC;
  return {(float)sumADC, (float)std::sqrt(sumADC)};
}

//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::SummedADCStdDev()
{
  if (NInputHits() < 2) return {0.F};

  algo.GetAverages();
  return {(float)algo.GetParams().rms_ADC};
}

//------------------------------------------------------------------------------
size_t
cluster::StandardClusterParamsAlg::NHits()
{
  if (NInputHits() < 2) return NInputHits();

  algo.GetAverages();
  return (size_t)algo.GetParams().N_Hits;
}

//------------------------------------------------------------------------------
float
cluster::StandardClusterParamsAlg::MultipleHitDensity()
{
  if (NInputHits() < 2) return 0.0F;

  algo.GetAverages();
  return algo.GetParams().N_Wires ? algo.GetParams().multi_hit_wires / algo.GetParams().N_Wires :
                                    0.;
}

//------------------------------------------------------------------------------
float
cluster::StandardClusterParamsAlg::Width(util::GeometryUtilities const& gser)
{
  if (NInputHits() < 3) return 0.0F;

  algo.GetProfileInfo(gser);
  return algo.GetParams().width;
}

//------------------------------------------------------------------------------
size_t
cluster::StandardClusterParamsAlg::NInputHits() const
{
  return algo.GetNHits();
}
