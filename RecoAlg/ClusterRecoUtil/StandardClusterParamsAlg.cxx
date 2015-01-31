/** ****************************************************************************
 * @file   StandardClusterParamsAlg.cxx
 * @brief  Implementation of interface to class computing cluster parameters
 * @author petrillo@fnal.gov
 * @date   January 22, 2015
 * @see    StandardClusterParamsAlg.h
 * 
 * ****************************************************************************/

// C/C++ standard library
#include <vector>

// LArSoft libraries
#include "SimpleTypesAndConstants/PhysicalConstants.h" // util::DegreesToRadians()
#include "Utilities/PxHitConverter.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"


//==============================================================================
//===  cluster::StandardClusterParamsAlg
//===  

cluster::StandardClusterParamsAlg::StandardClusterParamsAlg() {
  SetVerbose(0);
} // StandardClusterParamsAlg::StandardClusterParamsAlg()


//------------------------------------------------------------------------------
void cluster::StandardClusterParamsAlg::SetVerbose(int level /* = 1 */) {
  ClusterParamsAlgBase::SetVerbose(level);
  algo.SetVerbose(level > 0);
} // StandardClusterParamsAlg::SetVerbose()


//------------------------------------------------------------------------------
void cluster::StandardClusterParamsAlg::Clear() {
  algo.Initialize();
} // StandardClusterParamsAlg::Clear()


//------------------------------------------------------------------------------
void cluster::StandardClusterParamsAlg::SetHits
  (std::vector<recob::Hit const*> const& hits)
{
  Clear();
  util::PxHitConverter pxhitconverter;
  algo.SetHits(pxhitconverter.ToPxHitVector(hits));
} // StandardClusterParamsAlg::SetHits()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::StartCharge()
{
  return { (float) algo.StartCharge() };
} // StandardClusterParamsAlg::StartCharge()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::EndCharge()
{
  return { (float) algo.EndCharge() };
} // StandardClusterParamsAlg::EndCharge()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::StartAngle()
{
  // compute the rough direction and related quantities
  algo.GetRoughAxis();
  // return the relevant information, no uncertainty
  return { (float) util::DegreesToRadians(algo.GetParams().cluster_angle_2d) };
} // StandardClusterParamsAlg::StartAngle()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::EndAngle()
{
  return StartAngle();
} // StandardClusterParamsAlg::EndAngle()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::StartOpeningAngle()
{
  // compute the direction and related quantities
  algo.RefineDirection();
  // return the relevant information, no uncertainty
  return { (float) algo.GetParams().opening_angle_charge_wgt };
} // StandardClusterParamsAlg::StartOpeningAngle()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::EndOpeningAngle()
{
  // compute the direction and related quantities
  algo.RefineDirection();
  // return the relevant information, no uncertainty
  return { (float) algo.GetParams().closing_angle_charge_wgt };
} // StandardClusterParamsAlg::EndOpeningAngle()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::Integral()
{
  // compute all the averages
  algo.GetAverages();
  // return the relevant information, no uncertainty
  return { (float) algo.GetParams().sum_charge };
} // StandardClusterParamsAlg::Integral()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::IntegralStdDev()
{
  // compute all the averages
  algo.GetAverages();
  // return the relevant information, no uncertainty
  return { (float) algo.GetParams().rms_charge };
} // StandardClusterParamsAlg::IntegralStdDev()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::SummedADC()
{
  // compute all the averages
  algo.GetAverages();
  double sumADC = algo.GetParams().sum_ADC;
  // return the relevant information, no uncertainty
  return { (float) sumADC, (float) std::sqrt(sumADC) };
} // StandardClusterParamsAlg::SummedADC()


//------------------------------------------------------------------------------
cluster::StandardClusterParamsAlg::Measure_t
cluster::StandardClusterParamsAlg::SummedADCStdDev()
{
  // compute all the averages
  algo.GetAverages();
  // return the relevant information, no uncertainty
  return { (float) algo.GetParams().rms_ADC };
} // StandardClusterParamsAlg::SummedADCStdDev()


//------------------------------------------------------------------------------
size_t cluster::StandardClusterParamsAlg::NHits() {
  // compute all the averages
  algo.GetAverages();
  // return the relevant information, no uncertainty
  return { (size_t) algo.GetParams().N_Hits };
} // StandardClusterParamsAlg::NHits()


//------------------------------------------------------------------------------
float cluster::StandardClusterParamsAlg::MultipleHitWires() {
  // compute all the averages
  algo.GetAverages();
  // return the relevant information
  return algo.GetParams().N_Wires?
    algo.GetParams().multi_hit_wires / algo.GetParams().N_Wires: 0.;
} // StandardClusterParamsAlg::MultipleHitWires()
    

//------------------------------------------------------------------------------
float cluster::StandardClusterParamsAlg::Width() {
  // compute all the shower profile information
  algo.GetProfileInfo();
  // return the relevant information, no uncertainty
  return algo.GetParams().width;
} // StandardClusterParamsAlg::Width()


//------------------------------------------------------------------------------
