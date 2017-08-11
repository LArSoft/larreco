///////////////////////////////////////////////////////////////////////
///
/// \file   IWaveformAlgs.h
///
/// \brief  This provides an art tool interface to a package of
///         algorithms that provide handy functions on waveforms
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IWaveformAlgs_H
#define IWaveformAlgs_H

#include "fhiclcpp/ParameterSet.h"

namespace reco_tool
{
    class IWaveformAlgs
    {
    public:
        virtual ~IWaveformAlgs() noexcept = default;
        
        // Define standard art tool interface
        virtual void configure(const fhicl::ParameterSet& pset) = 0;
        
        // Rebin the input vector
        virtual void reBin(const std::vector<float>&, std::vector<float>&, size_t) const = 0;
        
        // Do a running average of the input vector
        virtual void doBinAverage(const std::vector<float>&, std::vector<float>&, size_t) const = 0;
        
        // Returned smoothed version of input vector ("triangle smooth")
        virtual void getSmoothVec(const std::vector<float>&,std::vector<float>&) const = 0;
        
        // Return the derivative of the input waveform
        virtual void getDerivativeVec(const std::vector<float>&,std::vector<float>&) const = 0;
        
        // Return the smoothed derivative of the input waveform
        virtual void getSmoothDerivativeVec(const std::vector<float>&,std::vector<float>&) const = 0;
    };
}

#endif
