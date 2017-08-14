////////////////////////////////////////////////////////////////////////
/// \file   WaveformAlgs.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/IWaveformAlgs.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include <cmath>
#include <fstream>

namespace reco_tool
{

class WaveformAlgs : IWaveformAlgs
{
public:
    explicit WaveformAlgs(const fhicl::ParameterSet& pset);
    
    ~WaveformAlgs();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    // Rebin the input vector
    void reBin(const std::vector<float>&, std::vector<float>&, size_t) const override;
    
    // Do a running average of the input vector
    void doBinAverage(const std::vector<float>&, std::vector<float>&, size_t) const override;
    
    // Returned smoothed version of input vector ("triangle smooth")
    void getSmoothVec(const std::vector<float>&,std::vector<float>&) const override;
    
    // Return the derivative of the input waveform
    void getDerivativeVec(const std::vector<float>&,std::vector<float>&) const override;
    
    // Return the smoothed derivative of the input waveform
    void getSmoothDerivativeVec(const std::vector<float>&,std::vector<float>&) const override;
    
private:
};
    
//----------------------------------------------------------------------
// Constructor.
WaveformAlgs::WaveformAlgs(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
WaveformAlgs::~WaveformAlgs()
{
}
    
void WaveformAlgs::configure(const fhicl::ParameterSet& pset)
{
    return;
}
    
void WaveformAlgs::doBinAverage(const std::vector<float>& inputVec,
                                      std::vector<float>& outputVec,
                                      size_t              binsToAverage) const
{
    size_t halfBinsToAverage(binsToAverage/2);
    
    float runningSum(0.);
    
    for(size_t idx = 0; idx < halfBinsToAverage; idx++) runningSum += inputVec[idx];
    
    outputVec.resize(inputVec.size());
    std::vector<float>::iterator outputVecItr = outputVec.begin();
    
    // First pass through to build the erosion vector
    for(std::vector<float>::const_iterator inputItr = inputVec.begin(); inputItr != inputVec.end(); inputItr++)
    {
        size_t startOffset = std::distance(inputVec.begin(),inputItr);
        size_t stopOffset  = std::distance(inputItr,inputVec.end());
        size_t count       = std::min(2 * halfBinsToAverage, std::min(startOffset + halfBinsToAverage + 1, halfBinsToAverage + stopOffset - 1));
        
        if (startOffset >= halfBinsToAverage) runningSum -= *(inputItr - halfBinsToAverage);
        if (stopOffset  >  halfBinsToAverage) runningSum += *(inputItr + halfBinsToAverage);
        
        *outputVecItr++ = runningSum / float(count);
    }
    
    return;
}
    
void WaveformAlgs::reBin(const std::vector<float>& inputVec,
                               std::vector<float>& outputVec,
                               size_t              nBinsToCombine) const
{
    size_t nNewBins = inputVec.size() / nBinsToCombine;
    
    if (inputVec.size() % nBinsToCombine > 0) nNewBins++;
    
    outputVec.resize(nNewBins, 0.);
    
    size_t outputBin = 0;
    
    for(size_t inputIdx = 0; inputIdx < inputVec.size();)
    {
        outputVec[outputBin] += inputVec[inputIdx++];
        
        if (inputIdx % nBinsToCombine == 0) outputBin++;
        
        if (outputBin > outputVec.size())
        {
            std::cout << "***** DISASTER!!! ****** outputBin: " << outputBin << ", inputIdx = " << inputIdx << std::endl;
            break;
        }
    }
    
    return;
}
    
void WaveformAlgs::getSmoothVec(const std::vector<float>& roiSignalVec, std::vector<float>& smoothVec) const
{
    // We don't want the smoothing procedure to factor into the output
    // So start by making a local copy of the input vector
    std::vector<float> tempVec = roiSignalVec;
    
    // Now run the "triangle" smoothing operation
    for(size_t idx = 2; idx < tempVec.size() - 2; idx++)
        smoothVec.at(idx) = (tempVec.at(idx-2) + 2.*tempVec.at(idx-1) + 3.*tempVec.at(idx) + 2.*tempVec.at(idx+1) + tempVec.at(idx+2))/9.;
    
    return;
}

void WaveformAlgs::getDerivativeVec(const std::vector<float>& roiSignalVec, std::vector<float>& derivativeVec) const
{
    derivativeVec.resize(roiSignalVec.size(),0.);
    
    for(size_t idx = 1; idx < roiSignalVec.size()-1; idx++)
        derivativeVec.at(idx) = 0.5 * (roiSignalVec.at(idx+1) - roiSignalVec.at(idx-1));
    
    return;
}
    
void WaveformAlgs::getSmoothDerivativeVec(const std::vector<float>& roiSignalVec, std::vector<float>& derivativeVec) const
{
    // First task is to compute the derivative of the input signal vector
    std::vector<float> tempVec;
    
    tempVec.resize(roiSignalVec.size(),0.);
    
    for(size_t idx = 1; idx < roiSignalVec.size()-1; idx++)
        tempVec.at(idx) = 0.5 * (roiSignalVec.at(idx+1) - roiSignalVec.at(idx-1));
    
    // Now smooth the derivative vector
    derivativeVec.resize(roiSignalVec.size(),0.);
    
    for(size_t idx = 2; idx < tempVec.size() - 2; idx++)
        derivativeVec.at(idx) = (tempVec.at(idx-2) + 2.*tempVec.at(idx-1) + 3.*tempVec.at(idx) + 2.*tempVec.at(idx+1) + tempVec.at(idx+2))/9.;
    
    // Simply copy the two unsmoothed elements to fill out vector
    if (derivativeVec.size() > 2)
    {
        derivativeVec.at(0)                      = tempVec.at(1);
        derivativeVec.at(1)                      = tempVec.at(1);
        derivativeVec.at(derivativeVec.size()-2) = tempVec.at(derivativeVec.size()-2);
        derivativeVec.at(derivativeVec.size()-1) = tempVec.at(derivativeVec.size()-2);
    }
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(WaveformAlgs)
}
