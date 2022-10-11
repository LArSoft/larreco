/**
 * \file CFAlgoTimeOverlap.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CFAlgoTimeOverlap
 *
 * @author ariana hackenburg
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CFALGOTIMEOVERLAP_H
#define RECOTOOL_CFALGOTIMEOVERLAP_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CFloatAlgoBase.h"

namespace cmtool {
  /**
     \class CFAlgoTimeOverlap
     User implementation for CFloatAlgoBase class
     doxygen documentation!
  */
  class CFAlgoTimeOverlap : public CFloatAlgoBase {

  public:
    /// Default constructor
    CFAlgoTimeOverlap();

    /**This algorithm calculates the difference between start and end times for merged clusters,
                and compares across planes to form matches.
    */
    float Float(util::GeometryUtilities const&,
                const std::vector<const cluster::ClusterParamsAlg*>& clusters) override;

    void SetStartTimeCut(float start_time) { _start_time_cut = start_time; }

    void SetRatioCut(float ratio) { _time_ratio_cut = ratio; }

    //Order the theta, phi, hits per plane to make cuts convenient
    /*
    virtual void SetMaxMiddleMin(const double first, const double second, const double third,
                                 double &most, double &middle, double &least) ;
    */
    void SetDebug(bool debug) { _debug = debug; }

    void SetVerbose(bool verbose) override { _verbose = verbose; }

    void RequireThreePlanes(bool doit) { _require_3planes = doit; }

    void Report() override;

    void Reset() override;

  protected:
    float _time_ratio_cut;
    float _start_time_cut;
    bool _debug;
    bool _verbose;
    bool _require_3planes;
  };
}
#endif
/** @} */ // end of doxygen group
