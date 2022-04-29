/**
 * \file CMatchManager.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CMatchManager
 *
 * @author kazuhiro
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CMATCHMANAGER_H
#define RECOTOOL_CMATCHMANAGER_H

#include <stddef.h>

#include "larreco/RecoAlg/CMTool/CMToolBase/CMManagerBase.h"
#include "larreco/RecoAlg/CMTool/CMToolBase/CMatchBookKeeper.h"

namespace cmtool {

  class CFloatAlgoBase;
  /**
     \class CMatchManager
     A class that instantiates merging algorithm(s) and run.
     The book-keeping of merged cluster sets are done by CMatchBookKeeper.
  */
  class CMatchManager : public CMManagerBase {

  private:
    /// Default constructor is private because we need an argument to configure w/ # planes in the detector
    CMatchManager();

  public:
    CMatchManager(size_t nplanes);

    /// Default destructor
    virtual ~CMatchManager() {}

    /// Method to reset itself
    virtual void Reset();

    /// A simple method to add an algorithm for merging
    void
    AddMatchAlgo(CFloatAlgoBase* algo)
    {
      _match_algo = algo;
    }

    /// A method to obtain book keeper
    const CMatchBookKeeper&
    GetBookKeeper() const
    {
      return _book_keeper;
    }

  protected:
    //
    // FMWK functions override
    //

    /// FMWK function called @ beginning of Process()
    virtual void EventBegin();

    /// FMWK function called @ beginning of iterative loop inside Process()
    virtual void IterationBegin();

    /// FMWK function called @ iterative loop inside Process()
    virtual bool IterationProcess(util::GeometryUtilities const& gser);

    /// FMWK function called @ end of iterative loop inside Process()
    virtual void IterationEnd();

    /// FMWK function called @ end of Process()
    virtual void EventEnd();

  protected:
    /// Book keeper instance
    CMatchBookKeeper _book_keeper;

    /// Merging algorithm
    ::cmtool::CFloatAlgoBase* _match_algo;

    /// Number of planes
    size_t _nplanes;
  };
}

#endif
/** @} */ // end of doxygen group
