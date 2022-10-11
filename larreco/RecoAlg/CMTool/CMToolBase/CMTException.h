/**
 * \file CMTException.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for exception classes in CMTException
 *
 * @author kazuhiro
 */

/** \addtogroup CMTException

    @{*/
#ifndef RECOTOOL_CMTEXCEPTION_H
#define RECOTOOL_CMTEXCEPTION_H

#include <exception>
#include <string>

namespace cmtool {
  /**
     \class CMTException
     Generic (base) exception class
  */
  class CMTException : public std::exception {

  public:
    CMTException(std::string msg = "") : _msg(msg) {}

    const char* what() const noexcept override { return _msg.c_str(); }

  private:
    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group
