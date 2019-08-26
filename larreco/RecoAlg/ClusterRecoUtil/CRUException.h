/**
 * \file CRUException.h
 *
 * \ingroup ClusterRecoUtil
 *
 * \brief Class def header for exception classes in ClusterRecoUtil package
 *
 * @author kazuhiro
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CRUEXCEPTION_H
#define CRUEXCEPTION_H

#include <string>
#include <exception>

namespace cluster {
  /**
     \class CRUException
     Generic (base) exception class
  */
  class CRUException : public std::exception{

  public:

    CRUException(std::string msg="") : _msg(msg)
    {}

    const char* what() const noexcept override
    {return _msg.c_str(); }

  private:

    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group
