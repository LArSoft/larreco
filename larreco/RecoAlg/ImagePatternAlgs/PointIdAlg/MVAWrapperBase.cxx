//////////////////////////////////////////////////////////////////////////////
// \version 
//
// \brief Helper functions for MVAReader and MVAWriter wrappers
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/MVAWrapperBase.h"

#include <cxxabi.h>
#include <algorithm>

namespace anab {

std::string MVAWrapperBase::getProductName(std::type_info const & ti) const
{
    char* realname;
    int status;

    realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
    std::string pname(realname);
    free(realname);

    pname.erase( std::remove(pname.begin(), pname.end(), ' '), pname.end() );
    pname.erase( std::remove(pname.begin(), pname.end(), ':'), pname.end() );

    return pname;
}

} // namespace anab
