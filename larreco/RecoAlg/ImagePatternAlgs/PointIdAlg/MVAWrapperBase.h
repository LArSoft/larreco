//////////////////////////////////////////////////////////////////////////////
// \version 
//
// \brief Helper functions for MVAReader and MVAWriter wrappers
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAWRAPPERBASE_H
#define ANAB_MVAWRAPPERBASE_H

#include <typeinfo>
#include <string>

namespace anab {

/// Helper functions for MVAReader and MVAWriter wrappers.
class MVAWrapperBase {
public:

protected:

    std::string getProductName(std::type_info const & ti) const;

};

} // namespace anab
#endif //ANAB_MVAWRAPPERBASE

