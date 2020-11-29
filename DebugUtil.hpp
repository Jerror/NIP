#ifndef INC_DebugUtil_hpp
#define INC_DebugUtil_hpp

#ifdef env_MPI
#include <mpi.h>
#define ABORT(code) MPI_Abort(MPI_COMM_WORLD, code)
#else
#include <cstdlib>
#define ABORT(code) abort()
#endif

#ifndef __GNUG__
#define __PRETTY_FUNCTION__ __func__
#endif

// formatting follows https://www.gnu.org/prep/standards/html_node/Errors.html
#define HERE \
    __FILE__<<":"<<__LINE__<<": "

#define ERROR2(what, code) \
    std::cout<<HERE<<" Error in call to "<<__PRETTY_FUNCTION__<<std::endl<<what<<std::endl; \
    ABORT(code);

// The utility of the error code is made mostly redundant by HERE, so we're
// usually happy with an arbitrary default value.
#define ERROR(what) ERROR2(what, 9999)

#define WARNING(what) \
    std::cout<<HERE<<" Warning in call to "<<__PRETTY_FUNCTION__<<std::endl<<what<<std::endl;

#define DIAGNOSTIC(what) \
    std::cout<<what<<std::endl;

#define DIAGNOSTIC_V(what) \
    std::cout<<HERE<<what<<std::endl;

#define DIAGNOSTIC_VV(what) \
    std::cout<<HERE<<__PRETTY_FUNCTION__<<": "<<what<<std::endl;

#endif // INC_DebugUtil_hpp
