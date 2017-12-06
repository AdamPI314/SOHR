#ifndef __GLOBAL_MACROS_H_
#define __GLOBAL_MACROS_H_

#ifdef _WIN64
    #define __WINDOWS_
	#define __NO_USE_CANTERA_
	#define __NO_USE_MPI_
#elif _WIN32
    #define __WINDOWS_
	#define __NO_USE_CANTERA_
	#define __NO_USE_MPI_

#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_OS_IPHONE && TARGET_IPHONE_SIMULATOR
        // define something for simulator   
    #elif TARGET_OS_IPHONE
        // define something for iphone  
    #else
        #define TARGET_OS_OSX 1
        // define something for OSX
    #endif

#elif __linux
    #define __LINUX_
	#define __CANTERA_AVAILABLE_
	#define __USE_CANTERA_

	#define __CHEMKIN_AVAILABLE_
	#define __LSODE_AVAILABLE_
	#define __USE_MPI_
#elif __unix 
    // all unices not caught above
#elif __posix
    // POSIX
#endif


#endif
