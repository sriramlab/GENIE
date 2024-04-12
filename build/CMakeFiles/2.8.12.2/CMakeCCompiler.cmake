set(CMAKE_C_COMPILER "/u/local/compilers/gcc/7.5.0/bin/gcc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "7.5.0")
set(CMAKE_C_PLATFORM_ID "Linux")

set(CMAKE_AR "/u/local/compilers/gcc/7.5.0/bin/ar")
set(CMAKE_RANLIB "/u/local/compilers/gcc/7.5.0/bin/ranlib")
set(CMAKE_LINKER "/u/local/compilers/gcc/7.5.0/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()




set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/u/local/compilers/gcc/7.5.0/lib64;/u/local/compilers/gcc/7.5.0/lib/gcc/x86_64-pc-linux-gnu/7.5.0;/lib64;/usr/lib64;/u/local/compilers/gcc/7.5.0/lib;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.8;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/daal/lib/intel64_lin;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib/release;/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/lib;/u/local/compilers/intel/2020.4/debugger_2020/libipt/intel64/lib;/u/local/compilers/intel/2020.4/debugger_2020/python/intel64/lib;/u/local/compilers/gcc/7.5.0/x86_64-pc-linux-gnu/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")



