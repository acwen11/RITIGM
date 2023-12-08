include_directories("/home1/apps/intel19/impi19_0/fftw3/3.3.10/include")
include_directories("/opt/apps/intel19/gsl/2.6/include")
include_directories("/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/include")
include_directories("/home1/apps/gcc9_1/boost/1.83.0/include")


set (PGPLOT_LIBRARIES "-L/home1/09319/acw6923/pgplot/lib/libpgplot.so -lgfortran")
set (GSL_LIBRARIES "/opt/apps/intel19/gsl/2.6/lib/libgsl.so")
set (SCALAPACK_LIBRARIES "-L/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")
#set (BLACS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblacs-openmpi.so")
set (FFTW_LIBRARIES "/home1/apps/intel19/impi19_0/fftw3/3.3.10/lib/libfftw3.so")
#set (BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so")
#set (LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so")
