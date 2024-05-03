Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/bhm/myheader.h>

All headers should include `IMP/bhm/bhm_config.h` as their
first include and surround all code with `IMPBHM_BEGIN_NAMESPACE`
and `IMPBHM_END_NAMESPACE` to put it in the IMP::bhm namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
