Place the private header files in this directory. They will be
available to your code with

     #include <IMP/bhm/internal/myheader.h>

All headers should include `IMP/bhm/bhm_config.h` as their
first include and surround all code with `IMPBHM_BEGIN_INTERNAL_NAMESPACE`
and `IMPBHM_END_INTERNAL_NAMESPACE` to put it in the
IMP::bhm::internal namespace and manage compiler warnings.
