/* -*- c++ -*-
 * Copyright (c) 2012-2016 by the GalSim developers team on GitHub
 * https://github.com/GalSim-developers
 *
 * This file is part of GalSim: The modular galaxy image simulation toolkit.
 * https://github.com/GalSim-developers/GalSim
 *
 * GalSim is free software: redistribution and use in source and binary forms,
 * with or without modification, are permitted provided that the following
 * conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions, and the disclaimer given in the accompanying LICENSE
 *    file.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions, and the disclaimer given in the documentation
 *    and/or other materials provided with the distribution.
 */

#include "IgnoreWarnings.h"

#include "Python.h"

#define BOOST_NO_CXX11_SMART_PTR
#include "boost/python.hpp"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL GALSIM_ARRAY_API
// This is the only one that doesn't have NO_IMPORT_ARRAY.
#include "numpy/arrayobject.h"

static int doImportNumpy() {
    import_array1(0);
    return 0;
}

namespace galsim {

    namespace hsm {
        void pyExportHSM();
    } // namespace hsm


} // namespace galsim

BOOST_PYTHON_MODULE(_hsm) {
    doImportNumpy();

    galsim::hsm::pyExportHSM();
}
