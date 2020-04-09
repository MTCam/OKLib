#ifndef __CHEMISTRY_TEST_FIXTURES_H__
#define __CHEMISTRY_TEST_FIXTURES_H__

#include <cmath>

namespace testfixtures {
  namespace chemistrytest {

    typedef plascom2::application::domainvector        DomainVector;
    typedef simulation::grid::parallel_blockstructured pbsgrid_t;
    typedef simulation::state::base                    state_t;


#ifdef HAVE_CANTERA
    int GenerateUniformChemistryExact(pbsgrid_t &exactGrid, state_t &exactState);
    int GenerateUniformChemistryPlasmaExact(pbsgrid_t &exactGrid, state_t &exactState);
#endif

  }
}
#endif
