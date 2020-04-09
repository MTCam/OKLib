#ifndef __TEST_UTIL_H__
#define __TEST_UTIL_H__

namespace testutil {

  // To avoid error-prone manual cleanup of things (PETSc, Overkit contexts, etc.) when exiting
  // tests early
  template <typename T, typename F> class handle {
  public:
    handle(T &object, F destroy):
      objectPtr_(&object),
      destroy_(destroy)
    {}
    ~handle() {
      destroy_(*objectPtr_);
    }
  private:
    T *objectPtr_;
    F destroy_;
  };

}

#endif
