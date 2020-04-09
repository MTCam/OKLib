///
/// @file
/// @ingroup plascom2_group
/// @brief Testing utilities for PlasCom2.
///
#ifndef __PC2_TWENKE_H__
#define __PC2_TWENKE_H__

#include "ComLine.H"

namespace twenke {
  
  ///
  /// ComLineObject for PlasCom2 testing command-line interface.
  ///
  class comline : public ix::util::ComLineObject
  {
  public:
    comline()
      : ix::util::ComLineObject()
    {};
    comline(const char *args[])
      : ix::util::ComLineObject(args)
    {};
    void Initialize(){
      AddOption('h',"help");
      AddHelp("help","Print out long version of help and exit.");
      AddOption('v',"verblevel",1,"level");
      AddHelp("verblevel","Set the verbosity level. (default = 0)");
      AddOption('c',"config",2,"filename");
      AddHelp("config","Set the configuration file to <filename>. (default = twenke.config)");
      _description.assign(std::string("Command-line interface to test PlasCom2 kernels and WENO."));
    };
  };
};
#endif
