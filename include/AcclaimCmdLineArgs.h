#ifndef ACCLAIM_CMD_LINE_ARGS
#define ACCLAIM_CMD_LINE_ARGS

#include "TString.h"

namespace Acclaim{


/** 
 * @class CmdLineArgs 
 * @brief A simple command line option parser.
 * 
 * Warning! Could prevent a program that uses it
 * from running if it doesn't like what it gets
 * 
 * int main(int argc, char* argv[]){
 *   Acclaim::CmdLineArgs args(argc, argv);
 *   args.run; // this is the run that was passed 
 *   ...
 */
class CmdLineArgs {

 public:
  CmdLineArgs(int argc, char* argv[]);
  void printHelp(const char* argv0);
  void getArgs(int argc, char* argv[]);
  void checkArgs(const char* argv0);

  int event_selection;
  int run;
  int numdivisions;
  int division;
  int anitaversion;
  TString settings_filename;
  TString output_filename;  
  
};



}

#endif
