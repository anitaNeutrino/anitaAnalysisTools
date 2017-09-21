#include "AcclaimCmdLineArgs.h"

#include <unistd.h>
#include <getopt.h>
#include <iostream>

#include "AnitaVersion.h"
#include "AnalysisFlow.h"
#include "RootTools.h"

// Values used if value not specified on command line
const int default_run = 352;
const int default_numdivisions = 1;
const int default_subdivision = 0;
const int default_anitaversion = 3;
const char* default_settings_filename = "";

/** 
 * Prints available options and default values
 * 
 * @param argv0 should be the program name, i.e. argv[0]
 */
void Acclaim::CmdLineArgs::printHelp(const char* argv0){

  // should probably use printf to get formatting better... oh well
  std::cerr << std::endl;
  std::cerr << argv0 << " options:" << std::endl;
  std::cerr << "-h     or --help            : print this help message" << std::endl;
  std::cerr << "-a " << default_anitaversion << "   or --anitaversion=" << default_anitaversion << "  : sets the initial anita version" << std::endl;
  std::cerr << "-r " << default_run << " or --run=" << default_run << "         : which run to process" << std::endl;
  std::cerr << "-m " << default_subdivision << "   or --subdivision=" << default_subdivision << "   : This job processes only part of the run, subdivision m of n divisions (where 0 <= m < n)" << std::endl;
  std::cerr << "-n " << default_numdivisions << "   or --numdivisions=" << default_numdivisions << "  : This job processes only part of the run, subdivision m of n divisions (where 0 <= m < n)" << std::endl;
  std::cerr << "-s " << default_settings_filename << "   or --settings=" << default_settings_filename << "  : Path to an AcclaimSettings.conf file (or directory containing), see AnalysisSettings class for default locations searched if not specified." << std::endl;
  std::cerr << "-o " << output_filename << "   or --settings=" << output_filename << "  : Path to an AcclaimSettings.conf file (or directory containing), see AnalysisSettings class for default locations searched if not specified." << std::endl;    
  std::cerr << std::endl;
}



/** 
 * Leverages getopt.h to actually do the parsing
 * 
 * @param argc passed from main (via constructor)
 * @param argv passed from main (via constructor)
 */
void Acclaim::CmdLineArgs::getArgs(int argc, char* argv[]){

  int c = 0;

  run = default_run;
  numdivisions = default_numdivisions;
  division = default_subdivision;
  anitaversion = default_anitaversion;
  settings_filename = default_settings_filename;
  event_selection = AnalysisFlow::kDecimated;

  std::vector<TString> argv0Tokens;
  RootTools::tokenize(argv0Tokens, argv[0], "/");
  output_filename = argv0Tokens.back();
  
  while (true){

    const int numEventSelectionFlags = 3; // should all go at the front
    static struct option long_options[] = {{"all",          no_argument,       &event_selection, AnalysisFlow::kAll},
                                           {"decimated",    no_argument,       &event_selection, AnalysisFlow::kDecimated},
                                           {"wais",         no_argument,       &event_selection, AnalysisFlow::kWaisPulser},
                                           {"help",         no_argument,       NULL,              'h'},
                                           {"anitaversion", required_argument, NULL,              'a'},
                                           {"run",          required_argument, NULL,              'r'},
                                           {"numdivisions", required_argument, NULL,              'n'},
                                           {"subdivision",  required_argument, NULL,              'm'},
                                           {"settings",     required_argument, NULL,              's'},
                                           {"output",       required_argument, NULL,              'o'},
                                           {0, 0, 0, 0}};
    
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "ha:r:n:d:s:",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1){
      break;
    }
    
    switch (c){
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0){
          if(option_index < numEventSelectionFlags){
            std::cout << "Info in " << __PRETTY_FUNCTION__ << ", got " << long_options[option_index].name <<  ", so set event_selection = " << event_selection << std::endl;
          }          
          break;
        }
        
        // std::cout << "option " << long_options[option_index].name << std::endl;
        // if (optarg){
        //   std::cout  << " with arg " <<  optarg << std::endl;
        // }
        // std::cout << std::endl;
        // break;

      case 'h':
        printHelp(argv[0]);
        exit(0);
        break;

      case 'a':
        anitaversion = atoi(optarg);
        std::cout << "Info in " << __PRETTY_FUNCTION__ << ", set anitaversion = " << anitaversion << std::endl;
        AnitaVersion::set(anitaversion);
        break;

      case 'r':
        run = atoi(optarg);
        std::cout << "Info in " << __PRETTY_FUNCTION__ << ", set run = " << run << std::endl;
        break;

      case 'n':
        numdivisions = atoi(optarg);
        std::cout << "Info in " << __PRETTY_FUNCTION__ << ", set numdivisions = " << numdivisions << std::endl;
        break;

      case 'm':
        division = atoi(optarg);        
        std::cout << "Info in " << __PRETTY_FUNCTION__ << ", set division = " << division << std::endl;
        break;

      case 's':
        settings_filename = TString::Format("%s", optarg);
        std::cout << "Info in " << __PRETTY_FUNCTION__ << ", set settings_filename = " << settings_filename << std::endl;

      case '?':
        // printHelp(argv[0]);
        /* getopt_long already printed an error message. */
        break;

      default:
        abort();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    std::cerr << "non-option ARGV-elements: ";
    while (optind < argc){
      std::cerr <<  argv[optind++] << " ";
    }
    std::cerr << std::endl;
  }
}


/** 
 * Constructor
 * 
 * @param argc passed from main
 * @param argv passed from main
 */
Acclaim::CmdLineArgs::CmdLineArgs(int argc, char* argv[]){

  getArgs(argc, argv);
  checkArgs(argv[0]);
  
}


void Acclaim::CmdLineArgs::checkArgs(const char* argv0){

  // if(anitaversion!=3){
  //   std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", currently only checking arguments for ANITA-3" << std::endl;
  // }
  
  Bool_t mustDie = false; // will set to true if have nonsense, program will then die

  if(numdivisions <= 0){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", require numdivisions to be a positive integer" << std::endl;
    mustDie = true;
  }
  
  if(division < 0 || division >= numdivisions){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", require 0 <= division < numdivisions" << std::endl;
    mustDie = true;
  }

  if(mustDie){
    printHelp(argv0);
    exit(1);
  }

  
  
}

