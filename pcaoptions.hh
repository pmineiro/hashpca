#ifndef __PCAOPTIONS_HH__
#define __PCAOPTIONS_HH__

#include <cstdlib>
#include <iostream>

namespace hashpca 
{
  const char* help="Usage: pca [options] input [input2]\nAvailable options:\n\
                    -b <int>  : size of feature hash (default: 65536)\n\
                    -k <int>  : rank of approximation (default: 40)\n\
                    -m <name> : model file\n\
                    -t        : projection mode\n\
                    -f        : flush after each output line\n\
                    -e        : evidence normalize\n\
                    -s        : tanh(0.85)ify\n\
                    -c        : center data\n\
                    -w        : do not whiten projection\n\
                    -a        : hash all features (including integers)\n\
                    -q ab     : pair features from a and b\n";

  struct PcaOptions
    {
      unsigned int hashsize;
      unsigned int rank;
      const char* model;
      bool project;
      bool flush;
      bool normalize;
      bool tanhify;
      bool center;
      bool whiten;
      bool hashall;
      const char* dashq;
  
      PcaOptions () : hashsize (65536),
                      rank (40),
                      model (0),
                      project (false),
                      flush (false),
                      normalize (false),
                      tanhify (false),
                      center (false),
                      whiten (true),
                      hashall (false),
                      dashq (0)
        {
        }
    };

  int
  parse_int (int   argc,
             char* argv[])
    {
      if (argc < 1)
        {
          std::cerr << "ERROR: missing expected integer argument" << std::endl;
          std::cerr << help << std::endl;
          exit (1);
        }

      char* endptr;
      int rv = ::strtol (argv[0], &endptr, 0);

      if (endptr == argv[0] || (*endptr != '\0' && ! ::isspace (*endptr)))
        {
          std::cerr << "ERROR: invalid integer argument '" 
                    << argv[0] << "'" << std::endl;
          std::cerr << help << std::endl;
          exit (1);
        }

      return rv;
    }

  char*
  parse_string (int   argc,
                char* argv[])
    {
      if (argc < 1)
        {
          std::cerr << "ERROR: missing expected string argument" << std::endl;
          std::cerr << help << std::endl;
          exit (1);
        }

      return argv[0];
    }

  PcaOptions
  parse_pca_options (int&    argc,
                     char**& argv)
    {
      PcaOptions options;

      --argc;
      ++argv;

      while (argc > 0 && argv[0][0] == '-')
        {
          switch (argv[0][1])
            {
              case '-':
                --argc;
                ++argv;
                return options;
              case 'b':
                --argc;
                ++argv;
                options.hashsize = parse_int (argc, argv);
                break;
              case 'k':
                --argc;
                ++argv;
                options.rank = parse_int (argc, argv);
                break;
              case 'e':
                options.normalize = true;
                break;
              case 's':
                options.tanhify = true;
                break;
              case 't':
                options.project = true;
                break;
              case 'f':
                options.flush = true;
                break;
              case 'c':
                options.center = true;
                break;
              case 'w':
                options.whiten = false;
                break;
              case 'a':
                options.hashall = true;
                break;
              case 'm':
                --argc;
                ++argv;
                options.model = parse_string (argc, argv);
                break;
              case 'q':
                --argc;
                ++argv;
                options.dashq = parse_string (argc, argv);
                break;
              default:
                std::cerr << "ERROR: unrecognized switch " << argv[0] << std::endl;
                std::cerr << help << std::endl;
                exit (1);
                break;
            }

          --argc;
          ++argv;
        }

      return options;
    }
}

#endif // __PCAOPTIONS_HH__
