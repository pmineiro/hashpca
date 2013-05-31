/*=====================================================================*
 *                   Copyright (C) 2013 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

#ifndef __PCAOPTIONS_HH__
#define __PCAOPTIONS_HH__

#include <cstdlib>
#include <iostream>

namespace hashpca 
{
  const char* help="Usage: pca [options] input\nAvailable options:\n\
                    -b <int>  : size of feature hash (default: 65536)\n\
                    -k <int>  : rank of approximation (default: 40)\n\
                    -m <name> : model file\n\
                    -t        : projection mode\n\
                    -e        : evidence normalize\n\
                    -s        : tanh(0.85)ify\n\
                    -c        : center data\n\
                    -q ab     : pair features from a and b\n";

  struct PcaOptions
    {
      unsigned int hashsize;
      unsigned int rank;
      const char* model;
      bool project;
      bool normalize;
      bool tanhify;
      bool center;
      const char* dashq;
  
      PcaOptions () : hashsize (65536),
                      rank (40),
                      model (0),
                      project (false),
                      normalize (false),
                      tanhify (false),
                      center (false),
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
              case 'c':
                options.center = true;
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
