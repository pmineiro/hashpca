/*=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
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

#ifndef __PARSE_HH__
#define __PARSE_HH__

#include <cstdlib>
#include <cstring>
#include <istream>
#include <utility>
#include <vector>

#include "example.hh"
#include "MurmurHash3.hh"

namespace veedubparse
{
  struct HashAll
    {
      uint64_t
      operator() (const char* s,
                  uint64_t    seed)
        {
          uint64_t hv[2];

          MurmurHash3_x64_128 (s, strlen (s), seed, hv);

          return hv[0];
        }
    };

  struct HashString
    {
      uint64_t
      operator() (const char* s,
                  uint64_t    seed)
        {
          uint64_t hv[2];

          const char* p;
          for (hv[0] = 0, p = s; isdigit (*p); ++p)
            {
              hv[0] *= 10;
              hv[0] += *p - '0';
            }

          if (*p == '\0')
            {
              return hv[0] + seed;
            }

          MurmurHash3_x64_128 (s, strlen (s), seed, hv);

          return hv[0];
        }
    };

  inline bool
  iswhite (char c)
    {
      return c == ' ' || c == '\t';
    }

  inline char*
  strrchr_back (char* start,
                char* end,
                char  value)
    {
      while (end > start)
        {
          --end;
          if (*end == value)
            {
              return end;
            }
        }

      return 0;
    }

  template<typename Hash>
  class StandardParse
    {
      private:
        Hash hash;

        class Problem { };

        char*
        parse_label (char*& buf)
          {
            if (iswhite (*buf) || *buf == '\0') { return 0; }

            char* rv = buf;
            while (! iswhite (*buf) && *buf != '\0') { ++buf; }
            if (*buf != '\0') { *buf++ = '\0'; }
            return rv;
          }

        std::pair<char*, char*>
        parse_importance_tag (char*& buf)
          {
            typedef std::pair<char*, char*> CharPtrPair;

            if (*buf == '|')
              {
                ++buf;
                return CharPtrPair (0, 0);
              }
            else if (*buf == '\0')
              {
                return CharPtrPair (0, 0);
              }
            else
              {
                char* first = buf;
                while (! iswhite (*buf) && *buf != '|' && *buf != '\0') 
                  { 
                    ++buf; 
                  }

                if (*buf == '\0') 
                  { 
                    return CharPtrPair (0, first);
                  }

                if (*buf == '|')
                  {
                    *buf++ = '\0';
                    return CharPtrPair (0, first);
                  }
                else // iswhite (*buf)
                  {
                    *buf++ = '\0';
                    while (iswhite (*buf) && *buf != '\0') { ++buf; }

                    if (*buf == '\0') 
                      {
                        return CharPtrPair (first, 0);
                      }
  
                    if (*buf == '|')
                      {
                        ++buf;
                        return CharPtrPair (first, 0);
                      }
                    else // ! iswhite (*buf)
                      {
                        char* second = buf;
                        while (*buf != '|' && *buf != '\0') { ++buf; }
                        if (*buf == '\0') 
                          { 
                            return CharPtrPair (first, second);
                          }

                        *buf++ = '\0';

                        return CharPtrPair (first, second);
                      }
                  }
              }
          }
  
      public:
        StandardParse () { }

        template<typename F>
        void
        operator() (std::istream& in,
                    const F&      f)
          {
            GeneralExample example;

            example.buf.reset (new std::string);

            while (in.good ())
              {
                try
                  {
                    example.reset ();

                    std::getline (in, *(example.buf));
                    if (! in.good ()) { return; }

                    char* buf = &(*example.buf)[0];
                    if (*buf == '\0') { throw Problem (); }

                    example.label = parse_label (buf);

                    while (iswhite (*buf) && *buf != '\0') { ++buf; }

                    std::pair<char*, char*> i_t (parse_importance_tag (buf));

                    example.importance = i_t.first;
                    example.tag = i_t.second;

                    while (*buf != '\0')
                      {
                        const char* ns = " ";
                        float nsweight = 1;
                        
                        if (! iswhite (*buf))
                          {
                            ns = buf;
                            while (! iswhite (*buf) && *buf != '\0') { ++buf; }
                            if (*buf == '\0') { break; }
                            *buf++ = '\0';

                            char* colon = strrchr_back (const_cast<char*> (ns), buf, ':');

                            if (colon)
                              {
                                char* endptr;

                                *colon++ = '\0';
                                nsweight = ::strtof (colon, &endptr);

                                if (*colon == '\0' || *endptr != '\0')
                                  {
                                    throw Problem ();
                                  }
                              }
                          }

                        uint64_t seed = hash (ns, 0);

                        while (iswhite (*buf) && *buf != '\0') { ++buf; }
                        if (*buf == '\0') { break; }
          
                        while (*buf != '|' && *buf != '\0')
                          {
                            char* feat = buf;
                            float weight = 1;

                            while (! iswhite (*buf) && *buf != '\0') { ++buf; }
                            if (*buf != '\0') { *buf++ = '\0'; }

                            char* colon = strrchr_back (feat, buf, ':');

                            if (colon)
                              {
                                char* endptr;

                                *colon++ = '\0';
                                weight = ::strtof (colon, &endptr);

                                if (*colon == '\0' || *endptr != '\0')
                                  {
                                    throw Problem ();
                                  }
                              }

                            uint64_t hv = hash (feat, seed);
                            float xi = (hv & (1ULL<<63)) ? 1 : -1;
                            unsigned int n = static_cast<unsigned char> (ns[0]);

                            example.f[n].push_back (Feature (hv, xi * weight * nsweight));
                          }

                        if (*buf == '|')
                          {
                            ++buf;
                          }
                      }

                    f (true, example);
                  }
                catch (Problem&)
                  {
                    f (false, example);
                  }
              }
          }
      };
}

#endif // __PARSE_HH__
