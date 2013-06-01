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

#ifndef __PCA_HH__
#define __PCA_HH__

#include <cstdlib>
#include <Eigen/Dense>
#include "iterator.hh"
#include "parse.hh"

namespace hashpca 
{
  template<typename Hash,
           typename Iterator>
  std::pair<double, uint64_t>
  pca_accumulate (std::istream&            in,
                  hashpca::MatrixXd&       Y,
                  const hashpca::MatrixXd& Omega,
                  Iterator&                iterator,
                  bool                     center,
                  Eigen::VectorXd&         sum)
    {
      std::pair<double, uint64_t> rv (0.0, 0);
      veedubparse::StandardParse<Hash> parse;
      SugaryVectorXd xOmega (Omega.cols ());

      if (center)
        sum.setZero (Omega.rows ());

      parse (in,
             [ & ] (bool ok,
                    const veedubparse::GeneralExample& ex) {
               double imp;

               if (ok &&
                  [] (char* s, double& z) { 
                    if (! s) { z = 1; return true; }
                    char* ep; z = ::strtod (s, &ep);
                    return ep != s && (*ep == '\0' || isspace (*ep));
                  } (ex.importance, imp))
                 {
                   auto x = iterator (ex);

                   xOmega = x.transpose () * Omega;
                   Y += imp * x * xOmega;
                   if (center)
                     sum += imp * x.transpose ();
                   rv.first += imp;
                 }
                 else {
                   ++rv.second;
                 }
              });

      double invn = 1.0 / rv.first;

      if (center)
        {
          Eigen::VectorXd c (Omega.cols ());
          c = sum.transpose () * Omega;
          Y.noalias () -= invn * sum * c;
        }

      Y *= invn;

      return rv;
    }

  template<typename Hash,
           typename Iterator>
  int 
  computeu (std::istream&            in,
            std::ostream&            out,
            const hashpca::MatrixXd& V,
            const Eigen::VectorXd&   mean,
            const Eigen::VectorXd&   sinv,
            Iterator                 iterator,
            bool                     tanhify,
            bool                     normalize,
            bool                     flush)
    {
      veedubparse::StandardParse<Hash> parse;
      SugaryVectorXd Vtx (V.cols ());
      Eigen::VectorXd u (V.cols ());

      parse (in,
             [ & ] (bool ok,
                    const veedubparse::GeneralExample& ex) {
               if (ok)
                 {
                   auto x = iterator (ex);

                   Vtx = x.transpose () * V;
                   u = sinv.asDiagonal () * (Vtx - mean);

                   if (tanhify)
                     u = u.unaryExpr ([] (double x) { return tanh (0.85 * x); });

                   if (normalize)
                     u.normalize ();

                   if (ex.tag)
                     out << ex.tag;

                   for (unsigned int i = 0; i < V.cols (); ++i)
                     out << " " << (i+1) << ":" << u (i);

                   out << std::endl;
                   
                   if (flush)
                     out << std::flush;
                 }
              });

      return 0;
    }
}

#endif // __PCA_HH__
