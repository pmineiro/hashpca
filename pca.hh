#ifndef __PCA_HH__
#define __PCA_HH__

#include <cstdlib>
#include <Eigen/Dense>
#include "iterator.hh"
#include "parse.hh"
#include "pcaoptions.hh"

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
            PcaOptions               options)
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
                   if (options.whiten)
                     u = sinv.asDiagonal () * (Vtx - mean);
                   else
                     u = Vtx - mean;

                   if (options.tanhify)
                     u = u.unaryExpr ([] (double x) { return tanh (0.85 * x); });

                   if (options.normalize)
                     u.normalize ();

                   if (ex.tag)
                     out << ex.tag;

                   for (unsigned int i = 0; i < V.cols (); ++i)
                     out << " " << (i+1) << ":" << u (i);

                   out << std::endl;
                   
                   if (options.flush)
                     out << std::flush;
                 }
              });

      return 0;
    }
}

#endif // __PCA_HH__
