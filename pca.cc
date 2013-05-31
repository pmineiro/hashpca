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

#define NDEBUG
#include "pca.hh"
#include "pcaoptions.hh"
#include <cctype>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <Eigen/Eigenvalues>

#include <fenv.h>

#define FUDGE 10

namespace
{
  using namespace hashpca;

  // http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf

  /* Implementation of a 32-bit KISS generator which uses no multiply instructions */
  unsigned int JKISS32 ()
  {
    static uint32_t x=123456789,y=234567891,z=345678912,w=456789123,c=0;
    int32_t t;
    y ^= (y<<5); y ^= (y>>7); y ^= (y<<22); 
    t = z+w+c; z = w; c = t < 0; w = t&2147483647;
    x += 1411392427;
    return x + y + w;
  }

  double urandom ()
    {
      return .0000000002328306089594001 * JKISS32() + FLT_EPSILON;
    }

  //This is faster than Box-Muller but not benchmarked against inverse cdf
  double
  gensample (double)
    {
      double u,v,x,x2;

      do
        {
          u = urandom();
          v = 1.71552776992141359294*urandom()-.85776388496070679647;
          x = v/u;
          x2 = x*x;
        }
      while(x2>6-8*u+2*u*u && (x2>2/u-2*u || x2>-4*log(u)));

      return x;
    }

  void
  orthonormalize (Eigen::MatrixXd& Y)
    {
      // // http://forum.kde.org/viewtopic.php?f=74&t=91271
      // // Unfortunately this appears memory intensive :(
      // Eigen::HouseholderQR<Eigen::MatrixXd> qr (Y);
      // Y = qr.householderQ () * Eigen::MatrixXd::Identity (Y.rows (), Y.cols ());

      // Gram-Schmidt has no space overhead
      for (unsigned int j = 0; j < Y.cols (); ++j)
        {
          for (unsigned int i = 0; i < j; ++i)
            Y.col (j) -= Y.col (i).dot (Y.col (j)) * Y.col (i);

          for (unsigned int i = 0; i < j; ++i)
            Y.col (j) -= Y.col (i).dot (Y.col (j)) * Y.col (i);

          Y.col (j).normalize ();
        }
    }

  Eigen::VectorXd
  smallsvd (const Eigen::MatrixXd& Z,
            Eigen::MatrixXd&       V)
    {
      Eigen::MatrixXd ZtZ = Z.transpose () * Z;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es (ZtZ);
      Eigen::VectorXd s = es.eigenvalues ();
      Eigen::MatrixXd Upsilon = es.eigenvectors ();

      s = s.unaryExpr ([] (double x) { return x <= 0.0 ? 0.0 : sqrt (x); });

      Eigen::VectorXd pinv = s.unaryExpr ([] (double x) { return 1.0 / (1e-6 + x); });

      V.noalias () = Z * Upsilon * pinv.asDiagonal ();

      s = s.unaryExpr ([] (double x) { return x <= 0.0 ? 0.0 : sqrt (x); });

      return s;
    }

  int  
  writemodel (PcaOptions             options,
              std::ostream&          model,
              const Eigen::MatrixXd& V,
              const Eigen::VectorXd& s,
              const Eigen::VectorXd& sum,
              double                 n)
    {
      std::cerr << "writing model ... ";

      model << V.rows () << " " << V.cols () - FUDGE << std::endl;
      model 
        << (options.dashq ? static_cast<unsigned char> (options.dashq[0]) : -1) 
        << " "
        << (options.dashq ? static_cast<unsigned char> (options.dashq[1]) : -1) 
        << std::endl;

      model << s.reverse ().topRows (s.rows () - FUDGE) << std::endl;

      Eigen::VectorXd mean;

      if (options.center)
        mean = (1.0 / n) * V.transpose () * sum;
      else
        mean = Eigen::VectorXd::Zero (V.cols ());

      mean = mean.reverse ();

      // http://forum.kde.org/viewtopic.php?f=74&t=107161

      model.write ((char*) & mean (0), sizeof (double) * (mean.rows () - FUDGE));

      for (unsigned int i = 0; i < V.rows (); ++i)
        {
          Eigen::VectorXd tmp = V.row (i).reverse ();
          model.write ((char*) & tmp (0),
                       sizeof (double) * (V.cols () - FUDGE));
        }

      std::cerr << "done." << std::endl;

      return 0;
    }

  char* 
  make_dashq (int a,
              int b)
    {
      static char buf[3];

      buf[0] = a;
      buf[1] = b;
      buf[2] = 0;

      return buf;
    }

  int  
  readmodel (PcaOptions       options,
             std::istream&    model,
             Eigen::MatrixXd& V,
             Eigen::VectorXd& s,
             Eigen::VectorXd& mean)
    {
      model >> options.hashsize >> options.rank;

      int a;
      int b;

      model >> a >> b;

      if (a >= 0 && b >= 0)
        options.dashq = make_dashq (a, b);
      else
        options.dashq = 0;

      s.resize (options.rank);
      for (unsigned int i = 0; i < options.rank; ++i)
        model >> s (i);

      char newline;
      model.read (&newline, 1);

      mean.resize (options.rank);
      model.read ((char*) & mean (0), sizeof (double) * options.rank);

      Eigen::VectorXd tmp (options.rank);
      V.resize (options.hashsize, options.rank);

      for (unsigned int i = 0; i < V.rows (); ++i)
        {
          model.read ((char*) & tmp (0), sizeof (double) * options.rank);
          V.row (i) = tmp;
        }

      return 0;
    }

  template<typename Iterator>
  int 
  computeu (PcaOptions options,
            std::istream&          in,
            std::ostream&          out,
            const Eigen::MatrixXd& V,
            const Eigen::VectorXd& mean,
            const Eigen::VectorXd& sinv,
            Iterator               iterator)
    {
      veedubparse::StandardParse<veedubparse::HashString> parse;
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

                   if (options.tanhify)
                     u = u.unaryExpr ([] (double x) { return tanh (0.85 * x); });

                   if (options.normalize)
                     u.normalize ();

                   if (ex.tag)
                     out << ex.tag;

                   for (unsigned int i = 0; i < options.rank; ++i)
                     out << " " << (i+1) << ":" << u (i);

                   out << std::endl << std::flush;
                 }
              });

      return 0;
    }

  int
  do_project (PcaOptions    options,
              int           argc,
              char*         argv[])
    {
      (void) argc;

      std::ifstream in (argv[0]);

      if (! in.good ())
        {
          std::cerr << "ERROR: failed to open file '" 
                    << argv[0] << "' for reading: " 
                    << strerror (errno) << std::endl;
          return 1;
        }

      std::ifstream model (options.model);

      if (! model.good ())
        {
          std::cerr << "ERROR: failed to open file '" 
                    << options.model << "' for reading: " 
                    << strerror (errno) << std::endl;
          return 1;
        }

      Eigen::MatrixXd V;
      Eigen::VectorXd s;
      Eigen::VectorXd mean;

      if (readmodel (options, model, V, s, mean))
        {
          std::cerr << "ERROR: model file is corrupt or invalid." << std::endl;
          return 1;
        }

      Eigen::VectorXd sinv = s.unaryExpr ([] (double x) { return (x < 1e-6) ? 0 : 1.0 / x; });

      if (options.dashq)
        {
          BiLinearIterator it (options.hashsize,
                               static_cast<unsigned char> (options.dashq[0]),
                               static_cast<unsigned char> (options.dashq[1]));

          return computeu (options, in, std::cout, V, mean, sinv, it);
        }
      else
        {
          LinearIterator it (options.hashsize);

          return computeu (options, in, std::cout, V, mean, sinv, it);
        }
    }

  int
  do_pca (PcaOptions    options,
          int           argc,
          char*         argv[])
    {
      std::ifstream in (argv[0]);

      if (! in.good ())
        {
          std::cerr << "ERROR: failed to open file '" 
                    << argv[0] << "' for reading: " 
                    << strerror (errno) << std::endl;
          return 1;
        }

      std::ifstream in2 (argv[argc > 1 ? 1 : 0]);

      if (! in2.good ())
        {
          std::cerr << "ERROR: failed to open file '" 
                    << argv[argc > 1 ? 1 : 0] << "' for reading: " 
                    << strerror (errno) << std::endl;
          return 1;
        }

      std::ofstream model (options.model);

      if (! model.good ())
        {
          std::cerr << "ERROR: failed to open file '" 
                    << options.model << "' for writing: " 
                    << strerror (errno) << std::endl;
          return 1;
        }

      options.rank += FUDGE;

      Eigen::MatrixXd Y; Y.setZero (options.hashsize, options.rank);
      Eigen::MatrixXd Omega (options.hashsize, options.rank);
      Eigen::VectorXd sum;

      Omega = Omega.unaryExpr (std::ptr_fun (gensample));

      std::cerr << "Processing examples ... ";
      std::pair<double, uint64_t> lines1;

      if (options.dashq)
        {
          BiLinearIterator it (options.hashsize,
                               static_cast<unsigned char> (options.dashq[0]),
                               static_cast<unsigned char> (options.dashq[1]));

          lines1 = pca_accumulate (in, Y, Omega, it, options.center, sum);
        }
      else
        {
          LinearIterator it (options.hashsize);

          lines1 = pca_accumulate (in, Y, Omega, it, options.center, sum);
        }

      std::cerr << "(" << lines1.first << ", " 
                << lines1.second << ") done." << std::endl;

      std::cerr << "Orthogonalizing ... ";
      orthonormalize (Y);
      std::cerr << "done." << std::endl;

      std::cerr << "Reprocessing examples ... ";
      Eigen::MatrixXd& Z (Omega);
      Z.setZero (Z.rows (), Z.cols ());

      std::pair<double, uint64_t> lines2;

      if (options.dashq)
        {
          BiLinearIterator it (options.hashsize,
                               static_cast<unsigned char> (options.dashq[0]),
                               static_cast<unsigned char> (options.dashq[1]));

          lines2 = pca_accumulate (in2, Z, Y, it, options.center, sum);
        }
      else
        {
          LinearIterator it (options.hashsize);

          lines2 = pca_accumulate (in2, Z, Y, it, options.center, sum);
        }

      std::cerr << "(" << lines2.first << ", " 
                << lines2.second << ") done." << std::endl;

      std::cerr << "Small svding + postprocessing ... ";
      Eigen::MatrixXd& V (Y);
      Eigen::VectorXd s = smallsvd (Z, V);
      std::cerr << "done." << std::endl;

      options.rank -= FUDGE;
      return writemodel (options, model, V, s, sum, lines2.first);
    }
}

int 
main (int   argc,
      char* argv[])
{
  std::ios::sync_with_stdio (false);

  if (argc < 2)
    {
      std::cerr << help << std::endl;
      return 1;
    }

  PcaOptions options = parse_pca_options (argc, argv);

  if (argc < 1)
    {
      std::cerr << "ERROR: did not specify input file\n" << std::endl;
      std::cerr << help << std::endl;
      return 1;
    }

  if (! options.model)
    {
      std::cerr << "ERROR: did not specify model file\n" << std::endl;
      std::cerr << help << std::endl;
      return 1;
    }

  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  return options.project ? do_project (options, argc, argv)
                         : do_pca (options, argc, argv);
}
