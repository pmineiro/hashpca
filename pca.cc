#define NDEBUG
#include "pca.hh"
#include "pcaoptions.hh"
#include <cctype>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <Eigen/Eigenvalues>

#define FUDGE 10

namespace
{
  using namespace hashpca;
  using veedubparse::HashAll;
  using veedubparse::HashString;

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
  orthonormalize (hashpca::MatrixXd& Y)
    {
      // // http://forum.kde.org/viewtopic.php?f=74&t=91271
      // // Unfortunately this appears memory intensive :(
      // Eigen::HouseholderQR<hashpca::MatrixXd> qr (Y);
      // Y = qr.householderQ () * hashpca::MatrixXd::Identity (Y.rows (), Y.cols ());

      Y.transposeInPlace ();

      // Gram-Schmidt has no space overhead
      for (unsigned int j = 0; j < Y.rows (); ++j)
        {
          for (unsigned int i = 0; i < j; ++i)
            Y.row (j) -= Y.row (i).dot (Y.row (j)) * Y.row (i);

          for (unsigned int i = 0; i < j; ++i)
            Y.row (j) -= Y.row (i).dot (Y.row (j)) * Y.row (i);

          Y.row (j).normalize ();
        }

      Y.transposeInPlace ();
    }

  Eigen::VectorXd
  smallsvd (const hashpca::MatrixXd& Z,
            hashpca::MatrixXd&       V)
    {
      hashpca::MatrixXd ZtZ = Z.transpose () * Z;
      Eigen::SelfAdjointEigenSolver<hashpca::MatrixXd> es (ZtZ);
      Eigen::VectorXd s = es.eigenvalues ();
      hashpca::MatrixXd Upsilon = es.eigenvectors ();

      s = s.unaryExpr ([] (double x) { return x <= 0.0 ? 0.0 : sqrt (x); });

      Eigen::VectorXd pinv = s.unaryExpr ([] (double x) { return 1.0 / (1e-6 + x); });

      V.noalias () = Z * Upsilon * pinv.asDiagonal ();

      s = s.unaryExpr ([] (double x) { return x <= 0.0 ? 0.0 : sqrt (x); });

      return s;
    }

  int  
  writemodel (PcaOptions                options,
              std::ostream&             model,
              const hashpca::MatrixXd&  V,
              const Eigen::VectorXd&    s,
              const Eigen::VectorXd&    sum,
              double                    n)
    {
      std::cerr << "Writing model ... ";

      model << V.rows () << " " << V.cols () - FUDGE << " " << options.hashall << std::endl;
      model 
        << (options.dashq ? static_cast<unsigned char> (options.dashq[0]) : -1) 
        << " "
        << (options.dashq ? static_cast<unsigned char> (options.dashq[1]) : -1) 
        << std::endl;

      Eigen::IOFormat HeavyFmt (Eigen::FullPrecision);
      model << s.reverse ().topRows (s.rows () - FUDGE).format (HeavyFmt) << std::endl;

      Eigen::VectorXd mean;

      if (options.center)
        mean = (1.0 / n) * V.transpose () * sum;
      else
        mean = Eigen::VectorXd::Zero (V.cols ());

      mean = mean.reverse ().eval ();

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
  readmodel (PcaOptions&        options,
             std::istream&      model,
             hashpca::MatrixXd& V,
             Eigen::VectorXd&   s,
             Eigen::VectorXd&   mean)
    {
      model >> options.hashsize >> options.rank >> options.hashall;

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

      hashpca::MatrixXd V;
      Eigen::VectorXd s;
      Eigen::VectorXd mean;

      if (readmodel (options, model, V, s, mean))
        {
          std::cerr << "ERROR: model file is corrupt or invalid." << std::endl;
          return 1;
        }

      Eigen::VectorXd pinv = s.unaryExpr ([] (double x) { return 1.0 / (1e-6 + x); });

      if (options.dashq)
        {
          BiLinearIterator it (options.hashsize,
                               static_cast<unsigned char> (options.dashq[0]),
                               static_cast<unsigned char> (options.dashq[1]));

          return (options.hashall) 
            ? computeu<HashAll> (in, std::cout, V, mean, pinv, it, options.tanhify, options.normalize, options.flush)
            : computeu<HashString> (in, std::cout, V, mean, pinv, it, options.tanhify, options.normalize, options.flush);
        }
      else
        {
          LinearIterator it (options.hashsize);

          return (options.hashall)
            ? computeu<HashAll> (in, std::cout, V, mean, pinv, it, options.tanhify, options.normalize, options.flush)
            : computeu<HashString> (in, std::cout, V, mean, pinv, it, options.tanhify, options.normalize, options.flush);
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

      hashpca::MatrixXd Y; Y.setZero (options.hashsize, options.rank);
      hashpca::MatrixXd Omega =
        hashpca::MatrixXd::Zero (options.hashsize, options.rank)
          .unaryExpr (std::ptr_fun (gensample));
      Eigen::VectorXd sum;

      std::cerr << "Processing examples ... ";
      std::pair<double, uint64_t> lines1;

      if (options.dashq)
        {
          BiLinearIterator it (options.hashsize,
                               static_cast<unsigned char> (options.dashq[0]),
                               static_cast<unsigned char> (options.dashq[1]));

          lines1 = (options.hashall)
            ? pca_accumulate<HashAll> (in, Y, Omega, it, options.center, sum)
            : pca_accumulate<HashString> (in, Y, Omega, it, options.center, sum);
        }
      else
        {
          LinearIterator it (options.hashsize);

          lines1 = (options.hashall)
            ? pca_accumulate<HashAll> (in, Y, Omega, it, options.center, sum)
            : pca_accumulate<HashString> (in, Y, Omega, it, options.center, sum);
        }

      std::cerr << "(" << lines1.first << ", " 
                << lines1.second << ") done." << std::endl;

      std::cerr << "Orthogonalizing ... ";
      orthonormalize (Y);
      std::cerr << "done." << std::endl;

      std::cerr << "Reprocessing examples ... ";
      hashpca::MatrixXd& Z (Omega);
      Z.setZero (Z.rows (), Z.cols ());

      std::pair<double, uint64_t> lines2;

      if (options.dashq)
        {
          BiLinearIterator it (options.hashsize,
                               static_cast<unsigned char> (options.dashq[0]),
                               static_cast<unsigned char> (options.dashq[1]));

          lines2 = (options.hashall)
            ? pca_accumulate<HashAll> (in2, Z, Y, it, options.center, sum)
            : pca_accumulate<HashString> (in2, Z, Y, it, options.center, sum);
        }
      else
        {
          LinearIterator it (options.hashsize);

          lines2 = (options.hashall)
            ? pca_accumulate<HashAll> (in2, Z, Y, it, options.center, sum)
            : pca_accumulate<HashString> (in2, Z, Y, it, options.center, sum);
        }

      std::cerr << "(" << lines2.first << ", " 
                << lines2.second << ") done." << std::endl;

      std::cerr << "Small svding + postprocessing ... ";
      hashpca::MatrixXd& V (Y);
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

  return options.project ? do_project (options, argc, argv)
                         : do_pca (options, argc, argv);
}
