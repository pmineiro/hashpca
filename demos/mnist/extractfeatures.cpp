#include <arpa/inet.h>
#include <cassert>
#include <cmath>
#include <cstring>
#include <stdint.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

namespace {
  using std::setprecision;
  using std::stringstream;

  template<typename T>
  T
  square (T x)
    {
      return x * x;
    }

  char* numbers[256];

  void
  make_numbers ()
    {
      for (unsigned int i = 0; i < 256; ++i)
	{
	  stringstream ss;

	  ss << ":" << setprecision (8) << static_cast<double> (i)/256.0;

	  numbers[i] = ::strdup (ss.str ().c_str ());
	}
    }
}

int 
main (void)
{
  using std::cin;
  using std::cout;
  using std::endl;
  using std::setprecision;

  make_numbers ();

  uint32_t magic;
  cin.read (reinterpret_cast<char*> (&magic), sizeof (uint32_t));
  magic = ntohl (magic);
  assert (magic == 2051);

  uint32_t n_images;
  cin.read (reinterpret_cast<char*> (&n_images), sizeof (uint32_t));
  n_images = ntohl (n_images);

  uint32_t n_rows;
  cin.read (reinterpret_cast<char*> (&n_rows), sizeof (uint32_t));
  n_rows = ntohl (n_rows);

  uint32_t n_columns;
  cin.read (reinterpret_cast<char*> (&n_columns), sizeof (uint32_t));
  n_columns = ntohl (n_columns);

  uint32_t rc = n_rows * n_columns;
  unsigned char buf[rc];

  size_t n = 0;

  for (cin.read (reinterpret_cast<char*> (buf), rc); 
       ! cin.eof (); 
       cin.read (reinterpret_cast<char*> (buf), rc))
    {
      ++n;

//      double sumsq = 0;
//
//      for (unsigned int p = 0; p < n_rows * n_columns; ++p)
//        {
//          sumsq += square (buf[p]);
//        }
//
//      sumsq = sqrt (sumsq);

      for (unsigned int p = 0; p < n_rows * n_columns; ++p)
        {
          if (buf[p])
            cout << " " << p+1 << numbers[buf[p]];
//            cout << " " << p+1 << ":" << static_cast<double> (buf[p]) / sumsq;
        }

      cout << endl;
    }

  return 0;
}
