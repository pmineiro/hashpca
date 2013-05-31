#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>

#include "../parse.hh"
#include "../MurmurHash3.cpp"

namespace
{
using namespace veedubparse;

bool
not_space (char p)
{
  return ! iswhite (p);
}

bool
not_space_or_colon (char p)
{
  return ! iswhite (p) && p != ':';
}

bool
not_space_or_bar (char p)
{
  return ! iswhite (p) && p != '|';
}

char*
random_string (bool (*test) (char))
{
  size_t n = 1 + lrand48 () % 16;
  char* buf = new char[1 + n];

  for (size_t i = 0; i < n; ++i)
    {
      do
        {
          buf[i] = lrand48 ();
//          buf[i] = ' ' + ('~' - ' ') * drand48 ();
        }
      while (buf[i] == '\n' || buf[i] == '\0' || ! test (buf[i]));
    }

  buf[n] = '\0';

  return buf;
}

std::ostream&
operator<< (std::ostream& o,
            const Feature& v)
{
  o << "{ " << v.id << " , " << v.value << " }";
  return o;
}

template<typename T>
std::ostream&
operator<< (std::ostream&         o,
            const std::vector<T>& v)
{
  for (auto x = v.begin (); x != v.end (); ++x)
    {
      o << *x << ", ";
    }

  return o;
}

typedef std::vector<std::pair<std::string, std::pair<bool, float> > > WeightedFeatureList;
typedef std::map<std::pair<std::string, std::pair<bool, float> >, WeightedFeatureList> AllFeatures;

std::ostream&
operator<< (std::ostream&       o,
            const AllFeatures&  af)
{
  for (AllFeatures::const_iterator x = af.begin (); x != af.end (); ++x)
    {
      o << "|" << x->first.first;

      if (x->first.second.first)
        {
          o << ":" << x->first.second.second;
        }

      o << " ";

      for (WeightedFeatureList::const_iterator y = x->second.begin ();
           y != x->second.end ();
           ++y)
        {
          o << y->first;
          if (y->second.first)
            {
              o << ":" << y->second.second;
            }
          o << " ";
        }
    }

  return o;
}

WeightedFeatureList
random_weighted_feature_list (void)
{
  WeightedFeatureList wf;
  size_t num_features = lrand48 () % 16;

  for (size_t i = 0; i < num_features; ++i)
    {
      char* feature = random_string (not_space_or_colon);
      while (feature[0] == '|')
        {
          delete[] feature;
          feature = random_string (not_space_or_colon);
        }

      std::string fs (feature);
      std::pair<bool, float> value;

      value.first = drand48 () < 0.5;
      value.second = drand48 ();

      wf.push_back (std::make_pair (fs, value));

      delete[] feature;
    }

  return wf;
}

AllFeatures
random_features (void)
{
  AllFeatures af;
  size_t num_ns = lrand48 () % 16;

  for (size_t i = 0; i < num_ns; ++i)
    {
      char* space = const_cast<char*> (" ");
      // TODO: add values to ns
      char* ns = (drand48 () < 0.2) ? space : random_string (not_space_or_colon);

      while (ns[0] == '|')
        {
          delete[] ns;
          ns = (drand48 () < 0.2) ? space : random_string (not_space_or_colon);
        }

      std::pair<bool, float> value;

      value.first = ns != space && drand48 () < 0.5;
      value.second = drand48 ();

      WeightedFeatureList fl = random_weighted_feature_list ();

      af[std::make_pair (std::string (ns), value)].insert
        (af[std::make_pair (std::string (ns), value)].end (),
         fl.begin (),
         fl.end ());

      if (ns != space) { delete[] ns; }
    }

  return af;
}

struct IndexedFeatureSet
{
  std::vector<Feature> spaces[256];

  const std::vector<Feature>&
  operator[] (const int i) const
    {
      return spaces[i];
    }

  std::vector<Feature>&
  operator[] (const int i)
    {
      return spaces[i];
    }
};

IndexedFeatureSet
get_all_features (const AllFeatures& af)
{
  IndexedFeatureSet f;

  for (AllFeatures::const_iterator x = af.begin (); x != af.end (); ++x)
    {
      const char* ns = x->first.first.c_str ();
      float nsvalue = x->first.second.first ? x->first.second.second : 1;
      char* endptr;
      uint64_t hv[2];
      uint64_t seed = strtoull (ns, &endptr, 10);

      if (! isdigit (*ns) || *endptr != '\0')
        {
          MurmurHash3_x64_128 (ns, strlen (ns), 0, hv);
          seed = hv[0];
        }

      for (WeightedFeatureList::const_iterator y = x->second.begin ();
           y != x->second.end ();
           ++y)
        {
          const char* feat = y->first.c_str ();
          float value = y->second.first ? y->second.second : 1;

          hv[0] = strtoull (feat, &endptr, 10);
          hv[0] += seed;
          if (! isdigit (*feat) || *endptr != '\0')
            {
              MurmurHash3_x64_128 (feat, strlen (feat), seed, hv);
            }

          float xi = hv[0] & (1ULL << 63) ? 1 : -1;

          unsigned int n = static_cast<unsigned char> (ns[0]);

          f[n].push_back (Feature (hv[0], xi * value * nsvalue));
        }
    }

  return f;
}

std::ostream&
operator<< (std::ostream&               o,
            const IndexedFeatureSet&    f)
{
  for (unsigned int n = 0; n < 256; ++n)
    {
      if (f[n].size ())
        {
          o << " (" << n << ") " << f[n] << " ";
        }
    }

  o << std::endl;

  return o;
}

bool
fuzzy_eq (const std::vector<Feature>& a,
          const std::vector<Feature>& b)
{
  if (a.size () != b.size ())
    {
      std::cerr << "fuzzy_eq feature set size mismatch " << a.size () << " ?= " << b.size () << std::endl;
      return false;
    }

  for (std::vector<Feature>::const_iterator x = a.begin (), y = b.begin ();
       x != a.end () && y != b.end ();
       ++x, ++y)
    {
      if (x->id != y->id)
        {
          std::cerr << "fuzzy_eq feature difference " << x->id << " ?= " << y->id << std::endl;
          return false;
        }

      float x_xi = x->id & (1ULL << 63) ? 1 : -1;
      float y_xi = y->id & (1ULL << 63) ? 1 : -1;

      if (fabs (x_xi * x->value - y_xi * y->value) 
          > 1e-4 * (1.0 + fabs (x->value) + fabs (y->value)))
        {
          std::cerr << "fuzzy_eq value difference " << x->value << " ?= " << y->value << std::endl;
          return false;
        }
    }

  return true;
}

bool
fuzzy_eq (const IndexedFeatureSet& a,
          const std::vector<Feature>*const b)
{
  for (unsigned int n = 0; n < 256; ++n)
    {
      if (! fuzzy_eq (a[n], b[n]))
        {
          return false;
        }
    }

  return true;
}

void
test_parse_once (void)
{
  GeneralExample result;
  result.label = (drand48 () < 0.3) ? 0 : random_string (not_space);
  result.importance = (drand48 () < 0.3) ? 0 : random_string (not_space_or_bar);
  result.tag = (drand48 () < 0.3) ? 0 : random_string (not_space_or_bar);

  AllFeatures af = random_features ();

  std::stringstream ss;

  if (result.label) { ss << result.label; }
  ss << " ";
  if (result.importance) { ss << result.importance << " "; }
  if (result.tag) { ss << result.tag; } 
  ss << af;
  ss << std::endl;

//  std::cout << "test line is " << ss.str (); 

  ss.seekg (0, std::ios_base::beg);

  StandardParse<HashString> parse;

  parse (ss, 
         [ af, result, & ss ] (bool ok, const GeneralExample& ex) {
    assert (ok);
    assert ((! ex.label && ! result.label) ||
            (ex.label && result.label && 
             strcmp (ex.label, result.label) == 0) || 
            ((std::cerr << "generated label is '" << (result.label ? result.label : "")
                        << "' (" << static_cast<void*> (result.label) << ") "
                        << " parsed label is '" << (ex.label ? ex.label : "")
                        << "' (" << static_cast<void*> (ex.label) << ") "
                        << std::endl), false));
  
    assert ((! ex.importance && ! result.importance) ||
            (ex.importance && result.importance && 
             strcmp (ex.importance, result.importance) == 0) || 
            ((std::cerr << "generated importance is '" << (result.importance ? result.importance : "")
                        << "' (" << static_cast<void*> (result.importance) << ") "
                        << " parsed importance is '" << (ex.importance ? ex.importance : "")
                        << "' (" << static_cast<void*> (ex.importance) << ") "
                        << std::endl), false));
    assert ((! ex.tag && ! result.tag) ||
            (ex.tag && result.tag && 
             strcmp (ex.tag, result.tag) == 0) || 
            ((std::cerr << "generated tag is '" << (result.tag ? result.tag : "")
                        << "' (" << static_cast<void*> (result.tag) << ") "
                        << " parsed tag is '" << (ex.tag ? ex.tag : "")
                        << "' (" << static_cast<void*> (ex.tag) << ") "
                        << std::endl), false));
  
    assert (ex.f);
  
    IndexedFeatureSet ifs = get_all_features (af);
  
    assert (fuzzy_eq (ifs, ex.f) || ((std::cerr << "ifs       = " << ifs << std::endl << "example.f = " << ex.f << std::endl << "test line is " << ss.str ()), false));
  });

  delete[] result.label; 
  delete[] result.importance;
  delete[] result.tag;
}

void
test_parse (void)
{
  for (size_t i = 0; i < 10000; ++i)
    {
      test_parse_once ();
    }
}

}

int 
main (void)
{
  srand48 (69);
  test_parse ();

//  std::cerr << "all tests passed" << std::endl;

  return 0;
}
