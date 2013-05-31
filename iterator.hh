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

#ifndef __ITERATOR_HH__
#define __ITERATOR_HH__

#include <cstdlib>
#include <Eigen/Dense>
#include "parse.hh"

namespace hashpca 
{
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          MatrixXd;

  template<typename S,
           typename T>
  struct Times
    {
      S left;
      T right;

      Times (S _left, T _right) : left (_left), right (_right) { }
    };

  template<typename T>
  struct Transposed
    {
      T t;

      Transposed (T _t) : t (_t) { }

      template<typename Right>
      Times<Transposed, Right&>
      operator* (Right& r)
        {
          return Times<Transposed, Right&> (*this, r);
        }
    };

  template<typename T>
  Times<double, Transposed<T> >
  operator* (double             x,
             Transposed<T>      t)
    {
      return Times<double, Transposed<T> > (x, t);
    }

  class LinearIterator
    {
      private:
        uint64_t hashsize;

      public:
        LinearIterator (uint64_t _hashsize) : hashsize (_hashsize) { }

        struct Deferred
          {
            uint64_t hashsize;
            const veedubparse::GeneralExample& ex;

            Deferred (uint64_t                           _hashsize,
                      const veedubparse::GeneralExample& _ex)
              : hashsize (_hashsize),
                ex (_ex)
              {
              }

            template<typename F>
            void
            foreach (F f) const
              {
                for (unsigned int n = 0; n < 256; ++n)
                  {
                    const veedubparse::Feature* fptr = &ex.f[n][0];

                    for (unsigned int j = 0; j < ex.f[n].size (); ++j, ++fptr)
                      {
                        veedubparse::Feature x (fptr->id % hashsize, fptr->value);

                        f (x);
                      }
                  }
              }

            Transposed<Deferred>
            transpose ()
              {
                return Transposed<Deferred> (*this);
              }
          };
  
        Deferred
        operator() (const veedubparse::GeneralExample& ex)
          {
            return Deferred (hashsize, ex);
          }
    };

  Times<double, LinearIterator::Deferred>
  operator* (double                     x,
             LinearIterator::Deferred   d)
    {
      return Times<double, LinearIterator::Deferred> (x, d);
    }

  class BiLinearIterator
    {
      private:
        uint64_t        hashsize;
        unsigned char   a;
        unsigned char   b;

        static const uint64_t quadratic_constant = 27942141;

      public:
        BiLinearIterator (uint64_t      _hashsize,
                          unsigned char _a,
                          unsigned char _b)
          : hashsize (_hashsize),
            a (_a),
            b (_b)
          {
          }

        struct Deferred
          {
            uint64_t hashsize;
            unsigned char a;
            unsigned char b;
            const veedubparse::GeneralExample& ex;

            Deferred (uint64_t                           _hashsize,
                      unsigned char                      _a,
                      unsigned char                      _b,
                      const veedubparse::GeneralExample& _ex)
              : hashsize (_hashsize),
                a (_a),
                b (_b),
                ex (_ex)
              {
              }

            template<typename F>
            void
            foreach (F f) const
              {
                const veedubparse::Feature* iptr = &ex.f[a][0];

                for (unsigned int i = 0; i < ex.f[a].size (); ++i, ++iptr)
                  {
                    const veedubparse::Feature* jptr = &ex.f[b][0];

                    for (unsigned int j = 0; j < ex.f[b].size (); ++j, ++jptr)
                      {
                        veedubparse::Feature 
                          x ((quadratic_constant * iptr->id + jptr->id) % hashsize,
                             iptr->value * jptr->value);

                        f (x);
                      }
                  }
              }

            Transposed<Deferred>
            transpose ()
              {
                return Transposed<Deferred> (*this);
              }
          };
  
        Deferred
        operator() (const veedubparse::GeneralExample& ex)
          {
            return Deferred (hashsize, a, b, ex);
          }
    };

  Times<double, BiLinearIterator::Deferred>
  operator* (double                     x,
             BiLinearIterator::Deferred d)
    {
      return Times<double, BiLinearIterator::Deferred> (x, d);
    }

  struct SugaryVectorXd : public Eigen::VectorXd
    {
      SugaryVectorXd (unsigned int cols) : Eigen::VectorXd (cols) { };

      template<typename Iterator>
      SugaryVectorXd&
      operator= (Times<Transposed<Iterator>, const hashpca::MatrixXd&> op)
        {
          setZero (size ());

          op.left.t.foreach ([ &op, this ] (const veedubparse::Feature& f) {
            *this += f.value * op.right.row (f.id);
          });

          return *this;
        }
    };

  template<typename T>
  Times<Times<double, T>, SugaryVectorXd&>
  operator* (Times<double, T>  left,
             SugaryVectorXd&   right)
    {
      return Times<Times<double, T>, SugaryVectorXd&> (left, right);
    }

  template<typename Iterator>
  hashpca::MatrixXd&
  operator+= (hashpca::MatrixXd&                                x,
              Times<Times<double, Iterator>, SugaryVectorXd&>   op)
    {
      op.left.right.foreach ([&] (const veedubparse::Feature& f) {
        x.row (f.id) += op.left.left * f.value * op.right;
      });

      return x;
    }

  template<typename Iterator>
  Eigen::VectorXd&
  operator+= (Eigen::VectorXd&                          x,
              Times<double, Transposed<Iterator> >      op)
    {
      op.right.t.foreach ([&] (const veedubparse::Feature& f) {
        x (f.id) += op.left * f.value;
      });

      return x;
    }
}

#endif // __ITERATOR_HH__
