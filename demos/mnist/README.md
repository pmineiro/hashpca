mnist pca demo
--------------

In this demo PCA is applied recursively to build a low-rank primal 
approximation to a (3rd-order) polynomial kernel.  This is used 
with l1-regularized logistic regression on (permutation-invariant) mnist
yielding 176 test errors with vowpal wabbit 7.2.  liblinear can achieve 150 
test errors with the same design matrix, but it requires >64G of ram so the 
demo here uses vowpal wabbit.

To execute the demo, type `make demo`.  This will:

  * Download the mnist data sets, which are fairly modest (circa 10Mb).
  * PCA the training set down to 40 dimensions.
  * Interact the 40 dimensional projection with the original pixels and 
    PCA the result down to 50 dimensions.
  * Train an L1-regularized logistic regression with interaction features
    between raw pixels and the 50 dimensional second stage projection.
  * Test the resulting model and report the test errors and confusion matrix.
