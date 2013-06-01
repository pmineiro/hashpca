mnist pca demo
--------------

__This demo is very disk intensive__.

In this demo PCA is applied to a matrix defined by the following
relation on Twitter.  This results in an embedding of any Twitter 
account into a 30 dimensional space based upon the set of users
followed by the account.  

To execute the demo, type `make demo`.  This will:

  * Download the __4.5 Gigabytes__ compressed Twitter graph.
  * Transpose the graph and convert it to vowpal wabbit format.
        * This involves an expensive invocation of sort,
          which uses about 20 Gb of temporary disk space,
          and takes about an hour on a fast desktop machine.
  * PCA the output of the previous step.
        * As is typical in machine learning, this is faster
          than the time it takes to prepare the input format.
