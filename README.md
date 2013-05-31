hashpca
=======

Scalable PCA via Hashing

Description
-----------

This software provides a highly scalable SVD/PCA implementation using
a combination of hashing and randomized linear algebra techniques.

The memory usage is independent of the number of examples so you will
never run out of memory.  

Directions
-----------

0. The basic workflow looks like:
	1. Build a model

        > pca -c -k 100 -m testmodel inputdata
	2. Use a model

        > pca -t -m testmodel inputdata  
1. The software expects input in [vowpal wabbit format](https://github.com/JohnLangford/vowpal_wabbit/wiki/Input-format).  
	1. The label is not used but is required, so basically start all your lines with whitespace.
	2. The importance weight is used but if omitted defaults to 1.
	3. The tag if present will be passed through to the output in projection mode.  (It is ignored when constructing a model).
2. The software will do 2 streaming passes over the data to build a model.  Here are two techniques you can use to avoid intermediate materialization of the input data: 
	1.  You can supply two arguments when building a model and the software will open the first one for the first pass and the second one for the second pass.  This is  useful when composed with [zsh process substitution](http://zsh.sourceforge.net/Intro/intro_7.html).  For example, with compressed data:

        >  pca -k 100 -m testmodel <(zcat compressed.gz) <(zcat compressed.gz)
	2. If you supply one argument when building a model the software will open that argument twice and access the data from each opened version sequentially without seeking.  This is useful when composed with [named pipes](https://en.wikipedia.org/wiki/Named_pipe).
