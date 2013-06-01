SHELL=/bin/bash

.SECONDARY:

all:
	@cat README.md

demo: embedding.out

clean:

%.check:
	@test -x "$$(which $*)" || {				\
	  echo "ERROR: you need to install $*" 1>&2;		\
	  exit 1;						\
	}

TWITTER2SVD='do { print " $$l|0 ", join " ", @ids    		\
	    	    if $$\#ids >= 6 && $$\#ids < 300; 		\
		  @ids = (); } 					\
	       if $$l && $$l != $$F[1];				\
	     push @ids, $$F[0];					\
	     $$l = $$F[1]; } {					\
	       do { print " $$l|0 ", join " ", @ids    		\
		      if $$\#ids >= 6 && $$\#ids < 300;		\
		    @ids = (); }'

twitter_rv.zip: wget.check
	wget http://an.kaist.ac.kr/~haewoon/release/twitter_social_graph/twitter_rv.zip

twitter.pca.gz: twitter_rv.zip perl.check unzip.check sort.check lzop.check gzip.check
	unzip -p $< | sort --compress-program lzop -S70% -t'	' -k2,2n -k1,1n \
	            | perl -lane $(TWITTER2SVD) 			\
		    | gzip > $@

pca.components: twitter.pca.gz ../../pca.check zcat.check
	time ../../pca -b 1299709 -c -k 30 -m $@ <(zcat $<) <(zcat $<)

embedding.out: pca.components twitter.pca.gz ../../pca.check zcat.check
	../../pca -t -m $(word 1,$^) <(zcat $(word 2,$^)) | head -1

