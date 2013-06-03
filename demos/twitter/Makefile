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
	time ../../pca -a -b 1299709 -c -k 30 -m $@ <(zcat $<) <(zcat $<)

embedding.out: pca.components twitter.pca.gz ../../pca.check zcat.check
	@echo "component scores for @biz (Twitter id 13):"
	@../../pca -t -m $(word 1,$^) <(zcat $(word 2,$^)) | head -1

allusers: twitter_rv.zip perl.check unzip.check
	unzip -p $<							\
	| perl -lane 'do { print "$$l\t$$n"; $$n = 0; } if 		\
	                   1 && $$l && $$l != $$F[0];			\
		      1; ++$$n; $$l = $$F[0]; } { print "$$l\t$$n";'	\
	| sort -k2,2rn -S30% > $@

topusers.pca.out: allusers twitter.pca.gz pca.components ../../pca.check perl.check
	../../pca -t -s -m $(word 3, $^)				\
	  <(perl -lane 'print " $$F[0]|0 $$F[0]";' $(word 1,$^))	\
	| head -10000 > $@

data.dat: topusers.pca.out topuser2tsne perl.check
	./topuser2tsne <(head -2000 $<) > $@

result.dat: data.dat tSNE_linux.check
	tSNE_linux $< 

coords: result.dat topusers.pca.out
	@./jointsne							\
	  <(./tsneresultdecode < $(word 1,$^)) 				\
	  <(cut -f1 -d' ' $(word 2,$^)) |				\
	env NO_LOOKUP=1 ./join-twitter twitter.cache |			\
	grep -v -e '_unknown_' -e '_notfamous_' |			\
	./diverse-coords | 						\
	( echo -n "coords = {";						\
	  perl -ane 'print "," if $$first++; 				\
	             1; print " { \"$$F[0]\", $$F[1], $$F[2], $$F[3], $$F[4] }"';	\
	  echo " };" ) > $@
