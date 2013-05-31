.SECONDARY:

CXXFLAGS=-O3 -Wall
LABELS=train-labels-idx1-ubyte.gz
IMAGES=train-images-idx3-ubyte.gz
TESTLABELS=t10k-labels-idx1-ubyte.gz
TESTIMAGES=t10k-images-idx3-ubyte.gz

all: 
	@cat README.md

demo: svd.test.confusion

clean:
	rm -f $(wildcard svd.feats sqsvd.feats extractfeatures *.cache* *.predictions)

%.check:
	@test -x "$$(which $*)" || {				\
	  echo "ERROR: you need to install $*" 1>&2;		\
	  exit 1;						\
	}

train-%.gz: wget.check
	wget http://yann.lecun.com/exdb/mnist/train-$*.gz

t10k-%.gz: wget.check
	wget http://yann.lecun.com/exdb/mnist/t10k-$*.gz

extractfeatures: extractfeatures.cpp

svd.feats: $(LABELS) $(IMAGES) extractfeatures ../../pca.check
	@./do-svd-feats $(word 1,$^) $(word 2,$^) $@

svd.model: svd.feats $(LABELS) $(IMAGES) zsh.check perl.check paste.check vw.check
	@./do-train $(word 1,$^) $(word 2,$^) $(word 3,$^) $@

svd.test.predictions: svd.model $(TESTLABELS) $(TESTIMAGES)
	@./do-test svd $(word 2,$^) $(word 3,$^) $@

CONFUSION='++$$n;                                                       \
           ++$$c if $$F[0] != $$F[1];                                   \
           ++$$m{"@{[int($$F[0])]}:$$F[1]"}; } {                        \
           $$c ||= 0;                                                   \
           print "$* test errors: $$c out of $$n";                      \
           print "confusion matrix (rows = truth, columns = prediction):"; \
           foreach $$true (1 .. 10) {                                   \
             print join "\t", map { $$m{"$$true:$$_"} || 0 } (1 .. 10); \
           }'

%.confusion: %.predictions 
	@perl -lane $(CONFUSION) $*.predictions
