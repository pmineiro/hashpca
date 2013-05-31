all: check

clean:
	rm -f data s u v testmodel

check: testsvd.ok testpca.ok

%.ok: % 
	@echo "=== running test $(word 1,$^) ===" 1>&2
	@./$(word 1,$^)
	@echo "=== $(word 1,$^) passed ===" 1>&2
