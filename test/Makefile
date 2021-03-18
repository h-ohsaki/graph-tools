all:
	@for i in 0*.py 40*.py [8-9]*.py; \
	do \
	  echo $$i:; \
	  python3 $$i; \
	  echo; \
	done | tee log

	@npass=`grep ' ok ' log | wc -l`; \
	nfail=`grep ' ng ' log | wc -l`; \
	echo "$$npass passed, $$nfail failed"
