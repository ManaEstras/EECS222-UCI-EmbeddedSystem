# Makefile: Canny Edge Detector model in SpecC

# --- settings ---

VIDEO	= EngPlaza
FRAMES	= $(VIDEO)[0-9][0-9][0-9]_edges.pgm

SCC	= scc
SCCOPT	= -g -ww -vv

# --- targets ---

all:	Canny

test:	Canny
	./Canny
	set -e;	\
		for f in video/$(FRAMES); do \
		diff `basename $$f` $$f; \
		done

clean:
	rm -f *~ *.bak *.BAK
	rm -f *.o *.h *.cc *.si *.sir
	rm -f Canny
	rm -f $(FRAMES)

# --- compile the example ---

Canny: Canny.sc
	$(SCC) $@ $(SCCOPT) -i $<

# EOF

