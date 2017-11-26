# ------------------------------------------------------------------------
# Makefile for simple examples
# ------------------------------------------------------------------------
#
# Modifications: (most recent first)
#
# RD 09/24/14	added Adder2.sc example
# RD 09/28/01	added Exception.sc example
# RD 05/18/01	added Handshaking1 and Handshaking2 examples
# RD 04/05/01	added Pipeline2.sc example with new 'pipe' syntax
# RD 03/23/99	added Pipeline.sc example
# RD 12/30/98	added Timing.sc example
# RD 10/01/98	added -g option so that examples are easy debuggable
# RD 09/14/98	made "make test" sensitive to failures
# RD 09/08/98	added Adder.sc example
# RD 06/02/98	first version


# --- macros -------------------------------------------------------------

ALL	= Adder Behaviors BitVectors Callback DataTypes FSM \
	  HelloWorld Notes Pipeline Pipeline2 Timing \
	  Handshaking1 Handshaking2 Exception Adder2

SCC	= scc
SCCOPT	= -vv -w -g -d
RM	= rm -f


# --- SpecC rules --------------------------------------------------------

.SUFFIXES:
.SUFFIXES:	.sc .cc .o

.sc.cc:
	$(SCC) $* -sc2cch $(SCCOPT)

.cc.o:
	$(SCC) $* -cc2o $(SCCOPT)

.o:
	$(SCC) $* -o2out $(SCCOPT)

.cc:
	$(SCC) $* -cc2out $(SCCOPT)

.sc:
	$(SCC) $* -sc2out $(SCCOPT)


# --- targets ------------------------------------------------------------

all:	$(ALL)

clean:
	-$(RM) *.bak *.BAK
	-$(RM) *.si *.sir *.cc *.h *.o
	-$(RM) $(ALL) *.exe

test:	$(ALL)
	set -e;					\
	for file in $(ALL); do ./$$file ; done


# --- EOF ----------------------------------------------------------------
