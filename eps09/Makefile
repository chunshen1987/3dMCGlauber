#=================================================================================
#
#  File: Makefile
#  Author: Thomas Ullrich (thomas.ullrich@bnl.gov)
#  Last modified: June 16, 2009
#
#  Compile and link eps09query, a command line tool to EPS09 
#
#=================================================================================
CCX = g++
CCXOPT = -O -W -Wall

eps09query: eps09.o eps09.h eps09query.o
	$(CCX) $(CCXOPT) eps09.o eps09query.o -o eps09query

eps09.o: eps09.cxx eps09.h
	$(CCX) $(CCXOPT) -c eps09.cxx

eps09query.o: eps09query.cxx eps09.h
	$(CCX) $(CCXOPT) -c eps09query.cxx

clean:
	rm -rf *.o eps09query eps09query.dSYM

