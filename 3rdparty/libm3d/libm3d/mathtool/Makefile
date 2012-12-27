# $Id: Makefile,v 1.1 2002/03/15 04:22:04 neilien Exp $
#------------------------------------------------------
# Make your own and use it or use the one in BASELINE
#------------------------------------------------------

#-----------------------------------------------------------------------------
# Choose a compiler & its options
#--------------------------------------------------------------------------
 CC = g++

 CXX = g++

 OPTS = -O3

 INCLUDE = -I/usr/local/include

 CFLAGS     = $(OPTS) $(INCLUDE) -I.
 CXXFLAGS   = $(CFLAGS) 

all: nprmUtil

#--------------------------------------------------------------------
.SUFFIXES: .cpp
.cpp.o:
	$(CC) ${CFLAGS} -c $<

% : %.o
	$(CXX) $(CXXFLAGS) -o $@ $@.cpp
#--------------------------------------------------------------------
OBJS         =  EulerAngle.o  Gauss.o  Quaternion.o Basic.o

objs: ${OBJS}
#--------------------------------------------------------------------
nprmUtil: objs
	ar cru libmathtool.a $(OBJS)
clean:
	rm -f *.o libmathtool.a
