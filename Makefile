#---------------------------------------------------------------------------------------------------
# Compiler selection
#---------------------------------------------------------------------------------------------------

MPICXX	= mpicxx
CXX = g++

#---------------------------------------------------------------------------------------------------
# Directories
#---------------------------------------------------------------------------------------------------

OBJDIR = build
SRCDIR = src

#---------------------------------------------------------------------------------------------------
# Executables
#---------------------------------------------------------------------------------------------------

EXE = sync-greedy

#---------------------------------------------------------------------------------------------------
# Object files
#---------------------------------------------------------------------------------------------------

SYNCOBJ			= ConfigParser.o ExprsData.o Parallel.o GreedyController.o GreedyWorker.o SolPool.o \
							Timer.o

#---------------------------------------------------------------------------------------------------
# Compiler options
#---------------------------------------------------------------------------------------------------

CXXFLAGS = -O3 -Wall -fPIC -fexceptions -DIL_STD -std=c++11 -fno-strict-aliasing

#---------------------------------------------------------------------------------------------------
# Link options and libraries
#---------------------------------------------------------------------------------------------------

OPENMPI		   	 = <OMPI_LIB_PATH>

CXXLNDIRS      = -L$(OPENMPI)
CXXLNFLAGS     = -lm -lpthread -ldl

#---------------------------------------------------------------------------------------------------

all: CXXFLAGS += -DNDEBUG
all: $(EXE)

debug: CXXFLAGS += -g
debug: $(EXE)


sync-greedy: $(OBJDIR)/main.o
	$(MPICXX) -o $@ $(addprefix $(OBJDIR)/, $(SYNCOBJ) main.o)

$(OBJDIR)/main.o:	$(addprefix $(SRCDIR)/, main.cpp) \
									$(addprefix $(OBJDIR)/, $(SYNCOBJ) ) 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/ConfigParser.o: $(addprefix $(SRCDIR)/, ConfigParser.cpp ConfigParser.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/ExprsData.o:	$(addprefix $(SRCDIR)/, ExprsData.cpp ExprsData.h) \
												$(addprefix $(OBJDIR)/, ConfigParser.o) 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Parallel.o:	$(addprefix $(SRCDIR)/, Parallel.cpp Parallel.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/GreedyController.o:	$(addprefix $(SRCDIR)/, GreedyController.cpp GreedyController.h) \
									$(addprefix $(OBJDIR)/, Parallel.o ExprsData.o ConfigParser.o SolPool.o Timer.o) 
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/GreedyWorker.o:	$(addprefix $(SRCDIR)/, GreedyWorker.cpp GreedyWorker.h) \
									$(addprefix $(OBJDIR)/, Parallel.o ExprsData.o ConfigParser.o SolPool.o) 
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/SolPool.o:	$(addprefix $(SRCDIR)/, SolPool.cpp SolPool.h Utils.h) \
											$(addprefix $(OBJDIR)/, ExprsData.o) 
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Timer.o:	$(addprefix $(SRCDIR)/, Timer.cpp Timer.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<


#---------------------------------------------------------------------------------------------------
.PHONY: clean cleanest
clean:
	/bin/rm -f $(OBJDIR)/*.o

cleanest:
	/bin/rm -f $(OBJDIR)/*.o *.log *.cuts *.lp $(EXE)