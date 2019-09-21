# --- SYSTEM ---

#linux
#SYSTEM = x86-64_linux
#mac
SYSTEM = x86-64_osx
LIBFORMAT  = static_pic

# --- DIRECTORIES ---

CCC = g++

BASISILOG = /Users/margarita/Applications/IBM/ILOG/CPLEX_Studio127
CPOPTDIR   = $(BASISILOG)/cpoptimizer
CONCERTDIR = $(BASISILOG)/concert
CPLEXDIR   = $(BASISILOG)/cplex
BOOSTDIR   = ./include/boost/include 

# --- FLAGS ---

CCOPT = -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -Wno-long-long -ffloat-store -m64 -DILOUSEMT -D_REENTRANT -DILM_REENTRANT
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPOPTLIBDIR = $(CPOPTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CPOPTINCDIR   = $(CPOPTDIR)/include

# --- OPTIMIZATION FLAGS ---

DEBUG_OPT = -DNDEBUG -O3 -g
#DEBUG_OPT = -g3 -O0
#PROF = -pg
#PROF =

CFLAGS = -DIL_STD $(DEBUG_OPT)  -I$(CPOPTDIR)/include -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -Wno-long-long -ffloat-store -m64 -DILOUSEMT -D_REENTRANT -DILM_REENTRANT -I$(BOOSTDIR) -c

LDFLAGS = -m64 -L$(CPOPTLIBDIR) -lcp -L$(CPLEXLIBDIR)  -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert  -lpthread -lm 

# ---- COMPILE  ----
SRC_DIR   := src
OBJ_DIR   := obj

SRC_DIRS  := $(shell find $(SRC_DIR) -type d)
OBJ_DIRS  := $(addprefix $(OBJ_DIR)/,$(SRC_DIRS))

SOURCES   := $(shell find $(SRC_DIR) -name '*.cpp')
OBJ_FILES := $(addprefix $(OBJ_DIR)/, $(SOURCES:.cpp=.o))

vpath %.cpp $(SRC_DIRS)

# ---- TARGETS ----

EXECUTABLE = mdd_mPDTSP

all: $(EXECUTABLE)

$(EXECUTABLE): makedir $(SOURCES) $(OBJ_FILES) 
	$(CCC) $(OBJ_FILES) $(LDFLAGS) $(PROF) -o $@

$(OBJ_DIR)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@

makedir: $(OBJ_DIRS)

$(OBJ_DIRS):
	@mkdir -p $@

clean:
	@rm -rf obj 
	@rm -rf $(EXECUTABLE)
