#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= 
SYSLIBS         =
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}


ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR) -lAnitaEvent -lAnitaCorrelator
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
ANITA_UTIL_CALIB_DIR=$(ANITA_UTIL_INSTALL_DIR)/share/anitaCalib
else
ANITA_UTIL_LIB_DIR=/usr/local/lib
ANITA_UTIL_INC_DIR=/usr/local/include
ANITA_UTIL_CALIB_DIR=/usr/local/share/anitaCalib
ifdef EVENT_READER_DIR
LD_ANITA_UTIL=-L$(EVENT_READER_DIR)  -lAnitaEvent
INC_ANITA_UTIL=-I$(EVENT_READER_DIR)
endif
endif

#Toggles the FFT functions on and off
USE_FFT_TOOLS=1

ifdef USE_FFT_TOOLS
FFTLIBS = -L/usr/local/lib -lRootFftwWrapper -lfftw3
FFTFLAG = -DUSE_FFT_TOOLS
else
FFTLIBS =
FFTFLAG =
endif

#Generic and Site Specific Flags
CXXFLAGS     = -std=c++11 $(ROOTCFLAGS) $(FFTFLAG) $(SYSINCLUDES) $(INC_ANITA_UTIL) 
LDFLAGS      = -g $(ROOTLDFLAGS) 


LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) $(LD_ANITA_UTIL) $(FFTLIBS) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Toggles google performance profile functionality on and off
USE_GPERFTOOLS=1

ifdef USE_GPERFTOOLS
LDFLAGS	+= #-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
LIBS += -lprofiler #-ltcmalloc
endif


HARDCORE_MODE=1
ifdef HARDCORE_MODE
CXXFLAGS += -Wextra -Wextra -Wshadow -Werror -Wpedantic
endif




#ROOT stuff
ROOT_LIBRARY = libCorrelator.${DLLSUF}
LIB_OBJS = CrossCorrelator.o FancyTTreeInterpolator.o 
CLASS_HEADERS = CrossCorrelator.h FancyTTreeInterpolator.h

#Now the bits we're actually compiling
all: $(ROOT_LIBRARY) testCorrelator testFancyTTreeInterpolator phaseCentre waisPulseReconstruction prioritizerdPerformance distanceFromWaisDivide commit

testCorrelator: $(ROOT_LIBRARY) testCorrelator.$(SRCSUF)
	@echo "<**Compiling**> "
	$(LD) $(CXXFLAGS) $(LDFLAGS) testCorrelator.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
	@if test $$? == 0; then git add testCorrelator.$(SRCSUF); fi

testFancyTTreeInterpolator: $(ROOT_LIBRARY) testFancyTTreeInterpolator.$(SRCSUF)
	@echo "<**Compiling**> "
	$(LD) $(CXXFLAGS) $(LDFLAGS) testFancyTTreeInterpolator.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
	@if test $$? == 0; then git add testFancyTTreeInterpolator.$(SRCSUF); fi

phaseCentre: $(ROOT_LIBRARY) phaseCentre.$(SRCSUF)
	@echo "<**Compiling**> "
	$(LD) $(CXXFLAGS) $(LDFLAGS) phaseCentre.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
	@if test $$? == 0; then git add phaseCentre.$(SRCSUF); fi

waisPulseReconstruction: $(ROOT_LIBRARY) waisPulseReconstruction.$(SRCSUF)
	@echo "<**Compiling**> "
	$(LD) $(CXXFLAGS) $(LDFLAGS) waisPulseReconstruction.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
	@if test $$? == 0; then git add waisPulseReconstruction.$(SRCSUF); fi

distanceFromWaisDivide: $(ROOT_LIBRARY) distanceFromWaisDivide.$(SRCSUF)
	@echo "<**Compiling**> "
	$(LD) $(CXXFLAGS) $(LDFLAGS) distanceFromWaisDivide.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
	@if test $$? == 0; then git add distanceFromWaisDivide.$(SRCSUF); fi


prioritizerdPerformance: $(ROOT_LIBRARY) prioritizerdPerformance.$(SRCSUF)
	@echo "<**Compiling**> "
	$(LD) $(CXXFLAGS) $(LDFLAGS) prioritizerdPerformance.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
	@if test $$? == 0; then git add prioritizerdPerformance.$(SRCSUF); fi


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .
	$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	$(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF) %.h
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@if test $$? == 0; then git add $^; fi

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@
	@if test $$? == 0; then git add $%.C; fi

correlatorDict.C : $(CLASS_HEADERS)
		@echo "<**And here's the dictionary...**>" $<
		@rm -f *Dict*
#		rootcint $@ -c $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h
		rootcint $@ -c $(INC_ANITA_UTIL) FancyTTreeInterpolator.h LinkDef.h

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f testCorrelator testFancyTTreeInterpolator phaseCentre waisPulseReconstruction prioritizerdPerformance distanceFromWaisDivide

commit: 
	@git add Makefile
	@git commit
