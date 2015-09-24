#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= -I/usr/local/include
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
CXXFLAGS     = -g -fPIC $(ROOTCFLAGS) $(FFTFLAG) $(SYSINCLUDES) $(INC_ANITA_UTIL) #-std=c++11 
LDFLAGS      = -g $(ROOTLDFLAGS) 

LIBS          = $(ROOTLIBS) -lMathMore -lMinuit -lGraf $(SYSLIBS) $(LD_ANITA_UTIL) $(FFTLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Toggles google performance profile functionality on and off
#USE_GPERFTOOLS=1

ifdef USE_GPERFTOOLS
LDFLAGS	+= -Wl,-no_pie -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
LIBS += -lprofiler -ltcmalloc
endif


# I use git svn and am too infrequent a commiter by nature...
# Use this flag to prompt a local commit after every compile
# (You probably don't actually want to commit after every compile, but the prompt is nice)
#FORCE_GIT=1 


# For those who like really bloody pedantic compiler warnings... like me
HARDCORE_MODE=1
ifdef HARDCORE_MODE
CXXFLAGS += -Wall -Wextra -Wshadow -Werror #-Wpedantic
endif


#ROOT stuff
ROOT_LIBRARY = libBensAnitaTools.${DLLSUF}
DICT = benToolsDict
LIB_OBJS = $(DICT).o CrossCorrelator.o FancyTTreeInterpolator.o RootTools.o FancyFFTsWisdomManager.o FancyFFTs.o ProgressBar.o
CLASS_HEADERS = CrossCorrelator.h FancyTTreeInterpolator.h FancyFFTsWisdomManager.h FancyFFTs.h RootTools.h ProgressBar.h
BINARIES = testCorrelator testFancyTTreeInterpolator testFancyFFTs testDeltaTsSpherical testProgressBar

#Now the bits we're actually compiling
all: $(ROOT_LIBRARY) $(BINARIES) commit

.PHONY: install commit clean docs

$(BINARIES): %: %.$(SRCSUF) $(ROOT_LIBRARY) 
	@echo "<**Compiling**> "
	@echo $<
	$(LD) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $< $(ROOT_LIBRARY) -o $@
ifdef FORCE_GIT
	-@if test $$? == 0; then git add $<; fi
endif

docs: Doxyfile
	doxygen Doxyfile

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .
	$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ 
	$(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF) %.h
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o $@
ifdef FORCE_GIT	
	@if test $$? == 0; then git add $^; fi
endif
%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


$(DICT).C : $(CLASS_HEADERS)
		@echo "<**And here's the dictionary...**>" $<
		@rm -f *Dict*
#		rootcint $@ -c -p $(CXXFLAGS) $(CLASS_HEADERS) LinkDef.h
		rootcint $@ -c -p $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(BINARIES) 

commit: 
ifdef FORCE_GIT
	-@git add Makefile
	-@git commit
endif

install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
#	install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
	-@install -c -m 755 $(DICT)_rdict.pcm $(ANITA_UTIL_LIB_DIR)

else
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	install -c -m 644  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)


