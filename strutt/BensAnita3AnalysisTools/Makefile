#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#If you have the 64-bit version of fftw installed, try this to help CINT out.
CINTFLAGS=-DFFTW_64_BIT

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
ANITA_UTIL_ETC_DIR=$(ANITA_UTIL_INSTALL_DIR)/etc
else
ANITA_UTIL_LIB_DIR=/usr/local/lib
ANITA_UTIL_INC_DIR=/usr/local/include
ANITA_UTIL_ETC_DIR=/usr/local/etc
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




# For the first time, I've come across something that changed across ROOT versions...
# So now we have to hack in this compile time macro.
ifeq ($(shell test $(shell root-config --version | cut -c1) -ge 6; echo $$?) ,  0) 
CXXFLAGS += -DIS_ROOT_6
endif


# For those who like really bloody pedantic compiler warnings... like me
HARDCORE_MODE=1
ifdef HARDCORE_MODE
CXXFLAGS += -Wall -Wextra -Wshadow -Werror #-Wpedantic
endif

#ROOT stuff
ROOT_LIBRARY = libBensAnitaTools.${DLLSUF}
DICT = benToolsDict
LIB_OBJS = $(DICT).o CrossCorrelator.o FancyTTreeInterpolator.o RootTools.o FancyFFTs.o ProgressBar.o OutputConvention.o AveragePowerSpectrum.o AntarcticaMapPlotter.o AnitaAveragePowerSpectrum.o
CLASS_HEADERS = CrossCorrelator.h FancyTTreeInterpolator.h FancyFFTs.h RootTools.h ProgressBar.h OutputConvention.h AveragePowerSpectrum.h AntarcticaMapPlotter.h AnitaAveragePowerSpectrum.h
BINARIES = testCorrelator testFancyTTreeInterpolator testFancyFFTs testDeltaTsSpherical testProgressBar testOutputConvention testAveragePowerSpectrum

#Now the bits we're actually compiling
all: $(ROOT_LIBRARY) $(BINARIES) 

.PHONY: install clean docs

$(BINARIES): %: %.$(SRCSUF) $(ROOT_LIBRARY) 
	@echo "<**Compiling**> "
	@echo $<
	$(LD) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $< $(ROOT_LIBRARY) -o $@

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
%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


$(DICT).C : $(CLASS_HEADERS)
		@echo "<**And here's the dictionary...**>" $<
		@rm -f *Dict*
#		rootcint $@ -c -p $(CXXFLAGS) $(CLASS_HEADERS) LinkDef.h
		rootcint $@ -c -p $(CINTFLAGS) $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(BINARIES) 

install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
#	install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
	-@install -c -m 755 $(DICT)_rdict.pcm $(ANITA_UTIL_LIB_DIR)

else
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	install -c -m 644  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)

	install -d $(ANITA_UTIL_ETC_DIR)

	install -d $(ANITA_UTIL_INSTALL_DIR)/share/anitaMap
	for file in anitaMap/*.png; do install -c -m 644 "$${file}" $(ANITA_UTIL_INSTALL_DIR)/share/anitaMap; done
