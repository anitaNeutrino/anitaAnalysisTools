#ifndef ANALYSIS_RECO_H
#define ANALYSIS_RECO_H

#include "AnitaEventReconstructor.h"
#include "InterferometricMap.h"
#include "InterferometryCache.h"
#include "CrossCorrelator.h"
#include "AnalysisSettings.h"
#include <list>

namespace AnitaResponse {
  class DeconvolutionMethod;
  class ResponseManager;
}
class NoiseMonitor;
class TruthAnitaEvent;

namespace Acclaim
{

  class InterferometryCache;
  
  /**
   * @brief Does the event reconstruction, and produces a summary of it
   *
   * The main workhorse function is process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL, TruthAnitaEvent* truth = NULL).
   * This takes in a FilteredAnitaEvent, produced a bunch of InterferometricMap objects, finds their peak directions, creates a coherently averaged waveform
   * (and deconvolved version) for each, and produces some summary numbers from the maps and waveforms and stores them in an AnitaEventSummary.
   */
  class AnalysisReco : public TObject {

  public:
    friend class InterferometryCache; ///< For accessing the antenna position arrays #fRArray, #fPhiArrayDeg, #fZArray

    enum DrawDomain{
      kTimeDomain = 0,
      kFreqDomain = 1
    };

    /**
     * @brief Constructor
     *
     * There are no constructor paramters, but a bunch of things are configurable with config files.
     * @see AnalysisSettings
     */
    AnalysisReco();

    /**
     * @brief Destructor
     */
    virtual ~AnalysisReco();

    /**
     * Reconstructs the FilteredAnitaEvent and fills all the reconstruction relative stuff in the passed AnitaEventSummary
     *
     * @param ev is the FilteredAnitaEvent to reconstruct
     * @param summary is the AnitaEventSummary to fill during the recontruction
     * @param noiseMonitor is an optional parameter, which points to a NoiseMonitor, which tracks the RMS of MinBias events (used in SNR calculation)
     * @param truth is an optional pointer to the TruthAnitaEvent class generated with MC events
     */
    void process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL, TruthAnitaEvent* truth = NULL);

    /**
     * @brief Get the expected delay between antenna pairs for a given direction
     *
     * @param pol is the polarisation
     * @param ant1 is the first antenna
     * @param ant2 is the second antenna
     * @param phiWave is the incoming plane wave direction in radians in payload coordinate relative to ADU5 aft-fore
     * @param thetaWave is the incoming plane wave direction in radians (theta=0 is horizontal, +ve theta is up)
     *
     * @return the expected time difference in nano-seconds
     */
    Double_t getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave) const;

    /**
     * @brief Get the off axis delay between two antennas for a given phi angle
     *
     * @param pol is the polarisation
     * @param ant1 is the first antenna
     * @param ant2 is the second antennas
     * @param phiDeg is the position in degrees in payload coordinates relative to ADU5 aft-fore
     *
     * @return the relative off axis delay
     */
    Double_t relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiDeg) const;

    /**
     * Get the off boresight delay for a single antenna
     *
     * @param deltaPhiDeg is the phi angle from boresight
     *
     * @return the delay
     */
    Double_t singleAntennaOffAxisDelay(Double_t deltaPhiDeg) const;

    /**
     * @brief Draw everything interesting onto a TPad.
     *
     * Having this here greatly simplifies the interface with MagicDisplay
     *
     * @param pad is the pad do draw on if it already exists (makes a new Canvas if passed NULL)
     * @param pol is the polarization to draw
     */
    void drawSummary(TPad* pad, AnitaPol::AnitaPol_t pol);

    /**
     * @brief Get a pointer to the coarsely binned interferometric map stored in memory, once called, you own this InterferometricMap and must delete it.
     *
     * The InterferometricMap is from the last event passed to process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL, TruthAnitaEvent* truth = NULL)
     * For speed, once this function is called the internal pointer is set to NULL and responsibility for deletion is given to the function caller.
     *
     * @param pol is the polarisation of the map to get.
     *
     * @return the InterferometricMap pointer previously stored in coarseMaps[AnitaPol::kNotAPol]
     */
    InterferometricMap* getMap(AnitaPol::AnitaPol_t pol);

    /**
     * @brief Get a pointer to a finely binned interferometric map stored in memory, once called, you own this InterferometricMap and must delete it.
     *
     * The InterferometricMap is from the last event passed to process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL, TruthAnitaEvent* truth = NULL)
     * For speed, once this function is called the internal pointer is set to NULL and responsibility for deletion is given to the function caller.
     *
     * @param pol is the polarisation of the map to get.
     * @param peakInd is the index of the finely binned peak, corresponding to the local maximum in the coarsely binned map (see )
     *
     * @return the InterferometricMap pointer previously stored in #fineMaps[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol]
     */
    InterferometricMap* getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd=0);

    /**
     * @brief produce and store the coarsely binned interferometric map, (overwrites #coarseMaps)
     *
     * @param pol is the polarisation
     * @param pat is a pointer to the GPS data
     */
    void reconstruct(AnitaPol::AnitaPol_t pol, const Adu5Pat* pat = NULL);

    /**
     * @brief produce and store a finely binned interferometric map, (overwrites #fineMaps)
     *
     * @param pol is the polarisation
     * @param peakIndex is the index of the coarsely binned maximum (from 0 to #fNumPeaks-1)
     * @param zoomCenterPhiDeg is the centre of the histogram in phi (azimuth, degrees)
     * @param zoomCenterThetaDeg is the centre of the histogram in theta (elevation, degrees)
     * @param pat is a pointer to the GPS data for the eventf
     */
    void reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg, const Adu5Pat* pat = NULL);


    /**
     * @brief Make a coherently summed waveform for a given polarization, set of antennas in a particular direction
     *
     * @param fEv is the FilteredAnitaEvent from which to draw the waveforms to coherently sum
     * @param pol is the polarisation of the waveforms in the FilteredAnitaEvent
     * @param theAnts is the set of antennas to use in the coherent sum
     * @param peakPhiDeg is the incoming phi-direction in degrees in payload coordinates relative to ADU5 aft-fore
     * @param peakThetaDeg is the elevation in degrees (theta=0 is horizontal, +ve theta is up) in payload coordinates relative to ADU5 aft-fore
     * @param biggestPeakToPeak if non-NULL stores the largest peak-to-peak in the set of channels passed to the function
     *
     * @return pointer to a newly created AnalysisWaveform, produces by coherently summing the parts
     */
    AnalysisWaveform* coherentlySum(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol, const std::vector<Int_t>& theAnts, Double_t peakPhiDeg, Double_t peakThetaDeg, Double_t* biggestPeakToPeak=NULL);


    /**
     * @brief Coherenty sum a set of AnalysisWaveforms with a set of dts.
     *
     * @param waves is a set of waveforms to coherently sum
     * @param dts is the set of time offsets (ns) to apply to those waveforms when summing
     *
     * @return a pointer to a newly created AnalysisWaveform, produced by coherently averaging the input waveforms
     */
    AnalysisWaveform* coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts);



    /** 
     * Directly insert some geometry generated by Linda.
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @warning This class overwrites the "extra cable delays" in AnitaEventCalibrator with zero! 
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * Probably don't use this, or if you do, use with caution.
     * 
     * @param pathToLindasFile points to the file containing the antenna positions/cable delays.
     * @param pol is the polarisation to overwrite.
     * 
     * @return 1 on error, 0 if successful.
     */
    static Int_t directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol);

    /** 
     * @brief Inserts the photogrammetry geometry from AnitaGeomTool into this classes copy of the antenna position. 
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @warning This class overwrites the "extra cable delays" in AnitaEventCalibrator with zero! 
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * Probably don't use this, or if you do, use with caution.
     */    
    void insertPhotogrammetryGeometry();

    /**
     * @brief Coherently summed filtered (un-deconvolved) waveform accessor for external processes.
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @param pol is the polarisation
     * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
     * @param xPol direct or cross polarisation compared to pol.
     *
     * @return pointer to the AnalysisWaveform
     */
    AnalysisWaveform* getCoherentFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);


    /**
     * @brief Coherently summed (un-filtered, un-deconvolved) waveform accessor for external processes.
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @param pol is the polarisation
     * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
     * @param xPol direct or cross polarisation compared to pol.
     *
     * @return pointer to the AnalysisWaveform
     */
    AnalysisWaveform* getCoherent(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);


    /**
     * @brief Coherently summed (un-filtered) deconvolved waveform accessor for external processes.
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @param pol is the polarisation
     * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
     * @param xPol direct or cross polarisation compared to pol.
     *
     * @return pointer to the AnalysisWaveform
     */
    AnalysisWaveform* getDeconvolved(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);


    /**
     * @brief Coherently summed filtered deconvolved waveform accessor for external processes.
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @param pol is the polarisation
     * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
     * @param xPol direct or cross polarisation compared to pol.
     *
     * @return pointer to the AnalysisWaveform
     */
    AnalysisWaveform* getDeconvolvedFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd=0, bool xPol=false);


    /**
     * Access for internally produce minimally filtered event
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @return The minimally filtered version of the event
     */
    FilteredAnitaEvent* getEvMin();


    /**
     * Access for internally produced minimally filtered, deconvolved event
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @return The minimally filtered, deconvolved version of the event
     */
    FilteredAnitaEvent* getEvMinDeco();


    /**
     * Access for internally produced filtered, deconvolved event
     * Only works once per event processed as ownership is transferred to the function caller.
     *
     * @return The deconvolved version of the event
     */
    FilteredAnitaEvent* getEvDeco();


    /**
     * @brief From a list of antennas, get a set of dts relative to the first antenna
     *
     * This function modularizes one of the particular steps in making the coherently summed waveforms
     *
     * @param theAnts are the antenna indices
     * @param pol is the polarisation
     * @param peakPhiDeg is the phi angle relative to ADU5 aft-fore (degrees)
     * @param peakThetaDeg is the theta angle (0 = horizontal, +ve theta is up)
     * @param dts is a vector which gets filled with the time offsets
     */
    void directionAndAntennasToDeltaTs(const std::vector<Int_t>& theAnts, AnitaPol::AnitaPol_t pol,
				       Double_t peakPhiDeg, Double_t peakThetaDeg, std::vector<double>& dts);



    /**
     * @brief Calculate and fill the info the AnitaEventSummary relating the compatibility of the hardware trigger and the interferometric peak
     *
     * @param header is the event header
     * @param pol is the polarization of interest
     * @param peakPhiSector is the phi-sector of the peak in the interferometric map
     * @param peak is a reference to the relevant part of the AnitaEventSummary
     */
    static void setTriggerInfoFromPeakPhi(const RawAnitaHeader* header, AnitaPol::AnitaPol_t pol,
					  Int_t peakPhiSector, AnitaEventSummary::PointingHypothesis& peak);



    /**
     * @brief Calculate and fill the numbers for Peng's ChannelInfo object
     *
     * Whether or not this actually gets done can be controlled with the #fFillChannelInfo setting
     * @see SetFillChannelInfo
     *
     * @param fEv is the filtered event being processed
     * @param sum is the pointer to the AnitaEventSummary
     */
    static void fillChannelInfo(const FilteredAnitaEvent* fEv, AnitaEventSummary* sum);


    /**
     * Fill the power related values in the AnitaEventSummary::EventFlags
     *
     * @param fEv is a pointer to the FilteredAnitaEvent being processed
     * @param flags is a reference to the appropriat bit of the AnitaEventSummary for the processed event
     */
    void fillPowerFlags(const FilteredAnitaEvent* fEv, AnitaEventSummary::EventFlags& flags);


    /**
     * @brief Get a list of antenna indices (in a std::vector) for a coherently averaged waveform based from the phi-sector
     *
     * Depending on where the peak in your InterferometricMap lies you will want to use a subset of antennas to make your coherently summed waveform.
     * This function gives you that subset, from the phi sector, which is an index from 0 to #NUM_PHI - 1.
     *
     * @param peakPhiSector is the phi-sector index.
     *
     * @return a const reference to the vector of antenna indices.
     */
    const std::vector<int>& phiSectorToCoherentAnts(int peakPhiSector){return fPhiSectorToAnts[peakPhiSector];}


    /**
     * @brief Generate a new set of TGraphAligned such that the peaks are aligned,
     *
     * This is only used for plotting during debugging to make sure everything is working.
     *
     * @param waves is a set of waveforms to coherently sum
     * @param dts is the set of time offsets (ns) to apply to those waveforms when summing
     * @param grs are the same as the input AnalysisWaveform graphs but time shifted by dts so their peaks are aligned if plotted together. For debugging.
     *
     * @return a newly created AnalysisWaveform, created by coherently summing the input waves
     */
    void wavesInCoherent(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts, std::vector<TGraphAligned*>& grs);


    /**
     * Get a pointer to the CrossCorrelator, #fCrossCorr
     *
     * @return The #fCrossCorr pointer
     */
    CrossCorrelator* getCrossCorrelator(){return fCrossCorr;}

    /**
     * Get a const reference to the AnitaEventSummary generated from the last event given to process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary, NoiseMonitor* noiseMonitor=NULL, TruthAnitaEvent* truth = NULL)
     *
     * @return a const reference to protected #fSummary
     */
    const AnitaEventSummary& lastSummary(){return fSummary;}

  protected:

    /**
     * @brief Sets default values and zeros pointers for dynamically initialised heap members
     */
    void initializeInternals();

    InterferometricMap* coarseMaps[AnitaPol::kNotAPol]; ///< The coarsely binned InterferometricMap from the most recently processed event
    InterferometricMap* fineMaps[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol]; ///< The finely binned InterferometricMap from the most recently processed event


    /**
     * @brief Bounds checking for coherent averaging. Checks input vectors are the same length and zero pads/trims the vector in the case they're not.
     *
     * @param waves the waveforms to sum
     * @param dts the time offsets to apply
     *
     * @return the size of both vectors
     */
    size_t checkWavesAndDtsMatch(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts);

    /**
     * Does the calculations and stores the data for the AnitaEventSummary::WaveformInfo object
     *
     * @param pol is the polarisation
     * @param info is the selected WaveformInfo in the ANITA event summary
     * @param fEv is the FilteredAnitaEvent from which we are coherently summing
     * @param waveStore is the internal storage in the AnalysisReco class (used to save a copy for MagicDisplay)
     * @param h is the InterferometricMap which contains the peak direction in which we want to coherently sum
     * @param noiseMonitor Contains the min bias RMS values as a function of time, must be non-NULL to fill SNR values
     */
    void fillWaveformInfo(AnitaPol::AnitaPol_t pol, AnitaEventSummary::WaveformInfo& info, const FilteredAnitaEvent* fEv,
			  AnalysisWaveform** waveStore, InterferometricMap* h, NoiseMonitor* noiseMonitor);


    /**
     * @brief Generate a set of antennas to use to generate the coherent waveform.
     *
     * @param coherentDeltaPhi is the number of neighbouring phi-sectors to use.
     */
    void chooseAntennasForCoherentlySumming(int coherentDeltaPhi);

    /**
     * @brief does NULL pointer checking deletion on the internally generated FilteredAnitaEvents, #fEvMin, #fEvMinDeco, #fEvDeco
     */
    void nicelyDeleteInternalFilteredEvents();


    std::vector<Double_t> fRArray[AnitaPol::kNotAPol]; ///< Local copies of the antenna radial positions (metres) from AnitaGeomTool
    std::vector<Double_t> fPhiArrayDeg[AnitaPol::kNotAPol]; ///< Local copies of the antenna azimuth positions (in degrees) from AnitaGeomTool
    std::vector<Double_t> fZArray[AnitaPol::kNotAPol]; ///< Local copies of the antenna z positions (metres) from AnitaGeomTool

    AnalysisWaveform* fCoherentFiltered[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol]; ///< Copies of the coherently summed waveforms from the most recently processed event
    AnalysisWaveform* fCoherent[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol] ; ///< Copies of the unfiltered, coherently summed waveform from the most recently processed event
    AnalysisWaveform* fDeconvolved[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol]; ///< Copies of the dedispersed, unfiltered coherently summed waveform from the most recently processed event
    AnalysisWaveform* fDeconvolvedFiltered[AnitaPol::kNotAPol][AnitaEventSummary::maxDirectionsPerPol][AnitaPol::kNotAPol]; ///< Copies of the dedispersed, filtered coherently summed waveform from the most recently processed event
    std::vector<Int_t> fPhiSectorToAnts[NUM_PHI]; ///< Which antennas to use to make coherently summed waveforms depending on the phi-sector of the peak.

    AnitaEventSummary fSummary; ///< A copy of the AnitaEventSummary from the most recently processed event.
    CrossCorrelator* fCrossCorr; ///< CrossCorrelator generates the set of cross-correlations required to make an InterferometricMap
    bool fSpawnedCrossCorrelator; ///< Did I initialize a CrossCorrelator for #fCrossCorr, or was I given a CrossCorrelator initialized by someone else?
    InterferometryCache fDtCache; ///< Caches antenna delays as a function of incoming angle, for quickly making an Interferometric map

    FilterStrategy* fMinFilter; ///< Minimum set of filters (no filers, or just ALFA filters for ANITA-3)
    FilterStrategy* fMinDecoFilter; ///< Minimum filter with an appended deconvolution filter

    FilteredAnitaEvent* fEvMin; ///< Filtered event produced with the #fMinFilter
    FilteredAnitaEvent* fEvMinDeco; ///< Filtered event produced with the #fMinDecoFilter
    FilteredAnitaEvent* fEvDeco; ///< FilteredEvent produced by appending a deconvolution filter to the filter strategy passed to process()

    UInt_t fCurrentEventNumber; ///< Assigned at the start of process(), helpful for printing warning/info messages
    Int_t fCurrentRun; ///< Assigned at the start of process(), helpful for printing warning/info messages

    /**
     * @fn GetDebug
     * @brief Get the value of #fDebug
     */
    /**
     * @fn SetDebug
     * @brief Set the value of #fDebug
     */
    /**
     * @var fDebug
     * @brief Controls the printing of debug info, feel free to wrap noisy stderr messages with if(fDebug),
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, Debug);

    /**
     * @fn GetCoherentDeltaPhi
     * @brief Get the value of #fCoherentDeltaPhi
     */
    /**
     * @fn SetCoherentDeltaPhi
     * @brief Set the value of #fCoherentDeltaPhi
     */
    /**
     * @var fCoherentDeltaPhi
     * @brief The +/- range of phi-sectors around the map peak phi-sector to use when making a coherently averaged waveform, i.e. antennas from the nearest 1+2*#fCoherentDeltaPhi phi-sectors are used.
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, CoherentDeltaPhi);
    Int_t fLastCoherentDeltaPhi; ///< Used to figure out if we need to recalculate the cached values related to coherently summing waveforms


    /**
     * @fn GetWhichResponseDir
     * @brief Get the value of #fWhichResponseDir
     */
    /**
     * @fn SetWhichResponseDir
     * @brief Set the value of #fWhichResponseDir
     */
    /**
     * @var fWhichResponseDir
     * @brief Which version of the antenna response to use?
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, WhichResponseDir);

    /**
     * @fn GetUseOffAxisDelay
     * @brief Get the value of #fUseOffAxisDelay
     */
    /**
     * @fn SetUseOffAxisDelay
     * @brief Set the value of #fUseOffAxisDelay
     */
    /**
     * @var fUseOffAxisDelay
     * @brief Whether or not to use the derived off-axis delay values when making an InterferometricMap
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, UseOffAxisDelay);


    /**
     * @fn GetResponseNPad
     * @brief Get the value of #fResponseNPad
     */
    /**
     * @fn SetResponseNPad
     * @brief Set the value of #fResponseNPad
     */
    /**
     * @var fResponseNPad
     * @brief How up to how many samples do we need to pad the antenna response functions?
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, ResponseNPad);

    /**
     * @fn GetNumPeaks
     * @brief Get the value of #fNumPeaks
     */
    /**
     * @fn SetNumPeaks
     * @brief Set the value of #fNumPeaks
     */
    /**
     * @var fNumPeaks
     * @brief How many local maxima in the coarsely binned InterferometricMap should we make a finely binned InterferometricMap? (per polarization).
     * 
     * Must be less than or equal to AnitaEventSummary::maxDirectionsPerPol. Default is 3.
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, NumPeaks);

    /**
     * @fn GetCoherentDtNs
     * @brief Get the value of #fCoherentDtNs
     */
    /**
     * @fn SetCoherentDtNs
     * @brief Set the value of #fCoherentDtNs
     */
    /**
     * @var fCoherentDtNs
     * @brief Time in nano-seconds to interpolate to when coherently averaging waveforms
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Double_t, CoherentDtNs);

    /**
     * @fn GetSlopeFitStartFreqGHz
     * @brief Get the value of #fSlopeFitStartFreqGHz
     */
    /**
     * @fn SetSlopeFitStartFreqGHz
     * @brief Set the value of #fSlopeFitStartFreqGHz
     */
    /**
     * @var fSlopeFitStartFreqGHz
     * @brief Frequency at which to start the linear fit to the power spectrum of the coherently averaged waveforms
     * @see AnalysisSettings
     */    
    ANALYSIS_SETTING(Double_t, SlopeFitStartFreqGHz);

    /**
     * @fn GetSlopeFitEndFreqGHz
     * @brief Get the value of #fSlopeFitEndFreqGHz
     */
    /**
     * @fn SetSlopeFitEndFreqGHz
     * @brief Set the value of #fSlopeFitEndFreqGHz
     */
    /**
     * @var fSlopeFitEndFreqGHz
     * @brief Frequency at which to stop the linear fit to the power spectrum of the coherently averaged waveforms
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Double_t, SlopeFitEndFreqGHz);

    
    /**
     * @fn GetMeanPowerFlagLowFreqGHz
     * @brief Get the value of #fMeanPowerFlagLowFreqGHz
     */
    /**
     * @fn SetMeanPowerFlagLowFreqGHz
     * @brief Set the value of #fMeanPowerFlagLowFreqGHz
     */
    /**
     * @var fMeanPowerFlagLowFreqGHz
     * @brief Frequency at which to start the mean power calculation
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Double_t, MeanPowerFlagLowFreqGHz);

    /**
     * @fn GetMeanPowerFlagHighFreqGHz
     * @brief Get the value of #fMeanPowerFlagHighFreqGHz
     */
    /**
     * @fn SetMeanPowerFlagHighFreqGHz
     * @brief Set the value of #fMeanPowerFlagHighFreqGHz
     */
    /**
     * @var fMeanPowerFlagHighFreqGHz
     * @brief Frequency at which to stop the mean power calculation
     * @see AnalysisSettings
     */    
    ANALYSIS_SETTING(Double_t, MeanPowerFlagHighFreqGHz);

    /**
     * @fn GetFillChannelInfo
     * @brief Get the value of #fFillChannelInfo
     */
    /**
     * @fn SetFillChannelInfo
     * @brief Set the value of #fFillChannelInfo
     */
    /**
     * @var fFillChannelInfo
     * @brief Whether or not to fill Peng's AnitaEventSummary::ChannelInfo (smaller output files from better compression if you don't)
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, FillChannelInfo);

    /**
     * @fn GetFillSpectrumInfo
     * @brief Get the value of #fFillSpectrumInfo
     */
    /**
     * @fn SetFillSpectrumInfo
     * @brief Set the value of #fFillSpectrumInfo
     */
    /**
     * @var fFillSpectrumInfo
     * @brief Whether or not to fill the spectrum info (smaller output files from better compression if you don't)
     * @see AnalysisSettings
     */    
    ANALYSIS_SETTING(Int_t, FillSpectrumInfo);

    /**
     * @fn GetFillUnfiltered
     * @brief Get the value of #fFillUnfiltered
     */
    /**
     * @fn SetFillUnfiltered
     * @brief Set the value of #fFillUnfiltered
     */
    /**
     * @var fFillUnfiltered
     * @brief  Whether or not to fill the AnitaEventSummary::WaveformInfo objects for the unfiltered versions of the AnalysisWaveforms (smaller output files from better compression if you don't)
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, FillUnfiltered);












    /**
     * He we get variables and functions to tweak the GUI options for #DrawSummary
     */




    /**
     * @fn GetDrawNPeaks
     * @brief Get the value of #fDrawNPeaks
     */

    /**
     * @fn SetDrawNPeaks
     * @brief Set the value of fDrawNPeaks
     */
    /**
     * @var fDrawNPeaks
     * @brief How many peaks to draw in the summary, must be less than or equal to fNumPeaks
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, DrawNPeaks)


    /**
     * @fn GetEnumDrawDomain
     * @brief Get the value of #fDrawDomain as AnalysisFlow::selection enum
     * @see AnalysisSettings
     */
    /**
     * @fn GetDrawDomain
     * @brief Get the value of #fDrawDomain as an integer
     * @see AnalysisSettings
     */
    /**
     * @fn SetDrawDomain(Int_t val)
     * @brief Set the value of #fDrawDomain from an integer value
     * @see AnalysisSettings
     */
    /**
     * @fn SetDrawDomain(AnitaDataset::BlindingStrategy val)
     * @brief Set the value of #fDrawDomain from an enum value
     * @see AnalysisSettings
     */
    /**
     * @var fDrawDomain
     * @brief Set as 0 to draw time domain, or 1 to draw frequency domain
     * @see AnalysisSettings
     */
    ENUM_ANALYSIS_SETTING(DrawDomain, DrawDomain)


    /**
     * @fn GetDrawCoherent
     * @brief Get the value of #fDrawCoherent
     */

    /**
     * @fn SetDrawCoherent
     * @brief Set the value of #fDrawCoherent
     */
    /**
     * @var fDrawCoherent
     * @brief
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, DrawCoherent)


    /**
     * @fn GetDrawDedispersed
     * @brief Get the value of #fDrawDedispersed
     */

    /**
     * @fn SetDrawDedispersed
     * @brief Set the value of #fDrawDedispersed
     */
    /**
     * @var fDrawDedispersed
     * @brief
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, DrawDedispersed)


    /**
     * @fn GetDrawXPol
     * @brief Get the value of #fDrawXPol
     */

    /**
     * @fn SetDrawXPol
     * @brief Set the value of #fDrawXPol
     */
    /**
     * @var fDrawXPol
     * @brief
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, DrawXPol)


    /**
     * @fn GetDrawXPolDedispersed
     * @brief Get the value of #fDrawXPolDedispersed
     */

    /**
     * @fn SetDrawXPolDedispersed
     * @brief Set the value of #fDrawXPolDedispersed
     */
    /**
     * @var fDrawXPolDedispersed
     * @brief
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, DrawXPolDedispersed)





    ClassDef(AnalysisReco, 0)

  };

}


#endif
