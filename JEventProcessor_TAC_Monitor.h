/*
 * JEventProcessor_TAC_Monitor.h
 *
 *  Created on: Mar 30, 2017
 *      Author: hovanes
 */

#ifndef JEVENTPROCESSOR_TAC_MONITOR_H_
#define JEVENTPROCESSOR_TAC_MONITOR_H_

#include <iostream>
#include <map>
#include <vector>
#include <iterator>
#include <algorithm>

#include <TH1.h>

#include <JANA/JEventProcessor.h>
#include <DANA/DApplication.h>
#include <TRIGGER/DL1Trigger.h>
#include <DAQ/Df250WindowRawData.h>

#include "CompressionTester.h"

class JEventProcessor_TAC_Monitor: public jana::JEventProcessor {
protected:

	// Map of all histograms for this monitoring plugin. the first index identifies the
	// name of the histogram, the second index (inner) identifies the trigger bit.
	std::map<std::string, std::map<unsigned,TH1*> > histoMap;

	// ROOT file name
	std::string rootFileName = "tac_monitor.root";

	// ROOT directory pointer
	TDirectory* rootDir = nullptr;

	CompressionTester* dataCompressor = nullptr;

	// Mask indicating which trigger bits this class cares for.
	static uint32_t triggerMask;
//	// Mask that specifies the TAC trigger bit
//	static uint32_t tacTriggerBit;
	// Maximum numbr of trigger bits considered in this plugin
	static uint32_t numberOfTriggerBits;
	// FADC value starting from which the FADC is considered overflowed
	static unsigned maxPulseValue;
	// Minimum peak value
	static unsigned minPulseValue;
	// Value that assigned to the peak in cases of overlfow
	static unsigned overflowPulseValue;

	// Timing cut value between the TAGH and TAC coincidence
	static double timeCutValue_TAGH;
	// Timing cut width between the TAGH and TAC coincidence
	static double timeCutWidth_TAGH;

	// Timing cut value between the TAGH and TAC coincidence
	static double timeCutValue_TAGM;
	// Timing cut width between the TAGH and TAC coincidence
	static double timeCutWidth_TAGM;

	// Threshold that will define when the TAC hit occurred
	static unsigned tacThreshold;

	// Time units for the timing from the Digi bank
	static double fadc250DigiTimeScale;
	// Time units for the timing from raw FADC
	static double fadc250RawTimeScale;

	virtual jerror_t init(void);          ///< Called once at program start.
	virtual jerror_t brun(jana::JEventLoop *eventLoop, int32_t runNumber);          ///< Called everytime a new run number is detected.
	virtual jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventNumber);          ///< Called every event.
	virtual jerror_t erun(void);          ///< Called every time run number changes, provided brun has been called.
	virtual jerror_t fini(void);          ///< Called after last event of last event source has been processed.

	// Method where the histograms are created
	virtual void createHistograms();
	// Fill raw data histograms (the ones related to waveforms
	virtual jerror_t fillRawDataHistograms(jana::JEventLoop* eventLoop,
			uint32_t trigBit);
	// Fill pulse data histograms
	virtual jerror_t fillPulseDataHitograms(jana::JEventLoop* eventLoop,
			uint32_t trigBit);

	// Fill F1TDC related histograms
	virtual jerror_t fillTDCHistograms( jana::JEventLoop* eventLoop, uint32_t trigBit );

	// Fill tagger related histos
	template<typename DATA_TYPE, typename CounterID, typename TIME_CUT>
	jerror_t fillTaggerRelatedHistograms(
			std::vector<const DATA_TYPE*>& digiHitVector, uint32_t trigBit,
			std::string detComp, std::string tacMethod, double tacPeak,
			double tacTime, CounterID idFunctor, TIME_CUT timeCut);

	// Return a pair giving the peak location (first) and the peak value (second)
	virtual std::pair<unsigned,unsigned> getPeakLocationAndValue( const Df250WindowRawData* data );
	// Return the time where the signal crosses the threshold
	virtual unsigned getPulseTime( const Df250WindowRawData* data, unsigned threshold );

	template<typename TH1_TYPE>
	jerror_t createHisto(unsigned trigBit, std::string key,
			std::string titlePrefix, std::string xTitle, int nBins, double xMin,
			double xMax);
	template<typename TH2_TYPE>
	jerror_t createHisto(unsigned trigBit,
			std::string histKey, std::string xTitlePrefix, std::string xTitle,
			std::string yTitle, int nBinsX, double xMin, double xMax,
			int nBinsY, double yMin, double yMax);

	// Write histograms into the file
	virtual jerror_t writeHistograms();

	// Check the file compression by writing out some files.
	static jerror_t writeRawData( const Df250WindowRawData* tacRawData );

	// Check if the trigger bits for the event are useful
	static bool triggerIsUseful(const DL1Trigger* trigWords) {
		if (trigWords) {
			return triggerIsUseful( trigWords->trig_mask );
		}
		return false;
	}
	// Check if the trigger bits for the event are useful
	static bool triggerIsUseful(unsigned trigBits) {
		if ((trigBits & triggerMask) == 0) {
			return false;
		}
		return true;
	}

public:

	JEventProcessor_TAC_Monitor(){}
	virtual ~JEventProcessor_TAC_Monitor(){}

//	const char* className(void) {
//		return "JEventProcessor_TAC_Monitor";
//	}

	static uint32_t getTriggerMask() {
		return triggerMask;
	}

	static void setTriggerMask(uint32_t mask) {
		triggerMask = mask;
	}

	const TDirectory* getRootDir() {
		return rootDir;
	}

	void setRootDir(TDirectory* rootDir) {
		this->rootDir = rootDir;
	}
};

#endif /* JEVENTPROCESSOR_TAC_MONITOR_H_ */
