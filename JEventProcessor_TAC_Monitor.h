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

#include <TH1.h>

#include "TRIGGER/DL1Trigger.h"
#include <JANA/JEventProcessor.h>

class JEventProcessor_TAC_Monitor: public jana::JEventProcessor {
protected:

	// Map of all histograms for this monitoring plugin. the first index identifies the
	// name of the histogram, the second index (inner) identifies the trigger bit.
	std::map<std::string, std::map<unsigned,TH1*> > histoMap;

	// ROOT file name
	std::string rootFileName = "tac_monitor.root";

	// Mask indicating which trigger bits this class cares for.
	static uint32_t triggerMask;
	static unsigned maxPulseValue;
	static unsigned minPulseValue;
	static unsigned overflowPulseValue;


	virtual jerror_t init(void);          ///< Called once at program start.
	virtual jerror_t brun(jana::JEventLoop *eventLoop, int32_t runNumber);          ///< Called everytime a new run number is detected.
	virtual jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventNumber);          ///< Called every event.
	virtual jerror_t erun(void);          ///< Called everytime run number changes, provided brun has been called.
	virtual jerror_t fini(void);          ///< Called after last event of last event source has been processed.


	// Method where the histograms are created
	virtual void createHistograms();
	// Fill raw data histograms (the ones related to waveforms
	virtual jerror_t fillRawDataHistograms( jana::JEventLoop* eventLoop, uint32_t trigBits );
	// Fill pulse data histograms
	virtual jerror_t fillPulseDataHitograms( jana::JEventLoop* eventLoop, uint32_t trigBits );

	// Write histograms into the file
	virtual jerror_t writeHistograms();


	// Check if the trigger bits for the event are useful
	static bool triggerIsUseful(const DL1Trigger* trigWords) {
		if (trigWords) {
			if ((trigWords->trig_mask & triggerMask) == 0)
				return false;
		} else {
			return false;
		}
		return true;
	}
public:
	JEventProcessor_TAC_Monitor(){}
	virtual ~JEventProcessor_TAC_Monitor(){}

	const char* className(void) {
		return "JEventProcessor_TAC_Monitor";
	}

};

#endif /* JEVENTPROCESSOR_TAC_MONITOR_H_ */
