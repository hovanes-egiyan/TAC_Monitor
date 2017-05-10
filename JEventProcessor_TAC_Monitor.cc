/*
 * JEventProcessor_TAC_Monitor.cc
 *
 *  Created on: Mar 30, 2017
 *      Author: hovanes
 */

#include <iostream>
#include <map>
#include <vector>
#include <sstream>

#include "TApplication.h"  // needed to display canvas
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"

#include <JANA/JApplication.h>

#include "TRIGGER/DL1Trigger.h"
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseData.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/DF1TDCHit.h>

#include <DAQ/DEPICSvalue.h>
#include <TAC/DTACDigiHit.h>
#include <TAC/DTACTDCDigiHit.h>
#include <TAC/DTACHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHTDCDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>
#include <TTAB/DTTabUtilities.h>

//#include "libRootSpy/DRootSpy.h"

#include "ReadWriteLock.h"
#include "JEventProcessor_TAC_Monitor.h"

using namespace jana;
using namespace std;

// Routine used to create our JEventProcessor
extern "C" {
void InitPlugin(JApplication *app) {
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_TAC_Monitor());
}
}

// Mask that specifies the bits of interest for TAC runs
uint32_t JEventProcessor_TAC_Monitor::triggerMask = 0b00000010;
// Maximum numbr of trigger bits considered in this plugin
uint32_t JEventProcessor_TAC_Monitor::numberOfTriggerBits = 16;
//// Mask that specifies the TAC trigger bit
//uint32_t JEventProcessor_TAC_Monitor::tacTriggerBit = 1 << 1;
// Mimimum ADC values that is used
unsigned JEventProcessor_TAC_Monitor::minPulseValue = 400;
// ADC value which indicates overflow (maximum+1)
unsigned JEventProcessor_TAC_Monitor::maxPulseValue = 4095;
// ADC value assigned for overlfows
unsigned JEventProcessor_TAC_Monitor::overflowPulseValue = 4500;

// Timing cut value between the TAGH and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutValue_TAGH = 100.0;
// Timing cut width between the TAGH and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutWidth_TAGH = 20;

// Timing cut value between the TAGM and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutValue_TAGM = 90.0;
// Timing cut width between the TAGM and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutWidth_TAGM = 20;

// Threshold that will define when the TAC hit occurred
unsigned JEventProcessor_TAC_Monitor::tacThreshold = 200;

// Time units for the timing from the Digi bank
double JEventProcessor_TAC_Monitor::fadc250DigiTimeScale = (1./16.0);
// Time units for the timing from raw FADC
double JEventProcessor_TAC_Monitor::fadc250RawTimeScale = 4.0;


jerror_t JEventProcessor_TAC_Monitor::init(void) {
	volatile WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());

	// Create parameters and assign values
	gPARMS->SetDefaultParameter<string,double>( "TAC:TAGH_FADC_MEAN_TIME", timeCutValue_TAGH );
	gPARMS->GetParameter( "TAC:TAGH_FADC_MEAN_TIME" )->GetValue( timeCutValue_TAGH );
	gPARMS->SetDefaultParameter<string,double>( "TAC:TAGM_FADC_MEAN_TIME", timeCutValue_TAGM );
	gPARMS->GetParameter( "TAC:TAGM_FADC_MEAN_TIME" )->GetValue( timeCutValue_TAGM );
	gPARMS->SetDefaultParameter<string,unsigned>( "TAC:TAC_FADC_THRESHOLD", tacThreshold );
	gPARMS->GetParameter( "TAC:TAC_FADC_THRESHOLD" )->GetValue( tacThreshold );

	// Create TAC directory and the histograms
	TDirectory *mainDir = gDirectory;
	rootDir = gDirectory->mkdir("TAC");
	rootDir->cd();
	createHistograms();
	mainDir->cd();
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::brun(jana::JEventLoop* eventLoop,
		int32_t runnumber) {
	stringstream fileNameStream;
	fileNameStream << "tac_monitor_" << runnumber << ".root" ;
	rootFileName = fileNameStream.str();
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::evnt(jana::JEventLoop* eventLoop,
		uint64_t eventNumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// do not do anything if the application is not DApplication
	if (dynamic_cast<DApplication*>(japp) == nullptr)
		return NOERROR;

	// Get First Trigger Type
	const DL1Trigger *trigWords = nullptr;
	try {
		eventLoop->GetSingle(trigWords);
	} catch (...) {
	};

	// Decide if to continue considering this event based on the trigger bit pattern
	if (!triggerIsUseful(trigWords))
		return NOERROR;
	uint32_t usefulTriggerBits = triggerMask & trigWords->trig_mask;

	// Here we fill the raw waveforms
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned singleBit = 1 << trigBit;
		if ((singleBit & usefulTriggerBits) != 0) {
			this->fillRawDataHistograms(eventLoop, trigBit);
			this->fillPulseDataHitograms(eventLoop, trigBit);
			this->fillTDCHistograms(eventLoop, trigBit);
		}
	}
	// Write histograms into ROOT file once in a while
	if( eventNumber % 200000 == 0 ) {
		this->writeHistograms();
	}

	return NOERROR;
}

// Handle histograms with FADC250 raw data
jerror_t JEventProcessor_TAC_Monitor::fillRawDataHistograms(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {
	vector<const Df250WindowRawData*> rawDataVector;
	eventLoop->Get(rawDataVector);

	unsigned tacDataCounter = 0;
	// Pick the TAC raw data
	const Df250WindowRawData* tacRawData = nullptr;
	for (auto& rawData : rawDataVector) {
		if (rawData) {
			// Hard-coded the location of the TAC FADC signal, pick the TAC raw data
			if (rawData->rocid == 14 && rawData->slot == 20
					&& rawData->channel == 0) {
				tacRawData = rawData;
				tacDataCounter++;
			}
		}
	}

	unsigned tacTDCDataCounter = 0;
	vector<const DCAEN1290TDCHit*> rawTDCDataVector;
	eventLoop->Get(rawTDCDataVector);

	const DDAQAddress* tacTDCRawData = nullptr;
	for (auto& rawTDCData : rawTDCDataVector) {
		if (rawTDCData) {
//			cout << rawTDCData->rocid << " , " << rawTDCData->slot << "  , " << rawTDCData->channel << endl;
			// Hard-coded the location of the TAC FADC signal, pick the TAC raw data
			if (rawTDCData->rocid == 78 && rawTDCData->slot == 8
					&& rawTDCData->channel == 18) {
				tacTDCRawData = rawTDCData;
				tacTDCDataCounter++;
//				cout << "Found TAC TDC" << endl;
			}
		}
	}
//	if( tacTDCDataCounter > 0 )
//		cout << "Found TAC TDC data " << endl;


	if (tacDataCounter < 1) {
		cout << "Too few TAC raw hits: " << tacDataCounter << endl;
	}

	if (tacDataCounter > 1) {
		cout << "Too many TAC raw hits: " << tacDataCounter << endl;
	}
	// Only analyze events with a single hit in the FADC
	if (tacDataCounter != 1)
		return NOERROR;

	// Fill the waveform histograms
	{
		volatile WriteLock rootRWLock(
				*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
		int binNumber = 0;
		for (auto& rawDataValue : tacRawData->samples) {
			binNumber++;
			if (binNumber <= histoMap["TACFADCRAW"][trigBit]->GetNbinsX()) {
				histoMap["TACFADCRAW"][trigBit]->SetBinContent(binNumber,
						rawDataValue);
				histoMap["TACFADCRAW_ENTRIES"][trigBit]->SetBinContent(
						binNumber,
						histoMap["TACFADCRAW_ENTRIES"][trigBit]->GetBinContent(
								binNumber) + 1.0);
				histoMap["TACFADCRAW_SUM"][trigBit]->SetBinContent(binNumber,
						histoMap["TACFADCRAW_SUM"][trigBit]->GetBinContent(
								binNumber) + rawDataValue);
			}
		}
		histoMap["TACFADCRAW_AVG"][trigBit]->Divide(
				histoMap["TACFADCRAW_SUM"][trigBit],
				histoMap["TACFADCRAW_ENTRIES"][trigBit], 1, 1);
	}

	// Find the maximum by going through the raw data and comparing samples
	pair<unsigned, unsigned> maxInfoPair = this->getPeakLocationAndValue(
			tacRawData);
	double maxValue = maxInfoPair.second;
	// Find the TAC signal time by finding the bin where the signal is above threshold
	double tacPeakTime = this->getPulseTime( tacRawData, tacThreshold );
	// Assign a bigger values for cases with overflows
	if (maxValue >= maxPulseValue)
		maxValue = overflowPulseValue;
	{
		volatile WriteLock rootRWLock(
				*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
		histoMap["TACAmpWAVE"][trigBit]->Fill(maxValue);
		histoMap["TACTimeWAVE"][trigBit]->Fill(tacPeakTime*fadc250RawTimeScale);
	}

	// Call methods to fill tagger (TAGH and TAGM) related histograms
	vector<const DTAGHDigiHit*> taghDigiHitVector;
	eventLoop->Get(taghDigiHitVector);
	fillTaggerRelatedHistograms(taghDigiHitVector, trigBit, "TAGH", "WAVE", maxValue,
			tacPeakTime, [&]( const DTAGHDigiHit* hit ) {return hit->counter_id;},
			[&](const DTAGHDigiHit* hit ) -> bool {return fabs( hit->pulse_time*fadc250DigiTimeScale - timeCutValue_TAGH ) < timeCutWidth_TAGH;});
	vector<const DTAGMDigiHit*> tagmDigiHitVector;
	eventLoop->Get(tagmDigiHitVector);
	fillTaggerRelatedHistograms(tagmDigiHitVector, trigBit, "TAGM", "WAVE", maxValue,
			tacPeakTime, [&]( const DTAGMDigiHit* hit ) {return hit->column;},
			[&](const DTAGMDigiHit* hit ) -> bool {return fabs( hit->pulse_time*fadc250DigiTimeScale - timeCutValue_TAGM ) < timeCutWidth_TAGM;});

	return NOERROR;
}

// Handle histogram from FADC250 pulse data
jerror_t JEventProcessor_TAC_Monitor::fillPulseDataHitograms(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {
	vector<const DTACDigiHit*> tacDigiHitVector;
	eventLoop->Get(tacDigiHitVector);

	// Pick the only peak value from the TACDigiHit object
	double pulsePeak = 0;
	double pulseTime = 0;
	{
		volatile WriteLock rootRWLock(
				*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
		histoMap["TAC_NHITS"][trigBit]->Fill(tacDigiHitVector.size());
	}
	// Find the digi hit with the largest pulse and use its height and time
	for (auto tacDigiHit : tacDigiHitVector) {
		if (tacDigiHit) {
			double currentPeak = double(tacDigiHit->getPulsePeak());
			if (currentPeak > pulsePeak) {
				pulsePeak = currentPeak;
				pulseTime = double(tacDigiHit->getPulseTime())
						* fadc250DigiTimeScale;
			}
		}
	}
	// Assign larger value when overflow is detected in FADC
	if (pulsePeak >= maxPulseValue)
		pulsePeak = overflowPulseValue;
	{
		volatile WriteLock rootRWLock(
				*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
		histoMap["TACAmpPULSE"][trigBit]->Fill(pulsePeak);
		histoMap["TACTimePULSE"][trigBit]->Fill(pulseTime);
	}
	vector<const DTAGHDigiHit*> taghDigiHitVector;
	eventLoop->Get(taghDigiHitVector);
	fillTaggerRelatedHistograms(taghDigiHitVector, trigBit, "TAGH", "PULSE",
			pulsePeak, pulseTime,
			[&]( const DTAGHDigiHit* hit ) {return hit->counter_id;},
			[&](const DTAGHDigiHit* hit ) -> bool {return fabs( hit->pulse_time*fadc250DigiTimeScale - timeCutValue_TAGH ) < timeCutWidth_TAGH;});

	vector<const DTAGMDigiHit*> tagmDigiHitVector;
	eventLoop->Get(tagmDigiHitVector);
	fillTaggerRelatedHistograms(tagmDigiHitVector, trigBit, "TAGM", "PULSE",
			pulsePeak, pulseTime,
			[&]( const DTAGMDigiHit* hit ) {return hit->column;},
			[&](const DTAGMDigiHit* hit ) -> bool {return fabs( hit->pulse_time*fadc250DigiTimeScale - timeCutValue_TAGM ) < timeCutWidth_TAGM;});
	return NOERROR;
}


jerror_t JEventProcessor_TAC_Monitor::fillTDCHistograms(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {
	vector<const DTACTDCDigiHit*> tacTDCDigiHits;
	eventLoop->Get(tacTDCDigiHits);

	const DTTabUtilities* ttabUtilities = nullptr;
	eventLoop->GetSingle(ttabUtilities);

	for (auto& tacTDCDigiHit : tacTDCDigiHits) {
		if (tacTDCDigiHit) {
			const DF1TDCHit* tacF1RawHit = nullptr;
			tacTDCDigiHit->GetSingleT(tacF1RawHit);
			if (tacTDCDigiHit) {
				double tacTDCTime = ttabUtilities->Convert_DigiTimeToNs_F1TDC(
						tacF1RawHit);
				histoMap["TAC_TDCTIME"][trigBit]->Fill(double(tacTDCTime));
			}
		}
	}
	return NOERROR;
}


jerror_t JEventProcessor_TAC_Monitor::erun(void) {
	this->writeHistograms();
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::fini(void) {
//	TDirectory* oldDir = gDirectory;
//	rootDir->cd();
//	for (auto& histMapIter : histoMap) {
//		auto& histMap = histMapIter.second;
//		for (auto& histIter : histMap) {
//			auto histPtr = histIter.second;
////			if (histPtr != 0)
////				histPtr->Write();
//		}
//	}
//	oldDir->cd();
	return NOERROR;
}

void JEventProcessor_TAC_Monitor::createHistograms() {
	cout << "Creating TAC histos" << endl;
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned trigPattern = 1 << trigBit;
		if (triggerIsUseful(trigPattern)) {
			// Create TAC FADc raw data
			createHisto<TH1D>(trigBit, "TACFADCRAW", "Single TAC FADC waveform ",
					"FlashADC sample number [#]", 100, 0., 100.);
			// Create TAC summed FADc raw data
			createHisto<TH1D>(trigBit, "TACFADCRAW_SUM",
					"Summed TAC FADC waveform", "FlashADC sample number [#]",
					100, 0., 100.);
			// Create TAC FADc raw data for entries
			createHisto<TH1D>(trigBit, "TACFADCRAW_ENTRIES",
					"Entries in TAC FADC waveform bit ",
					"FlashADC sample number [#]", 100, 0., 100.);
			// Create TAC averaged FADc raw data
			createHisto<TH1D>(trigBit, "TACFADCRAW_AVG",
					"Averaged TAC FADC waveform", "FlashADC sample number [#]",
					100, 0., 100.);

			// Create TAC number of hits histogram
			createHisto<TH1D>(trigBit, "TAC_NHITS",
					"Number of hits in TAC", "number of hits from FADC FPGA [#]",
					7, 0., 7.);

			// Create TAC TDC hit time
			createHisto<TH1D>(trigBit, "TAC_TDCTIME",
					"TDC time in TAC", "TDC time [ns]",
					500, 0., 500.);

			// Create TAC amplitude histos
			createHisto<TH1D>(trigBit, "TACAmpPULSE",
					"TAC Signal Amplitude for Trigger ", "TAC Amplitude", 500,
					0., 5000.);
			// Create TAC amplitude histos for going through the data and picking the highest bin
			createHisto<TH1D>(trigBit, "TACAmpWAVE",
					"TAC Signal Maximum from Raw for Trigger ", "TAC Amplitude",
					500, 0., 5000.);
			// Create TAC signal time histo
			createHisto<TH1D>(trigBit, "TACTimePULSE",
					"TAC Signal time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);
			// Create TAC signal time based on raw data histo
			createHisto<TH1D>(trigBit, "TACTimeWAVE",
					"TAC Signal based on raw data time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);

			// Create TAGH Hits detector ID
			createHisto<TH1D>(trigBit, "TAGH_ID",
					"TAGH Hits Detector ID for Trigger ", "Tagger Hodoscope Det. Number [#]",
					320, 0., 320.);
			// Create TAGH Hits detector ID
			createHisto<TH1D>(trigBit, "TAGH_ID_MATCHEDPULSE",
					"Matched TAGH Hits Detector ID for Trigger ", "Tagger Hodoscope Det. Number [#]",
					320, 0., 320.);
			createHisto<TH1D>(trigBit, "TAGH_ID_MATCHEDWAVE",
					"Matched TAGH Hits Detector ID for Trigger ", "Tagger Hodoscope Det. Number [#]",
					320, 0., 320.);
			// Create TAGH signal time histo
			createHisto<TH1D>(trigBit, "TAGHSigTime",
					"TAGH Signal time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);
			// Create TAC time vs TAGH FADC time histo
			createHisto<TH2D>(trigBit, "TACTIMEPULSEvsTAGHTIME",
					"TAC time vs TAGH time ", "FlashADC peak time for TAGH (ns)", "FlashADC peak time for TAC (ns)",
					400, 0., 400., 400, 0., 400. );
			createHisto<TH2D>(trigBit, "TACTIMEWAVEvsTAGHTIME",
					"TAC time vs TAGH time ", "FlashADC peak time for TAGH (ns)", "FlashADC peak time for TAC (ns)",
					400, 0., 400., 400, 0., 400. );
			// Create TAC amplitude vs TAGH ID histo
			createHisto<TH2D>(trigBit, "TACAMPPULSEvsTAGHID",
					"TAC FADC Amplitude vs TAGH ID ", "Tagger Hodoscope Det. Number [#]", "FlashADC peak for TAC",
					320, 0., 320., 1000, 10., 5000. );
			createHisto<TH2D>(trigBit, "TACAMPWAVEvsTAGHID",
					"TAC FADC Amplitude vs TAGH ID ", "Tagger Hodoscope Det. Number [#]", "FlashADC peak for TAC",
					320, 0., 320., 1000, 10., 5000. );
			// Create TAGH time vs TAGH ID histo
			createHisto<TH2D>(trigBit, "TAGHTIMEvsTAGHID",
					"TAGH Time vs TAGH ID ", "Tagger Hodoscope Det. Number [#]", "TAGH time",
					320, 0., 320., 400, 0., 400. );

			// Create TAGM Hits detector ID
			createHisto<TH1D>(trigBit, "TAGM_ID",
					"TAGM Hits Detector ID for Trigger ", "Tagger Microscope Det. Number [#]",
					110, 0., 110.);
			// Create TAGM Hits detector ID
			createHisto<TH1D>(trigBit, "TAGM_ID_MATCHEDPULSE",
					"Matched TAGM Hits Detector ID for Trigger ", "Tagger Microscope Det. Number [#]",
					110, 0., 110.);
			createHisto<TH1D>(trigBit, "TAGM_ID_MATCHEDWAVE",
					"Matched TAGM Hits Detector ID for Trigger ", "Tagger Microscope Det. Number [#]",
					110, 0., 110.);
			// Create TAGM signal time histo
			createHisto<TH1D>(trigBit, "TAGMSigTime",
					"TAGM Signal time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);
			createHisto<TH2D>(trigBit, "TACTIMEPULSEvsTAGMTIME",
					"TAC time vs TAGM time ", "FlashADC peak time for TAGM (ns)", "FlashADC peak time for TAC (ns)",
					400, 0., 400., 400, 0., 400. );
			createHisto<TH2D>(trigBit, "TACTIMEWAVEvsTAGMTIME",
					"TAC time vs TAGM time ", "FlashADC peak time for TAGM (ns)", "FlashADC peak time for TAC (ns)",
					400, 0., 400., 400, 0., 400. );
			// Create TAC amplitude vs TAGM ID histo
			createHisto<TH2D>(trigBit, "TACAMPPULSEvsTAGMID",
					"TAC FADC Amplitude vs TAGM ID ", "Tagger Microscope Det. Number [#]", "FlashADC peak for TAC",
					110, 0., 110., 1000, 10., 5000. );
			createHisto<TH2D>(trigBit, "TACAMPWAVEvsTAGMID",
					"TAC FADC Amplitude vs TAGM ID ", "Tagger Microscope Det. Number [#]", "FlashADC peak for TAC",
					110, 0., 110., 1000, 10., 5000. );
			// Create TAGH time vs TAGH ID histo
			createHisto<TH2D>(trigBit, "TAGMTIMEvsTAGMID",
					"TAGM Time vs TAGM ID ", "Tagger Microscope Det. Number [#]", "TAGM time",
					110, 0., 110., 400, 0., 400. );
		}
	}
}

// Create a 1D histogram of type TH1_TYPE and assign it to the histogram map based on the argument valeus
// provided in the function call.
template<typename TH1_TYPE>
jerror_t JEventProcessor_TAC_Monitor::createHisto(unsigned trigBit,
		string histKey, string titlePrefix, string xTitle, int nBins,
		double xMin, double xMax) {
	static_assert(std::is_base_of<TH1, TH1_TYPE>::value,
	              "TH1_TYPE must be derived from TH1");
	stringstream histName;
	stringstream histTitle;
	histName << histKey << "_" << trigBit;
	histTitle << titlePrefix << trigBit;
	if (histoMap.count(histName.str()) == 0)
		histoMap[histName.str()] = map<unsigned, TH1*>();
	histoMap[histKey][trigBit] = new TH1_TYPE(histName.str().c_str(),
			histTitle.str().c_str(), nBins, xMin, xMax);
	histoMap[histKey][trigBit]->GetXaxis()->SetTitle(xTitle.c_str());
	return NOERROR;
}

// Create a 2D histogram of type TH2_TYPE and assign it to the histogram map based on the argument valeus
// provided in the function call.
template<typename TH2_TYPE>
jerror_t JEventProcessor_TAC_Monitor::createHisto(unsigned trigBit,
		string histKey, string titlePrefix, string xTitle, string yTitle,
		int nBinsX, double xMin, double xMax, int nBinsY, double yMin,
		double yMax) {
	static_assert(std::is_base_of<TH2, TH2_TYPE>::value,
	              "TH2_TYPE must be derived from TH2");
	stringstream histName;
	stringstream histTitle;
	histName << histKey << "_" << trigBit;
	histTitle << titlePrefix << trigBit;
	if (histoMap.count(histName.str()) == 0)
		histoMap[histName.str()] = map<unsigned, TH1*>();
	histoMap[histKey][trigBit] = new TH2_TYPE(histName.str().c_str(),
			histTitle.str().c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax);
	histoMap[histKey][trigBit]->GetXaxis()->SetTitle(xTitle.c_str());
	histoMap[histKey][trigBit]->GetYaxis()->SetTitle(yTitle.c_str());
	return NOERROR;
}


jerror_t JEventProcessor_TAC_Monitor::writeHistograms() {
	volatile WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());

	TDirectory* oldDir = gDirectory;
	TFile outFile( rootFileName.c_str(), "RECREATE" );
	outFile.cd();
	for( auto& histNameIter : histoMap ) {
//		auto& histName =  histNameIter.first;
		auto& histMapTrig = histNameIter.second;
		for( auto& histTrigIter : histMapTrig ) {
			auto histPointer = histTrigIter.second;
			histPointer->Write();
		}
	}
	outFile.Write();
	outFile.Close();
	oldDir->cd();

	return NOERROR;
}

// Return a pair giving the peak location (first) and the peak value (second)
inline pair<unsigned, unsigned> JEventProcessor_TAC_Monitor::getPeakLocationAndValue(
		const Df250WindowRawData* tacRawData) {
	// Find the maximum by going through the raw data and comparing samples
	pair<double, double> maxInfo(0, 0);
	if (tacRawData != nullptr) {
		maxInfo.first = std::distance(tacRawData->samples.begin(),
				std::max_element(tacRawData->samples.begin(),
						tacRawData->samples.end()));
		maxInfo.second = *std::max_element(tacRawData->samples.begin(),
				tacRawData->samples.end());
	}
	return maxInfo;
}

// Return the time where the signal crosses the threshold
inline unsigned JEventProcessor_TAC_Monitor::getPulseTime(
		const Df250WindowRawData* tacRawData, unsigned threshold) {
	if (tacRawData) {
		// Find the iterator where the value in the sample is greater than the threshold
		auto iterFound = std::find_if(tacRawData->samples.begin(),
				tacRawData->samples.end(),
				[&threshold]( const uint16_t& val) {return (val > threshold);});
		// If the iterator found was meaningful then return the sample number
		if (iterFound != tacRawData->samples.end()) {
			return distance(tacRawData->samples.begin(), iterFound);
		}
	}
	return 0;
}

// Fill Tagger-related histograms
template<typename TAG_TYPE, typename CounterID, typename TIME_CUT>
jerror_t JEventProcessor_TAC_Monitor::fillTaggerRelatedHistograms(
		vector<const TAG_TYPE*>& digiHitVector, uint32_t trigBit,
		string detComp, string tacMethod, double tacPeak, double tacTime,
		CounterID idFunctor, TIME_CUT timeCut) {
	for (auto digiHit : digiHitVector) {
		if (digiHit != nullptr) {
			double tagTime = double(digiHit->pulse_time) * fadc250DigiTimeScale;
			double detID = idFunctor(digiHit);
			bool match = timeCut(digiHit);

			volatile WriteLock rootRWLock(
					*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
			histoMap[detComp + "_ID"][trigBit]->Fill(detID);
			histoMap[detComp + "SigTime"][trigBit]->Fill(tagTime);
			histoMap["TACTIME" + tacMethod + "vs" + detComp + "TIME"][trigBit]->Fill(
					tagTime, tacTime);
			histoMap[detComp + "TIMEvs" + detComp + "ID"][trigBit]->Fill(detID,
					tagTime);
			if (match) {
				histoMap[detComp + "_ID_MATCHED" + tacMethod][trigBit]->Fill(
						detID);
				histoMap["TACAMP" + tacMethod + "vs" + detComp + "ID"][trigBit]->Fill(
						detID, tacPeak);
			}
		}
	}
	return NOERROR;
}

