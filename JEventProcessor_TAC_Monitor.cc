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
#include <DAQ/DEPICSvalue.h>
#include <TAC/DTACDigiHit.h>
#include <TAC/DTACHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>

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
double JEventProcessor_TAC_Monitor::timeCutValue_TAGH = 115.0;
// Timing cut width between the TAGH and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutWidth_TAGH = 30;

// Timing cut value between the TAGM and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutValue_TAGM = 90.0;
// Timing cut width between the TAGM and TAC coincidence in ns
double JEventProcessor_TAC_Monitor::timeCutWidth_TAGM = 20;

jerror_t JEventProcessor_TAC_Monitor::init(void) {
	WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
	TDirectory *mainDir = gDirectory;
	rootDir = gDirectory->mkdir("TAC");
	rootDir->cd();
	createHistograms();
	mainDir->cd();
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::brun(jana::JEventLoop* eventLoop,
		int32_t runnumber) {
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

	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned singleBit = 1 << trigBit;
		if ((singleBit & usefulTriggerBits) != 0) {
			this->fillRawDataHistograms(eventLoop, trigBit);
			this->fillPulseDataHitograms(eventLoop, trigBit);
		}
	}

//	const Df250WindowRawData* TACRaw = NULL;

	if( eventNumber % 200000 == 0 ) {
		this->writeHistograms();
	}

	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::fillRawDataHistograms(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {
	vector<const Df250WindowRawData*> rawDataVector;
	eventLoop->Get(rawDataVector);

	unsigned tacDataCounter = 0;
	// Pick the TAC raw data
	const Df250WindowRawData* tacRawData = nullptr;
	for (auto& rawData : rawDataVector) {
		if (rawData->rocid == 14 && rawData->slot == 20
				&& rawData->channel == 0) {
			tacRawData = rawData;
			tacDataCounter++;
		}
	}

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
		WriteLock rootRWLock(
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

	std::vector<const Df250PulsePedestal*> pulsePedestals;
	tacRawData->Get(pulsePedestals);
	if (pulsePedestals.size() > 1)
		return NOERROR;

	std::vector<const Df250PulseIntegral*> pulseIntegral;
	tacRawData->Get(pulseIntegral);
	if (pulseIntegral.size() > 1)
		return NOERROR;

	std::vector<const Df250PulseData*> pulseData;
	tacRawData->Get(pulseData);
	if (pulseData.size() > 1)
		return NOERROR;

	// Find the maximum by going through the raw data and comparing samples
	unsigned currentMaxValue = 0;
	for (auto& rawDataValue : tacRawData->samples) {
		if (rawDataValue > currentMaxValue) {
			currentMaxValue = rawDataValue;
		}
	}
	if (currentMaxValue >= maxPulseValue)
		currentMaxValue = overflowPulseValue;
	{
		stringstream histKey;
		histKey << "TACMax";
		WriteLock rootRWLock(
				*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
		histoMap[histKey.str()][trigBit]->Fill(currentMaxValue);
	}
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::fillPulseDataHitograms(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {
	vector<const DTACDigiHit*> tacDigiHitVector;
	eventLoop->Get(tacDigiHitVector);
	// pick only events with one digi object
	if (tacDigiHitVector.size() == 0 || tacDigiHitVector.size() > 1)
		return NOERROR;

	// Pick the only peak value from the TACDigiHit object
	const DTACDigiHit* tacDigiHit = tacDigiHitVector[0];
	double pulsePeak = double(tacDigiHit->getPulsePeak());
	double pulseTime = double(tacDigiHit->getPulseTime()) / 16.0;
//			cout << "Peak is at " << pulsePeak << " time is at " << pulseTime << endl;
	if (pulsePeak >= maxPulseValue)
		pulsePeak = overflowPulseValue;
	WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
	histoMap["TACAmp"][trigBit]->Fill(pulsePeak);
	histoMap["TACSigTime"][trigBit]->Fill(pulseTime);

	vector<const DTAGHDigiHit*> taghDigiHitVector;
	eventLoop->Get(taghDigiHitVector);
	for (auto taghDigiHit : taghDigiHitVector) {
		double taghTime = double(taghDigiHit->pulse_time) / 16.0;
		histoMap["TAGH_ID"][trigBit]->Fill(double(taghDigiHit->counter_id));
		histoMap["TAGHSigTime"][trigBit]->Fill(taghTime);
		histoMap["TACTIMEvsTAGHTIME"][trigBit]->Fill( taghTime, pulseTime );
		histoMap["TAGHTIMEvsTAGHID"][trigBit]->Fill( double(taghDigiHit->counter_id), taghTime );

		if( fabs( taghTime - timeCutValue_TAGH ) < timeCutWidth_TAGH ) {
			histoMap["TAGH_ID_MATCHED"][trigBit]->Fill(double(taghDigiHit->counter_id));
			histoMap["TACAMPvsTAGHID"][trigBit]->Fill(double(taghDigiHit->counter_id), pulsePeak );
		}
	}
	vector<const DTAGMDigiHit*> tagmDigiHitVector;
	eventLoop->Get(tagmDigiHitVector);
	for (auto tagmDigiHit : tagmDigiHitVector) {
		double tagmTime = double(tagmDigiHit->pulse_time) / 16.0;
		histoMap["TAGM_ID"][trigBit]->Fill(double(tagmDigiHit->column));
		histoMap["TAGMSigTime"][trigBit]->Fill(tagmTime);
		histoMap["TACTIMEvsTAGMTIME"][trigBit]->Fill( tagmTime, pulseTime );

		if( fabs( tagmTime - timeCutValue_TAGM ) < timeCutWidth_TAGM ) {
			histoMap["TAGM_ID_MATCHED"][trigBit]->Fill(double(tagmDigiHit->column));
			histoMap["TACAMPvsTAGMID"][trigBit]->Fill(double(tagmDigiHit->column), pulsePeak );
			histoMap["TAGMTIMEvsTAGMID"][trigBit]->Fill(  double(tagmDigiHit->column), tagmTime );
		}
	}


	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::erun(void) {
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
	this->writeHistograms();
	return NOERROR;
}

void JEventProcessor_TAC_Monitor::createHistograms() {
	cout << "Creating TAC histos" << endl;
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned trigPattern = 1 << trigBit;
		if (triggerIsUseful(trigPattern)) {
			// Create TAC FADc raw data
			createHisto<TH1D>(trigBit, "TACFADCRAW", "TAC FADC waveform ",
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
			// Create TAC amplitude histos
			createHisto<TH1D>(trigBit, "TACAmp",
					"TAC Signal Amplitude for Trigger ", "TAC Amplitude", 500,
					0., 5000.);
			// Create TAC amplitude histos for going through the data and picking the highest bin
			createHisto<TH1D>(trigBit, "TACMax",
					"TAC Signal Maximum from Raw for Trigger ", "TAC Amplitude",
					500, 0., 5000.);
			// Create TAC signal time histo
			createHisto<TH1D>(trigBit, "TACSigTime",
					"TAC Signal time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);
			// Create TAGH Hits detector ID
			createHisto<TH1D>(trigBit, "TAGH_ID",
					"TAGH Hits Detector ID for Trigger ", "Tagger Hodoscope Det. Number [#]",
					320, 0., 320.);
			// Create TAGH Hits detector ID
			createHisto<TH1D>(trigBit, "TAGH_ID_MATCHED",
					"Matched TAGH Hits Detector ID for Trigger ", "Tagger Hodoscope Det. Number [#]",
					320, 0., 320.);
			// Create TAGH signal time histo
			createHisto<TH1D>(trigBit, "TAGHSigTime",
					"TAGH Signal time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);
			// Create TAC time vs TAGH FADC time histo
			createHisto<TH2D>(trigBit, "TACTIMEvsTAGHTIME",
					"TAC time vs TAGH time ", "FlashADC peak time for TAGH (ns)", "FlashADC peak time for TAC (ns)",
					400, 0., 400., 400, 0., 400. );
			// Create TAC amplitude vs TAGH ID histo
			createHisto<TH2D>(trigBit, "TACAMPvsTAGHID",
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
			createHisto<TH1D>(trigBit, "TAGM_ID_MATCHED",
					"Matched TAGM Hits Detector ID for Trigger ", "Tagger Microscope Det. Number [#]",
					110, 0., 110.);
			// Create TAGM signal time histo
			createHisto<TH1D>(trigBit, "TAGMSigTime",
					"TAGM Signal time for Trigger ", "FlashADC peak time (ns)",
					400, 0., 400.);
			createHisto<TH2D>(trigBit, "TACTIMEvsTAGMTIME",
					"TAC time vs TAGM time ", "FlashADC peak time for TAGM (ns)", "FlashADC peak time for TAC (ns)",
					400, 0., 400., 400, 0., 400. );
			// Create TAC amplitude vs TAGM ID histo
			createHisto<TH2D>(trigBit, "TACAMPvsTAGMID",
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
	WriteLock rootRWLock(
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
