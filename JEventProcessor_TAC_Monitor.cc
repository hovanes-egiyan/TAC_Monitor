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
	// Force RootSpy lock to be set to the global RootReadWrite lock the moment this
	// plugin is loaded and this init function is called, unless it is non-zero.
	// This will interfere with all other plugins that use DRootSpy class object without
	// actual locking RootSpy plugin.
//	if (gROOTSPY_RW_LOCK == nullptr) {
//		if (dynamic_cast<DApplication*>(japp)) {
//			gROOTSPY_RW_LOCK =
//					dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock();
//			gROOTSPY->Initialize( (dynamic_cast<DApplication*>(japp))->GetRootReadWriteLock(), "<default>");
//		}
//	}
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

jerror_t JEventProcessor_TAC_Monitor::init(void) {
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
	if( dynamic_cast<DApplication*>(japp) == nullptr )
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

//	if( eventNumber % 200000 ) {
//		this->writeHistograms();
//	}

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
//	cout << "Window length is " << tacRawData->samples.size() << endl;

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
	if (tacDigiHitVector.size() > 1)
		return NOERROR;

	// Find the peak value from the TACDigiHit object
	if (tacDigiHitVector.size() > 0) {
		for (auto tacDigiHit : tacDigiHitVector) {
			double pulsePeak = tacDigiHit->getPulsePeak();
			if (pulsePeak >= maxPulseValue)
				pulsePeak = overflowPulseValue;
			stringstream histKey;
			histKey << "TACAmp";
			WriteLock rootRWLock(
					*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
			histoMap[histKey.str()][trigBit]->Fill(pulsePeak);
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
//	this->writeHistograms();
	return NOERROR;
}

void JEventProcessor_TAC_Monitor::createHistograms() {
	cout << "Creating TAC histos" << endl;
	WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned trigPattern = 1 << trigBit;
		if (triggerIsUseful(trigPattern)) {
			// Create TAC FADc raw data
			{
				stringstream histKey;
				stringstream histName;
				stringstream histTitle;
				histKey << "TACFADCRAW";
				histName << histKey.str();
				histTitle << "TAC FADC waveform ";
				if (histoMap.count(histName.str()) == 0)
					histoMap[histName.str()] = map<unsigned, TH1*>();
				histoMap[histKey.str()][trigBit] = new TH1F(
						histName.str().c_str(), histTitle.str().c_str(), 100,
						0., 100.);
				histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle(
						"FlashADC sample number [#]");
			}
			// Create TAC summed FADc raw data
			{
				stringstream histKey;
				stringstream histName;
				stringstream histTitle;
				histKey << "TACFADCRAW_SUM";
				histName << histKey.str();
				histTitle << "Summed TAC FADC waveform" << trigBit;
				if (histoMap.count(histName.str()) == 0)
					histoMap[histName.str()] = map<unsigned, TH1*>();
				histoMap[histKey.str()][trigBit] = new TH1F(
						histName.str().c_str(), histTitle.str().c_str(), 100,
						0., 100.);
				histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle(
						"FlashADC sample number [#]");
			}
			// Create TAC FADc raw data for entries
			{
				stringstream histKey;
				stringstream histName;
				stringstream histTitle;
				histKey << "TACFADCRAW_ENTRIES";
				histName << histKey.str();
				histTitle << "Entries in TAC FADC waveform bit " << trigBit;
				if (histoMap.count(histName.str()) == 0)
					histoMap[histName.str()] = map<unsigned, TH1*>();
				histoMap[histKey.str()][trigBit] = new TH1F(
						histName.str().c_str(), histTitle.str().c_str(), 100,
						0., 100.);
				histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle(
						"FlashADC sample number [#]");
			}

			// Create TAC averaged FADc raw data
			{
				stringstream histKey;
				stringstream histName;
				stringstream histTitle;
				histKey << "TACFADCRAW_AVG";
				histName << histKey.str();
				histTitle << "Averaged TAC FADC waveform" << trigBit;
				if (histoMap.count(histName.str()) == 0)
					histoMap[histName.str()] = map<unsigned, TH1*>();
				histoMap[histKey.str()][trigBit] = new TH1F(
						histName.str().c_str(), histTitle.str().c_str(), 100,
						0., 100.);
				histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle(
						"FlashADC sample number [#]");
			}

		}
	}
	// Create TAC amplitude histos
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned trigPattern = 1 << trigBit;
		if (triggerIsUseful(trigPattern)) {
			stringstream histKey;
			stringstream histName;
			stringstream histTitle;
			histKey << "TACAmp";
			histName << histKey.str() << trigBit;
			histTitle << "TAC Signal Amplitude for Trigger " << trigBit;
			cout << "Reserving histo with name " << histName.str()
					<< " and trigger bit " << trigBit << endl;
			if (histoMap.count(histName.str()) == 0)
				histoMap[histName.str()] = map<unsigned, TH1*>();
			histoMap[histKey.str()][trigBit] = new TH1D(histName.str().c_str(),
					histTitle.str().c_str(), 500, 0., 5000.);
			histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle(
					"TAC Amplitude");
		}
	}
	// Create TAC amplitude histos for going through the data and picking the highest bin
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned trigPattern = 1 << trigBit;
		if (triggerIsUseful(trigPattern)) {
			stringstream histKey;
			stringstream histName;
			stringstream histTitle;
			histKey << "TACMax";
			histName << histKey.str() << trigBit;
			histTitle << "TAC Signal Maximum from Raw for Trigger " << trigBit;
			cout << "Reserving histo with name " << histName.str()
					<< " and trigger bit " << trigBit << endl;
			if (histoMap.count(histName.str()) == 0)
				histoMap[histName.str()] = map<unsigned, TH1*>();
			histoMap[histKey.str()][trigBit] = new TH1D(histName.str().c_str(),
					histTitle.str().c_str(), 500, 0., 5000.);
			histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle(
					"TAC Amplitude");
		}
	}
}

//jerror_t JEventProcessor_TAC_Monitor::writeHistograms() {
//	return NOERROR;
//	japp->RootReadLock();
//
//	TDirectory* oldDir = gDirectory;
//	TFile outFile( rootFileName.c_str(), "RECREATE" );
//	outFile.cd();
//	for( auto& histNameIter : histoMap ) {
////		auto& histName =  histNameIter.first;
//		auto& histMapTrig = histNameIter.second;
//		for( auto& histTrigIter : histMapTrig ) {
//			auto histPointer = histTrigIter.second;
//			histPointer->Write();
//		}
//	}
//
//	outFile.Write();
//	outFile.Close();
//	oldDir->cd();
//
//	japp->RootUnLock();
//
//	return NOERROR;
//}
