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

uint32_t JEventProcessor_TAC_Monitor::triggerMask = 0x2;
unsigned JEventProcessor_TAC_Monitor::minPulseValue = 400;
unsigned JEventProcessor_TAC_Monitor::maxPulseValue = 4095;
unsigned JEventProcessor_TAC_Monitor::overflowPulseValue = 4500;

jerror_t JEventProcessor_TAC_Monitor::init(void) {
	TDirectory *mainDir = gDirectory;
	TDirectory *tacDir = gDirectory->mkdir("TAC");
	tacDir->cd();
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

	// Get First Trigger Type
	const DL1Trigger *trigWords = nullptr;
	uint32_t trigBits = 0;
	try {
		eventLoop->GetSingle( trigWords );
	} catch ( ... ) {
	};
	if ( trigWords ) {
		trigBits = trigWords->trig_mask;
		if( ( trigBits & triggerMask ) == 0 ) return NOERROR;
	} else {
		trigBits = 0;
		return NOERROR;
	}

	// Decide if to continue considering this event based on the trigger bit pattern
	if( ! triggerIsUseful( trigWords ) ) return NOERROR;
	trigBits = trigWords->trig_mask;



	this->fillRawDataHistograms( eventLoop, trigBits );
	this->fillPulseDataHitograms( eventLoop, trigBits );


//	const Df250WindowRawData* TACRaw = NULL;

	if( eventNumber % 200000 ) {
		this->writeHistograms();
	}

	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::fillRawDataHistograms(
		jana::JEventLoop* eventLoop, uint32_t trigBits) {
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
	int binNumber = 0;
	for (auto& rawDataValue : tacRawData->samples) {
		binNumber++;
		unsigned tacBit = 0xFFFF;
		histoMap["TACFADCRAW"][tacBit]->SetBinContent( binNumber, rawDataValue );
		histoMap["TACFADCRAW_ENTRIES"][tacBit]->SetBinContent( binNumber,  histoMap["TACFADCRAW_ENTRIES"][tacBit]->GetBinContent(binNumber) + 1.0 );
		histoMap["TACFADCRAW_SUM"][tacBit]->SetBinContent( binNumber,  histoMap["TACFADCRAW_SUM"][tacBit]->GetBinContent(binNumber) + rawDataValue );
		histoMap["TACFADCRAW_AVG"][tacBit]->Divide( histoMap["TACFADCRAW_SUM"][tacBit], histoMap["TACFADCRAW_ENTRIES"][tacBit], 1, 1 );
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
	if (currentMaxValue > maxPulseValue)
		currentMaxValue = overflowPulseValue;
	{
		stringstream histKey;
		histKey << "TACMax";
		japp->RootFillLock(this);
		histoMap[histKey.str()][trigBits]->Fill(currentMaxValue);
		japp->RootFillUnLock(this);
	}
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::fillPulseDataHitograms(
		jana::JEventLoop* eventLoop, uint32_t trigBits) {
	vector<const DTACDigiHit*> tacDigiHitVector;
	eventLoop->Get(tacDigiHitVector);
	// pick only events with one digi object
	if (tacDigiHitVector.size() > 1)
		return NOERROR;

	// Find the peak value from the TACDigiHit object
	if (tacDigiHitVector.size() > 0) {
		for (auto tacDigiHit : tacDigiHitVector) {
			double pulsePeak = tacDigiHit->getPulsePeak();
			if (pulsePeak > maxPulseValue)
				pulsePeak = overflowPulseValue;
			stringstream histKey;
			histKey << "TACAmp";
			japp->RootFillLock(this);
			histoMap[histKey.str()][trigBits]->Fill(pulsePeak);
			japp->RootFillUnLock(this);
		}
	}
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::erun(void) {
	return NOERROR;
}

jerror_t JEventProcessor_TAC_Monitor::fini(void) {
	this->writeHistograms();
	return NOERROR;
}



void JEventProcessor_TAC_Monitor::createHistograms() {
	cout << "Creating TAC histos" << endl;
	japp->RootWriteLock();

	// Create TAC FADc raw data
	{
		stringstream histKey;
		stringstream histName;
		stringstream histTitle;
		unsigned trigBit = 0xFFFF;
		histKey << "TACFADCRAW";
		histName << histKey.str();
		histTitle << "TAC FADC waveform " ;
		if( histoMap.count( histName.str() ) == 0 ) histoMap[histName.str()] = map<unsigned, TH1*>();
		histoMap[histKey.str()][trigBit] = new TH1F( histName.str().c_str(), histTitle.str().c_str(), 100, 0., 100. );
		histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle( "FlashADC sample number [#]" );
	}
	// Create TAC summed FADc raw data
	{
		stringstream histKey;
		stringstream histName;
		stringstream histTitle;
		unsigned trigBit = 0xFFFF;
		histKey << "TACFADCRAW_SUM";
		histName << histKey.str();
		histTitle << "Summed TAC FADC waveform" ;
		if( histoMap.count( histName.str() ) == 0 ) histoMap[histName.str()] = map<unsigned, TH1*>();
		histoMap[histKey.str()][trigBit] = new TH1F( histName.str().c_str(), histTitle.str().c_str(), 100, 0., 100. );
		histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle( "FlashADC sample number [#]" );
	}
	// Create TAC FADc raw data for entries
	{
		stringstream histKey;
		stringstream histName;
		stringstream histTitle;
		unsigned trigBit = 0xFFFF;
		histKey << "TACFADCRAW_ENTRIES";
		histName << histKey.str();
		histTitle << "Entries in TAC FADC waveform" ;
		if( histoMap.count( histName.str() ) == 0 ) histoMap[histName.str()] = map<unsigned, TH1*>();
		histoMap[histKey.str()][trigBit] = new TH1F( histName.str().c_str(), histTitle.str().c_str(), 100, 0., 100. );
		histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle( "FlashADC sample number [#]" );
	}

	// Create TAC averaged FADc raw data
	{
		stringstream histKey;
		stringstream histName;
		stringstream histTitle;
		unsigned trigBit = 0xFFFF;
		histKey << "TACFADCRAW_AVG";
		histName << histKey.str();
		histTitle << "Averaged TAC FADC waveform" ;
		if( histoMap.count( histName.str() ) == 0 ) histoMap[histName.str()] = map<unsigned, TH1*>();
		histoMap[histKey.str()][trigBit] = new TH1F( histName.str().c_str(), histTitle.str().c_str(), 100, 0., 100. );
		histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle( "FlashADC sample number [#]" );
	}

	// Create TAC amplitude histos
	for (unsigned trigBit = 0; trigBit < 16; trigBit++) {
		stringstream histKey;
		stringstream histName;
		stringstream histTitle;
		histKey << "TACAmp";
		histName << histKey.str() << trigBit;
		histTitle << "TAC Signal Amplitude for Trigger " << trigBit;
		cout << "Reserving histo with name " <<  histName.str() << " and trigger bit " << trigBit << endl;
		if( histoMap.count( histName.str() ) == 0 ) histoMap[histName.str()] = map<unsigned, TH1*>();
		histoMap[histKey.str()][trigBit] = new TH1D(histName.str().c_str(),
				histTitle.str().c_str(), 500, 0., 5000.);
		histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle("TAC Amplitude");
	}
	// Create TAC amplitude histos for going through the data and picking the highest bin
	for (unsigned trigBit = 0; trigBit < 16; trigBit++) {
		stringstream histKey;
		stringstream histName;
		stringstream histTitle;
		histKey << "TACMax";
		histName << histKey.str() << trigBit;
		histTitle << "TAC Signal Maximum from Raw for Trigger " << trigBit;
		cout << "Reserving histo with name " <<  histName.str() << " and trigger bit " << trigBit << endl;
		if( histoMap.count( histName.str() ) == 0 ) histoMap[histName.str()] = map<unsigned, TH1*>();
		histoMap[histKey.str()][trigBit] = new TH1D(histName.str().c_str(),
				histTitle.str().c_str(), 500, 0., 5000.);
		histoMap[histKey.str()][trigBit]->GetXaxis()->SetTitle("TAC Amplitude");
	}

	japp->RootUnLock();
}


jerror_t JEventProcessor_TAC_Monitor::writeHistograms() {
	TFile outFile( rootFileName.c_str(), "RECREATE" );
	outFile.cd();
	japp->RootReadLock();
	for( auto& histNameIter : histoMap ) {
//		auto& histName =  histNameIter.first;
		auto& histMapTrig = histNameIter.second;
		for( auto& histTrigIter : histMapTrig ) {
			auto histPointer = histTrigIter.second;
			histPointer->Write();
		}
	}
	japp->RootUnLock();

	outFile.Write();
	outFile.Close();

	return NOERROR;
}
