// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /TAC/TACTIMEPULSEvsTAGHTIME_1
// hnamepath: /TAC/TACAMPPULSEvsTAGHID_1
// hnamepath: /TAC/TAGHTIMEvsTAGHID_1
// hnamepath: /TAC/TACTIMEPULSEvsTAGMTIME_1
// hnamepath: /TAC/TACAMPPULSEvsTAGMID_1
// hnamepath: /TAC/TAGMTIMEvsTAGMID_1
//
// e-mail: davidl@jlab.org
// e-mail: hovanes@jlab.org
//

/*
 * TAC_2D_hist.C
 *
 *  Created on: Dec 6, 2017
 *      Author: Hovanes Egiyan
 */

#include <iostream>
#include <map>

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TPad.h"

{
//void main() {

	std::map<std::string,TH1*> histoMap;
	std::map<std::string, TLegend*> legendMap;
	std::map<std::string, std::string> drawOptMap;

	std::vector<TH1*> histOrderVector;

	TDirectory *savedir = gDirectory;
	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("TAC");
//	std::cout << "Directory is at " << dir << std::endl;
	if(dir) dir->cd();

	histOrderVector.push_back( histoMap["TACTIMEPULSEvsTAGHTIME_1"] = (TH2D*)gDirectory->FindObjectAny("TACTIMEPULSEvsTAGHTIME_1") );
	histOrderVector.push_back( histoMap["TACAMPPULSEvsTAGHID_1"] = (TH2D*)gDirectory->FindObjectAny("TACAMPPULSEvsTAGHID_1") );
	histOrderVector.push_back( histoMap["TAGHTIMEvsTAGHID_1"] = (TH2D*)gDirectory->FindObjectAny("TAGHTIMEvsTAGHID_1") );
	histOrderVector.push_back( histoMap["TACTIMEPULSEvsTAGMTIME_1"] = (TH2D*)gDirectory->FindObjectAny("TACTIMEPULSEvsTAGMTIME_1") );
	histOrderVector.push_back( histoMap["TACAMPPULSEvsTAGMID_1"] = (TH2D*)gDirectory->FindObjectAny("TACAMPPULSEvsTAGMID_1" ) );
	histOrderVector.push_back( histoMap["TAGMTIMEvsTAGMID_1"] = (TH2D*)gDirectory->FindObjectAny("TAGMTIMEvsTAGMID_1") );

	drawOptMap["TACTIMEPULSEvsTAGHTIME_1"] = "COLZ";
	drawOptMap["TACAMPPULSEvsTAGHID_1"] = "COLZ";
	drawOptMap["TAGHTIMEvsTAGHID_1"] = "COLZ";
	drawOptMap["TACTIMEPULSEvsTAGMTIME_1"] = "COLZ";
	drawOptMap["TACAMPPULSEvsTAGMID_1"] = "COLZ";
	drawOptMap["TAGMTIMEvsTAGMID_1"] = "COLZ";

	for (auto& histoPair : histoMap) {
		auto& histoName = histoPair.first;
		auto& histoPtr = histoPair.second;
		if (histoPtr) {
			legendMap[histoName] = new TLegend(0.5, 0.8, 0.9, 0.9);
			if (dynamic_cast<TH2D*>(histoPtr) == nullptr) {
				histoPtr->SetBarWidth(0.5);
				histoPtr->SetBarOffset(0);
				histoPtr->SetFillColor(kGreen);
			}
			histoPtr->SetStats(0);
			histoPtr->SetTitleSize(0.05, "X");
			histoPtr->GetXaxis()->CenterTitle();
			histoPtr->SetTitleSize(0.05, "Y");
			histoPtr->GetYaxis()->CenterTitle();
			if (dynamic_cast<TH2D*>(histoPtr) == nullptr)
				histoPtr->GetYaxis()->SetRangeUser(0.0, histoPtr->GetMaximum());
			legendMap[histoName]->AddEntry(histoPtr, histoPtr->GetTitle(), "f");
		}
	}

	if(!gPad) {savedir->cd(); return;}

	TCanvas *c1 = gPad->GetCanvas();
	if( c1 == nullptr ) return;
	c1->cd(0);
	c1->Clear();

	c1->Divide(3,2, 0.01, 0.01);

	unsigned iPad = 0;
	for (auto histoPtr : histOrderVector) {
		if (histoPtr) {
			c1->cd(++iPad);
			gPad->SetTicks();
			gPad->SetGridx();
			gPad->SetGridy();
			gPad->SetLogz();
			if (histoPtr) {
				histoPtr->DrawCopy(drawOptMap[histoPtr->GetName()].c_str());
//			legendMap[histoPtr->GetName()]->Draw();
			}
			gPad->Update();
		}
	}

	c1->cd(0);
	c1->Update();

}
