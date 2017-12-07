// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /TAC/TAC_NHITS_1
// hnamepath: /TAC/TACFADCRAW_1
// hnamepath: /TAC/TACAmpPULSE_1
// hnamepath: /TAC/TACIntegral_1
// hnamepath: /TAC/TACTimePULSE_1
// hnamepath: /TAC/TAC_NTDCHITS_1
// hnamepath: /TAC/TAC_TDCTIME_1
// hnamepath: /TAC/TAC_TDCADCTIME_1
// hnamepath: /TAC/TAGHSigTime_1
// hnamepath: /TAC/TAGH_ID_MATCHEDPULSE_1
// hnamepath: /TAC/TAGMSigTime_1
// hnamepath: /TAC/TAGM_ID_MATCHEDPULSE_1

//
// e-mail: davidl@jlab.org
// e-mail: hovanes@jlab.org
//
/*
 * TAC_1D_hits.C
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

	histOrderVector.push_back( histoMap["TAC_NHITS_1"] = (TH1D*)gDirectory->FindObjectAny("TAC_NHITS_1") );
	histOrderVector.push_back( histoMap["TACFADCRAW_1"] = (TH1D*)gDirectory->FindObjectAny("TACFADCRAW_1") );
	histOrderVector.push_back( histoMap["TACAmpPULSE_1"] = (TH1D*)gDirectory->FindObjectAny("TACAmpPULSE_1") );
	histOrderVector.push_back( histoMap["TACIntegral_1"] = (TH1D*)gDirectory->FindObjectAny("TACIntegral_1") );
	histOrderVector.push_back( histoMap["TACTimePULSE_1"] = (TH1D*)gDirectory->FindObjectAny("TACTimePULSE_1") );
	histOrderVector.push_back( histoMap["TAC_NTDCHITS_1"] = (TH1D*)gDirectory->FindObjectAny("TAC_NTDCHITS_1") );
	histOrderVector.push_back( histoMap["TAC_TDCTIME_1"] = (TH1D*)gDirectory->FindObjectAny("TAC_TDCTIME_1") );
	histOrderVector.push_back( histoMap["TAC_TDCADCTIME_1"] = (TH1D*)gDirectory->FindObjectAny("TAC_TDCADCTIME_1") );
	histOrderVector.push_back( histoMap["TAGHSigTime_1"] = (TH1D*)gDirectory->FindObjectAny("TAGHSigTime_1") );
	histOrderVector.push_back( histoMap["TAGH_ID_MATCHEDPULSE_1"] = (TH1D*)gDirectory->FindObjectAny("TAGH_ID_MATCHEDPULSE_1") );
	histOrderVector.push_back( histoMap["TAGMSigTime_1"] = (TH1D*)gDirectory->FindObjectAny("TAGMSigTime_1") );
	histOrderVector.push_back( histoMap["TAGM_ID_MATCHEDPULSE_1"] = (TH1D*)gDirectory->FindObjectAny("TAGM_ID_MATCHEDPULSE_1") );


	drawOptMap["TAC_NHITS_1"] = "BAR";
	drawOptMap["TACFADCRAW_1"] = "L";
	drawOptMap["TACAmpPULSE_1"] = "L";
	drawOptMap["TACIntegral_1"] = "L";
	drawOptMap["TACTimePULSE_1"] = "BAR";
	drawOptMap["TAC_TDCTIME_1"] = "BAR";
	drawOptMap["TAC_NTDCHITS_1"] = "BAR";
	drawOptMap["TAC_TDCTIME_1"] = "BAR";
	drawOptMap["TAC_TDCADCTIME_1"] = "BAR";
	drawOptMap["TAGHSigTime_1"] = "BAR";
	drawOptMap["TAGH_ID_MATCHEDPULSE_1"] = "BAR";
	drawOptMap["TAGMSigTime_1"] = "BAR";
	drawOptMap["TAGM_ID_MATCHEDPULSE_1"] = "BAR";



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

	c1->Divide(4,3, 0.01, 0.01);

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
