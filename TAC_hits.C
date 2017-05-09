// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /TAC/TACFADCRAW_1
// hnamepath: /TAC/TAC_NHITS_1
// hnamepath: /TAC/TACAmpPULSE_1
// hnamepath: /TAC/TAC_TDCTIME_1
// hnamepath: /TAC/TACTIMEPULSEvsTAGHTIME_1
// hnamepath: /TAC/TACTIMEPULSEvsTAGMTIME_1
// hnamepath: /TAC/TACAMPPULSEvsTAGHID_1
// hnamepath: /TAC/TACAMPPULSEvsTAGMID_1
//
// e-mail: davidl@jlab.org
// e-mail: hovanes@jlab.org
//

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

	histOrderVector.push_back( histoMap["TACFADCRAW_1"] = (TH1D*)gDirectory->FindObjectAny("TACFADCRAW_1") );
	histOrderVector.push_back( histoMap["TAC_NHITS_1"] = (TH1D*)gDirectory->FindObjectAny("TAC_NHITS_1") );
	histOrderVector.push_back( histoMap["TACAmpPULSE_1"] = (TH1D*)gDirectory->FindObjectAny("TACAmpPULSE_1") );
	histOrderVector.push_back( histoMap["TAC_TDCTIME_1"] = (TH1D*)gDirectory->FindObjectAny("TAC_TDCTIME_1") );

	histOrderVector.push_back( histoMap["TACTIMEPULSEvsTAGHTIME_1"] = (TH2D*)gDirectory->FindObjectAny("TACTIMEPULSEvsTAGHTIME_1") );
	histOrderVector.push_back( histoMap["TACTIMEPULSEvsTAGMTIME_1"] = (TH2D*)gDirectory->FindObjectAny("TACTIMEPULSEvsTAGMTIME_1") );
	histOrderVector.push_back( histoMap["TACAMPPULSEvsTAGHID_1"] = (TH2D*)gDirectory->FindObjectAny("TACAMPPULSEvsTAGHID_1") );
	histOrderVector.push_back( histoMap["TACAMPPULSEvsTAGMID_1"] = (TH2D*)gDirectory->FindObjectAny("TACAMPPULSEvsTAGMID_1" ) );

	drawOptMap["TACFADCRAW_1"] = "L";
	drawOptMap["TAC_NHITS_1"] = "BAR";
	drawOptMap["TACAmpPULSE_1"] = "L";
	drawOptMap["TAC_TDCTIME_1"] = "BAR";

	drawOptMap["TACTIMEPULSEvsTAGHTIME_1"] = "COLZ";
	drawOptMap["TACTIMEPULSEvsTAGMTIME_1"] = "COLZ";
	drawOptMap["TACAMPPULSEvsTAGHID_1"] = "COLZ";
	drawOptMap["TACAMPPULSEvsTAGMID_1"] = "COLZ";

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

	// Just for testing
	if(gPad == nullptr){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}

	if(!gPad) {savedir->cd(); return;}

	TCanvas *c1 = gPad->GetCanvas();
	if( c1 == nullptr ) return;
	c1->cd(0);
	c1->Clear();

	c1->Divide(4,2, 0.01, 0.01);

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
