/*
 * CompressionTester.cc
 *
 *  Created on: Jun 12, 2017
 *      Author: hovanes
 */

#include "CompressionTester.h"

using namespace std;

std::mutex CompressionTester::fileAccessMutex;

void CompressionTester::writeData(const vector<uint16_t>& inputData) {
	std::vector<int> intVector(inputData.size());
	std::transform(inputData.begin(), inputData.end(), intVector.begin(),
			[](uint16_t u) {return static_cast<int>(u);});

	data::data decoder;

	std::vector<uint16_t> low, high;
	std::vector<char> encoded;
	std::vector<char> encodedLossy;

	std::vector<char> vec16(2 * intVector.size());
	for (unsigned p = 0; p < intVector.size(); p++) {
		uint16_t *ptr = reinterpret_cast<uint16_t*>(&vec16[p * 2]);
		uint16_t value = (uint16_t) intVector[p];
		*ptr = value;
	}

	decoder.decompose(intVector, low, high);
	decoder.encode(intVector, encoded);
	decoder.encodeLossy(intVector, encodedLossy);

	std::vector<int> vecSUB = decoder.getSubtracted(intVector);
	std::vector<int> vecREM = decoder.getReiman(vecSUB);

//    decoder.print(vecSUB);
//    decoder.print(vecREM);
//    decoder.print(intVector);
//    decoder.print(low);
//    decoder.print(high);

	std::lock_guard<std::mutex> guard(fileAccessMutex);

	losslessStream.write(&encoded[0], encoded.size());
	lossyStream.write(&encodedLossy[0], encodedLossy.size());
	rawStream.write(&vec16[0], vec16.size());

	for (unsigned iSample = 0; iSample < inputData.size(); iSample++) {
		asciiStream << std::setw(5) << inputData[iSample];
		if (iSample < inputData.size() - 1)
			asciiStream << " , ";
	}
	asciiStream << endl;

	return;
}
