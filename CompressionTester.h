/*
 * CompressionTester.h
 *
 *  Created on: Jun 12, 2017
 *      Author: hovanes
 */

#ifndef COMPRESSIONTESTER_H_
#define COMPRESSIONTESTER_H_

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <ios>
#include <iosfwd>
#include <stdint.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <thread>
#include <mutex>
#include <iomanip>
#include "data.h"

class CompressionTester {
protected:
	std::string rawFileName;
	std::string losslessFileName;
	std::string lossyFileName;
	std::string asciiFileName;

	std::ofstream rawStream;
	std::ofstream losslessStream;
	std::ofstream lossyStream;

	std::ofstream asciiStream;

	static std::mutex fileAccessMutex;

public:
	CompressionTester(std::string prefix = "tacCompression") :
			rawFileName(prefix + "_Raw.bin"), losslessFileName(
					prefix + "_Lossless.bin"), lossyFileName(
					prefix + "_Lossy.bin"), asciiFileName(prefix+"_Samples.txt") {
		rawStream.open(rawFileName, std::ios::out | std::ios::binary);
		losslessStream.open(losslessFileName, std::ios::out | std::ios::binary);
		lossyStream.open(lossyFileName, std::ios::out | std::ios::binary);
		asciiStream.open(asciiFileName, std::ios::out);
	}
	virtual ~CompressionTester() {
		std::cout << "Closing files" << std::endl;
		rawStream.close();
		losslessStream.close();
		lossyStream.close();
		asciiStream.close();
	}

	virtual void writeData(const std::vector<uint16_t>& inputData);

	static const std::mutex& getFileAccessMutex() {
		return fileAccessMutex;
	}

	const std::string& getLosslessFileName() const {
		return losslessFileName;
	}

	void setLosslessFileName(const std::string& losslessFileName) {
		this->losslessFileName = losslessFileName;
	}

	const std::ofstream& getLosslessStream() const {
		return losslessStream;
	}

	const std::string& getLossyFileName() const {
		return lossyFileName;
	}

	void setLossyFileName(const std::string& lossyFileName) {
		this->lossyFileName = lossyFileName;
	}

	const std::ofstream& getLossyStream() const {
		return lossyStream;
	}

	const std::string& getRawFileName() const {
		return rawFileName;
	}

	void setRawFileName(const std::string& rawFileName) {
		this->rawFileName = rawFileName;
	}

	const std::ofstream& getRawStream() const {
		return rawStream;
	}

	const std::string& getAsciiFileName() const {
		return asciiFileName;
	}

	void setAsciiFileName(const std::string& asciiFileName) {
		this->asciiFileName = asciiFileName;
	}

	const std::ofstream& getAsciiStream() const {
		return asciiStream;
	}
};

#endif /* COMPRESSIONTESTER_H_ */
