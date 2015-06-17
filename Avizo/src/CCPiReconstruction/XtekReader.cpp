#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "XtekReader.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include "exlibtiff.h"
#include "tiffio.h"
#include <hxcore/HxMessage.h>
#define Message theMsg->stream()

//Constructor
XtekReader::XtekReader(std::string filename)
{
	test = true;
	Message<<"Starting to read the .xtekct file"<<std::endl;
	if(readXtekCTFile(filename))
	{
		Message<<"Completed reading .xtekct file"<<std::endl;
		//Split and get the file path and append the _ctdata.txt to read angles data
		boost::filesystem::path xtekfile(filename);
		boost::filesystem::path ctname("_ctdata.txt");
		boost::filesystem::path full_ctname = xtekfile.parent_path() / ctname;
		if(!readCTDataFile(full_ctname.string()))
		{
			Message<<"Error reading the _ctdata.txt file"<<std::endl;
			//Try ang file
			boost::filesystem::path full_ang = xtekfile.replace_extension("ang");
			Message<<"Reading ang file instead"<<full_ang.string()<<std::endl;
			if(!readAngFile(full_ang.string()))
				return;
		}
		Message<<"Completed reading the _ctdata.txt file"<<std::endl;
		//allocate memory for the images
		unsigned long long width = detectorPixels[0];
		unsigned long long height = detectorPixels[1];
		std::cout<<projections<<"x"<<width<<"x"<<height<<std::endl;
		imageData = new float[projections * width * height];
		Message<<"Allocated the image data ("<<projections*width*height<<")"<<std::endl;
		boost::filesystem::path imgname(name); //prefix of the images
		boost::filesystem::path prefix_image = xtekfile.parent_path() / imgname;
		std::cout<<"Tiff file prefix: "<<prefix_image.string()<<std::endl;
		std::cout<<"File prefix: "<<imgname.string()<<std::endl;
		std::cout<<"File name: "<<name<<std::endl;
		//read the images
		readImageFiles(prefix_image.string(), projections, imageData);
	}
}

//Destructor
XtekReader::~XtekReader()
{
}

bool XtekReader::readXtekCTFile(std::string filename)
{
	//Open the xxx.xtetct file
	std::ifstream ctfile(filename.c_str());
	if(!ctfile.is_open()) return false;
	//read the header
	std::string line;
	getline(ctfile, line);
	boost::algorithm::trim(line);
	if(line.compare("[XTekCT]")!=0) return false; // return on invalid header
	while( getline(ctfile, line))
	{
		boost::algorithm::trim(line);
		std::vector<std::string> keyvalue;
		boost::algorithm::split(keyvalue, line, boost::is_any_of("="));
		if(keyvalue[0].compare("Name")==0)
		{
			name = keyvalue[1];
		}else if(keyvalue[0].compare("VoxelsX")==0)
		{
			voxels[0] = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("VoxelsY")==0)
		{
			voxels[1] = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("VoxelsZ")==0)
		{
			voxels[2] = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("VoxelSizeX")==0)
		{
			voxelSize[0] = strtod(keyvalue[1].c_str(),NULL);
		}else if(keyvalue[0].compare("VoxelSizeY")==0)
		{
			voxelSize[1] = strtod(keyvalue[1].c_str(),NULL);
		}else if(keyvalue[0].compare("VoxelSizeZ")==0)
		{
			voxelSize[2] = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("OffsetX")==0)
		{
			offset[0] = strtod(keyvalue[1].c_str(),NULL);
		}else if(keyvalue[0].compare("OffsetY")==0)
		{
			offset[1] = strtod(keyvalue[1].c_str(), NULL);
		}else if (keyvalue[0].compare("OffsetZ")==0)
		{
			offset[2] = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("SrcToObject")==0)
		{
			srcToObject = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("SrcToDetector")==0)
		{
			srcToDetector = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("MaskRadius")==0)
		{
			maskRadius = strtod(keyvalue[1].c_str(),NULL);
		}else if(keyvalue[0].compare("DetectorPixelsX")==0)
		{
			detectorPixels[0]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("DetectorPixelsY")==0)
		{
			detectorPixels[1]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("DetectorPixelSizeX")==0)
		{
			detectorPixelSize[0]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("DetectorPixelSizeY")==0)
		{
			detectorPixelSize[1]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("DetectorOffsetX")==0)
		{
			detectorOffset[0]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("DetectorOffsetY")==0)
		{
			detectorOffset[1]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CentreOfRotation")==0)
		{
			centreOfRotation[0]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CentreOfRotationTop")==0)
		{
			centreOfRotation[1]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CentreOfRotationBottom")==0)
		{
			centreOfRotation[2]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("WhiteLevel")==0)
		{
			whiteLevel= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Scattering")==0)
		{
			scattering= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CoefX4")==0)
		{
			coefX[4]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CoefX3")==0)
		{
			coefX[3]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CoefX2")==0)
		{
			coefX[2]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CoefX1")==0)
		{
			coefX[1]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("CoefX0")==0)
		{
			coefX[0]= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Scale")==0)
		{
			scale= strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("RegionStartX")==0)
		{
			regionStart[0] = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("RegionStartY")==0)
		{
			regionStart[1] = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("RegionPixelsX")==0)
		{
			regionPixels[0] = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("RegionPixelsY")==0)
		{
			regionPixels[1] = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Projections")==0)
		{
			projections = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("InitialAngle")==0)
		{
			initialAngle = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("AngularStep")==0)
		{
			angularStep = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("FilterType")==0)
		{
			filterType = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("CutOffFrequency")==0)
		{
			cutOffFrequency = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Exponent")==0)
		{
			exponent = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Normalisation")==0)
		{
			normalisation = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("InterpolationType")==0)
		{
			interpolationType = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("Scaling")==0)
		{
			scaling = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("OutputType")==0)
		{
			outputType = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("ImportConversion")==0)
		{
			importConversion = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("TIFFScaling")==0)
		{
			tiffScaling = atoi(keyvalue[1].c_str());
		}else if(keyvalue[0].compare("XraykV")==0)
		{
			xrayKV = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("XrayuA")==0)
		{
			xrayuA = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Filter_ThicknessMM")==0)
		{
			filterThicknessMM = strtod(keyvalue[1].c_str(), NULL);
		}else if(keyvalue[0].compare("Filter_Material")==0)
		{
			filterMaterial = keyvalue[1];
		}else if(keyvalue[0].compare("Shuttling")==0)
		{
			boost::algorithm::to_lower(keyvalue[1]);
			std::istringstream is(keyvalue[1]);
			is >> std::boolalpha >> shuttling;
		}else if(keyvalue[0].compare("InputSeparator") == 0) 
		{
			inputSeperator = keyvalue[1];
		}
	}
	return true;
}

bool XtekReader::readCTDataFile(std::string filename)
{
	//Open the _ctdata.txt file
	std::ifstream ctfile(filename.c_str());
	if(!ctfile.is_open()) return false;
	//read the header
	std::string line;
	getline(ctfile, line);
	boost::algorithm::trim(line);
	if(line.compare("Projs	Frames	Exposure(ms)")!=0) return false; // return on invalid header
	//Read number of projections
	getline(ctfile, line);
	boost::algorithm::trim(line);
	std::vector<std::string> values;
	boost::algorithm::split(values, line, boost::is_any_of("\t "), boost::token_compress_on);
	if (atoi(values[0].c_str())!= projections) return false; // number of projections don't match
	framesPerProjection = atoi(values[1].c_str());
	exposureTime = strtod(values[2].c_str(), NULL);

	//Read the next header
	getline(ctfile, line);
	//Allocate angles memory
	angles = new float[projections];
	angleTime = new double[projections];
	while( getline(ctfile, line))
	{
		boost::algorithm::trim(line);
		std::vector<std::string> values;
		boost::algorithm::split(values, line, boost::is_any_of("\t "), boost::token_compress_on);
		int proj = atoi(values[0].c_str()) - 1;
		if(proj > projections) return false; // projection index cannot be greater than number projections
		angles[proj] = strtof(values[1].c_str(), NULL); 
		angleTime[proj] = strtod(values[2].c_str(), NULL);
	}
	return true;
}

bool XtekReader::readAngFile(std::string filename)
{
	//Open the xxx.ang file
	std::ifstream ctfile(filename.c_str());
	if(!ctfile.is_open()) return false;
	//read the header
	std::string line;
	getline(ctfile, line);
	boost::algorithm::trim(line);
	if(line.compare("Proj	Angle(deg)")!=0) return false; // return on invalid header
	angles = new float[projections];
	while( getline(ctfile, line))
	{
		boost::algorithm::trim(line);
		std::vector<std::string> values;
		boost::algorithm::split(values, line, boost::is_any_of("\t "), boost::token_compress_on);
		int proj = atoi(values[0].c_str());
		if(proj > projections) return false; // projection index cannot be greater than number projections
		angles[projections - proj] = strtof(values[1].c_str(), NULL); 
	}
	return true;
}

bool XtekReader::readImageFiles(std::string filePrefix, int numberOfAngles, float *data)
{
	unsigned long long width,height;
	width = detectorPixels[0];
	height = detectorPixels[1];
	for (unsigned long long i=0;i<numberOfAngles;i++) {
		std::stringstream fullpath;
		fullpath << boost::format("%s%s%|04d|.tif")%filePrefix %inputSeperator %(i+1);
		bool ok = readTiffFile(fullpath.str(), width, height, data+i*width*height);
		if(!ok)
		{
			std::cout<<"Error reading file: "<<fullpath.str()<<std::endl;
			return false;
		}
	}
	return true;
}

bool XtekReader::readTiffFile(std::string filename, int width, int height, float *data)
{
	TIFF *tif = TIFFOpen(filename.c_str(), "r");
	if (tif == 0) {
		std::cout<<"Error opening the file"<<std::endl;
		return false;
	} else {
		if (TIFFIsTiled(tif)) {
			std::cout<<"Error: Doesn't support tiled tiff"<<std::endl;
			return false; // Doesn't support tiled tiff
		} else {
			uint16 bps;
			if (TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps) == 0) {
				std::cout<<"Error: reading bits per sample"<<std::endl;
				return false; //TIFF error reading bits per sample
			} else if (bps != 16) {
				std::cout<<"Error: TIFF is not 16bit data"<<std::endl;
				return false; //TIFF is not 16bit data
			} else {
				uint16 photo;
				if (TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photo) == 0) {
					std::cout<<"Error: getting TIFF photometric info"<<std::endl;
					return false; //Error getting TIFF photometric info
				} else if (photo != PHOTOMETRIC_MINISBLACK) {
					std::cout<<"Error: TIFF photometric type not supported by XTek Reader"<<std::endl;
					return false;//TIFF photometric type not supported by XTek reader
				} else {
					if (TIFFIsMSB2LSB(tif) == 0) {
						std::cout<<"Error: TIFF is not MSB to LSB"<<std::endl; 
						return false;//TIFF is not MSB to LSB
					} else {
						// I'm assuming orientation would be 1, (0,0) in top left corner
						// but XTek tiffs don't usually work with TIFFTAG_ORIENTATION
						uint32 w;
						if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w) == 0) {
							std::cout<<"Error: getting TIFF width"<<std::endl;
							return false;//Error getting TIFF width
						} else if ((int)w != width) {
							std::cout<<"Error: Width don't match"<<std::endl;
							return false; //Image width mismatch
						} else {
							tsize_t strip_size = TIFFStripSize(tif);
							int nstrips = TIFFNumberOfStrips(tif);
							char *buf = (char *)_TIFFmalloc(strip_size * nstrips);
							int offset = 0;
							int result;
							for (int count = 0; (count < nstrips); count++) {
								if ((result = TIFFReadEncodedStrip(tif, count, buf + offset,
									strip_size)) == -1) {
										std::cout<<"Error: read error in file"<<std::endl;
										return false;//Read error in file
								} else
									offset += result;
							}
							int extra = offset % 2;
							offset /= 2;
							if (extra == 1) {
								std::cout<<"Error: Odd number of bytes for 16bit image"<<std::endl;
								return false;//Odd number of bytes for 16 bit image
							} else if (offset != width * height) {
								std::cout<<"Error: Image size mismatch"<<std::endl;
								return false;//Image size mismatch
							} else {
								//std::memcpy(data, buf, sizeof(uint16)*width*height);
								uint16 *b = (uint16 *)buf;
								for(unsigned long i=0;i<width;i++)
									for(unsigned long j=0;j<height;j++)
										data[i*height+j] = float(b[i*height+j]);		
							}
							_TIFFfree(buf);
						}
					}
				}
			}
		}
		TIFFClose(tif);
	}
	return true;
}
