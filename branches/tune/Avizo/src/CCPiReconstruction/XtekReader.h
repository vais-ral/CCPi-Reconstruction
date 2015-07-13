/**
 * This is a Xtek file format reader
 * Author : Mr. Srikanth Nagella
 */
#ifndef XTEKREADER_H
#define XTEKREADER_H
#include <string>

class XtekReader
{
private:
	//.xtekct parameters
	std::string name;
	std::string inputSeperator;
	int	voxels[3];
	double voxelSize[3];
	double offset[3];
	double srcToObject;
	double srcToDetector;
	double maskRadius;
	double detectorPixels[2];
	double detectorPixelSize[2];
	double detectorOffset[2];
	double centreOfRotation[2]; //0- center, 1 -- top, 2 -- bottom
	double whiteLevel;
	double scattering;
	double coefX[5]; 
	double scale;
	int regionStart[2];
	int regionPixels[2];
	int projections;
	double initialAngle;
	double angularStep;
	int filterType;
	double cutOffFrequency;
	double exponent;
	double normalisation;
	int interpolationType;
	double scaling;
	int outputType;
	int importConversion;
	int tiffScaling;
	double xrayKV;
	double xrayuA;
	double filterThicknessMM;
	std::string filterMaterial;
	bool shuttling;
	//_ctdata.txt parameters
	int framesPerProjection;
	double exposureTime;
	float *angles;
	double *angleTime;
	//Image data
	float* imageData;

	//Read xtekct file
	bool readXtekCTFile(std::string filename);
	bool readCTDataFile(std::string filename);
	bool readAngFile(std::string filename);
	bool readImageFiles(std::string filePrefix, int numberOfAngles, float *data);
	bool readTiffFile(std::string filename,int width, int height, float *data);

	//test
	bool test;
public:
	//Constructor
	XtekReader(std::string filename);
	//Destructor
	~XtekReader();
	//return data
	float* getImageData(){return imageData;}
	int getNumberOfProjections(){return projections;}
	int getImageWidth(){return detectorPixels[0];}
	int getImageHeight(){return detectorPixels[1];}
	int getFramesPerProjection(){return framesPerProjection;}
	double getExposureTime(){return exposureTime;}
	float* getAngles(){return angles;}
	double* getVoxelSize(){return voxelSize;}
	double* getOffset(){return offset;}
	double getSourceToObject(){return srcToObject;}
	double getSourceToDetector(){return srcToDetector;}
	double* getDetectorPixelSize(){return detectorPixelSize;}
	double getInitialAngle(){return initialAngle;}
	std::string getName(){return name;}
	double getWhiteLevel(){return whiteLevel;}
	double getScattering(){return scattering;}
	double* getCoefX(){return coefX;}
	double getScale(){return scale;}
	double getMaskRadius(){return maskRadius;}
};
#endif
