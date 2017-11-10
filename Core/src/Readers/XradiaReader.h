/*
  This is class hold the data for Xradia instruments (Zeiss CT). This uses pole library for reading txrm files(Encoded in OLE)
  Author: Mr. Srikanth Nagella
*/

#ifndef XRADIAREADER_H
#define XRADIAREADER_H
#include "CCPiDefines.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

#include "CCPiUserApplicationInterface.h"
namespace CCPi
{
	typedef boost::multi_array<float,3> float_pixel;
	typedef boost::multi_array<short,3> short_pixel;

	class CCPI_EXPORT XradiaReader
	{
	public:
		XradiaReader(std::string filename, CCPiUserApplicationInterface *uai=0);
		~XradiaReader();
		enum DataType { FLOAT, UINT_16};
		int getImageWidth(); //Returns the image width
		int getImageHeight(); //Returns the image height
		DataType getDataType(); //Return the data type of the image data
		int getNumberOfImages();//returns the number of Images
		std::vector<float> getAngles();//Return the Angles
		std::vector<float> getXShifts();//Return the xShift values
		std::vector<float> getYShifts();//Return the yShift values
		boost::shared_ptr< short_pixel > getImageDataInShort(); //Return the projection data in short datatype
		boost::shared_ptr< float_pixel > getImageDataInFloat(); //Returns the projection data in float
		float getDetectorToObject(){return detToObject;}
		float getSourceToObject(){return srcToObject;}
 	private:
		std::string filename; // file name to use to read the data
		int imageWidth; //Image Width
		int imageHeight; //Image Height
		int numberOfImages; // Number of Images/Angles
		std::vector<float> angles; // Array of angles
		std::vector<float> xSamplePosition; //Sample position in X
		std::vector<float> ySamplePosition; //Sample position in Y
		std::vector<float> zSamplePosition; //Sample position in Z
		std::vector<float> xShifts; //X Shift
		std::vector<float> yShifts;//Y Shift
		DataType type;
		boost::shared_ptr< float_pixel > floatImageData;
		boost::shared_ptr< short_pixel > shortImageData;
		float srcToObject;
		float detToObject;
		bool localUAI;
		CCPiUserApplicationInterface *uai;
		//Private methods
		void extractMetaDataAll();
		int extractMetaDataInt(std::string streamName);
		std::vector<float> extractMetaDataFloatArray(std::string streamName,int numberOfValues);
		void extractImageData();
	};
}

#endif // XRADIAREADER_H
