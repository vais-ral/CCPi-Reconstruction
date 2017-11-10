#include "XradiaReader.h"

#include<iostream>
#include <cstdint>
#ifdef WIN32
#include "Winsock2.h"
#endif
//#include <boost/endian/conversion.hpp>
#include <boost/make_shared.hpp>
#include <boost/format.hpp>

#include "pole.h"
#include "CCPiConsoleUserInterface.h"

namespace CCPi{

	XradiaReader::XradiaReader(std::string filename, CCPiUserApplicationInterface *uai):localUAI(false)
	{
		this->filename = filename;
		type = UINT_16;
		if(uai==NULL)
		{
			this->uai = new CCPiConsoleUserInterface();
			localUAI = true;
		}else{
			this->uai = uai;
		}
		extractMetaDataAll();
	}

	XradiaReader::~XradiaReader()
	{
		if(localUAI)
			delete uai;
	}

	int XradiaReader::getImageWidth()
	{
		return imageWidth;
	}

	int XradiaReader::getImageHeight()
	{
		return imageHeight;
	}

	XradiaReader::DataType XradiaReader::getDataType()
	{
		return type;
	}

	std::vector<float> XradiaReader::getAngles()
	{
		return angles;
	}

	int XradiaReader::getNumberOfImages()
	{
		return numberOfImages;
	}

	std::vector<float> XradiaReader::getXShifts()
	{
		return xShifts;
	}
	std::vector<float> XradiaReader::getYShifts()
	{
		return yShifts;
	}

	boost::shared_ptr< short_pixel > XradiaReader::getImageDataInShort()
	{
		return shortImageData;
	}
	boost::shared_ptr< float_pixel > XradiaReader::getImageDataInFloat()
	{
		return floatImageData;
	}

	void XradiaReader::extractMetaDataAll()
	{
		imageWidth = extractMetaDataInt("ImageInfo/ImageWidth");
		imageHeight = extractMetaDataInt("ImageInfo/ImageHeight");
		int dtype = extractMetaDataInt("ImageInfo/DataType");
		if(dtype==10)
			type = FLOAT;
		else if(dtype==5)
			type = UINT_16;
		else
			uai->LogMessage("Error data type not read properly");
		numberOfImages = extractMetaDataInt("ImageInfo/ImagesTaken");
		xShifts = extractMetaDataFloatArray("Alignment/X-Shifts", numberOfImages);
		yShifts = extractMetaDataFloatArray("Alignment/Y-Shifts", numberOfImages);
		xSamplePosition = extractMetaDataFloatArray("ImageInfo/XPosition", numberOfImages);
		ySamplePosition = extractMetaDataFloatArray("ImageInfo/YPosition", numberOfImages);
		zSamplePosition = extractMetaDataFloatArray("ImageInfo/ZPosition", numberOfImages);
		srcToObject = extractMetaDataFloatArray("ImageInfo/StoRADistance", numberOfImages)[0]; //Just use the first value
		detToObject = extractMetaDataFloatArray("ImageInfo/DtoRADistance", numberOfImages)[0]; //Just use the first value
		angles = extractMetaDataFloatArray("ImageInfo/Angles", numberOfImages);
		if(type==UINT_16){
			shortImageData = boost::make_shared< short_pixel >(boost::extents[numberOfImages][imageWidth][imageHeight]);
		}else if(type==FLOAT){
			floatImageData = boost::make_shared< float_pixel >(boost::extents[numberOfImages][imageWidth][imageHeight]);
		}
		extractImageData();
	}

	int XradiaReader::extractMetaDataInt(std::string streamName)
	{
		POLE::Storage* storage = new POLE::Storage( filename.c_str() );
		storage->open(); // Open the storage
		if(storage->result() != POLE::Storage::Ok)
		{
			uai->LogMessage("Error: opening the file: "+filename);
			return -1;
		}
		POLE::Stream* stream = new POLE::Stream(storage, streamName);
		POLE::uint64 datasize = stream->size();
		void *data = new char[datasize];
		stream->read( (unsigned char*)data, datasize);
		int* iData = (int*)data;
		int retval = iData[0];//boost::endian::little_to_native(iData[0]);
		delete[] data;
		delete stream;
		delete storage;
		return retval;
	}

	std::vector<float> XradiaReader::extractMetaDataFloatArray(std::string streamName,int numberOfValues)
	{
		std::vector<float> retval;
		POLE::Storage* storage = new POLE::Storage( filename.c_str() );
		storage->open(); // Open the storage
		if(storage->result() != POLE::Storage::Ok)
		{
			uai->LogMessage("Error: opening the file: " + filename);
			return retval;
		}
		POLE::Stream* stream = new POLE::Stream(storage, streamName);
		POLE::uint64 datasize = stream->size();
		void *data = new char[datasize];
		stream->read( (unsigned char*)data, datasize);
		float* iData = (float*)data;
		for(int i=0;i<numberOfValues;i++)
		{
			float val = iData[i];//boost::endian::little_to_native(iData[i]);
			retval.push_back(val);
		}
		delete[] data;
		delete stream;
		delete storage;
		return retval;
	}

	void XradiaReader::extractImageData()
	{
		POLE::Storage* storage = new POLE::Storage( filename.c_str() );
		storage->open(); // Open the storage
		if(storage->result() != POLE::Storage::Ok)
		{
			uai->LogMessage("Error: opening the file: "+filename);
			return;
		}
		uai->SetProgressValue(0.0);
		if(type==FLOAT)
		{
			float val;
			for(int imageId=1;imageId<=numberOfImages;imageId++)
			{
				int pageId = (int)ceil(imageId/100.0);
				std::stringstream streamName;
				streamName<<"ImageData"<<pageId+1<<"/Image"<<imageId;
				uai->SetStatusMessage("Reading Image:"+streamName.str());
				uai->SetProgressValue(imageId/(numberOfImages*1.0f));
				POLE::Stream* stream = new POLE::Stream(storage, streamName.str().c_str());
				POLE::uint64 datasize = stream->size();
				void *data = new char[datasize];
				stream->read( (unsigned char*)data, datasize);
				float* iData = (float*)data;
				for(int i=0;i<imageWidth;i++)
				{
					for(int j=0;j<imageHeight;j++)
					{
						val = iData[i*imageHeight+j];//boost::endian::little_to_native(iData[i*imageHeight+j]);
						(*floatImageData)[imageId][i][j]=val;
					}
				}
				delete[] data;
				delete stream;
			}
		}else if(type=UINT_16){
			short val;
			void *data = new char[imageHeight*imageWidth*4];
			for(int imageId=0;imageId<numberOfImages;imageId++)
			{
				int pageId = (int)ceil((imageId+1)/100.0);
				std::stringstream streamName;
				streamName<<"ImageData"<<pageId<<"/Image"<<boost::format("%i")%(imageId+1);
				uai->SetStatusMessage("Reading Image:"+streamName.str());
				uai->SetProgressValue(imageId/(numberOfImages*1.0f));
				POLE::Stream* stream = new POLE::Stream(storage, streamName.str().c_str());
				POLE::uint64 datasize = stream->size();
				stream->read( (unsigned char*)data, datasize);
				short* iData = (short*)data;
				for(int i=0;i<imageWidth;i++)
				{
					for(int j=0;j<imageHeight;j++)
					{
						val = iData[i*imageHeight+j];//boost::endian::little_to_native(iData[i*imageHeight+j]);
						(*shortImageData)[imageId][i][j]=val;
					}
				}
				delete stream;
			}
			delete[] data;
		}
		delete storage;
		uai->SetProgressValue(1.0);
	}

}
