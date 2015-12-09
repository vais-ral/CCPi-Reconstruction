#include <base_types.hpp>
#include <instruments.hpp>
#include <algorithms.hpp>
#include <results.hpp>

namespace CCPi{	

	class ReconstructionImpl
	{
	public:
		ReconstructionImpl(void);
		~ReconstructionImpl(void);
		void setDeviceId(devices id);
		void setAlgorithmId(algorithms id);
		void setNumberOfProcessors(int nprocs);
		bool run();
		bool saveResults(std::string outputFilename, CCPi::output_format outputFormat);
		void enableHyperThreads();
		void disableHyperThreads();
		void setResolution(int resolution);
		void setNumberOfIterations(int iterations);
		void setRegularisation(double value);
		void enableBeamHardening();
		void disableBeamHardening();
		void setBeamHardening(bool value);
		void setFilename(std::string name);
	private:
		devices deviceId;
		int numberOfProcessors;
		algorithms algorithmId;
		bool bHyperThreads;
		int resolution;
		int numberOfIterations;
		double regularise;
		bool bBeamHardening;
		voxel_data *voxels;
		real voxel_origin[3];
		real voxel_size[3];
		std::string filename;
	};

}