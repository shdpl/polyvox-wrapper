#include <stdint.h>


#define HANDLE int
#define NULL 0
#ifndef __cplusplus
#define bool char
#endif

typedef HANDLE PVoxRegion;
typedef HANDLE PVoxAOCalculator;
typedef HANDLE PVoxVolume;
typedef HANDLE PVoxAStarPathfinder;
typedef HANDLE PVoxSurfaceExtractor;
typedef HANDLE PVoxSurfaceMesh;
typedef HANDLE PVoxArray;
typedef HANDLE PVoxMeshDecimator;
typedef HANDLE PVoxRaycast;


//
#define MATERIALTYPE					uint16_t
//#define DENSITYTYPE					uint8_t
#define VOXELTYPE						VOXELTYPE_MATERIAL

#ifdef MATERIALTYPE
	#ifdef DENSITYTYPE
		#define VOXELCLASS				MaterialDensityPair< MATERIALTYPE, DENSITYTYPE >
	#else
		#define VOXELCLASS				Material< MATERIALTYPE >
	#endif
#else
	#define VOXELCLASS					Density< DENSITYTYPE >
#endif
	

//
#define VOLUMETYPE						VOLUMETYPE_SIMPLE

#define VOLUMETYPE_SIMPLE				SimpleVolume< VOXELCLASS >
#define VOLUMETYPE_LARGE				LargeVolume< VOXELCLASS >
#define VOLUMETYPE_RAW					RawVolume< VOXELCLASS >


//
#define VERTEXTYPE						VERTEXTYPE_POSMAT

#define VERTEXTYPE_POSMATNOR			PositionMaterialNormal
#define VERTEXTYPE_POSMAT				PositionMaterial


//
#define EXTRACTORTYPE					EXTRACTORTYPE_CUBIC

#if VERTEXTYPE == VERTEXTYPE_POSMAT
	#define EXTRACTORTYPE_CUBIC			CubicSurfaceExtractor< VOLUMETYPE >
#else
	#define EXTRACTORTYPE_CUBIC			CubicSurfaceExtractorWithNormals< VOLUMETYPE >
#endif
#define EXTRACTORTYPE_MARCHINGCUBES		MarchingCubesSurfaceExtractor< VOLUMETYPE >


struct PVoxVector3DInt32
{
	int32_t x;
	int32_t y;
	int32_t z;
};
typedef struct PVoxVector3DInt32 PVoxVector3DInt32;

struct PVoxVector3DUint32
{
	uint32_t x;
	uint32_t y;
	uint32_t z;
};
typedef struct PVoxVector3DUint32 PVoxVector3DUint32;

struct PVoxVector3DFloat
{
	float x;
	float y;
	float z;
};
typedef struct PVoxVector3DFloat PVoxVector3DFloat;

struct PVoxVectorNDUint32
{
	uint32_t count;
	PVoxVector3DUint32* nodes; 
};
typedef struct PVoxVectorNDUint32 PVoxVectorNDUint32;

struct PVoxRaycastResult
{
	bool intersectionFound;
	PVoxVector3DUint32 intersectionVoxel;
	PVoxVector3DUint32 previousVoxel;
};
typedef struct PVoxRaycastResult PVoxRaycastResult;


struct PVoxVoxel
{
	uint32_t	material;
};
typedef struct PVoxVoxel PVoxVoxel;


struct PVoxVertex
{
	float material;
	PVoxVector3DFloat position;
/*#if VERTEXCLASS == VERTEXCLASS_POSMATNOR
	PVoxVector3DFloat normal;
#endif*/
};
typedef struct PVoxVertex PVoxVertex;

#define APIEXPORT __declspec(dllexport)

APIEXPORT PVoxAOCalculator pvoxAOCalculatorAdd
(
	PVoxVolume volume,
	uint8_t (*result)[3],
	PVoxRegion region,
	float fRayLength,
	uint8_t uNoOfSamplesPerOutputElement,
	bool (*funcIsTransparent)(PVoxVoxel voxel)
);
APIEXPORT void pvoxAOCalculatorRemove(PVoxAOCalculator handle);
APIEXPORT void pvoxAOCalculatorExecute(PVoxAOCalculator handle);


APIEXPORT PVoxAStarPathfinder pvoxASPathfinderAdd(
	PVoxVolume volData,
	PVoxVector3DUint32* v3dStart,
	PVoxVector3DUint32* v3dEnd,
	PVoxVectorNDUint32* listResult,
	float fHBias/*=1.0*/,
	uint32_t uMaxNoOfNodes/*=10000*/,
	//Connectivity connectivity=TwentySixConnected,
	bool funcIsVoxelValidForPath(const PVoxVolume, const PVoxVector3DInt32*)/*=NULL*/,
	void funcProgressCallback(float)/*=NULL*/
);
APIEXPORT void pvoxASPathfinderRemove(PVoxAStarPathfinder handle);
APIEXPORT void pvoxASPathfinderExecute(PVoxAStarPathfinder handle);


APIEXPORT PVoxVolume pvoxVolumeAdd(PVoxRegion region);
//
APIEXPORT PVoxVoxel pvoxVolumeGetBorderValue(PVoxVolume volume);
APIEXPORT void pvoxVolumeSetBorderValue(PVoxVolume volume, PVoxVoxel tBorder);
//
APIEXPORT int32_t pvoxVolumeGetWidth(PVoxVolume volume);
//
APIEXPORT int32_t pvoxVolumeGetHeight(PVoxVolume volume);
//
APIEXPORT int32_t pvoxVolumeGetLongestSideLength(PVoxVolume volume);
APIEXPORT int32_t pvoxVolumeGetShortestSideLength(PVoxVolume volume);
//
APIEXPORT float pvoxVolumeGetDiagonalLength(PVoxVolume volume);
//
APIEXPORT PVoxVoxel pvoxVolumeGetVoxelAt(PVoxVolume volume, PVoxVector3DUint32 vec);
APIEXPORT bool pvoxVolumeSetVoxelAt(PVoxVolume volume, PVoxVector3DUint32 vec, PVoxVoxel tValue);
//
APIEXPORT uint32_t pvoxVolumeCalculateSizeInBytes(PVoxVolume volume);
//
APIEXPORT void pvoxVolumeRemove(PVoxVolume volume);


APIEXPORT PVoxSurfaceExtractor pvoxSurfaceExtractorAdd
(
	PVoxVolume volData,
	PVoxRegion region,
	PVoxSurfaceMesh result,
	bool bMergeQuads/*=true*/
	//isQuadNeeded
);
APIEXPORT void pvoxSurfaceExtractorRemove(PVoxSurfaceExtractor handle);
//
APIEXPORT void pvoxSurfaceExtractorExecute(PVoxSurfaceExtractor handle);

APIEXPORT PVoxSurfaceMesh pvoxSurfaceMeshAdd();
APIEXPORT void pvoxSurfaceMeshRemove(PVoxSurfaceMesh handle);
//
APIEXPORT uint32_t pvoxSurfaceMeshGetNoOfIndices(PVoxSurfaceMesh handle);
APIEXPORT void pvoxSurfaceMeshGetIndices(PVoxSurfaceMesh handle, uint32_t result[]);
//
APIEXPORT uint32_t pvoxSurfaceMeshGetNoOfVertices(PVoxSurfaceMesh handle);
APIEXPORT void pvoxSurfaceMeshGetVertices(PVoxSurfaceMesh handle, PVoxVertex result[]);
//
APIEXPORT uint32_t pvoxSurfaceMeshGetNoOfNonUniformTrianges(PVoxSurfaceMesh handle);
APIEXPORT uint32_t pvoxSurfaceMeshGetNoOfUniformTrianges(PVoxSurfaceMesh handle);
//TODO: getRawVertexData (void)


APIEXPORT PVoxMeshDecimator pvoxMeshDecimatorAdd(PVoxSurfaceMesh input, PVoxSurfaceMesh output, float fEdgeCollapseThreshold/*=0.95f*/);
APIEXPORT void pvoxMeshDecimatorRemove(PVoxMeshDecimator handle);
//
APIEXPORT void pvoxMeshDecimatorExecute(PVoxMeshDecimator handle);


APIEXPORT PVoxRaycast pvoxRaycastAdd(PVoxVolume volume, PVoxVector3DFloat* start, PVoxVector3DFloat* direction, PVoxRaycastResult* result);
APIEXPORT void pvoxRaycastRemove(PVoxRaycast handle);
//
APIEXPORT void pvoxRaycastExecute(PVoxRaycast handle);


APIEXPORT PVoxRegion pvoxRegionAdd();
APIEXPORT void pvoxRegionRemove(PVoxRegion handle);
//
APIEXPORT PVoxVector3DInt32 pvoxRegionGetLowerCorner(PVoxRegion handle);
APIEXPORT PVoxVector3DInt32 pvoxRegionGetUpperCorner(PVoxRegion handle);
//
APIEXPORT void pvoxRegionSetUpperCorner(PVoxRegion handle, PVoxVector3DUint32* upperCorner);
APIEXPORT void pvoxRegionSetLowerCorner(PVoxRegion handle, PVoxVector3DUint32* lowerCorner);
//
APIEXPORT bool pvoxRegionContainsPointF(PVoxRegion handle, PVoxVector3DFloat* pos, float boundary/*=0.0f*/);
APIEXPORT bool pvoxRegionContainsPointI(PVoxRegion handle, PVoxVector3DInt32* pos, uint8_t boundary/*=0.0f*/);
//
APIEXPORT bool pvoxRegionContainsPointInXF(PVoxRegion handle, float x, float boundary/*=0.0f*/);
APIEXPORT bool pvoxRegionContainsPointInXI(PVoxRegion handle, uint32_t x, uint8_t boundary/*=0.0f*/);
//
APIEXPORT bool pvoxRegionContainsPointInYF(PVoxRegion handle, float y, float boundary/*=0.0f*/);
APIEXPORT bool pvoxRegionContainsPointInYI(PVoxRegion handle, uint32_t y, uint8_t boundary/*=0.0f*/);
//
APIEXPORT bool pvoxRegionContainsPointInZF(PVoxRegion handle, float z, float boundary/*=0.0f*/);
APIEXPORT bool pvoxRegionContainsPointInZI(PVoxRegion handle, uint32_t z, uint8_t boundary/*=0.0f*/);
//
APIEXPORT void pvoxRegionCropTo(PVoxRegion handle, PVoxRegion other);
APIEXPORT void pvoxRegionShift(PVoxRegion handle, PVoxVector3DInt32* vec);
//
APIEXPORT void pvoxRegionShiftLowerCorner(PVoxRegion handle, PVoxVector3DInt32* vec);
APIEXPORT void pvoxRegionShiftUpperCorner(PVoxRegion handle, PVoxVector3DInt32* vec);