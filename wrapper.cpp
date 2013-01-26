#include <tuple>
#include <map>

#include "PolyVoxCore/MaterialDensityPair.h"
#include "PolyVoxCore/CubicSurfaceExtractorWithNormals.h"
#include "PolyVoxCore/CubicSurfaceExtractor.h"
#include "PolyVoxCore/SurfaceMesh.h"
#include "PolyVoxCore/SimpleVolume.h"
#include "PolyVoxCore/AmbientOcclusionCalculator.h"
#include "PolyVoxCore/AStarPathfinder.h"
#include "PolyVoxCore/MeshDecimator.h"
#include "PolyVoxCore/Raycast.h"

extern "C"
{
#include "wrapper.h"
}

using namespace PolyVox;

std::map< PVoxAOCalculator,std::tuple< AmbientOcclusionCalculator< VOLUMETYPE >*, Array3DUint8*, uint8_t(*)[3] > > _AOCalculator;
uint32_t _AOCalculatorSeq = 0;

std::map< PVoxVolume,VOLUMETYPE* > _Volume;
uint32_t _VolumeSeq;

std::map< PVoxRegion,Region* > _Region;
uint32_t _RegionSeq;

std::map< PVoxAStarPathfinder,std::tuple< AStarPathfinder< VOLUMETYPE >*,std::list< Vector3DInt32 >*,PVoxVectorNDUint32* > > _ASPathfinder;
uint32_t _ASPathfinderSeq;

std::map< PVoxSurfaceExtractor,EXTRACTORTYPE* > _SurfaceExtractor;
uint32_t _SurfaceExtractorSeq;

std::map< PVoxSurfaceMesh,SurfaceMesh< VERTEXTYPE >* > _SurfaceMesh;
uint32_t _SurfaceMeshSeq;

std::map< PVoxMeshDecimator,MeshDecimator< VERTEXTYPE >* > _MeshDecimator;
uint32_t _MeshDecimatorSeq;

std::map< PVoxRaycast,std::tuple< Raycast< VOLUMETYPE >*,PVoxRaycastResult*,RaycastResult* > > _Raycast;
uint32_t _RaycastSeq;

PVoxVoxel convertVoxelCppToC(const VOXELCLASS& voxel)
{ //TODO
	PVoxVoxel ret;
	ret.material = voxel.getMaterial();
	return ret;
}

const VOXELCLASS convertVoxelCToCpp(PVoxVoxel& voxel)
{ //TODO
	return VOXELCLASS(voxel.material);
}

extern "C"
{
	PVoxAOCalculator pvoxAOCalculatorAdd
	(
		PVoxVolume volume,
		uint8_t (*result)[3],
		PVoxRegion region,
		float fRayLength,
		uint8_t uNoOfSamplesPerOutputElement,
		bool (*funcIsTransparent)(PVoxVoxel voxel)
	)
	{
		Array3DUint8 *arr = new Array3DUint8();
		_AOCalculator[++_AOCalculatorSeq] = std::make_tuple< AmbientOcclusionCalculator< VOLUMETYPE >*, Array3DUint8*, uint8_t(*)[3] >
			(
			new AmbientOcclusionCalculator< VOLUMETYPE >(
					_Volume[volume], arr, *_Region[region], fRayLength,
					uNoOfSamplesPerOutputElement, [&](const VOXELCLASS& voxel) { return funcIsTransparent(convertVoxelCppToC(voxel)); } ),
				arr,
				result
			);
		return _AOCalculatorSeq;
	}

	void pvoxAOCalculatorRemove(PVoxAOCalculator handle)
	{
		std::tuple< AmbientOcclusionCalculator< VOLUMETYPE >*, Array3DUint8*, uint8_t(*)[3] > tuple = _AOCalculator[handle];
		delete std::get<0>(tuple);
		delete std::get<1>(tuple);
	}

	void pvoxAOCalculatorExecute(PVoxAOCalculator handle)
	{
		std::tuple< AmbientOcclusionCalculator< VOLUMETYPE >*, Array3DUint8*, uint8_t(*)[3] > tuple = _AOCalculator[handle];
		std::get<0>(tuple)->execute();

		uint8_t* rd = std::get<1>(tuple)->getRawData();
		(*std::get<2>(tuple))[0] = rd[0];
		(*std::get<2>(tuple))[1] = rd[1];
		(*std::get<2>(tuple))[2] = rd[2];
	}

	PVoxAStarPathfinder pvoxASPathfinderAdd(
		PVoxVolume volData,
		PVoxVector3DUint32* v3dStart,
		PVoxVector3DUint32* v3dEnd,
		PVoxVectorNDUint32* listResult,
		float fHBias/*=1.0*/,
		uint32_t uMaxNoOfNodes/*=10000*/,
		//Connectivity connectivity=TwentySixConnected,
		bool funcIsVoxelValidForPath(const PVoxVolume, const PVoxVector3DInt32*)/*=NULL*/,
		void funcProgressCallback(float)/*=NULL*/
	)
	{
		std::list< Vector3DInt32 >* list = new std::list< Vector3DInt32 >();

		polyvox_function<bool (const VOLUMETYPE*, const Vector3DInt32&)> _isVoxelValidForPath =
			[&](const VOLUMETYPE* volume, const Vector3DInt32& vec_cpp) -> bool
			{	//TODO: birectional map or something
				PVoxVector3DInt32 vec_c;
				vec_c.x = vec_cpp.getX();
				vec_c.y = vec_cpp.getY();
				vec_c.z = vec_cpp.getZ();
				
				for (std::map< PVoxVolume,VOLUMETYPE* >::iterator it = _Volume.begin(); it != _Volume.end(); it++)
				{
					if ( it->second == volume )
						return funcIsVoxelValidForPath(it->first, &vec_c);
				}
				assert(0);
			};

		_ASPathfinder[++_ASPathfinderSeq] = std::make_tuple< AStarPathfinder< VOLUMETYPE >*,std::list< Vector3DInt32 >*,PVoxVectorNDUint32* >(
			new AStarPathfinder< VOLUMETYPE >(
				AStarPathfinderParams< VOLUMETYPE >(
					_Volume[volData],
					Vector3DInt32(v3dStart->x, v3dStart->y, v3dStart->z),
					Vector3DInt32(v3dEnd->x, v3dEnd->y, v3dEnd->z),
					list, fHBias, uMaxNoOfNodes, TwentySixConnected,
					funcIsVoxelValidForPath == NULL ? NULL : _isVoxelValidForPath,
					funcProgressCallback
				)
			),
			list,
			listResult
		);
		return _ASPathfinderSeq;
	}

	void pvoxASPathfinderRemove(PVoxAStarPathfinder handle)
	{
		std::tuple< AStarPathfinder< VOLUMETYPE >*,std::list< Vector3DInt32 >*,PVoxVectorNDUint32* > tuple = _ASPathfinder[handle];
		delete std::get<0>(tuple);
		delete std::get<1>(tuple);
		delete[] std::get<2>(tuple)->nodes;
	}

	void pvoxASPathfinderExecute(PVoxAStarPathfinder handle)
	{
		std::tuple< AStarPathfinder< VOLUMETYPE >*,std::list< Vector3DInt32 >*,PVoxVectorNDUint32* > tuple = _ASPathfinder[handle];
		std::get<0>(tuple)->execute();

		std::list< Vector3DInt32 > *list = std::get<1>(tuple);
		PVoxVectorNDUint32* path = std::get<2>(tuple);
		path->count = list->size();
		path->nodes = new PVoxVector3DUint32[path->count];
		int i = 0;
		for (std::list< Vector3DInt32 >::iterator it = list->begin(); it != list->end(); it++, i++)
		{
			PVoxVector3DUint32* node = path->nodes+i;
			node->x = it->getX();
			node->y = it->getY();
			node->z = it->getZ();
		}
	}

	//
	
	PVoxVolume pvoxVolumeAdd(PVoxRegion region)
	{
		_Volume[++_VolumeSeq] = new VOLUMETYPE(*_Region[region]);
		return _VolumeSeq;
	}

	PVoxVoxel pvoxVolumeGetBorderValue(PVoxVolume volume)
	{
		return convertVoxelCppToC(_Volume[volume]->getBorderValue());
	}

	int32_t pvoxVolumeGetWidth(PVoxVolume volume)
	{
		return _Volume[volume]->getWidth();
	}

	int32_t pvoxVolumeGetHeight(PVoxVolume volume)
	{
		return _Volume[volume]->getHeight();
	}

	int32_t pvoxVolumeGetLongestSideLength(PVoxVolume volume)
	{
		return _Volume[volume]->getLongestSideLength();
	}

	int32_t pvoxVolumeGetShortestSideLength(PVoxVolume volume)
	{
		return _Volume[volume]->getShortestSideLength();
	}

	float pvoxVolumeGetDiagonalLength(PVoxVolume volume)
	{
		return _Volume[volume]->getDiagonalLength();
	}

	PVoxVoxel pvoxVolumeGetVoxelAt(PVoxVolume volume, PVoxVector3DUint32 vec)
	{
		return convertVoxelCppToC(_Volume[volume]->getVoxelAt(Vector3DInt32(vec.x, vec.y, vec.z)));
	}

	void pvoxVolumeSetBorderValue(PVoxVolume volume, PVoxVoxel tBorder)
	{
		_Volume[volume]->setBorderValue(convertVoxelCToCpp(tBorder));
	}

	bool pvoxVolumeSetVoxelAt(PVoxVolume volume, PVoxVector3DUint32 vec, PVoxVoxel tValue)
	{
		return _Volume[volume]->setVoxelAt(Vector3DInt32(vec.x, vec.y, vec.z), convertVoxelCToCpp(tValue));
	}

	uint32_t pvoxVolumeCalculateSizeInBytes(PVoxVolume volume)
	{
		return _Volume[volume]->calculateSizeInBytes();
	} 

	void pvoxVolumeRemove(PVoxVolume volume)
	{
		delete _Volume[volume];
	}

//#if VOLUMETYPE == VOLUMETYPE_LARGE
//	//
//#endif

	//

	PVoxSurfaceExtractor pvoxSurfaceExtractorAdd
	(
		PVoxVolume volData,
		PVoxRegion region,
		PVoxSurfaceMesh result,
		bool bMergeQuads=true
		//isQuadNeeded
	)
	{
		_SurfaceExtractor[++_SurfaceExtractorSeq] = new EXTRACTORTYPE(_Volume[volData], *_Region[region], _SurfaceMesh[result], bMergeQuads);
		return _SurfaceExtractorSeq;
	}

	void pvoxSurfaceExtractorExecute(PVoxSurfaceExtractor handle)
	{
		_SurfaceExtractor[handle]->execute();
	}

	void pvoxSurfaceExtractorRemove(PVoxSurfaceExtractor handle)
	{
		delete _SurfaceExtractor[handle];
	}

	//

	PVoxSurfaceMesh pvoxSurfaceMeshAdd()
	{
		_SurfaceMesh[++_SurfaceMeshSeq] = new SurfaceMesh< VERTEXTYPE >();
		return _SurfaceMeshSeq;
	}

	void pvoxSurfaceMeshRemove(PVoxSurfaceMesh handle)
	{
		delete _SurfaceMesh[handle];
	}

	uint32_t pvoxSurfaceMeshGetNoOfIndices(PVoxSurfaceMesh handle)
	{
		return _SurfaceMesh[handle]->getNoOfIndices();
	}

	void pvoxSurfaceMeshGetIndices(PVoxSurfaceMesh handle, uint32_t result[])
	{
		std::vector< uint32_t > indices = _SurfaceMesh[handle]->getIndices();
		size_t i = 0;
		for( std::vector< uint32_t >::const_iterator it = indices.begin(); it != indices.end(); it++)
		{
			result[i++] = *it;
		}
	}

	uint32_t pvoxSurfaceMeshGetNoOfVertices(PVoxSurfaceMesh handle)
	{
		return _SurfaceMesh[handle]->getNoOfVertices();
	}

	void pvoxSurfaceMeshGetVertices(PVoxSurfaceMesh handle, PVoxVertex result[])
	{
		std::vector< VERTEXTYPE > vertices = _SurfaceMesh[handle]->getVertices();
		size_t i = 0;
		for( std::vector< VERTEXTYPE >::const_iterator it = vertices.begin(); it != vertices.end(); it++, i++)
		{
			result[i].material = it->getMaterial();
			Vector3DFloat pos = it->getPosition();
			result[i].position.x = pos.getX();
			result[i].position.y = pos.getY();
			result[i].position.z = pos.getZ();
//#if VERTEXTYPE == VERTEXTYPE_POSMATNOR
//			result[i].normal = it->getNormal();
//#endif
		}
	}

	uint32_t pvoxSurfaceMeshGetNoOfNonUniformTrianges(PVoxSurfaceMesh handle)
	{
		return _SurfaceMesh[handle]->getNoOfNonUniformTrianges();
	}

	uint32_t pvoxSurfaceMeshGetNoOfUniformTrianges(PVoxSurfaceMesh handle)
	{
		return _SurfaceMesh[handle]->getNoOfUniformTrianges();
	}

	//TODO: getRawVertexData (void)

	//

	PVoxMeshDecimator pvoxMeshDecimatorAdd(PVoxSurfaceMesh input, PVoxSurfaceMesh output, float fEdgeCollapseThreshold=0.95f)
	{
		_MeshDecimator[++_MeshDecimatorSeq] = new MeshDecimator< VERTEXTYPE >(_SurfaceMesh[input], _SurfaceMesh[output], fEdgeCollapseThreshold);
		return _MeshDecimatorSeq;
	}

	void pvoxMeshDecimatorExecute(PVoxMeshDecimator handle)
	{
		_MeshDecimator[handle]->execute();
	}

	void pvoxMeshDecimatorRemove(PVoxMeshDecimator handle)
	{
		delete _MeshDecimator[handle];
	}

	//

	bool funcIsPassable(const VOLUMETYPE::Sampler& sampler)
	{
		return true;
	}

	PVoxRaycast pvoxRaycastAdd(PVoxVolume volume, PVoxVector3DFloat* start, PVoxVector3DFloat* direction, PVoxRaycastResult* result)
	{
		RaycastResult *res = new RaycastResult;
		_Raycast[++_RaycastSeq] = std::make_tuple< Raycast< VOLUMETYPE >*,PVoxRaycastResult*,RaycastResult* >(
			new Raycast< VOLUMETYPE >(
				_Volume[volume],
				Vector3DFloat(start->x,start->y,start->z),
				Vector3DFloat(direction->x,direction->y,direction->z),
				*res,
				funcIsPassable
			),
			result,
			res);
		return _RaycastSeq;
	}

	void pvoxRaycastExecute(PVoxRaycast handle)
	{
		std::tuple< Raycast< VOLUMETYPE >*,PVoxRaycastResult*,RaycastResult* > tuple = _Raycast[handle];
		std::get<0>(tuple)->execute();
		PVoxRaycastResult* out = std::get<1>(tuple);
		RaycastResult* in = std::get<2>(tuple);
		out->intersectionFound = in->foundIntersection;
		out->intersectionVoxel.x = in->intersectionVoxel.getX();
		out->intersectionVoxel.y = in->intersectionVoxel.getY();
		out->intersectionVoxel.z = in->intersectionVoxel.getZ();
		out->previousVoxel.x = in->previousVoxel.getX();
		out->previousVoxel.y = in->previousVoxel.getY();
		out->previousVoxel.z = in->previousVoxel.getZ();
	}

	void pvoxRaycastRemove(PVoxRaycast handle)
	{
		std::tuple< Raycast< VOLUMETYPE >*,PVoxRaycastResult*,RaycastResult* > tuple = _Raycast[handle];
		delete std::get<0>(tuple);
		delete std::get<2>(tuple);
	}

	//TODO: raycast with callback

	PVoxRegion pvoxRegionAdd()
	{
		_Region[++_RegionSeq] = new Region;
		return _RegionSeq;
	}

	PVoxVector3DInt32 pvoxRegionGetLowerCorner(PVoxRegion handle)
	{
		Vector3DInt32 corner = _Region[handle]->getLowerCorner();
		PVoxVector3DInt32 ret;
		ret.x = corner.getX();
		ret.y = corner.getY();
		ret.z = corner.getZ();
		return ret;
	}

	PVoxVector3DInt32 pvoxRegionGetUpperCorner(PVoxRegion handle)
	{
		Vector3DInt32 corner = _Region[handle]->getUpperCorner();
		PVoxVector3DInt32 ret;
		ret.x = corner.getX();
		ret.y = corner.getY();
		ret.z = corner.getZ();
		return ret;
	}

	void pvoxRegionSetUpperCorner(PVoxRegion handle, PVoxVector3DUint32* upperCorner)
	{
		_Region[handle]->setUpperCorner(Vector3DInt32(upperCorner->x,upperCorner->y,upperCorner->z));
	}

	void pvoxRegionSetLowerCorner(PVoxRegion handle, PVoxVector3DUint32* lowerCorner)
	{
		_Region[handle]->setLowerCorner(Vector3DInt32(lowerCorner->x,lowerCorner->y,lowerCorner->z));
	}

	bool pvoxRegionContainsPointF(PVoxRegion handle, PVoxVector3DFloat* pos, float boundary)
	{
		return _Region[handle]->containsPoint(Vector3DFloat(pos->x, pos->y, pos->z), boundary);
	}

	bool pvoxRegionContainsPointI(PVoxRegion handle, PVoxVector3DInt32* pos, uint8_t boundary)
	{
		return _Region[handle]->containsPoint(Vector3DInt32(pos->x, pos->y, pos->z), boundary);
	}

	bool pvoxRegionContainsPointInXF(PVoxRegion handle, float x, float boundary)
	{
		return _Region[handle]->containsPointInX(x, boundary);
	}

	bool pvoxRegionContainsPointInXI(PVoxRegion handle, uint32_t x, uint8_t boundary)
	{
		return _Region[handle]->containsPointInX(x, boundary);
	}

	bool pvoxRegionContainsPointInYF(PVoxRegion handle, float y, float boundary)
	{
		return _Region[handle]->containsPointInY(y, boundary);
	}

	bool pvoxRegionContainsPointInYI(PVoxRegion handle, uint32_t y, uint8_t boundary)
	{
		return _Region[handle]->containsPointInY(y, boundary);
	}

	bool pvoxRegionContainsPointInZF(PVoxRegion handle, float z, float boundary)
	{
		return _Region[handle]->containsPointInZ(z, boundary);
	}

	bool pvoxRegionContainsPointInZI(PVoxRegion handle, uint32_t z, uint8_t boundary)
	{
		return _Region[handle]->containsPointInZ(z, boundary);
	}

	void pvoxRegionCropTo(PVoxRegion handle, PVoxRegion other)
	{
		_Region[handle]->cropTo(*_Region[other]);
	}

	void pvoxRegionShift(PVoxRegion handle, PVoxVector3DInt32* vec)
	{
		_Region[handle]->shift(Vector3DInt32(vec->x, vec->y, vec->z));
	}

	void pvoxRegionShiftLowerCorner(PVoxRegion handle, PVoxVector3DInt32* vec)
	{
		_Region[handle]->shiftLowerCorner(Vector3DInt32(vec->x, vec->y, vec->z));
	}

	void pvoxRegionShiftUpperCorner(PVoxRegion handle, PVoxVector3DInt32* vec)
	{
		_Region[handle]->shiftUpperCorner(Vector3DInt32(vec->x, vec->y, vec->z));
	}

	void pvoxRegionRemove(PVoxRegion handle)
	{
		delete _Region[handle];
	}

	//


}