#pragma once
class Pair
{
public:
	int x, y;
	Pair(){ x = -1; y = -1; }
	Pair(int xx, int yy){ x = xx; y = yy; }
	Pair reverse(){ return (Pair(this->y, this->x)); }
	bool operator==(const Pair& arg) const
	{
		if ((this->x == arg.x) && (this->y == arg.y))
			return (true);
		else
			return (false);
	}
	bool operator<(const Pair& arg) const
	{
		if (this->x < arg.x)
			return (true);
		
		if (this->x == arg.x)
			return( this->y < arg.y );

		return (false);
	}
	~Pair(){}
};

class PointMap
{
	//std::vector < std::vector<EPoint_2> > pMap;
	//std::vector < Pair > indexMap;
	std::vector<EPoint_2> nullVector;
	std::vector<Point_3> nullVec;

	std::vector<EPoint_2> newPoints;
	std::vector <int> newIndices;

	std::map< Pair, std::vector<EPoint_2> > pairToVectorMap;
	std::map< int, EPoint_2 > indexToPointMap;
	std::map < EPoint_2, int > invIndexToPointMap;

	std::map< int, Point_3 > indexToNewUV;
	std::map < Point_3, int > invIndexToNewUV;

	std::map< Pair, std::vector<Point_3> > EdgesUV;
	std::map< Pair, std::vector<Point_3> > edgeToMesh;

public:
	std::map< Pair, int > edgeToFace;
	PointMap(){ nullVector.clear(); nullVec.clear(); }
	~PointMap(){}
	std::vector<EPoint_2>& operator[](Pair key)
	{ 
		/*int val = -1;
		for (int i = 0; i < (int)this->indexMap.size(); ++i)
			if (indexMap[i] == key)
			{
				val = i;
				break;
			}
		if (val != -1)
			return (pMap[val]);
		else
			return (nullVector);*/
		if (pairToVectorMap.count(key) == 0)
			return (nullVector);
		return (pairToVectorMap[key]);
	}

	void insert(Pair key, std::vector<EPoint_2>& val)
	{
		pairToVectorMap[key] = val;
		//indexMap.push_back(key);
		//this->pMap.push_back(val);

		Pair revKey = key.reverse();
		std::vector<EPoint_2> revVal(val.rbegin(),val.rend());
		//std::vector<CGAL::Point_2<Kernel>> revPointsVec(pointsVec.rbegin(), pointsVec.rend());
		//for (int i = (int)val.size()-1; i >= 0; --i)
		//	revVal.push_back(val[i]);

		pairToVectorMap[revKey] = revVal;
		//indexMap.push_back(revKey);
		//this->pMap.push_back(revVal);
	}

	void insertUV(Pair key, std::vector<Point_3>& val)
	{
		EdgesUV[key] = val;
		Pair revKey = key.reverse();
		std::vector<Point_3> revVal(val.rbegin(), val.rend());
		EdgesUV[revKey] = revVal;
	}

	void insertNewMeshPoint(Pair key, std::vector<Point_3>& val)
	{
		edgeToMesh[key] = val;
		Pair revKey = key.reverse();
		std::vector<Point_3> revVal(val.rbegin(), val.rend());
		edgeToMesh[revKey] = revVal;
	}

	std::vector<Point_3>& getEdgeUVs(Pair edge)
	{
		if (EdgesUV.count(edge) == 0)
			return (nullVec);
		return (EdgesUV[edge]);
	}

	std::vector<Point_3>& getNewMeshPointsVec(Pair edge)
	{
		return (edgeToMesh[edge]);
	}

	void updateNewPoint(int newIndex, EPoint_2& newPoint)
	{
		indexToPointMap[newIndex] = newPoint;
		invIndexToPointMap[newPoint] = newIndex;
		//this->newPoints.push_back(newPoint);
		//this->newIndices.push_back(newIndex);
	}

	void updateNewUV(int newIndex, Point_3& newUV)
	{
		indexToNewUV[newIndex] = newUV;
		invIndexToNewUV[newUV] = newIndex;
	}

	int getNewIndex(EPoint_2& newPoint)
	{
		return (invIndexToPointMap[newPoint]);
		/*for (int i = 0; i < (int)newPoints.size(); ++i)
			if (newPoints[i] == newPoint)
				return (newIndices[i]);*/
	}

	EPoint_2& getNewPoint(int index)
	{
		return (indexToPointMap[index]);
		//for (int i = 0; i < (int)newIndices.size(); ++i)
		//	if (newIndices[i] == index)
		//		return (newPoints[i]);
	}

	int numOfNewPoints()
	{
		return((int)indexToPointMap.size());
	}

	void test()
	{
		/*for (int i = 0; i < (int)pMap.size(); ++i)
		{
			std::cout << "Pair - (" << indexMap[i].x << "," << indexMap[i].y << ")\n";
			for (int j = 0; j < (int)pMap[i].size(); ++j)
			{
				std::cout << "(" << pMap[i][j].x() << "," << pMap[i][j].y() << ")\n";
			}
		}*/
	
	}
};


void loadSourceMesh( Mesh &source_mesh , std::vector<Kernel::Point_3> &pVec , std::vector<int> &fVec );
void addPointsToTarget( Polygon_2 &poly , int numOfBorder , double avg_arc );
void HarmonicFlattening(Mesh &source_mesh, GMMSparseRowMatrix &u, GMMSparseRowMatrix &weightsMat, bool harmonic = true);
void getPointsFromFace( const Arrangement_2::Face_const_handle& face, std::vector<EPoint_2>& points , std::vector<int>& indicesOrder);
void getPointsFromFace_Mesh( Mesh& targetMesh/*const Mesh::Face_const_handle& face*/, std::vector<EPoint_2>& points , std::vector<int>& indicesOrder );
void BuildArrangement(Arrangement_2& arr, Landmarks_pl& trap, const std::vector<EPoint_2>& vertices, const std::vector<int>& faces, /*Face_index_observer& obs,*/ Mesh &source_mesh);
int findTarget(const Landmarks_pl& target, const EPoint_2& point, int& type, Arrangement_2::Face_const_handle& targetFace);
ARRNumberType crossProduct ( EVector_2 v1 , EVector_2 v2 );
void barycentricCord( const std::vector<EPoint_2>& points, EPoint_2 point, ARRTraits_2::Point_3 &res );
EPoint_2 reverseBarycentric ( const std::vector<EPoint_2>& points, ARRTraits_2::Point_3 bar );
void updateMeshUV( Mesh& mesh, int index , Mesh::Point_3 point );
std::vector<int> updateUVs(Mesh& sourceMesh, Mesh& targetMesh, const Arrangement_2& source, const Arrangement_2& target, const Landmarks_pl& sourceLandMark, const Landmarks_pl& targetLandMark, /*Face_index_observer& sourceObs, Face_index_observer& targetObs,*/ std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<int>& fVec, std::vector<Point_3>& uvVector);
void setBoundaryUV( Mesh &source_mesh, Mesh &target_mesh, std::vector<Point_3>& uvVector );

void matchPointsIndices(Arrangement_2& arrSource, const std::vector<EPoint_2>& sourceHarmonicMapPoints, const Landmarks_pl& sourceLandMark, Mesh& sourceMesh);
void matchEdges(Arrangement_2& arr, Mesh& source_mesh);
void matchFaces(Arrangement_2& arr, const std::vector<EPoint_2>& mapPoints, const std::vector<int> &fVec, Landmarks_pl& trap, Mesh& sourceMesh);

int refine(std::vector<int>& neg, Mesh& sourceMesh, Mesh& targetMesh, Arrangement_2& source, Arrangement_2& target, const Landmarks_pl& sourceLandMark, const Landmarks_pl& targetLandMark, /*Face_index_observer& targetObs,*/ std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<Kernel::Point_3> &pVec, std::vector<int>& fVec, std::vector<Point_3>& uvVector);
void findIntersection(Arrangement_2& target, const Landmarks_pl& targetLandMark, ESegment_2 seg, std::vector<EPoint_2>& intersectionPoints);
void inverseList(std::vector<EPoint_2>& intersectionPoints);
Point_3 calcNewUV(EPoint_2 tempP, Mesh& targetMesh, const Landmarks_pl& targetLandMark, /*Face_index_observer& targetObs,*/ std::vector<Point_3>& uvVector);
Point_3 calcNewMeshPoint(EPoint_2 tempP, int faceIndex, std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<Kernel::Point_3> &pVec, std::vector<int>& fVec);
void extractIndicesFromPair(Pair& p, PointMap& pMap, std::vector<int>& polygonIndices, std::map<Point_3, int>& newUVtoIndicesMap);
void extractIndicesFromPair(Pair& e1, Pair& e2, Pair& e3, PointMap& pMap, std::vector<int>& polygonIndices, std::map<Point_3, int>& newUVtoIndicesMap);

void extractIndicesFromPair(Pair& p, PointMap& pMap, std::vector<int>& polygonIndices);
void extractIndicesFromPair(Pair& e1, Pair& e2, Pair& e3, PointMap& pMap, std::vector<int>& polygonIndices);

void triangulateNewPolygon(std::vector<int>& polygonIndices, std::vector<Point_3>& uvVector, std::vector<int>& fVec);

void findNegativeNeighbors(std::vector<int>& neg, std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<int>& fVec, Arrangement_2& source, const Landmarks_pl& sourceLandMark);
void triangulateNeighbors(std::vector<int>& neighTri, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, PointMap& pMap, std::vector<Point_3>& uvVector);
bool triangulateNewPolygon(std::vector<int>& polygonIndices, std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<int>& fVec, PointMap& pMap, std::vector<Point_3>& uvVector);

bool checkIfBoundaryEdge(ESegment_2& seg, const Landmarks_pl& sourceLandMark);

void refineTriangle(int triIndex, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, PointMap& pMap, std::vector<Point_3>& uvVector, Mesh& targetMesh, Arrangement_2& target, const Landmarks_pl& targetLandMark, const Landmarks_pl& sourceLandMark, /*Face_index_observer& targetObs,*/ std::vector<Kernel::Point_3> &pVec, int& updatedNumOfPoints);
void findTriangleNeighbors(int triIndex, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, const Landmarks_pl& sourceLandMark, std::vector<int>& neighTri, std::vector<bool>& inTheList);
bool triangulateNeighbor(int triIndex, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, PointMap& pMap, std::vector<Point_3>& uvVector, std::map<Point_3, int>& newUVtoIndicesMap);

bool checkIfSimple(int triIndex, std::vector<int>& fVec, PointMap& pMap, std::vector<Point_3>& uvVector);

void simplify(PointMap& pMap, std::vector<int>& fVec, std::vector<Point_3>& uvVector);

void barycentricCord(const std::vector<Point_3>& points, Point_3 point, Point_3 &res);
double crossProduct(Vector_3 v1, Vector_3 v2);

void meanValueWeights(Mesh &source_mesh, GMMSparseRowMatrix &u, GMMSparseRowMatrix &weightsMat);

int countFoldovers(const Mesh& mesh, bool printReport = false);