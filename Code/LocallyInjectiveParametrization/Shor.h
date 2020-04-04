#pragma once
class Shor
{
private:
	int N ;
	int** Q , **K ;
	int *rArray;
	// ------------for ear clipping----------------
	std::vector<int> rVector;
	std::vector<unsigned int> triangles;
	Polygon_2 simpPoly;
	std::map<int, int> simpToOriginalIndices;
	Polygon_2 tempForSwap;
	// --------------------------------------------
	Angle** firstAng , **lastAng ;
	Polygon_2 poly;
	Polygon_2 bPoly;
	int numOfTriangles;
	//std::vector<int> tri_indices;
	std::map<int , int> polyIndicesMap;
	std::map<int , int> invMap;

	double sourceBoundaryMinArc;
	double sourceArea;
	static bool _first;
	bool isSimple;
public:
	int numOfWantedTriangles;
	bool isTriangultae;
	std::vector<Triangle_2> triangulated;
	Mesh target_mesh;
	std::vector<Kernel::Point_3> pVec;
	std::vector<int> fVec;

	Shor();
	~Shor();
	void load_polygon ( Polygon_2 poly , Polygon_2 boundaryPoly );
	void load_rArray ( int *arr );
	void init_tables();
	void play();
	bool addToTable ( int i , int k , int j );
	void build_triangulation();
	void addTriangle( int i , int j );
	void simplify_triangulation();
	void sort_and_check( int *arr , Polygon_2 &res_poly ,  int* temp );
	void setSourceMinArc ( double min ){this->sourceBoundaryMinArc = min;}
	void setSourceArea ( double a ){this->sourceArea = a;}
	void createSimplePolygonList( std::vector<Polygon_2> &simple_list );
	void findSharedEdgesAndAddPoints(std::vector< std::vector<CGAL::Point_2<Kernel>> >& sharedVec);
	void addSharedEdgeToSimplePoly(Polygon_2& toTheListPoly, std::vector< std::vector<CGAL::Point_2<Kernel>> >& sharedVec, CGAL::Point_2<Kernel> start, CGAL::Point_2<Kernel> end);
};

