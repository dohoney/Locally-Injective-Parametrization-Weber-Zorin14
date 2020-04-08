#pragma once


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of an extension to the ear clipping algorithm by David Eberly for self-overlapping polygons.
// Note: the algorithm only simplify the polygon by removing some ears but not all of them.
// Once the algorithm is done, we can run our triangulation algorithm for self-overlapping polygons on the simplified polygon which
// is expected to run faster
// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include "STL_Macros.h"
#include "CGAL_Mesh.h"
//#include "CGAL_Macros.h"


class EarClipper
{
	class Vertex
	{
	public:
		Vertex();
		Complex mPosition;
		bool mActive; //inactive vertex is a vertex that is not part of the polygon
		bool mIsEar;
		bool mIsEdgeIntersecting; //true if the edge between this vertex and the next one intersect any other polygon edge or vertex
		int mPrevVertexIndex, mNextVertexIndex; // links to the previous and next active vertices in the polygon
		int mPrevEarIndex, mNextEarIndex; // links to the previous and next ear vertices in the polygon
		int mFullRotationIndex;
	};


public:

	EarClipper(const std::vector<Complex>& polygon, const std::vector<int>& fullRotationIndices);
	~EarClipper();
	void clipAllEars(std::vector<unsigned int>& triangles, std::vector<unsigned int>& simplifiedPolygonIndices);

protected:
	
	void insertEarToEarsList(int i);
	void removeEarFromEarsList(int i);
	void removeVertexAndItsEar(int i);
	bool isEar(int i);
	bool doTwoSegmentsIntersect(const Complex& segmentA_start, const Complex& segmentA_end, const Complex& segmentB_start, const Complex& segmentB_end);
	void initialize();


protected:

	std::vector<Vertex> mPolygon; //the vertices of the input polygon ordered in counter-clockwise direction (the list is cyclic)
	int mFirstEarVertex;
	int mNumActiveVertices;
};

