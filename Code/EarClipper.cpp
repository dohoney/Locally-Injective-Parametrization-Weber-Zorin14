#include "stdafx.h"

#include "EarClipper.h"



EarClipper::EarClipper(const std::vector<Complex>& polygon, const std::vector<int>& fullRotationIndices)
{
	mNumActiveVertices = polygon.size();
	assert(mNumActiveVertices >= 3);
	assert(mNumActiveVertices == fullRotationIndices.size());

	mFirstEarVertex = -1;

	mPolygon.resize(mNumActiveVertices);

	for(int i = 0; i < mNumActiveVertices; i++)
	{
		Vertex& v_i = mPolygon[i];
		v_i.mFullRotationIndex = fullRotationIndices[i];
		v_i.mPosition = polygon[i];
		v_i.mActive = true;
		v_i.mIsEar = false;
		v_i.mIsEdgeIntersecting = false;

		v_i.mPrevVertexIndex = i > 0 ? i - 1 : mNumActiveVertices - 1;
		v_i.mNextVertexIndex = i < mNumActiveVertices - 1 ? i +1 : 0;
		v_i.mPrevEarIndex = -1;
		v_i.mNextEarIndex = -1;
	}

	initialize();
}


EarClipper::~EarClipper()
{

}


//remove as many as possible ears from the polygon as long as the polygon is not degenerate (means it has only two vertices)
//triangles - a list of triplets of (global. i.e. original) indices
//simplifiedPolygonIndices - the size of this vector equals the size of the simplified polygon. This vector contains the indices of the original (unsimplified) polygon that take part in the simplified polygon
void EarClipper::clipAllEars(std::vector<unsigned int>& triangles, std::vector<unsigned int>& simplifiedPolygonIndices)
{
	triangles.clear();
	triangles.reserve(3*mNumActiveVertices);

	int numEarsRemoved = 0;

	//run as long as there are ears and the polygon is not degenerated
	//while(mFirstEarVertex != -1 && mNumActiveVertices > 2)
	while(mFirstEarVertex != -1 && mNumActiveVertices > 3)
	{
		int prevVertexIndex = mPolygon[mFirstEarVertex].mPrevVertexIndex;
		int nextVertexIndex = mPolygon[mFirstEarVertex].mNextVertexIndex;

		assert(prevVertexIndex >= 0 && mFirstEarVertex >= 0 && nextVertexIndex >= 0);
		
		triangles.push_back(prevVertexIndex);
		triangles.push_back(mFirstEarVertex);
		triangles.push_back(nextVertexIndex);

		removeVertexAndItsEar(mFirstEarVertex);
		numEarsRemoved++;

		bool isPrevEar = isEar(prevVertexIndex);
		bool isNextEar = isEar(nextVertexIndex);

		if(mPolygon[prevVertexIndex].mIsEar && !isPrevEar)
		{
			removeEarFromEarsList(prevVertexIndex);
		}
		if(!mPolygon[prevVertexIndex].mIsEar && isPrevEar)
		{
			insertEarToEarsList(prevVertexIndex);
		}
		if(mPolygon[nextVertexIndex].mIsEar && !isNextEar)
		{
			removeEarFromEarsList(nextVertexIndex);
		}
		if(!mPolygon[nextVertexIndex].mIsEar && isNextEar)
		{
			insertEarToEarsList(nextVertexIndex);
		}
	}

	simplifiedPolygonIndices.clear();
	
	if(mNumActiveVertices >= 3)
	{
		simplifiedPolygonIndices.reserve(mNumActiveVertices);

		int sizeOfOriginalPolygon = mPolygon.size();

		for(int i = 0; i < sizeOfOriginalPolygon; i++)
		{
			if(mPolygon[i].mActive)
			{
				simplifiedPolygonIndices.push_back(i);
			}
		}
	}
}


bool EarClipper::doTwoSegmentsIntersect(const Complex& segmentA_start, const Complex& segmentA_end, const Complex& segmentB_start, const Complex& segmentB_end)
{

	double test1 = CounterClockWise(segmentA_start, segmentA_end, segmentB_start);
	double test2 = CounterClockWise(segmentA_start, segmentA_end, segmentB_end);

	if((test1 > 0.0 && test2 > 0.0) || (test1 < 0.0 && test2 < 0.0))
	{
		return false;
	}

	std::vector<ESegment_2> segments;

	segments.push_back(ESegment_2(EPoint_2(segmentA_start.real(), segmentA_start.imag()), EPoint_2(segmentA_end.real(), segmentA_end.imag())));
	segments.push_back(ESegment_2(EPoint_2(segmentB_start.real(), segmentB_start.imag()), EPoint_2(segmentB_end.real(), segmentB_end.imag())));

	return CGAL::do_curves_intersect(segments.begin(), segments.end());
}


bool EarClipper::isEar(int i)
{
	assert(i >= 0);
	assert(i < mPolygon.size());
	assert(mPolygon[i].mActive);

	int prev_i = mPolygon[i].mPrevVertexIndex;
	int next_i = mPolygon[i].mNextVertexIndex;

	if(mPolygon[i].mIsEdgeIntersecting)
	{
		return false;
	}

	if(mPolygon[prev_i].mIsEdgeIntersecting)
	{
		return false;
	}

	double cross = CounterClockWise(mPolygon[prev_i].mPosition, mPolygon[i].mPosition, mPolygon[next_i].mPosition);

	if(cross <= 0.0)
	{
		return false;
	}

	if(mPolygon[i].mFullRotationIndex != 0)
	{
		return false;
	}
	if(mPolygon[prev_i].mFullRotationIndex != 0)
	{
		return false;
	}
	if(mPolygon[next_i].mFullRotationIndex != 0)
	{
		return false;
	}

	if(mNumActiveVertices == 3) //if this is a triangle (it will also be oriented correctly), all three vertices are ear tips
	{
		return true;
	}

	int index = i;

	for(int j = 0; j < mNumActiveVertices; j++, index = mPolygon[index].mNextVertexIndex)
	{
		if(index == i || index == prev_i || index == next_i)
		{
			continue;
		}

		if(CounterClockWise(mPolygon[i].mPosition, mPolygon[next_i].mPosition, mPolygon[index].mPosition) >= 0.0 &&
		   CounterClockWise(mPolygon[prev_i].mPosition, mPolygon[i].mPosition, mPolygon[index].mPosition) >= 0.0 &&
		   CounterClockWise(mPolygon[next_i].mPosition, mPolygon[prev_i].mPosition, mPolygon[index].mPosition) >= 0.0)
		{
			return false;
		}
	}

	return true; //if we got here it means that we passed all tests and the vertex is an ear tip
}


void EarClipper::initialize()
{
	for(int i = 0; i < mNumActiveVertices; i++)
	{
		mPolygon[i].mIsEdgeIntersecting = false;
	}

	//initialize edge intersection states
	for(int i = 0; i < mNumActiveVertices; i++)
	{
		if(mPolygon[i].mIsEdgeIntersecting)
		{
			continue; //if the current segment already intersects another segment, there is no need to look for intersections
		}
		int next_i = mPolygon[i].mNextVertexIndex;

		int index = i;

		for(int j = 0; j < mNumActiveVertices; j++)
		{
			if(i == j)
			{
				continue; //there is no point to check whether a segment intersect itself
			}
			int next_j = mPolygon[j].mNextVertexIndex;

			if(doTwoSegmentsIntersect(mPolygon[i].mPosition, mPolygon[next_i].mPosition, mPolygon[j].mPosition, mPolygon[next_j].mPosition))
			{
				mPolygon[i].mIsEdgeIntersecting = true;
				mPolygon[j].mIsEdgeIntersecting = true;
				break; //no need to continue looking for intersections since one intersection is enough to turn on the flag
			}
		}
	}

	mFirstEarVertex = -1; //there are no ears yet

	//initialize ears
	for(int i = 0; i < mNumActiveVertices; i++)
	{
		bool isThisEar = isEar(i);
		if(isThisEar)
		{
			assert(!mPolygon[i].mIsEar);
			insertEarToEarsList(i);
			mPolygon[i].mIsEar = true;
		}
	}
}


void EarClipper::insertEarToEarsList(int i)
{
	assert(mPolygon[i].mActive);

	mPolygon[i].mPrevEarIndex = -1; //we always insert at the beginning of the list, hence the prev of the newly inserted ear is always -1

	if(mFirstEarVertex == -1) //this is the first ear found
	{
		mPolygon[i].mNextEarIndex = -1;
	}
	else //the list of ears already contain some ears
	{
		mPolygon[i].mNextEarIndex = mFirstEarVertex;
		mPolygon[mFirstEarVertex].mPrevEarIndex = i;
	}
	mFirstEarVertex = i;
	mPolygon[i].mIsEar = true; //mark this vertex as ear tip
}


void EarClipper::removeEarFromEarsList(int i)
{
	assert(mPolygon[i].mIsEar); //if we want to make this vertex "NOT ear", it has to be an ear first 
	assert(mFirstEarVertex != -1);

	int prevEarIndex = mPolygon[i].mPrevEarIndex;
	int nextEarIndex = mPolygon[i].mNextEarIndex;

	if(prevEarIndex == -1) //if this ear doesn't have a prev, it means that it is the first ear in the ears list
	{
		assert(i == mFirstEarVertex); //this should always be the case
		mFirstEarVertex = nextEarIndex;
	}
	else //this ear is not the first element in the ears list
	{
		mPolygon[prevEarIndex].mNextEarIndex = nextEarIndex;
	}

	if(nextEarIndex == -1) //if this ear doesn't have a next, it means that it is the last ear in the ears list
	{
		//do nothing
	}
	else
	{
		mPolygon[nextEarIndex].mPrevEarIndex = prevEarIndex;
	}

	mPolygon[i].mIsEar = false; //not an ear anymore
}


void EarClipper::removeVertexAndItsEar(int i)
{
	assert(mNumActiveVertices >= 3); //assuming the polygon is not degenerated
	assert(mPolygon[i].mActive); //only active vertex can be removed
	assert(mPolygon[i].mIsEar); //we only remove a vertex from the polygon if it is an ear tip

	int prevVertexIndex = mPolygon[i].mPrevVertexIndex;
	int nextVertexIndex = mPolygon[i].mNextVertexIndex;

	mPolygon[prevVertexIndex].mNextVertexIndex = nextVertexIndex;
	mPolygon[nextVertexIndex].mPrevVertexIndex = prevVertexIndex;

	//this should always be the case since the vertex we remove is an ear tip, it means that the diagonal from v_prev to v_next does not intersect anything else
	assert(mPolygon[prevVertexIndex].mIsEdgeIntersecting == false);

	removeEarFromEarsList(i);

	mPolygon[i].mActive = false;

	mNumActiveVertices--;
}


EarClipper::Vertex::Vertex() : 
	mPosition(0.0, 0.0),
	mActive(false),
	mIsEar(false),
	mIsEdgeIntersecting(false),
	mPrevVertexIndex(-1),
	mNextVertexIndex(-1),
	mPrevEarIndex(-1),
	mNextEarIndex(-1),
	mFullRotationIndex(0)
{

}

