// last update - stop counting
#include "stdafx.h"
#include "EarClipper.h"
#define FINAL 1		// for debug (search me in the simplify function)

double complexLength ( std::complex<double> t1  )
{
	return ( std::sqrt (t1.real()*t1.real() + t1.imag()*t1.imag() ) );
}

double complexLengthVector ( std::complex<double> t1 , std::complex<double> t2 )
{
	std::complex<double> t3;
	t3.real(t2.real()-t1.real());
	t3.imag(t2.imag()-t1.imag());
	return ( complexLength(t3) );
}

bool searchForShare( std::complex<double> &t1, std::complex<double> &t2, std::vector<std::complex<double>> &in_points, int &index)
{
	for ( int i = 0; i < (int)in_points.size(); i=i+2 )
	{
		if ( (t1 == in_points[i])&&(t2 == in_points[i+1]) || (t2 == in_points[i])&&(t1 == in_points[i+1]) )
		{
			index = i;
			return (true);
		}	
	}
	return (false);
}

bool Shor::_first = true;

Shor::Shor()
{
	this->N = 0;
	this->Q = NULL;
	this->K = NULL;
	this->firstAng = NULL;
	this->lastAng = NULL;
	this->rArray = NULL;
	this->isSimple = false;
	this->sourceBoundaryMinArc = -1;
	this->numOfWantedTriangles = -1;
}

Shor::~Shor()
{
	if ( this -> N == 0 )
		return;

	for ( int i = 0 ; i < N ; ++i )
	{
		delete []Q[i];
		delete []K[i];
		delete []firstAng[i];
		delete []lastAng[i];
	}
	delete []Q;
	delete []K;
	delete []firstAng;
	delete []lastAng;

	if ( rArray != NULL )
		delete []rArray;
}

void Shor::load_polygon( Polygon_2 poly , Polygon_2 boundaryPoly )
{
	this -> poly = poly ;
	this -> N = poly.size();
	this -> bPoly = boundaryPoly;

	Q = new int*[N];
	K = new int*[N];
	firstAng = new Angle*[N];
	lastAng = new Angle*[N];
	int j = 0, bN = (int)boundaryPoly.size();

	for ( int i = 0 ; i < N ; ++i )
	{
		Q[i] = new int[N]();
		K[i] = new int[N]();
		firstAng[i] = new Angle[N];
		lastAng[i] = new Angle[N];

		Kernel::Point_3 point( poly[i][0], poly[i][1], 0 );
		this->pVec.push_back(point);

		this->polyIndicesMap[i] = j;
		this->invMap[j] = i;

		while ( poly[(i+1)%N] != boundaryPoly[j] )
		{
			j = (j+1)%bN;
		}
	}

	for ( int i = 0 ; i < N ; ++i )
		for ( int j = 0 ; j < N ; ++j )
			K[i][j] = -1;
}

void Shor::load_rArray ( int *arr )
{
	rArray = new int [ N ]();
	rVector.resize(N);
	if (arr != NULL)
	{
		for (int i = 0; i < N; i++)
		{
			rArray[i] = arr[i];
			rVector[i] = arr[i];
		}
	}
	else
	{
		for (int i = 0; i < N; i++)
		{
			rArray[i] = 1;
			rVector[i] = 1;
		}
	}
}

void Shor::init_tables()
{
	Point_2 a;
	Point_2 b;
	Point_2 c;
	N = poly.size();
	for ( int i = 0 ; i < N ; ++i ) // init
	{
		a = poly[i];
		b = poly[(i+1)%N];
		c = poly[(i+2)%N];

		Q[i][(i+1)%N] = 1;
		firstAng[i][(i+1)%N] = Angle ( b , a , b );
		lastAng[i][(i+1)%N] = Angle ( a , b , a );
		
		Triangle_2 tri(a , b , c );

		if ( tri.orientation() > 0 )
		{
			Q[i][(i+2)%N] = 1;
			K[i][(i+2)%N] = (i + 1)%N;
			firstAng[i][(i+2)%N] = Angle( b , a , c ); //[k , i , j]
			lastAng[i][(i+2)%N] = Angle ( a , c , b ); //[i , j , k]
		}
		else
			Q[i][(i+2)%N] = 0;
	}
}

void Shor::play()
{
	//if (_first)
	//	std::cout << "Running Shor algoritem...\n";
	int i , j , k , d; // i - start index, k- spliting index , j- end index, d- distance

	if (this->poly.is_simple())	// if target polygon is simple we dont need Shor algorithem
	{
		this->isSimple = true;
		return;
	}

	//------------ear clipping------------------------------------------------------
	std::vector < std::complex<double> > polygonP;
	for (int i = 0; i < N; ++i)
		polygonP.push_back(std::complex<double>(poly[i][0],poly[i][1]));

	EarClipper earClipper(polygonP, rVector);
	std::vector<unsigned int> simplifiedPolygonIndices;
	earClipper.clipAllEars(triangles, simplifiedPolygonIndices);
	int numVerticesInSimplifiedPolygon = simplifiedPolygonIndices.size();
	std::vector<int> simplifiedSOPfullRotationIndices;
	for (int i = 0; i < numVerticesInSimplifiedPolygon; i++)
	{
		simpToOriginalIndices[i] = simplifiedPolygonIndices[i];
		simpPoly.push_back( Point_2( polygonP [ simplifiedPolygonIndices[i] ].real(), polygonP [ simplifiedPolygonIndices[i] ].imag()) );
		simplifiedSOPfullRotationIndices.push_back(rVector[simplifiedPolygonIndices[i]]);
	}
	rVector = simplifiedSOPfullRotationIndices;
	tempForSwap = poly;
	poly = simpPoly;
	//-------------------------------------------------------------------------------
	this->init_tables();
	N = poly.size();
	bool stop = false;
	for ( d = 3 ; d < N ; ++d )			//sweep the distances
	{
		for ( i = 0 ; i < N ; ++i )		//sweep the vertices
		{
			j = (i + d)%N;
			k = (i + 1)%N;
			while ( k != j )			//sweep from v_i to v_j
			{
				if ( Q[i][k] == 1 && Q[k][j] == 1 )
				{
					if  ( this->addToTable( i , k , j ) )
					{
						Q[i][j] = 1;
						K[i][j] = k;

						if ( i == (j + 1)%N )
						{
							stop = true;
							break;
						}
					}
				}
				k = (k + 1)%N;
			}
			if ( stop )
				break;
		}
		if ( stop )
			break;
	}
	//if (_first)
	//{
	//	std::cout << "Done!\n";
	//	_first = false;
	//}
}

bool Shor::addToTable ( int i , int k , int j )
{
	Point_2 a,b,c,helpP;
	Angle temp1,temp2,tempF,tempL;
	N = poly.size();
	a = poly[i];
	b = poly[k];
	c = poly[j];

	Triangle_2 tri(a , b , c );

	Angle k1 ( poly[ j ] , poly[ k ] , poly[ i ] );
	Angle k2 ( poly[ i ] , poly[ k ] , poly[ ( k-1+N )%N ] );
	Angle k3 ( poly[ ( k-1+N )%N ] , poly[ k ] , poly[ ( k+1 )%N ] );
	Angle k4 ( poly[ ( k+1 )%N ] , poly[ k ] , poly[ j ] );
	Angle kSum = k1 + k2 + k3 + k4;

	bool triFlag = true;
	helpP = poly[(i+1)%N];
	Triangle_2 triI1(a , helpP , b );//[i,i+1,k]
	Triangle_2 triI2(a , c , helpP );//[i,j,i+1]
	if ( triI1.orientation() < 0 )
		if ( triI2.orientation() < 0 )
			triFlag = false;

	helpP = poly[(j-1+N)%N];
	Triangle_2 triJ1(c , b , helpP );//[j,k,j-1]
	Triangle_2 triJ2(a , c , helpP );//[i,j,j-1]
	if ( triJ1.orientation() < 0 )
		if ( triJ2.orientation() < 0 )
			triFlag = false;

	// extend algoritem: check the rotation of the angles
	temp1 = Angle ( b , a , c ); // [k,i,j]
	temp2 = Angle ( a , c , b ); // [i,j,k]
	tempF = firstAng[i][k] + temp1;
	tempL = temp2 + lastAng[k][j];

	//update angle table
	firstAng[i][j] = tempF;	
	lastAng[i][j] = tempL;

	bool kAngleFlag = (kSum.getR() == 0);

	if ( ( tempF.getR() <= rVector[i] ) && ( tempL.getR() <= rVector[j] ) && ( (k4+k1+k2).getR() == rVector[k] ) )
						kAngleFlag = true;
					/*else
						triFlag = false;*/
	if  ( ( tri.orientation() > 0 ) && ( kAngleFlag ) && ( triFlag ) )
		return true;
	else
		return false;
}

void Shor::build_triangulation()
{
	if (this->isSimple)
	{
		std::vector<std::complex<double> > polygonPoints;
		std::vector<unsigned int> triangleIndices;

		for (int i = 0; i < (int)poly.size(); ++i)
			polygonPoints.push_back(std::complex<double>(poly[i].x(), poly[i].y()));

		triangulatePolygonWithoutAddingVertices( polygonPoints , triangleIndices );
		for (int i = 0; i < (int)triangleIndices.size(); ++i)
			this->fVec.push_back(triangleIndices[i]);

		this->isTriangultae = true;
		//---------------
		MeshBuilder<Mesh::HalfedgeDS, Kernel> meshBuilder(&pVec, &fVec);
		this->target_mesh.clear();
		this->target_mesh.delegate(meshBuilder);
		this->target_mesh.updateAllGlobalIndices();
		return;
	}

	int i;
	bool build = false;
	this->numOfTriangles = 0;
	N = poly.size();
	for ( i = 0 ; i < N ; i++)
	{
		if ( Q[i][ (i-1+N)%N ] == 1 )
		{
			build = true;
			break;
		}
	}

	if (build)
	{
		this->addTriangle( i , (i - 1 + N)%N );
		this->isTriangultae = true;
		for (int i = 0; i < fVec.size(); ++i)
			fVec[i] = simpToOriginalIndices[fVec[i]];
		fVec.insert(fVec.end(), triangles.begin(), triangles.end());
		
		/*/---------------
		//DUBUG in matlab
		GMMDenseColMatrix pp(pVec.size(), 2);
		GMMDenseColMatrix ff(fVec.size()/3, 3);
		for (int i = 0; i < (int)pVec.size(); ++i)
		{
			pp(i, 0) = pVec[i][0];
			pp(i, 1) = pVec[i][1];
		}
		int index = 0;
		for (int i = 0; i < (int)fVec.size(); i=i+3)
		{
			ff(index, 0) = fVec[i];
			ff(index, 1) = fVec[i+1];
			ff(index, 2) = fVec[i+2];
			index++;
		}
		MatlabGMMDataExchange::SetEngineDenseMatrix("pp", pp);
		MatlabGMMDataExchange::SetEngineDenseMatrix("ff", ff);*/

		MeshBuilder<Mesh::HalfedgeDS,Kernel> meshBuilder( &pVec, &fVec );
		this->target_mesh.clear();
		this->target_mesh.delegate( meshBuilder );
		this->target_mesh.updateAllGlobalIndices();

		poly = tempForSwap;
		N = poly.size();
		//----------------
	}
	else
		this->isTriangultae = false;
}

void Shor::addTriangle( int i , int j )
{
	if ( K[i][j] == -1 )	//( ( std::abs( i - j ) < 2 ) && ( i != (j+1)%N ) )
		return;

	int k = K[i][j];

	triangulated.push_back( Triangle_2( poly[i] , poly[k] , poly[j] ) );

	//tri_indices.push_back ( i );
	//tri_indices.push_back ( k );
	//tri_indices.push_back ( j );
	this->fVec.push_back(/*this->polyIndicesMap[i]*/i);
	this->fVec.push_back(/*this->polyIndicesMap[k]*/k);
	this->fVec.push_back(/*this->polyIndicesMap[j]*/j);

	this->numOfTriangles++;
	addTriangle( i , k );
	addTriangle( k , j );
}

void simplify_mesh(Mesh& target_mesh)
{
		bool it_flag = true;
		std::vector<int> bad_he;
		Mesh::Halfedge_iterator heIt,split_edge;
		Mesh::Halfedge_around_facet_circulator circHE;
		Polygon_2 temp_poly;
		heIt = target_mesh.halfedges_begin();
		bool /*resetItFlag = true ,*/ badFlag = false; 
		while ( heIt != target_mesh.halfedges_end() )
		{
			if (!badFlag)
			{
				heIt = target_mesh.halfedges_begin();
				badFlag = false;
			}
			temp_poly.clear();
			while ( ( heIt->is_border() ) || ( heIt->opposite()->is_border() ) )
			{
				if ( heIt == target_mesh.halfedges_end() )
				{
					it_flag = false;
					break;
				}
				heIt++;
			}
			if (!it_flag)
				break;

			for ( int i = 0; i < bad_he.size(); ++i)
				if ( heIt->index() == bad_he[i] )
					badFlag = true;

			if (badFlag)
			{
				heIt++;
				continue;
			}
			auto g = heIt->prev_on_vertex();
			bad_he.push_back(heIt->index());
			split_edge = target_mesh.join_facet( heIt );
			auto face = split_edge->facet();
			circHE = face->facet_begin();
			auto end = circHE;
			temp_poly.push_back( CGAL::Point_2<Kernel>(circHE->vertex()->point()[0] , circHE->vertex()->point()[1]  ) );
			circHE++;
			while ( circHE != end)
			{
				temp_poly.push_back( CGAL::Point_2<Kernel>(circHE->vertex()->point()[0] , circHE->vertex()->point()[1]  ) );
				circHE++;
			}
			if ( !temp_poly.is_simple() )
			{
				target_mesh.split_facet(split_edge,g);
				heIt++;
				//resetItFlag = false;
			}
			else
				bad_he.pop_back();

		}

}

void Shor::createSimplePolygonList( std::vector<Polygon_2> &simple_list )
{
	Polygon_2 temp_poly , toTheListPoly;
	Mesh::Face_iterator f_it = target_mesh.facets_begin();
	Mesh::Halfedge_around_facet_circulator circHE;

	std::vector< std::vector<CGAL::Point_2<Kernel>> > sharedVec;
	findSharedEdgesAndAddPoints(sharedVec);

	while ( f_it != target_mesh.facets_end() )
	{
		temp_poly.clear();
		toTheListPoly.clear();
		circHE = f_it->facet_begin();
		auto end = circHE;
		temp_poly.push_back( CGAL::Point_2<Kernel>(circHE->vertex()->point()[0] , circHE->vertex()->point()[1]  ) );
		circHE++;
		while ( circHE != end)
		{
			temp_poly.push_back( CGAL::Point_2<Kernel>(circHE->vertex()->point()[0] , circHE->vertex()->point()[1]  ) );
			circHE++;
		}
		
		int i = 0 , j = 0 , tempSize = (int)temp_poly.size() ;
		while ( temp_poly[0] != this->bPoly[i] )
			i++;
		for ( int k = 0; k < tempSize; ++k )
		{
			if ( temp_poly[(k+1)%tempSize] == this->bPoly[ this->polyIndicesMap [( (this->invMap[i])+1 ) % this->N ] ] )	//if this is a boundary edge - take all the points on the boundary
			{
				while ( this->bPoly[i] != temp_poly[(k+1)%tempSize] )
				{
					toTheListPoly.push_back ( this->bPoly[i] );
					i = (i+1)%this->bPoly.size();
				}
			}
			else			// this is shared edge
			{
				//toTheListPoly.push_back ( temp_poly[k] );
				addSharedEdgeToSimplePoly(toTheListPoly, sharedVec, temp_poly[k], temp_poly[(k + 1) % tempSize]);
				//maybe need to add code that insert more points here --- done!

				i = 0;
				while ( temp_poly[(k+1)%tempSize] != this->bPoly[i] )
					i = (i+1)%this->bPoly.size();
			}

		}


		simple_list.push_back( toTheListPoly );
		f_it++;
	}
}

void Shor::findSharedEdgesAndAddPoints(std::vector< std::vector<CGAL::Point_2<Kernel>> >& sharedVec)
{
	/*double minTargetEdge = (bPoly[1] - bPoly[0]).squared_length();
	int size = bPoly.size();
	for (int i = 1; i < size; ++i)
	{
		if (minTargetEdge >(bPoly[i] - bPoly[(i + 1) % size]).squared_length())
			minTargetEdge = (bPoly[i] - bPoly[(i + 1) % size]).squared_length();
	}
	minTargetEdge = std::sqrt(minTargetEdge);*/
	// instead of "min" i calculate avarage
	double minTargetEdge = std::sqrt( (bPoly[1] - bPoly[0]).squared_length() ),temp;
	int size = bPoly.size();
	for (int i = 1; i < size; ++i)
	{
		temp = std::sqrt( (bPoly[(i + 1) % size] - bPoly[i]).squared_length() );
		minTargetEdge += temp;
	}
	minTargetEdge = minTargetEdge / size;

	Mesh::Edge_iterator e_it = target_mesh.edges_begin();
	std::vector<CGAL::Point_2<Kernel>> pointsVec;

	while (e_it != target_mesh.edges_end())
	{
		if (!e_it->is_border_edge())
		{
			pointsVec.clear();
			CGAL::Point_2<Kernel> start =	CGAL::Point_2<Kernel>(e_it->vertex()->point()[0], e_it->vertex()->point()[1]);
			CGAL::Point_2<Kernel> end	=	CGAL::Point_2<Kernel>(e_it->opposite()->vertex()->point()[0], e_it->opposite()->vertex()->point()[1]);

			pointsVec.push_back(start);

			double len = (end - start).squared_length();
			len = std::sqrt(len);
			len = len / minTargetEdge;
			int num = len;

			for (int i = 1; i < num; ++i)
			{
				pointsVec.push_back(CGAL::Point_2<Kernel>(start[0] + (end[0] - start[0]) * i/num, start[1] + (end[1] - start[1]) * i/num));
			}
			pointsVec.push_back(end);

			std::vector<CGAL::Point_2<Kernel>> revPointsVec(pointsVec.rbegin(), pointsVec.rend());
			sharedVec.push_back(pointsVec);
			sharedVec.push_back(revPointsVec);

		}
		e_it++;
	}

}

void Shor::addSharedEdgeToSimplePoly(Polygon_2& toTheListPoly, std::vector< std::vector<CGAL::Point_2<Kernel>> >& sharedVec, CGAL::Point_2<Kernel> start, CGAL::Point_2<Kernel> end)
{
	int vectorSize = sharedVec.size();
	for (int i = 0; i < vectorSize; ++i)
	{
		int L = sharedVec[i].size();
		if ((start == sharedVec[i][0]) && (end == sharedVec[i][L - 1]))
		{
			for (int j = 0; j < L - 1; ++j)
				toTheListPoly.push_back(sharedVec[i][j]);
			return;
		}
	}
	assert(0);	// if we got here its mean the shared edge is not at the shared vector


}

void sendToMatlab(std::vector<std::complex<double>>& meshVertices, std::vector<unsigned int>& triangleIndices)
{
	//DUBUG in matlab
	GMMDenseColMatrix pp(meshVertices.size(), 2);
	GMMDenseColMatrix ff(triangleIndices.size() / 3, 3);

	for (int i = 0; i < (int)meshVertices.size(); ++i)
	{
		pp(i, 0) = meshVertices[i].real();
		pp(i, 1) = meshVertices[i].imag();
	}

	int index = 0;
	for (int i = 0; i < (int)triangleIndices.size(); i = i + 3)
	{
		ff(index, 0) = triangleIndices[i];
		ff(index, 1) = triangleIndices[i+1];
		ff(index, 2) = triangleIndices[i+2];
		index++;
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("pp", pp);
	MatlabGMMDataExchange::SetEngineDenseMatrix("ff", ff);
	MatlabInterface::GetEngine().Eval("ff=ff+1;trimesh( ff , pp(:,1) , pp(:,2) );hold on");

}

void Shor::simplify_triangulation()
{
	if ( !this->isTriangultae )
		return;
	std::cout << "Building mesh from target polygon...\n";
	if (!this->isSimple)
		simplify_mesh(target_mesh);

	//stage 2 of simplification : use the triangle alogoritem//

	std::vector<Polygon_2> simple_list;
	if (this->isSimple)
		simple_list.push_back(bPoly);
	else
		createSimplePolygonList( simple_list );


	/*Mesh::Edge_iterator e_it = target_mesh.edges_begin();
	std::vector<std::complex<double>> in_points;
	while ( e_it != target_mesh.edges_end() )
	{
		if ( !e_it->is_border_edge() )
		{
			in_points.push_back (std::complex<double> ( e_it->vertex()->point()[0] , e_it->vertex()->point()[1] ) );
			in_points.push_back (std::complex<double> (e_it->opposite()->vertex()->point()[0],e_it->opposite()->vertex()->point()[1] ));
		}
		e_it++;
	}
	std::vector<int> shareEdgesFlag;
	for ( int i = 0; i < (int)in_points.size()/2; ++i)
		shareEdgesFlag.push_back(-1);
	std::vector < std::vector<std::complex<double>> > shareVertices;*/

	//int maxPolygonPoints=0 ,maxTriIndices=0;
	std::vector < std::vector<std::complex<double>> >  polygonPoints_vec , meshVertices_vec;
	std::vector < std::vector<std::pair<int, std::complex<double> > > > boundaryVertices_vec;
	std::vector < std::vector<unsigned int> > triangleIndices_vec;

	std::vector<std::complex<double>> polygonPoints;
	std::vector<std::pair<int, std::complex<double> > > boundaryVertices;
		
	//std::complex<double> t1,t2,t_temp;
	//double L , delta , dis;
	//int add = 0,index;
	double targetArea = 0;//simple_list[0].area();

	// add points on boundery //
	/*for ( int i = 0 ; i < (int)simple_list.size() ; i++ )
	{
		targetArea += simple_list[i].area();
		for ( int j = 0; j < (int)simple_list[i].size() ; ++j )
		{	
			polygonPoints.push_back ( std::complex<double> ( simple_list[i][j][0] , simple_list[i][j][1] ) );
			boundaryVertices.push_back ( std::pair<int,std::complex<double>> ( j+add , polygonPoints[j] ) );
			t1 = std::complex<double> ( simple_list[i][j][0] , simple_list[i][j][1] );
			t2 = std::complex<double> ( simple_list[i][(j+1)%simple_list[i].size()][0] , simple_list[i][(j+1)%simple_list[i].size()][1] );
			if ( searchForShare(t1,t2,in_points,index) )
			{
				if ( shareEdgesFlag[index/2] == -1 )
				{
					std::vector<std::complex<double>> temp_vec;
					L = complexLengthVector(t1,t2);
					delta = this->sourceBoundaryMinArc / L;
					dis = complexLengthVector(t1,t2);
					int plus = 1;
					temp_vec.push_back(t1);
					while ( dis > this->sourceBoundaryMinArc )
					{
						t_temp.real( t1.real() - t1.real()*delta*plus + t2.real()*delta*plus );
						t_temp.imag( t1.imag() - t1.imag()*delta*plus + t2.imag()*delta*plus );
						add++;
						plus++;
						polygonPoints.push_back ( t_temp );
						temp_vec.push_back(t_temp);
						boundaryVertices.push_back ( std::pair<int,std::complex<double>> ( j+add , t_temp ) );
						dis = complexLengthVector(t_temp,t2);
					}
					temp_vec.push_back(t2);
					shareEdgesFlag[index/2] = (int)shareVertices.size();
					shareVertices.push_back(temp_vec);	
					
				}
				else
				{
					if ( t1 == shareVertices[ shareEdgesFlag[index/2] ][0] )
					{
						for ( int ii = 1; ii < (int)shareVertices[ shareEdgesFlag[index/2] ].size()-1; ++ii)
						{
							add++;
							polygonPoints.push_back ( shareVertices[ shareEdgesFlag[index/2] ][ii] );
							boundaryVertices.push_back ( std::pair<int,std::complex<double>> ( j+add , shareVertices[ shareEdgesFlag[index/2] ][ii] ) );
						}
					}
					else
					{
						for ( int ii = (int)shareVertices[ shareEdgesFlag[index/2] ].size()-2; ii > 0; --ii )
						{
							add++;
							polygonPoints.push_back ( shareVertices[ shareEdgesFlag[index/2] ][ii] );
							boundaryVertices.push_back ( std::pair<int,std::complex<double>> ( j+add , shareVertices[ shareEdgesFlag[index/2] ][ii] ) );
						}
					}
					
				}
			}
			else
			{
				L = complexLengthVector(t1,t2);
				delta = this->sourceBoundaryMinArc / L;
				dis = complexLengthVector(t1,t2);
				int plus = 1;
				while ( dis > this->sourceBoundaryMinArc )
				{
					t_temp.real( t1.real() - t1.real()*delta*plus + t2.real()*delta*plus );
					t_temp.imag( t1.imag() - t1.imag()*delta*plus + t2.imag()*delta*plus );
					add++;
					plus++;
					polygonPoints.push_back ( t_temp );
					boundaryVertices.push_back ( std::pair<int,std::complex<double>> ( j+add , t_temp ) );
					dis = complexLengthVector(t_temp,t2);
				}
			}
		}

		polygonPoints_vec.push_back(polygonPoints);
		boundaryVertices_vec.push_back(boundaryVertices);
		polygonPoints.clear();
		boundaryVertices.clear();
		add = 0;
	}*/

	//*******************************
	for ( int i = 0 ; i < (int)simple_list.size() ; i++ )
	{
		if (simple_list[i].area() > 0)
			targetArea += simple_list[i].area();
		else
			targetArea += (-1*simple_list[i].area());
		for ( int j = 0; j < (int)simple_list[i].size() ; ++j )
		{
			polygonPoints.push_back ( std::complex<double> ( simple_list[i][j][0] , simple_list[i][j][1] ) );
			boundaryVertices.push_back ( std::pair<int,std::complex<double>> ( j , polygonPoints[j] ) );

		}

		polygonPoints_vec.push_back(polygonPoints);
		boundaryVertices_vec.push_back(boundaryVertices);
		polygonPoints.clear();
		boundaryVertices.clear();
	}
	//*******************************
	
	double avgTriArea = targetArea / this->numOfWantedTriangles;//std::pow(this->sourceBoundaryMinArc,2) * std::sqrt((double)3) / 4;
	double maxTriangleArea = 1.7*avgTriArea;//0.01 * avgTriArea * this->sourceArea;// / 100;// targetArea;
	//std::complex<double> temp_p;
	std::vector<std::complex<double>> meshVertices;
	std::vector<unsigned int> triangleIndices;

	for ( int i = 0; i < (int)polygonPoints_vec.size(); ++i )
	{
		triangulatePolygon(polygonPoints_vec[i],meshVertices,triangleIndices,boundaryVertices_vec[i],maxTriangleArea,false);
		meshVertices_vec.push_back(meshVertices);

		//sendToMatlab(meshVertices, triangleIndices);

		//if ( meshVertices.size() > maxPolygonPoints )
		//	maxPolygonPoints = meshVertices.size();

		meshVertices.clear();
		triangleIndices_vec.push_back(triangleIndices);
		//if ( triangleIndices.size() > maxTriIndices )
		//	maxTriIndices = triangleIndices.size();
		triangleIndices.clear();
	}		
		

	pVec.clear();
	fVec.clear();
		
	/*for ( int i = 0; i < (int)meshVertices_vec[0].size(); ++i )
		pVec.push_back(Kernel::Point_3(meshVertices_vec[0][i].real(),meshVertices_vec[0][i].imag(),0));
	for ( int j = 0; j < (int)triangleIndices_vec[0].size(); ++j )
		fVec.push_back(triangleIndices_vec[0][j]);

		
	bool sFlag = false;
	for ( int i = 1; i < (int)meshVertices_vec.size(); ++i )
	{
		std::vector<int> temp_indices;
		for ( int a =0; a<triangleIndices_vec[i].size(); ++a)
			temp_indices.push_back( (int)triangleIndices_vec[i][a] );
		int checkSize = pVec.size();

		for ( int j = 0; j < (int)meshVertices_vec[i].size(); ++j)
		{
				
			sFlag = false;
			for ( int k = 0; k < checkSize; ++k )
			{
				if ( pVec[k] == Kernel::Point_3(meshVertices_vec[i][j].real(),meshVertices_vec[i][j].imag(),0) )
				{
					sFlag = true;
					for ( int l = 0; l < (int)temp_indices.size(); ++l )
						if ( temp_indices[l] == j )
							temp_indices[l] = (k+1)*-1;
					break;	
				}
			}
			if ( !sFlag )
			{
				pVec.push_back(Kernel::Point_3(meshVertices_vec[i][j].real(),meshVertices_vec[i][j].imag(),0));
				//mesh_mat(j,i) = meshVertices_vec[i][j];
				//jump++;
				int newIndex = pVec.size() - 1;
				for ( int m = 0; m < (int)temp_indices.size(); ++m )
				{
					if (temp_indices[m] == j)
						temp_indices[m] = newIndex;	//pVec.size()-1;
				}
			}
		}

		// update fVec
		for ( int m = 0; m < (int)temp_indices.size(); ++m )
		{
			if (temp_indices[m]<0)
				temp_indices[m] = (temp_indices[m]+1)*-1;
			fVec.push_back( temp_indices[m] );
		}
	}*/

	std::map<Kernel::Point_3, int> pointToRealIndex;
	/*for (int i = 0; i < (int)meshVertices_vec[0].size(); ++i)
	{
		pVec.push_back(Kernel::Point_3(meshVertices_vec[0][i].real(), meshVertices_vec[0][i].imag(), 0));
		pointToRealIndex[pVec[i]] = i;
	}
	for (int j = 0; j < (int)triangleIndices_vec[0].size(); ++j)
		fVec.push_back(triangleIndices_vec[0][j]);*/
	for (int i = 0; i < (int)this->bPoly.size(); ++i)
	{
		pVec.push_back(Kernel::Point_3(bPoly[i][0], bPoly[i][1], 0));
		pointToRealIndex[pVec[i]] = i;
	}

	for (int i = 0; i < (int)meshVertices_vec.size(); ++i)
	{
		std::map<int, int> localToRealIndices;
		for (int j = 0; j < (int)meshVertices_vec[i].size(); ++j)
		{
			Kernel::Point_3 p(meshVertices_vec[i][j].real(), meshVertices_vec[i][j].imag(), 0);
			if (pointToRealIndex.count(p) == 1)
				localToRealIndices[j] = pointToRealIndex[p];
			else
			{
				pVec.push_back(p);
				pointToRealIndex[p] = pVec.size() - 1;
				localToRealIndices[j] = pVec.size() - 1;
			}
		}
		for (int j = 0; j < triangleIndices_vec[i].size(); ++j)
			fVec.push_back(localToRealIndices[ triangleIndices_vec[i][j] ]);
	}



	GMMDenseComplexColMatrix mesh_mat(pVec.size(),1);
	GMMDenseColMatrix tri_indices(fVec.size(),1);
	for ( int i = 0; i < (int)pVec.size(); ++i )
		mesh_mat(i,0) = std::complex<double> ( pVec[i][0] , pVec[i][1] );
	for ( int i = 0; i < (int)fVec.size(); ++i)
		tri_indices(i,0) = fVec[i];

	MatlabGMMDataExchange::SetEngineDenseMatrix( "matrix" , mesh_mat );
	MatlabGMMDataExchange::SetEngineDenseMatrix( "tri_indices" , tri_indices );

	MeshBuilder<Mesh::HalfedgeDS,Kernel> meshBuilder( &pVec, &fVec );
	this->target_mesh.clear();
	this->target_mesh.delegate( meshBuilder );
	this->target_mesh.updateAllGlobalIndices();
		


	/*auto a = target_mesh.size_of_border_edges();
	GMMDenseColMatrix pp(a, 2);
	auto bheit = target_mesh.border_edges_begin();
	for (int i = 0; i < a; ++i)
	{
		pp(i,0) = bheit->vertex()->point()[0];
		pp(i,1) = bheit->vertex()->point()[1];
		bheit++;
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("pp", pp);
	MatlabInterface::GetEngine().Eval("impoly(gca,pp)");*/


	std::cout << "Done!\n";
}



