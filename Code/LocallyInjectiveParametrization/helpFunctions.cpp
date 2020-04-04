#include "stdafx.h"
#define RESET_NUM -1000.123
#define DEBUG_MATLAB 1

CGAL::Vector_3<Kernel> normalizeVector(CGAL::Vector_3<Kernel>& v)
{
	double len = std::sqrt(v.squared_length());
	return (CGAL::Vector_3<Kernel>( v[0] / len , v[1] / len ,v[2] / len) );
}

void matchPointsIndices(Arrangement_2& arrSource, const std::vector<EPoint_2>& sourceHarmonicMapPoints, const Landmarks_pl& sourceLandMark, Mesh& sourceMesh)
{
	/*auto it = arrSource.vertices_begin();
	while ( it!= arrSource.vertices_end() )
	{
		for ( int i = 0; i < (int)sourceHarmonicMapPoints.size(); ++i )
		{
			if ( it->point() == sourceHarmonicMapPoints[i] )
			{
				//it->data()->index() = i ;
				it->data()->userIndex() = i;
				break;
			}
		}
		it++;	
	}*/

	/*for (int i = 0; i < (int)sourceHarmonicMapPoints.size(); ++i)
	{
		auto res = sourceLandMark.locate(sourceHarmonicMapPoints[i]);
		Arrangement_2::Vertex_const_iterator vertex;
		CGAL::assign(vertex, res);
		if (vertex.ptr() == nullptr)
			assert(0);
		vertex->data()->userIndex() = i;
		//CGAL::Arrangement_on_surface_2<ARRTraits_2, Dcel>::non_const_handle(vertex);
		//Arrangement_2::non_const_handle(vertex);
		
	}*/
	auto vIt = sourceMesh.vertices_begin() , vEnd = sourceMesh.vertices_end();
	int index;
	while (vIt != vEnd)
	{
		index = vIt->index();
		vIt->userIndex() = index;
		auto res = sourceLandMark.locate(sourceHarmonicMapPoints[index]);
		Arrangement_2::Vertex_const_iterator vertex;
		CGAL::assign(vertex, res);
		if (vertex.ptr() == nullptr)
			assert(0);
		//Arrangement_2::Vertex_handle v;
		(arrSource.non_const_handle(vertex))->set_data(vIt);
		
		vIt++;
	}
}

void matchEdges(Arrangement_2& arr, Mesh& source_mesh)
{
	Arrangement_2::Halfedge_handle halfEdge = arr.halfedges_begin();
	const Arrangement_2::Halfedge_handle halfEdgeEnd = arr.halfedges_end();

	//go over the halfedges of the arrangement
	while(halfEdge != halfEdgeEnd)
	{
		Mesh::Vertex_handle v1 = halfEdge->source()->data();
		Mesh::Vertex_handle v2 = halfEdge->target()->data();

		Mesh::Halfedge_around_vertex_circulator h = v2->vertex_begin();
		const Mesh::Halfedge_around_vertex_circulator hEnd = h;

		bool foundIt = false;
		do
		{
			if (v1 == h->prev()->vertex())
			{
				foundIt = true;
				break;
			}
			h++;
		} while (h != hEnd);
		assert(foundIt);
		halfEdge->set_data(h);

		halfEdge++;
	}

}

void matchFaces(Arrangement_2& arr, const std::vector<EPoint_2>& mapPoints, const std::vector<int> &fVec, Landmarks_pl& trap, Mesh& sourceMesh)
{
	/*Arrangement_2::Face_handle face = arr.faces_begin();
	const Arrangement_2::Face_handle faceEnd = arr.faces_end();
	
	//go over the faces of the arrangement
	while (face != faceEnd)
	{
		if (face->is_unbounded())
		{
			face++;
			continue;
		}
		face->set_data(face->outer_ccb()->data()->face());
		face++;
	}


	bool res1 = arr.is_valid();
	bool res2 = arr.number_of_unbounded_faces() == 1;
	bool res3 = sourceMesh.size_of_facets() == arr.number_of_faces() - arr.number_of_unbounded_faces();
	bool res4 = arr.number_of_isolated_vertices() == 0;

	bool allResults = res1 && res2 && res3 && res4;

	if (!allResults)
		assert(0);*/
	//update indices
	auto fIt = sourceMesh.facets_begin();
	int index = 0;
	for (int i = 0; i < arr.number_of_faces() - 1; ++i)
	{
		ARRNumberType x = mapPoints[fVec[index]].x() + mapPoints[fVec[index + 1]].x() + mapPoints[fVec[index + 2]].x();
		x = x / 3;
		ARRNumberType y = mapPoints[fVec[index]].y() + mapPoints[fVec[index + 1]].y() + mapPoints[fVec[index + 2]].y();
		y = y / 3;
		EPoint_2 p(x, y);

		auto res = trap.locate(p);
		Arrangement_2::Face_const_handle face;
		CGAL::assign(face, res);
		if (face.ptr() == nullptr)
			assert(0);
		if (face->is_unbounded())
			assert(0);
		//Face_handle f = Arrangement_2::non_const_handle(face);
		//Mesh::Facet_handle handle = new Mesh::Facet;
		//handle->index() = i;
		fIt->index() = i;
		arr.non_const_handle(face)->set_data(fIt);
		index = index + 3;
		fIt++;
	}

}

void loadSourceMesh( Mesh &source_mesh , std::vector<Kernel::Point_3> &pVec , std::vector<int> &fVec  )
{
	MatlabInterface::GetEngine().Eval( "nis" );
	//---------------load source mesh----------------------------
	Wavefront_obj objParser;
	const int strMaxLen = 10000;
	OPENFILENAME ofn = {0};
	TCHAR fileStr[strMaxLen] = {0};

	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFile = fileStr;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(fileStr)/sizeof(TCHAR) - 1;

	GetOpenFileName(&ofn);
	std::cout << "Loading source mesh...\n";
	std::wstring str;
	
	for ( int i = 0; i < (int)strlen(fileStr); ++i)
		str.push_back( fileStr[i] );
	
	objParser.load_file( str );
	
	int p_size = (int)objParser.m_points.size();
	GMMDenseColMatrix m_points(p_size,3);
	for ( int i = 0 ; i < p_size ; ++i )
	{
		Kernel::Point_3 point( objParser.m_points[i][0], objParser.m_points[i][1], objParser.m_points[i][2] );
		pVec.push_back( point );
		m_points(i,0) = objParser.m_points[i][0];
		m_points(i,1) = objParser.m_points[i][1];
		m_points(i,2) = objParser.m_points[i][2];
	}
	
	int f_size = (int)objParser.m_faces.size();
	GMMDenseColMatrix m_faces(f_size,3);
	for ( int i = 0 ; i < f_size ; ++i )
	{
		fVec.push_back( objParser.m_faces[i].v[0]);
		fVec.push_back( objParser.m_faces[i].v[1]);
		fVec.push_back( objParser.m_faces[i].v[2]);
		m_faces(i,0) = objParser.m_faces[i].v[0];
		m_faces(i,1) = objParser.m_faces[i].v[1];
		m_faces(i,2) = objParser.m_faces[i].v[2];
	}

	int t_size = (int)objParser.m_textureCoordinates.size();
	GMMDenseColMatrix t_points(t_size, 2);
	for (int i = 0; i < t_size; ++i)
	{
		t_points(i, 0) = objParser.m_textureCoordinates[i][0];
		t_points(i, 1) = objParser.m_textureCoordinates[i][1];
	}

	//-----------------finish loading--------------------------------
	//-----------------build mesh-------------------------------
	MeshBuilder<Mesh::HalfedgeDS,Kernel> meshBuilder( &pVec, &fVec );
	source_mesh.clear();
	source_mesh.delegate( meshBuilder );
	source_mesh.updateAllGlobalIndices();

	auto sourceVertices = source_mesh.vertices_begin();
	while (sourceVertices != source_mesh.vertices_end())
	{
		int index = sourceVertices->index();
		sourceVertices->uv() = Point_3(objParser.m_textureCoordinates[index][0], objParser.m_textureCoordinates[index][1], 0);
		sourceVertices++;
	}
	//---------------pass matlab the mesh----------------------
	MatlabGMMDataExchange::SetEngineDenseMatrix( "m_points" , m_points );
	MatlabGMMDataExchange::SetEngineDenseMatrix( "m_faces" , m_faces );
	MatlabGMMDataExchange::SetEngineDenseMatrix("t_points", t_points);
	MatlabInterface::GetEngine().Eval("m_faces=m_faces+1");
	std::cout << "Done!\n";
}

void HarmonicFlattening(Mesh &source_mesh, GMMSparseRowMatrix &u, GMMSparseRowMatrix &weightsMat, bool harmonic)
{
			//Harmonic flattening


		//***** Extracting Boundray Lengths & Boundray Vertices
		// Input:
		// std::vector<Mesh::Halfedge_iterator> borderHDS
		//
		// Output:
		// std::vector<int> verticesIndices
		// std::vector<double> partialLengths

		static bool isFirst = true;
		static std::vector<double> sourceBoundary;
		static Point_3 firstBoundaryVertex(RESET_NUM, RESET_NUM, RESET_NUM), secondBoundaryVertex(RESET_NUM, RESET_NUM, RESET_NUM);

		std::vector<Mesh::Halfedge_iterator> borderHDS;
		std::vector<int> verticesIndices;
		std::vector<double> partialLengths;
		double totalLength = 0, currLength = 0;
		Mesh::Vertex_const_handle currentV;
		Mesh::Halfedge_const_handle currH;

		source_mesh.getBorderHalfEdges(borderHDS);
		if (isFirst)
		{
			firstBoundaryVertex = borderHDS[0]->vertex()->uv();
			secondBoundaryVertex = borderHDS[1]->vertex()->uv();
		}

		//GMMDenseColMatrix dMat(borderHDS.size(),2);

		for ( int i = 0 ; i < (int)borderHDS.size() ; ++i )
		{

			//dMat(i, 0) = borderHDS[i]->vertex()->point().x();
			//dMat(i, 1) = borderHDS[i]->vertex()->point().y();

			currH = borderHDS[i];
			currentV = currH->vertex();
			currLength = currH->length();
			totalLength += currLength;
			verticesIndices.push_back(currentV->index());
			partialLengths.push_back(totalLength);
		}

		//MatlabGMMDataExchange::SetEngineDenseMatrix("dMat", dMat);

	
		int sourceMeshSize = source_mesh.size_of_vertices();
		int count = 0;
		int targetIndex = 0;
		if (!isFirst)
		{
			targetIndex = -1;
			for (int i = 0; i < (int)borderHDS.size(); ++i)
				if (borderHDS[i]->vertex()->point() == firstBoundaryVertex)
				{
					targetIndex = i;
					break;
				}
			assert(targetIndex != -1);
			if (borderHDS[(targetIndex + 2) % borderHDS.size()]->vertex()->point() != secondBoundaryVertex)
				assert(0);
		}
		//targetIndex = 0;
		//u(verticesIndices[0],0) = 1;
		u(verticesIndices[targetIndex], 0) = 1;
		int i = 1;
		int NN = borderHDS.size();

		while (i < NN)
		{
			if (isFirst)
			{
				u(verticesIndices[i], 0) = std::cos(2 * M_PI * partialLengths[i - 1] / totalLength);
				u(verticesIndices[i], 1) = std::sin(2 * M_PI * partialLengths[i - 1] / totalLength);
				sourceBoundary.push_back(std::cos(2 * M_PI * partialLengths[i - 1] / totalLength));
				sourceBoundary.push_back(std::sin(2 * M_PI * partialLengths[i - 1] / totalLength));
			}
			else
			{
				if (i  != 1)
				{
					if (i + 1 != NN)
					{
						u(verticesIndices[(targetIndex + i) % NN], 0) = 0.5*sourceBoundary[count - 2] + 0.5*sourceBoundary[count];
						u(verticesIndices[(targetIndex + i) % NN], 1) = 0.5*sourceBoundary[count - 1] + 0.5*sourceBoundary[count + 1];
					}
					else
					{
						u(verticesIndices[(targetIndex + i) % NN], 0) = 0.5*sourceBoundary[count - 2] + 0.5;
						u(verticesIndices[(targetIndex + i) % NN], 1) = 0.5*sourceBoundary[count - 1];
					}
				}
				else		// the first 'u' is (1,0)
				{
					u(verticesIndices[(targetIndex + i) % NN], 0) = 0.5*sourceBoundary[count] + 0.5;
					u(verticesIndices[(targetIndex + i) % NN], 1) = 0.5*sourceBoundary[count + 1] ;
				}

				i++;
				if (i == NN)
					break;
				u(verticesIndices[(targetIndex+i)%NN], 0) = sourceBoundary[count];
				count++;
				u(verticesIndices[(targetIndex+i)%NN], 1) = sourceBoundary[count];
				count++;
			}
			++i;
		}


	

		//************* Compute Weights **************
		// Input:
		// Mesh::Vertex_const_iterator currV; --begining of the vertices
		// Mesh::Vertex_const_iterator endV;  --ending of the vertices
		//
		// Output:
		// WeightsMat
		if (harmonic)
		{
			int numOfNegativeWeights = 0;
			Mesh::Vertex_const_iterator currV, endV;
			currV = source_mesh.vertices_begin();
			endV = source_mesh.vertices_end();

			Mesh::Halfedge_around_vertex_const_circulator hdgAroundV;

			int vi, vj;
			double sumCot = 0, cot1, cot2;
			Mesh::Point_3 p1, p2, pi, pj;
			CGAL::Vector_3<Kernel> vv1, vv2, vv3, vv4;

			//std::vector<std::vector<double>> weightsMat; //replace with GMM matrix and resize before sending to function


			// ********* remove code after replacing matrix ***********
			//weightsMat.resize( source_mesh.size_of_vertices() );
			//for ( int i = 0 ; i < (int)weightsMat.size() ; ++i )
			//	weightsMat[i].resize( mesh.size_of_vertices() );
			// ********************************************************

			for (; currV != endV; ++currV) // traverse through all mesh vertices
			{
				vi = currV->index();
				pi = currV->point();
				if (currV->is_border())
				{
					weightsMat(vi, vi) = 1;
					continue;
				}

				hdgAroundV = currV->vertex_begin();

				do // traverse through connected edges
				{
					auto oppositeV = hdgAroundV->opposite()->vertex();
					vj = oppositeV->index(); // opposite vertex index for curr edge
					pj = oppositeV->point();
					p1 = hdgAroundV->next()->vertex()->point(); // 3rd point from first triangle
					p2 = hdgAroundV->opposite()->next()->vertex()->point(); // 3rd point from second triangle

					/*vv1 = p1 - pi;
					vv2 = p1 - pj;
					vv3 = p2 - pi;
					vv4 = p2 - pj; */

					vv1 = pi - p1;
					vv2 = pj - p1;
					vv3 = pi - p2;
					vv4 = pj - p2;

					vv1 = normalizeVector(vv1); vv2 = normalizeVector(vv2);
					vv3 = normalizeVector(vv3); vv4 = normalizeVector(vv4);

					cot1 = (vv1 * vv2) / std::sqrt(CGAL::cross_product(vv1, vv2).squared_length());
					cot2 = (vv3 * vv4) / std::sqrt(CGAL::cross_product(vv3, vv4).squared_length());
					weightsMat(vi, vj) = (cot1 + cot2) / 2;
					sumCot += weightsMat(vi, vj);

					if (weightsMat(vi, vj) < 0)
						numOfNegativeWeights++;
					hdgAroundV++;
				} while (hdgAroundV != currV->vertex_begin());

				weightsMat(vi, vi) = -1 * sumCot;
				sumCot = 0;
			}
			/*if (isFirst)
				std::cout << "# of negative weights in source map: " << numOfNegativeWeights << "\n";
			else
				std::cout << "# of negative weights in target map: " << numOfNegativeWeights << "\n";*/
		}
		else	//mean value
		{

			Mesh::Vertex_const_iterator currV, endV;
			currV = source_mesh.vertices_begin();
			endV = source_mesh.vertices_end();

			Mesh::Halfedge_around_vertex_const_circulator hdgAroundV;

			int vi, vj;
			double sumCot = 0;// , t1, t2;
			Mesh::Point_3 p1, p2, pi, pj;
			CGAL::Vector_3<Kernel> vv1, vv2, vv3, vv4;

			for (; currV != endV; ++currV) // traverse through all mesh vertices
			{
				vi = currV->index();
				pi = currV->point();
				if (currV->is_border())
				{
					weightsMat(vi, vi) = 1;
					continue;
				}

				hdgAroundV = currV->vertex_begin();

				do // traverse through connected edges
				{
					/*auto oppositeV = hdgAroundV->opposite()->vertex();
					vj = oppositeV->index(); // opposite vertex index for curr edge
					pj = oppositeV->point();
					p1 = hdgAroundV->next()->vertex()->point(); // 3rd point from first triangle
					p2 = hdgAroundV->opposite()->next()->vertex()->point(); // 3rd point from second triangle

					vv1 = p1 - pi;
					vv2 = pj - pi;
					vv3 = p2 - pi;
					//vv4 = p2 - pj;

					vv1 = normalizeVector(vv1); vv2 = normalizeVector(vv2);
					vv3 = normalizeVector(vv3); //vv4 = normalizeVector(vv4);
					double y_ij = std::acos(vv2*vv1);
					double l_ij = std::acos(vv3*vv2);
					if ((y_ij > 3.15) || (l_ij > 3.15))
						assert(0);
					y_ij = y_ij / 2;
					l_ij = l_ij / 2;
					t1 = std::tan(y_ij);//  std::sqrt(CGAL::cross_product(vv1, vv2).squared_length());
					t2 = std::tan(l_ij);// (vv3 * vv4) / std::sqrt(CGAL::cross_product(vv3, vv4).squared_length());
					if ((t1 < 0) || (t2 < 0))
						assert(0);
						
					weightsMat(vi, vj) = (t1 + t2) / (2 * std::sqrt(vv2.squared_length()));
					sumCot += weightsMat(vi, vj);

					hdgAroundV++;
						*/


					auto oppositeV = hdgAroundV->opposite()->vertex();
					vj = oppositeV->index(); // opposite vertex index for curr edge
					Point_3 p0 = pi;
					Point_3 p1 = oppositeV->point();
					Point_3 p2 = hdgAroundV->next()->vertex()->point(); // 3rd point from first triangle
					Point_3 p3 = hdgAroundV->opposite()->next()->vertex()->point(); // 3rd point from second triangle

					Vector_3 u = p1 - p0;
					Vector_3 v = p2 - p0;
					Vector_3 w = p3 - p0;

					double u_length = sqrt(u.squared_length());
					double v_length = sqrt(v.squared_length());
					double w_length = sqrt(w.squared_length());

					double cross_uv = sqrt(cross_product(u, v).squared_length());
					double cross_wu = sqrt(cross_product(w, u).squared_length());

					double alphaUV_half = 0.0;
					double alphaWU_half = 0.0;

					
					alphaUV_half = (u_length*v_length - u*v) / cross_uv;
					alphaWU_half = (w_length*u_length - w*u) / cross_wu;
					
					double weight = (alphaUV_half + alphaWU_half) / u_length;
					assert(weight > 0.0); //mean value weights are supposed to be positive

					weightsMat(vi, vj) = weight;
					sumCot += weightsMat(vi, vj);

					hdgAroundV++;
				} while (hdgAroundV != currV->vertex_begin());

				weightsMat(vi, vi) = -1 * sumCot;
				sumCot = 0;
			}
		}
		isFirst = false;
}

void addPointsToTarget( Polygon_2 &poly , int numOfBorder , double avg_arc )
{
	int n = (int)poly.size();
	numOfBorder -= n;
	int i = 0;

	while ( numOfBorder > 0 )
	{
		if ( i == n )
		{
			i = 0;
			avg_arc = avg_arc / 2;
		}

		if ( std::sqrt ( poly.edge(i).squared_length() ) > avg_arc )
		{
			Point_2 p( 0.5 * poly[i][0] + 0.5 * poly[(i+1)%n][0] , 0.5 * poly[i][1] + 0.5 * poly[(i+1)%n][1] );
			poly.insert( poly.vertices_begin() + 1 + i , p );
			numOfBorder--;
			n++;
			i = i+2;
		}
		else
			i++;			
	}
	Polygon_2 tempPoly = poly;
	n = (int)tempPoly.size();
	poly.clear();
	for (int j = 0; j < n; ++j)
	{
		poly.push_back(tempPoly[j]);
		Point_2 p(0.5 * tempPoly[j][0] + 0.5 * tempPoly[(j + 1) % n][0], 0.5 * tempPoly[j][1] + 0.5 * tempPoly[(j + 1) % n][1]);
		poly.push_back(p);
	}
	/*/ debug in matlab
	GMMDenseColMatrix targetP((int)poly.size(),2);
	for ( int j = 0; j < (int)poly.size(); ++j )
	{
		targetP(j,0) = poly[j][0];
		targetP(j,1) = poly[j][1];
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix( "targetP" , targetP );*/
	
}

	void getPointsFromFace( const Arrangement_2::Face_const_handle& face, std::vector<EPoint_2>& points , std::vector<int>& indicesOrder)
	{
		Arrangement_2::Ccb_halfedge_const_circulator it = face->outer_ccb();
		points.clear();
		indicesOrder.clear();
		points.resize(3);
		indicesOrder.resize(3);
		int i = 0;
		do 
		{
			points[i] = it->target()->point();
			indicesOrder[i] = it->target()->data()->userIndex();
			i++;
			it++;
		} while ( it != face->outer_ccb() );
	}

	int orderOfIndex(const int index , std::vector<int>& indicesOrder)
	{
		if ( indicesOrder[0] == index )
			return 0;
		if ( indicesOrder[1] == index )
			return 1;
		else
			return 2;
	}

	void getPointsFromFace_Mesh( Mesh& targetMesh /*const Mesh::Face_const_handle& face*/, std::vector<EPoint_2>& points , std::vector<int>& indicesOrder )
	{/*
		Mesh::Halfedge_around_facet_const_circulator it = face->facet_begin();
		points.clear();
		points.resize(3);
		int i = 0;
		do 
		{
			int j = orderOfIndex(it->index(),indicesOrder);
			points[j] = EPoint_2(it->vertex()->point().x() , it->vertex()->point().y());
			i++;
			it++;
		} while ( it != face->facet_begin() );*/
		points.clear();
		points.resize(3);
		points[0] = EPoint_2 (targetMesh.vertex(indicesOrder[0])->point().x(),targetMesh.vertex(indicesOrder[0])->point().y());
		points[1] = EPoint_2 (targetMesh.vertex(indicesOrder[1])->point().x(),targetMesh.vertex(indicesOrder[1])->point().y());
		points[2] = EPoint_2 (targetMesh.vertex(indicesOrder[2])->point().x(),targetMesh.vertex(indicesOrder[2])->point().y());

	}


	void BuildArrangement(Arrangement_2& arr, Landmarks_pl& trap, const std::vector<EPoint_2>& vertices, const std::vector<int>& faces, /*Face_index_observer& obs,*/ Mesh &source_mesh)
	{
		static bool firstTime = true;
		if (firstTime)
		{
			firstTime = false;
			std::cout << "Building arrangement from source unit disk map...\n";
		}
		else
			std::cout << "Building arrangement from target unit disk map...\n";
/*
		std::vector<ESegment_2>    segments;
		//std::vector<EPoint_2> Evertices;
		arr.clear();

		segments.resize(faces.size());
		//Evertices.resize(vertices.size());
		//for ( int i = 0; i < (int)vertices.size(); ++i )
		//	Evertices[i] = EPoint_2 ( vertices[i][0] , vertices[i][1] );
		EPoint_2 p1, p2, p3;
		int i1,i2,i3;
		for ( int i = 0 ; i < (int)faces.size() - 2; i +=3 )
		{
			i1 = faces[i];
			i2 = faces[i+1];
			i3 = faces[i+2];
			p1 = vertices[i1];
			p2 = vertices[i2];;
			p3 = vertices[i3];;
			//segments.push_back ( ESegment_2 (p1, p2) );
			//segments.push_back ( ESegment_2 (p2, p3) );
			//segments.push_back ( ESegment_2 (p3, p1) );
			segments[i]   =	ESegment_2 (p1, p2) ;
			segments[(i+1)] =  ESegment_2 (p2, p3) ;
			segments[(i+2)] =  ESegment_2 (p3, p1);
		}
	//	insert (arr, segments.begin(), segments.end() , trap); //doesn't compile
		CGAL::insert(arr, segments.begin(), segments.end());
		//CGAL::insert_non_intersecting_curves(arr, segments.begin(), segments.end());
		trap.attach(arr);
		matchPointsIndices(arr,vertices);
		obs.matchFaces(arr,vertices,faces);
		// DEBUG //
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		//		std::cout<<"num of faces: "<<arr.number_of_faces()<<"\n";
		//		std::cout<<"num of vertices: "<<arr.number_of_vertices()<<"\n\n\n";
		//		std::cout<<"number_of_unbounded_faces: "<<arr.number_of_unbounded_faces()<<"\n\n\n\n";
		//		int Nfaces = arr.number_of_faces();
		//		int NUnBoundFaces = arr.number_of_unbounded_faces();
		//		int NumOfVertices = arr.number_of_vertices();
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/



		Mesh::Halfedge_iterator he = source_mesh.halfedges_begin();
		const Mesh::Halfedge_iterator heEnd = source_mesh.halfedges_end();

		std::list<ESegment_2> segmentsList;

		int i = 0;
		while ( he != heEnd )
		{
			if (he->is_border())
			{
				he++;
				continue;
			}
			if (he->opposite()->is_border() || he->index() < he->opposite()->index())
			{
				Point_3 p1 = he->prev()->vertex()->uv();
				Point_3 p2 = he->vertex()->uv();
				assert(p1.z() == 0.0 && p2.z() == 0.0);
				ESegment_2 segment(EPoint_2(p1.x(), p1.y()), EPoint_2(p2.x(), p2.y()));
				segmentsList.push_back(segment);
			}
			i++;
			he++;
		}

		CGAL::insert_non_intersecting_curves(arr, segmentsList.begin(), segmentsList.end());

		trap.attach(arr);
		matchPointsIndices(arr, vertices, trap, source_mesh);
		matchEdges(arr, source_mesh);
		matchFaces(arr , vertices, faces, trap, source_mesh);

		//obs.matchFaces(arr, vertices, faces, trap);

		std::cout << "Done!\n";
	}

	int findTarget(const Landmarks_pl& target, const EPoint_2& point, int& type, Arrangement_2::Face_const_handle& targetFace)
	{
		Landmarks_pl::result_type result = target.locate( point );
		Arrangement_2::Face_const_iterator face;
		Arrangement_2::Vertex_const_iterator vertex;
		Arrangement_2::Halfedge_const_iterator halfedge;

		CGAL::assign( vertex, result);
		CGAL::assign( face, result);
		CGAL::assign( halfedge, result);

		if (nullptr != face.ptr())
		{
			type = 1;
			if (face->is_unbounded())
				assert(0);
			targetFace = face;
			//return face->data()->index();
			return 1;
		}

		if( nullptr != vertex.ptr() )
		{
			type = 0;
			face = vertex->incident_halfedges()->face();
			if ( face->is_unbounded() )
				face = vertex->incident_halfedges()->twin()->face();

			targetFace = face;
			//return face->data()->index();
			return 1;
		}

		if( nullptr != halfedge.ptr() )
		{
			type = 2;
			face = halfedge->face();
			if ( face->is_unbounded() )
				face = halfedge->twin()->face();
			targetFace = face;
			//return face->data()->index();
			return 1;
		}

		
		return -1;
	}

	ARRNumberType crossProduct ( EVector_2 v1 , EVector_2 v2 )
	{
		ARRNumberType temp = v1.x() * v2.y() - v1.y() * v2.x() ;
		if ( temp < 0 ) 
			temp = temp*-1;
		return (temp);
	}

	void barycentricCord( const std::vector<EPoint_2>& points, EPoint_2 point, ARRTraits_2::Point_3 &res )
	{
		/*if ( 3 != points.size() )
			return;*/
		EPoint_2 p1 = points[0], p2 = points[1], p3 = points[2];

		EVector_2 vec1 = p1 - point;
		EVector_2 vec2 = p2 - point;
		EVector_2 vec3 = p3 - point;

		
		ARRNumberType s1 = crossProduct(vec3,vec2)/2;
		ARRNumberType s2 = crossProduct(vec3,vec1)/2;
		ARRNumberType s3 = crossProduct(vec1,vec2)/2;
		ARRNumberType s = s1+s2+s3;

		ARRTraits_2::Point_3 resp( s1/s,s2/s,s3/s);
		res = resp;

		//std::cout<<"("<<res[0]<<","<<res[1]<<","<<res[2]<<")\t= "<<res[0]+res[1]+res[2]<<"\n";

		return;

	}

	EPoint_2 reverseBarycentric ( const std::vector<EPoint_2>& points, ARRTraits_2::Point_3 bar )
	{
		EPoint_2 p1 = points[0], p2 = points[1], p3 = points[2]; 
		return ( EPoint_2( p1.x()*bar[0] + p2.x()*bar[1] + p3.x()*bar[2] ,  p1.y()*bar[0] +  p2.y()*bar[1] + p3.y()*bar[2] ) );
	}

	Point_3 reverseBarycentric(const std::vector<Point_3>& points, ARRTraits_2::Point_3 bar)
	{
		Point_3 p0 = points[0], p1 = points[1], p2 = points[2];
		Point_3 B = Mesh::Point_3(CGAL::to_double<ARRNumberType>(bar[0]), CGAL::to_double<ARRNumberType>(bar[1]), CGAL::to_double<ARRNumberType>(bar[2]));
		return (Point_3(p0.x()*B[0] + p1.x()*B[1] + p2.x()*B[2], p0.y()*B[0] + p1.y()*B[1] + p2.y()*B[2], p0.z()*B[0] + p1.z()*B[1] + p2.z()*B[2]));
	}

	void updateMeshUV( Mesh& mesh, int index , Mesh::Point_3 point )
	{
		if ( mesh.vertex(index)->uv() == Point_3(RESET_NUM, RESET_NUM, RESET_NUM) )
			mesh.vertex(index)->uv() = point;
		/*////////////
		else
		{
			if (mesh.vertex(index)->uv() != point)
			{
				std::cout << "old: (" << mesh.vertex(index)->uv().x() << "," << mesh.vertex(index)->uv().y() << ")\n";
				std::cout << "new: (" << point.x() << "," << point.y() << ")\n";

			}
		}
		/////////////*/
	}

	std::vector<int> updateUVs(Mesh& sourceMesh, Mesh& targetMesh, const Arrangement_2& source, const Arrangement_2& target, const Landmarks_pl& sourceLandMark, const Landmarks_pl& targetLandMark, /*Face_index_observer& sourceObs, Face_index_observer& targetObs,*/ std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<int>& fVec, std::vector<Point_3>& uvVector)
	{
		/*GMMDenseColMatrix newFvec(source.number_of_faces()-1,3);
		int count=0;
		EPoint_2 tempP;
		std::vector<int> indicesOrder;

		std::vector<int> negativeOrientationTriangles;
		std::vector<EPoint_2> potentialUV;
		for ( Arrangement_2::Face_const_iterator face = source.faces_begin(); face != source.faces_end(); ++face )
		{
			if ( face->is_unbounded() )
				continue;
			Arrangement_2::Face_const_handle targetFace;
			//Mesh::Face_const_handle targetMeshFace;
			std::vector<EPoint_2> points;
			std::vector<EPoint_2> targetMeshPoints;
			ARRTraits_2::Point_3 barPoint;
			Arrangement_2::Ccb_halfedge_const_circulator halfedge = face->outer_ccb();
			std::vector<int> faceVerticesIndexes;

			faceVerticesIndexes.clear();
			potentialUV.clear();

			int type, index;
			do 
			{
				int verIndex = halfedge->target()->data()->index();
				tempP = halfedge->target()->point();
				faceVerticesIndexes.push_back(verIndex); // vertex index
				if ( sourceMesh.vertex( verIndex )->uv() != Point_3 ( RESET_NUM , RESET_NUM , RESET_NUM ) )
				{
					potentialUV.push_back ( EPoint_2 ( sourceMesh.vertex( verIndex )->uv()[0] , sourceMesh.vertex( verIndex )->uv()[1] ) ); 
					
				}
				else
				{
				index = findTarget( targetLandMark, tempP , type ); // face index
				targetFace = targetObs.getFace(index);
				getPointsFromFace( targetFace, points , indicesOrder);
				barycentricCord( points, tempP , barPoint );
				//targetMeshFace = targetMesh.face(index);
				getPointsFromFace_Mesh( targetMesh , targetMeshPoints , indicesOrder);
				potentialUV.push_back(reverseBarycentric( targetMeshPoints, barPoint ));
				}
				halfedge++;

			} while ( halfedge != face->outer_ccb() ); // three potential UVs
			
			if ( potentialUV.size() != 3 )
				assert(0); // more or less than 3 vertices in the face

			ARRKernel::Triangle_2 tri( potentialUV[2], potentialUV[1], potentialUV[0]);

			//if ( tri.orientation() > 0 )
			//{
				updateMeshUV( sourceMesh, faceVerticesIndexes[0], Mesh::Point_3( CGAL::to_double<ARRNumberType>(potentialUV[0].x() ) , CGAL::to_double<ARRNumberType>(potentialUV[0].y() ) , 0 ) );
				updateMeshUV( sourceMesh, faceVerticesIndexes[1], Mesh::Point_3( CGAL::to_double<ARRNumberType>(potentialUV[1].x() ) , CGAL::to_double<ARRNumberType>(potentialUV[1].y() ) , 0 ) );
				updateMeshUV( sourceMesh, faceVerticesIndexes[2], Mesh::Point_3( CGAL::to_double<ARRNumberType>(potentialUV[2].x() ) , CGAL::to_double<ARRNumberType>(potentialUV[2].y() ) , 0 ) );
			//}
			//else
			if ( tri.orientation() != 1 )
			{
				negativeOrientationTriangles.push_back(faceVerticesIndexes[2]);
				negativeOrientationTriangles.push_back(faceVerticesIndexes[1]);
				negativeOrientationTriangles.push_back(faceVerticesIndexes[0]);
			}
			else
			{
				newFvec(count, 0) = faceVerticesIndexes[0];
				newFvec(count, 1) = faceVerticesIndexes[1];
				newFvec(count, 2) = faceVerticesIndexes[2];
				count++;
			}
		}
		MatlabGMMDataExchange::SetEngineDenseMatrix( "newFvec" , newFvec );
		return negativeOrientationTriangles;*/
		int numOfTri = (int)fVec.size() / 3 , index,type;
		EPoint_2 potentialUV[3];
		std::vector<int> indicesOrder ;
		std::vector<EPoint_2> points;
		std::vector<EPoint_2> targetMeshPoints;
		ARRTraits_2::Point_3 barPoint;
		Arrangement_2::Face_const_handle targetFace;
		std::vector<int> negativeOrientationTriangles;

		for (int i = 0; i < numOfTri; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (uvVector[fVec[3*i + j]].z() != RESET_NUM)
				{
					potentialUV[j] = EPoint_2(uvVector[fVec[3*i + j]].x(), uvVector[fVec[3*i + j]].y());	
				}
				else
				{
					EPoint_2 tempP = sourceHarmonicMapPoints[fVec[3 * i + j]];
					indicesOrder.clear();
					points.clear();
					targetMeshPoints.clear();

					index = findTarget(targetLandMark, tempP, type, targetFace); // face index
					//targetFace = targetMesh.face(index);//targetObs.getFace(index);
					getPointsFromFace(targetFace, points, indicesOrder);
	
					barycentricCord(points, tempP, barPoint);
					getPointsFromFace_Mesh(targetMesh, targetMeshPoints, indicesOrder);
					potentialUV[j] = (reverseBarycentric(targetMeshPoints, barPoint));
					uvVector[fVec[3*i + j]] = Mesh::Point_3(CGAL::to_double<ARRNumberType>(potentialUV[j].x()), CGAL::to_double<ARRNumberType>(potentialUV[j].y()), 0);
				}
			}

			ARRKernel::Triangle_2 tri(potentialUV[0], potentialUV[1], potentialUV[2]);
			if (tri.orientation() != 1)
				negativeOrientationTriangles.push_back(i);

			
		}
		return (negativeOrientationTriangles);
	}

	void setBoundaryUV(Mesh &source_mesh, Mesh &target_mesh, std::vector<Point_3>& uvVector)
	{
		std::vector<Mesh::Halfedge_iterator> borderSourceHDS, borderTargetHDS;
		source_mesh.getBorderHalfEdges(borderSourceHDS);
		target_mesh.getBorderHalfEdges(borderTargetHDS);

		auto firstBoundaryVertex = borderSourceHDS[0]->vertex()->uv();
		auto secondBoundaryVertex = borderSourceHDS[1]->vertex()->uv();

		int targetIndex = -1;
		for (int i = 0; i < (int)borderTargetHDS.size(); ++i)
		{
			if (borderTargetHDS[i]->vertex()->uv() == firstBoundaryVertex)
			{
				targetIndex = i;
				break;
			}
		}
		assert(targetIndex != -1);
		if (secondBoundaryVertex != borderTargetHDS[(targetIndex + 2) % borderTargetHDS.size()]->vertex()->uv())
			assert(0);
	
		int jj = 0;
		for (auto i = source_mesh.vertices_begin(); i != source_mesh.vertices_end(); ++i)	//reset all uv's
		{
			i->uv() = Point_3(RESET_NUM, RESET_NUM, RESET_NUM);
			uvVector[jj] = i->uv();
			jj++;
		}

		int N = borderTargetHDS.size();

		for ( int i = 0 ; i < (int)borderSourceHDS.size() ; ++i )
		{
			//std::cout<<"borderS index: "<<borderSourceHDS[i]->vertex()->index()<<"\t"<<"borderT index: "<<borderTargetHDS[i]->vertex()->index()<<"\n";
			//std::cout << "borderS: (" << borderSourceHDS[i]->vertex()->point().x() << "," << borderSourceHDS[i]->vertex()->point().y() << ")\t" << "borderT: (" << borderTargetHDS[i]->vertex()->point().x() << "," << borderTargetHDS[i]->vertex()->point().y()<< ")\n";
			source_mesh.vertex(borderSourceHDS[i]->vertex()->index())->uv() = borderTargetHDS[(targetIndex+2*i)%N]->vertex()->point();
			uvVector[borderSourceHDS[i]->vertex()->index()] = Point_3(borderTargetHDS[(targetIndex+2*i)%N]->vertex()->point().x(), borderTargetHDS[(targetIndex+2*i)%N]->vertex()->point().y(), 1);	// the z=1 is to mark that this is boundary vertex
	
			//std::cout << borderSourceHDS[i]->vertex()->index() << ":\t(" << uvVector[borderSourceHDS[i]->vertex()->index()].x() << "," << uvVector[borderSourceHDS[i]->vertex()->index()].y() << ")\n";
		}
	}

	int refine(std::vector<int>& neg, Mesh& sourceMesh, Mesh& targetMesh, Arrangement_2& source, Arrangement_2& target, const Landmarks_pl& sourceLandMark, const Landmarks_pl& targetLandMark, /*Face_index_observer& targetObs,*/ std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<Kernel::Point_3> &pVec, std::vector<int>& fVec, std::vector<Point_3>& uvVector)
	{
		int size = neg.size();

		if (size == 0)
		{
			std::cout << "No need to refine.\n";
			return 0;
		}
		else
			std::cout << "# of negative triangles: " << size << "\nrefine...\n";

		// --- find all the neighbors of the triangles in 'neg' ---// we found them inside the loop
		//findNegativeNeighbors(neg, sourceHarmonicMapPoints, fVec, source, sourceLandMark);
		int startSize = pVec.size();
		std::vector<Point_3> tempUV = uvVector;
		PointMap pMap;
		int updatedNumOfPoints = uvVector.size();

		std::vector<bool> isRefined(fVec.size() / 3), inTheList(fVec.size() / 3), isTriangulated(fVec.size() / 3);

		// ----init the boolean vectors-------------
		for (int i = 0; i < isRefined.size(); ++i)
		{
			isRefined[i] = false;
			inTheList[i] = false;
			isTriangulated[i] = false;
		}

		for (int i = 0; i < size; ++i)
			inTheList[neg[i]] = true;
		//-----------------------------------------

		for (int i = 0; i < neg.size(); ++i)
		{
			if (i < size)
			{
				findTriangleNeighbors(neg[i], fVec, sourceHarmonicMapPoints, sourceLandMark, neg, inTheList);
				refineTriangle(neg[i], fVec, sourceHarmonicMapPoints, pMap, uvVector, targetMesh, target, targetLandMark, sourceLandMark, /*targetObs,*/ pVec, updatedNumOfPoints);
				isRefined[ neg[i] ] = true;
				continue;
			}

			if (isRefined[neg[i]])
				continue;

			if (checkIfSimple(neg[i], fVec, pMap, uvVector))
			{
				inTheList[neg[i]] = false;
				continue;
			}
			else
			{
				findTriangleNeighbors(neg[i], fVec, sourceHarmonicMapPoints, sourceLandMark, neg, inTheList);
				refineTriangle(neg[i], fVec, sourceHarmonicMapPoints, pMap, uvVector, targetMesh, target, targetLandMark, sourceLandMark, /*targetObs,*/ pVec, updatedNumOfPoints);
				isRefined[neg[i]] = true;
			}
			/*if ( !triangulateNeighbor(neg[i], fVec,  sourceHarmonicMapPoints,  pMap, uvVector) )
			{
				findTriangleNeighbors(neg[i], fVec, sourceHarmonicMapPoints, sourceLandMark, neg);
				refineTriangle(neg[i], fVec, sourceHarmonicMapPoints, pMap, uvVector, targetMesh, target, targetLandMark, sourceLandMark, targetObs, pVec, updatedNumOfPoints);
			}*/
		}

		std::cout << "Done!\nSimplify - trying to reduce number of new points...\n";

		// need to add here simplify function
		CGAL::Timer simplifyTimer;
		simplifyTimer.start();
		simplify(pMap, fVec, uvVector);
		simplifyTimer.stop();
		
		std::cout << "Done!\n" << "Total time to simplify : " << simplifyTimer.time() << " seconds\nTriangulating new polygons...\n";

		int updateSize = tempUV.size();
		std::map<Point_3, int> newUVtoIndicesMap;
		std::set<Pair> visitedPairs;
		auto mapIt = pMap.edgeToFace.begin();
		auto mapEnd = pMap.edgeToFace.end();
		while (mapIt != mapEnd)
		{
			Pair e = mapIt->first;
			if (visitedPairs.count(e) != 0)
			{
				mapIt++;
				continue;
			}
			visitedPairs.insert(e);
			visitedPairs.insert(e.reverse());
			std::vector<Point_3> &vecOfPointsP1 = pMap.getEdgeUVs(e);
			std::vector<Point_3> & meshVecP1 = pMap.getNewMeshPointsVec(e);
			for (int i = 1; i < vecOfPointsP1.size() - 1; ++i)
			{
				newUVtoIndicesMap[vecOfPointsP1[i]] = updateSize;
				tempUV.push_back(vecOfPointsP1[i]);
				pVec.push_back(meshVecP1[i]);
				updateSize++;
			}
			mapIt++;
		}

		bool result;
		int len = neg.size();
		for (int i = 0; i < len; ++i)
		{
			if (isTriangulated[neg[i]])
				continue;
			result = triangulateNeighbor(neg[i], fVec, sourceHarmonicMapPoints, pMap, tempUV, newUVtoIndicesMap);
			assert(result);
			isTriangulated[neg[i]] = true;
		}
		
		//////////////////////////////////////////////////
		// -- need to triangulate all neighbors faces --//
		//////////////////////////////////////////////////

		//triangulateNeighbors(neighTri, fVec, sourceHarmonicMapPoints,pMap, uvVector);
		int numOfNewPoints = pVec.size() - startSize;
		std::cout << "Done!\n# of new points: " << numOfNewPoints << "\n";
		uvVector = tempUV;
		return (numOfNewPoints);
	}

	void findIntersection( Arrangement_2& target, const Landmarks_pl& targetLandMark, ESegment_2 seg , std::vector<EPoint_2>& intersectionPoints)
	{
		std::list<CGAL::Object> intersectionList;
		Arrangement_2::Halfedge_handle halfedge;
		Arrangement_2::Vertex_handle ver;
		Arrangement_2::Face_handle face;
		EPoint_2 p1, p2, pp1(RESET_NUM, RESET_NUM), pp2(RESET_NUM, RESET_NUM);

		CGAL::zone(target, seg, std::back_inserter(intersectionList), targetLandMark);

		//std::vector<Arrangement_2::Halfedge_handle> temp;
		//int k = 0;
		//std::cout << "from:\t(" << seg.source().x() << "," << seg.source().y() << ")\n";
		if (intersectionList.size() < 2)
			return;
		auto it = intersectionList.begin();
		auto itEnd = intersectionList.end();
		it++;
		itEnd--;

		for (; it != itEnd; it++)
		//for (int i = 1; i < (int)intersectionList.size()-1; ++i)
		{
			if (assign(face, *it /*intersectionList[i]*/))
			{
				if (face->is_unbounded())
					assert(0);
				continue;
			}

			if (assign(halfedge, *it /*intersectionList[i]*/))
			{
				if (halfedge->next()->next()->next()->data()->index() != halfedge->data()->index())
				{
					halfedge = halfedge->twin();
					//continue;
				}

				p1 = halfedge->source()->point();
				p2 = halfedge->target()->point();
				if ( ( (p1 == pp1) && (p2 == pp2) ) || ( (p1 == pp2) && (p2 == pp1) ) )
					continue;
				pp1 = p1;
				pp2 = p2;

				//std::cout << "p1 = (" << p1.x() << "," << p1.y() << ")\tp2 = (" << p2.x() << "," << p2.y() << ")\n";

				std::vector<EPoint_2> localIntersectionPoints;
				std::vector<ESegment_2> pp(2);
				pp[0] = seg;
				pp[1] = ESegment_2(p1, p2);				
				CGAL::compute_intersection_points(pp.begin(), pp.end(), std::back_inserter(localIntersectionPoints));

				if (localIntersectionPoints.size() == 1)
					intersectionPoints.push_back(localIntersectionPoints[0]);

				continue;
				//std::cout << k << ":\t(" << intersectionPoints[k].x() << "," << intersectionPoints[k].y() << ")\n";
				//k++;

			}

			if (assign(ver, *it /*intersectionList[i]*/))
			{
				intersectionPoints.push_back(ver->point());
				//std::cout << k << ":\t(" << intersectionPoints[k].x() << "," << intersectionPoints[k].y() << ")\n";
				//k++;
				//std::cout << "I'm vertex!\n";
			}
			
		}


		//std::cout << "to:\t(" << seg.target().x() << "," << seg.target().y() << ")\n\n\n";
	}

	void inverseList(std::vector<EPoint_2>& intersectionPoints)
	{
		std::reverse(intersectionPoints.begin(), intersectionPoints.end());
		/*std::vector<EPoint_2> revList;
		for (int i = (int)intersectionPoints.size() - 1; i >= 0; --i)
			revList.push_back(intersectionPoints[i]);
		for (int i = 0; i < (int)intersectionPoints.size(); ++i)
			intersectionPoints[i] = revList[i];*/
		
	}

	Point_3 calcNewUV(EPoint_2 tempP, Mesh& targetMesh, const Landmarks_pl& targetLandMark, /*Face_index_observer& targetObs,*/ std::vector<Point_3>& uvVector)
	{
		std::vector<int> indicesOrder;
		std::vector<EPoint_2> points;
		std::vector<EPoint_2> targetMeshPoints;
		ARRTraits_2::Point_3 barPoint;
		Arrangement_2::Face_const_handle targetFace;

		int type;
		int index = findTarget(targetLandMark, tempP, type, targetFace); // face index
		//targetFace = targetObs.getFace(index);
		getPointsFromFace(targetFace, points, indicesOrder);

		barycentricCord(points, tempP, barPoint);
		getPointsFromFace_Mesh(targetMesh, targetMeshPoints, indicesOrder);
		auto res = reverseBarycentric(targetMeshPoints, barPoint);
		uvVector.push_back(Mesh::Point_3(CGAL::to_double<ARRNumberType>(res.x()), CGAL::to_double<ARRNumberType>(res.y()), 0));
		return ( uvVector[uvVector.size()-1] );
	}

	Point_3 calcNewMeshPoint(EPoint_2 tempP, int faceIndex, std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<Kernel::Point_3> &pVec, std::vector<int>& fVec)
	{
		std::vector<Point_3> meshFacePoints;
		std::vector<int> indicesOrder;
		std::vector<EPoint_2> points;
		ARRTraits_2::Point_3 barPoint;
	
		points.resize(3);
		meshFacePoints.resize(3);
		for (int i = 0; i < 3; ++i)
		{
			points[i] = sourceHarmonicMapPoints[fVec[3 * faceIndex + i]];
			meshFacePoints[i] = pVec[fVec[3 * faceIndex + i]];
		}

		barycentricCord(points, tempP, barPoint);
		Point_3 p = reverseBarycentric(meshFacePoints, barPoint);
		
		return(p);
		//pVec.push_back(p);
	}

	void extractIndicesFromPair(Pair& p, PointMap& pMap, std::vector<int>& polygonIndices, std::map<Point_3, int>& newUVtoIndicesMap)
	{
		polygonIndices.push_back(p.x);
		auto vec = pMap.getEdgeUVs(p);
		for (int i = 1; i < (int)vec.size()-1; ++i)
			polygonIndices.push_back( newUVtoIndicesMap[vec[i]] );
	}

	void extractIndicesFromPair(Pair& e1, Pair& e2, Pair& e3, PointMap& pMap, std::vector<int>& polygonIndices, std::map<Point_3, int>& newUVtoIndicesMap)
	{
		extractIndicesFromPair(e1, pMap, polygonIndices, newUVtoIndicesMap);
		extractIndicesFromPair(e2, pMap, polygonIndices, newUVtoIndicesMap);
		extractIndicesFromPair(e3, pMap, polygonIndices, newUVtoIndicesMap);
	}

	void triangulateNewPolygon(std::vector<int>& polygonIndices, std::vector<Point_3>& uvVector, std::vector<int>& fVec)
	{
		std::map <int, int> mapToOriginalIndices;
		/*std::vector<std::complex<double> > polygonPoints;
		std::vector<unsigned int> triangleIndices;*/
		Polygon_2 localPoly;
		int N = (int)polygonIndices.size();

		//GMMDenseColMatrix testPoly(N, 2);

		for (int i = 0; i < N; ++i)
		{
			mapToOriginalIndices[i] = polygonIndices[i];
			//polygonPoints.push_back(std::complex<double>(uvVector[polygonIndices[i]].x(), uvVector[polygonIndices[i]].y()));
			localPoly.push_back(Point_2(uvVector[polygonIndices[i]].x(), uvVector[polygonIndices[i]].y()));

			//testPoly(i, 0) = uvVector[polygonIndices[i]].x();
			//testPoly(i, 1) = uvVector[polygonIndices[i]].y();

		}

					/*#ifdef DEBUG_MATLAB	////GMMDenseColMatrix testPoly(N, 2);
						MatlabGMMDataExchange::SetEngineDenseMatrix("testPoly", testPoly);
						MatlabInterface::GetEngine().Eval("figure");
						MatlabInterface::GetEngine().Eval("impoly(gca,testPoly)");
					#endif*/
		Shor localShor;
		localShor.load_polygon(localPoly, localPoly);
		localShor.load_rArray(NULL);
		localShor.play();
		localShor.build_triangulation();

		if (localShor.fVec.size() == 0)
			assert(0);

		for (int i = 0; i < (int)localShor.fVec.size(); ++i)
			fVec.push_back(mapToOriginalIndices[localShor.fVec[i]]);
		/*triangulatePolygonWithoutAddingVertices( polygonPoints , triangleIndices );
		for (int i = 0; i < (int)triangleIndices.size(); ++i)
			fVec.push_back(mapToOriginalIndices[triangleIndices[i]]);*/
	
	}



	void findNegativeNeighbors(std::vector<int>& neg, std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<int>& fVec, Arrangement_2& source, const Landmarks_pl& sourceLandMark)
	{
		// i don't use this function anymore!
		int N = neg.size();
		for (int i = 0; i < N; ++i)
			std::cout << i;
			//findTriangleNeighbors(neg[i], fVec, sourceHarmonicMapPoints, sourceLandMark, neg);
	
	}

	void findTriangleNeighbors(int triIndex, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, const Landmarks_pl& sourceLandMark, std::vector<int>& neighTri, std::vector<bool>& inTheList)
	{
		for (int j = 0; j < 3; ++j)
		{
			int v1 = 3 * triIndex + j;
			int v2 = 3 * triIndex + ((j + 1) % 3);
			EPoint_2 p1 = sourceHarmonicMapPoints[fVec[v1]];
			EPoint_2 p2 = sourceHarmonicMapPoints[fVec[v2]];
			ARRNumberType x = p1.x() + p2.x();
			ARRNumberType y = p1.y() + p2.y();
			EPoint_2 p(x / 2, y / 2);
			auto res = sourceLandMark.locate(p);
			Arrangement_2::Halfedge_const_iterator halfedge;
			if (!CGAL::assign(halfedge, res))
				assert(0);
			if ((halfedge->face()->is_unbounded()) || (halfedge->twin()->face()->is_unbounded()))
				continue;
			if (halfedge->face()->data()->index() != triIndex)
			{
				if (inTheList[halfedge->face()->data()->index()] == false)
				{
					inTheList[halfedge->face()->data()->index()] = true;
					neighTri.push_back(halfedge->face()->data()->index());
				}
			}
			else
			{
				if (inTheList[halfedge->twin()->face()->data()->index()] == false)
				{
					inTheList[halfedge->twin()->face()->data()->index()] = true;
					neighTri.push_back(halfedge->twin()->face()->data()->index());
				}
			}
		}
	}

	void refineTriangle(int triIndex, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, PointMap& pMap, std::vector<Point_3>& uvVector, Mesh& targetMesh, Arrangement_2& target, const Landmarks_pl& targetLandMark, const Landmarks_pl& sourceLandMark, /*Face_index_observer& targetObs,*/ std::vector<Kernel::Point_3> &pVec, int& updatedNumOfPoints)
	{
		for (int j = 0; j < 3; ++j)	//for each edge in the triangle find the intersections
		{
			int v1 = 3 * triIndex + j;
			int v2 = 3 * triIndex + ((j + 1) % 3);
			std::vector<EPoint_2> intersectionPoints;
			Pair e(fVec[v1], fVec[v2]);


			if (((pMap[e]).size() != 0) || ((pMap[e.reverse()]).size() != 0))	//check if we already find the intersections
				continue;

			ESegment_2 seg = ESegment_2(sourceHarmonicMapPoints[fVec[v1]], sourceHarmonicMapPoints[fVec[v2]]);
			if (checkIfBoundaryEdge(seg, sourceLandMark))
				continue;

			findIntersection(target, targetLandMark, seg, intersectionPoints);
			if (intersectionPoints.size() == 0)
				continue;

			std::vector<Point_3> edgeUV, newMeshVec;
			edgeUV.push_back(uvVector[fVec[v1]]);
			newMeshVec.push_back(pVec[fVec[v1]]);

			if (intersectionPoints.size()>1)
			{
				ESegment_2 s1(sourceHarmonicMapPoints[fVec[v1]], intersectionPoints[0]);
				ESegment_2 s2(sourceHarmonicMapPoints[fVec[v1]], intersectionPoints[1]);
				if (s1.line().to_vector().squared_length() > s2.line().to_vector().squared_length())
					inverseList(intersectionPoints);
			}

			pMap.insert(e, intersectionPoints);	//update the map

			for (int k = 0; k < (int)intersectionPoints.size(); ++k)	// for each new point find the uv cord , and the point to refine in the original source mesh
			{
				pMap.updateNewPoint(updatedNumOfPoints, intersectionPoints[k]);
				updatedNumOfPoints++;

				edgeUV.push_back ( calcNewUV(intersectionPoints[k], targetMesh, targetLandMark, /*targetObs,*/ uvVector) );
				newMeshVec.push_back( calcNewMeshPoint(intersectionPoints[k], triIndex, sourceHarmonicMapPoints, pVec, fVec) );
				//pVec.push_back(newMeshVec[newMeshVec.size() - 1]);	//need to remove it after the simplify function work!
			}
			edgeUV.push_back(uvVector[fVec[v2]]);
			newMeshVec.push_back(pVec[fVec[v2]]);
			pMap.insertUV(e, edgeUV);
			pMap.insertNewMeshPoint(e, newMeshVec);

			//now we found the face index of the edge and the opposite edge
			pMap.edgeToFace[e] = triIndex;

			EPoint_2 p1 = sourceHarmonicMapPoints[fVec[v1]];
			EPoint_2 p2 = sourceHarmonicMapPoints[fVec[v2]];
			ARRNumberType x = p1.x() + p2.x();
			ARRNumberType y = p1.y() + p2.y();
			EPoint_2 p(x / 2, y / 2);
			auto res = sourceLandMark.locate(p);
			Arrangement_2::Halfedge_const_iterator halfedge;
			if (!CGAL::assign(halfedge, res))
				assert(0);
			if (halfedge->face()->data()->index() != triIndex)
				pMap.edgeToFace[e.reverse()] = halfedge->face()->data()->index();
			else
				pMap.edgeToFace[e.reverse()] = halfedge->twin()->face()->data()->index();
			

		}
		// now we need to triangulate the polygon with the new points

		/*
		Pair e1(fVec[3 * triIndex + 0], fVec[3 * triIndex + 1]);
		Pair e2(fVec[3 * triIndex + 1], fVec[3 * triIndex + 2]);
		Pair e3(fVec[3 * triIndex + 2], fVec[3 * triIndex + 0]);
		std::vector<int> polygonIndices;
		extractIndicesFromPair(e1, e2, e3, pMap, polygonIndices);
		if (polygonIndices.size() == 3)	//no points added
			assert(0);
		triangulateNewPolygon(polygonIndices, uvVector, fVec);

		fVec[3 * triIndex + 0] = -1;
		fVec[3 * triIndex + 1] = -1;
		fVec[3 * triIndex + 2] = -1;*/

	}


	void triangulateNeighbors(std::vector<int>& neighTri, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, PointMap& pMap, std::vector<Point_3>& uvVector)
	{
		int N = neighTri.size();
		for (int i = 0; i < N; ++i)
		{
			Pair e1(fVec[3 * neighTri[i] + 0], fVec[3 * neighTri[i] + 1]);
			Pair e2(fVec[3 * neighTri[i] + 1], fVec[3 * neighTri[i] + 2]);
			Pair e3(fVec[3 * neighTri[i] + 2], fVec[3 * neighTri[i] + 0]);

			if (e1.x == -1)
				continue;

			std::vector<int> polygonIndices;
			extractIndicesFromPair(e1, e2, e3, pMap, polygonIndices);
			
			if (polygonIndices.size() == 3)
				continue;

			triangulateNewPolygon(polygonIndices, sourceHarmonicMapPoints, fVec, pMap, uvVector);

			fVec[3 * neighTri[i]] = -1;
			fVec[3 * neighTri[i] + 1] = -1;
			fVec[3 * neighTri[i] + 2] = -1;

		}
	
	}

	bool triangulateNewPolygon(std::vector<int>& polygonIndices, std::vector<EPoint_2>& sourceHarmonicMapPoints, std::vector<int>& fVec, PointMap& pMap, std::vector<Point_3>& uvVector)
	{

		std::map <int, int> mapToOriginalIndices;
		Polygon_2 localPoly;
		int N = (int)polygonIndices.size() , len = sourceHarmonicMapPoints.size();

		for (int i = 0; i < N; ++i)
		{
			mapToOriginalIndices[i] = polygonIndices[i];
			localPoly.push_back(Point_2(uvVector[polygonIndices[i]].x(), uvVector[polygonIndices[i]].y()));
		}

#ifdef DEBUG_MATLAB1
		GMMDenseColMatrix testPoly(N, 2);
		for (int i = 0; i < N; ++i)
		{
			testPoly(i, 0) = localPoly[i].x();
			testPoly(i, 1) = localPoly[i].y();
		}
		MatlabGMMDataExchange::SetEngineDenseMatrix("testPoly", testPoly);
		MatlabInterface::GetEngine().Eval("figure");
		MatlabInterface::GetEngine().Eval("impoly(gca,testPoly)");
#endif

		Shor localShor;
		localShor.load_polygon(localPoly, localPoly);
		localShor.load_rArray(NULL);
		localShor.play();
		localShor.build_triangulation();

		if (localShor.fVec.size() == 0)
			return (false);

		for (int i = 0; i < (int)localShor.fVec.size(); ++i)
			fVec.push_back(mapToOriginalIndices[localShor.fVec[i]]);
		
		return(true);
	}
	


	bool checkIfBoundaryEdge(ESegment_2& seg, const Landmarks_pl& sourceLandMark)
	{
		EPoint_2 p1 = seg.source();
		EPoint_2 p2 = seg.target();
		ARRNumberType x = p1.x() + p2.x();
		ARRNumberType y = p1.y() + p2.y();
		EPoint_2 p(x / 2, y / 2);
		auto res = sourceLandMark.locate(p);
		Arrangement_2::Halfedge_const_iterator halfedge;
		if (!CGAL::assign(halfedge, res))
			assert(0);

		if ((halfedge->face()->is_unbounded()) || (halfedge->twin()->face()->is_unbounded()))
			return (true);
		else
			return (false);
	}

	


	bool triangulateNeighbor(int triIndex, std::vector<int>& fVec, std::vector<EPoint_2>& sourceHarmonicMapPoints, PointMap& pMap, std::vector<Point_3>& uvVector, std::map<Point_3, int>& newUVtoIndicesMap)
	{
		Pair e1(fVec[3 * triIndex + 0], fVec[3 * triIndex + 1]);
		Pair e2(fVec[3 * triIndex + 1], fVec[3 * triIndex + 2]);
		Pair e3(fVec[3 * triIndex + 2], fVec[3 * triIndex + 0]);

		if (e1.x == -1)
			return (true);	// already triangulated

		std::vector<int> polygonIndices;
		extractIndicesFromPair(e1, e2, e3, pMap, polygonIndices, newUVtoIndicesMap);

		if (polygonIndices.size() == 3)
			return (true); // no need to triangulate

		if (!triangulateNewPolygon(polygonIndices, sourceHarmonicMapPoints, fVec, pMap, uvVector))
			return (false);

		fVec[3 * triIndex + 0] = -1;
		fVec[3 * triIndex + 1] = -1;
		fVec[3 * triIndex + 2] = -1;
	
		return (true);
	}


	bool checkIfSimple(int triIndex, std::vector<int>& fVec,  PointMap& pMap, std::vector<Point_3>& uvVector)
	{
		Pair e1(fVec[3 * triIndex + 0], fVec[3 * triIndex + 1]);
		Pair e2(fVec[3 * triIndex + 1], fVec[3 * triIndex + 2]);
		Pair e3(fVec[3 * triIndex + 2], fVec[3 * triIndex + 0]);

		if (e1.x == -1)
			assert(0);	//if the triangle is already triangulated we dont suppose to be here!

		std::vector<int> polygonIndices;
		extractIndicesFromPair(e1, e2, e3, pMap, polygonIndices);
		
		Polygon_2 localPoly;
		int N = (int)polygonIndices.size();

		for (int i = 0; i < N; ++i)
			localPoly.push_back(Point_2(uvVector[polygonIndices[i]].x(), uvVector[polygonIndices[i]].y()));
		

		if (localPoly.is_simple())
			return (true);
		else
			return (false);
	}


	//zone(Arrangement_2& arr, const ESegment_2& seg, std::list<CGAL::Object>& intersectionList, Landmarks_pl& pointLocatorStructure)
	/*
	std::vector<CGAL::Object> intersectionList;
	ESegment_2 seg(EPoint_2(-0.2, -0.66), EPoint_2(-0.4, -0.5));
	Arrangement_2::Halfedge_handle halfedge;
	CGAL::zone(target, seg, std::back_inserter(intersectionList), targetLandMark);
	std::vector<EPoint_2> intersectionPoints;
	for (int i = 0; i < (int)intersectionList.size(); ++i)
	{
	assign(halfedge, intersectionList[i]);
	if (nullptr != halfedge.ptr())
	{
	EPoint_2 p1 = halfedge->source()->point();
	EPoint_2 p2 = halfedge->target()->point();
	int te = halfedge->face()->data()->index();
	int te2 = halfedge->twin()->face()->data()->index();
	std::vector<ESegment_2> pp;

	pp.push_back(ESegment_2 (p1, p2));
	pp.push_back(seg);

	CGAL::compute_intersection_points(pp.begin(), pp.end(), std::back_inserter(intersectionPoints));
	for (int j = 0; j < (int)intersectionPoints.size(); ++j)
	std::cout << i << ":\t" << "(" << intersectionPoints[j].x() << "," << intersectionPoints[j].y() << ")\n";// ---> (" << p2.x() << ", " << p2.y() << ")\n";
	}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	std::vector<int> indicesOrder;
					std::vector<EPoint_2> points;
					std::vector<EPoint_2> targetMeshPoints;
					ARRTraits_2::Point_3 barPoint;
					Arrangement_2::Face_const_handle targetFace;

					EPoint_2 tempP = intersectionPoints[i];
					int type;
					int index = findTarget(targetLandMark, tempP, type); // face index
					targetFace = targetObs.getFace(index);
					getPointsFromFace(targetFace, points, indicesOrder);

					barycentricCord(points, tempP, barPoint);
					getPointsFromFace_Mesh(targetMesh, targetMeshPoints, indicesOrder);
					auto res = reverseBarycentric(targetMeshPoints, barPoint);
					uvVector.push_back ( Mesh::Point_3(CGAL::to_double<ARRNumberType>(res.x()), CGAL::to_double<ARRNumberType>(res.y()), 0) );
	
	
	
	
	
	
	
	
	
	
	*/
	void setOrderOfEdges(Pair *p, int face, std::vector<int>& fVec)
	{
		Pair e1(fVec[3 * face + 0], fVec[3 * face + 1]);
		Pair e2(fVec[3 * face + 1], fVec[3 * face + 2]);
		Pair e3(fVec[3 * face + 2], fVec[3 * face + 0]);

		if (p[0] == e1)
		{
			p[1] = e2;
			p[2] = e3;
			return;
		}
		if (p[0] == e2)
		{
			p[1] = e3;
			p[2] = e1;
			return;
		}
		if (p[0] == e3)
		{
			p[1] = e1;
			p[2] = e2;
			return;
		}
		assert(0); //problem
	}

	void checkForNull(std::vector<Point_3>& v, Pair e, std::vector<Point_3>& uvVector, std::vector<Point_3>& refVec)
	{
		if (v.size() > 0)
		{
			if (v[0] == refVec[refVec.size() - 1])
				return;
			assert(0);
		}
		if (uvVector[e.x] == refVec[refVec.size() - 1]) 
		{
			v.push_back(uvVector[e.x]);
			v.push_back(uvVector[e.y]);
		}
		else
		{
			v.push_back(uvVector[e.y]);
			v.push_back(uvVector[e.x]);
		}
		
	}

	void simplify(PointMap& pMap, std::vector<int>& fVec, std::vector<Point_3>& uvVector)
	{
		bool stopFlag = true;
		Pair p1[3], p2[3];	//represent the polygons we want to reduce in thier shared edges the number of points
		std::vector<Point_3> tempVec, tempVec2;	//the polygons that we check if they are simple will be here
		while (stopFlag)
		{
			stopFlag = false;
			auto mapIt = pMap.edgeToFace.begin();
			auto mapEnd = pMap.edgeToFace.end();
			while (mapIt != mapEnd)
			{
				p1[0] = Pair((mapIt->first).x, (mapIt->first).y);
				p2[0] = p1[0].reverse();
				setOrderOfEdges(p1, mapIt->second, fVec);
				setOrderOfEdges(p2, pMap.edgeToFace[p2[0]], fVec);

				std::vector<Point_3> &vecOfPointsP1 = pMap.getEdgeUVs(p1[0]);
				std::vector<Point_3> &vecOfPointsP2 = pMap.getEdgeUVs(p2[0]);
				std::vector<Point_3> & meshVecP1 = pMap.getNewMeshPointsVec(p1[0]);
				std::vector<Point_3> & meshVecP2 = pMap.getNewMeshPointsVec(p2[0]);

				tempVec.clear();
				tempVec2.clear();

				tempVec = vecOfPointsP1;
				tempVec2 = vecOfPointsP2;

				if (tempVec.size() == 0)
					assert(0);

				std::vector<Point_3>  vec21 = pMap.getEdgeUVs(p1[1]), vec31 = pMap.getEdgeUVs(p1[2]);	//p1 edges
				std::vector<Point_3>  vec22 = pMap.getEdgeUVs(p2[1]), vec32 = pMap.getEdgeUVs(p2[2]);	//p2 edges

				checkForNull(vec21, p1[1], uvVector, tempVec);
				checkForNull(vec22, p2[1], uvVector, tempVec2);
				tempVec.pop_back();
				tempVec2.pop_back();

				tempVec.insert(tempVec.end(), vec21.begin(), vec21.end());
				tempVec2.insert(tempVec2.end(), vec22.begin(), vec22.end());

				checkForNull(vec31, p1[2], uvVector, tempVec);
				checkForNull(vec32, p2[2], uvVector, tempVec2);
				tempVec.pop_back();
				tempVec2.pop_back();

				tempVec.insert(tempVec.end(), vec31.begin(), vec31.end());
				tempVec2.insert(tempVec2.end(), vec32.begin(), vec32.end());
				tempVec.pop_back();
				tempVec2.pop_back();

				Polygon_2 poly1, poly2;
				for (int j = 0; j < (int)tempVec.size(); ++j)
					poly1.push_back(Point_2(tempVec[j][0], tempVec[j][1]));
				for (int j = 0; j < (int)tempVec2.size(); ++j)
					poly2.push_back(Point_2(tempVec2[j][0], tempVec2[j][1]));

				int len = vecOfPointsP1.size();
				for (int i = 1; i < len - 1; ++i)
				{
					auto save = poly1[i];
					if (save != Point_2(vecOfPointsP1[i][0], vecOfPointsP1[i][1]))
						assert(0);
					poly1.erase(poly1.vertices_begin() + i);
					poly2.erase(poly2.vertices_begin() + (len - 1 - i));
					//vecOfPointsP1.erase(vecOfPointsP1.begin() + i);
					//vecOfPointsP2.erase(vecOfPointsP2.begin() + (len-1-i));
					
					/*/GMMDenseColMatrix testPoly(tempVec.size(), 2) , testPoly2(tempVec2.size(),2);

					

						//testPoly(j, 0) = tempVec[j][0];
						//testPoly(j, 1) = tempVec[j][1];
					//}
					

						//testPoly2(j, 0) = tempVec2[j][0];
						//testPoly2(j, 1) = tempVec2[j][1];
					//}

					//MatlabGMMDataExchange::SetEngineDenseMatrix("testPoly", testPoly);
					//MatlabGMMDataExchange::SetEngineDenseMatrix("testPoly2", testPoly2);
					//MatlabInterface::GetEngine().Eval("figure");
					//MatlabInterface::GetEngine().Eval("subplot(2,1,1);impoly(gca,testPoly);subplot(2,1,2);impoly(gca,testPoly2)");*/



					if ((poly1.is_simple()) && (poly2.is_simple()) && (poly1.orientation() == 1) && (poly2.orientation() == 1))
					{
						vecOfPointsP1.erase(vecOfPointsP1.begin() + i);
						vecOfPointsP2.erase(vecOfPointsP2.begin() + (len-1-i));
						len--;
						i--;
						meshVecP1.erase(meshVecP1.begin()+i);
						meshVecP2.erase(meshVecP2.begin() + (len-1-i) );
						stopFlag = true;
					}
					else
					{
						poly1.insert(poly1.vertices_begin() + i, save);
						poly2.insert(poly2.vertices_begin() + (len - 1 - i), save);
						//vecOfPointsP1.insert(vecOfPointsP1.begin() + i , save);
						//vecOfPointsP2.insert(vecOfPointsP2.begin() + (len - 1 - i) , save);
					}
				}

				mapIt++;
			}

		}
	
	}


	void extractIndicesFromPair(Pair& p, PointMap& pMap, std::vector<int>& polygonIndices)
	{
		polygonIndices.push_back(p.x);
		for (int i = 0; i < (int)pMap[p].size(); ++i)
			polygonIndices.push_back(pMap.getNewIndex(pMap[p][i]));
	}

	void extractIndicesFromPair(Pair& e1, Pair& e2, Pair& e3, PointMap& pMap, std::vector<int>& polygonIndices)
	{
		extractIndicesFromPair(e1, pMap, polygonIndices);
		extractIndicesFromPair(e2, pMap, polygonIndices);
		extractIndicesFromPair(e3, pMap, polygonIndices);
	}




	void barycentricCord(const std::vector<Point_3>& points, Point_3 point, Point_3 &res)
	{
		/*if ( 3 != points.size() )
		return;*/
		Point_3 p1 = points[0], p2 = points[1], p3 = points[2];

		Vector_3 vec1 = p1 - point;
		Vector_3 vec2 = p2 - point;
		Vector_3 vec3 = p3 - point;

		double s1 = crossProduct(vec3, vec2) / 2;
		double s2 = crossProduct(vec3, vec1) / 2;
		double s3 = crossProduct(vec1, vec2) / 2;
		double s = s1 + s2 + s3;

		Point_3 resp(s1 / s, s2 / s, s3 / s);
		res = resp;

		//std::cout<<"("<<res[0]<<","<<res[1]<<","<<res[2]<<")\t= "<<res[0]+res[1]+res[2]<<"\n";

		return;

	}

	double crossProduct(Vector_3 v1, Vector_3 v2)
	{
		double temp1 = v1.y() * v2.z() - v1.z() * v2.y();
		double temp2 = v1.x() * v2.z() - v1.z() * v2.x();
		temp2 = temp2*-1;
		double temp3 = v1.x() * v2.y() - v1.y() * v2.x();
		double temp = temp1*temp1 + temp2*temp2 + temp3*temp3;
		temp = std::sqrt(temp);

		return (temp);
	}









	void meanValueWeights(Mesh &source_mesh, GMMSparseRowMatrix &u, GMMSparseRowMatrix &weightsMat)
	{
		static bool isFirst = true;
		static std::vector<double> sourceBoundary;
		static Point_3 firstBoundaryVertex(RESET_NUM, RESET_NUM, RESET_NUM), secondBoundaryVertex(RESET_NUM, RESET_NUM, RESET_NUM);

		std::vector<Mesh::Halfedge_iterator> borderHDS;
		std::vector<int> verticesIndices;
		std::vector<double> partialLengths;
		double totalLength = 0, currLength = 0;
		Mesh::Vertex_const_handle currentV;
		Mesh::Halfedge_const_handle currH;

		source_mesh.getBorderHalfEdges(borderHDS);
		if (isFirst)
		{
			firstBoundaryVertex = borderHDS[0]->vertex()->uv();
			secondBoundaryVertex = borderHDS[1]->vertex()->uv();
		}

		for (int i = 0; i < (int)borderHDS.size(); ++i)
		{
			currH = borderHDS[i];
			currentV = currH->vertex();
			currLength = currH->length();
			totalLength += currLength;
			verticesIndices.push_back(currentV->index());
			partialLengths.push_back(totalLength);
		}

		int sourceMeshSize = source_mesh.size_of_vertices();
		int count = 0;
		int targetIndex = 0;
		if (!isFirst)
		{
			targetIndex = -1;
			for (int i = 0; i < (int)borderHDS.size(); ++i)
				if (borderHDS[i]->vertex()->point() == firstBoundaryVertex)
				{
				targetIndex = i;
				break;
				}
			assert(targetIndex != -1);
			if (borderHDS[(targetIndex + 2) % borderHDS.size()]->vertex()->point() != secondBoundaryVertex)
				assert(0);
		}

		u(verticesIndices[targetIndex], 0) = 1;
		int i = 1;
		int NN = borderHDS.size();

		while (i < NN)
		{
			if (isFirst)
			{
				u(verticesIndices[i], 0) = std::cos(2 * M_PI * partialLengths[i - 1] / totalLength);
				u(verticesIndices[i], 1) = std::sin(2 * M_PI * partialLengths[i - 1] / totalLength);
				sourceBoundary.push_back(std::cos(2 * M_PI * partialLengths[i - 1] / totalLength));
				sourceBoundary.push_back(std::sin(2 * M_PI * partialLengths[i - 1] / totalLength));
			}
			else
			{
				if (i != 1)
				{
					if (i + 1 != NN)
					{
						u(verticesIndices[(targetIndex + i) % NN], 0) = 0.5*sourceBoundary[count - 2] + 0.5*sourceBoundary[count];
						u(verticesIndices[(targetIndex + i) % NN], 1) = 0.5*sourceBoundary[count - 1] + 0.5*sourceBoundary[count + 1];
					}
					else
					{
						u(verticesIndices[(targetIndex + i) % NN], 0) = 0.5*sourceBoundary[count - 2] + 0.5;
						u(verticesIndices[(targetIndex + i) % NN], 1) = 0.5*sourceBoundary[count - 1];
					}
				}
				else		// the first 'u' is (1,0)
				{
					u(verticesIndices[(targetIndex + i) % NN], 0) = 0.5*sourceBoundary[count] + 0.5;
					u(verticesIndices[(targetIndex + i) % NN], 1) = 0.5*sourceBoundary[count + 1];
				}

				i++;
				if (i == NN)
					break;
				u(verticesIndices[(targetIndex + i) % NN], 0) = sourceBoundary[count];
				count++;
				u(verticesIndices[(targetIndex + i) % NN], 1) = sourceBoundary[count];
				count++;
			}
			++i;
		}

		// ********************************************************

		Mesh::Vertex_const_iterator currV, endV;
		currV = source_mesh.vertices_begin();
		endV = source_mesh.vertices_end();

		Mesh::Halfedge_around_vertex_const_circulator hdgAroundV;

		int vi, vj;
		double sumCot = 0, t1, t2;
		Mesh::Point_3 p1, p2, pi, pj;
		CGAL::Vector_3<Kernel> vv1, vv2, vv3, vv4;

		for (; currV != endV; ++currV) // traverse through all mesh vertices
		{
			vi = currV->index();
			pi = currV->point();
			if (currV->is_border())
			{
				weightsMat(vi, vi) = 1;
				continue;
			}

			hdgAroundV = currV->vertex_begin();

			do // traverse through connected edges
			{
				auto oppositeV = hdgAroundV->opposite()->vertex();
				vj = oppositeV->index(); // opposite vertex index for curr edge
				pj = oppositeV->point();
				p1 = hdgAroundV->next()->vertex()->point(); // 3rd point from first triangle
				p2 = hdgAroundV->opposite()->next()->vertex()->point(); // 3rd point from second triangle

				vv1 = p1 - pi;
				vv2 = pj - pi;
				vv3 = p2 - pi;
				//vv4 = p2 - pj;

				vv1 = normalizeVector(vv1); vv2 = normalizeVector(vv2);
				vv3 = normalizeVector(vv3); //vv4 = normalizeVector(vv4);
				double y_ij = std::acos(vv2*vv1);
				double l_ij = std::acos(vv3*vv2);
				y_ij = y_ij / 2;
				l_ij = l_ij / 2;
				t1 = std::tan(y_ij);//  std::sqrt(CGAL::cross_product(vv1, vv2).squared_length());
				t2 = std::tan(l_ij);// (vv3 * vv4) / std::sqrt(CGAL::cross_product(vv3, vv4).squared_length());
				if ((t1 < 0) || (t2 < 0))
					assert(0);
				weightsMat(vi, vj) = (t1 + t2) / ( 2*std::sqrt( vv2.squared_length() ) );
				sumCot += weightsMat(vi, vj);

				hdgAroundV++;
			} while (hdgAroundV != currV->vertex_begin());

			weightsMat(vi, vi) = -1 * sumCot;
			sumCot = 0;
		}


		isFirst = false;
	}


	double exactFaceOrientation(Facet_const_handle& face)
	{
		Vertex_const_handle p[3];
		face->getVertices(p);

		std::complex<double> p1(p[0]->uv().x(), p[0]->uv().y());
		std::complex<double> p2(p[1]->uv().x(), p[1]->uv().y());
		std::complex<double> p3(p[2]->uv().x(), p[2]->uv().y());

		double cross = CounterClockWise(p1, p2, p3);

		return cross;
	}

	int countFoldovers(const Mesh& mesh, bool printReport)
	{
		assert(mesh.is_pure_triangle());
		int numFlipped = 0;
		int numGood = 0;

		for_each_const_facet(face, mesh)
		{
			double cross = exactFaceOrientation(face);

			if (cross <= 0.0)
			{
				numFlipped++;
			}
			else
			{
				numGood++;
			}
		}
		int numTotal = numFlipped + numGood;

		if (printReport)
		{
			cout << "There are: " << numFlipped << " flipped triangles out of " << numTotal << " in total" << endl;
		}

		return numFlipped;
	}







