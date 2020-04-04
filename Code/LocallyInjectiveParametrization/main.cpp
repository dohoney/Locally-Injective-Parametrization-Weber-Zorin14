#include "stdafx.h"

void run();
const std::string currentDateTime();

void main()
{
	//std::ofstream out("log.txt");
	//std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	//std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
	std::cout << "****************\nProgram start at: " << currentDateTime() <<"\n";
	run();
	std::cout << "****************";
}

void run()
{
	ofstream logFile;
	logFile.open("log.txt");
	double sumTime = 0;
	std::vector<Kernel::Point_3> pVec;
	std::vector<int> fVec;
	Mesh source_mesh;
	loadSourceMesh( source_mesh ,pVec , fVec);

	logFile << "Mesh loaded successfully.\n# of vertices: " << pVec.size() << "\n# of faces: " << fVec.size()/3 << "\n\n" ;

	std::vector<Mesh::Halfedge_iterator> border;
	source_mesh.getBorderHalfEdges( border );
	int numOfBorder = (int)border.size();

	//--------------get data from matlab about the target polygon-----------
	MatlabInterface::GetEngine().Eval("nis2");
	GMMDenseColMatrix target_size;
	MatlabGMMDataExchange::GetEngineDenseMatrix("n_bSize" , target_size);
	GMMDenseColMatrix targetVertices((int)target_size(0, 0), 2), rotIndices(1, (int)target_size(0, 0));
	MatlabGMMDataExchange::GetEngineDenseMatrix("n_b" , targetVertices);
	MatlabGMMDataExchange::GetEngineDenseMatrix("rotIndices", rotIndices);

	// set full rotation indices array
	int *rArr = new int[(int)target_size(0, 0)];
	for (int i = 0; i < (int)target_size(0, 0); ++i)
		rArr[i] = rotIndices(0, i);
	
	// calculate avg length of source edges on border
	double avg_arc = 0;
	
	for ( int i = 0; i < numOfBorder; ++i )
		avg_arc += border[i]->length();
	
	avg_arc = avg_arc / numOfBorder;
	//----------------------------------------------------------------------------------------------------------------------
	Polygon_2 poly,bPoly;

	for ( int i = 0; i < target_size(0,0); ++i )
		poly.push_back( Point_2( targetVertices(i,0) , targetVertices(i,1) ) );

	bPoly = poly;
	addPointsToTarget( bPoly , numOfBorder , avg_arc );

	CGAL::Timer triangulateTimer;
	triangulateTimer.start();
	std::cout << "Triangulate target polygon...\n";
	Shor shor;
	shor.setSourceMinArc(avg_arc);
	shor.setSourceArea( source_mesh.area() );
	shor.numOfWantedTriangles = 2*source_mesh.size_of_facets();
	shor.load_polygon (poly , bPoly);
	shor.load_rArray(rArr);
	//shor.load_rArray(NULL);
	shor.play();
	shor.build_triangulation();
	std::cout << "Done!\n";
	shor.simplify_triangulation();
	if (!shor.isTriangultae)	//fail to triangulate target polygon
	{
		std::cout << "Error: the target polygon is not self-overlapping polygon! \n";
		logFile << "Error: the target polygon is not self-overlapping polygon! \nIf you think it's indeed SOP, try to choose the 'reverse boundary orientation' option, or change the rotation indices.\n";
		logFile.close();
		return;		
	}

	triangulateTimer.stop();
	std::cout << "Total time to generate mesh from the target polygon: " << triangulateTimer.time() << " seconds\n";
	logFile << "Total time to generate mesh from the target polygon: " << triangulateTimer.time() << " seconds\n";
	sumTime += triangulateTimer.time();
	//if we got here that means the target mesh is set
	//***********************************************

	MatlabInterface::GetEngine().Eval("nis3");
	GMMDenseColMatrix weightsSelect(1, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("weightsSelect", weightsSelect);
	bool isSourceHarmonic = weightsSelect(0, 0) == 1;
	bool isTargetHarmonic = weightsSelect(0, 1) == 1;

	int sourceMeshSize = source_mesh.size_of_vertices();
	int targetMeshSize = shor.target_mesh.size_of_vertices();
	GMMSparseRowMatrix uSource(sourceMeshSize,2),uTarget(targetMeshSize,2);
	GMMSparseRowMatrix weightsMatSource(sourceMeshSize,sourceMeshSize), weightsMatTarget(targetMeshSize,targetMeshSize);

	std::cout << "Mapping source and target meshes to the unit disk... \n";
	CGAL::Timer harmonicTimer;
	harmonicTimer.start();
	HarmonicFlattening(source_mesh, uSource, weightsMatSource, isSourceHarmonic);
	HarmonicFlattening(shor.target_mesh, uTarget, weightsMatTarget, isTargetHarmonic);
	//meanValueWeights(source_mesh, uSource, weightsMatSource);
	//meanValueWeights(shor.target_mesh, uTarget, weightsMatTarget);
	harmonicTimer.stop();
	std::cout << "Done!\n";
	
	MatlabGMMDataExchange::SetEngineSparseMatrix( "uSource" , uSource );
	MatlabGMMDataExchange::SetEngineSparseMatrix( "weightsMatSource" , weightsMatSource );
	MatlabGMMDataExchange::SetEngineSparseMatrix( "uTarget" , uTarget );
	MatlabGMMDataExchange::SetEngineSparseMatrix( "weightsMatTarget" , weightsMatTarget );
	MatlabInterface::GetEngine().Eval("nis4");
	//*************************************************
	GMMDenseColMatrix sTime(1, 1), tTime(1, 1);
	GMMDenseColMatrix sourceMap(sourceMeshSize, 2), targetMap(targetMeshSize, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("outSource", sourceMap);
	MatlabGMMDataExchange::GetEngineDenseMatrix("outTarget", targetMap);
	MatlabGMMDataExchange::GetEngineDenseMatrix("sTime", sTime);
	MatlabGMMDataExchange::GetEngineDenseMatrix("tTime", tTime);
	std::cout << "Total time to construct the 2 harmonic maps: " << harmonicTimer.time() + sTime(0, 0) + tTime(0, 0) << " seconds\n";
	logFile << "Total time to construct the 2 harmonic maps (to the unit disk): " << harmonicTimer.time() + sTime(0, 0) + tTime(0, 0) << " seconds\n";
	sumTime += harmonicTimer.time() + sTime(0, 0) + tTime(0, 0);

	Arrangement_2 arrSource,arrTarget;
	Landmarks_pl sourceLm,targetLm;
	//Face_index_observer sourceObs(arrSource),targetObs(arrTarget);

	std::vector<EPoint_2> sourceHarmonicMapPoints,targetHarmonicMapPoints;
	sourceHarmonicMapPoints.resize(sourceMeshSize);
	targetHarmonicMapPoints.resize(targetMeshSize);
	for (int i = 0; i < sourceMeshSize; ++i)
	{
		sourceHarmonicMapPoints[i] = EPoint_2(sourceMap(i, 0), sourceMap(i, 1));
		//std::cout << "(" << sourceHarmonicMapPoints[i].x() << "," << sourceHarmonicMapPoints[i].y() << ")\n";
	}
	for ( int i = 0; i < targetMeshSize; ++i )
		targetHarmonicMapPoints[i] = EPoint_2(targetMap(i,0),targetMap(i,1)) ;

	auto vItSource = source_mesh.vertices_begin();
	while (vItSource != source_mesh.vertices_end())
	{
		int i = vItSource->index();
		vItSource->uv() = Point_3(sourceMap(i, 0) , sourceMap(i, 1) , 0);
		vItSource++;
	}
	auto vItTarget = shor.target_mesh.vertices_begin();
	while (vItTarget != shor.target_mesh.vertices_end())
	{
		int i = vItTarget->index();
		vItTarget->uv() = Point_3(targetMap(i, 0), targetMap(i, 1), 0);
		vItTarget++;
	}


	CGAL::Timer arrangementBuildTimer;
	arrangementBuildTimer.start();
	BuildArrangement ( arrSource , sourceLm , sourceHarmonicMapPoints , fVec , /*sourceObs ,*/ source_mesh);
	BuildArrangement ( arrTarget , targetLm , targetHarmonicMapPoints , shor.fVec , /*targetObs ,*/ shor.target_mesh);
	arrangementBuildTimer.stop();
	std::cout << "Total time to build the CGAL Arrangements: " << arrangementBuildTimer.time() << " seconds\n";
	logFile << "Total time to build the CGAL Arrangements: " << arrangementBuildTimer.time() << " seconds\n";
	sumTime += arrangementBuildTimer.time();
	//matchPointsIndices( arrSource , sourceMap );
	//matchPointsIndices( arrTarget , targetMap );
	
	CGAL::Timer buildMapTimer;
	buildMapTimer.start();

	std::vector<Point_3> uvVector;
	uvVector.resize(sourceMeshSize);
	setBoundaryUV(source_mesh, shor.target_mesh, uvVector);
	std::cout << "Calculating new UV's... \n";
	std::vector<int> neg = updateUVs(source_mesh, shor.target_mesh, arrSource, arrTarget, sourceLm, targetLm, /*sourceObs, targetObs,*/ sourceHarmonicMapPoints, fVec, uvVector);
	std::cout << "Done!\n";

	int aa = refine(neg, source_mesh, shor.target_mesh, arrSource, arrTarget, sourceLm, targetLm, /*targetObs,*/ sourceHarmonicMapPoints, pVec, fVec, uvVector);
	
	buildMapTimer.stop();
	std::cout << "Total time of composition and refinement: " << buildMapTimer.time() << " seconds\n";
	logFile << "Total time of composition and refinement: " << buildMapTimer.time() << " seconds\n";
	sumTime += buildMapTimer.time();
	logFile << "\n# of new points: " << aa << "\n\nTotal run time: " << sumTime <<"\n";
	

	GMMDenseColMatrix finalOut(uvVector.size(), 2);
	//auto it = source_mesh.vertices_begin();
	//int indexUV;
	/*while ( it!= source_mesh.vertices_end() )
	{
		indexUV = it->index();
		finalOut(indexUV,0) = it->uv().x();
		finalOut(indexUV,1) = it->uv().y();
		it++;
	}*/
	for (int i = 0; i < (int)uvVector.size(); ++i)
	{
		finalOut(i, 0) = uvVector[i].x();
		finalOut(i, 1) = uvVector[i].y();
	}
	
	std::vector<int> newFvec;
	for (int i = 0; i < fVec.size(); i=i+3)
	{
		if ((fVec[i] == -1) && (fVec[i + 1] == -1) && (fVec[i + 2] == -1))
			continue;
		newFvec.push_back(fVec[i]);
		newFvec.push_back(fVec[i+1]);
		newFvec.push_back(fVec[i+2]);
	}

	GMMDenseColMatrix finalFvec(newFvec.size()/3, 3);
	int j = 0;
	for (int i = 0; i < newFvec.size(); i=i+3)
	{
		finalFvec(j, 0) = newFvec[i];
		finalFvec(j, 1) = newFvec[i+1];
		finalFvec(j, 2) = newFvec[i+2];
		j++;
	}

	GMMDenseColMatrix finalPvec(pVec.size(), 3);
	for (int i = 0; i < (int)pVec.size(); ++i)
	{
		finalPvec(i, 0) = pVec[i].x();
		finalPvec(i, 1) = pVec[i].y();
		finalPvec(i, 2) = pVec[i].z();
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix( "finalOut" , finalOut );
	MatlabGMMDataExchange::SetEngineDenseMatrix("finalFvec", finalFvec);
	MatlabGMMDataExchange::SetEngineDenseMatrix("finalPvec", finalPvec);
	MatlabInterface::GetEngine().Eval("finalFvec = finalFvec +1");
	MatlabInterface::GetEngine().Eval("resMap");

	delete[] rArr;
	logFile.close();
}


const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}