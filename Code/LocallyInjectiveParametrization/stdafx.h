#pragma  once
#define M_PI 3.1415926535897932384626433832795

// CGAL include
#include <CGAL/Timer.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Gmpq.h>
#include <CGAL/gmp.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Gmp_coercion_traits.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/number_utils.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include "CGAL_Vertex.h"
#include "CGAL_Halfedge.h"
#include "CGAL_Face.h"
#include "CGAL_Mesh.h"

#include <CGAL/Arrangement_on_surface_2.h>



// CGAL typedef
typedef CGAL::Simple_cartesian<double>	Kernel;

typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef CGAL::Polygon_2<Kernel>    Polygon_2;//Contour;
typedef CGAL::Triangle_2<Kernel> Triangle_2;

typedef CGAL::Polyhedron_3<Kernel>		pHedron;
typedef pHedron::HalfedgeDS				HDS;

typedef CGAL::Exact_predicates_exact_constructions_kernel		ARRKernel;
typedef ARRKernel::FT											ARRNumberType;
typedef CGAL::Arr_segment_traits_2<ARRKernel>					ARRTraits_2;
typedef ARRTraits_2::Point_2									EPoint_2;
typedef ARRKernel::Vector_2										EVector_2;
typedef ARRTraits_2::X_monotone_curve_2							ESegment_2;
typedef ARRTraits_2::Ray_2										ERay_2;
typedef CGAL::Arr_extended_dcel<ARRTraits_2,
								Mesh::Vertex_handle,
								Mesh::Halfedge_handle,
								Mesh::Facet_handle>				Dcel;
typedef CGAL::Arrangement_2<ARRTraits_2, Dcel>					Arrangement_2;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2>		Landmarks_pl;
typedef CGAL::Arr_default_overlay_traits<Arrangement_2>			Overlay_traits;
typedef CGAL::Polygon_2<ARRKernel>								EPolygon_2;




// other include

//#include "Face_index_observer.h"

#include "gmm/gmm.h"
#include "MatlabGMMDataExchange.h"
#include "MatlabInterface.h"
#include "GMM_Macros.h"

#include "Angle.h"
#include "Shor.h"


#include "wavefront_obj.h"


// general include 
#include <iostream>
using namespace std;
#include <string>
#include <map>
#include <windows.h>
#include <vector>
#include <complex>

#include "triangulation.h"
#include "cassert"

extern "C"
{
	#include "triangle.h"
}

#include "helpFunctions.h"

#include <CGAL/Sweep_line_2_algorithms.h>