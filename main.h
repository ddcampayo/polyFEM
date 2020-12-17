//#define DEBUG
#undef DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/spatial_sort.h>

using std::vector;
using std::list;
using std::endl;
using std::cout;
using std::cin;


#include <fstream>

//using namespace CGAL;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//typedef CGAL::Periodic_2_triangulation_filtered_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;

typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT>    VbDS;

//typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT>    VbDS;

typedef CGAL::Periodic_2_triangulation_face_base_2<GT>      FbDS;


typedef K::FT FT;

typedef CGAL::Vector_2<K>       Vector_2;
typedef CGAL::Segment_2<K>      Segment;
typedef CGAL::Triangle_2<K>     Triangle;
typedef CGAL::Point_2<K>        Point;

typedef CGAL::Creator_uniform_2<FT,Point> Creator;


#include"vrtx_info.h"


typedef My_vertex_base<GT, VbDS>                             Vb;
typedef My_face_base<GT, FbDS>                               Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb>         Tds;

//#include"my_Tr.h"
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds>   Triangulation;


typedef Triangulation::Face_handle    Face_handle;
typedef Triangulation::Face           Face;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Vertex         Vertex;
typedef Triangulation::Edge           Edge;
//typedef Triangulation::Facet          Facet;
typedef Triangulation::Locate_type    Locate_type;
//typedef Triangulation::Point          Point;

typedef Triangulation::Finite_vertices_iterator   F_v_it;
typedef Triangulation::Finite_edges_iterator      F_e_it;
typedef Triangulation::Finite_faces_iterator      F_f_it;
typedef Triangulation::Edge_circulator            Edge_circulator;
typedef Triangulation::Vertex_circulator          Vertex_circulator;
typedef Triangulation::Face_circulator            Face_circulator;
typedef Triangulation::Iso_rectangle     Iso_rectangle;


typedef Triangulation::Periodic_point           Periodic_point;
typedef Triangulation::Periodic_point_iterator  Periodic_point_iterator;

typedef Triangulation::Periodic_segment           Periodic_segment;
typedef Triangulation::Periodic_segment_iterator  Periodic_segment_iterator;

typedef Triangulation::Periodic_triangle           Periodic_triangle;
typedef Triangulation::Periodic_triangle_iterator  Periodic_triangle_iterator;

typedef Triangulation::Iterator_type               Iterator_type;

const Iterator_type stored=Triangulation::STORED;
const Iterator_type unique=Triangulation::UNIQUE;
const Iterator_type stored_cover=Triangulation::STORED_COVER_DOMAIN;
const Iterator_type unique_cover=Triangulation::UNIQUE_COVER_DOMAIN;


typedef CGAL::Periodic_2_offset_2 Offset;

void draw(Triangulation& T,  const std::string file_name , const bool setup ) ;
void number(Triangulation& T);
void quad_coeffs( Triangulation& T , bool ) ;
void matrices(void);
FT solve_linear(std::vector<Periodic_point>& vP,
		  std::vector<FT>& AA );

void prune(std::vector<Vertex_handle>& v);
void set_fields_Zalesak(void);
void set_vels(void);
void curvature(void) ;
FT move(Triangulation& T, FT, FT& ) ;
void move_info(Triangulation& T) ;
void u_star(Triangulation& T, FT, bool) ;
void u_new(Triangulation& T, const FT) ;
//void u_star_new(Triangulation& T, FT) ;
void update_half_velocity( Triangulation& T , const bool );
void update_half_alpha( Triangulation& T );
void nabla(Triangulation& T) ;
void Delta(Triangulation& T);

void setup_v(void);


void areas(Triangulation& T);
void volumes(Triangulation& T, bool);
void integrals(Triangulation& T, std::ofstream& log_file);
void fidelity(Triangulation& T, std::ofstream& log_file );

// TODO: don't have these around as global vars
//extern Triangulation Tm;
//extern Triangulation Tp;

//#include <Eigen/Dense>
// #include"linear.h"


