template <  class type >
class field {
 public:
  field() : f_()  {}
  type f() const {return f_;}
  type val() const {return f_;}
  type operator()() const {return f_;}

  void set(const type& ff) {f_ = ff;}
  void reset() {f_= type();}

  void operator =(const type& ff) {f_=ff;}
  type operator +(const type& ff) const {return f_ +ff;}
  type operator += (const type& ff) {type f2=  f_ + ff ; f_=f2 ; return f2;}
  type operator/=(const FT& ff) {return (f_ = (1.0)/ff * f_ );}
  type operator*=(const FT& ff) {return (f_ =       ff * f_ );}

 private:
  type f_;
};

//template<class K>
class autoVector_2 :
public
//CGAL::Vector_2<K> {
Vector_2 {
 public:

 autoVector_2() : Vector_2(CGAL::NULL_VECTOR) {}
 autoVector_2(const Vector_2& v) : Vector_2(v) {}

  //  operator Vector_2() {return Vector_2( this->x(), this->y()  ); }

};


namespace kind {
  enum f {P , PSTAR, VOL , GRADP , LAPLP , 
	  U, USTAR, UOLD, DIVU, LAPLU, ALPHA};
};

template < class Gt, class Vb >
class My_vertex_base
  : public  Vb
{
  typedef Vb                              Base;
public:
  typedef typename Vb::Vertex_handle      Vertex_handle;
  typedef typename Vb::Face_handle        Face_handle;
  typedef typename Vb::Point              Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef My_vertex_base<Gt,Vb2>                           Other;
  };

  typedef std::map<Vertex_handle,      FT>      scalar_link; 
  typedef std::map<Vertex_handle,autoVector_2>  vector_link; 
  //  typedef std::map<int,int> vector_link; 


private:
  scalar_link mass_; // vertices linked by coeffs of mass
  scalar_link stiff_;  // idem stiffness 
  scalar_link Delta_;  // idem Delta
  vector_link nabla_;  // idem nabla (vector)

public:
 My_vertex_base() :
  Base() {reset_fields(); }
 My_vertex_base(const Point & p) :
  Base(p) {reset_fields(); }
 My_vertex_base(const Point & p, Face_handle f) :
  Base(f,p) {reset_fields(); }
 My_vertex_base(Face_handle f) : 
  Base(f) {reset_fields(); }

  FT mass(const Vertex_handle& v, const FT& mm) {
    mass_[v] += mm;
    return mm;
  };

  FT stiff(const Vertex_handle& v, const FT& ss) {
    return stiff_[v] += ss;
  };


  FT Delta(const Vertex_handle& v, const FT& dd) {
    return Delta_[v] += dd;
  };

  Vector_2 nabla(const Vertex_handle& v, const Vector_2& vv) {
    Vector_2 v0= nabla_[v];
    return     nabla_[v] = v0+vv;
  };

  vector_link& nabla(void) {
    return nabla_;
  };

  void get_nabla(vector_link& nn) {  nn=nabla_;  };

  scalar_link mass() {
    return mass_;
  };

  scalar_link& stiff() {
    return stiff_;
  };

  scalar_link& Delta() {
    return Delta_;
  };

  // FT& sf(field::f w) { // Scalar field
  //   if(w==field::P) return p_;
  // }

  // Vector_2& vf(field::f w) { // Vector field
  //   if(w==field::GRADP) return gradp_;
  // }
  

  typedef field<FT>            scalar_field;
  typedef field<autoVector_2>  vector_field; // TODO: make this autoVector_2
  typedef field<Point>         point_field; // TODO: make this autoVector_2
  scalar_field p;          //scalar field (e.g pressure)
  scalar_field pstar;
  scalar_field vol;
  scalar_field fvol;
  scalar_field alpha;

  vector_field U;
  vector_field Uold;
  vector_field Ustar;
  vector_field disp;

  point_field rold;

  vector_field gradp;          //vector field (e.g grad of pressure)

  scalar_field laplp;          //scalar field (e.g Laplacian of pressure)
  scalar_field divU;

  vector_field laplU;

  field<int> idx;       // an index

  field<bool> moved;       // whether it's moved or not

  // aux fields:
  scalar_field p_1;
  vector_field Uold_1;

  void reset_fields(void) { // set fields to zero
    p.reset();
    pstar.reset();
    vol.reset();
    fvol.reset();
    gradp.reset();
    laplp.reset();
    U.reset();
    Uold.reset();
    Ustar.reset();
    divU.reset();
    laplU.reset();
  }


  scalar_field& sf(const kind::f a) {
    if(a==kind::P)
      return p;
    else if(a==kind::VOL)
      return vol;
    else if(a==kind::PSTAR)
      return pstar;
    else if(a==kind::LAPLP)
      return laplp;
    else if(a==kind::DIVU)
      return divU;
    else if(a==kind::ALPHA)
      return alpha;
  }

  vector_field& vf(const kind::f a) {
    if(a==kind::GRADP)
      return gradp;
    else if(a==kind::U)
      return U;
    else if(a==kind::USTAR)
      return Ustar;
    else if(a==kind::UOLD)
      return Uold;
    else if(a==kind::LAPLU)
      return laplU;
  }
 

  
};





template < typename Gt, typename Fb >
  class My_face_base
  : public Fb
{
 public:
 typedef Gt                                           Geom_traits;
 typedef typename Fb::Vertex_handle                   Vertex_handle;
 typedef typename Fb::Face_handle                     Face_handle;

 template < typename TDS2 >
 struct Rebind_TDS {
   typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
   typedef My_face_base<Gt, Fb2>             Other;
 };

 public:
 My_face_base()
 : Fb() {}

 My_face_base(Vertex_handle v0, 
			   Vertex_handle v1, 
			   Vertex_handle v2)
 : Fb(v0,v1,v2) {}

 My_face_base(Vertex_handle v0, 
			   Vertex_handle v1, 
			   Vertex_handle v2,
			   Face_handle n0, 
			   Face_handle n1, 
			   Face_handle n2)
 : Fb(v0,v1,v2,n0,n1,n2) {}

 static int ccw(int i) {return CGAL::Triangulation_cw_ccw_2::ccw(i);}
 static int  cw(int i) {return CGAL::Triangulation_cw_ccw_2::cw(i);}

#ifndef CGAL_NO_DEPRECATED_CODE
 Vertex_handle mirror_vertex(int i) const;
 int mirror_index(int i) const;
#endif


 FT area;

 // typedef std::vector<int> v_idx;

 typedef std::vector<Vertex_handle> v_vrtcs;
//   typedef std::pair<int,double> coef;

 typedef  std::vector<v_vrtcs> v_connects;

 typedef std::vector<FT> v_A;

 typedef std::vector<v_A> v_coeffs;

 // std::vector<v_vrtcs>
 v_connects connects;

 v_coeffs  coeffs;
 
 std::vector<FT> ll;

//   
//   typedef v_coef edge_info;

// // typedef std::map<Vertex_handle,vrtx_info> vrtxmap;

//   typedef std::map<Edge,edge_info> edgemap;
// }

// edge::edgemap ed_map;


};
