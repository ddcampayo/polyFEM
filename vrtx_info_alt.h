// Fixes errors on some compilers:

//typedef Vector_2 autoVector_2;

template<class K>
//class autoVector_2 : Vector_2 {
class autoVector_2 : public CGAL::Vector_2<K> {

 public:
//
 autoVector_2() : CGAL::Vector_2<K>(CGAL::NULL_VECTOR) {}
 autoVector_2(const Vector_2& v) : Vector_2(v) {}
//
//  //  operator Vector_2() {return Vector_2( this->x(), this->y()  ); }
//
};


