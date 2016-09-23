#include"main.h"
#include"periodic.h"

//extern const FT LL;

FT per_dist(const FT& x1, const FT& x2) {
  FT dx=x2-x1;

  if(fabs(dx)>LL/2.0)
    if(dx>0)
      dx -= LL;
    else
      dx += LL;
  return dx;
}

Point per_point(const Point& p) {
  return per_point( p.x(),p.y() );
}

Point per_point(const FT& x10,const  FT& x20) {

  FT x1=x10;
  FT x2=x20;

  if (x1 > LL/2.0)
    x1 -= LL;
  else if (x1 < -LL/2.0)
    x1 += LL;

  if (x2 > LL/2.0)
    x2 -= LL;
  else if(x2 < -LL/2.0)
    x2 += LL;

  return Point(x1,x2);
}


Vector_2 per_vect(const Point& P1, const Point& P2) {
  FT dx=per_dist(P1.x(),P2.x());
  FT dy=per_dist(P1.y(),P2.y());

  return Vector_2(dx,dy);
}

FT per_dist2(const Point& P1, const Point& P2) {
  FT dx=per_dist(P1.x(),P2.x());
  FT dy=per_dist(P1.y(),P2.y());

  return dx*dx+dy*dy;
}

FT per_area(const Point& P0, const Point& P1, const Point& P2) {
  Vector_2 v1=per_vect(P0,P1);
  Vector_2 v2=per_vect(P0,P2);

  return fabs((
	       v1.x()*v2.y()-
	       v1.y()*v2.x())/2.0);
}







