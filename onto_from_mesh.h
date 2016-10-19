void onto_mesh_delta( Triangulation& Tfrom, Triangulation& Tto, const kind::f);
void onto_mesh_delta_v( Triangulation& Tfrom, Triangulation& Tto, const kind::f);
//void onto_mesh_lumped(void);

void onto_mesh_full( Triangulation& Tfrom, Triangulation& Tto, linear& algebra, const kind::f scalarf);
void onto_mesh_full_v( Triangulation& Tfrom, Triangulation& Tto, linear& algebra, const kind::f vectorf);
void onto_mesh_flip(Triangulation& Tpart, Triangulation& Tmesh, bool FEM, const kind::f scalarf );
void onto_mesh_flip_v(Triangulation& Tpart, Triangulation& Tmesh, bool FEM, const kind::f vectorf );

void from_mesh( Triangulation& Tfrom, Triangulation& Tto,const kind::f ) ;
void from_mesh_v( Triangulation& Tfrom, Triangulation& Tto,const kind::f ) ;

void from_mesh_full( Triangulation& Tfrom, Triangulation& Tto,linear& algebra_p, const kind::f scalarf);
void from_mesh_full_v( Triangulation& Tfrom, Triangulation& Tto,linear& algebra_p, const kind::f vectorf);
void from_mesh_lumped( Triangulation& Tfrom, Triangulation& Tto,const kind::f scalarf);
void from_mesh_lumped_v( Triangulation& Tfrom, Triangulation& Tto,const kind::f vectorf);

void flip_volumes(Triangulation& Tpart, Triangulation& Tmesh, bool FEM);
