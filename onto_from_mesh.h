void onto_mesh_delta(const kind::f);
void onto_mesh_delta_v(const kind::f);
//void onto_mesh_lumped(void);

void onto_mesh_full(linear& algebra, const kind::f scalarf);
void onto_mesh_full_v(linear& algebra, const kind::f vectorf);
void onto_mesh_flip(Triangulation& Tpart, Triangulation& Tmesh, bool FEM, const kind::f scalarf );

void from_mesh(const kind::f ) ;
void from_mesh_v(const kind::f ) ;

void from_mesh_full(linear& algebra_p, const kind::f scalarf);
void from_mesh_lumped(const kind::f scalarf);

void flip_volumes(Triangulation& Tpart, Triangulation& Tmesh, bool FEM);
