int ii;

// 2
main_data >> ii;
vh->idx.set(ii);

FT ff;

//vol 3
main_data >> ff;
vh->vol.set(ii);


// alpha 4
main_data >> ff;
vh->alpha.set(ii);

// p 5

main_data >> ff;
vh->p.set(ff);

//grad p 6, 7
main_data >> x;
main_data >> y;

// U 8, 9
main_data >> x;
main_data >> y;
vh->U.set(Vector_2(x,y));

// div U , 10
main_data >> ff;

// lapl U 11, 12
main_data >> x;
main_data >> y;

// Ustar 13,14 
main_data >> x;
main_data >> y;

// Uold 15,16 

main_data >> x;
main_data >> y;
vh->Uold.set(Vector_2(x,y));

// p star 17
main_data >> ff;

