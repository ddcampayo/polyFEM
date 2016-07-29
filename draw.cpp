#include"main.h"
#include"sim_pars.h"

extern sim_pars simu;



void draw(void) {

  cout << "draw on step  "<< simu.current_step() << endl;

  std::stringstream  mkdir;
  std::stringstream  dirname;

  dirname << simu.current_step();
  mkdir << "mkdir -p "  << dirname.str();
  cout << "running : " << mkdir.str() << endl;
  system(mkdir.str().c_str());

  std::stringstream  cp_cfg;
  cp_cfg << "cp -f simu.cfg " << dirname.str();
  system(cp_cfg.str().c_str());

  {
  std::stringstream  namefile;
  namefile << simu.current_step() << '/' << "mesh.dat";

  cout << "writing on file : " << namefile.str() << endl;
  std::ofstream main_data;
  main_data.open(namefile.str().c_str() );

  for(F_v_it vit=Tm.vertices_begin();
      vit != Tm.vertices_end();
      vit++) {

    Periodic_point pp=Tm.periodic_point(vit);
    Point p=Tm.point(pp);

    #include"printout.h"

  }

  main_data << endl;

  main_data.close();
  }
  {
  std::stringstream  namefile;
  namefile << simu.current_step() << '/' << "particles.dat";

  cout << "writing on file : " << namefile.str() << endl;
  std::ofstream main_data;
  main_data.open(namefile.str().c_str() );

  for(F_v_it vit=Tp.vertices_begin();
      vit != Tp.vertices_end();
      vit++) {

    Periodic_point pp=Tp.periodic_point(vit);
    Point p=Tp.point(pp);

    #include"printout.h"

  }

  main_data << endl;

  main_data.close();
  }


  std::stringstream  vtkfile;
  vtkfile << simu.current_step() << '/' << "points.vtu";

  cout << "writing on file : " << vtkfile.str() << endl;


//
//  write_c2t3_to_vtk_xml_file.h
//
//  Created by David Bernstein on 5/1/13.
//


  std::ofstream  vtk_file(vtkfile.str().c_str() );
  
  // header
  vtk_file << "<VTKFile type=\"UnstructuredGrid\" ";
  vtk_file << "version=\"0.1\" ";
  vtk_file << "byte_order=\"BigEndian\">" << std::endl;
        
  int indent_size = 2;
  std::string indent_unit(indent_size, ' ');
  std::string indent = indent_unit;
  vtk_file << indent + "<UnstructuredGrid>" << std::endl;
        
  int num_vertices = Tp.number_of_vertices();
  int num_cells = Tp.number_of_faces();
        
  indent += indent_unit;
  vtk_file << indent + "<Piece NumberOfPoints=\"" << num_vertices << "\" ";
  vtk_file << "NumberOfCells=\"" << num_cells << "\">" << std::endl;
        
  // Write vertices
  indent += indent_unit;
  vtk_file << indent + "<Points>" << std::endl;
        
  indent += indent_unit;
  vtk_file << indent;
  vtk_file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
        
  std::map<Point, int> V;
  int i=0;
  indent += indent_unit;
        
  for (F_v_it it=Tp.finite_vertices_begin(); it != Tp.finite_vertices_end(); ++it)
    {
      vtk_file << indent;
      vtk_file << it->point().x() << " " << it->point().y() << " " << 0 << std::endl;
      V[it->point()] = i;
      ++i;
    }
        
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</DataArray>" << std::endl;

  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</Points>" << std::endl;

  // Write tetrahedra
  vtk_file << indent << "<Cells>" << std::endl;
   
  indent += indent_unit;
  vtk_file << indent;
  vtk_file << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">";
  vtk_file << std::endl;
        
  indent += indent_unit;
  F_f_it it;

  for (it = Tp.finite_faces_begin(); it != Tp.finite_faces_end(); ++it)
    {
	vtk_file << indent;
	vtk_file << V[it->vertex(0)->point()] << " ";
	vtk_file << V[it->vertex(1)->point()] << " ";
	vtk_file << V[it->vertex(2)->point()] << " " << std::endl;
    }
        
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</DataArray>" << std::endl;

  // offsets
  // every element is a three node triangle so all offsets are multiples of 3
  vtk_file << indent;
  vtk_file << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">";
  vtk_file << std::endl;
  i = 3;
  indent += indent_unit;
  for (int j = 0; j < num_cells; ++j)
    {
      vtk_file << indent << i << std::endl;
      i += 3;
    }
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</DataArray>" << std::endl;
        
  // cell types (type 5 is a three node triangle)
  vtk_file << indent;
  vtk_file << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">";
  vtk_file << std::endl;
  indent += indent_unit;
  for (int j = 0; j < num_cells; ++j)
    {
      vtk_file << indent << "5" << std::endl;
    }
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</DataArray>" << std::endl;
  
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</Cells>" << std::endl;
  
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</Piece>" << std::endl;
        
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << indent + "</UnstructuredGrid>" << std::endl;
        
  indent.erase(indent.length()-indent_size, indent_size);
  vtk_file << "</VTKFile>" << std::endl;

  vtk_file.close();
        
  return;

}
