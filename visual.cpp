#include "helper.hpp"
#include "visual.hpp"
#include <stdio.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTuple.h>
#include <vtkStructuredGridWriter.h>
#include "vector"

void write_gpu_vtkFile(const char *szProblem,
		 SimulationContext& sim_ctx,
         GpuGrid& grid) {
  
  int i,j;
  int il = 0;
  int ir = grid._imax - 2;
  int jb = 0;
  int jt = grid._jmax - 2;
  matrix2dRM<float> U(grid._imax,grid._jmax, 0), V(grid._imax, grid._jmax, 0), P(grid._imax, grid._jmax, 0), T(grid._imax, grid._jmax, 0), D(grid._imax, grid._jmax, 0);
  matrix2dRM<std::pair<float, float>> uv(grid._imax, grid._jmax, {});
  std::vector<matrix2dRM<float>> s(sim_ctx.species.size(), matrix2dRM<float>(grid._imax, grid._jmax, 0));
  grid.temperature.get_cur_image(T.data(), T.size() * sizeof(float));
  grid.pressure.get_image(P.data(), P.size() * sizeof(float));
  grid.velocity.get_image(uv.data(), uv.size() * sizeof(std::pair<float, float>));
  grid.density.get_image(D.data(), D.size() * sizeof(float));
  for (int i = 0; i < sim_ctx.species.size(); ++i) {
      grid.species_concentrations.get_cur_image(i, s[i].data(), s[i].size() * sizeof(float));
  }
  for (int w = 0; w < uv.width(); ++w) {
      for (int h = 0; h < uv.height(); ++h) {
          U(w, h) = uv(w, h).first;
          V(w, h) = uv(w, h).second;
      }
  }
  char szFileName[150];
  FILE *fp=NULL;
  sprintf( szFileName, "%s.%i.vtk", szProblem, sim_ctx.output_count );
  fp = fopen( szFileName, "w");
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
  }

  write_vtkHeader( fp, il, ir, jb, jt, sim_ctx.dx, sim_ctx.dy);
  write_vtkPointCoordinates(fp, il, ir, jb, jt, 0, 0, sim_ctx.dx, sim_ctx.dy);

  fprintf(fp,"POINT_DATA %i \n", (ir - il + 1)*(jt - jb + 1) );
	
  fprintf(fp,"\n");
  fprintf(fp, "VECTORS velocity float\n");
  for(j = jb; j < jt+1; j++) {
    for(i = il; i < ir+1; i++) {
      if(grid.cell_flagg(i,j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i,j)>9)
         fprintf(fp, "%f %f 0\n", (U(i,j) + U(i,j+1)) * 0.5, (V(i,j) + V(i+1,j)) * 0.5 );
	  else
		 fprintf(fp, "0 0 0\n");
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"CELL_DATA %i \n", ((ir - il)*(jt - jb)) );
  fprintf(fp, "SCALARS pressure float 1 \n"); 
  fprintf(fp, "LOOKUP_TABLE default \n");
  for(j = jb + 1; j <= jt; j++) {
    for(i = il + 1; i <= ir; i++) {
		if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i,j)>9)
			fprintf(fp, "%f\n", P(i,j));
		else
			fprintf(fp, "0\n");
    }
  }

  fprintf(fp, "\n");
  fprintf(fp, "SCALARS cell_type int 1 \n");
  fprintf(fp, "LOOKUP_TABLE default \n");
  for (j = jb + 1; j <= jt; j++) {
      for (i = il + 1; i <= ir; i++) {
          if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW|| (int)grid.cell_flagg(i,j)>9)
            fprintf(fp, "%d\n", grid.cell_flagg(i, j));
          else
          	fprintf(fp, "0\n");
      }
  }

  fprintf(fp, "\n");
  fprintf(fp, "SCALARS temperature float 1 \n");
  fprintf(fp, "LOOKUP_TABLE default \n");
  for (j = jb + 1; j <= jt; j++) {
	for (i = il + 1; i <= ir; i++) {
		if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW|| (int)grid.cell_flagg(i,j)>9)
			fprintf(fp, "%f\n", T(i,j));
		else
			fprintf(fp, "0\n");
	}
  }

  fprintf(fp, "\n");
  fprintf(fp, "SCALARS density float 1 \n");
  fprintf(fp, "LOOKUP_TABLE default \n");
  for (j = jb + 1; j <= jt; j++) {
      for (i = il + 1; i <= ir; i++) {
          if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW|| (int)grid.cell_flagg(i,j)>9)
              fprintf(fp, "%f\n", D(i, j));
          else
              fprintf(fp, "0\n");
      }
  }

  for (int c = 0; c < s.size(); ++c) {
      fprintf(fp, "\n");
      fprintf(fp, "SCALARS ");
      fprintf(fp, (sim_ctx.species[c].name + " float 1 \n").c_str());
      fprintf(fp, "LOOKUP_TABLE default \n");
      for (j = jb + 1; j <= jt; j++) {
          for (i = il + 1; i <= ir; i++) {
              if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i,j)>9)
                  fprintf(fp, "%f\n", s[c](i, j));
              else
                  fprintf(fp, "0\n");
          }
      }
  }

  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    ERROR( szBuff );
  }
}

void write_vtkFile(const char* szProblem,
    SimulationContext& sim_ctx,
    Grid& grid) {

    int i, j;
    int il = (sim_ctx.omg_i) ? 1 : 1;
    int ir = (sim_ctx.omg_i == sim_ctx.iproc - 1) ? grid.pressure.width() - 2 : grid.pressure.width() - 2;
    int jb = (sim_ctx.omg_j) ? 1 : 1;
    int jt = (sim_ctx.omg_j == sim_ctx.jproc - 1) ? grid.pressure.height() - 2 : grid.pressure.height() - 2;
    matrix2d<double>& U = grid.velocity[0], & V = grid.velocity[1], & P = grid.pressure, & T = grid.temperature;
    char szFileName[150];
    FILE* fp = NULL;
    sprintf(szFileName, "%s_%i_%i.%i.vtk", szProblem, sim_ctx.omg_i, sim_ctx.omg_j, sim_ctx.output_count);
    fp = fopen(szFileName, "w");
    if (fp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }

    write_vtkHeader(fp, il, ir, jb, jt, sim_ctx.dx, sim_ctx.dy);
    write_vtkPointCoordinates(fp, il, ir, jb, jt, sim_ctx.il - 2 * (sim_ctx.omg_i != 0), sim_ctx.jb - 2 * (sim_ctx.omg_j != 0), sim_ctx.dx, sim_ctx.dy);

    fprintf(fp, "POINT_DATA %i \n", (ir - il + 1) * (jt - jb + 1));

    fprintf(fp, "\n");
    fprintf(fp, "VECTORS velocity float\n");
    for (j = jb; j < jt + 1; j++) {
        for (i = il; i < ir + 1; i++) {
            if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i, j) > 9)
                fprintf(fp, "%f %f 0\n", (U(i, j) + U(i, j + 1)) * 0.5, (V(i, j) + V(i + 1, j)) * 0.5);
            else
                fprintf(fp, "0 0 0\n");
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "CELL_DATA %i \n", ((ir - il) * (jt - jb)));
    fprintf(fp, "SCALARS pressure float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (j = jb + 1; j <= jt; j++) {
        for (i = il + 1; i <= ir; i++) {
            if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i, j) > 9)
                fprintf(fp, "%f\n", P(i, j));
            else
                fprintf(fp, "0\n");
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS temperature float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (j = jb + 1; j <= jt; j++) {
        for (i = il + 1; i <= ir; i++) {
            if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i, j) > 9)
                fprintf(fp, "%f\n", T(i, j));
            else
                fprintf(fp, "0\n");
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS density float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (j = jb + 1; j <= jt; j++) {
        for (i = il + 1; i <= ir; i++) {
            if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i, j) > 9)
                fprintf(fp, "%f\n", grid.density(i, j));
            else
                fprintf(fp, "0\n");
        }
    }

    std::vector<matrix2d<double>>& s = grid.species_concentrations;
    for (int c = 0; c < s.size(); ++c) {
        fprintf(fp, "\n");
        fprintf(fp, "SCALARS ");
        fprintf(fp, (grid.species[c].name + " float 1 \n").c_str());
        fprintf(fp, "LOOKUP_TABLE default \n");
        for (j = jb + 1; j <= jt; j++) {
            for (i = il + 1; i <= ir; i++) {
                if (grid.cell_flagg(i, j) == cell_flag::FLUID || grid.cell_flagg(i, j) == cell_flag::INFLOW || grid.cell_flagg(i, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i, j) > 9)
                    fprintf(fp, "%f\n", s[c](i, j));
                else
                    fprintf(fp, "0\n");
            }
        }
    }

    if (fclose(fp))
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}



void write_vtkHeader( FILE *fp, int il, int ir, int jb, int jt,
                      double dx, double dy) {
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Null pointer in write_vtkHeader" );
    ERROR( szBuff );
    return;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated by CFD-lab course output of group D \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"\n");	
  fprintf(fp,"DATASET STRUCTURED_GRID\n");
  fprintf(fp,"DIMENSIONS  %i %i 1 \n", (ir - il)+1, (jt - jb)+1);
  fprintf(fp,"POINTS %i float\n", ((ir - il) +1)*((jt - jb) +1) );
  fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int il, int ir, int jb, int jt, int ibase, int jbase,
                      double dx, double dy) {
  double originX = 0.0;  
  double originY = 0.0;

  int i = 0;
  int j = 0;

  for(j = jb; j < jt+1; j++) {
    for(i = il; i < ir+1; i++) {
      fprintf(fp, "%f %f 0\n", originX+((i+ibase)*dx), originY+((j+jbase)*dy) );
    }
  }
}


void VTKHelper::printVTKFile(GpuGrid grid, double dx, double dy,
        std::string casename, std::string outputdir, int timestep) {

    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Lower Left Corner
    //double x(0), y(0), z(0);
    //for(int col=0;col < grid.jmaxb()-1;col++) {
    //    // Reset x again
    //    x = 0;
    //    for(int row=0;row < grid.imaxb()-1;row++){
    //        points->InsertNextPoint(x, y, z);
    //        x += dx;
    //    }
    //    y += dy;
    //}

    // Specify the dimensions of the grid
    //structuredGrid->SetDimensions(grid.imaxb()-1,grid.jmaxb()-1,1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkDoubleArray* Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray* Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Set number of tuples
    //std::vector<std::vector<Cell>> cells;
    std::vector<std::vector<double>> pressure;
    //grid.pressure(pressure);
    //
    //// Print pressure from bottom to top
    //for(int j = 1; j < grid.jmaxb()-1; j++) {
    //    for(int i = 1; i < grid.imaxb()-1; i++) {
    //        Pressure->InsertNextTuple(&pressure.at(i).at(j));
    //    }
    //}

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Get Velocity
    std::vector<std::vector<double>> velocity_U;
    std::vector<std::vector<double>> velocity_V;
    //grid.velocity(velocity_U, velocity_type::U);
    //grid.velocity(velocity_V, velocity_type::V);
    //
    //// Print Velocity from bottom to top
    //for(int j = 0; j < grid.jmaxb()-1; j++) {
    //        for(int i = 0; i < grid.imaxb()-1; i++) {
    //        vel[0] = (velocity_U.at(i).at(j) + velocity_U.at(i).at(j+1)) * 0.5;
    //        vel[1] = (velocity_V.at(i).at(j) + velocity_V.at(i+1).at(j)) * 0.5;
    //        Velocity->InsertNextTuple(vel);
    //    }
    //}

    // Add Pressure to Structured GpuGrid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured GpuGrid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write GpuGrid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname = outputdir + '/' + casename + "." + std::to_string(timestep) + ".vtk";
    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}