//
// Created by Moritz Gnisia on 30.01.20.
//

#include "Reader.hpp"
// define LOG

Reader::Reader(std::filesystem::path file) {
    file_vtk = file;
    vtkSmartPointer<vtkStructuredGridReader> reader = vtkSmartPointer<vtkStructuredGridReader>::New();
    #ifdef LOG
        std::cout << "Reading: " << file.c_str() << std::endl;
    #endif
    reader->SetFileName(file.c_str());
    reader->SetReadAllScalars(true);
    reader->SetReadAllVectors(true);
    reader->Update();
    sg = reader -> GetOutput();
    pd = sg -> GetPointData();
    cd = sg -> GetCellData();
}
