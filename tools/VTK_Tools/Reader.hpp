//
// Created by Moritz Gnisia on 30.01.20.
//

#ifndef STRUCTUREDGRIDREADER_READER_H
#define STRUCTUREDGRIDREADER_READER_H
#include <vtkStructuredGridReader.h>
#include <filesystem>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
class Reader {        // The class
public:          // Access specifier
    std::filesystem::path file_vtk;  // Attribute
    vtkSmartPointer<vtkStructuredGrid> sg;
    vtkSmartPointer<vtkPointData> pd;
    vtkSmartPointer<vtkCellData> cd;
    Reader(std::filesystem::path file);
};


#endif //STRUCTUREDGRIDREADER_READER_H
