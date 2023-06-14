
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkBoundingBox.h>
#include <vtkCellLocator.h>
#include <vtkNew.h>
#include <vtkCell.h>
#include <vtkTuple.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkRenderWindow.h>
#include <vtkCellData.h>
#include <vtkRenderer.h>
#include <vtkDataArray.h>
#include <vtkObject.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTypedDataArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPoints.h>
#include <cmath>
#include "FileIterator.hpp"
#include <iterator>     // std::advance
#include <algorithm>
#include "filesystem"
#include <regex>
#include <set>
#include "cxxopts.hpp"
#include <fmt/format.h>
#include "helper.hpp"

int main(int argc, char** argv) {

    cxxopts::Options options("Merge", "Merging vtk files for different subdomains back into one main vtk files");

    std::string prefix;
    std::string suffix;
    std::string origin;
    std::string dest;
    std::regex regex_folder("[\\S]+\\/");
    std::regex regex_suffix("[a-zA-Z]+");
    
    options.add_options()
        ("o,origin", "Origin Folder of the vtk files e.g. folder/", cxxopts::value<std::string>())
        ("d,destination", "Destination Folder, if no destination folder is given the files will be written into the root of folder of this executable e.g. dest/", cxxopts::value<std::string>()->default_value(""))
        ("p,prefix", "prefix of the vtk files e.g. cavity_", cxxopts::value<std::string>()->default_value("cavity_"))
        ("s,suffix", "prefix of the vtk files e.g. .vtk, default value vtk", cxxopts::value<std::string>()->default_value("vtk"))
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);
    std::string source;

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        std::cout << print_nomenclature() << std::endl;
        exit(0);
    }
    if (!result.count("origin")) {
        std::cout << "Please provide the origin location of the subdomain vtk files" << std::endl;
        std::cout << "run ./VTK_Merge -help for more information" << std::endl;
        exit(1);
    }
    else {
        // Test if the folder contains a / in the end
        bool regex_match = std::regex_match(result["origin"].as<std::string>(), regex_folder);
        if (!regex_match) {
            std::cout << "The provided path contains either whitespaces or doesn't have a trailing slash" << std::endl;
            exit(1);
        }
        else {
            origin = result["origin"].as<std::string>();
        }
    dest = result["destination"].as<std::string>();
    // Create Directory if this doesn't exists
    if (dest.length() > 0) {
        if (!std::filesystem::is_directory(dest)) {
            std::filesystem::create_directory(dest);
        }
    }

    prefix = result["prefix"].as<std::string>();
    suffix = result["suffix"].as<std::string>();
    if (!std::regex_match(suffix, regex_suffix)) {
        std::cout << "The provided suffix contains a non alphabetic value e.g.: .;/;" << std::endl;
        exit(1);
    }
    }

    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // File Reader
    vtkSmartPointer<vtkStructuredGridReader> reader = vtkSmartPointer<vtkStructuredGridReader>::New();

    // Variables for finding closest cell and points
    vtkSmartPointer<vtkCellLocator> cellLocator =  vtkSmartPointer<vtkCellLocator>::New();
    vtkSmartPointer<vtkKdTreePointLocator> PointLocator =  vtkSmartPointer<vtkKdTreePointLocator>::New();
    double center[3];
    double closestPoint[3];//the coordinates of the closest point will be returned here
    double closestPointDist2; //the squared distance to the closest point will be returned here
    vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
    int subId; //this is rarely used (in triangle strips only, I believe)

    // Directory Iterator
    std::filesystem::directory_iterator true_sol_iter(origin);

    // Create File Iterator
    FileIterator file_iter_true(true_sol_iter, prefix, true);
    vector_path* vec_true = file_iter_true.get_dir_vec();

    // Check if files exists
    std::string file = (*vec_true)[0];
    std::cout << "Origin " << origin << std::endl;
    std::cout << "prefix " << prefix << std::endl;
    std::cout << "suffix " << suffix << std::endl;
    // Example filename: merged/output_0.0.vtk
    /* Formating Overview   
    merged/output_0.0.vtk   Table       
    <-----><----->| | {4}   {0} Origin  
        {0}    {1}  |       {1} Prefix  
                {2} {3}     {2} Rank    
                            {3} Timestep
                            {4} Suffix  
    */
    std::string file_fmt = fmt::format("{0}{1}{2}.{3}.{4}", origin, prefix, 0, 0, suffix);
    if (file.compare(file_fmt) != 0) {
        std::cout << "The first file in the target directory is: " << file << std::endl;
        std::cout << "The formated file string is: " << file_fmt << std::endl;
        std::cout << "Please change the formating in Line " << __LINE__ - 4 << std::endl;
        exit(1);
    }   
    // Get Unique number of timesteps and ranks
    std::vector<int> timesteps = get_unique(vec_true, get_timestamp);
    std::vector<int> ranks = get_unique(vec_true, get_rank);

    // Create Sample Batch
    std::vector<std::string> files;
    for (int rank: ranks) {
        files.push_back(fmt::format("{0}{1}{2}.{3}.{4}", origin, prefix, rank, timesteps[0], suffix));
    }


    // Init Global Dimension, (!) Assumption Dim z = 0
    int dim_x(0), dim_y(0);

    // Iterate over files to and dimensions, Assumption: Merge two files in x and y direction
    for (int i = 0; i < 2; ++i) {

        // Read files
        reader->SetFileName(files[i].c_str());
        reader->ReadAllVectorsOn();
        reader->ReadAllScalarsOn();
        reader->Update();

        // local pressure
        vtkStructuredGrid* grid_local = reader->GetOutput();

        // Print Extent
        dim_x  += grid_local -> GetExtent()[1];
        dim_y  += grid_local -> GetExtent()[3];
    }

    std::cout << "Dim x: " << dim_x << std::endl;
    std::cout << "Dim y: " << dim_y << std::endl;

    // Get Spacing
    reader->SetFileName(files[0].c_str());
    reader->ReadAllVectorsOn();
    reader->ReadAllScalarsOn();
    reader->Update();

    double* bounds;
    float spacing_x(0), spacing_y(0);
    // Bounds are defined here: https://vtk.org/doc/nightly/html/classvtkDataSet.html#a6a25018c1af03c2fc1068ea4a0221281
    // Structure: (xmin,xmax, ymin,ymax, zmin,zmax)
    bounds = reader->GetOutput()->GetCell(0) -> GetBounds();
    // Spacing x
    spacing_x =  bounds[1] - bounds[0];
    // Spacing y
    spacing_y =  bounds[3] - bounds[2];

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    float x(0.0), y(0.0), z(0.0);
    for(size_t j = 0; j < dim_y+1; j++) {
        x = 0;
        for(size_t i = 0; i < dim_x+1; i++) {
            points->InsertNextPoint(x, y, z);
            x += spacing_x;
        }
        y += 0.02;
    }

    
    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(dim_x+1,dim_x+1,1);
    structuredGrid->SetPoints(points);

    // new global pressure
    vtkFloatArray* Pressure = vtkFloatArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Set number of tuples
    Pressure->SetNumberOfTuples(structuredGrid -> GetNumberOfCells());
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Velocity
    vtkFloatArray* Velocity = vtkFloatArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);
    // Set Number of Tuples
    Velocity -> SetNumberOfTuples(structuredGrid -> GetNumberOfPoints());
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Locator for larger Grid
    // Celllocator is required for the pressure values
    cellLocator->SetDataSet(structuredGrid);
    cellLocator->BuildLocator();
    // Celllocator is required for the velocity values
    PointLocator->SetDataSet(structuredGrid);
    PointLocator->BuildLocator();
    double coords[3];
    // std::vector<int> subset_timesteps = subset(timesteps, 40);
    for (int step: timesteps) {
        std::cout << "Merging TS: " << step << std::endl;
        files.erase(files.begin(), files.end());
        for (int rank: ranks) {
            files.push_back(fmt::format("{0}{1}{2}.{3}.{4}", origin, prefix, rank, timesteps[0], suffix));
        }
        
        for(std::string str: files) {
            reader->SetFileName(str.c_str());
            reader->ReadAllVectorsOn();
            reader->ReadAllScalarsOn();
            reader->Update();

            // local pressure
            vtkStructuredGrid* grid_local = reader->GetOutput();
            int cells = grid_local -> GetNumberOfCells();
            int points = grid_local -> GetNumberOfPoints();

            for (int k = 0; k < cells; ++k) {
                // using the bbox of every cell to find the corresponding cell id in the larger grid
                vtkBoundingBox bbox_test(grid_local -> GetCell(k) -> GetBounds());
                bbox_test.GetCenter(center);
                cellLocator->FindClosestPoint(center, closestPoint, cellId, subId, closestPointDist2);
                Pressure -> InsertTuple(cellId , grid_local -> GetCellData() -> GetScalars("pressure") -> GetTuple(k));
            }
            for (int k = 0; k < points; ++k) {
                grid_local -> GetPoint(k, coords);
                // for the velocity we need the pointlocator
                Velocity -> InsertTuple(PointLocator->FindClosestPoint(coords) , grid_local -> GetPointData() -> GetVectors("velocity") -> GetTuple(k));
            }

        }

        // New writer
        vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
        std::string output = fmt::format("{0}{1}{2}.{3}", dest, prefix, step, suffix);
        writer->SetFileName(output.c_str());
        writer->SetInputData(structuredGrid);
        // writer->SetInputData(dataset);
        writer->Write();

    }
    
}
