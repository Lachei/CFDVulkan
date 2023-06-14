//
// Created by Moritz Gnisia on 08.02.20.
//

#include "FileIterator.hpp"
#define logging


FileIterator::FileIterator(dir_iter& dir, std::string& prefix, bool merge) {
    // Definition of lambda function to check weather the prefix is in the file or not
    auto func = [prefix](std::filesystem::path i){
        return !(i.filename().string().find(prefix) == std::string::npos);
};
    // copy the directory content to the directory vector dir_vec and ensure that the prefix is in the filename
    std::copy_if(std::filesystem::begin(dir), std::filesystem::end(dir), std::back_inserter(dir_vec), func);
    // sort the vector
    std::sort(dir_vec.begin(), dir_vec.end());
    // Definition of index variable
    std::size_t index;
    std::regex file_structure;
    // Iterate over all files
    for (vector_path::iterator it(dir_vec.begin()), it_end(dir_vec.end()); it != it_end; ++it) {
        // check that the filename consists out of three elements: prefix, timestep, vtk
        if (merge) {
            // Filename cavity_3.9998.vtk
            file_structure = std::regex("[a-zA-Z]+_\\d.\\d+.vtk");
        }
        else {
            // Filename cavity.9998.vtk
            file_structure = std::regex("[a-zA-Z]+\\.\\d+\\.vtk");
        }
        if (!std::regex_match(it->filename().string(), file_structure)) {
            std::cout << it->filename().string() << std::endl;
            std::cout << "did not match file type" << std::endl;
        }
        index++;
    }
}



bool FileIterator::compare(vector_path& vec_compare) {
    for (std::pair<vector_path::iterator, vector_path::iterator> files(FileIterator::dir_vec.begin(), vec_compare.begin());
     files.first != FileIterator::dir_vec.end(); ++files.first, ++files.second) {
        if (files.first->filename() != files.second->filename()) {
            // return false if the ordering of the files for the directories is not identical
            return false;
        }
    }
    // return true if the ordering of the files for the directories is identical
    return true;
}


vector_path* FileIterator::get_dir_vec() {
    return &dir_vec;
}
