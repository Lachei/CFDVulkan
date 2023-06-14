//
// Created by Moritz Gnisia on 08.02.20.
//

#ifndef VTKREADER_FILEITERATOR_H
#define VTKREADER_FILEITERATOR_H

#include <vector>
#include <filesystem>
#include <iostream>
#include <regex>
typedef std::vector<std::filesystem::path> vector_path;
typedef std::filesystem::directory_iterator dir_iter;

class FileIterator {

private:
    dir_iter dir;
    vector_path  dir_vec;

public:
    FileIterator(dir_iter& dir, std::string& prefix, bool merge);
    bool compare(vector_path& vec_compare);
    vector_path* get_dir_vec();
};


#endif //VTKREADER_FILEITERATOR_H
