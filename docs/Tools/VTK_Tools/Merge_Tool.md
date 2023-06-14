## VTK Merge Tool

This small CLI tool allows you to merge the VTK files from different subdomains into one main vtk file for different timesteps.

### Requirements

- gcc9
- vtk 7.1.1
- fmt
- cmake

### CLI Arguments

The following CLI Arguments can be used:

This also visible by running `./Merge -h`

```shell
Merging vtk files for different subdomains back into one main vtk files
Usage:
  Merge [OPTION...]

  -o, --origin arg       Origin Folder of the vtk files e.g. folder/
  -d, --destination arg  Destination Folder, if no destination folder is
                         given the files will be written into the root of folder
                         of this executable (default: "")
  -p, --prefix arg       prefix of the vtk files (default: cavity_)
  -s, --suffix arg       prefix of the vtk files (default: vtk)
  -h, --help             Print usage
```

### Example Run

```
./VTK_Merge -o ../Misc/VTK_Subdomains/ -d ../Misc/Merged/ -p output_ 
```

### Tests

You can run the tests with `ctest` in the build directory. Currently they only check different prefixes and timesteps.