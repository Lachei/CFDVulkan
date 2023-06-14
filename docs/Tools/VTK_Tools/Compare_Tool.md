# VTK Compare Tool

## Requirements

- gcc9
- vtk 7.1.1
- cmake

## Building


1. `cd tools`
2. `mkdir build && cd build`
3. `cmake .. && make`

Oneliner:

```
cd tools && mkdir build && cd build && cmake .. && make
```

## Usage

```shell
Comparing two folders of VTK Files
Usage:
  VTK Compare Tool [OPTION...]

  -o, --origin arg   Origin folder name of the vtk files e.g. folder/; this
                     is the reference solution
  -c, --current arg  Folder name which contains the vtk files that should be
                     compared e.g. dest/
  -f, --first arg    prefix of the vtk files in the origin (first) given
                     folder e.g. cavity
  -s, --second arg   prefix of the vtk files in the current (second) given
                     folder e.g. cavity
  -h, --help         Print usage
```

## Example Usage

```shell
# executed within the build directory
./VTK_Comp -o ../Misc/VTK_Origin/ -c ../Misc/VTK_Compare/ -f output -s output
```

!!! warning
    This example is intended to fail