Model
=========

folder for the model code in C

## Installation of C code
Download the code as a zip or clone (or checkout) the repository.

## Compilation of C code (using command-line)
Guide to compile the files and generate the executable.
  1. Change into directory `phylosim/model/`
  2. Execute the  command `g++ -Wall -std=c++11 -c -o File.o File.cpp` for every **File** in the directory 
  3. Compile the files in the directory *Utils* by executing the command `g++ -Wall -std=c++11 -c -o File.o File.cpp` for every **File** in the directory *Utils*
  4. (Compile third party software) Compile the file from pugi-xml using the commands above e.g. 
    1. Change directory to `phylosim/model/XML/ThirdPartySoftware/PugiXML/src/`
    2. Execute the  command `g++ -Wall -std=c++11 -c -o pugixml.o pugixml.cpp`
  5. To generate the executable file you execute the command `g++ -Wall -std=c++11 -o PhyloSim.exe File.o ` where **File** must again be replaced by the name of the every file from the directory  `phylosim/model/`,  `phylosim/model/Utils/` and `phylosim/model/XML/ThirdPartySoftware/PugiXML/src/`
