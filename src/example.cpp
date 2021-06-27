#include <iostream>
#include <cstring>

#include "four_russians.h"

// simple program showing how to use the implementation

int main(int argc, char ** argv)
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " string1 string2" << std::endl;
    return -1;
  }
  
  // parameters small enough to facilitate fast response
  const int tm = 2, tn = 3;
  
  // set the parameters tm and tn by creating a new object of the class Four_Russians
  FR::Four_russians<tm, tn> four_russians;
  // run the preprocessing step
  four_russians.preprocessing_step();
  
  // and run the computation step by supplying the input strings and their sizes
  long result = four_russians.computation_step(
    argv[1], std::strlen(argv[1]), argv[2], std::strlen(argv[2]));
  
  std::cout << "The edit distance between \"" << argv[1] << "\" and \""
    << argv[2] << "\" is " << result << ".\n";
}
