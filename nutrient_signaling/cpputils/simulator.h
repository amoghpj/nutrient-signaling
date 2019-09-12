#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <string>

extern int ode(float*, float, float, float*,
               float,
               float* (*function) (float,float*,float*,std::map<std::string,float>),
               std::string,
               bool, bool , std::string,
               int, char**);
#endif
  


