#include "simulator.h"
#include "model.h"

float *rk2(float x[], float y[],float t, std::map<std::string, float> Plist,float* (*slope)(float ,float*,float*, std::map<std::string,float>),float step){
  float k1[50];
  float k2[50];
  float xhalfstep[50];
  
  y = (*slope)(t,x,y,Plist);
  
  for (int i=0; i<25;i++){
    k1[i] = y[i];
    xhalfstep[i] = x[i] + step*k1[i];
  }

  y = (*slope)(t+step,xhalfstep,y,Plist);

  for (int i=0; i<25;i++){
    k2[i] = y[i];
    x[i] = x[i] + step*(k1[i] + k2[i])*0.5;
  }
  return x;
  
}



float *rk4(float x[], float y[],float t, std::map<std::string, float> Plist,float* (*slope)(float ,float*,float*, std::map<std::string,float>),float step){
  float k1[50];
  float k2[50];
  float k3[50];
  float k4[50];
  float xhalfstep[50];
  float xfullstep[50];

  y = (*slope)(t,x,y,Plist);
  
  for (int i=0; i<25;i++){
    k1[i] = y[i];
    xhalfstep[i] = x[i]+0.5*step*k1[i];
  }
  
  y = (*slope)(t+step/2,xhalfstep,y,Plist);
  
  for (int i=0; i<25;i++){
    k2[i] = y[i];
    xhalfstep[i] = x[i]+0.5*step*k2[i];
  }

  y = (*slope)(t + step/2,xhalfstep,y,Plist);
  for (int i=0; i<25;i++){
    k3[i] = y[i];
    xfullstep[i] = x[i]+ k3[i]*step;
  }

  y = (*slope)(t + step,xfullstep,y,Plist);
  
  for (int i=0; i<25;i++){
    k4[i] = y[i];
  }
 

  for (int i=0; i<25;i++){
    x[i] = x[i] + 1.0/6.0*(k1[i]+ 2*k2[i] + 2*k3[i] + k4[i])*step;
  }

  return x;}


float *euler(float x[], float y[],float t, std::map<std::string, float> Plist,float* (*slope)(float ,float*,float*, std::map<std::string,float>),float step){

  y = (*slope)(t,x,y,Plist);
  for (int i=0; i<25;i++){
    x[i] = y[i]*step + x[i];
  }

  return x;    
  
}

int ode(float x[], float tstart, float tend, float y[], float step, 
        float* (*function) (float,float*,float*,std::map<std::string,float>),
        std::string solver,
        bool verb,
        bool outflag, std::string outname,
        int argc=0, char** argv= (char) 0){

  // Defined in model.cpp
  // Generated by pydstool2cpp.py because
  // var and par list might be model dependent
  std::map<std::string,float> Plist;
  std::map<std::string,float> Vlist;
  
  Plist = initializeParamMap(Plist, argc, argv, verb);
  

  Vlist = initializeICSMap(Vlist, argc, argv, verb);
  
  int counter = 0;
  
  for (int i=0; i<26;i++){
    x[counter] = Vlist[names[i]];
    counter++;}
  
  std::ofstream values;
  std::string defaultoutname = "./values.dat";
  if (outflag == true){
    defaultoutname = outname;}
  if (verb == true){
    std::cout<<"\nOutput file name is "<<defaultoutname;}
  
  values.open(defaultoutname.c_str());

  values<<"t";
  for (int i=0; i<26; i++){
    values<<"\t"<<names[i];}
  values<<"\n";
  
  float tcurr = tstart;

  values<<tcurr;
  for (int j=0;j<26;j++){
    values<<"\t"<<x[j];
  }
  
  for (int i=0;i<tend/step;i++){
    if (solver=="euler"){
      x=euler(x,y, tcurr, Plist, function,step);}
    if (solver=="rk4"){
      x=rk4(x,y, tcurr, Plist, function,step);}
    if (solver=="rk2"){
      x=rk2(x,y, tcurr, Plist, function,step);}
    values<<"\n";
    
    tcurr = tcurr + step;
    
    values<<tcurr;
    for (int j=0;j<26;j++){
      values<<"\t"<<x[j];
    }
    

  }
  values.close();
  
  return 0;
}
