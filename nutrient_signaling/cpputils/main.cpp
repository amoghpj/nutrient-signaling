#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <cmath>
#include <map>
#include <time.h>
#include <string>
#include "model.h"
//#include "simulator.h"

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;

static void show_help(std::string name){
 std:cerr<<"Usage:"<<name<<" <option(s)> "
         <<"\nOptions:\n"
         <<"\t--help\tShow this help message\n"
         <<"\t--ode\tRun the ODE solver. If no other\n\t\t options are passed, default parameters are used\n\tIf ODE, additional options are:\n"
         <<"\t\t--step\tset step size for solver [Default = 0.01]\n"
         <<"\t\t--tend\tset tend for solver [Default = 90]\n"
         <<"\t\t--pars\tSpecify parameter name SPACE value\n"
         <<"\t\t--ics\tSpecify variable name SPACE initial value\n"
         <<"\t\t--solver\tSpecify solver [rk4,EULER]\n"
         <<"\tMiscellaneous:\n"
         <<"\t--plot\tplot results using plotfile.py [Default = FALSE]\n"
         <<"\t--verb\tVerbose output [Default = FALSE]\n";
}


// Original observer
// struct push_back_state_and_time{
//     std::vector< state_type >& m_states;
//     std::vector< double >& m_times;

//     push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
//     : m_states( states ) , m_times( times ) { }

//     void operator()( const state_type &x , double t )
//     {
//         m_states.push_back( x );
//         m_times.push_back( t );
//     }
// };

int main(int argc, char** argv){
  
  
  time_t timerstart;
  time_t timerstop;
  bool ODEflag = false;
  //float ics[50];
  
  bool userDef_step = false;
  double userDef_stepval = 0.0;
  
  bool userDef_tend = false;
  double userDef_tendval = 0.0;

  bool userDef_out = false;
  std::string userDef_outval;
  
  bool verb = false;

  bool plotresults = false;
  // std::map<string,float> Plist;
  // std::map<string,float> Vlist;

  std::string solver="euler";

  
  if (argc==1){
    show_help(argv[0]);
    return 1;}

  else{   
    for (int i=1;i<argc;i++){
      std::string arg = argv[i];

      if (arg == "--ode"){
        ODEflag = true;
      }
      if (arg == "--step"){
        if (i+1 < argc){
          userDef_step = true;
          userDef_stepval = atof(argv[i+1]);
        }
        else{
          cout<<"Please pass value of step\n";
          return 1;}
        }
      if (arg == "--solver"){
        if (i+1 < argc){
          solver = argv[i+1];
        }

        if ((solver != "euler") && (solver != "rk2") &&  (solver != "rk4")){
          cout<<"\nINVALID SOLVER "<<solver<<"\n";
          return 1;}
        
      }
      if (arg == "--tend"){
        if (i+1 < argc){
          userDef_tend = true;
          userDef_tendval = atof(argv[i+1]);
        }
        else{
          cout<<"Please pass value of tend\n";
          return 1;}
        }

      if (arg == "--out"){
        if (i+1 < argc){
          userDef_out = true;
          userDef_outval = argv[i+1];
        }
        else{
          cout<<"Please pass value of tend\n";
          return 1;}
        }
      
      
      if (arg == "--plot"){
        plotresults = true;}
      
      if (arg == "--verb"){
        verb = true;}

      
    }
  }
  
  if (ODEflag == true ){
    //   timerstart = time(NULL);
    //   float x[50];
    //   float y[50];
    double tend;
    tend = 90.0;
    //   float tstart = 0.0;
    double step = 0.01;
  //   float r=1.0;
    
    if (userDef_step == true){
      step = userDef_stepval;}
    
    if (userDef_tend == true){
      tend = userDef_tendval;}
    //   if (verb ==true){
    //     cout<<"\nRunning ODE solver!";}
    //   // function call
    
    //   ode(x, tstart, tend,y, step, NutSig, solver, verb, userDef_out, userDef_outval, argc, argv);
    
    //   if (verb ==true){
    //     cout<<"\nDone!\n";}
    //   }
    // timerstop = time(NULL);
    // double seconds;
    // seconds = difftime(timerstop,timerstart);
    


    state_type x(25);
    std::map<std::string, float> Vlist;
    
    bool verb;
    verb = true;
    
    Vlist = initializeICSMap(Vlist,argc, argv, verb);
    for (int i = 0; i<25; i++){
      x[i]=Vlist[names[i]];
      //cout<<x[i];
    }
    vector<state_type> x_vec;
    vector<double> times;
    
    std::map<std::string, float> Plist;
    Plist = initializeParamMap(Plist,argc, argv, verb);
    
    NutSig nsm(Plist);

    std::ofstream outvals;
    std::string outfilename;
    
     outfilename = "values.dat";
     if (userDef_out = true){
       outfilename = userDef_outval;}
     
    ////////////////////////////////////////////////////////    
    //size_t steps = integrate(nsm, x, 0.0, tend, step, push_back_state_and_time( x_vec , times ) );


     ////////////////////////////////////////
     // Use a stepper:
     size_t imax = 100000;
     size_t i = 0;
     
     auto stepper = make_dense_output(1.0e-6, 1.0e-6, runge_kutta_dopri5< state_type>());
     
     stepper.initialize(x,0,step);
     
     ////////////////////////////////////////
     outvals.open(outfilename.c_str());
     
     outvals<<"t"; 
     for (int j=0;j<25;j++){
       outvals<<"\t"<<names[j];}
     outvals<<"\n";
     
     while((stepper.current_time() < tend) && (i < imax)){
       outvals<<stepper.current_time()<<'\t';
       for (int j=0;j<25;j++){
         outvals<<stepper.current_state()[j]<<"\t";}
       outvals<<"\n";
       stepper.do_step(nsm);
       ++i;
     }
     ///////////////////////////////////////
     // Old
     // for (int i=0;i<steps;i++){
     //   outvals<<times[i]<<"\t";
     //   for (int j=0;j<25;j++){
     //     outvals<<x_vec[i][j]<<"\t";}
     //   outvals<<"\n";
     // }
    ////////////////////////////////////////////////////////
    // RK4!!
     // outvals.open(outfilename.c_str());
     
     // runge_kutta4< state_type > stepper;
     // integrate_const(stepper,nsm, x, 0.0, tend, step);

     // outvals<<"t"; 
     // for (int j=0;j<25;j++){
     //   outvals<<"\t"<<names[j];}
     // outvals<<"\n";
     
     // const double dt = step;
     // for (double t=0.0;t<tend;t+=dt){
     //   stepper.do_step(nsm,x,t,dt);
       
     //   outvals<<t<<"\t";
       
     //   for (int j=1;j<26;j++){
     //     outvals<<x[j]<<"\t";
         
     //   }
     //   outvals<<"\n";
     // }
     ////////////////////////////////////////////////////////

     

    outvals.close();
    if (plotresults == true){
      system("python plotfile.py");}
    // if (verb == true){
    //   cout<<"This took "<<seconds<<" s for 1 simulation\n";}
    
    return 0;
  }
}

