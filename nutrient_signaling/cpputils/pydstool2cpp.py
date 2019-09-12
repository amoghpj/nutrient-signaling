import os
import sys

class PyDSTool2CPP:
    def __init__(self, modelpath):
        self.modelpath = modelpath
        self.writepath = './src/'
        if not os.path.isdir(self.modelpath):
            print("Path to model defintion must be a folder")
            sys.exit()

            
    def setwritepath(self, writepath):
        self.writepath = writepath

    def getwritepath(self):
        return(self.writepath)    
        
    def checkstring(self, inputfile,line):
        if line[0]=='#':
            return 1
        elif line[0]=='\n':
            return 1
        else:
            K,V=line.split('\n')[0].split('\t')
            if '#' in V:
                V=V.split('#')[0]
            if inputfile!='variables':
                return((K,float(V)))
            else:
                return((K,V))
    
    def readinput(self):
        modeldefinition={
            'variables':{},
            'parameters':{},
            'initialconditions':{}
        }
        
        for inputfile in modeldefinition.keys():
            print("Reading " + inputfile)
            with open(self.modelpath +'/'+inputfile+".txt",'r') as infile:
                for line in infile.readlines():
                    if self.checkstring(inputfile,line) == 1:
                        continue
                    else:
                        key,value = self.checkstring(inputfile,
                                                     line)
                        modeldefinition[inputfile][key] = value
        return(modeldefinition)

    def writecpp(self):
        if not os.path.exists(self.writepath):
            os.mkdir(self.writepath)        
        ModelDef = self.readinput()
        with open(self.writepath + "model.h",'w') as outfile:
            outfile.write("#ifndef MODEL_H\n")
            outfile.write("#define MODEL_H\n\n")
            outfile.write("#include <iostream>\n")
            outfile.write("#include <cmath>\n")
            outfile.write("#include <map>\n")
            outfile.write("#include <string>\n")
            outfile.write("#include <algorithm>\n")
            outfile.write("#include <vector>\n")
        
            outfile.write("extern std::map<std::string,float> initializeParamMap(std::map<std::string,float> Plist,int,char**, bool);\n")
            outfile.write("extern std::map<std::string,float> initializeICSMap(std::map<std::string,float> Vlist,int ,char** , bool);\n")
        
            outfile.write("extern char names[25][25];\n")
            outfile.write("typedef std::vector <double> state_type;\n")
            ###########################################################################
            #### Class Nutsig
            outfile.write("class NutSig{\n")
            outfile.write("\n//Parameters\n")
            outfile.write("float ")
            i = 0 
            for parameter in ModelDef['parameters'].keys():
                if parameter != 'index':
                    outfile.write(parameter)
                    i +=1
                    if i < len(ModelDef['parameters'].keys()):
                        outfile.write(",\n")
                    else:
                        outfile.write(";\n\n")
            
        
            outfile.write("public:\n")
        
            outfile.write("NutSig(std::map<std::string, float> Plist){\n")
            s = ''
            for parameter in ModelDef['parameters'].keys():
                if parameter != 'index':        
                    s = parameter + " = Plist[\"" + parameter+"\"]"
                    i += 1
                    if i < len(ModelDef['parameters'].keys()):
                        outfile.write(s + ",\n")
                    else:
                        outfile.write(s + ";\n")
            
            outfile.write("}\n")
            outfile.write("float shs(float sigma, float omega){\nreturn 1/(1+exp(-sigma*omega));}\n\n")
            outfile.write("float tRNA(float tRNA_tot,float AmAc){\nreturn std::min(tRNA_tot, AmAc);}\n\n")
            outfile.write("float pRib(float rib_comp, float init_factor){\nreturn std::min(rib_comp,init_factor);}\n\n")
        
            outfile.write("void operator() (const state_type &x, state_type &xdot, const double ){\n")
            
            outfile.write("\n//Variables\n")
            outfile.write("float ")
            i = 0
            s = ''
            for variable in ModelDef['variables'].keys():
                s = variable + " = " + "x[" + str(i) +"]"
                i += 1
                if i < len(ModelDef['variables'].keys()):
                    outfile.write(s + ",\n")
                else:
                    outfile.write(s + ";\n")
        
            i = 0
            outfile.write("\n//Equations\n")
            for variable in ModelDef['variables'].keys():
                outfile.write("xdot[" +str(i) + "] = " +str(ModelDef['variables'][variable]).replace('*min','*std::min')+ ";\n")
                i += 1
        
            outfile.write("}};\n")
        
            ###########################################################################
            outfile.write("#endif")
            
        with open(self.writepath + "model.cpp",'w') as outfile:
            outfile.write("#include \"model.h\"\n")
        
            outfile.write("char names[25][25] = {")
            i=0
            for variable in ModelDef['variables'].keys():
                outfile.write("\"" + str(variable) + "\"")
                i+=1
                if i<len(ModelDef['variables'].keys()):
                    outfile.write(',')
            outfile.write("};\n\n")
            
            ## initializeParamsMap
            
            outfile.write("std::map<std::string,float> initializeParamMap(std::map<std::string,float> Plist, int argc,char** argv, bool verb){\n")
            outfile.write("bool paramflag = false;\n")
            outfile.write("bool icsflag = false;\n")
            outfile.write("char ListOfPars[200][25] = {\n")
            i=0
            for parameter in ModelDef['parameters'].keys():
                outfile.write("\"" + parameter + "\"")
                if i < len(ModelDef['parameters'].keys()):
                    outfile.write(",\n")
                i += 1
        
            outfile.write("};\n\n")
                
            for parameter in ModelDef['parameters'].keys():
                s = "Plist[\"" + parameter+ "\"] = " + str(ModelDef['parameters'][parameter]) +";\n"
        
                outfile.write(s)
            outfile.write(" if (verb == true){\n")
            outfile.write("   std::cout<<\"\\nIn initializeParamMap()\";")
            outfile.write("   std::cout<<\"\\nWe have \"<<argc<<\" arguments\\n\";}\n")
            outfile.write(" \n")
            outfile.write(" for (int i=1;i<argc;i++){\n")
            outfile.write("    std::string arg = argv[i];\n")
            outfile.write("    if (arg == \"--pars\"){\n")
            outfile.write("      paramflag = true;\n")
            outfile.write("      icsflag = false;}\n")
            outfile.write("    if (i+1 < argc){\n")
            outfile.write("      if ((paramflag==true) && (icsflag==false)){\n")
            outfile.write("        arg = argv[i];\n")
            outfile.write("        for (int j = 0; j<200;j++){\n")
            outfile.write("          if (arg == ListOfPars[j]){\n")
            outfile.write("            if (verb ==true){\n")
            outfile.write("            std::cout<<arg<<\" found in ListofPars!\\n\";\n")
            outfile.write("            std::cout<<\"It will be assigned the value=\"<<atof(argv[i+1])<<\"\\n\";}\n")
            outfile.write("            Plist[arg] = atof(argv[i+1]);}\n")
            outfile.write("        }}}}\n")
            outfile.write("return Plist;}\n")
            outfile.write("\n")
            # ## initializeICSMap
            #outfile.write("std::map<string,float> Vlist;\n")
        
            outfile.write("std::map<std::string,float> initializeICSMap(std::map<std::string,float> Vlist, int argc,char** argv, bool verb){\n")
            outfile.write("  bool paramflag =false;\n")
            outfile.write("  bool icsflag =false;\n")
            outfile.write("  char ListOfVars[200][25] = {\n")
            i=0
            for var in ModelDef['variables'].keys():
                outfile.write("\"" + var + "\"")
                if i < len(ModelDef['variables'].keys()):
                    outfile.write(",\n")
                i += 1
            outfile.write("};\n")
        
        
            outfile.write("\n")
                
            for variable in ModelDef['variables'].keys():
                s = "Vlist[\"" + variable +"\"] = " + str(ModelDef['initialconditions'][variable]) +";\n"
                outfile.write(s)
            outfile.write(" if (verb ==true){\n")
            outfile.write("   std::cout<<\"In initializeICSMap()\";")
            outfile.write("   std::cout<<\"\\nWe have \"<<argc<<\" arguments\";}\n")
            outfile.write("\n")
            outfile.write(" \n")
            outfile.write(" for (int i=1;i<argc;i++){\n")
            outfile.write("    std::string arg = argv[i];\n")
            outfile.write("    if (arg == \"--ics\"){\n")
            outfile.write("      paramflag = false;\n")
            outfile.write("      icsflag = true;}\n")
            outfile.write("    if (i+1 < argc){\n")
            outfile.write("      if ((paramflag==false) && (icsflag==true)){\n")
            outfile.write("        arg = argv[i];\n")
            outfile.write("\n")
            outfile.write("        for (int j = 0; j<200;j++){\n")
            outfile.write("          if (arg == ListOfVars[j]){\n")
            outfile.write("\n")
            outfile.write("            if (verb == true){\n")
            outfile.write("            std::cout<<arg<<\" found in ListofVars!\\n\";\n")
            outfile.write("            std::cout<<\"It will be assigned the value=\"<<atof(argv[i+1])<<\"\\n\";}\n")
            outfile.write("\n")
            outfile.write("            Vlist[arg] = atof(argv[i+1]);\n")
            outfile.write("\n")
            outfile.write("          }\n")
            outfile.write("          \n")
            outfile.write("        }}}}\n")
            outfile.write("    return Vlist;\n")
            outfile.write("}\n")
