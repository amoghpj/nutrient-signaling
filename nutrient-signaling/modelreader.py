'''
Author: Amogh Jalihal
Date: 2018-02-14
'''

'''
TODO Add parameter to main to get simulation duration from user from command line
'''

def checkstring(inputfile,line):
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

def readinput(PATH):

    modeldefinition={
    'variables':{},
    'parameters':{},
    'initialconditions':{}
    }
    
    for inputfile in modeldefinition.keys():
        print("Reading " + inputfile)
        with open(PATH+'/'+inputfile+".txt",'r') as infile:
            for line in infile.readlines():
                if checkstring(inputfile,line)==1:
                    continue
                else:
                    key,value=checkstring(inputfile,line)
                    modeldefinition[inputfile][key]=value
    return(modeldefinition)

def main(inputpath):
    return(readinput(inputpath))

if __name__=='__main__':

    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option('-p','--path',dest='path',default=False,
                  help='Path to model definition')
    (options,args)=parser.parse_args()
    
    main(options.path)
