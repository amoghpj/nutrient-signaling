# Author: Amogh Jalihal
import os
import sys
import pandas as pd
import numpy as np
from pyDOE import lhs
import copy
import time
import multiprocessing as mp
import logging
mpl = mp.log_to_stderr()
mpl.setLevel(logging.INFO)
# local
from nutrient_signaling import cost
from nutrient_signaling import modelreader as md

class Settings:
    def __init__(self):
        # model and data definitions
        self.modelpath = "./data/2018-9-26-12-3-no-sigma"
        self.perturbpath = "./data/yaml/perturbation-data-rebuttal.yaml"
        self.timecoursepath = "./data/yaml/time-course-data-rebuttal.yaml"
        # Output paths
        self.datapath = "./tmp-1/"        
        self.hessianpath =  self.datapath + "Hessians/"
        self.write_psetspath = self.datapath + "Generated-Parameter-Sets/"
        self.lhspath = self.write_psetspath + 'lhs.txt'
        for p in [self.datapath, self.hessianpath, self.write_psetspath]:
            if not os.path.exists(p):
                os.makedirs(p)
        # Run settings
        self.runname = "corrected_data_full_plist"    
        self.simulator = "cpp"
        self.dest = "automatic-iterations/"
        self.debug = True
        self.pathToCombinedPsets = "./Output/2020-11-11/full_guided.txt"
        self.mineig = 0.1
        self.costmultiple = 3.
        self.ReferenceCost = 0.0
        self.lhs_range = 0.025
        self.numPsetsPerIter = 20# 15000
        self.number_of_lhs_sets = 20
        self.startiter = 0
        self.num_iters = 4
        self.exclude_params = exclude_params = [
##################################################    
    # Do NOT Modify this section
           'Carbon','Tps1_T',
           'Glutamine_ext','Sak_T',
           'NH4','Snf1_T',
           'ATP','Sch9_T',
           'Proline','Gln3_T',
           'Cyr1_T','Rtg13_T'
           'PDE_T','Gis1_T',
           'Ras_T','Gcn4_T',
           'PKA_T','Dot6_T',
           'Trehalase_T','Gcn2_T',
           'TORC1_T','Gln1_T',
           'EGO_T','eIF_T',
           'EGOGAP_T','Rib_T',
           'Mig1_T',
           ## The sigma values
           'sigma_gln','sigma_gap',
           'sigma_gln1','sigma_ego',
           'sigma_rtg','sigma_tor',
           'sigma_gis1','sigma_sak',
           'sigma_gcn4','sigma_tps',
           'sigma_dot','sigma_trehalase',
           'sigma_gcn2','sigma_pka',
           'sigma_rib','sigma_ras',
           'sigma_eif','sigma_pde',
           'sigma_cyr',       
           'sigma_mig1',
           'sigma_snf',
           'sigma_sch9',
           ##################################################
           ## The gamma values
           # 'gammacyr', 'gammagap',   
           # 'gammapde', 'gammasynth', 
           # 'gammaras', 'gammaeif',      
           # 'gammapka', 'gamma_rib',  
           # 'gammatre', 'gamma_gcn2', 
           # 'gammatps', 'gamma_mig',  
           # 'gammator', 'gammasnf',   
           # 'gammaego', 'gammasch9',  
           # 'gammagln1',  
           # 'gammagln3',  
       ]
    def setReferenceCost(self, cc):
        c, ct, cp = cc.compute()
        self.ReferenceCost = c
    
def makePars(iternumber, settings):
    runname = settings.runname
    mineig = settings.mineig
    genpestPath = settings.genpestPath
    pathToCombinedPsets = settings.pathToCombinedPsets
    numPsetsPerIter = settings.numPsetsPerIter
    numPsetsPerMachine = settings.numPsetsPerMachine
    costmultiple = settings.costmultiple
    ReferenceCost = settings.ReferenceCost
    
    HessianName = \
        runname + '_' +\
        'hessian' + '_' +\
        str(iternumber) + '.txt'
    #mineig = 0.1
    
    hessian_path = hessianPath + HessianName
    numPSetsForHessian = getNumPsetsLessThanCutoff(pathToCombinedPsets,
                                                   costmultiple,
                                                   ReferenceCost)
    P = np.arange(numPsetsPerMachine,\
                  numPsetsPerIter + numPsetsPerMachine, numPsetsPerMachine)
    for p in P:
        psets_suffix = runname + '_' +\
            str(iternumber) + '_' + str(p)
        guidedPsetGeneration(
                False,             # generate_hessian         
                HessianName,      # hessian_name
                mineig,            # mineig                
                pathToCombinedPsets,       # path_psets_for_hessian
                numPSetsForHessian,# num_psets_for_hessian 
                costmultiple,      # min_cost_multiplier   
                True,              # generate_psets        
                HessianName,      # path_to_hessian       
                numPsetsPerMachine,# number_to_generate    
                psets_suffix,      # psets_suffix          
                False)             # debug          

def get_lowest_cost(psetpath=None):
    df = pd.read_csv(psetpath,index_col=None)
    return(df.cost.min())

def makeHessian(iternumber, settings):
    runname = settings.runname
    hessianpath = settings.hessianpath
    costmultiple = settings.costmultiple
    mineig = settings.mineig
    ReferenceCost = settings.ReferenceCost

    if iternumber == 0:
        pathToCombinedPsets = settings.lhspath
        
    numPSetsForHessian = getNumPsetsLessThanCutoff(pathToCombinedPsets,
                                                   costmultiple,
                                                   ReferenceCost)

    minCostMultiplier = costmultiple
    HessianName = \
        runname + '_' +\
        'hessian' + '_' +\
        str(iternumber) + '.txt'
    hessian_path = hessianPath + HessianName 
    guidedPsetGeneration(
            True,              # generate_hessian         
            HessianName,#hessian_path,      # hessian_path
            mineig,            # mineig                
            pathToCombinedPsets,       # path_psets_for_hessian
            numPSetsForHessian,# num_psets_for_hessian 
            False,             # generate_psets        
            '',                # path_to_hessian       
            0,                 # number_to_generate    
            '',                # psets_suffix          
            True,              # baumann               
            False,             # debug
        settings)
    
def getNumPsetsLessThanCutoff(pathToCombinedPsets,costmultiple, ReferenceCost):
    DF = pd.read_csv(pathToCombinedPsets,sep='\t')
    return DF.shape[0]

def guidedPsetGeneration(generate_hessian, hessian_name,
                         mineig, path_psets_for_hessian,
                         num_psets_for_hessian,
                         generate_psets, path_to_hessian,
                         number_to_generate, psets_suffix, debug,
                         settings):
    """
    Computes Hessian, writes to file, and generates guided parameter sets.
    """
    lowest_cost = get_lowest_cost() #0.0261908040909428

    if not os.path.exists(path_psets_for_hessian):
        print("Parameter Set file to compute Hessian not found")
        sys.exit()
    
    # Read results of cost evaluation carried out on
    # LHS sample
    psetsFile = Path(path_psets_for_hessian).stem
    print(psetsFile)
    
    psetsFile_noext = psetsFile.split('.txt')[0]
    psetsDF = pd.read_csv(path_psets_for_hessian, sep='\t')

    # cleanup file
    if 'index' in list(psetsDF.columns):
        psetsDF = psetsDF.drop(['index'],axis=1)
        
    psetsDF = psetsDF.drop(['c_t','c_p'],axis=1)

    if is_string_dtype(psetsDF['cost']):
        psetsDF = psetsDF.loc[psetsDF['cost'] != 'NAN']
        psetsDF.reset_index(inplace=True,drop=True)
        psetsDF['cost'] = pd.to_numeric(psetsDF['cost'])

    cutoff = settings.costmultiple*lowest_cost
    psetsDF_cutoff = psetsDF.loc[psetsDF['cost'] < cutoff].reset_index(drop=True)
    
    logmessage(path_psets_for_hessian + "\tCutoff="+str(cutoff)+" produces a set containing "\
               +str(psetsDF_cutoff.shape[0]) + ' rows', settings)
    
    min_num_psets_for_hessian = 10
    if generate_hessian:
        if psetsDF_cutoff.shape[0] < min_num_psets_for_hessian:
            m = path_psets_for_hessian + '\t' + hessian_name + '\t ABORTED: NOT ENOUGH PSETS'
            logmessage(m, settings)
            sys.exit()

    # The cutoff only serves to choose which psets build the hessian.
    # This is wrt lowest_cost because lowest_cost establishes a good enough fit
    # The actual min cost in a given pset should contribute to refining
    ## the hessian.
    min_row = psetsDF_cutoff.loc[psetsDF_cutoff['cost'] == min(psetsDF_cutoff['cost'])]
    min_row.reset_index(inplace=True,drop=True)
    min_row = min_row.to_dict(orient='index')[0]
    if generate_hessian:
        if len(hessian_name) == 0:
            print("Please Specify Hessian File Name")
            sys.exit()
        logmessage(path_psets_for_hessian + '\tC_min='+str(min_row['cost']),settings)
        if min_row['cost'] < lowest_cost:
            logmessage(path_psets_for_hessian + "\t Found lower cost than lowest_cost ****", settings)
        Hpath = settings.hessianpath + hessian_name #'H_' + psetsFile_noext + '_'+ hessian_suffix + '.txt'
        ## Compute the Hessian
        ### computeHessian will first check if a hessian file
        ### already exists for a given input file name
        computeHessian(psetsDF_cutoff, psetsFile_noext,
                       Hpath, min_row,  debug)
        print("Wrote Hessian to file")
    if generate_psets:
        if not os.path.exists(settings.hessianpath + path_to_hessian):
            print("Hessian file with that name not found at " + settings.hessianpath)
            sys.exit()
            
        H_mat = pd.read_csv(settings.hessianpath + path_to_hessian, sep='\t',header=None)
        H = H_mat.values
        parnamesDF = pd.read_csv(settings.hessianpath + path_to_hessian.split('.txt')[0] + '_par_names.txt',sep='\t')
        par_names = list(parnamesDF['parnames'])
        DeltaC = 0.5
        # Compute the matrix of eigenvectors and the diagonal matrix containing eigenvalues
        num_pars = H_mat.shape[0]
        W, V = np.linalg.eig(H)

        # This is from Jignesh's script'
        W = abs(W)
        # The numpy way of doing things
        
        W_corrected = W.copy()

        small_eig_count = len(list(np.where(W_corrected < mineig)[0]))
        
        W_corrected[ np.where(W_corrected) < mineig] = mineig
        eig_val_matrix = np.diag(W_corrected)
        eig_vector_matrix = V
        
        logmessage(path_psets_for_hessian + "\tNumber eigenvalues < " + str(mineig) + "= " + str(small_eig_count))

        p_new = {p:[] for p in par_names}
        
        # This is the version used in the paper
        q_vec_skeleton = (DeltaC**(0.5))*np.matmul(eig_vector_matrix,\
                                                   np.linalg.inv(slinalg.sqrtm(eig_val_matrix)))
        
        print(q_vec_skeleton.shape)
        
        for i in tqdm(range(0,number_to_generate)):
            alpha_vector = create_alphavector(num_pars)
            q_vec = list(np.matmul(q_vec_skeleton,alpha_vector)*np.random.normal(0,1))
                  
            for q, p in zip(q_vec,par_names):
                # p*.e^(DeltaC^(0.5)V.Lambda^(-0.5))
                p_new[p].append(float(min_row[p])*np.exp(float(q.real)))
        dirs = path_to_hessian.split('/')
        HessianFile = ''
        for n in dirs:
            if '.txt' in n:
                HessianFile = n
                break
        print(p_new.keys())
        NewPDF = pd.DataFrame(p_new)
        NewPDF.to_csv(settings.write_psetspath + '/guided_'+ psets_suffix +'.txt',index=None,sep='\t')         

def computeHessian(data, fname, Hpath, min_row,  debug):
    print("Computing Hessian using Dr. Baumann's method")
    num_psets = data.shape[0]
    par_names = list(data.columns)
    par_names.remove('cost')

    for excludep in exclude_params:
        if excludep in par_names:
            par_names.remove(excludep)

    par_names_filename = Hpath.split('.txt')[0] + '_par_names.txt'
    with open(par_names_filename,'w') as parfile:
        parfile.write('parnames\n')
        for p in par_names:
            parfile.write(p + '\n')

    num_pars = len(par_names)
    nphv = int(num_pars*(num_pars+1.0)/2.0) # from Baumann's script
    vecV = np.zeros((nphv, 1))
    Q = np.zeros((nphv,nphv))
    L_mat = loadL(num_pars)
    for i, row in tqdm(data.iterrows()):
        q_k = np.zeros((num_pars,1))
        for j in range(0,num_pars):
            denom = min_row[par_names[j]]
            q_k[j] = np.log(row[par_names[j]]) - np.log(denom)
        qqT = np.outer(q_k,q_k.transpose()) # size num
        qqThv = np.matmul(L_mat,qqT.reshape((num_pars*num_pars,1)))
        vecV = vecV + float((row['cost']-min_row['cost']))*qqThv
        Q = Q + np.outer(qqThv,qqThv.transpose())
    print('Q' + str(Q.shape))
    print('vecV' +str(vecV.shape))
    D_mat = loadD(num_pars)
    print('Dmat' + str(D_mat.shape))
    vecH = np.linalg.solve(np.matmul(Q,np.matmul(D_mat.transpose(),D_mat)),vecV)
    print(vecH.shape)
    H = np.reshape(np.matmul(D_mat, vecH), (num_pars,num_pars))
    H_DF = pd.DataFrame(H)
    print('Writing Hessian to file')
    print(Hpath)
    H_DF.to_csv(Hpath,sep='\t',header=False, index=False)
    return(num_pars, par_names)

def loadD(num_pars):
    print('Checking if D_mat exists...')
    if os.path.exists('./Input/duplicator_' + str(num_pars) + '.txt'):
        print('D_mat of dimension ' + str(num_pars) + ' exists')
        print('Reading D_mat')
        D_mat_DF = pd.read_csv('./Input/duplicator_' + str(num_pars) +'.txt',sep='\t',header=None)
        D_mat = D_mat_DF.values
        return D_mat
    else:
        print('Please generate D_mat of size ' + str(num_pars) +\
              ' by calling D(' + str(num_pars) + ') in eliminator_duplicator.py')
        sys.exit()
        
def loadL(num_pars):
    print('Checking if L_mat exists...')
    if os.path.exists('./Input/eliminator_' + str(num_pars) + '.txt'):
        print('L_mat of dimension ' + str(num_pars) + ' exists')
        print('Reading L_mat')
        L_mat_DF = pd.read_csv('./Input/eliminator_' + str(num_pars) +'.txt',sep='\t',header=None)
        L_mat = L_mat_DF.values
        print('Done')
        return L_mat
    else:
        print('Please generate L_mat of size ' + str(num_pars) +\
              ' by calling L(' + str(num_pars) + ') in eliminator_duplicator.py')
        sys.exit()

def logmessage(message, settings):
    print(message)
    with open(settings.datapath + 'LOGS','a') as outfile:
        today = datetime.datetime.today()
        date = str(today.year) + '-' + \
               str(today.month) + '-' + str(today.day) + \
                 '-' + str(today.hour) + '-' + str(today.minute)
        
        outfile.write(date + '\t' + message + "\n")
       

    
def createExpansion(iternumber, runname, costmultiple,
                    datapath, dest, ReferenceCost):

    # Read combine file
    Combine = pd.read_csv(datapath + dest + runname\
                          + "/iter" + str(iternumber) +\
                          "/combine_iter" + str(iternumber) + ".txt",sep='\t')

    if is_string_dtype(Combine['cost']):
        Combine = Combine.loc[Combine['cost'] != 'NAN']
        Combine.reset_index(inplace=True,drop=True)
        Combine['cost'] = pd.to_numeric(Combine['cost'])

    cmin = Combine.cost.min()

    if cmin > ReferenceCost:
        cmin = ReferenceCost

    extract = Combine.loc[Combine.cost <= costmultiple*cmin]

    expansion_path =  datapath + dest + runname +\
        "/iter" + str(iternumber) +\
        "/expansion_iter" + str(iternumber) + ".txt"

    if iternumber == 0:
        # Create expansion_0 as subset of combine_0
        extract.to_csv(expansion_path, sep='\t', index=False)
    else:
        list_ = []
        list_.append(extract)
        prev_expansion = pd.read_csv(datapath + dest + runname\
                                     + "/iter" + str(iternumber - 1) +\
                                     "/expansion_iter" + str(iternumber - 1) +\
                                     ".txt",sep='\t')

        list_.append(prev_expansion)

        D = pd.concat(list_)

        D.to_csv(expansion_path,sep='\t',index=False)
    return expansion_path

def rundone(iternumber, runname, numPsetsPerMachine,datapath,dest):
    doneflag = True
    P = np.arange(numPsetsPerMachine, numPsetsPerIter + numPsetsPerMachine,\
                  numPsetsPerMachine)
    for m,p in zip(runlist, P):
        psets_suffix = runname + '_' +\
            str(iternumber) + '_' + str(p)
        DF = pd.read_csv(datapath + dest + runname +"/iter" + str(iternumber) +\
                         '/eval_' + psets_suffix + '.txt')
        if DF.shape[0] != numPsetsPerMachine:
            doneflag = False
    return doneflag

def combinefiles(iternumber, runname, datapath, dest):
    path = datapath + dest + runname + '/iter' + str(iternumber)
    list_ = []
    P = np.arange(numPsetsPerMachine, numPsetsPerIter + numPsetsPerMachine,\
                  numPsetsPerMachine)

    for m,p in zip(runlist, P):
        psets_suffix = runname + '_' +\
            str(iternumber) + '_' + str(p)

        f = "/eval_" + psets_suffix + ".txt"
        df = pd.read_csv(path + f, sep='\t')
        list_.append(df)
        D = pd.concat(list_)

    D.to_csv(path + "/combine_iter" + str(iternumber) +".txt",sep='\t',index=False)

def LHS(model, settings):
    """
    The LHS function returns a value between 0 and 1
    In order to use this to perturb the parameter value between +-0.05
    we transform the LHS output as follows:
    P_new = P_gset*[1+(0.5-r)*0.01] where r is the value from LHS()
    """
    number_of_sets = settings.number_of_lhs_sets
    lhs_range = settings.lhs_range

    print("You have requested " + str(number_of_sets) + " sets." )
    num_factors = len(model['parameters'].keys())
    
    ## Function Call
    LHSOutput = lhs(num_factors,samples=number_of_sets)
    ###
    PSetDict = {k:[] for k in model['parameters'].keys()}
    # The following will be used from the specified run and not the hand tuned
    for lhsVec in LHSOutput:
        for l, parname, parval in zip(lhsVec,
                                      model['parameters'].keys(),
                                      model['parameters'].values()):
            if parname in settings.exclude_params:
                # Use the default version
                PSetDict[parname].append(parval)
            else:
                ## Use the lhs sample version
                PSetDict[parname].append(parval*\
                                   (1 + ( l*2.0*lhs_range - lhs_range)))
    DF = pd.DataFrame(data=PSetDict)
    print("Writing to file...")
    
    DF.to_csv(settings.lhspath, sep='\t', index=False)
    
def startCostEvaluation(iternumber, cc, settings):
    runname = settings.runname
    numPsetsPerIter = settings.numPsetsPerIter

    # Distribute parameter sets to be evaluated across cores
    if mp.cpu_count() < 5:
        denom = mp.cpu_count()
    else:
        denom = ((mp.cpu_count() // 5) * 5)
    numPsetsPerCore = int(numPsetsPerIter/denom)

    PsetSplits = np.arange(numPsetsPerCore, numPsetsPerIter + numPsetsPerCore,\
                  numPsetsPerCore)
    startp = 0
    if iternumber == 0:                                                        #
        path_to_psets = settings.lhspath                                       #    
    # Debug
    # for pid, endp  in enumerate(PsetSplits):
    #     args = {'cc':copy.copy(cc),
    #             'settings':copy.copy(settings),
    #             'pid':pid,
    #             'startp':int(startp),
    #             'endp':int(endp),
    #             'iternumber':iternumber,
    #     }
    #     startp = int(endp)
    #     compute_cost_wrapper(args)
    ##############################################################################
    # Same problem as async

    # arglist = []                                                               #
    # for pid, endp in enumerate(PsetSplits):                                    #
    #     parsdf = pd.read_csv(path_to_psets,sep='\t')                           #
    #     parsoi = parsdf.loc[startp:endp]                                       #
    #     args = {                                                               #
    #         #'cc':copy.copy(cc),                                               #
    #         'parsoi':parsoi,                                                   #
    #         'outlist':list(),                                                  #
    #         'settings':copy.copy(settings),                                    #
    #         'pid':pid,                                                         #
    #         'startp':int(startp),                                              #
    #         'endp': int(endp),                                                 #
    #         'iternumber':iternumber                                            #
    #     }                                                                      #
    #     startp = int(endp)                                                     #
    #     arglist.append(args)                                                   #
    # processes = []                                                             #
    #     for arg in arglist:                                                    #
    #     processes.append(mp.Process(target=compute_cost_wrapper, args=(arg,))) #
    # for p in processes:                                                        #
    #     p.start()                                                              #
    # for p in processes:                                                        #
    #     p.join()                                                               #
    ##############################################################################

    with mp.Pool(3) as pool:
        jobs = []
        startp = 0
        arglist = []
        for pid, endp in enumerate(PsetSplits):
            parsdf = pd.read_csv(path_to_psets,sep='\t')
            parsoi = parsdf.loc[startp:endp]
            outlist = []            
            args = {
                'parsoi':parsoi,
                'outlist':list(),
                'settings':copy.copy(settings),
                'pid':pid,
                'startp':int(startp),
                'endp': int(endp),
                'iternumber':iternumber
            }
            # Sequential, blocking
            # job = pool.apply(compute_cost_wrapper, args=(args,))
            job = pool.apply_async(compute_cost_wrapper, args=(args,))            
            startp = int(endp)
            jobs.append(job)
        for job in jobs:
            job.wait()            
        
    # This works, at least with 'py'
    # for pid, endp in enumerate(PsetSplits):
    #     parsdf = pd.read_csv(path_to_psets,sep='\t')
    #     parsoi = parsdf.loc[startp:endp]
    #     args = {
    #         #'cc':copy.copy(cc),
    #         'parsoi':parsoi,
    #         'outlist':list(),
    #         'settings':copy.copy(settings),
    #         'pid':pid,
    #         'startp':int(startp),
    #         'endp': int(endp),
    #         'iternumber':iternumber
    #     }
    #     startp = int(endp)
    #     arglist.append(args)
    # with mp.Pool() as pool:
    #     results = pool.map(compute_cost_wrapper, arglist)
    #     pool.close()
    #     pool.join()
        
def compute_cost_wrapper(args):
    settings = args['settings']
    cc = cost.ComputeCost(settings.modelpath, settings.simulator,
                          perturbpath=settings.perturbpath, 
                          timecoursepath=settings.timecoursepath)
    parsoi = args['parsoi']
    pid = args['pid']    

    iternumber = args['iternumber']
    executable = 'main_' + str(pid) + '.o'
    simfilename = 'guided_' + str(pid) + '.dat'
    outlist = args['outlist']
    for i, paramrow in parsoi.iterrows():
        paramset = paramrow.to_dict()        
        args = {'simulator':'cpp',
                'executable':executable,
                'simfilename':simfilename}
        try:
            c, ct, cp = cc.compute(paramset, args)
            outlist.append(dict(paramset,cost=c, c_t=ct, c_p=cp))
            print("success")
        except Exception as e:
            print(pid, "fail", e)
    outdf = pd.DataFrame(outlist)
    outdf.to_csv(f"%s/iter-%d-%d.csv" % (settings.write_psetspath, iternumber, pid))
    print(pid, os.getpid())    
        
def main():
    settings = Settings()
    print("setting up cost obj")
    cc = cost.ComputeCost(settings.modelpath, settings.simulator,
                          perturbpath=settings.perturbpath, 
                          timecoursepath=settings.timecoursepath)
    print("computing ref cost")
    settings.setReferenceCost(cc)
    model = md.readinput(settings.modelpath) 
    parnames = []
    parvalues = []
    for k in model['parameters'].keys():
        if k not in settings.exclude_params:
            parnames.append(k)
            parvalues.append(model['parameters'][k]) 

    for i in range(settings.startiter,settings.num_iters):
        print("In Iter " + str(i))
        if not os.path.exists(settings.datapath +\
                              settings.dest +\
                              settings.runname + "/iter" + str(i)):
            os.makedirs(settings.datapath +\
                        settings.dest +\
                        settings.runname + "/iter" + str(i))
            print("Created folder")

        if i == 0 :
            print("starting lhs gen")
            LHS(model, settings)
            print("starting lhs eval")
            startCostEvaluation(i, cc, settings)
        sys.exit()
        print("\tStarting makeHessian()...")
        makeHessian(i, settings)

        print("\tStarting makePars()")
        makePars(i, settings)
        sys.exit()
        
        print("\tStarting startCostEvaluation()...")
        startCostEvaluation(i, runname, runlist, genpsetsPath,\
                            numPsetsPerIter, numPsetsPerMachine,datapath,dest)
        
        print("\tParameter Sets are being evaluated")
        time.sleep(30)
    
        while not rundone(i, runname, numPsetsPerMachine,datapath,dest):
            time.sleep(wait)
            print("Next check in " + str(wait) + "s")
            
        print("\tEvaluation Complete, combining files...")
        combinefiles(i, runname, datapath, dest)
        
        print("\tExtracting psets to =expansion=")
        pathToCombinedPsets = createExpansion(i, runname, costmultiple,
                                              datapath, dest, ReferenceCost)
        
    if not debug:
        killscreens(runlist)
        
if __name__ == '__main__':
    main()
