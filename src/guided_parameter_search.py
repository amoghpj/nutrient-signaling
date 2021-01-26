# Author: Amogh Jalihal
import os
import sys
import pandas as pd
import scipy.linalg as slinalg
from pandas.api.types import is_string_dtype
import numpy as np
np.seterr(all='raise')
from pyDOE import lhs
import copy
import time
import multiprocessing as mp
import logging
import datetime
from tqdm import tqdm
from pathlib import Path
# mpl = mp.log_to_stderr()
# mpl.setLevel(logging.INFO)
# local
from nutrient_signaling import cost
from nutrient_signaling import utils
from nutrient_signaling import modelreader as md

class Settings:
    def __init__(self):
        # model and data definitions
        self.modelpath = "./data/2018-9-26-12-3-no-sigma"
        self.perturbpath = "./data/yaml/perturbation-data-rebuttal.yaml"
        self.timecoursepath = "./data/yaml/time-course-data-rebuttal.yaml"
        # Output paths
        self.datapath = "./search-redone-2/"        
        self.hessianpath =  self.datapath + "Hessians/"
        self.write_psetspath = self.datapath + "Generated-Parameter-Sets/"
        self.lhspath = self.write_psetspath + 'lhs.csv'
        for p in [self.datapath, self.hessianpath, self.write_psetspath]:
            if not os.path.exists(p):
                os.makedirs(p)
        # Run settings
        self.simulator = "cpp"
        self.dest = "automatic-iterations/"
        self.debug = True
        self.current_pset = ""        
        self.min_num_psets_for_hessian = 100
        self.mineig = 0.001
        self.costmultiple = 3. #10
        self.ReferenceCost = 0.0
        self.lhs_range = 0.025
        ####
        self.numPsetsPerIter = 15000
        self.numproc = 0
        self._get_numproc()
        self.numPsetsPerCore = int(self.numPsetsPerIter/self.numproc)
        self.PsetSplits = np.arange(self.numPsetsPerCore,
                                    self.numPsetsPerIter + self.numPsetsPerCore,\
                                    self.numPsetsPerCore)
        self.number_of_lhs_sets = self.numPsetsPerIter
        self.startiter = 0
        self.num_iters = 6        
        self.expansion_pset =  f"%s/expansion_iter-%d.csv" % (self.write_psetspath, self.startiter)        
        ####
        self.exclude_params =  [
##################################################    
    # Do NOT Modify this section
            'Carbon','Tps1_T','PDE_T','Rtg13_T'
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

    def _get_numproc(self):
        # Distribute parameter sets to be evaluated across cores
        if mp.cpu_count() < 5:
            self.numproc = mp.cpu_count() - 1
        else:
            self.numproc = ((mp.cpu_count() // 5) * 5)

def makePars(iternumber, settings):
    guidedPsetGeneration(
        iternumber, settings,
        generate_hessian=False,             # generate_hessian         
        generate_pars=True,              # generate_pars       
        debug=False)             # debug          

def makeHessian(iternumber, settings):
    guidedPsetGeneration(
        iternumber,settings,
        generate_hessian=True,              # generate_hessian         
        generate_pars=False,             # generate_pars
        debug=False)
    
def get_lowest_cost(iternumber, settings):
    df = pd.read_csv(settings.expansion_pset,index_col=None)
    return(df.cost.min())

def create_alphavector(num_pars):    
    alpha_vector = np.random.normal(0,1, num_pars)
    alpha_vector = alpha_vector/np.linalg.norm(alpha_vector)
    return alpha_vector

def guidedPsetGeneration(iternumber,settings,
                         generate_hessian=False,
                         generate_pars=False,
                         debug=True):
    """
    Computes Hessian, writes to file, and generates guided parameter sets.
    """    

    lowest_cost = get_lowest_cost(iternumber, settings) #0.0261908040909428
    cutoff = settings.costmultiple*lowest_cost
    # Read results of cost evaluation carried out in 
    # previous iteration    
    if not os.path.exists(settings.expansion_pset):
        print("Parameter Set file to compute Hessian not found")
        sys.exit()
    psetsFile = str(Path(settings.expansion_pset).stem)
    print(psetsFile)
    psetsFile_noext = psetsFile.split('.txt')[0]
    psetsDF = pd.read_csv(settings.expansion_pset, index_col=None)
    psetsDF = psetsDF.drop(['c_t','c_p'],axis=1)    
    hessian_name = f"hessian_%d.txt" % (iternumber)        
    if is_string_dtype(psetsDF['cost']):
        psetsDF = psetsDF.loc[psetsDF['cost'] != 'NAN']
        psetsDF.reset_index(inplace=True,drop=True)
        psetsDF['cost'] = pd.to_numeric(psetsDF['cost'])
        
    numPSetsForHessian = psetsDF.shape[0]
    # This is important: Filter out parameter sets in case lowest cost has decreased
    psetsDF_cutoff = psetsDF.loc[psetsDF['cost'] < cutoff].reset_index(drop=True)
    # The cutoff only serves to choose which psets build the hessian.
    # This is wrt lowest_cost because lowest_cost establishes a good enough fit
    # The actual min cost in a given pset should contribute to refining
    ## the hessian.
    min_row = psetsDF_cutoff.loc[psetsDF_cutoff['cost'] == min(psetsDF_cutoff['cost'])]
    min_row.reset_index(inplace=True,drop=True)
    min_row = min_row.to_dict(orient='index')[0]
    Hpath = settings.hessianpath + hessian_name    
    if generate_hessian:
        if psetsDF_cutoff.shape[0] < settings.min_num_psets_for_hessian:
            m = path_psets_for_hessian + '\t' + hessian_name + '\t ABORTED: NOT ENOUGH PSETS'
            logmessage(m, settings)
            sys.exit()
        if min_row['cost'] < lowest_cost:
            logmessage(path_psets_for_hessian + "\t Found lower cost than lowest_cost ****", settings)
        ## Compute the Hessian
        ### computeHessian will first check if a hessian file
        ### already exists for a given input file name
        computeHessian(psetsDF_cutoff,
                       psetsFile_noext,
                       Hpath, min_row,
                       debug, settings)
        print("Wrote Hessian to file")
    if generate_pars:
        if not os.path.exists(Hpath):
            print("Hessian file with that name not found at " + settings.hessianpath)
            sys.exit()
        H_mat = pd.read_csv(settings.hessianpath + hessian_name, header=None, index_col=None)
        H = H_mat.values
        parnamesDF = pd.read_csv(Hpath.split('.txt')[0] + '_par_names.csv',index_col=None)
        par_names = list(parnamesDF['parnames'])
        DeltaC = 0.5
        # Compute the matrix of eigenvectors and the diagonal matrix containing eigenvalues
        num_pars = H_mat.shape[0]
        W, V = np.linalg.eig(H)

        # This is from Jignesh's script'
        W = abs(W)
        # The numpy way of doing things
        W_corrected = W.copy()
        mineig = settings.mineig
        small_eig_count = len(list(np.where(W_corrected < mineig)[0]))
        W_corrected[ np.where(W_corrected <  mineig)] = mineig
        eig_val_matrix = np.diag(W_corrected)
        eig_vector_matrix = V
        # logmessage(path_psets_for_hessian +\
        #            "\tNumber eigenvalues < " + str(mineig) + "= " + str(small_eig_count))
        p_new = {p:[] for p in par_names}
        # This is the version used in the paper
        q_vec_skeleton = (DeltaC**(0.5))*np.matmul(eig_vector_matrix,\
                                                   np.linalg.inv(slinalg.sqrtm(eig_val_matrix)))
        print(q_vec_skeleton.shape)
        for i in tqdm(range(0,settings.numPsetsPerIter)):
            alpha_vector = create_alphavector(num_pars)
            q_vec = list(np.matmul(q_vec_skeleton,alpha_vector)*np.random.normal(0,1))
            for q, p in zip(q_vec,par_names):
                # p*.e^(DeltaC^(0.5)V.Lambda^(-0.5))
                p_new[p].append(float(min_row[p])*np.exp(float(q.real)))
        HessianFile = Path(Hpath).stem
        print(p_new.keys())
        generatedParsDF = pd.DataFrame(p_new)
        generated_pset_path = f"%s/guided_%d.csv" % (settings.write_psetspath, iternumber)
        generatedParsDF.to_csv(generated_pset_path,index=None)
        settings.current_pset = generated_pset_path

def computeHessian(psetsDF_cutoff,
                   psetsFile_noext,
                   Hpath, min_row,
                   debug, settings):
    print("Computing Hessian using Dr. Baumann's method")
    num_psets = psetsDF_cutoff.shape[0]
    par_names = list(psetsDF_cutoff.columns)
    par_names.remove('cost')

    for excludep in settings.exclude_params:
        if excludep in par_names:
            par_names.remove(excludep)
            
    ## Debug
    par_names_filename = Hpath.split('.txt')[0] + '_par_names.csv'
    with open(par_names_filename,'w') as parfile:
        parfile.write('parnames\n')
        for p in par_names:
            parfile.write(p + '\n')
    ## -- 
    num_pars = len(par_names)
    nphv = int(num_pars*(num_pars+1.0)/2.0) # from Baumann's script
    vecV = np.zeros((nphv, 1))
    Q = np.zeros((nphv,nphv))
    L_mat = loadL(num_pars, settings)
    for i, row in tqdm(psetsDF_cutoff.iterrows()):
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
    D_mat = loadD(num_pars, settings)
    print('Dmat' + str(D_mat.shape))
    vecH = np.linalg.solve(np.matmul(Q,np.matmul(D_mat.transpose(),D_mat)),vecV)
    print(vecH.shape)
    H = np.reshape(np.matmul(D_mat, vecH), (num_pars,num_pars))
    H_DF = pd.DataFrame(H)
    print('Writing Hessian to file')
    print(Hpath)
    H_DF.to_csv(Hpath,
                header=False,
                index=False)
    return(num_pars, par_names)

def loadD(num_pars, settings):
    print('Checking if D_mat exists...')
    path_to_D = settings.write_psetspath + '/duplicator_' + str(num_pars) + '.txt'
    if os.path.exists(path_to_D):
        print('D_mat of dimension ' + str(num_pars) + ' exists')
        print('Reading D_mat')
        D_mat_DF = pd.read_csv(path_to_D, header=None, index_col=None)
        D_mat = D_mat_DF.values
        print('Done')
        return D_mat
    else:
        d = utils.D(num_pars)
        pd.DataFrame(d).to_csv(path_to_D, header=False, index=False)
        return d
        
def loadL(num_pars, settings):
    print('Checking if L_mat exists...')
    path_to_L = settings.write_psetspath + '/eliminator_' + str(num_pars) + '.txt'
    if os.path.exists(path_to_L):
        print('L_mat of dimension ' + str(num_pars) + ' exists')
        print('Reading L_mat')
        L_mat_DF = pd.read_csv(path_to_L, header=None, index_col=None)
        L_mat = L_mat_DF.values
        print('Done')
        return L_mat
    else:
        l = utils.L(num_pars)
        pd.DataFrame(l).to_csv(path_to_L, header=False, index=False)
        return l
    
def logmessage(message, settings):
    print(message)
    with open(settings.datapath + 'LOGS','a') as outfile:
        today = datetime.datetime.today()
        date = str(today.year) + '-' + \
               str(today.month) + '-' + str(today.day) + \
                 '-' + str(today.hour) + '-' + str(today.minute)
        
        outfile.write(date + '\t' + message + "\n")
           
def createExpansion(iternumber, settings):
    # Read combine file
    combinedf = pd.read_csv(f"%s/combine_iter-%d.csv" % (settings.write_psetspath, iternumber),
                          index_col=False)

    if is_string_dtype(combinedf['cost']):
        combinedf = combinedf.loc[combinedf['cost'] != 'NAN']
        combinedf.reset_index(inplace=True,drop=True)
        combinedf['cost'] = pd.to_numeric(combinedf['cost'])

    cmin = combinedf.cost.min()
    if cmin > settings.ReferenceCost:
        cmin = settings.ReferenceCost
    cutoffdf = combinedf.loc[combinedf.cost <= settings.costmultiple*cmin]
    expansion_pset = f"%s/expansion_iter-%d.csv" % (settings.write_psetspath, iternumber)
    if iternumber < 2:
        # Create expansion_0 as subset of combine_0
        cutoffdf.to_csv(expansion_pset, index=False)
    else:
        list_ = []
        list_.append(cutoffdf)
        prev_expansion = pd.read_csv(f"%s/expansion_iter-%d.csv" % (settings.write_psetspath,
                                                                    iternumber - 1),
                                     index_col=None)
        list_.append(prev_expansion)
        expansiondf = pd.concat(list_)
        expansiondf.to_csv(expansion_pset,index=False)
    # update path to parameter sets
    settings.expansion_pset = expansion_pset

# def rundone(iternumber, runname, numPsetsPerMachine,datapath,dest):
#     doneflag = True
#     P = np.arange(numPsetsPerMachine, numPsetsPerIter + numPsetsPerMachine,\
#                   numPsetsPerMachine)
#     for m,p in zip(runlist, P):
#         psets_suffix = runname + '_' +\
#             str(iternumber) + '_' + str(p)
#         DF = pd.read_csv(datapath + dest + runname +"/iter" + str(iternumber) +\
#                          '/eval_' + psets_suffix + '.txt')
#         if DF.shape[0] != numPsetsPerMachine:
#             doneflag = False
#     return doneflag

def combinefiles(iternumber, settings):
    list_ = []
    for pid, endp in enumerate(settings.PsetSplits):
        path = f"%s/iter-%d-%d.csv" % (settings.write_psetspath, iternumber, pid)
        print(path)
        df = pd.read_csv(path, index_col=False)
        list_.append(df.reset_index(drop=True))
    D = pd.concat(list_)
    D.to_csv(f"%s/combine_iter-%d.csv" % (settings.write_psetspath, iternumber),
             index=False)
    
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
    
    DF.to_csv(settings.lhspath, index=False)
    
def startCostEvaluation(iternumber, cc, settings):
    numPsetsPerIter = settings.numPsetsPerIter
    numPsetsPerCore = settings.numPsetsPerCore
    parsdf = pd.read_csv(settings.current_pset, index_col=False)            
    # Debug
    # startp = 0    
    # for pid, endp  in enumerate(settings.PsetSplits):
    #     args = {'cc':copy.copy(cc),
    #             'settings':copy.copy(settings),
    #             'pid':pid,
    #             'startp':int(startp),
    #             'endp':int(endp),
    #             'iternumber':iternumber,
    #     }
    #     startp = int(endp)
    #     compute_cost_wrapper(args)
    with mp.Pool(settings.numproc) as pool:
        jobs = []
        startp = 0
        arglist = []
        for pid, endp in enumerate(settings.PsetSplits):
            parsoi = parsdf.loc[startp:endp-1]
            outlist = []            
            args = {
                'parsoi':parsoi,
                'outlist':list(),
                'settings':copy.copy(settings),
                'pid':pid,
                'startp':int(startp),
                'endp': int(endp),
                'iternumber':iternumber}
            job = pool.apply_async(compute_cost_wrapper, args=(args,))            
            startp = int(endp)
            jobs.append(job)
        for job in jobs:
            job.wait()            

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
        except Exception as e:
            print(pid, "fail", e)
    outdf = pd.DataFrame(outlist)
    outdf.to_csv(f"%s/iter-%d-%d.csv" % (settings.write_psetspath, iternumber, pid), index=False)
    # print(pid, os.getpid())    
        
def main():
    start  = time.time()
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
                              "/iter" + str(i)):
            os.makedirs(settings.datapath +\
                        settings.dest +\
                        "/iter" + str(i))
            print("Created folder")

        if i == 0 :
            settings.current_pset = settings.lhspath
            print("starting lhs gen")
            LHS(model, settings)
        else:
            # Get started on making Hessian from parameters generated
            # from previous iteration
            print("\tStarting makeHessian()...")
            makeHessian(i, settings)
            print("\tStarting makePars()")
            makePars(i, settings)
        
        print("\tStarting startCostEvaluation()...")
        startCostEvaluation(i, cc, settings)        
        
        print("\tEvaluation Complete, combining files...")
        combinefiles(i, settings)
        
        print("\tExtracting psets to =expansion=")
        createExpansion(i, settings)

    print(f"This took %0.2fs" % (time.time() - start))
        
if __name__ == '__main__':
    main()
