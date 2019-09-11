"""Collection of utilities for simulating and analyzing model
"""

import PyDSTool as dst
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size':8})
import datetime
import numpy as np
import scipy.linalg as slinalg
import numpy.random as random
#random.seed(0)
#np.random.seed()
import copy
import ast
import pandas as pd
from pandas.api.types import is_string_dtype
import os
import time
import sys
import elimination_duplication as LD
from pyDOE import lhs
from tqdm import tqdm

def readMutantPhenotypeTable():
    """
    Read the file `mutant-phenotype-table.tsv'
    :return MutantDF: Table of mutants
    :type MutantDF: pandas DataFrame        
    """
    mutant_path = "./mutant-phenotype-table.tsv"
    MutantDF = pd.read_csv(mutant_path, sep='\t')
    MutantDF.set_index('sim_id', inplace=True)
    return(MutantDF)


def readTimeSeriesTable():
    """
    :return TimeSeriesDF: Reads `time-series-table.tsv`. Will error out if this file is not  present
    :type TimeSeriesDF: pandas DataFrame
    """
    PATH = "./time-series-table.tsv"
    TimeSeriesDF = pd.read_csv(PATH, sep='\t')
    return(TimeSeriesDF)


def minmaxnorm(input_list):
    """
    Min-Max scale the values in a list. If there is no sufficient 
    difference between min and max in the list, write to log in home directory.

    :param input_list: list of values
    :type input_list: list
    :return output_list: Min-Max scaled values from input_list
    :type output_list: list
    """
    mi = min(input_list)
    mx = max(input_list)
    output_list = []
    for ele in input_list:
        if (mx-mi) < 1e-5:
            home = os.path.expanduser('~')
            with open(home + 'SIMULATION_LOGS','a') as outfile:
                today = datetime.datetime.today()
                outfile.write("\nDIVIDE BY ZERO ENCOUNTERED IN MINMAXNORM on " + str(today.year) + '-' + \
                              str(today.month) + '-' + str(today.day) + \
                              '-' + str(today.hour) + '-' +str(today.minute))
            denom = 1
        else:
            denom = mx - mi
        output_list.append((ele - mi) / (denom))
        
    return(output_list)


def minmaxnorm_special(input_list, mi, mx):
    """
    Like minmaxnorm, but with specified min (mi) and max (mx)

    :param input_list: list of values
    :type input_list: list
    :param mi: user defined min
    :type mi: float
    :param mx: user defined max
    :type mx: float
    :return output_list: List of scaled values
    :type output_list: list
    """
    output_list = []
    for ele in input_list:
        output_list.append((ele - mi) / (mx - mi))
    return(output_list)


def processliteral(string_specification):
    """
    Wraps the `literal_eval` function from python's 
    Abstract Syntax Tree. Used in reading dictionaries 
    storing mutant specifications.

    :param string_specification: String representing a python object
    :type string_specification: str
    :return  S: Object interpreted by the `literal_eval` function.
    """
    S = ast.literal_eval(string_specification)
    return(S)


def createModelObject(model):
    """
    Takes dictionary object as input, and
    returns a PyDSTool Model object

    :param model: Model stored as dictionary. Should contain the keys `variables`, `parameters`, and `initiaconditions`
    :type model: dict
    :return ModelDS: a PyDSTool object that can be simulated to obtain the simulated trajectory
    :type ModelDS: PyDSTool Object
    """
    ModelArgs = dst.args(
        name='test',
        varspecs=model['variables'],
        pars=model['parameters'],
        ics=model['initialconditions'],
        tdata=[0, 60],
        ## Define inline functions that are specific to the NutSig modelx
        ## tRNA() computes min(tRNA_total, Amino acid)
        ## pRib() computes min(Rib, eIF)
        ## shs() is the soft-heaviside function
        ## shsm() is the soft heaviside function with a maximum specified.
        fnspecs={'tRNA': (['tRNA_tot', 'AmAc'], 'min(tRNA_tot, AmAc)'),
                 'pRib': (['rib_comp', 'init_factor'],
                          'min(rib_comp,init_factor)'),
                 'shs': (['sig', 'summation'],
                         '1/(1+e^(-sig*summation))'),
                 'shsm': (['sig', 'm', 'summation'],
                          'm/(1+e^(-sig*summation))')}
    )
    ModelDS = dst.Vode_ODEsystem(ModelArgs)
    return(ModelDS)


def simulateModel(model):
    """
    Takes as input PyDSTool Model object
    and returns PyDSTool Pointset with
    default dt=0.01.

    :param model: PyDSTool Model object
    :return pts: Solution of Model
    :type pts: dict
    """
    pts = model.compute('test').sample()
    return(pts)


def plotTrajectory(points, variable='TORC1'):
    """
    Plots trajectory of single variable in time.
    Default is set to TORC1
    """
    plt.plot(points['t'], points[variable])
    plt.show()


def modelAnalysis(model, perturb, timeseries, robust,
                  num_sets,
                  ts_fit_error_cutoff,
                  allowable_mutant_prediction_error,
                  normaldist,
                  SearchRange,
                  normalsd,
                  vis_only,
                  normalizepred,
                  globalfit,
                  usegoldenset,
                  uselhs,
                  outfilename,
                  debug):
    """
    Calls analysis functions for fitting time series,
    predicting mutants.

    :param model: PyDSTool model object constructed from file
    :type model: PyDSTool object
    :param perturb: Specifies if model comparison to perturbation data should be carried out
    :type perturb: bool
    :param timeseries: Specifies if model comparison to time series data should be carried out
    :type timeseries: bool
    :param robust: Specifies if parameter robustness analysis should be carried out
    :type robust: bool
    :param num_sets: If robustness analysis, how many parameter sets should be simulated?
    :type num_sets: int
    :param ts_fit_error_cutoff: Error cutoff while constructing quadratic cost function
    :type ts_fit_error_cutoff: float
    :param allowable_mutant_prediction_error: Number of mutants that can be wrongly predicted. Obsolete
    :type allowable_mutant_prediction_error: int
    :param normaldist: Sampling distribution from which new parameter is drawn
    :type normaldist: bool
    :param SearchRange: Specifies lower and upper bounds of range to search if using uniform dist.
    :type SearchRange: list
    :param normalsd: Obsolete
    :type normalsd: bool
    :param vis_only: Display plot from perturbation/time series comparisons
    :type vis_only: bool 
    :param normalizepred: scale predicted values between 0 and 1
    :type normalizepred: bool 
    :param globalfit: Obsolete
    :type globalfit: bool 
    :param usegoldenset: Specify whether to use initial parameter set for acceptance criteria. Obsolete
    :type usegoldenset: bool
    :param uselhs: Specify if Latin Hypercube Sampling is to be carried out
    :type uselhs: bool
    :param outfilename: Prefix name for output files
    :type outfilename: str
    :param debug: Write debug logs to file, print verbose messages
    :type debug: bool
    """
    # CHANGELOG:
    # 20180305 Added key 'readout' to experimental_data_dictionary,
    #          this variable is now used for comparison with predictions
    # 20180312 Reading directly from mutant-phenotype-table.tsv instead
    #          of the hardcoded experimental_data_dictionary.
    #          Revamping the structure of the loop to accomodate for the
    #          additional fields
    # 20180319 Now accepts mutlitple types of simulations
    # 20180320 Fixed pandas file reading. Renamed function to modelAnalysis()
    # 20180321 Moved contents to PerturbationAnalysis, added option to
    #          choose between perturbation analysis and time series
    # analysis
    # 20180607 Added option to do robustness analysis

    start=time.clock()
    today = datetime.datetime.today()
    foldername = str(today.year) + '-' + \
                 str(today.month) + '-' + str(today.day) + \
                 '-' + str(today.hour)
    
    if not os.path.exists('./Output/' + foldername):
        os.makedirs('./Output/' + foldername)
    folderpath = './Output/' + foldername

    OutputString=''
    
    if timeseries:
        print("Plotting time series data")
        OutputString += TimeSeriesAnalysis(model, folderpath, normalizepred)
        print("Done!")

    if perturb:
        print("\nStarting perturbation analysis...")
        OutputString += PerturbationAnalysis(model, folderpath,
                                             vis_only, normalizepred)
        print("\nDone!...")

    if robust:
        print("Please use the --cpp option")
        sys.exit()

    writeResultsToFile(OutputString, folderpath)
    print("This analysis took " + str(time.clock() - start) + 's')

"""
Defines functions for processing Mig1 and Rib time courses
"""
customFuncDef = {'mig1': (lambda mig_ts: [np.log10((m)/(1e-3+1-m)) for m in mig_ts]), # ensures denom is not 0
                 'rib': (lambda rib_ts: [rib_ts[i]/(1e-2+rib_ts[0]) for i in range(0,len(rib_ts))])} #ensures denom is not 0

customylabel = {'mig1' : 'log(nMig1/cMig1)',
                    'rib':'Rel Rib'}

def TimeSeriesAnalysis(model, folderpath, normalizepred):
    """
    Read time series data from file, carry out model simulations
    reflecting the wt and mutant strains used in experiments,
    and plot time courses for comparison to data.
    
    :param model: Model defined in file
    :type model: PyDSTool object
    :param folderpath: Path to write output plots
    :type folderpath: str
    :param normalizepred: Specify whether to normalize the predicted time series
    :type normalizepred: bool
    :return OutputString: Write output to org-file
    :type OutputString: str
    """
    TimeSeriesDF = readTimeSeriesTable()

    OutputString = ''
    for i, exp_data in TimeSeriesDF.iterrows():
        print("\t\t Simulating " + str(exp_data.description))
        plt.figure()


        T = processliteral(exp_data.time)
        if exp_data.tunits == 's':
            T = [t / 60.0 for t in T]
            
        if exp_data.ts_type == 'ts':
            if exp_data.lower == 'Nan' and exp_data.upper == 'Nan':
                #print(exp_data['values'])
                V = minmaxnorm(processliteral(exp_data['values']))
            else:
                V = minmaxnorm_special(processliteral(exp_data['values']),
                                       float(exp_data.lower),
                                       float(exp_data.upper))
    
        elif 'func' in exp_data.ts_type:
            V = processliteral(exp_data['values'])
        
        plt.plot(T, V, 'k.')
        plt.xlabel('Time (min)')

    
        model_copy = copy.deepcopy(model)
        
        TimeSeriesDS = createModelObject(model_copy)

        pre_shift_spec = processliteral(exp_data.pre_shift)
        post_shift_spec = processliteral(exp_data.post_shift)
        
        TimeSeriesDS.set(pars=pre_shift_spec['parameters'],
                         ics=pre_shift_spec['inconds'],
                         tdata=[0, 90])
        
        new_ics = get_ss(TimeSeriesDS)

        TimeSeriesDS.set(ics=new_ics)
        TimeSeriesDS.set(pars=post_shift_spec['parameters'],
                         ics=post_shift_spec['inconds'],
                         tdata=[0, np.ceil(T[-1])])

        simRes = simulateModel(TimeSeriesDS)

        if exp_data.ts_type == 'ts':
            if exp_data.readout == 'cAMP' or exp_data.readout == 'Glutamine':
                plt.plot(simRes['t'], simRes[exp_data.readout], 'k--')
            else:
                if not normalizepred:
                    plt.plot(simRes['t'], simRes[exp_data.readout], 'k--')
                else:
                    plt.plot(simRes['t'], minmaxnorm(simRes[exp_data.readout]), 'k--')
            plt.ylabel(exp_data.readout)
            
        elif 'func' in exp_data.ts_type:
            customFuncID = exp_data.ts_type.split(':')[1]
            CustomTimeSeries = customFuncDef[customFuncID](simRes[exp_data.readout])
            plt.plot(simRes['t'], CustomTimeSeries, 'k--')
            plt.ylabel(customylabel[customFuncID])
        #plt.ylim([0,1.0])

        plt.title(exp_data.title)
        plt.savefig(folderpath + "/" +
                    exp_data.description.replace(':', '-') + ".pdf")
        OutputString += '#+ATTR_LATEX: width 0.5\\textwidth\n'
        OutputString += '#+CAPTION: ' + exp_data.description.replace(':', '-')  + '\n'
        OutputString += '[[./' + exp_data.description.replace(':', '-') + '.pdf]]\n\n'

    return(OutputString)


def PerturbationAnalysis(model, folderpath,
                         vis_only, normalizepred):
    """
    Read perturbation data from file, carry out model simulations
    reflecting the wt and mutant strains used in experiments,
    and plot spans for comparison to data.
    
    :param model: Model defined in file
    :type model: PyDSTool object
    :param folderpath: Path to write output plots
    :type folderpath: str
    :param vis_only: Specify whether to display plot instead of writing to file
    :type vis_only: bool
    :param normalizepred: Specify whether to normalize the predicted time series
    :type normalizepred: bool
    :return OutputString: Write output to org-file
    :type OutputString: str
    """
    
    ExperimentalDataDF = readMutantPhenotypeTable()
    Predictions = {}

    Predictions = SimulateRowsInExpDatDF(ExperimentalDataDF, model)
    l, b = ExperimentalDataDF.shape
    f,ax = plt.figure(1,2,figsize=(1., 2.*float(l)*4.))    
    ax_left = ax[0]#f.subplot(121)
    ax_right = ax[1]#plt.subplot(122)
    if vis_only is True:
        ax_list = visualizePerturbAnalysis(ExperimentalDataDF, Predictions,
                                 folderpath, vis_only, normalizepred,[ax_left, ax_right])
        f.subplots_adjust(wspace=0, hspace=0)
        plt.suptitle("High Nutrient -> Low Nutrient")
        plt.tight_layout()
        plt.show()
        return('')
    else:
        imgname = visualizePerturbAnalysis(ExperimentalDataDF,
                                           Predictions, folderpath,
                                           vis_only, normalizepred,[ax_left, ax_right])
        OutputString = ''
        OutputString = writePerturbResults(ExperimentalDataDF,
                                           Predictions,
                                           imgname,
                                           folderpath)
        return(OutputString)

def SimulateRowsInExpDatDF(ExperimentalDataDF, model, parameterSet=None):
    """
    Construct a simulation by reading the configuration specified
    in each line of ExperimentalDataDF

    :param ExperimentalDataDF: Perturbation data read from file
    :type ExperimentalDataDF: pandas DataFrame
    :param model: Model Object read from file
    :type model: PyDSTool Object
    :param parameterSet: Parameter set to use when simulating each row in ExperimentalDataDF
    :type parameterSet: dict
    :return Predictons: Dictionary of predicted steady states
    :type Predictions: dict
    """
    Predictions = {}
    for i, mutant_def in ExperimentalDataDF.iterrows():
        if mutant_def.simulate == 'y':
            description = str(mutant_def['description'])
            sim_id = mutant_def._name
            print("-----------------------------------")
            print("\tSimulation ID: " + str(int(sim_id)))
            print("Description: " + description)

            if mutant_def.sim_type == 'shift':
                Predictions[sim_id] = simulateShifts(mutant_def,
                                                     model, parameterSet)

            elif mutant_def.sim_type == 'treat':
                Predictions[sim_id] = simulateTreatment(mutant_def,
                                                        model, parameterSet)
                
            elif 'func' in mutant_def.sim_type:
                Predictions[sim_id] = simulateFunction(mutant_def,
                                                        model, parameterSet)
    return Predictions
    
def simulateShifts(mutant_def, model, parameterSet):
    """Utility function to simulate perturbations of type `shift`

    :param mutant_def: Dictionary specifying pre- and post-shift states
    :param model: Model read from fil
    :param parameterSet: Parameter Set to be used to simulate shift 
    :return ModelPredictions: Pre and post-shift predictions by model
    :rtype: dict
    """
    #  Make copy of model object
    model_temp = copy.deepcopy(model)

    # Empty dictionary
    ModelPredictions = {}
    # Extract data from DF
    readout = mutant_def.readout
    print("Readout: " + readout)

    wt_pre = mutant_def.wt_pre
    wt_post = mutant_def.wt_post
    mut_pre = mutant_def.mut_pre
    mut_post = mutant_def.mut_post

    normalized_exp_values = minmaxnorm([
        wt_pre,
        wt_post,
        mut_pre,
        mut_post])

    e_wt_pre, e_wt_post, e_mut_pre, e_mut_post = normalized_exp_values

    perturb_spec = processliteral(mutant_def.parameter_change)
    tend_wt = int(mutant_def.tend_wt)
    tend_mut = int(mutant_def.tend_mut)
    print("Simulating Wildtype...")

    WTModelDS = createModelObject(model_temp)

    if parameterSet is not None:
        WTModelDS.set(pars=parameterSet)

    print("Setting sim pre specification")
    WTModelDS.set(ics=get_ss(WTModelDS),
                  tdata=[0, tend_wt], pars=processliteral(mutant_def.pre_spec))
    print("Done!")
    
    ModelPredictions['wt_pre'] = get_ss(WTModelDS)[readout]
    
    whichcondition = mutant_def.whichcondition

    WTModelDS.set(tdata=[0, tend_wt], pars=processliteral(mutant_def.post_spec))

    # Simulate WT with the set initial conditions
    p_wt_post = get_ss(WTModelDS)[readout]

    ModelPredictions['wt_post'] = p_wt_post

    print("Simulating Perturbation...")
    MutantModelDS=createModelObject(model_temp)

    if parameterSet is not None:
        MutantModelDS.set(pars=parameterSet)

    MutantModelDS.set(ics=get_ss(MutantModelDS))
    
    MutantModelDS.set(pars=perturb_spec['parameters'],
                      ics=perturb_spec['inconds'])

    MutantModelDS.set(ics=get_ss(MutantModelDS),
                      tdata=[0, tend_mut],
                      pars=processliteral(mutant_def.pre_spec)
    )

    p_mut_pre = get_ss(MutantModelDS)[readout]
    ModelPredictions['mut_pre'] = p_mut_pre

    MutantModelDS.set(ics=get_ss(MutantModelDS),
                      tdata=[0,tend_mut],
                      pars=processliteral(mutant_def.post_spec))
    
    p_mut_post = get_ss(MutantModelDS)[readout]

    ModelPredictions['mut_post']=p_mut_post
    return(ModelPredictions)


def simulateTreatment(mutant_def, model, parameterSet):
    """Utility function to simulate perturbations of type `treat`

    :param mutant_def: Dictionary specifying pre- and post-shift states
    :param model: Model read from fil
    :param parameterSet: Parameter Set to be used to simulate shift 
    :return ModelPredictions: Pre and post-shift predictions by model
    :rtype: dict
    """    
    #  Make copy of model object
    model_temp = copy.deepcopy(model)

    # Empty dictionary
    ModelPredictions = {}
    # Extract data from DF
    readout = mutant_def.readout
    print("Readout: " + readout)

    wt_pre = mutant_def.wt_pre
    wt_post = mutant_def.wt_post
    mut_pre = mutant_def.mut_pre
    mut_post = mutant_def.mut_post
    
    tend_wt = int(mutant_def.tend_wt)
    tend_mut = int(mutant_def.tend_mut)
    
    normalized_exp_values = minmaxnorm([
        wt_pre,
        wt_post,
        mut_pre,
        mut_post])

    e_wt_pre, e_wt_post, e_mut_pre, e_mut_post = normalized_exp_values

    perturb_spec = processliteral(mutant_def.parameter_change)

    print("Simulating Wildtype...")

    WTModelDS = createModelObject(model_temp)
    if parameterSet is not None:
        WTModelDS.set(pars=parameterSet)

    
    WTModelDS.set(tdata=[0, tend_wt])
    WTModelDS.set(ics=get_ss(WTModelDS))
    WTModelDS.set(pars=perturb_spec[0]['parameters'], # Nutrient Condition
                  ics=perturb_spec[0]['inconds'])
    
    ModelPredictions['wt_pre'] = get_ss(WTModelDS)[readout]

    # Simulate WT with the treatment
    WTModelDS.set(tdata=[0, tend_wt],
                  pars=perturb_spec[0]['parameters'], # Nutrient Condition
                  ics=perturb_spec[0]['inconds'])
    WTModelDS.set(pars=perturb_spec[1]['parameters'], # Treatment
                  ics=perturb_spec[1]['inconds'])

    p_wt_post = get_ss(WTModelDS)[readout]
    
    ModelPredictions['wt_post'] = p_wt_post

    print("Simulating Perturbation+treatment...")

    MutantModelDS=createModelObject(model_temp)

    if parameterSet is not None:
        MutantModelDS.set(pars=parameterSet)


    MutantModelDS.set(pars=perturb_spec[0]['parameters'], # Nutrient Condition
                      ics=perturb_spec[0]['inconds'])
    MutantModelDS.set(pars=perturb_spec[2]['parameters'], # Perturbation
                      ics=perturb_spec[2]['inconds'])
    MutantModelDS.set(ics=get_ss(MutantModelDS),
                      tdata=[0,tend_mut])

    p_mut_pre = get_ss(MutantModelDS)[readout]
    ModelPredictions['mut_pre'] = p_mut_pre

    MutantModelDS.set(pars=perturb_spec[0]['parameters'],
                      ics=perturb_spec[0]['inconds'])
    MutantModelDS.set(pars=perturb_spec[2]['parameters'],
                      ics=perturb_spec[2]['inconds'])
    MutantModelDS.set(pars=perturb_spec[1]['parameters'], # Treatment
                      ics=perturb_spec[1]['inconds'])
    MutantModelDS.set(ics=get_ss(MutantModelDS),
                      tdata=[0,tend_mut])

    p_mut_post = get_ss(MutantModelDS)[readout]

    ModelPredictions['mut_post']=p_mut_post

    return(ModelPredictions)

def simulateFunction(mutant_def, model, parameterSet):
    """Utility function to simulate perturbations of type `func`

    :param mutant_def: Dictionary specifying pre- and post-shift states
    :param model: Model read from fil
    :param parameterSet: Parameter Set to be used to simulate shift 
    :return ModelPredictions: Pre and post-shift predictions by model
    :rtype: dict
    """    
    #  Make copy of model object
    model_temp = copy.deepcopy(model)

    ModelPredictions = {}
    
    readout = mutant_def.readout
    function_definition = mutant_def['sim_type'].split(':')[1]
    print("Special function definition: " + str())
    func_name_dict = {'eifratio': func_eifratio}

    # Function must return a ModelPrediction dictionary
    
    ModelPredictions = func_name_dict[function_definition](mutant_def, model, parameterSet)
                                              
    return(ModelPredictions)
    
def writeResultsToFile(OutputString, folderpath):
    """generate automated string describing summary of results to org file

    :param OutputString: String containing content formatted for org-markdown
    :param folderpath: Path to output org file
    """

    today=datetime.datetime.today()
    outfilename=str(today.year) + '-' + \
                 str(today.month) + '-' + str(today.day) + \
                 '-' + str(today.hour) + ':' + str(today.minute) + '.org'

    with open(folderpath + '/' + outfilename, 'w') as outfile:
        outfile.write("#+TITLE: Model fits to experimental steady states\n")
        outfile.write(r"#+LATEX_HEADER: \usepackage[margin=1in,top=1in,bottom=1in]{geometry}")
        outfile.write("\n")
        outfile.write(r"#+LATEX_HEADER: \usepackage{xcolor}")
        outfile.write("\n")
        outfile.write("#+AUTHOR: Amogh Jalihal\n")
        outfile.write("Notes:\n- values appearing in parenthesis are best fit nutrient input values for the specified condition\n")

        outfile.write("\n\n")

        outfile.write(OutputString)

        outfile.write("\n\nbibliographystyle:unsrt")
        outfile.write("\nbibliography:~/jalihal_projects/Research/references.bib")


def writePerturbResults(exdatDF, pred_dict, imgname, folderpath):
    OutputString=''
    OutputString+=("|*No.*|*Mutant*| *Readout* | *Condition*"
                   " |*WT: Pre to Post*| *Per: Pre to Post*|*Ref*|\n")
    count = 1

    for i, mutant_def in exdatDF.iterrows():
        if mutant_def.simulate=='y':
            
            description = mutant_def['description']
            
            sim_id=mutant_def._name
            
            citation=mutant_def.citation
            
            wt_high = mutant_def.wt_pre
            wt_low = mutant_def.wt_post
            mut_high = mutant_def.mut_pre
            mut_low = mutant_def.mut_post
            
            if mutant_def.sim_type == 'treat':
                condition='treat'
            elif mutant_def.sim_type == 'shift':
                condition = mutant_def.whichcondition
                
            normalized_exp_values = minmaxnorm([
                wt_high,
                wt_low,
                mut_high,
                mut_low])
            
            exp_wt_high, exp_wt_low, \
                exp_mut_high, exp_mut_low = normalized_exp_values
            
            pred_wt_high = pred_dict[sim_id]['wt_pre']
            pred_wt_low = pred_dict[sim_id]['wt_post']
            pred_mut_high = pred_dict[sim_id]['mut_pre']
            pred_mut_low = pred_dict[sim_id]['mut_post']
            
            exp_wt_string = ''
            # mut_wt_string = ''

            if exp_wt_low < exp_wt_high:
                exp_wt_string = 'Decr'
            else:
                exp_wt_string = 'Incr'
                
            if pred_wt_low < pred_wt_high:
                pred_wt_string = 'Decr'
            else:
                pred_wt_string = 'Incr'

            if abs(exp_mut_low - exp_mut_high) < 0.05:
                exp_mut_string = 'Same'
            elif exp_mut_low < exp_mut_high:
                exp_mut_string = 'Decr'
            else:
                exp_mut_string = 'Incr'

            if abs(pred_mut_high - pred_mut_low) < 0.05:
                pred_mut_string='Same'
            elif pred_mut_low < pred_mut_high:
                pred_mut_string='Decr'
            else:
                pred_mut_string='Incr'
                
            if exp_mut_string != pred_mut_string:
                OutputString+= "|" + str(sim_id) + "|" + \
                               description + "|" + \
                               mutant_def.readout + "|" + condition +\
                               "|" + exp_wt_string + "/" + \
                               pred_wt_string + "|\\color{red}" +\
                               exp_mut_string + "/" + pred_mut_string +\
                               "|" + str(citation) + "|\n"
            else:
                OutputString+="|" + str(sim_id) + "|" +\
                               description + "|" +\
                               mutant_def.readout + "|" +\
                               condition + "|" + exp_wt_string +\
                               "/" + pred_wt_string + "|" +\
                               exp_mut_string + "/" + pred_mut_string +\
                               "|" + str(citation) + "|\n"
    
            count+=1

    ## Custom function handling starts
    #count = 1
    writeFunctionDict = {'eifratio':write_eifratio}
    
    for i, mutant_def in exdatDF.iterrows():
        if mutant_def.simulate=='y':
            if 'func' in mutant_def.sim_type:
                customname = mutant_def.sim_type.split(':')[1]
                OutputString += '\n'
                OutputString += writeFunctionDict[customname]\
                                (mutant_def,
                                 pred_dict
                                )
    ## Custom function handling ends
    
    if imgname is not None:
        OutputString+="#+ATTR_LATEX: width: 0.9\\textwidth"
        OutputString+="\n[[./" + imgname + "]]"

    return(OutputString)


def visualizePerturbAnalysis(exdatDF, pred_dict, folderpath, visonly, normalizepred, ax_list):
    """
    Creates a two column span plot with the 'wt' mutant ranges
    on the left and the treatment/mutant ranges on the right

    :param exdatDF: 
    :type exdatDF:  pandas DataFrame
    :param folderpath: output path
    :type folderpath: str
    :param visonly: do not write to file
    :type visonly: boolean
    :param normalizepred: min-max normalize the simulation output
    :type normalizepred: boolean
    :param ax_list: list of axis objects
    :type ax_list: list
    """
    variable_count = 0

    BLACKLIST = []#[10, 12]
    variable_count=0

    conditionSpecificExp = []
    conditionSpecificReadouts = []

    """
    Separate the simulations based on the condition
    """
    classify_plots = {'C': [], 'N': [], 'treat': []}


    for i, mutant_def in exdatDF.iterrows():
        if i not in BLACKLIST and mutant_def.simulate=='y':
            sim_id = mutant_def._name
                
            if mutant_def.sim_type == 'shift':
                if mutant_def.whichcondition == 'C':
                    classify_plots['C'].append(sim_id)
                elif mutant_def.whichcondition == 'N':
                    classify_plots['N'].append(sim_id)
            elif mutant_def.sim_type == 'treat':
                    classify_plots['treat'].append(sim_id)

    # ax_left = plt.subplot(121)
    # ax_right = plt.subplot(122)
    ax_left, ax_right = ax_list
    
    Order_of_plots = ['C', 'N', 'treat']
    line_draw_dec = 0
    
    for oop in Order_of_plots:
        line_draw_dec+=1

        for s in classify_plots[oop]:
            mutant_def = exdatDF.loc[s]
            normalized_exp_values = minmaxnorm([mutant_def.wt_pre,
                                                mutant_def.wt_post,
                                                mutant_def.mut_pre,
                                                mutant_def.mut_post])
            
            e_wt_pre, e_wt_post, \
                e_mut_pre, e_mut_post = normalized_exp_values

            if not normalizepred:
            
                plotWTAutoFits([e_wt_post,
                                e_wt_pre,
                                pred_dict[s]['wt_post'],
                                pred_dict[s]['wt_pre']],
                               variable_count, ax_left)
            
                plotMUTAutoFits([e_mut_post,
                                 e_mut_pre,
                                 pred_dict[s]['mut_post'],
                                 pred_dict[s]['mut_pre']],
                                variable_count, ax_right)

            else:
                normalized_pred_values = minmaxnorm([pred_dict[s]['wt_pre'],
                                                pred_dict[s]['wt_post'],
                                                pred_dict[s]['mut_pre'],
                                                pred_dict[s]['mut_post']])
            
                p_wt_pre, p_wt_post, \
                    p_mut_pre, p_mut_post = normalized_pred_values
            
                plotWTAutoFits([e_wt_post,
                                e_wt_pre,
                                p_wt_post,
                                p_wt_pre],
                               variable_count, ax_left)
            
                plotMUTAutoFits([e_mut_post,
                                 e_mut_pre,
                                 p_mut_post,
                                 p_mut_pre],
                                variable_count, ax_right)
            if mutant_def.annotation == '0':
                annotation = ''
            else:
                annotation = ' ' + mutant_def.annotation
            
            conditionSpecificReadouts.append(str(mutant_def.readout) +
                                                 ' [' + oop + ']' + annotation)

            if mutant_def.tend_mut != mutant_def.tend_wt:
                conditionSpecificExp.append(str(mutant_def['description'])\
                                            + ' ' + str(mutant_def.tend_mut)\
                                            + ' min')
            else:
                conditionSpecificExp.append(str(mutant_def['description']))
            variable_count+=1
            
        ax_left.axhline(variable_count, color='black')
        ax_right.axhline(variable_count, color='black')
        variable_count+=1
        conditionSpecificReadouts.append('')
        conditionSpecificExp.append('')

    variable_count-=1
    ax_left.set_ylim([-1, variable_count])
    ax_left.set_yticks(np.arange(0, variable_count, 1.0))
    ax_left.set_yticklabels(conditionSpecificReadouts)
    ax_left.set_title("Wildtype")

    ax_right.set_ylim([-1, variable_count])
    ax_right.set_yticks(np.arange(0, variable_count, 1.0))
    ax_right.set_yticklabels(['' for v in conditionSpecificExp])
    ax_right_twin=ax_right.twinx()
    ax_right_twin.set_ylim(ax_right.get_ylim())
    ax_right_twin.set_yticks(ax_right.get_yticks())
    ax_right_twin.set_yticklabels(conditionSpecificExp)
    ax_right.set_title("Perturbation")
    
    # plt.subplots_adjust(wspace=0, hspace=0)
    # plt.suptitle("High Nutrient -> Low Nutrient")
    # plt.tight_layout()

    if visonly is True:
        return([ax_left,ax_right])
        
    else:
        today=datetime.datetime.today()
        outfilename = str(today.year) + '-' + \
                      str(today.month) + '-' + str(today.day) + \
                      '-' + str(today.hour) + '_' + str(today.minute) + \
                      '-Perturbation-Summary.pdf'
        imgname = outfilename
        plt.tight_layout()        
        plt.savefig(folderpath + '/' + imgname)
        return(imgname)


def plotWTAutoFits(simvals, variable_count, ax):
    """
    DATE: 2018-01-11
    WT: black
    High: circle
    Low: square
    Experimental: filled
    Predicted: empty
    """
    e_low, e_high, \
        p_low, p_high = simvals

    ax.plot([e_low, e_high],
            [variable_count + 0.15, variable_count + 0.15],
            'k')
    ax.plot([p_low, p_high],
            [variable_count - 0.15, variable_count - 0.15],
            'k--')

    if e_low < e_high:
        ax.plot(e_low, variable_count + 0.15, 'k', marker='<', markersize=10)
        ax.plot(e_high, variable_count + 0.15, 'k', marker='.', markersize=10)
    else:
        ax.plot(e_low, variable_count + 0.15, 'k', marker='>', markersize=10)
        ax.plot(e_high, variable_count + 0.15, 'k', marker='.', markersize=10)

    if p_low < p_high:
        ax.plot(p_low, variable_count - 0.15, 'k', marker='<',
                markersize=10, markerfacecolor='white',
                markeredgecolor='black')
        ax.plot(p_high, variable_count - 0.15, 'k', marker='.',
                markersize=10, markerfacecolor='white',
                markeredgecolor='black')

    else:
        ax.plot(p_low, variable_count - 0.15, 'k', marker='>',
                markersize=10, markerfacecolor='white',
                markeredgecolor='black')
        ax.plot(p_high, variable_count - 0.15, 'k', marker='.',
                markersize=10, markerfacecolor='white',
                markeredgecolor='black')


def plotMUTAutoFits(simvals, variable_count, ax):
    """
    DATE: 2018-01-11
    MUT: red
    High: circle
    Low: square
    Experimental: filled
    Predicted: empty
    """
    e_low, e_high, \
        p_low, p_high = simvals
    
    ax.plot([e_low, e_high],
            [variable_count + 0.15, variable_count + 0.15],
            'r')
    ax.plot([p_low, p_high],
            [variable_count - 0.15, variable_count - 0.15],
            'r--')

    if e_low<e_high:
        ax.plot(e_low, variable_count + 0.15,
                'r', marker='<', markersize=10)
        ax.plot(e_high, variable_count + 0.15,
                'r', marker='.', markersize=10)
    else:
        ax.plot(e_low, variable_count + 0.15,
                'r', marker='>', markersize=10)
        ax.plot(e_high, variable_count + 0.15,
                'r', marker='.', markersize=10)

    if p_low < p_high:
        ax.plot(p_low, variable_count - 0.15, 'r',
                marker='<', markersize=10,
                markerfacecolor='white', markeredgecolor='red')
        ax.plot(p_high, variable_count - 0.15, 'r',
                marker='.', markersize=10,
                markerfacecolor='white', markeredgecolor='red')

    else:
        ax.plot(p_low, variable_count - 0.15, 'r',
                marker='>', markersize=10,
                markerfacecolor='white', markeredgecolor='red')
        ax.plot(p_high, variable_count - 0.15, 'r',
                marker='.', markersize=10,
                markerfacecolor='white', markeredgecolor='red')


def searchPspace(model, variable, value_to_fit, whichcondition, tol=0.1):
    """Search Parameter Space

    :param model: 
    :param variable: 
    :param value_to_fit: 
    :param whichcondition: 
    :param tol: 
    :returns: 
    :rtype: 

    """


    """
    CHANGELOG
    20180105 Added third condition 'Bb' for both nitrogen and carbon
    20180307 Now returns the best fit input values
    """
    value=float(value_to_fit)
    found_best_fit=False
    best_fit_so_far=2.0  # Some large value that has to be minimized
    
    if whichcondition=='c' or whichcondition=='C':

        CarbonValues = np.linspace(0.00, 1.0, 15)
        for c in CarbonValues:
            M = dst.deepcopy(model)
            M.set(pars={"Carbon": c,
                        "ATP": c,
                        "Nitrogen": 1.0})
            SSPoints=get_ss(M)
            if abs(SSPoints[variable] - value) <= tol:

                print("Best fit parameter C=" + str(c))
                found_best_fit=True
                return(SSPoints[variable], c)
                break
            if abs(SSPoints[variable] - value) <= abs(best_fit_so_far - value):
                best_fit_so_far, best_fit_input=SSPoints[variable], c

    if whichcondition=='n' or whichcondition=='N':
        LOWER_LIM=1e-4
        UPPER_LIM=1.0

        NitrogenValues = np.logspace(np.log10(LOWER_LIM),
                                     np.log10(UPPER_LIM),
                                     15)
        for n in NitrogenValues:
            M=dst.deepcopy(model)
            M.set(pars={"Carbon": 1.0,
                        "ATP": 1.0,
                        "Nitrogen": n})
            SSPoints=get_ss(M)
            if abs(SSPoints[variable] - value) <= tol:

                print("Best fit parameter N=" + str(n))
                found_best_fit = True
                return(SSPoints[variable], n)
                break
            if abs(SSPoints[variable] - value) <= abs(best_fit_so_far - value):
                best_fit_so_far, best_fit_input = SSPoints[variable], n

    if whichcondition=='b' or whichcondition=='B':
        LOWER_LIM=1e-4
        UPPER_LIM=1.0
        print("Range of N values to check = [" +
              str(LOWER_LIM) + ", " + str(UPPER_LIM) + "]")
        NitrogenValues = np.logspace(np.log10(LOWER_LIM),
                                     np.log10(UPPER_LIM), 15)

        print("Range of C values to check = [0, 1]")
        CarbonValues=np.linspace(0.00, 1.0, 15)

        for n in NitrogenValues:
            for c in CarbonValues:
                M=dst.deepcopy(model)
                M.set(pars={"Carbon": c,
                            "ATP": c,
                            "Nitrogen": n})
                SSPoints=get_ss(M)
                if abs(SSPoints[variable] - value) <= tol:
    
                    print("Best fit parameter N, C=" + str(n) + ", " + str(c))
                    found_best_fit=True
                    return(SSPoints[variable], (c, n))
                    break
                if abs(SSPoints[variable] - value) <= abs(best_fit_so_far - value):
                    best_fit_so_far, best_fit_input=SSPoints[variable], (c, n)

    if found_best_fit is False:
        print("No best fit found within tolerance limit,  returning best fit ")
        return(best_fit_so_far, best_fit_input)


def get_ss(model):
    """
    Returns steady states of all variables in the model

    :return SSPoints: Dictionary containing steady state values
    """
    Points=simulateModel(model)
    SSPoints={}
    for k in Points.keys():
        SSPoints[k]=Points[k][-1]
    return(SSPoints)

def growth_rate(points):
    """
    Compute predicted growth rate Lambda according to
    Scott et al 2010, calibrated to data from Metzl-Raz et al, 2017
    """
    k_NegativeFeedBack=20.0
    PhiRmax=0.325
    PhiR0=0.08
    rho=0.76

    if isinstance(points['AApool'], float):
        kn = 1.0 / (1.0 + (points['AApool'] / k_NegativeFeedBack) ** 2)
        kt = points['TORC1'] * (1 - points['Snf1']) * (points['PKA'])
    else:
        kn = 1.0 / (1.0 + (points['AApool'][-1] / k_NegativeFeedBack) ** 2)
        kt = points['TORC1'][-1] * (1 - points['Snf1'][-1]) * (points['PKA'][-1])
        
    Lambda = 6.5 * ((PhiRmax - PhiR0) / rho) * kt * kn / (kt + kn)
    # Returns in minutes
    if Lambda != 0.0:
        Lambda_minutes=float(int(6000.0 / Lambda)) / 100.0
        return(Lambda_minutes)
    else:
        return(Lambda)

def sse(parvalues, expValues, model, parnames, varname):
    """FIXME! briefly describe function

    :param parvalues: 
    :param expValues: 
    :param model: 
    :param parnames: 
    :param varname: 
    :returns: 
    :rtype: 

    """
    print(parvalues)
    expY, expt = expValues
    model = model
    parameters = {}
    
    for i in range(0,len(parvalues)):
        parameters[parnames[i]] = parvalues[i]

    modelDS = copy.deepcopy(model)
    modelDS.set(pars=parameters)
    
    P = modelDS.compute('test').sample(dt=0.1)
    X = [x for x in P[varname]]
    T = P['t']
    errors = []
    for i in range(0, len(expt)):
        e = float( X[ findindex(T, expt[i]) ] - expY[i])
        errors.append(e)
    del modelDS
    return(errors)


def findindex(t, expt,tol=1e-1):
    """FIXME! briefly describe function

    :param t: 
    :param expt: 
    :param tol: 
    :returns: 
    :rtype: 

    """
    my_index = 0
    for T in t:
        if abs(T - expt) < tol:
            return my_index
        else:
            my_index += 1


def ParameterEstimation(model, simspec, InitialGuess):
    """FIXME! briefly describe function

    :param model: 
    :param simspec: 
    :param InitialGuess: 
    :returns: 
    :rtype: 

    """
    if InitialGuess is None:
        print("\033[1;31;40m !!Please specify initial parameter guess!!")
        sys.exit()

    if type(InitialGuess) == str:
        print(InitialGuess)
        InitialGuess = processliteral(InitialGuess)
        print(type(InitialGuess))


    print(InitialGuess)
    parnames = InitialGuess.keys()
    parvalues = InitialGuess.values()

    TimeSeriesDF = readTimeSeriesTable()

    if TimeSeriesDF.loc[TimeSeriesDF.description == simspec].empty:
        print("\033[1;31;40m !!Incorrect specification!!")
        sys.exit()
        
    data = TimeSeriesDF.loc[TimeSeriesDF.description == simspec]
    for i, d in data.iterrows():
        exp_dat = d
    varname = exp_dat['readout']    
    """
    Extract experimental values
    """
    expT = processliteral(exp_dat['time'])
    if exp_dat.tunits =='s':
        expT = [t/60.0 for t in expT]

    if exp_data.ts_type == 'ts':
        if exp_data.lower == 'Nan' and exp_data.upper == 'Nan':
            #print(exp_data['values'])
            expV = minmaxnorm(processliteral(exp_data['values']))
        else:
            expV = minmaxnorm_special(processliteral(exp_data['values']),
                                   float(exp_data.lower),
                                   float(exp_data.upper))
            
    elif 'func' in exp_data.ts_type:
        expV = processliteral(exp_data['values'])

    #expV = minmaxnorm(processliteral(exp_dat['values']))
    #expV = minmaxnorm_special(processliteral(exp_dat['values']), float(exp_dat.lower), float(exp_dat.upper))

    expValues = [expV, expT]
    
    plt.plot(expT, expV, 'r.')

    """
    Get simulation conditions from table
    """
    preSpec = processliteral(exp_dat['pre_shift'])
    postSpec = processliteral(exp_dat['post_shift'])
    
    modelDS = createModelObject(model)
    
    modelDS.set(pars=preSpec['parameters'],
                ics=preSpec['inconds'])

    NewICs = get_ss(modelDS)

    modelDS.set(ics=NewICs)
    modelDS.set(pars=postSpec['parameters'],
                ics=postSpec['inconds'],
                tdata=[0.0, expT[-1]])
    P = simulateModel(modelDS)
    if exp_data.ts_type == 'ts':
        if exp_data.readout == 'cAMP' or exp_data.readout == 'Glutamine':
            plt.plot(simRes['t'], simRes[exp_data.readout], 'k--')
        else:
            if not normalizepred:
                plt.plot(simRes['t'], simRes[exp_data.readout], 'k--')
            else:
                plt.plot(simRes['t'], minmaxnorm(simRes[exp_data.readout]), 'k--')
        plt.ylabel(exp_data.readout)
            
    elif 'func' in exp_data.ts_type:
        customFuncID = exp_data.ts_type.split(':')[1]
        CustomTimeSeries = customFuncDef[customFuncID](simRes[exp_data.readout])
        plt.plot(simRes['t'], CustomTimeSeries, 'k--')

    """
    Simulate the default parameters
    """
    #plt.plot(P['t'], P[varname])

    # Kinitial = {'w_pka_camp': 2.0, 'w_pka_sch9':3.0}


    BestFit = leastsq(sse,
                      parvalues,
                      args=(expValues,
                            copy.deepcopy(modelDS),
                            parnames,
                            varname),
                      full_output=True)


    KFit = BestFit[0]
    parameters = {}

    for i in range(0,len(KFit)):
        parameters[parnames[i]] = KFit[i]

    print(parameters)
        
    modelDS.set(pars=parameters, tdata=[0, max(expT)])

    P = simulateModel(modelDS)
    plt.plot(P['t'], P[varname], 'g')
    plt.show()



def mysse(parvalues, parnames, sim_spec, datadict, globalfit):
    """FIXME! briefly describe function

    :param parvalues: 
    :param parnames: 
    :param sim_spec: 
    :param datadict: 
    :param globalfit: 
    :returns: 
    :rtype: 

    """
    #print(parvalues)
    parameters = {}
    errors = []

    parameters = {pname:pvalue for pname,pvalue in zip(parnames,parvalues)}

    for condition in sim_spec:
        if not globalfit:
            errors.append([])
        """
        Get simulation conditions from table
        """
        expV = datadict[condition]['v']
        expT = datadict[condition]['time']
        varname = datadict[condition]['readout']
        # modelDS = copy.deepcopy(datadict[condition]['model'])
        modelDS = datadict[condition]['model'] ##PROBLEMATIC?
        
        modelDS.set(pars=parameters)

        P = modelDS.compute('test').sample(dt=0.1)
        if 'func' in datadict[condition]['ts_type']:
            X = customFuncDef[datadict[condition]['ts_type'].split(':')[1]](P[varname])
        else:
            X = [x for x in P[varname]]
        T = P['t']

        for i in range(0, len(expT)):
            e = float( X[ findindex(T, expT[i]) ] - expV[i] )
            if not globalfit:
                errors[-1].append(e)
            else:
                errors.append(e)
        del modelDS
    return(errors)


def MyParameterEstimation(model, specList, InitialGuess):
    if InitialGuess is None:
        print("\033[1;31;40m !!Please specify initial parameter guess!!")
        sys.exit()

    if type(InitialGuess) == str:
        print(InitialGuess)
        InitialGuess = processliteral(InitialGuess)
        print(type(InitialGuess))

    if type(specList) == str:
        specList = processliteral(specList)

    print(InitialGuess)
    parnames = InitialGuess.keys()
    parvalues = InitialGuess.values()
    parameters = {}

    # for i in range(0,len(parvalues)):
    #     parameters[parnames[i]] = parvalues[i]

    TimeSeries = readTimeSeriesTable()

    DataOfInterest = TimeSeries.loc[TimeSeries.description.isin(specList)]
    datadict = dict()

    for i, row in DataOfInterest.iterrows():
        print(row)
        datadict[row.description] = {'time': processliteral(row.time),
                                     'v': processliteral(row['values']),
                                     'ts_type':str(row['ts_type']),
                                     'pre_shift': processliteral(row.pre_shift),
                                     'post_shift': processliteral(row.post_shift),
                                     'readout': row.readout}

    Number_of_Conditions = len(specList)

    f, ax = plt.subplots(Number_of_Conditions, 1)
    #modelDS = createModelObject(model)

    i = 0
    for condition in specList:
        """
        Get simulation conditions from table
        """
        
        preSpec = datadict[condition]['pre_shift']
        postSpec = datadict[condition]['post_shift']

        expV = minmaxnorm(datadict[condition]['v'])
        expT = datadict[condition]['time']
        if Number_of_Conditions == 1:
            ax.plot(expT, expV, 'r.')
        else:
            ax[i].plot(expT, expV, 'r.')
        # Reinitialize
        modelDS = createModelObject(model)

        modelDS.set(pars=preSpec['parameters'],
                    ics=preSpec['inconds'])

        NewICs = get_ss(modelDS)

        modelDS.set(ics=NewICs)
        modelDS.set(pars=postSpec['parameters'],
                    ics=postSpec['inconds'],
                    tdata=[0.0, expT[-1]])
        P = simulateModel(modelDS)

        if len(datadict.keys())==1:
            ax.plot(P['t'], P[datadict[condition]['readout']], 'b')
        else:
            ax[i].plot(P['t'], P[datadict[condition]['readout']], 'b')
        i += 1         

        datadict[condition]['model'] = copy.deepcopy(modelDS)
        del modelDS
        
    BestFit = leastsq(mysse,
                      parvalues,
                      args=(parnames,
                            specList,
                            datadict,
                            False),
                      full_output=True)

    # BestFit = leastsq(sse,
    #                   parvalues,
    #                   args=([datadict[specList[0]]['v'],datadict[specList[0]]['time']],
    #                         modelDS,
    #                         parnames,
    #                         datadict[specList[0]]['readout'])
    # )
    KFit = BestFit[0]
    parameters = {}

    for i in range(0,len(KFit)):
        parameters[parnames[i]] = KFit[i]

    #print(parameters)

    i = 0
    for condition in specList:
        # Reinitialize
        modelDS = copy.deepcopy(datadict[condition]['model']) #createModelObject(model)
        # # modelDS.set(pars=model['parameters'],
        # #             ics=model['initialconditions'])

        preSpec = datadict[condition]['pre_shift']
        postSpec = datadict[condition]['post_shift']
        modelDS.set(pars=preSpec['parameters'],
                    ics=preSpec['inconds'])
        modelDS.set(ics=get_ss(modelDS))
        modelDS.set(pars=postSpec['parameters'],
                    ics=postSpec['inconds'])

        modelDS.set(pars=parameters, tdata=[0, max(expT)])

        P = simulateModel(modelDS)
        if Number_of_Conditions == 1:
            ax.plot(P['t'], P[datadict[condition]['readout']], 'g')
        else:
            ax[i].plot(P['t'], P[datadict[condition]['readout']], 'g')
        del modelDS
        i += 1
    plt.show()


def writeRobustnessParameters(num_sets,
                              ts_fit_error_cutoff,
                              SearchRange,
                              allowable_mutant_prediction_error,
                              normaldist,
                              normalsd):
    OutputString = ''
    OutputString += '\n* Robustness Analysis'
    OutputString += '\nParameters:\n'
    
    OutputString += '|Number of sets  | ' + str(num_sets) + '|\n'
    OutputString += '|Time Series Error cutoff  | ' + str(ts_fit_error_cutoff) + '|\n'
    OutputString += '|Search range | ' + str(SearchRange) + '|\n'
    OutputString += '|Allowable mutant prediction error | ' + str(allowable_mutant_prediction_error) + '|\n'
    #OutputString += '|Global Fit|' + str(globalfit) + '|\n'
    if normaldist:
        OutputString += '|Distribution|Normal, std='+str(normalsd)+'|\n'
    else:
        OutputString += '|Distribution| Uniform|\n'

    return OutputString
    
## Currently hard coded the list of time series data to compare against
specList = [
    'Snf1:glucose_repression',
    'wt:cAMP',
    'sch9Delta:cAMP',
    'Sch9P:100_glutamine',
    'Sch9P:1_glutamine',
    'Sch9:gtr1Delta',
    'Sch9:glucose_starve',
    'Sch9:glucose_relief',
    'Mig1:glucose_relief',
    'Rib:rap',
    'Gln3:rap'
]

exclude_params = [
##################################################    
    # Do NOT Modify this section
    'Carbon',
    'Glutamine_ext',
    'NH4',
    'ATP',
    'Proline',
    'Cyr1_T',
    'PDE_T',
    'Ras_T',
    'PKA_T',
    'Trehalase_T',
    'Tps1_T',
    'Sak_T',
    'Snf1_T',
    'Sch9_T',
    'Gln3_T',
    'Rtg13_T',
    'Gis1_T',
    'Gcn4_T',
    'Dot6_T',
    'Gcn2_T',
    'TORC1_T',
    'EGO_T',
    'EGOGAP_T',
    'Gln1_T',
    'eIF_T',
    'Rib_T',
    'Mig1_T',
    ##################################################
    ## The sigma values
    'sigma_gln',
    'sigma_gln1',
    'sigma_rtg',
    'sigma_gis1',
    'sigma_gcn4',
    'sigma_dot', 
    'sigma_gcn2',
    'sigma_rib', 
    'sigma_eif',
    'sigma_gap',
    'sigma_ego',  
    'sigma_tor',
    'sigma_sak', 
    'sigma_tps', 
    'sigma_trehalase', 
    'sigma_pka',
    'sigma_ras', 
    'sigma_pde',
    'sigma_cyr',
    'sigma_mig1',
    'sigma_snf',
    'sigma_sch9',
    ##################################################
    ## The gamma values
    # 'gammacyr',
    # 'gammapde',
    # 'gammaras',    
    # 'gammapka',
    # 'gammatre',
    # 'gammatps',
    # 'gammator',
    # 'gammaego',
    # 'gammagap',
    # 'gammasynth',
    # 'gammaeif',
    # 'gamma_rib',
    # 'gamma_gcn2',
    # 'gamma_mig',
    # 'gammasnf',
    # 'gammasch9',
    # 'gammagln1',
    # 'gammagln3',
    ##################################################
    # 'w_sng_glc',
    # 'w_snf_sak',
    # 'w_snf',
    
    # 'k_degr',
    # 'k_acc_glu',
    # 'k_acc_pro',
    # 'k_acc_nh4',
    
    # 'w_cyr',
    # 'w_cyr_glu',
    # 'w_cyr_snf',
    
    
    # 'w_pde',
    # 'w_pde_pka',
    
    # 'w_ras_pka',
    # 'w_ras_glu',

    # 'w_ras',
    
    
    # 'w_pka',
    # 'w_pka_camp',
    # 'w_pka_sch9',
    
    
    # 'w_tre_pka',
    # 'w_tre',
    
    # 'w_tps_pka',
    # 'w_tps',
    
    # 'w_sak',
    # 'w_sak_pka',
     
    # 'w_torc_ego',
    # 'w_torc_egoin',
    # 'w_torc_glut',
    # 'w_torc',
    # 'w_torc_snf',
    
    
    # 'w_ego',
    # 'w_ego_gap',
    # 'w_ego_basal',
    
    # 'w_gap_torc',
    # 'w_gap_N',
    
    
    # 'w_synthase_tor',
    # 'w_synthase',
    # 'k_pr',
    
    # 'w_eif',
    # 'w_eif_gcn2',
    
    # 'w_rib_dot',
    # 'w_rib',
    
    # 'w_gcn',
    # 'w_gcn_torc',
    # 'tRNA_total',
    
    # 'w_dot_sch_pka',
    # 'w_dot',
    
    # 'w_gcn4_gcn2_trna',
    # 'w_gcn4',
    # 'tRNA_sensitivity',
    
    # 'w_gis_pka',
    # 'w_gis_sch',
    # 'w_gis',
    
    # 'w_rtg_torc',
    # 'w_rtg',
    
    
    # 'w_gln_torc',
    # 'w_gln_snf',
    # 'w_gln_sit',
    # 'k_camp_cyr',
    # 'k_camp_pde',
    # 'k_camp_deg',
    
    # 'w_mig_pka',
    # 'w_mig_snf',
    # 'w_mig',
    
    # 'w_snf_glc',
    # 'w_snf_sak',
    # 'w_snf',
    
    # 'w_sch9_torc',
    # 'w_sch9',
    # 'k_transcription',
    # 'k_mRNA_degr',
    
    # ####################################################################################################
    # # # Experiment: exclude the negative (reverse reaction) omegas
    # ##################################################

    # 'k_degr',
    # # 'k_acc_glu',
    # # 'k_acc_pro',
    # # 'k_acc_nh4',

    # 'w_cyr',
    # #                   'w_cyr_glu',
    # #                   'w_cyr_snf',
    
    
    # 'w_pde',
    # #                   'w_pde_pka',
    
    # #                   'w_ras_pka',
    # #                   'w_ras_glu',
    
    # 'w_ras',
    
    
    # 'w_pka',
    # #                   'w_pka_camp',
    # #                   'w_pka_sch9',
    
    
    # #                   'w_tre_pka',
    # 'w_tre',
    
    # #                   'w_tps_pka',
    # 'w_tps',
    
    # 'w_sak',
    # #                   'w_sak_pka',
    
    # #                   'w_torc_ego',
    # #                   'w_torc_egoin',
    # #                   'w_torc_glut',
    # 'w_torc',
    # #                   'w_torc_snf',
    
    
    # #                   'w_ego',
    # #                   'w_ego_gap',
    # 'w_ego_basal',
    
    # 'w_gap_torc',
    # #                   'w_gap_N',
    
    
    # #                   'w_synthase_tor',
    # #                   'w_synthase',
    # 'k_pr',
    
    # 'w_eif',
    # #                   'w_eif_gcn2',
    
    # #                   'w_rib_dot',
    # #                   'w_rib',
    
    # 'w_gcn',
    # #                   'w_gcn_torc',
    # #                   'tRNA_total',
    
    # #                   'w_dot_sch_pka',
    # 'w_dot',
    
    # #                   'w_gcn4_gcn2_trna',
    # 'w_gcn4',
    # #                   'tRNA_sensitivity',
    
    # #                   'w_gis_pka',
    # #                   'w_gis_sch',
    # 'w_gis',
    
    # #                   'w_rtg_torc',
    # 'w_rtg',
    
    # 'w_gln3',
    # #                   'w_gln_torc',
    # #                   'w_gln_snf',
    # #                   'w_gln_sit',
    
    # 'k_camp_cyr',
    # #                   'k_camp_pde',
    # #                   'k_camp_deg',
    
    # # 'w_gln1_gln3',
    # 'w_gln1',
    # #                   'w_mig_pka',
    # #'w_mig_snf',
    # 'w_mig',
    # 'w_snf',    
    # #                   'w_snf_glc',
    # #                   'w_snf_sak',
    # #                   'w_sch9_torc',
    # 'w_sch9',
    # 'k_mRNA_degr'
    
]
# Include


def extractParsAndWriteHeaderText(foldername, outfilename_prefix, model,
                                  exclude_params,
                                  ExperimentalDataDF, uselhs):
    parnames = []
    parvalues = []
    mut_ids = []
    #HeaderText = '#Index\t'
    HeaderText = ''

    summary_file_path = foldername + '/' + outfilename_prefix + '.txt'

    if os.path.exists(summary_file_path):
        filemode = 'a'
    else:
        filemode = 'w'

    print(filemode)
    with open(summary_file_path,filemode) as outfile:
        for k in model['parameters'].keys():
            if k not in exclude_params:
                if filemode == 'w':
                    HeaderText += k + '\t'
                parnames.append(k)
                parvalues.append(model['parameters'][k])
                
        if filemode == 'w':
            if uselhs:
                HeaderText += 'c_t\tc_p\tcost\n'
            else:
                HeaderText += 'c_ts\tc_p\tcost\tdecision\n'
        outfile.write(HeaderText)
        
    return (parnames, parvalues, mut_ids)


def checkMutantPredictions(PSetPredictions, GoldenSetPredictions, mut_ids):
    # Return a decision based on the GoldenPrediction and the
    # Prediction made by the parameterSet
    mutantDecision = {}
    for mutantID in mut_ids:
        wpre_pset = PSetPredictions[mutantID]['wt_pre']
        wpost_pset = PSetPredictions[mutantID]['wt_post']
        mpre_pset = PSetPredictions[mutantID]['mut_pre']
        mpost_pset = PSetPredictions[mutantID]['mut_post']

        wpre_gold = GoldenSetPredictions[mutantID]['wt_pre']
        wpost_gold = GoldenSetPredictions[mutantID]['wt_post']
        mpre_gold = GoldenSetPredictions[mutantID]['mut_pre']
        mpost_gold = GoldenSetPredictions[mutantID]['mut_post']


        wt_decision = ((wpre_gold-wpost_gold) -\
                       (wpre_pset - wpost_pset))\
                       /((wpre_gold-wpost_gold))

        if (mpre_gold-mpost_gold) > 1e-2:
            m_decision = ((mpre_gold-mpost_gold) -\
                          (mpre_pset - mpost_pset))\
                          /((mpre_gold-mpost_gold))
        else:
            # Heuristic to take care of small numbers
            if (mpre_pset - mpost_pset) > 1e-1:
                m_decision = 0.5 # REJECT
            else:
                m_decision = 0.01 # ACCEPT
        
        
        if abs(wt_decision) < 0.1 and abs(m_decision) < 0.1:
            mutantDecision[mutantID]='PASS'
        else:
            mutantDecision[mutantID]='FAIL'

    return mutantDecision
        
def writePSetToFile(line, folderpath, outfilename_prefix):
    with open(folderpath + '/' + outfilename_prefix + '.txt','a') as outfile:
        outfile.write(line)
        
def makeDecision(decisionFuncArgs,
                 parnames, parvalues,
                 mut_ids,
                 allowable_mutant_prediction_error,
                 ts_fit_error_cutoff,
                 accepted_parameter_set_number,
                 folderpath,
                 outfilename_prefix):
    """FIXME! briefly describe function

    :param decisionFuncArgs: 
    :param parnames: 
    :param parvalues: 
    :param mut_ids: 
    :param allowable_mutant_prediction_error: 
    :param ts_fit_error_cutoff: 
    :param accepted_parameter_set_number: 
    :param folderpath: 
    :param outfilename_prefix: 
    :returns: 
    :rtype: 

    """
    
    line = '' #str(accepted_parameter_set_number) + '\t'
    
    # Write parameter values
    for i in range(0,len(parnames)):
        line += str(parvalues[i]) + '\t'

    printl = ''
    decisionFlag = False
    conditionalAccept = False
    gsetCost, psetCosts, LastAcceptedCost, uselhs = decisionFuncArgs
    
    # psetCosts in the plural
    psetCost, psetCost_splitup = psetCosts
    
    if uselhs:
        line += str(psetCost_splitup[0]) + '\t'+ str(psetCost_splitup[1]) + '\t' +str(psetCost) + '\n'
        printl = ''#str(round(LastAcceptedCost[0],4)) + "\t" + str(round(LastAcceptedCost[1],4)) + "\t"
        writePSetToFile(line, folderpath, outfilename_prefix)
        accepted_parameter_set_number += 1
        return (decisionFlag, accepted_parameter_set_number, LastAcceptedCost, printl)

    else:
        # Sample from exponential distribution

        DiffCost = psetCost - LastAcceptedCost
        #beta = 0.1#3.6
        #beta = 50.0
        beta = 1000.0        
        prob = min(1,np.exp(-beta*DiffCost))
        randS = random.rand()
        # DiffCost will always satisfy if it is negative.        
        if randS < prob:
            line += str(psetCost_splitup[0]) + '\t'+ str(psetCost_splitup[1]) + '\t' +str(psetCost) + "\t"
            if prob == 1:
                decisionFlag = True
                line += "ACCEPT\n"
                LastAcceptedCost = psetCost                                
                printl += str(round(LastAcceptedCost,4))
            else:
                conditionalAccept = True
                line += "WEAK ACCEPT\n"
                LastAcceptedCost = psetCost                                                

        if decisionFlag:
            writePSetToFile(line, folderpath, outfilename_prefix)
            printl += "\tACCEPT\t"
            accepted_parameter_set_number += 1
        else:
            if conditionalAccept:
                writePSetToFile(line, folderpath, outfilename_prefix)
                printl += "\tWEAK ACCEPT"
                accepted_parameter_set_number += 1

            else:
                printl += "\tREJECT"
        return (decisionFlag, accepted_parameter_set_number, LastAcceptedCost, printl)
            
def ALessThanEqualB(A, B):
    """Utility function to test magnitude of two quantities.

    :param A: float or list
    :param B: float or list
    :returns flag: decision after comparison of A and B
    :rtype: bool

    """
    flag = True
    if isinstance(A, float):
        if A > B:
            return False

    else:
        for a,b in zip(A, B):
            if a > b:
                return False
    
    return flag
            
def parameterSetFinder(LastAcceptedSet,
                       normaldist=False,
                       SearchRange=10,
                       normalsd=0.06):
    """TODO briefly describe function

    :param LastAcceptedSet: 
    :param normaldist: 
    :param SearchRange: 
    :param normalsd: 
    :returns: 
    :rtype: 

    """

    prevParamSet = LastAcceptedSet
    newset = dict()
    for k in prevParamSet.keys():
        if not normaldist:
            newset[k] = random.uniform((1.0 - SearchRange/200.0)*prevParamSet[k],\
                                       (1.0 + SearchRange/200.0)*prevParamSet[k])
        else:
            newset[k] = random.normal(prevParamSet[k],scale=normalsd*prevParamSet[k])
        
    return newset
    
        
def PCA(paramset_path,annotate):
    from sklearn.decomposition import PCA
    print("Starting Principal component analysis...")
    fileextension = paramset_path.split('.')[-1]
    print(fileextension)
    if fileextension != 'txt' and fileextension != 'tsv':
        print("Invalid input. Please provide path to a table in a txt or tsv.")
        sys.exit()
    print("Reading table")
    DF = pd.read_csv(paramset_path,
                     sep='\t',
                     index_col=None)

    # DF = DF.loc[DF['decision']=='ACCEPT']
    # DF = DF.reset_index(drop=True)
    #print(DF)
    
    today = datetime.datetime.today()
    foldername = str(today.year) + '-' + \
                 str(today.month) + '-' + str(today.day) + \
                 '-' + str(today.hour)
    
    if not os.path.exists('./Output/' + foldername):
        os.makedirs('./Output/' + foldername)
    folderpath = './Output/' + foldername

    colnames = list(DF.columns)
    non_parameters = ['cost','c_ts', 'c_p', 'decision']
    for np in non_parameters:
        colnames.remove(np)

    features = []
    for f in colnames:
        if f not in exclude_params:
            features.append(f)
    # features =  [
    #     # 'k_degr',
    #     # 'k_acc_glu',
    #     # 'k_acc_pro',
    #     # 'k_acc_nh4',
    #     # 'gammacyr',
    #     # 'sigma_cyr',
    #     # 'w_cyr',
    #     # 'w_cyr_glu',
    #     # 'w_cyr_snf',
    #     'k_camp_cyr',
    #     'k_camp_pde',
    #     'k_camp_deg',
    #     # 'sigma_pde',
    #     # 'gammapde',
    #     # 'w_pde',
    #     # 'w_pde_pka',
    #     # 'w_ras_pka',
    #     # 'w_ras_glu',
    #     # 'gammaras',
    #     # 'w_ras',
    #     # 'sigma_ras',
    #     # 'w_pka',
    #     # 'w_pka_camp',
    #     # 'w_pka_sch9',
    #     # 'gammapka',
    #     # 'sigma_pka',
    #     # 'gammatre',
    #     # 'sigma_trehalase',
    #     # 'w_tre_pka',
    #     # 'w_tre',
    #     # 'gammatps',
    #     # 'sigma_tps',
    #     # 'w_tps_pka',
    #     # 'w_tps',

    #     # 'w_sak',
    #     # 'w_sak_pka',
    #     # 'sigma_sak',

    #     'gammasnf',
    #     'sigma_snf',
    #     'w_snf_glc',
    #     'w_snf_sak',

    #     # 'gammator',
    #     # 'sigma_tor',
    #     # 'w_torc_ego',
    #     # 'w_torc_egoin',
    #     # 'w_torc_glut',
    #     # 'w_torc',
    #     # 'w_torc_snf',

    #     # 'gammaego',
    #     # 'sigma_ego',
    #     # 'w_ego',
    #     # 'w_ego_gap',
    #     # 'w_ego_basal',

    #     'gammasch9',
    #     'sigma_sch9',
    #     'w_sch9_torc',
    #     'w_sch9',

    #     # 'gammagap',
    #     # 'sigma_gap',
    #     # 'w_gap_torc',
    #     # 'w_gap_N',

    #     # 'gammasynth',
    #     # 'sigma_synthase',
    #     # 'w_synthase_tor',
    #     # 'w_synthase',
    #     # 'k_pr',

    #     # 'sigma_eif',
    #     # 'w_eif',
    #     # 'w_eif_gcn2',
    #     # 'gammaeif',

    #     # 'gamma_rib',
    #     # 'sigma_rib',
    #     # 'w_rib_dot',
    #     # 'w_rib',

    #     # 'gamma_gcn2',
    #     # 'sigma_gcn2',
    #     # 'w_gcn',
    #     # 'w_gcn_torc',
    #     # 'tRNA_total',

    #     # 'sigma_dot',
    #     # 'w_dot_sch_pka',
    #     # 'w_dot',

    #     # 'sigma_gcn4',
    #     # 'w_gcn4_gcn2_trna',
    #     # 'w_gcn4',
    #     # 'tRNA_sensitivity',

    #     # 'sigma_gis1',
    #     # 'w_gis_pka',
    #     # 'w_gis_sch',
    #     # 'w_gis',
    #     'sigma_mig1',
    #     'w_mig_pka',
    #     'w_mig_snf',
    #     'w_mig',
    #     'gamma_mig',

    #     # 'sigma_rtg',
    #     # 'w_rtg_torc',
    #     # 'w_rtg',

    #     # 'sigma_gln',
    #     # 'w_gln_torc',
    #     # 'w_gln_snf',
    #     # 'w_gln_sit'
    # ]
    
    # # 'Gln3_T',    
    # # 'Rtg13_T',        
    # # 'Cyr1_T',
    # # 'PDE_T',
    # # 'Ras_T',
    # # 'PKA_T',
    # # 'Trehalase_T',
    # # 'Tps1_T',
    # # 'Sak_T',
    # #'Snf1_T',
    # # 'TORC1_T',        
    # # 'EGO_T',        
    # # 'Sch9_T',        
    # # 'EGOGAP_T',        
    # # 'Synthase_T',        
    # # 'eIF_T',        
    # # 'Rib_T',        
    # # 'Gcn2_T',
    # # 'Dot6_T',
    # # 'Gis1_T',
    # # 'Gcn4_T',
        
    # #print(DF.columns)
    
    X = DF.loc[:, features].values
    
    #Y = DF.loc[:, ['decision']].values
    
    from sklearn.preprocessing import StandardScaler
    
    X_std = StandardScaler().fit_transform(X)
    
    pca = PCA(n_components=2)
    
    principalComponents = pca.fit_transform(X_std)
    #print(pca.components_)
    # print("PC1 largest weight")
    # print(features[list(pca.components_[0]).index(max(pca.components_[0]))])
    # print("PC1 smallest weight")
    # print(features[list(pca.components_[0]).index(min(pca.components_[0]))])

    # print("PC2 largest weight")
    # print(features[list(pca.components_[1]).index(max(pca.components_[1]))])
    # print("PC2 smallest weight")
    # print(features[list(pca.components_[1]).index(min(pca.components_[1]))])

    
    ## This gives us the PC1 and PC2 values for each row
    principalDF = pd.DataFrame(data=principalComponents, columns=\
                               ['PC1',
                                'PC2'])
                                # 'PC3',
                                # 'PC4',
                                # 'PC5'])
    ## To help with the plotting, we add the 'decision' col to
    ## the principalDF dataframe
    print(DF['decision'])
    finalDF = pd.concat([principalDF, DF[['decision']]], axis=1)
    #print(finalDF)
    targets = ['ACCEPT', 'WEAK ACCEPT']
    colors = ['#264f12', '#b2451a']
    alpha = [1, 0.3]
    plt.close()

    if not annotate:
    ## For each row, color the point using the decision
        for target, color, a in zip(targets, colors, alpha):
            indicestoKeep = finalDF['decision'] == target
            plt.scatter(finalDF.loc[indicestoKeep, 'PC1'],
                        finalDF.loc[indicestoKeep, 'PC2'],
                        c=color, alpha=a)
            # plt.annotate(indicestoKeep, (finalDF.loc[indicestoKeep, 'PC1'],
            #                              finalDF.loc[indicestoKeep, 'PC2'],))
            plt.xlabel('PC1')
            plt.ylabel('PC2')
    else:
        f = plt.figure()
        ax = f.add_subplot(111)

        print("plotting")
        for target, color in zip(targets, colors):
            indicestoKeep = pd.DataFrame(finalDF['decision'] == target)
            indicestoKeep.columns = ['Bool']
            indiceslist = indicestoKeep.index[indicestoKeep['Bool']].tolist()
            #print(indiceslist)
            for i in indiceslist:
                xcoord = float(finalDF['PC1'].loc[i])
                ycoord = float(finalDF['PC2'].loc[i])
                #print((xcoord,ycoord))
                ax.annotate(s=str(i),
                             xy=(xcoord,ycoord),
                             xycoords='data',color=color
                )
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_xlim([-10,10])
        ax.set_ylim([-10,10])


    # #ax.annotate('Hello', xy=(0,0),xycoords='data')
    # ax.annotate('There', xy=(0.5,0.5),xycoords='data')
    # ax.annotate('!', xy=(0,0.5),xycoords='data')
    plt.show()
    #plt.savefig(folderpath + '/' + foldername + '-PCA.png')


def func_eifratio(mutant_def, model, parameterSet):
    ## TODOS in this function
    ## eIF measurement from Urban, 2007 Fig7B
    ## - Calculate wt ss ratio of (eIF)/(1-eIF)
    ## - Calculate wt ss ratio after rap treatment
    ##   - this can be encoded in the table
    ## - Calculate rap ratio to wt ratio. approx 4
    ## - Repeat for sch9 delta
    ## Instead of repeating low level definitions, I have decided to
    ## simply use the output of simulateTreatment in this case.
    
    ## Notes:
    ## In Fig 4B of Urban 2007, effect of rapamycin treatment is
    ## reported as the ratio of eIF-P/eIF level in the treatment to
    ## that in the control
    ## wt_pre is thus 1
    ## wt_post is 4
    ## mut_pre is 5.6
    ## mut_post is 5.8

    MPFromSimulateTreatment = simulateTreatment(mutant_def,
                                                model,
                                                parameterSet)

    CorrectedModelPredictions = {} # This will do the ratio calculations
    eIF_wt = MPFromSimulateTreatment['wt_pre']

    eIFP_eIF_wt = (1 - eIF_wt)/eIF_wt
    
    CorrectedModelPredictions['wt_pre'] = round(eIFP_eIF_wt/\
                                                eIFP_eIF_wt,2)

    eIF = MPFromSimulateTreatment['wt_post']
    CorrectedModelPredictions['wt_post'] = round(((1-eIF)/eIF)/\
                                           eIFP_eIF_wt,2)

    eIF = MPFromSimulateTreatment['mut_pre']
    CorrectedModelPredictions['mut_pre'] = round(((1-eIF)/eIF)/\
                                           eIFP_eIF_wt,2)

    eIF = MPFromSimulateTreatment['mut_post']
    CorrectedModelPredictions['mut_post'] = round(((1-eIF)/eIF)/\
                                             eIFP_eIF_wt,2)

    # This should ideally work
    return CorrectedModelPredictions
    
def write_eifratio(mutant_def, pred_dict):
    ## Some form of the following:
    description = mutant_def['description']
    
    sim_id=mutant_def._name
    
    citation=mutant_def.citation

    OutputString = ''

    OutputString += description + '\n'
    OutputString += 'Ratio of eIF-P/eIF treatment to control cf. ' + mutant_def.source +\
                    ' Urban, 2007 ' +\
                    '[['  + citation + ']]\n'

    OutputString += '|Condition|Reported|Simulated|\n'
    OutputString += '|wt|' + str(mutant_def.wt_pre) +\
                    '|' + str(pred_dict[sim_id]['wt_pre']) + '|\n'
    
    OutputString += '|wt + rap|' + str(mutant_def.wt_post) +\
                    '|' + str(pred_dict[sim_id]['wt_post']) + '|\n' 

    OutputString += '|sch9 $\Delta$|' + str(mutant_def.mut_pre) +\
                    '|' + str(pred_dict[sim_id]['mut_pre']) + '|\n'
    
    OutputString += '|sch9 $\Delta$ + rap|' + str(mutant_def.mut_post) +\
                    '|' + str(pred_dict[sim_id]['mut_post']) + '|\n' 

    return OutputString
    
def LHS(number_of_sets, lhs_range, filename, ModelDefinition):
    """
    The LHS function returns a value between 0 and 1
    In order to use this to perturb the parameter value between +-0.05
    we transform the LHS output as follows:
    P_new = P_gset*[1+(0.5-r)*0.01] where r is the value from LHS()
    """
    print("You have requested " + str(number_of_sets) + " sets." )
    num_factors = len(ModelDefinition['parameters'].keys())
    
    ## Function Call
    LHSOutput = lhs(num_factors,samples=number_of_sets)
    ###
    PSetDict = {k:[] for k in ModelDefinition['parameters'].keys()}
    AllParNames = [k for k in ModelDefinition['parameters'].keys()]
    # The following will be used from the specified run and not the hand tuned
    for lhsVec in LHSOutput:
        for p,k in zip(lhsVec, AllParNames):
            if k in exclude_params:
                # Use the default version
                PSetDict[k].append(ModelDefinition['parameters'][k])
            else:
                ## Use the lhs sample version
                PSetDict[k].append(ModelDefinition['parameters'][k]*\
                                   (1 + ( p*2.0*lhs_range - lhs_range)))
    DF = pd.DataFrame(data=PSetDict)
    print("Writing to file...")
    
    DF.to_csv("./" + filename, sep='\t', index=False)



def writePsetToFile(ModelDefinition,Path_to_PSet):
    data = pd.read_csv(Path_to_PSet, sep='\t')
    foldername = raw_input('Please specify foldername for output: ')
    if not os.path.exists('./ParameterSets/' + foldername):
        os.makedirs('./ParameterSets/' + foldername)

    new_Pars_dict = ModelDefinition['parameters']
    pset_dict = data.loc[data['cost'] == min(list(data['cost']))].reset_index(drop=True).to_dict(orient='index')[0]
    print(pset_dict['cost'])
    #sys.exit()
    for c in ['cost','c_ts', 'c_p', 'decision']:
        del pset_dict[c]
    new_Pars_dict.update(pset_dict)
    print('Writing variables to file')
    with open('./ParameterSets/' + foldername + '/variables.txt','w') as outfile:
        for k,v in ModelDefinition['variables'].iteritems():
            outfile.write(k + '\t' + v + '\n')
    print('Done')
    print('Writing ICs to file')
    with open('./ParameterSets/' + foldername + '/initialconditions.txt','w') as outfile:
        for k,v in ModelDefinition['initialconditions'].iteritems():
            outfile.write(k + '\t' + str(float(v)) + '\n')
    print('Done')
    print('Writing parameters to file')
    with open('./ParameterSets/' + foldername + '/parameters.txt','w') as outfile:
        for k,v in new_Pars_dict.iteritems():
            outfile.write(k + '\t' + str(float(v)) + '\n')
    print('All Done!')
    
    
def guidedPsetGeneration(generate_hessian, hessian_name,
                         mineig, path_psets_for_hessian,
                         num_psets_for_hessian,
                         min_cost_multiplier,                         
                         generate_psets, path_to_hessian,
                         number_to_generate, psets_suffix,baumann, debug):
    """
    Computes Hessian, writes to file, and generates guided parameter sets.
    """
    # .. py:function:: guidedPsetGeneration(generate_hessian, hessian_suffix,
    #                      mineig, path_psets_for_hessian,
    #                      num_psets_for_hessian,
    #                      generate_psets, path_to_hessian,
    #                      number_to_generate, psets_suffix, debug)
    
    # Hard coded paths! Modify these by hand-------------------->
    HESSIANPATH = "/data/amogh-jalihal/Nutrient-Signaling/Hessians/"
    WRITE_PSETSPATH = "/data/amogh-jalihal/Nutrient-Signaling/Generated-Parameter-Sets/"
    # Hard coded min MCMC cost
    MCMCmin = 0.0261908040909428
    # <--------------------------------------------------
    if not os.path.exists(path_psets_for_hessian):
        print("Parameter Set file to compute Hessian not found")
        sys.exit()
    
    # Read results of cost evaluation carried out on
    # LHS sample
    dirs = path_psets_for_hessian.split('/')
    psetsFile = ''
    for n in dirs:
        if '.txt' in n:
            psetsFile = n
            break
    #print(psetsFile)
    
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


    # Filter to acceptable cost cutoff
    costs = sorted(list(psetsDF['cost']))
    ##################################################
    ##################################################
    ##################################################
    # Default 6
    cutoff = min_cost_multiplier*MCMCmin
    ##################################################
    ##################################################
    ##################################################    
    psetsDF_cutoff = psetsDF.loc[psetsDF['cost'] < cutoff].reset_index(drop=True)
    
    logmessage(path_psets_for_hessian + "\tCutoff="+str(cutoff)+" produces a set containing "\
          +str(psetsDF_cutoff.shape[0]) + ' rows')

    # Now, how do I handle cases where I don't have enough psets?
    # A: expand cutoff

    if generate_hessian:
        if psetsDF_cutoff.shape[0] < num_psets_for_hessian:
            print(path_psets_for_hessian + '\t' + hessian_name + '\t ABORTED: NOT ENOUGH PSETS')
            logmessage(path_psets_for_hessian + '\t' + hessian_name + '\t ABORTED: NOT ENOUGH PSETS')
            sys.exit()
    
        else:
            # Prune        
            psetsDF_cutoff = psetsDF_cutoff.iloc[0:num_psets_for_hessian]

    # else:
    #     psetsDF_cutoff is unchanged 

    
    # The cutoff only serves to choose which psets build the hessian.
    # This is wrt MCMCmin because MCMCmin establishes a good enough fit
    # The actual min cost in a given pset should contribute to refining
    ## the hessian.
        
    min_row = psetsDF_cutoff.loc[psetsDF_cutoff['cost'] == min(psetsDF_cutoff['cost'])]

    min_row.reset_index(inplace=True,drop=True)

    min_row = min_row.to_dict(orient='index')[0]
        
    if generate_hessian:
        
        if len(hessian_name) == 0:
            print("Please Specify Hessian File Name")
            sys.exit()
        
        logmessage(path_psets_for_hessian + '\tC_min='+str(min_row['cost']))
        if min_row['cost'] < MCMCmin:
            logmessage(path_psets_for_hessian + "\t Found lower cost than MCMCmin ****")
        
        Hpath = HESSIANPATH + hessian_name #'H_' + psetsFile_noext + '_'+ hessian_suffix + '.txt'
    
        ## Compute the Hessian
        # computeHessian will first check if a hessian file already exists for a given input file name
        

        if baumann:
            computeHessian_baumann(psetsDF_cutoff, psetsFile_noext,
                                   Hpath, min_row, num_psets_for_hessian, debug)
        else:
            computeHessian(psetsDF_cutoff, psetsFile_noext,
                           Hpath, min_row, num_psets_for_hessian, debug)            

        print("Wrote Hessian to file")
    
            
    if generate_psets:
        if not os.path.exists(HESSIANPATH + path_to_hessian):
            print("Hessian file with that name not found at " + HESSIANPATH)
            sys.exit()
            
        H_mat = pd.read_csv(HESSIANPATH + path_to_hessian, sep='\t',header=None)
        H = H_mat.values
        parnamesDF = pd.read_csv(HESSIANPATH + path_to_hessian.split('.txt')[0] + '_par_names.txt',sep='\t')
        par_names = list(parnamesDF['parnames'])
        #DeltaC = 2.0 # Used by Jignesh#1.0#3.0#4e-4
        DeltaC = 0.5 # Approx 20 times min.
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
        # else:
        #     # Dr. Baumann's implementation
        #     print("Coming soon")
        #     sys.exit()
        
        NewPDF = pd.DataFrame(p_new)
        NewPDF.to_csv(WRITE_PSETSPATH + '/guided_'+ psets_suffix +'.txt',index=None,sep='\t')

    
def create_alphavector(num_pars):    
    alpha_vector = [random.normal(0,1) for i in range(0,num_pars)]
    alpha_vector = np.array(alpha_vector)
    alpha_vector = alpha_vector/np.linalg.norm(alpha_vector)
    return alpha_vector

def computeHessian(data, fname, Hname, min_row, psetsForHessian, debug):
    # if data.shape[0] > psetsForHessian:
    #     data = data.iloc[0:psetsForHessian]
    #     #print(data.head())
    #     logmessage('Top ' + str(psetsForHessian)+ ' psets used to generate ' + Hname)

    print("Computing Hessian")
    num_psets = data.shape[0]
    par_names = list(data.columns)
            
    if 'index' in data.columns:
        print('shit')
        sys.exit()
    par_names.remove('cost')

    for excludep in exclude_params:
        if excludep in par_names:
            par_names.remove(excludep)
    print(len(par_names))
    par_names_filename = Hname.split('.txt')[0] + '_par_names.txt'
    with open(par_names_filename,'w') as parfile:
        parfile.write('parnames\n')
        for p in par_names:
            parfile.write(p + '\n')
    
    num_pars = len(par_names)
    if debug:
        print('num_pars='+str(num_pars))
        print('num_psets='+str(num_psets))

    # First check if D, L exist for the given dimension
    D_mat = loadD(num_pars)
    L_mat = loadL(num_pars)

    # Initialize the V matrix, nxn
    V = np.zeros((num_pars,num_pars))
    # Initialize the Q matrix, n^2xn^2
    Q = np.zeros((num_pars**2, num_pars**2))
    # Initialize matrix of q vectors, Kxn for K psets
    q_mat = np.zeros((num_psets,num_pars))
    
    COSTS = np.array((num_psets,1))
    
    # Identify the p_opt 

    # Fill the matrix of q vectors
    for i, row in tqdm(data.iterrows()):
        for j in range(0,num_pars):

            denom = min_row[par_names[j]]
            
            q_mat[i][j] = np.log(row[par_names[j]]) - np.log(denom)
    
    print('Calculating V, Q')
    for q_k in tqdm(q_mat):
        qqT = np.outer(q_k,q_k.transpose())
        V = V + float((row['cost']-min_row['cost']))*qqT
        Q = Q + np.outer(qqT,qqT)
    
    print('Done!')
    print('Calculating H:')
    print('Generating vec V')
    vecV = np.reshape(V, (num_pars**2,1))


    print('Computing L vec(V)')
    LvecV = np.matmul(L_mat,vecV)
    print(Q.shape)
    print(LvecV.shape)
    print('LQD')
    LQD = np.matmul(L_mat,np.matmul(Q,D_mat))
    print('(LQD)^-1')
    LQDinv = np.linalg.inv(LQD)
    print(LQDinv.shape)
    print('Putting it all together')
    vecH = np.matmul(D_mat,np.matmul(LQDinv,LvecV))
    print('Reshaping!')
    H = np.reshape(vecH,(num_pars,num_pars))
    print('Done!')
    H_DF = pd.DataFrame(H)
    print('Writing Hessian to file')
    print(Hname)
    H_DF.to_csv(Hname,sep='\t',header=False, index=False)
    return(num_pars, par_names)

def computeHessian_baumann(data, fname, Hpath, min_row, psetsForHessian, debug):
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

def logmessage(message):
    with open('LOGS','a') as outfile:
        today = datetime.datetime.today()
        date = str(today.year) + '-' + \
               str(today.month) + '-' + str(today.day) + \
                 '-' + str(today.hour) + '-' + str(today.minute)
        
        outfile.write(date + '\t' + message + "\n")
