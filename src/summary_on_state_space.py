import pandas as pd
import numpy as np
import confidence_state_space as css
import sys

def enclose_in_tabular(s, numcols):
    return('\\begin{tabular}' +
           '|'.join(['c' for _ in range(len(numcols))]) +
           '\n' +
           s +
           '\\end{tabular}\n')
def main():
    settings = css.Settings()
    
    ## State Predictions:
    ## This tab-separated file is generated using 
    ## src/compare-state-space-predictions.py
    predictionpath = 'output/global_space_24066.csv'
    predictiondf = pd.read_csv(predictionpath,
                                sep='\t', index_col=False,
                                dtype='str')
    ## Filter by cost, retain anything less than twice the minimum
    predictiondf = predictiondf.astype({'cost':'float'})
    cmin = predictiondf['cost'].min()
    predictiondf = predictiondf[predictiondf['cost'] <= 2.*cmin]

    # Curated evidence for states:
    ## This comma-separated file is manually curated,
    ## with references to the original publication
    evidencepath = 'data/csv/tf-state-experimental-evidence.csv'
    evidencedf = pd.read_csv(evidencepath, index_col=None)
    evidencedict = {row['strain']:row['state']                     # Create dict {strain:state}
                       for ind, row in                                # Loop over key, value of 
                       evidencedf.dropna(how='any').T.to_dict().items()} # evidence table converted to dict
    # Put wt on top
    colorder = [c for c in reversed(list(predictiondf.columns))]

    colorder.remove('cost')

    seen = set()
    allstrains = []
    for c in colorder:
        strain = c.split('_')[0]
        if strain not in seen:
            allstrains.append(strain)
            seen.add(strain)
            
    # make a nested dictionary to hold predictions
    confidence = [{'name':strain,
                   'states':{nstate:{'off':0.0,'on':0.0} for nstate in
                             settings.nutrientStates}}
                  for strain in allstrains]
    strainMapper = {conf['name']:i for i, conf in enumerate(confidence)}
    print(strainMapper)
    # Record predictions
    ## Each entry in predictiondf is a string of 0s and 1s
    ## which records the state of the TFs in the order given
    ## by settings.readouts.  The columns in this file take
    ## the form STRAIN-NAME_NUTRIENT-CONDITION.
    ## The function compare_states() reads the list containing
    ## predictions from each parameter set, and for each TF,
    ## counts the number of time the parameter sets predicted
    ## an ON vs an OFF. This is reported as a fraction of the total
    ## number of predictions, and is returned as a list. (See docs for visualize())
    
    for col in colorder:
        strain, nstate = col.split('_')
        strainid = strainMapper[strain]
        confidence[strainid]['states'][nstate] = css.compare_states(list(predictiondf[col]), settings)
    numpsets = predictiondf.shape[0]
    strains = [conf['name'] for conf in confidence]
    strainnames = []
    for strain in strains:
        if strain == 'wt':
            strainnames.append('wt')
        else:
            strainnames.append(settings.mapper[strain])
    summarydf = pd.DataFrame(columns=strainnames, index=pd.Index(settings.nutrientStates))
    revmapper = {s:r for r,s in settings.mapper.items()}
    revmapper['wt'] = 'wt'
    
    numcols = len(settings.readouts) + 1
    # ltstr = enclose_in_tabular("",numcols)
    
    for s in strainnames:
        current = None
        for c in confidence:
            if c['name'] == revmapper[s]:
                current = c
        for n in settings.nutrientStates:
            summarydf[s].loc[n] = ''
            for tf in current['states'][n]:
                summarydf[s].loc[n] += str(round(100*tf['on'],2)) + ','
    summarydf = summarydf.T
    summarydf = summarydf.iloc[::-1]    
    write_latex(summarydf, settings, evidencedict, revmapper)
    write_csv(summarydf, settings, evidencedict, revmapper)

def write_org(summarydf, settings):
    # OLD: write each nutrient state as separate table to a single org file
    lt= ""
    for ns in summarydf.columns:
        lt += '#+NAME: ' + ns + '\n'
        lt += '#+CAPTION: ' + ns + '\n'        
        lt += "|Strain|" + '|'.join(settings.readouts) + '|\n'
        for strain, row in summarydf.iterrows():
            vals = [float(s) for s in row[ns].split(',')]
            form = ''
            lt += '|' +  strain.replace('$',' ') + '|' + row[ns].replace(',', '|')  + '\n'
        lt += '\n\n'
    
    with open('summary_state_space_per_nutrient.org','w') as outfile:
        outfile.write(lt)

def write_csv(summarydf, settings, evidencedict, revmapper):
    columns = ["Strains"]
    for ns in summarydf.columns:
        for r in settings.readouts:
            columns.append(ns + '-' + r)
    expandeddf = pd.DataFrame(columns=columns, index=summarydf.index)
    for strain, row in summarydf.iterrows():
        for ns in summarydf.columns:
            vals = [float(v) for v in row[ns].split(',') if v != ""]
            for r, v in zip(settings.readouts, vals):
                expandeddf[ns + '-' + r].loc[strain] = v
    expandeddf.to_csv("summary_state_space.csv")
        
def write_latex(summarydf, settings, evidencedict, revmapper):
    lateximports = "\\documentclass[landscape]{article}\n"\
        "\\usepackage[a3paper, margin=1in]{geometry}\n"\
        "\\usepackage[table,dvipsnames]{xcolor}\n"\
        "\\usepackage{longfbox}\n"\
        "\\definecolor{robustgreen}{HTML}{D4EFDF}\n"\
        "\\definecolor{robustred}{HTML}{F5B7B1}\n"\
        "\\begin{document}\n"
    latexfooter = "\\end{document}"
    
    tableheader = lambda nsstart, nsend: "\\begin{table}\n\\centering\n"\
        "\\begin{tabular}{|" + '|'.join(['c' for c in range((nsend-nsstart)*6 + 1)])+ "|}\n" +\
        "\\hline\n" +\
        "\\textbf{Strain}  & " +\
        " & ".join(["\\multicolumn{6}{c|}{\\textbf{" + ns + "}}"\
                    for ns in settings.nutrientStates[nsstart:nsend]])+ "\\\\\n" +\
        "& " + " & ".join([" & ".join(["\\textbf{"+ r + "}"\
                                                        for r in settings.readouts])\
                           for _ in range(nsend-nsstart)]) + "\\\\\n"+\
        "\\hline \n"
    tableclose = "\\hline\n\\end{tabular}\n\\end{table}\n"
    print(settings.readouts)
    
    lt = lateximports 

    ##############################
    numcols = 4                #  Aesthetic blocks of 4 nutrient conditions 
    numtables = int(8/numcols) #  seems to look good for now
    start = 0                  #
    ##############################

    # main loop over split tables
    for tabid in range(numtables):
        pertable = summarydf[list(summarydf.columns[start:start + numcols])]
        lt += tableheader(start,start + numcols)
        for strain, row in pertable.iterrows():
            pernut = []            
            for ns in pertable.columns:
                vals = [int(np.ceil(float(s))) for s in row[ns].split(',') if s != ""]
                for counter, v in enumerate(vals):
                    cellcolor, textcolor = estimate_robustness(v)
                    borcol, bwidth, style = experiment_agreement(v, revmapper[strain], ns, counter, settings,evidencedict)
                    pernut.append("\\lfbox[border-color="+borcol+",border-width="+bwidth+",border-style="+style+"]{\\cellcolor{" + cellcolor+ "}\\textcolor{"+textcolor+"}{" + str(v) + "}}")
            lt += strain + "&" + " & ".join(pernut) + "\\\\\n"
        lt += tableclose
        start += numcols
    lt += latexfooter
    with open('summary_state_space.tex','w') as outfile:
        outfile.write(lt)
        
def experiment_agreement(v, strain, nstate, readout_counter, settings, evidencedict):
    strain_nstate_readout = strain + '_' + nstate + '_' + settings.readouts[readout_counter]
    if strain_nstate_readout in evidencedict.keys():
        gtstate = evidencedict[strain_nstate_readout]
        color = 'red'
        if gtstate == 'ON':
            color = 'green'
        predstate = 'ON'
        ## NOTE Model prediction matches the experiment if the
        ## majority prediction matches.
        if (100-v) > v:
            predstate = 'OFF'
        if gtstate == predstate:
            return(color, "3pt", "solid")
        else:
            return(color, "3pt", "dashed")
    return("black", "0pt", "solid")            
    
def estimate_robustness(v):
    if min(v, 100.-v) < 10.:                
        if v > 100. - v:
            cellcolor = "robustgreen"
        else:
            cellcolor = "robustred"
        textcolor = "black"
    else:
        if v > 100. - v:
            cellcolor = "ForestGreen"
        else:
            cellcolor = "red"
        textcolor = "white"            
    return(cellcolor, textcolor)

    
if __name__ == '__main__':
    main()

