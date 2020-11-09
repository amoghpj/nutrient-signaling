import pandas as pd
from decimal import Decimal

outfile = "output/parameter-latex-tables.tex"
parametersetpath = "data/2018-9-26-12-3-no-sigma/parameters.txt"

pars = pd.read_csv(parametersetpath, index_col=None, header=None,sep='\t')

tex = "\\begin{table}\n"

step = int(pars.shape[0]/3.0)

header = "\\begin{minipage}{0.33\\textwidth}\n"\
    "\\begin{tabular}{|p{2.5cm}|p{2cm}|}\n"\
    "\\hline\n"\
    "Name & Value \\\\\n"\
    "\\hline\n"

footer = "\\hline\n"\
    "\\end{tabular}\n"\
    "\\end{minipage}\n"\

tex += header

counter = 0
while counter < pars.shape[0]:
    name = pars.loc[counter, 0].replace("_","$\\_$")
    value = "%0.2e" % Decimal(pars.loc[counter, 1])
    tex +=  name + " & " + str(value) 
    if counter % step == 0 and counter >0:
        tex += "\n"        
        tex += footer
        tex += header

    else:
        tex +=  "\\\\\n"
    counter += 1
tex += "\\end{table}"
with open(outfile, "w") as out:
    out.write(tex)
