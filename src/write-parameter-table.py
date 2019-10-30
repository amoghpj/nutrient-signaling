from nutrient_signaling.utils import *
import nutrient_signaling.modelreader as md
import os
home = os.path.expanduser('~')

def make_parameter_values_file(modeldefinition, outpath):
    replace_things = ['(',')','*','^','+','-',' ','shs',',']
    modelinputs = ['Glutamine_ext','NH4','Carbon','Proline']
    model = modeldefinition
    sorted_variables = sorted(model['variables'].keys())
    parametertable = []
    for v in sorted_variables:
        s = model['variables'][v]
        for r in replace_things:
            s = s.replace(r, '|')
        s = s.split('|')
        s = [t for t in s \
             if t != '' \
             and t not in sorted_variables\
             and t not in modelinputs]
        found = [(p,round(val,2),v)\
                 for p, val in model['parameters'].items()\
                 if p in s]
        parametertable.extend(found)
        if v == 'Glutamine':
            parametertable.extend([(mi,model['parameters'][mi],'Glutamine') for mi in modelinputs])
    
    df = pd.DataFrame(parametertable, columns=['Parameter', 'Value','Variable'])
    df.to_csv(outpath)

model = md.readinput(home + '/jalihal_projects/Research/data/ModelAnalysis/ParameterSets/2018-9-26-12-3-no-sigma')
outpath = home + '/group/amogh-jalihal/papers/Nutrient-Signaling-Model/supplementary-files/parameter-values.csv'

# Function Call
make_parameter_values_file(model, outpath)
