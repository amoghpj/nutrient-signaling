import nutrient_signaling.modelreader as md
from optparse import OptionParser

def modelToLaTeX(modelpath, writepath="model.tex"):
    model = md.readinput(modelpath)
    latexstring = "\\begin{align}\n"
    variables = model['variables'].keys()
    parameters = model['parameters'].keys()
    for variable, equation in model['variables'].items():
        for v in variables:
            equation = equation.replace(v, "[\\text{" + v + "}]")
            equation = equation.replace(v+"_T", "[\\text{" + v + "_T}]")
        equation = equation.replace("shs", "\\mathcal{H}")
        equation = equation.replace("gamma_", "gamma")
        equation = equation.replace("w_", "\\omega_")
        equation = equation.replace("sigma_", "\\sigma_")                
        equation = equation.replace("*", " ")
        equation = equation.replace("+", " + ")
        equation = equation.replace("-", " - ")
        equation = equation.replace(",", " , ")                        
        equation = equation.replace("gamma", "\\gamma_")

        
        atoms = equation.split(" ")
        corrected = []
        for a in atoms:
            if '_' in a:
                a = fixUnderscore(a)

            corrected.append(a)
        equation = " ".join(corrected)
            
        latexstring += "\\frac{d[\\text{"+ variable+"}]}{dt} &= " + equation + "\\\\\n"
    latexstring += "\\end{align}"
    with open(writepath, 'w') as outfile:
        outfile.write(latexstring)

def fixUnderscore(a):
    asplit = a.split("_")
    afixed = [asplit[0]]
    afixed.extend(["{" +asplit[i+1] + "}" for i in range(len(asplit) - 1)])
    a = "_".join(afixed)
    return a

def main():
    parser = OptionParser()
    parser.add_option('-m','--model-path',dest='model',default="",
                  help='Path to model definition.')
    parser.add_option('-w','--write-path',dest='write',default="model.tex",
                  help='Path to write model.')    
    (options,args)=parser.parse_args()
    modelToLaTeX(options.model, writepath=options.write)
    
if __name__=='__main__':
    main()

    
