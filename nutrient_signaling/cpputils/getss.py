import os
import pandas as pd
import matplotlib.pyplot as plt
print(os.path.realpath('.'))
PATH = "./"
plotvar = "cAMP"
cmd = ''
cmd += PATH + 'main.o --ode --tend 90 --solver rk4 --out initial.dat --step 0.025 --pars  Carbon 1.0 Glutamine_ext 1.0 ATP 1.0'
so = os.popen(cmd).read()
print(so)
Res = pd.read_csv(PATH + 'initial.dat',index_col=False, sep='\t')
plt.plot(Res['t'], Res[plotvar],'b')
SS = Res.tail(1)
cmd = ''
cmd+='./main.o --ode --tend 90 --solver rk4 --step 0.025 --out final.dat --pars Carbon 0.0 Glutamine_ext 0.0 ATP 0.0 --ics '

for h in list(Res.columns.values):
    if h!='t':
        cmd += h
        cmd += ' '
        cmd += str(float(SS[h]))
        cmd += ' '


so = os.popen(cmd).read()
print(so)
Res1 = pd.read_csv(PATH + 'final.dat',index_col=False, sep='\t')
plt.plot(float(SS['t']) + Res1['t'], Res1[plotvar],'r')
plt.show()


