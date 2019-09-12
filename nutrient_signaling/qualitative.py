class QualExp:
    def __init__(self):
        self.readouts = ['Gis1',
                         'Mig1',
                         'Dot6',
                         'Gcn4',
                         'Rtg13',
                         'Gln3']
        self.cutoffs = {
            'Gis1':0.5,
            'Gcn4':0.28,
            'Rtg13':0.3,
            'Gln3':0.5,
            'Mig1':1.325,
            'Dot6':1.85}

    def read_data(self, path_to_yaml):
        with open(path_to_yaml,'r') as infile:
            self.data = yaml.safe_load(infile)
            
    def list_experiments(self):
        for exp in self.data['experiments']:
            print("%s\t%s" %(exp['exp_id'],exp['name']))

    def describe_experiment(self, exp_id):
        for exp in self.data['experiments']:
            if exp['exp_id'] == exp_id:
                expoi = exp
                break
            
        print('ID\t\t%s' % (exp.get('exp_id',"")))
        print('NAME\t\t%s' % (exp.get('name',"")))
        print('PRESHIFT ICS\t\t%s' % (str(exp.get('pre_ics',""))))
        print('POSTSHIFT_ICS\t\t%s' % (str(exp.get('post_ics',""))))
        print('PRESHIFT PARS\t\t%s' % (str(exp.get('pre_pars',""))))
        print('POSTSHIFT_PARS\t\t%s' % (str(exp.get('post_pars',""))))
        print('MUTANT_SPEC\t\t%s' % (str(exp.get('mutant',""))))                

    def simulate_experiment():
        """
        TODO In cases where parameter changes are relative to the 
        defaut values, handle string parsing. 
        Ex. in 14-ure2-2x-rich-rap, 'w_gln3'='2X' should 
        be implemented correctly.
        """
        return
    
    def tfStates(self, ss,readouts,cutoffs):
        for rd in readouts:
            if rd in ['Mig1','Dot6']:
                val = np.log10((ss[rd])/(1-ss[rd]))
            else:
                val =  ss[rd]
            if val > cutoffs[rd]:
                dec = 'On'
            else:
                dec = 'Off'
            print("%s\t%0.3f\t%s" %(rd, val,dec))

    def plotTfs(self, P, whatsimulation):                                     
        plt.close()                                                               
        f,ax = plt.subplots(1,3,figsize=(6,2))                                    
        for rd in self.readouts:                                                       
            if rd == 'Mig1':                                                      
                ax[1].plot(P['t'], [np.log10((p)/(1-p)) for p in P[rd]], label=rd)
            elif rd == 'Dot6':                                                    
                ax[2].plot(P['t'], [np.log10((p)/(1-p)) for p in P[rd]], label=rd)
            else:                                                                 
                ax[0].plot(P['t'], P[rd], label=rd)                               
        ax[1].set_ylim([1.1,1.6])                                                 
        ax[2].set_ylim([1.2,2.8])                                                 
        ax[0].set_ylim([0.0,1.0])                                                 
        ax[0].legend(fancybox=False, framealpha=0.0)                              
        ax[1].legend(fancybox=False, framealpha=0.0)                              
        ax[2].legend(fancybox=False, framealpha=0.0)                              
        plt.tight_layout()                                                        
        plt.suptitle(whatsimulation)                                              
        plt.savefig(whereami + '/img/' + whatsimulation + '.png')
