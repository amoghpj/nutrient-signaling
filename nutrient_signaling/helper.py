from optparse import OptionParser
import pandas as pd
import yaml

def read_data(path):
    """
    returns: pandas dataframe
    """
    df = pd.read_csv(path, index_col=0)
    return df


def convert_to_yaml(df):
    """
    returns: yaml object
    """
    yamlfile = df.to_dict(orient='index')
    return(yamlfile)


def main():
    parse = OptionParser()
    parse.add_option('-e', '--experiments', type='str',
                     help='Path to experimental data, csv file')
    parse.add_option('-o', '--out-prefix', type='str',
                     help='Output Prefix')
    opts, args = parse.parse_args()

    path = opts.experiments
    outPrefix = opts.out_prefix
    y = convert_to_yaml(read_data(path))
    
    with open(outPrefix + '-experiments.yaml', 'w') as outfile:
        outfile.write(yaml.dump(y))

if __name__ == '__main__':
    main()

