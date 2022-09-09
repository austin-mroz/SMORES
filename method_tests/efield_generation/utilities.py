import pandas as pd

def read_xyz(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    df = pd.DataFrame()
    for line in lines[2:]:
        line = list(filter(None, line.removesuffix('\n').split(' ')))
        nrow = pd.DataFrame([[str(line[0]), float(line[1]), float(line[2]), float(line[3])]])
        tempdf = pd.concat([df, nrow])
        df = pd.DataFrame(tempdf)
    df.columns = ['atom', 'x', 'y', 'z']
    
    return df


