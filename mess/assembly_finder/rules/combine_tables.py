import pandas as pd
"""
Main
"""
column = snakemake.params['column']
df_list = []
for file in snakemake.input:
    entry = file.split('/')[1].split('-filtered.tsv')[0]
    tb = pd.read_csv(file, sep='\t')
    tb.insert(loc=0, column=f'{column}', value=[entry]*len(tb))
    df_list.append(tb)
df = pd.concat(df_list, sort=False)
df.to_csv(snakemake.output[0], sep='\t', index=None)
