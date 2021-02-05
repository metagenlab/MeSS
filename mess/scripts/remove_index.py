import pandas as pd
checkpoint_output = snakemake.input[0]
table=pd.read_csv(checkpoint_output,delimiter='\t',index_col=1)
new_tb=table.drop(columns='AssemblyNames')
new_tb.to_csv(snakemake.output[0],sep='\t',header=None,index=True)


