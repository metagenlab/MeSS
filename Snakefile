
include : 'rules/simulate_reads.rules'
include: config['assembly_finder_Snakefile']

if config['read_status']=='single':
    sample_suffix=['1']
if config['read_status']=='paired':
    sample_suffix=['1','2']
    rule all_sim :
        input: expand('combined_{format}_shuffled/shuffled_samples_R{dir}.fq',format=config['read_status'],dir=sample_suffix)

