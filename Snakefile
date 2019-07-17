include : 'rules/summary_table/summary_table.rules'
#include : 'rules/download/download.rules'
include : 'rules/simulate_reads/simulate_reads.rules'

if config['read_status']=='single':
    sample_suffix=['1']
if config['read_status']=='paired':
    sample_suffix=['1','2']
    rule all :
        input: expand('combined_{format}_shuffled/shuffled_samples_R{dir}.fq',format=config['read_status'],dir=sample_suffix)


