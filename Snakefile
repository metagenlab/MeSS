include : 'rules/summary_table/summary_table.rules'
#include : 'rules/download/download.rules'
include : 'rules/simulate_reads/simulate_reads.rules'

rule all :
    input: 'simulated_reads/shuffled/simulated.fq'

