
include : "rules/simulate_reads.rules"

include: config['assembly_finder_Snakefile']

nb_rep=config["replicates"]

if config['read_status']=='single':
    read_direction=['1']
if config['read_status']=='paired':
    read_direction=['1','2']

def rep_expand(nb_rep):
    output_files=expand('simreads/{community}/simrep_{rep}/R{rd}.fq',rep=list(range(1,nb_rep+1)),\
                        community=config["community_name"],rd=read_direction)
    return output_files

rule all_sim :
    input: rep_expand(nb_rep)

