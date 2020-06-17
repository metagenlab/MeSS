include : "rules/simulate_reads.rules"
include : "assembly_finder/Snakefile"

nb_rep=config["replicates"]

if config['read_status']=='single':
    read_direction=['1']
if config['read_status']=='paired':
    read_direction=['1','2']

def rep_expand(nb_rep):
    output_files=expand('simreads/{community}/simrep-{rep}/{community}-{rep}_R{rd}.fq.gz',rep=list(range(1,nb_rep+1)),\
                        community=config["community_name"],rd=read_direction)
    return output_files



rule all_sim :
    input: rep_expand(nb_rep),
           expand('krona/{community}/simrep-{rep}.html',rep=list(range(1,nb_rep+1)),community=config["community_name"])