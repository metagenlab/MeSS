#!/bin/bash
####Functions
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
  }
##This function creates a file with random data and use it as an external random generator for shuf
##See https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html#Random-sources
seed_nb=$(( ( RANDOM % 1000 )  + 1 ))
####Parameters
input=$1
output=$2
echo $1 $2 $seed_seed_nb


####Main
awk '{OFS="\t"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}' $input | \
    shuf --random-source=<(get_seeded_random $seed_nb) | \
        awk '{OFS="\n"; print $1,$2,$3,$4}' > $output
#####