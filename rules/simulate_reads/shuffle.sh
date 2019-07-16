#!/bin/bash
####Functions

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
  }
####Parameters
input=$1
output=$2
echo $1 $2 $3



####Main
awk '{OFS="\t"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}' $input | \
    shuf --random-source=<(get_seeded_random $3) | \
        awk '{OFS="\n"; print $1,$2,$3,$4}' > $output
#####