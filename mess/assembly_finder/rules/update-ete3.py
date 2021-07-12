import sys
import os
from datetime import datetime
from ete3 import NCBITaxa
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    ncbi = NCBITaxa()
    sqldb = os.path.join(os.environ['HOME'], '.etetoolkit', 'taxa.sqlite')  # path to ete sql db
    db_modification_time = datetime.fromtimestamp(os.path.getctime(sqldb))
    database_age_days = abs((db_modification_time-datetime.now()).days)
    if database_age_days >= 10:
        ncbi.update_taxonomy_database()
        comment = f'taxa.sqlite is more than {database_age_days} days old, updating database'
    else:
        comment = 'taxa.sqlite is up to date'
    file = open(snakemake.output[0], 'w')
    file.write(comment)
    file.close()






