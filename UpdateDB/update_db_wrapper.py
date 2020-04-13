## update specific series and arrays

import subprocess as sbp
from pymongo import MongoClient
import logging, os
from datetime import datetime

now = datetime.now()
dir = 'logs'
if not os.path.exists(dir):
    os.makedirs(dir)

logging.basicConfig(filename='logs/%s_update_callset_variants.log' % now.strftime("%Y%m%d-%Hh%M"),level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
client = MongoClient()
db_names = ['arraymap','progenetix']
overwrite = "--overwrite" #''
for db_name in db_names:

    db = client[db_name]
    for item in db.callsets.find({'info.cnvstatistics':{"$exists":0}}):
        ser, arr = item['id'].split('::')[1:3]

        p = sbp.run('python3 /Users/Shared/ProgenetixTools/update_db.py -dbout {} -s {} -a {} {}'.format(db_name, ser, arr, overwrite), capture_output=True, shell = True)
        if p.returncode == 0:
            tmp = 'temporary' if overwrite == '' else 'original'
            logging.info('{} {} written to {} in {} collections'.format(ser, arr, db_name, tmp))

