import datetime
import os, sys
import json
from pymongo import MongoClient
import click
import logging

@click.command()
@click.option('-dbout', '--output_db', default='', multiple=True, help='The database to write into.')
@click.option('-s', '--sers', default=[], multiple=True, help='Series to process.')
@click.option('-a', '--arrs', default=[], multiple=True, help='Array to process. If this is provided, then only one series is accepted.')
@click.option('-r', '--rootdir', default = '/Volumes/arrayData/arraymap/grch38/', help='Provide root directory where variants.json and callset.json should be loaded, until the level grch38/.')
@click.option('--overwrite/--no-overwrite', default=False,  help='default write into temporary collections; if flagged, overwrite current collection')
def cli(output_db, sers, arrs, rootdir, overwrite):

    """
    This script is used to update variants and callsets in the db. The following collections are affected.

    biosamples.info.cnvstatistics (currently it does not, because info.cnvstatistics doesn't exist in json file)
    variants
    callsets

    """


    def initiate_vs_cs(ser, arr):

        ## variant collections
        with open(os.path.join(rootdir, ser, arr, 'variants.json'),'r') as json_data:
            variants_json = json.load(json_data)

        variant_obj = []

        for v in variants_json:

            v.pop('no', None)
            v['info']['cnv_value'] = v['info'].pop('value')
            v['info']['cnv_length'] = v['info'].pop('svlen')
            v['info'].pop('probes', None)
            v['variantset_id'] = 'AM_VS_GRCH38'

            variant_obj.append(v)


        ## callset collections
        with open(os.path.join(rootdir, ser, arr, 'callset.json'),'r') as json_data:
            callset_json = json.load(json_data)

        callset_json.pop('callset_id', None)
        callset_obj = callset_json

        return variant_obj, callset_obj


    def retrieveTime(return_type):
        time = datetime.datetime.now()#.replace(microsecond=0).isoformat()
#        t=os.path.getmtime(filename)
#        time2 = datetime.datetime.fromtimestamp(t).isoformat()
        if return_type == 'str':
            return time.isoformat()
        elif return_type == 'date':
            return time

    def run(ser, arr):

        ## loading common variables
        biosample_id = 'PGX_AM_BS_' + arr
        callset_id = '{}::{}::{}'.format('pgxcs',ser,arr)

        variants, callset  = initiate_vs_cs(ser, arr)

        ## variants
        for variant in variants:
            variant['callset_id'] = callset_id
            variant['biosample_id'] = biosample_id
            variant['updated'] = retrieveTime('str')


        ## callsets
        callset['id'] = callset_id
        callset['updated'] = retrieveTime('date')

        # for i, j in callset.items():
        #     if i == 'info':
        #         print(i, j)

        # the rest stay the same

        client = MongoClient()
        for db_name in output_db:
            db = client[db_name]

            # ### update biosamples

            # find_all = db.biosamples.find({'id': biosample_id})
            # count_find = 0
            # for i in find_all:
            #     count_find += 1

            # if count_find > 1:
            #     logging.error('duplicated biosample {} entry in {}'.format(biosample_id ,db_name))
            #     return

            # elif count_find == 0:
            #     logging.error('no biosample entry {} in {} for editting'.format(biosample_id ,db_name))
            #     return

            # else:
            #     current_biosample = db.biosamples.find_one({'id': biosample_id})

            # print(callset['info']['cnvstatistics'])
            # current_biosample['info']['cnvstatistics'] = callset['info']['cnvstatistics']
            # current_biosample['updated'] = retrieveTime('date')

            # if overwrite:
            #     db.biosamples.delete_one({'id': biosample_id})
            #     db.biosamples.insert_one(current_biosample)
            # else:
            #     current_biosample.pop('_id', None)
            #     print(current_biosample)
            #     db.biosamples_tmp.insert_one(current_biosample)

            ## update variants

            if variants: # can be normal sample with 0 variants
                if overwrite:
                    del_log = db.variants.delete_many({'callset_id': callset_id})
                    logging.debug('deleted {} variants for callset {}'.format(
                      str(del_log.deleted_count),
                      callset_id))
                    db.variants.insert_many(variants)
                else:
                    db.variants_tmp.insert_many(variants)

            ### update callsets
            find_all = db.callsets.find({'id': callset_id})
            count_find = 0
            for i in find_all:
                count_find += 1

            if count_find > 1:
                logging.error('duplicated callset {} entry in {}'.format(callset_id ,db_name))
                return

            elif count_find == 0:
                logging.error('no callset entry {} in {} for editing'.format(callset_id ,db_name))
                return

            else:
                current_callset = db.callsets.find_one({'id': callset_id})

            for key, val in callset.items():
                current_callset[key] = val
            if overwrite:
                db.callsets.delete_one({'id': callset_id})
                db.callsets.insert_one(current_callset)
            else:
                current_callset.pop('_id', None)
                db.callsets_tmp.insert_one(current_callset)

    ## logging set-up
    today = datetime.datetime.today()
    today = today.strftime('%Y-%m-%d')
    logging.basicConfig(level = 10, filename='%s_update.log'%(today), filemode='w', format='%(name)s - %(levelname)s - %(message)s')

    ## make list of queried arrays

    if arrs:
        for arr in arrs:
            if len(sers) > 1:
                    sys.exit('more than one series for provided arrays')
            run(sers[0], arr)
    else:
        for ser in sers:
            print(ser)
            datadir = os.path.join(rootdir, ser)
            arrs = [i for i in os.listdir(datadir) if i != '.DS_Store']
            # print(arrs)
            for arr in arrs:
               run(ser, arr)

if __name__ == '__main__':
	cli()
