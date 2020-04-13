import datetime
import sys
import json
import pandas as pd
from pymongo import MongoClient
import click
import random



@click.command()
@click.option('-dbout', '--output_db', default='', multiple=True, help='The database to write into.')
@click.option('-m', '--metapath', default='', help='The path to the metadata table.')
@click.option('-d', '--demo', default=0, type=click.IntRange(0, 10000), help='Enter a small number of entries to process as demo.')
def cli(output_db, metapath, demo):

    """
    This script takes metadata table and generates attributes for 4 collections to insert into the output database.

    biosamples
    variants
    callsets
    individuals

    """


    def initiate_vs_cs(ser, arr):
       
        ## variant collections
        with open('/Volumes/arrayMaster/arraymap/grch38/{0}/{1}/variants.json'.format(ser,arr),'r') as json_data:
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
        with open('/Volumes/arrayMaster/arraymap/grch38/{0}/{1}/callset.json'.format(ser,arr),'r') as json_data:
            callset_json = json.load(json_data)

        callset_json.pop('callset_id', None)
        callset_obj = callset_json

        return variant_obj, callset_obj

    def retrieveTime(filename, return_type):
        time = datetime.datetime.now()#.replace(microsecond=0).isoformat()
#        t=os.path.getmtime(filename)
#        time2 = datetime.datetime.fromtimestamp(t).isoformat()
        if return_type == 'str':
            return time.isoformat()
        elif return_type == 'date':
            return time


    def retrieveGEO(city):
        client = MongoClient()
        cl = client['progenetix'].geolocs
        geo_obj = cl.find_one({'CITY':city})
        if geo_obj:
            geoPrecision = 'city'
            country = geo_obj['COUNTRY']
            geoLabel = city + ', ' + country
            [geolong, geolat] = geo_obj['GEOLOC']['coordinates']
            geoLabel = city + ', ' + country

        else:
            geoLabel = None
            geoPrecision = None
            geolat = None
            geolong = None
            city = None
            country = None

        geo_info = { 'label': geoLabel,
                'precision': geoPrecision,
                'city': city,
                'country': country,
                'latitude': geolat,
                'longitude': geolong}

        return geo_info

    def retrievePlatformLabel(platformID):
        client = MongoClient()
        cl = client['arraymap'].platforms
        plf_obj = cl.find_one({'PLATFORMID':platformID})
        try:
            return plf_obj['PLATFORMLABEL']
        except:
            print('{} has no platform description in arraymap.platforms'.format(platformID))


    ### read in meta table
    mytable = pd.read_csv(metapath, sep = '\t', dtype = {'YEAR': str, 'PMID': str, 'SEERCODE': str})
    mytable = mytable.where((pd.notnull(mytable)), None) ## convert pd.nan to None
    no_row = mytable.shape[0]
    if demo != 0:
        bar_length = demo
        rdm_row = random.sample(range(no_row), demo)
        mytable = mytable.iloc[rdm_row, :]
    else:
        bar_length = no_row

    ### define list/dictionary of objects to insert in 4 collections
    variants_list = []
    callsets_dict = {}
    biosamples_dict = {}
    individuals_dict = {}

    ### find all existing ids in each output database and collections.
    exist_callset_id = {}
    exist_biosample_id = {}
    exist_individual_id  = {}

    client = MongoClient()
    for db_name in output_db:
        exist_callset_id[db_name] = client[db_name].callsets.distinct('id')
        exist_biosample_id[db_name] = client[db_name].biosamples.distinct('id')
        exist_individual_id[db_name] = client[db_name].individuals.distinct('id')

    with click.progressbar(mytable.itertuples(), label='Processing',
        fill_char=click.style('*', fg='green'), length=bar_length) as bar:

        for row in bar:
            ### read raw information from table
            ser = row.SERIESID
            arr = row.UID
            plf = row.PLATFORMID
            icdt = row.ICDTOPOGRAPHY
            icdtc = row.ICDTOPOGRAPHYCODE
            icdm = row.ICDMORPHOLOGY
            icdmc = row.ICDMORPHOLOGYCODE.replace('/','')
            diagnosis = row.DIAGNOSISTEXT
            if not diagnosis:
                diagnosis = ''
            age = row.AGE
            sex = row.SEX
            pmid = row.PMID
            tnm = row.TNM
            samplesource = row.SAMPLESOURCE
            ncitcode = row.NCITCODE
            if not ncitcode:
                ncitcode = 'None'
            ncitterm = row.NCITTERM
            seercode = row.SEERCODE
            city = row.CITY
            cellline_id = row.CELLLINEID
            cellosaurus_id = row.CELLOSAURUSID

            ### derived attributes that are shared by collections
            biosample_id = 'PGX_AM_BS_' + arr
            callset_id = '{}::{}::{}'.format('pgxcs',ser,arr)
            individual_id = 'PGX_IND_' + arr
            material_id = 'EFO:0009654' if icdmc == '00000' else 'EFO:0009656'
            material_label = 'reference sample' if icdmc == '00000' else 'neoplastic sample'
            geo_provenance = retrieveGEO(city)
            data_use = {
                    'label' : 'no restriction',
                    'id' : 'DUO:0000004'
                    }

            ############################
            ##   variants & callsets  ##
            ############################
            variants, callset  = initiate_vs_cs(ser, arr)

            ### check if callset_id exists already in the dababase and in the current process.

            check_callset = [callset_id in exist_callset_id[db_name] for db_name in output_db]
            if all(check_callset):
                continue ### if callset exists then the sample shouldn't be processed.

            ## variants
            for variant in variants:
                variant['callset_id'] = callset_id
                variant['biosample_id'] = biosample_id
                variant['updated'] = retrieveTime(metapath, 'str')

            variants_list.append(variants)


            ## callsets
            callset['id'] = callset_id
            callset['biosample_id'] = biosample_id
            callset['updated'] = retrieveTime(metapath, 'date')
            callset['description'] = retrievePlatformLabel(plf)
            callset['data_use_conditions'] = data_use

            callsets_dict[callset_id] = callset

        	############################
            ######   biosamples  #######
            ############################

            biosample= {
                    'id': biosample_id,
                    'description': diagnosis,
                    'biocharacteristics':[
                            {
                                    'description': diagnosis,
                                    'type': {
                                            'id':  'icdot:' + icdtc,
                                            'label': icdt
                                     }
                                    },
                            {
                                    'description': diagnosis,
                                    'type': {
                                            'id':  'icdom:' + icdmc,
                                            'label': icdm
                                     }
                                    },
                            {
                                    'description': diagnosis,
                                    'type': {
                                            'id':  'ncit:' + ncitcode.capitalize(),
                                            'label': ncitterm
                                     }
                                    }],
                    'updated': retrieveTime(metapath, 'date'),
                    'individual_id': individual_id,
                    'project_id': ser,
                    'age_at_collection':
                        {
                                'age_class': ### to be added later
                                    {
                                            'id': None,
                                            'label': None
                                        },
                                'age': age
                                    },
                    'external_references':[
                            {
                                    'type': {
                                            'id': 'geo:' + arr,
                                            'label': ''},
                                    'description': 'geo:gsm',
                                    'relation': 'denotes'
                                    },
                            {
                                    'type': {
                                            'id': 'geo:' + ser,
                                            'label': ''},
                                    'description': 'geo:gse',
                                    'relation': 'denotes'
                                    },
                            {
                                    'type': {
                                            'id': 'geo:' + plf,
                                            'label': ''},
                                    'description': 'geo:gpl',
                                    'relation': 'denotes'
                                    }
                            ],
                    'provenance':{
                            'material':{
                                    'description': diagnosis,
                                    'type':{
                                            'id': material_id,
                                            'label': material_label
                                            }
                                    },
                            'geo': geo_provenance
                            },
                    'info':{
                            'death': None,
                            'followup_months': None,
                            'samplesource': samplesource
                            },
                    'data_use_conditions': data_use

                    }

            if tnm:
                biosample['info']['tnm'] = tnm

            if seercode:
                biosample['info']['seer'] = seercode

            if pmid:
                biosample['external_references'].append(
                        {
                                    'type': {
                                            'id': 'pubmed:' + pmid,
                                            'label': ''},
                                    'description': 'pubmed',
                                    'relation': 'denotes'
                                    }
                        )

            if cellosaurus_id:
                biosample['external_references'].append(
                        {
                                    'type':{
                                            'id': 'cellosaurus:' + cellosaurus_id,
                                            'label': ''
                                            },
                                    'description': 'cellosaurus',
                                    'relation': 'denotes'

                                })

            if cellosaurus_id:
                biosample['info']['cell_line'] = cellline_id

            biosamples_dict[biosample_id] = biosample

            ############################
            ######   individuals  ######
            ############################

            if sex:
                sex_id = 'PATO:0020001' if sex == 'male' else 'PATO:0020002'
                sex_label = 'male genotypic sex' if sex == 'male' else 'female genotypic sex'
            else:
                sex_id = sex_label = None

            individual = {
                    'id': individual_id,
                    'description': None,
                    'biocharacteristics':[
                            {
                                    'description': sex,

                                    'type': {
                                            'id': sex_id,
                                            'label': sex_label}
                                    },
                            {
                                    'description': None,
                                    'type':{
                                            'id' : 'NCBITaxon:9606',
                                            'label' : 'Homo sapiens'
                                    }
                            }],
                     'data_use_conditions' : data_use,
                     'geo_provenance': geo_provenance,
                     'updated': retrieveTime(metapath, 'date')
                    }

            individuals_dict[individual_id] = individual

    ############################
    ###   database write-in  ###
    ############################

    click.echo()
    if click.confirm("""I have processed {} variants, {} callsets, {} biosamples and {} individuals for update.
Do you want to continue?""".format(sum([len(v) for v in variants_list]), len(callsets_dict), len(biosamples_dict),
                     len(individuals_dict)), default = True):
        client = MongoClient()
        for db_name in output_db:
            db = client[db_name]

            for variant_obj in variants_list:
                try:
                    db.variants.insert_many(variant_obj)
                except TypeError:
                    pass

            for callset_id, callset_obj in callsets_dict.items():
                if callset_id in exist_callset_id[db_name]:
                    continue
                db.callsets.insert_one(callset_obj)

            for biosample_id, biosample_obj in biosamples_dict.items():
                if biosample_id in exist_biosample_id[db_name]:
                    continue
                db.biosamples.insert_one(biosample_obj)

            for individual_id, individual_obj in individuals_dict.items():
                if individual_id in exist_individual_id[db_name]:
                    continue
                db.individuals.insert_one(individual_obj)

if __name__ == '__main__':
	cli()
