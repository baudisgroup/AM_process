import os,sys
import nltk

inputpath = sys.argv[1] # /Users/pgweb/arraydata/aroma/serarr_need_anno2.txt

def getWordFromPrefix (prefixes, tokens):
    matches = set()
    for w in tokens:
        for prefix in prefixes:
            if w.startswith(prefix):
                matches.add(w)
    return(matches)

def getBeforeAfter(inputword):
    flag = 1
    diag_before = ''
    previous = [all_words[offset-5:offset] for offset in c.offsets(inputword)]
    # print(previous)
    # print(len(previous))
    while isinstance(previous[0],list):
        previous=previous[0]

    for i in range(4,-1,-1):
        # print(previous[i])
        if previous[i] in puncs:
            flag = 0
        if flag:
            diag_before = previous[i] + ' ' + diag_before
    flag = 1
    diag_after = ''
    latter = [all_words[offset:offset+5] for offset in c.offsets(inputword)]
    while isinstance(latter[0],list):
        latter = latter[0]
    for i in range(5):
        if latter[i] in puncs:
            flag = 0
        if flag:
            diag_after += ' ' + latter[i]

    diagnosis = diag_before+' '+diag_after
    return(diagnosis)

def getDescriptionfromLine(lines, linestart):
    descriptionlist = []
    for l in lines:
        if l.startswith(linestart):
            descriptionlist.append(l.partition(' = ')[2].rstrip())
    description = ','.join(descriptionlist)
    return(description)


def parseStage(word):
    stage = None
    if word.startswith('I'):
        stage = 'T'+str(word.count('I'))
    if word[0] in list('1234'):
        stage = 'T'+word[0]
    if word.startswith('T'):
        try:
            if word[1] in list('1234'):
                stage = word
        except:
            pass
    return(stage)

def findGender(keywords):
    for keyword in keywords:
        words = [all_words[offset+2] for offset in c.offsets(keyword)]
        for word in words:
            if word in ['M','Male','male']:
                return('male')
            if word in ['F','Female','female']:
                return('female')

def findAge(keyword):
    age_final = None
    age = [all_words[offset+1:offset+6] for offset in c.offsets(keyword)]
    for found in range(len(age)):
        if not age_final:
            ## search next 5 words after age keyword
            for i in range(5):
                try:
                    ## if the word is a number > 1, presume it's an age, note the position
                    a = float(age[found][i])
                    if a > 1:
                        age_final = round(a)
                        found_final = found
                        i_final = i
                except:
                    continue
    if age_final:
        # print(age,found_final,i_final)
        try:
            if age[found_final][i_final+1] in ['months','month']:
                age_return = 'P'+str(age_final)+'M'
            else:
                age_return = 'P'+str(age_final)+'Y'
        except:
            age_return = 'P'+str(age_final)+'Y'
    else:
        age_return = ''
    return(age_return)

## get input text from geometa.soft by querying series array path
with open(inputpath,'r') as f:
    content = f.readlines()
with open('done.log','w+') as f:
    done = f.readlines()
with open('notFound.log','w+') as f:
    notFound = f.readlines()
doneserarrs = []
for l in done:
    doneserarrs.append(l.rstrip())
for l in notFound:
    doneserarrs.append(l.rstrip())

serarrs = []
for l in content:
    if l.startswith('GSE'):
        serarr = l.rstrip()
        if serarr not in doneserarrs:
            serarrs.append(serarr)

serarrs.sort()
for count,serarr in enumerate(serarrs):
    ### get general information from series metadata
    ## get IDs for the sample
    gseID = serarr.split('/')[0]
    gsmID = serarr.split('/')[1]
    try:
        with open(os.path.join('/Volumes/arraymapIncoming/GEOmeta/',gseID,'geometa.soft'),'r') as f:
            text = f.read()
    except:
        continue
    lines = nltk.line_tokenize(text)
    series_title = getDescriptionfromLine(lines,'!Series_title')
    bioproject = getDescriptionfromLine(lines, '!Series_relation = BioProject')
    bioproject = bioproject.split('/')[-1]

        # print(serarr)
    try:
        with open(os.path.join('/Volumes/arraymapIncoming/GEOmeta/',serarr,'geometa.soft'),'r') as f:
            if text:
                text += f.read()
            else:
                text = f.read()
    except:
        with open('notFound.log','a') as wf:
            wf.write(serarr+'\n')
        continue
    all_words = nltk.word_tokenize(text)
    c = nltk.ConcordanceIndex(all_words, key = lambda s: s.lower())
    all_words_str = [i.lower() for i in all_words]
    lines = nltk.line_tokenize(text)
    sample_title = getDescriptionfromLine(lines, '!Sample_title')
    source_name = getDescriptionfromLine(lines, '!Sample_source_name_ch1')

    sample_type = getDescriptionfromLine(lines, '!Sample_molecule_ch1')
    words_in_type = sample_type.split()
    if 'RNA'in words_in_type and gseID != 'GSE54504': ### wrong description in this series, used genomic DNA.
        with open('excluded.log','a') as wf:
            wf.write(serarr+'\t'+'RNA'+'\n')
        continue

    puncs = ['!', '#', '(', ')', ',', '.', ':', ';', '=', '@', '<', '>']

    ## get pubmedID
    try:
        with open(os.path.join('/Volumes/arraymapIncoming/GEOmeta/',gseID,'geometa.soft'),'r') as f:
            content = f.readlines()
        for l in content:
            if l.startswith('!Series_pubmed_id'):
                pubmedID = l.partition(' = ')[2].rstrip()
                break
        else:
            pubmedID = ''
    except:
        pass

    ## look for normal WORDS, if normal save in normal path

    ## tumor/normal
    normalwords = ("reference", "normal", "healthy", "unaffected","germline","control","peripheral","remission","non-tumoral","population","populations")
    linenames = ['!Series_title','!Sample_title','!Sample_source_name_ch1','!Sample_characteristics_ch1']
    descriptionlines = []

    tumorORnormal=None
    for l in lines:
        if not tumorORnormal:
            for linename in linenames:
                if l.startswith(linename):
                    descriptionlines.append(l.rstrip())
                    line_words = [i.lower() for i in l.rstrip().split()]
                    if any([w for w in normalwords if w in line_words]):
                        tumorORnormal = 'normal'
                        outputpath = 'normals.tsv'
                        diagnosis = 'normal tissue'
                        diagnosis += ',' + getDescriptionfromLine(lines,'!Sample_characteristics_ch1')
                        stage = ''
                        # print('im normal')
                        # print(outputpath)
                        break
    if not tumorORnormal:
        tumorORnormal = 'tumor'

    description_words = set([w.lower() for w in nltk.word_tokenize(' '.join(descriptionlines))])
    ## if think it's cancer, look for 'oma' and find offsets to determine type
    ## cancertype
    if tumorORnormal == 'tumor':
        outputpath = None
        omas = ['glioblastoma','lymphoma','myeloma','neuroblastoma','medulloblastoma','melanoma','meningioma','ependymoma','retinoblastoma']
        typewords = ['colon','bladder','breast','cervix','kidney','leukemia','liver','lung','ovary','pancreas','prostate','rectal','thyroid','oral','esophagus','gastric']
    ## if determine type set save path
        omatype=set()
        for i in all_words_str:
            if i.endswith('oma') or i in ['neoplasm','myelodysplastic','lynch']:
                omatype.add(i)

        matchtype = [i for i in omatype if i in omas ]
        diagnosis = ''
        if len(matchtype) > 0:
            diagnosis = matchtype[0]
            outputpath = '{0}.tsv'.format(matchtype[0])
            for i in matchtype:
                diagnosis += ','+getBeforeAfter(i)
        else:
            for cancertype in typewords:
                ## liver may have prefix as hepatic, heptocellular ...
                if cancertype == 'liver':
                    cancertype_pref = ['liver','hepat']
                if cancertype == 'ovary':
                    cancertype_pref = ['ovary','ovarian']
                if cancertype == 'colon':
                    cancertype_pref = ['colon','colorect']
                if cancertype == 'rectal':
                    cancertype_pref = ['recto','rectal']
                if cancertype == 'kidney':
                    cancertype_pref = ['renal','kidney']
                ###for the other organs, take the first 4 letters
                else:
                    cancertype_pref = [cancertype[:5]]
                words = getWordFromPrefix(cancertype_pref,all_words_str)
                if words:
                    diagnosis = cancertype
                    outputpath = '{0}.tsv'.format(cancertype)

                    for i in words:
                        diagnosis += ','+getBeforeAfter(i)

                    break

    ## else save in unknown tumor path
        if not outputpath:
            if len(omatype) > 0:
                outputpath = '{}.tsv'.format(omatype.pop())
            else:
                outputpath = 'unknown.tsv'

    ## tnm stage

        stages = [all_words[offset+1:offset+3] for offset in c.offsets('stage')]
        stages = sum(stages,[])
        for found in range(len(stages)):
                stage = parseStage(stages[found])
                if stage:
                    break
        else:
            stage=''

    ## sample source (primary/metastasis)
        if any([i for i in ['metastasis','metastatic'] if i in description_words]):
            source = 'metastasis'
        elif 'relapse' in description_words:
            source = 'relapse'
        elif 'xenograft' in description_words:
            source = 'xenograft'
        elif 'line' in description_words:
            source = 'cell line'
        else:
            source = 'primary'
    ## here is for normal samples
    else:
        stage=''
        source = 'non-cancer'


    ## look for platform, city, country, year, gender, age, tnm stage (if tumor)

    ## infochannels
    infochannel_1 = getDescriptionfromLine(lines, '!Sample_characteristics_ch1')
    infochannel_2 = getDescriptionfromLine(lines, '!Sample_characteristics_ch2')

    charchannel_1 = getDescriptionfromLine(lines, '!Sample_source_name_ch1')
    charchannel_2 = getDescriptionfromLine(lines, '!Sample_source_name_ch2')
    ## organism
    organism = getDescriptionfromLine(lines, '!Sample_organism_ch1')
    if organism not in ['','Homo sapiens','homo sapiens']:
        # print(organism)
        with open('excluded.log','a') as wf:
            wf.write(serarr+'\t'+'non human'+'\n')
        continue

    ## platform
    platID = [all_words[offset+2] for offset in c.offsets('Sample_platform_id')] [0]

    ## city
    city = [all_words[offset+2:offset+5] for offset in c.offsets('Sample_contact_city')]
    try:
        while isinstance(city[0],list) :
            city=city[0]
            # print(city)
        for i in range(3):
            ## ! marks the start of the next line
            if city[i] == '!':
                city = ' '.join(city[:i])
                break
        else:
            city = ' '.join(city)
    except:
        city=''

    ## country
    country = [all_words[offset+2:offset+5] for offset in c.offsets('Sample_contact_country')]
    try:
        while isinstance(country[0],list):
            country = country[0]
        # print(country)

        for i in range(3):
            ## ! marks the start of the next line
            if country[i] == '!':
                country = ' '.join(country[:i])
                break
        else:
            country = ' '.join(country)

    except:
        country=''

    ## year
    year = [all_words[offset+4] for offset in c.offsets('Sample_last_update_date')][0]

    ## age
    agewords = ['age','years','year']

    for keyword in agewords:
        age = findAge(keyword)
        if age:
            break
    if not age:
        age = ''

    ## gender

    gender = findGender(['gender','sex'])
    if not gender:
        if len([i for i in ['M','Male','male'] if i in description_words]):
            gender = 'male'
        elif len([i for i in ['F','Female','female'] if i in description_words]):
            gender = 'female'
        else:
            gender = ''


    if not os.path.isfile(outputpath):
        with open(outputpath,'w') as wf:
            wf.write('\t'.join(['STATUS','SERIESID','UID','PLATFORMID','DIAGNOSISTEXT','DIAGNOSISTEXT_original',
                'ICDMORPHOLOGY','ICDMORPHOLOGYCODE','ICDTOPOGRAPHY','ICDTOPOGRAPHYCODE','TECHNOLOGY','SERIES_TITLE', 
                'BIOPROJECT','SAMPLE_TITLE','TNM','CHANNEL_CANCER','INFORCHANNEL_1','INFORCHANNEL_2','CHAR_CHANNEL_1',
                'CHAR_CHANNEL_2','SAMPLESOURCE','CELLLINEID','CELLOSAURUSID','AGE','SEX','CITY','COUNTRY','YEAR',
                'PMID','ORGANISM'])+'\n')
    with open(outputpath,'a') as wf:
        wf.write('\t'.join(['',gseID,gsmID,platID,'',diagnosis,'','','','','in situ oligonucleotide',series_title,
                            bioproject,sample_title,stage,'',infochannel_1,infochannel_2,charchannel_1,
                            charchannel_2,source,'','',age,gender,city,country,year,pubmedID,organism])+'\n')
    with open('done.log','a') as wf:
        wf.write(serarr+'\n')

    if count % 500 == 0:
        print(count)
