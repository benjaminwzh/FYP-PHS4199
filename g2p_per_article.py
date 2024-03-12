from pymongo import MongoClient
from config import get_config
from datetime import datetime, timedelta
import pandas as pd
import pybiomart

# Authenticate into MongoDB
try:
    client = MongoClient('mongodb://localhost:27017/',
                         username=get_config("mongo_user"),
                         password=get_config("mongo_password"),
                         authSource=get_config("mongo_db"),
                         authMechanism='SCRAM-SHA-1')
    db_name = get_config("mongo_db")
    db = client[db_name]
    print('Authenticated into Mongo')
except Exception as e:
    print(f"could not authenticate into Mongo: {e}")
    client = None
    db = None

# Access collection
try:
    coll = db[get_config("mongo_article_collection")]
except Exception as e:
    print(f"could not use collection {get_config('mongo_collection')}: {e}")
    coll = None

# Pubmed id of random samples
with open('g2p_samples_210224', 'r') as f:
    ids = f.readlines()
    ids = list(map(int, ids))

# Human gene2pubmed entries
df = pd.read_csv('gene2pubmed', sep='\t')
df_human = df.loc[df['#tax_id'] == 9606]

# Set up Biomart, remove NaN HGNC ids
dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl',
                            host='http://www.ensembl.org')
gene_ids = dataset.query(attributes=['entrezgene_id', 'hgnc_id', 'hgnc_symbol'])
gene_ids = gene_ids.loc[gene_ids['HGNC ID'].notnull()]

# Maintain data structure for correct entities
tp_article = {'Reach': [],
       'Biobert': [],
       'consensus': []}
fp_article = {'Reach': [],
       'Biobert': [],
       'consensus': []}
fn_article = {'Reach': [],
       'Biobert': [],
       'consensus': []}
recognised_art_n = []
valid_art_n = []

### Calculate per article basis
# Assuming all extra entities can be removed, service did a correct recognition if
# ID annotated by G2P is present among the entities recognised by a service/consensus
for id in ids:
    # Get NCBI gene ID from g2p article; article is mapped to 1 gene only
    ncbi_geneid = df_human.loc[df_human['PubMed_ID'] == id, "GeneID"].item()
    # Map NCBI gene ID to HGNC id, may have multiple HGNC ids to one NCBI id
    hgnc_id = gene_ids.loc[gene_ids['NCBI gene (formerly Entrezgene) ID'] == ncbi_geneid, 'HGNC ID'].to_numpy()
    hgnc_id = list(map(lambda id: id.split(':')[1], hgnc_id))
    hgnc_id = hgnc_id[0] if hgnc_id else ''

    # Get recognised entities from GTT
    pmid = f'PMID:{str(id)}'
    date_to_query = datetime.now().strftime("%d/%m/%Y")
    in_mongo = coll.find_one({'_id': pmid})
    if in_mongo: valid_art_n.append(pmid)
    query = {'_id': pmid,
             'date': date_to_query}
    mongo_result = coll.find_one(query)

    if mongo_result:
        recognised_art_n.append(pmid)
        # Dictionary where key is "HGNCid_symbol" and value is Reach, Biobert, or consensus
        recognised = mongo_result.get('recognised')

        # Modify according to schema of collation in GTT and data structure of results of interest
        # Handle false negatives - where 1 service did not recognise any entities
        services = list(recognised.values())
        reach_count, biobert_count, consensus_count = services.count('Reach'), services.count('Biobert'), services.count('consensus')
        if consensus_count == 0:
            fn_article['consensus'].append(pmid)
        if reach_count == 0 and consensus_count == 0:
            fn_article['Reach'].append(pmid)
        if biobert_count == 0 and consensus_count == 0:
            fn_article['Biobert'].append(pmid)
        
        # Have a list of all entities recognised by services
        reach_ent = [k for k, v in recognised.items() if v == 'Reach']
        reach_ent = list(map(lambda k: k.split("_")[0], reach_ent))
        biobert_ent = [k for k, v in recognised.items() if v == 'Biobert']
        biobert_ent = list(map(lambda k: k.split("_")[0], biobert_ent))
        consens_ent = [k for k, v in recognised.items() if v == 'consensus']
        consens_ent = list(map(lambda k: k.split("_")[0], consens_ent))
        all_ent = [k for k in recognised.keys()]
        all_ent = list(map(lambda k: k.split("_")[0], all_ent))
        
        # Find if annotated G2P id is present among recognised entities
        if hgnc_id in all_ent:
            if hgnc_id in reach_ent:
                tp_article['Reach'].append(pmid)
                if biobert_ent: fp_article['Biobert'].append(pmid)
                if consens_ent: fp_article['consensus'].append(pmid)
            elif hgnc_id in biobert_ent:
                tp_article['Biobert'].append(pmid)
                if reach_ent: fp_article['Reach'].append(pmid)
                if consens_ent: fp_article['consensus'].append(pmid)
            elif hgnc_id in consens_ent:
                tp_article['consensus'].append(pmid)
                tp_article['Biobert'].append(pmid)
                tp_article['Reach'].append(pmid)
        else:
            if reach_ent: fp_article['Reach'].append(pmid)
            if biobert_ent: fp_article['Biobert'].append(pmid)
            if consens_ent: fp_article['consensus'].append(pmid)

valid = len(valid_art_n)
total = len(recognised_art_n)
print(f'Total valid articles in mongo: {valid}')
print(f'Total articles that passed through Recognise: {total}')

# Reach calculations
tp = len(tp_article['Reach'])
fp = len(fp_article['Reach'])
fn = len(fn_article['Reach'])
precision = tp / (tp + fp)
recall = tp / (tp + fn)
fscore = (2 * precision * recall) / (precision + recall)
print(f'Reach, perarticle: TP = {tp}, FP = {fp}, FN = {fn}, Precision = {precision}, Recall = {recall}, F = {fscore}')

# Biobert
tp = len(tp_article['Biobert'])
fp = len(fp_article['Biobert'])
fn = len(fn_article['Biobert'])
precision = tp / (tp + fp)
recall = tp / (tp + fn)
fscore = (2 * precision * recall) / (precision + recall)
print(f'Biobert, perarticle: TP = {tp}, FP = {fp}, FN = {fn}, Precision = {precision}, Recall = {recall}, F = {fscore}')

# consensus
tp = len(tp_article['consensus'])
fp = len(fp_article['consensus'])
fn = len(fn_article['consensus'])
precision = tp / (tp + fp)
recall = tp / (tp + fn)
fscore = (2 * precision * recall) / (precision + recall)
print(f'Consensus, perarticle: TP = {tp}, FP = {fp}, FN = {fn}, Precision = {precision}, Recall = {recall}, F = {fscore}')

'''
Results:
Total valid articles in mongo: 969
Total articles that passed through Recognise: 689
Reach, perarticle: TP = 555, FP = 113, FN = 20, Precision = 0.8308383233532934, Recall = 0.9652173913043478, F = 0.8930008045052292
Biobert, perarticle: TP = 10, FP = 322, FN = 353, Precision = 0.030120481927710843, Recall = 0.027548209366391185, F = 0.02877697841726619
Consensus, perarticle: TP = 10, FP = 17, FN = 662, Precision = 0.37037037037037035, Recall = 0.01488095238095238, F = 0.02861230329041488
'''