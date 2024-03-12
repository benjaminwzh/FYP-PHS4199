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
true_positive = {'Reach': [],
       'Biobert': [],
       'consensus': []}
false_positive = {'Reach': [],
       'Biobert': [],
       'consensus': []}
false_negative = {'Reach': [],
       'Biobert': [],
       'consensus': []}
recognised_articles = []
valid_articles = []


### Calculate per entity basis
for id in ids:
    # Get NCBI gene ID from g2p article; article is mapped to 1 gene only
    ncbi_geneid = df_human.loc[df_human['PubMed_ID'] == id, "GeneID"].item()
    # Map NCBI gene ID to HGNC id, may have multiple HGNC ids to one NCBI id
    hgnc_id = gene_ids.loc[gene_ids['NCBI gene (formerly Entrezgene) ID'] == ncbi_geneid, 'HGNC ID'].to_numpy()
    hgnc_id = list(map(lambda id: id.split(':')[1], hgnc_id))

    # Get recognised entities from GTT
    pmid = f'PMID:{str(id)}'
    date_to_query = datetime.now().strftime("%d/%m/%Y")
    in_mongo = coll.find_one({'_id': pmid})
    if in_mongo: valid_articles.append(pmid)
    query = {'_id': pmid,
             'date': date_to_query}
    mongo_result = coll.find_one(query)

    if mongo_result:
        recognised_articles.append(pmid)
        # Dictionary where key is "HGNCid_symbol" and value is Reach, Biobert, or consensus
        recognised = mongo_result.get('recognised')

        # Modify according to schema of collation in GTT and data structure of results of interest
        # Handle false negatives - where 1 service did not recognise any entities
        services = list(recognised.values())
        reach_count, biobert_count, consensus_count = services.count('Reach'), services.count('Biobert'), services.count('consensus')
        if consensus_count == 0:
            false_negative['consensus'].append(pmid)
        if reach_count == 0 and consensus_count == 0:
            false_negative['Reach'].append(pmid)
        if biobert_count == 0 and consensus_count == 0:
            false_negative['Biobert'].append(pmid)
        
        # Handle true and false positives - for each entity, whether they == g2p annotation
        for k, v in recognised.items():
            recognised_hgnc = k.split("_")[0]
            if recognised_hgnc in hgnc_id:
                # True positive case
                if v == 'Reach':
                    true_positive['Reach'].append(pmid)
                elif v == 'Biobert':
                    true_positive['Biobert'].append(pmid)
                elif v == 'consensus':
                    true_positive['consensus'].append(pmid)
                    true_positive['Reach'].append(pmid)
                    true_positive['Biobert'].append(pmid)
            else:
                # False positive case
                if v == 'Reach':
                    false_positive['Reach'].append(pmid)
                elif v == 'Biobert':
                    false_positive['Biobert'].append(pmid)
                elif v == 'consensus':
                    false_positive['consensus'].append(pmid)
                    false_positive['Reach'].append(pmid)
                    false_positive['Biobert'].append(pmid)

## Precision and Recall -> per entity basis confusion matrix
# P = TP / (TP + FP)
# R = TP / (TP + FN)
# TP = entities that match the g2p mapping for a pubmed article
# FP = entities that are recognised but are not the correct g2p mapping
# FN = x if an article has no recognised entities, but the article has x g2p mappings
valid = len(valid_articles)
total = len(recognised_articles)
print(f'Total valid articles in mongo: {valid}')
print(f'Total articles that passed through Recognise: {total}')

# Reach calculations
tp = len(true_positive['Reach'])
fp = len(false_positive['Reach'])
fn = len(false_negative['Reach'])
precision = tp / (tp + fp)
recall = tp / (tp + fn)
fscore = (2 * precision * recall) / (precision + recall)
print(f'Reach: TP = {tp}, FP = {fp}, FN = {fn}, Precision = {precision}, Recall = {recall}, F = {fscore}')

# Biobert
tp = len(true_positive['Biobert'])
fp = len(false_positive['Biobert'])
fn = len(false_negative['Biobert'])
precision = tp / (tp + fp)
recall = tp / (tp + fn)
fscore = (2 * precision * recall) / (precision + recall)
print(f'Biobert: TP = {tp}, FP = {fp}, FN = {fn}, Precision = {precision}, Recall = {recall}, F = {fscore}')

# consensus
tp = len(true_positive['consensus'])
fp = len(false_positive['consensus'])
fn = len(false_negative['consensus'])
precision = tp / (tp + fp)
recall = tp / (tp + fn)
fscore = (2 * precision * recall) / (precision + recall)
print(f'Consensus: TP = {tp}, FP = {fp}, FN = {fn}, Precision = {precision}, Recall = {recall}, F = {fscore}')

'''
Results:
Total valid articles in mongo: 969
Total articles that passed through Recognise: 689
Reach: TP = 555, FP = 1563, FN = 20, Precision = 0.26203966005665724, Recall = 0.9652173913043478, F = 0.4121797252135166
Biobert: TP = 10, FP = 597, FN = 353, Precision = 0.016474464579901153, Recall = 0.027548209366391185, F = 0.020618556701030927
Consensus: TP = 10, FP = 17, FN = 662, Precision = 0.37037037037037035, Recall = 0.01488095238095238, F = 0.02861230329041488
'''