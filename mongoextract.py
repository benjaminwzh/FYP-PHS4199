from pymongo import MongoClient
from config import get_config
from datetime import datetime, timedelta
import json
import random
import pandas as pd

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


## Query data from collection
## Based on articles added to Recognise at a specific date
date_to_query = datetime.now().strftime("%d/%m/%Y")
# With timedelta: 
# (datetime.now() - timedelta(days=(2))).strftime("%d/%m/%Y")

query = {'date': date_to_query}
result = list(coll.find(query))
print(f"queried {len(result)} articles")


# Counters for each recognised entity
none_found = 0
reach_only = 0
biobert_only = 0
consensus = 0

# Collection of collated articles
articles = {'no_entities': {},
            'consensus': {},
            'rest': {}}

# Counting each entity's recognised state, and allocating its article to the collection
for article in result:
    id = article['_id']
    recognised = article['recognised']
    if not recognised:
        none_found += 1
        articles['no_entities'][id] = article
    else:
        for service in recognised.values():
            if service == 'consensus':
                consensus += 1
                if not (id in articles['consensus']):
                    articles['consensus'][id] = article
            elif service == 'Reach':
                reach_only += 1
                if not (id in articles['rest']):
                    articles['rest'][id] = article
            elif service == 'Biobert':
                biobert_only += 1
                if not (id in articles['rest']):
                    articles['rest'][id] = article

# Random sample 50% of articles with consensus entities
n_none = len(articles['consensus'].keys())
samples = round(n_none * 0.5)
screen = random.sample(list(articles['consensus'].keys()), samples)

for key in screen:
    article = articles['consensus'][key]
    id = article['_id']
    title = article['article_data']['title']
    abstract = article['article_data']['abstract']
    recognised = article['recognised']
    recognised = [f"{k}, {v}" for k, v in recognised.items()]
    article_data = [id, title, abstract]
    article_data.extend(recognised)
    with open("consensus_screen.txt", 'a', encoding='utf8') as f:
        f.write('\n'.join(article_data))
        f.write('\n')
        f.write('\n')
print(f'set up {samples} articles for screening')

# Random sample 5% of articles with no entities found
n_none = len(articles['no_entities'].keys())
samples = round(n_none * 0.05)
screen = random.sample(list(articles['no_entities'].keys()), samples)

for key in screen:
    article = articles['no_entities'][key]
    id = article['_id']
    title = article['article_data']['title']
    abstract = article['article_data']['abstract']
    recognised = article['recognised']
    recognised = [f"{k}, {v}" for k, v in recognised.items()]
    article_data = [id, title, abstract]
    article_data.extend(recognised)
    with open("none_screen.txt", 'a', encoding='utf8') as f:
        f.write('\n'.join(article_data))
        f.write('\n')
        f.write('\n')
print(f'set up {samples} articles with no entities for screening')

# Random sample 5% of articles with disagreements
remaining = len(articles['rest'].keys())
samples = round(remaining * 0.05)
screen = random.sample(list(articles['rest'].keys()), samples)

for key in screen:
    article = articles['rest'][key]
    id = article['_id']
    title = article['article_data']['title']
    abstract = article['article_data']['abstract']
    recognised = article['recognised']
    recognised = [f"{k}, {v}" for k, v in recognised.items()]
    article_data = [id, title, abstract]
    article_data.extend(recognised)
    with open("rest_screen.txt", 'a', encoding='utf8') as f:
        f.write('\n'.join(article_data))
        f.write('\n')
        f.write('\n')
print(f'set up {samples} articles with disagreeing entities for screening')

# Print confusion matrix between Reach and BioBERT
table = [[consensus, biobert_only],
         [reach_only, none_found]]
df = pd.DataFrame(table,
                  columns=['Reach included', 'Reach excluded'],
                  index=['Biobert included', 'Biobert excluded'])
print(df)