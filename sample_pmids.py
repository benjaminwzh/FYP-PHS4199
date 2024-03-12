import pandas as pd
import random

df = pd.read_csv('gene2pubmed', sep='\t')
df_human = df.loc[df['#tax_id'] == 9606]

# Filter by articles with 1 gene
pmid_one_gene = df_human.groupby('PubMed_ID').filter(lambda x: len(x) == 1)

pubmed_ids = list(pmid_one_gene['PubMed_ID'])
print(len(pubmed_ids)) # 517595 articles
sample = random.sample(pubmed_ids, 1000)
sample = list(map(str, sample))

# with open('g2p_samples_210224', 'w') as f:
#     f.write('\n'.join(sample))