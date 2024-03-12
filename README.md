# FYP-PHS4199
Scripts and raw data used in NUS PHS4199 Final Year Project. For reference in the final report.

### Scripts
###### mongoextract.py
Queries MongoDB for articles that passed through Recognise. Maintain a collection of articles based on the service that recognised its entities. Samples each collection for articles for manual screening, by exporting their id, title, abstract, and entities into a .txt file.

###### sample_pmids.py
Filters the Gene2PubMed table for articles that are annotated with 1 human gene. Randomly samples 1000 articles from that collection, and outputs their PubMed IDs into a .txt file.

###### g2p_per_entity.py
Set up the Gene2PubMed gene IDs and Biomart for conversion of Entrez Gene IDs to HGNC Gene IDs. Queries MongoDB for articles sampled from Gene2PubMed, and for each article, checks each entity if it is the correct Gene2PubMed annotation. Calculates precision, recall, and F score.

###### g2p_per_article.py
Separate experiment that does g2p_per_entity.py, but checks whether the article contains a correctly annotated entity. Assumes all additional entities recognised can be removed safely, to see if precision can be improved.

### Raw data
###### g2p_samples_210224
A file that contains the PMIDs of the 1000 randomly sampled Gene2PubMed articles.

###### consensus_screen.txt, none_screen.txt, rest_screen.txt
These files contain the randomly sampled article data for manual screening. Contains each article's id, title, abstract, and recognised entities (if any).

#### References
Used packages from [pybiomart](https://github.com/jrderuiter/pybiomart), [pymongo](https://github.com/mongodb/mongo-python-driver), and [pandas](https://zenodo.org/record/3630805#.XjI9CyOIYdU).
