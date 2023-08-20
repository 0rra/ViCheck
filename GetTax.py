import pandas as pd
from ete3 import NCBITaxa

# function to add taxonomy information to a dataframe given column of taxids
# data =  dataframe with column containing taxid, column must be titled either "Subject_taxid" or "Query_taxid"
# subquery = query or subject
# desrank = desired rank, must be lower case, e.g.: "species","family",class","genus","order", "kingdom" etc

def get_tax(data, subquery, desrank):
    # gets ncbi taxonomy dump
    ncbi = NCBITaxa()
    # converts columns of taxids to list of taxids
    if subquery == 'subject':
        taxids = data['Subject_taxid'].tolist()
    if subquery == 'query':
        taxids = data['Query_taxid'].tolist()

    # finds name of desired taxonomy rank    
    def get_rank_info(desrank):
        # ranks and names
        taxinfo = {}
        for taxid in taxids:
            if ";" in str(taxid):
                taxid = taxid.split(';')[0]
            lineage = ncbi.get_lineage(taxid)
            names = ncbi.get_taxid_translator(lineage)
            # get lineage ranks
            ranksO = ncbi.get_rank(lineage)
            # use ranks to get code...
            # re-order keys and values in rank dictionary
            ranks = {rank: id for id, rank in ranksO.items()}
            try:
                # using rank code, get name
                rankcode = ranks[desrank]
                name = names[rankcode]
                # add taxonomy information to dictionary
                taxinfo.update({taxid: name})
            # if no name was found for the rank, assign name "not present"
            except KeyError:
                taxinfo.update({taxid: 'Not present'})
        return taxinfo
    # maps rank name to taxid using the taxonomy info dicionary, adding to new column 
    if subquery == 'subject':
        data[f'Subject_{desrank}'] = data['Subject_taxid'].map(get_rank_info(desrank))
    if subquery == 'query':
        data[f'Query_{desrank}'] = data['Query_taxid'].map(get_rank_info(desrank))
    return data
