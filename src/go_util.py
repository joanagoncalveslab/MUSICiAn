import os
import pickle as pkl
import biomart
from src.config import get_experiment_artifacts


def get_ensembl_mappings():
    filename = "{}/../gene_symbol_2_ids.pkl".format(get_experiment_artifacts())
    if os.path.exists(filename):
        return pkl.load(open(filename, 'rb'))


    # Set up connection to server                                               
    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')         
    mart = server.datasets['mmusculus_gene_ensembl']                            
                                                                                
    # List the types of data we want                                            
    attributes = ['ensembl_transcript_id', 'mgi_symbol', 
                  'ensembl_gene_id', 'entrezgene_id']
                                                                                
    # Get the mapping between the attributes                                    
    response = mart.search({'attributes': attributes})                          
    data = response.raw.data.decode('ascii')                                    
                                                                                
    gene_symbol_to_ensemble = {}                                                  
    # Store the data in a dict                                                  
    for line in data.splitlines():                                              
        line = line.split('\t')                                                 
        # The entries are in the same order as in the `attributes` variable
        transcript_id = line[0]                                                 
        gene_symbol = line[1]                                                   
        ensembl_gene = line[2]                                                  
        entrezgene_id = line[3] 
                                                                                
        # Some of these keys may be an empty string. If you want, you can 
        # avoid having a '' key in your dict by ensuring the 
        # transcript/gene/peptide ids have a nonzero length before
        # adding them to the dict
        gene_symbol_to_ensemble[gene_symbol] = {
            "transcript_id": transcript_id,
            "ensembl_gene": ensembl_gene,
            "entrezgene_id": entrezgene_id,
        }

    pkl.dump(gene_symbol_to_ensemble, open(filename, 'wb'))
                                                                                
    return gene_symbol_to_ensemble


if __name__ == "__main__":
    get_ensembl_mappings()