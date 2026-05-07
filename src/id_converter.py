from bioservices import UniProt
import time

def get_gene_names(uniprot_ids):
    """
    Queries UniProt using the bioservices library to retrieve gene names.
    
    Args:
        uniprot_ids (str or list): A single UniProt ID or a list of IDs.
        
    Returns:
        dict: A dictionary mapping UniProt IDs to Gene Names.
    """
    u = UniProt()
    
    if isinstance(uniprot_ids, str):
        ids_to_query = [uniprot_ids]
    else:
        ids_to_query = uniprot_ids
        
    mapping = {}
    print(f"Fetching gene names for {len(ids_to_query)} proteins using bioservices...")
    
    for count, uid in enumerate(ids_to_query, 1):
        try:
            # Query UniProt for the specific entry
            # We request 'gene_names' column specifically
            res = u.search(f"accession:{uid}", columns="id,gene_names")
            
            # bioservices.search returns a TSV string
            # Line 0: Header, Line 1: Data
            lines = res.strip().split('\n')
            if len(lines) > 1:
                data_parts = lines[1].split('\t')
                if len(data_parts) > 1:
                    # Get the primary gene name (usually the first one)
                    full_gene_names = data_parts[1].split()
                    gene_name = full_gene_names[0] if full_gene_names else "Unknown"
                    mapping[uid] = gene_name
                else:
                    mapping[uid] = "Unknown"
            else:
                mapping[uid] = "Unknown"
                
        except Exception as e:
            print(f"Error fetching ID {uid}: {e}")
            mapping[uid] = "Error"
            
        # Progress update
        if count % 5 == 0 or count == len(ids_to_query):
            print(f"Processed {count}/{len(ids_to_query)} IDs...")
            
        # Small delay to be polite to the UniProt servers
        time.sleep(0.1) 
        
    return mapping

if __name__ == "__main__":
    # Test with known IDs
    test_ids = ["P15056", "P08069", "O15111"]
    print("Running ID Converter Test (Bioservices)...\n")
    result_map = get_gene_names(test_ids)
    
    print("\n--- Test Results ---")
    for uid, gene in result_map.items():
        print(f"UniProt ID: {uid} --> Gene Name: {gene}")