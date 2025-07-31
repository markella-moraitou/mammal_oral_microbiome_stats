import pandas as pd
import sys
import os
import time
from external.omnicrobe import Taxon, OBT, Relation, version

## Define functions
def search_OBT(search_term, obt_type):
    ''' Gets a search term and OBT type and returns a list of OBTs that match the search term '''
    regex_term = "^" + search_term + "$"
    search_results = pd.DataFrame(columns=["search_term", "name", "id", "type"])
    search = OBT.search(s=regex_term, obt_type=obt_type)
    for o in search:
        new_row = pd.DataFrame([{"search_term": search_term, "name": o.name, "id": o.identifier, "type": o.obt_type.name}])
        search_results = pd.concat([search_results, new_row], ignore_index=True)
    return search_results

def search_OBT_list(search_term_list, obt_type_list):
    ''' Runs search_OBT for a list of search terms and OBT types and returns a dataframe with the results '''
    # search term list and obt_type list must have the same length
    if len(search_term_list) != len(obt_type_list):
        print("search term list and obt_type list must have the same length")
        return
    search_results = pd.DataFrame(columns=["search_term", "name", "id", "type"])
    for i in range(len(search_term_list)):
        term = search_term_list[i]
        obt_type = obt_type_list[i]
        print("Searching for term '" + term + "' in '" + obt_type +"'")
        results = search_OBT(term, obt_type)
        search_results = pd.concat([search_results, results], ignore_index=True)
    return search_results

def search_relations(taxon, OBT_search):
    ''' Gets a taxon and a list of OBTs and finds if the taxon is related to any of them '''
    relations = pd.DataFrame(columns=["taxon", "OBT", "id", "obt_type"])
    # Check if taxon causes error, in which case it was not created properly and will be skipped
    try:
        hasattr(taxon, 'name')
    except Exception as e:
        print(f"Error accessing taxon")
        return relations  # Return an empty DataFrame or handle as needed
    print("Searching relations for taxon: " + taxon.name)
    total_rels = 0
    for i in range(0, OBT_search.shape[0]):
        obt_id = OBT_search.iloc[i]["id"]
        obt_type = OBT_search.iloc[i]["type"]
        # Get OBT object
        obt = OBT(obt_id)
        # Check that OBT is of the expected type
        if obt.obt_type.name != obt_type:
            print("OBT " + obt.name + " is not the correct type")
            continue
        max_retries = 2
        for attempt in range(max_retries):
            try:
                rels = list(Relation.search(taxon=taxon, obt=obt))
                break  # If successful, exit the retry loop
            except Exception as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                time.sleep(10)  # Wait 10 second before retrying
                if attempt == max_retries - 1:
                    print("Max retries reached. Skipping this OBT.")
                    rels = []
        total_rels += len(rels)
        for r in rels:
            new_row = pd.DataFrame([{"taxon": taxon.name, "OBT": r.obt.name, "obt_type": obt_type, "id": r.obt.identifier}])
            relations = pd.concat([relations, new_row], ignore_index=True)
    print("Total relations found: " + str(total_rels))
    return relations

def simplify_results(relations):
    ''' Keeps only unique entries and adds a column with the number of times the relation was found '''
    simplified = relations.groupby(list(relations.columns.values)).size().reset_index(name='occurences')
    return simplified

def search_relations_taxa_list(taxon_list, search_terms_list, obt_type_list, output):
    ''' Gets a list of taxa and a list of search terms and obt types and finds if the taxa are related to any of the search terms '''
    all_relations = []
    # Identify OBTs to search
    obts = search_OBT_list(search_terms_list, obt_type_list)
    # The check the OBTs agaist each taxon
    for t in taxon_list:
        try:
            taxon_obj = Taxon(t)
            relations = search_relations(taxon_obj, obts)
            # Simplify results
            relations = simplify_results(relations)
            all_relations.append(relations)
            relations.to_csv(output, mode = 'a',
                             header=not os.path.exists(output),
                             index = False)
        except Exception as e:
            print(f"Error processing taxon {t}: {e}")

if __name__ == "__main__":
    print("Accessing omnicrobe database")
    print('Version')
    print(version())
    # Define search terms
    habitat_terms = ["laboratory equipment", "marine water", "deep sea", "soil",
                 "mammalian", "wild animal", "mammalian livestock",
                 "biofilm in natural environment", "host associated biofilm",
                 "mouth", "dental plaque", "gut", "rumen"]
    
    phenotype_terms = ["biofilm forming", "non-biofilm forming", 
                   "animal symbiont", "animal commensal", "animal pathogen", "plant hosted",
                   "aerobe", "aerotolerant", "anaerobe"]
    
    use_terms = ["antimicrobial activity", "amination activity", "fermentative", "lipolytic activity", "proteolytic activity"]
    
    # Read CLI arguments
    print("Reading arguments:")
    taxon_list_path = sys.argv[1]
    #taxon_list_path = "../output/community_analysis/taxids.tmp"
    
    print("Taxon list:" + taxon_list_path)
    
    outdir = sys.argv[2]
    #outdir = "../output/community_analysis/"
    print("Output directory: " + outdir)
    
    relation_type = sys.argv[3]
    
    # Read taxon list
    taxon_list = pd.read_csv(taxon_list_path, header=None)
    # Paste "ncbi:" to each taxon id in the list
    taxon_list = ["ncbi:" + str(taxon) for taxon in taxon_list[0]]
    
    if relation_type == "habitat":
        print("HABITATS:")
        habitat_relations = search_relations_taxa_list(taxon_list, habitat_terms, ["habitat"]*len(habitat_terms), os.path.join(outdir, "habitat_relations.csv"))
    
    if relation_type == "phenotype":
        print("PHENOTYPES:")
        phenotype_relations = search_relations_taxa_list(taxon_list, phenotype_terms, ["phenotype"]*len(phenotype_terms), os.path.join(outdir, "phenotype_relations.csv"))
    
    if relation_type == "use":
        print("USES:")
        use_relations = search_relations_taxa_list(taxon_list, use_terms, ["use"]*len(use_terms), os.path.join(outdir, "use_relations.csv"))
