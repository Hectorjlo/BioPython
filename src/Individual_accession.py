# Entrez verification
try:
    from Bio import Entrez
    ENTREZ_AVAILABLE = True
except ImportError:
    print("Warning: Biopython(Entrez) not installed.")
    ENTREZ_AVAILABLE = False



def entrez_credentials(entrez_email: str, entrez_api_key: str|None = None) -> tuple[str, str]:
    # Entrez credentials
    Entrez.email = entrez_email
    if entrez_api_key:
        Entrez.api_key = entrez_api_key
    
    return Entrez.email, Entrez.api_key

def get_project_IUD(project_acession: str) -> "Entrez.Parser.StringElement":
    
    # Handle
    handle = Entrez.esearch(db='bioproject', term=project_acession)
    # Retrieve
    entrez_search_results = Entrez.read(handle)
    # Close
    handle.close()

    # Get the UID
    project_uid = entrez_search_results['IdList'][0]
    print(f'[INFO] Bioproject IUD: {project_uid}')
    
    # print(type(project_uid))
    return project_uid
    
def summary_from_uid(project_iud):
    # Get a summary of the BioProject given the uid
    
    # Handle 
    handle = Entrez.esummary(db='bioproject', id=project_iud)
    # Read
    project_record = Entrez.read(handle)
    # Close
    handle.close()
    
    # Summary
    summary = project_record['DocumentSummarySet']['DocumentSummary'][0]

    # Show all keys of the record
    # print(f'[INFO] Record keys: {summary.keys()}')
    print('\n__________________________________\n')
    print(f'[INFO] Bioproject Info')
    print(f"[INFO] Access: {summary['Project_Acc']}")
    print(f"[INFO] Tittle: {summary['Project_Title']}")
    print(f"[INFO] Organism: {summary['Organism_Name']}")
    print(f"[INFO] Description: {summary.get('Project_Description')}")

def links_from_bioproject_to_GEO(project_iud):
    # Get db linked with the bioproject provieded 
    # Handle
    handle = Entrez.elink(dbfrom='bioproject', db='gds', id=project_iud)
    # Read
    links_of_bioproject = Entrez.read(handle)
    # Close
    handle.close()

    if links_from_bioproject_to_GEO[0]['LinkSetDb']:
        geo_ids = [link['Id'] for link in links_from_bioproject_to_GEO[0]['LinkSetDb'][0]['Link']]
        print(f'[INFO] The project with UID: {project_iud} is associated with the GEO data base')
        print(f'[INFO] UIDs with GEO associated: {geo_ids[:5]} ... total of: {len(geo_ids)}')
    else:
        print(f'[INFO] The project with UID: {project_iud} is NOT associated with the GEO')





if __name__ == '__main__':
    Entrez.email, Entrez.api_key = entrez_credentials(entrez_email='hectorjl@lcg.unam.mx', entrez_api_key='6416622762c8827ebf4f980b047bb0f25b08')
    project_iud = get_project_IUD(project_acession='PRJNA552284')
    project_summary = summary_from_uid(project_iud=project_iud)
    links_from_bioproject_to_GEO = links_from_bioproject_to_GEO(project_iud=project_iud)
    
    