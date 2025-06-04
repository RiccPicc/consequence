import re
def read_headers(fasta):
    with open(fasta, "r") as handle:
        lines = handle.readlines()
        headers = [line[1:] for line in lines if line.startswith(">")]
    return headers

def get_species(header):
    species_match = re.search(r'(?<=\.\d\s)(.+?)(?=\s(?:gene|isolate|strain|ribosomal|spacer|clone|partial|16S|23S|internal))', header)
    if species_match:
        species_name = species_match.group(1).strip()
        species_name += f" {header.split()[0][1:]}"
        return species_name
    return None

def map_id_species(fasta):
    headers = read_headers(fasta) 
    ids_species_map = {header.split()[0]: get_species(header) for header in headers}
    return ids_species_map

def rename_leafs(treefile, ids_species_map):
    with open(treefile, "r") as handle:
        tree = handle.read()
        for id, species in ids_species_map.items():
            tree = tree.replace(id, species)
    return tree