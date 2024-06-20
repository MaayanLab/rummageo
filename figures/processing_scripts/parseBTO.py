# Parse XML data from the Brenda Tissue Ontology to create a mapper of synomns
# for tissues and cell types, and cell lines (need better way to determine cell line entries) 
#http://purl.obolibrary.org/obo/bto.owl

import xml.etree.ElementTree as ET
import re
import json


def create_BTO_mappers(map_synomns=True):
    tree = ET.parse('data/bto.owl')
    root = tree.getroot()

    tissues_cell_types = {}
    cell_lines = {}
    synomn_mapper = {}

    patterns_to_remove = [r'\b(cell|line)\b', r'[-/]+']
    regex_patterns = re.compile('|'.join(patterns_to_remove), re.IGNORECASE)

    for child in root.iter('{http://www.w3.org/2002/07/owl#}Class'):
        try:
            bto = list(child.attrib.values())[0].split('/')[-1]
            name = child.find('{http://www.w3.org/2000/01/rdf-schema#}label').text.strip()
            synomns = []

            if map_synomns:
                for syn in child.iter('{http://www.geneontology.org/formats/oboInOwl#}hasRelatedSynonym'):
                    synomns.append(re.sub(r'[-_.]', ' ', syn.text.strip()))

            desc = child.find('{http://purl.obolibrary.org/obo/}IAO_0000115').text
            str_to_check = f"{desc} {name} {' '.join(synomns)}"

            # TODO: better segregation of cell lines
            if 'cell line' in str_to_check or 'cell lines' in str_to_check:
                stripped_name = re.sub(regex_patterns, '', name).strip()
                cell_lines[stripped_name] = bto
                if map_synomns:
                    for s in synomns:
                        synomn_mapper[s] = stripped_name
                        cell_lines[s] = bto
                    synomn_mapper[stripped_name] = stripped_name
            else:
                tissues_cell_types[re.sub(r'[-_.]', ' ', name).lower()] = bto
                if map_synomns:
                    for s in synomns:
                        synomn_mapper[s.lower()] = re.sub(r'[-_.]', ' ', name).lower()
                        tissues_cell_types[s.lower()] = bto
                    synomn_mapper[re.sub(r'[-_.]', ' ', name).lower()] = re.sub(r'[-_.]', ' ', name).lower()
        except:
            continue

    tissues_cell_types
    num_tissue_names = len(tissues_cell_types.keys())
    num_cell_line_names = len(cell_lines.keys())
    print('Extracted', num_tissue_names, 'tissue and cell type names and', num_cell_line_names, 'cell line names')

    with open('data/BTO_tissue_cell_types.json', 'w') as f:
        json.dump(tissues_cell_types, f)

    with open('data/BTO_cell_lines.json', 'w') as f:
        json.dump(cell_lines, f)



if __name__ == '__main__':
    create_BTO_mappers()