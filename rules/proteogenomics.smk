rule download_protein_xml:
    output: "human.protein.xml.gz"
    shell: "wget -O human.protein.xml.gz https://www.uniprot.org/uniprot/?query=organism:9606&format=rdf&compress=yes" # https://www.uniprot.org/help/api_queries
