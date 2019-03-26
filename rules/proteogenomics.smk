rule download_protein_xml:
    output: "human.protein.xml"
    script: "scripts/download_protein_xml.py" # this works, but it's slow. Might try out urllib2
