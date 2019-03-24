from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
fasta_sequences = SeqIO.parse(open("ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),'fasta')

ordered = []
chrn = range(1, 23)
chrs = []
x = ""
y = ""
m = ""
gl = []
ki = []
other = []
for seq in fasta_sequences:
    if seq.id.split(" ")[0] in chrn: chrs.append(seq)
    elif seq.id.split(" ")[0].startswith("X"): x = seq
    elif seq.id.split(" ")[0].startswith("Y"): y = seq
    elif seq.id.split(" ")[0].startswith("MT"): m = seq
    elif seq.id.split(" ")[0].startswith("GL"): gl.append(seq)
    elif seq.id.split(" ")[0].startswith("KI"): ki.append(seq)
    else: other.append(seq)
ordered.extend(chrs)
ordered.append(x)
ordered.append(y)
ordered.append(m)
ordered.extend(gl)
ordered.extend(ki)
ordered.extend(other)

with open("ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa","w") as out:
    for seq in ordered:
        out.write(seq.format("fasta"))
