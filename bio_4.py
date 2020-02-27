# Multiple Sequence Alignment

from Bio import AlignIO

alignment = AlignIO.read("C:/"
                         "Users/"
                         "SD NOH/"
                         "PycharmProjects/"
                         "First/"
                         "Bioinformatics_Biopython-master/"
                         "Bioinformatics_Biopython-master/"
                         "Section1/"
                         "Chap7/"
                         "example.aln", "clustal")

for record in alignment:
    print("%s - %s" % (record.seq, record.id))

for record in alignment:
    print("%s - %s" % (record.seq[:10], record.id))

# WebLogo로 보존 서열 알아보기

from Bio.motifs import Motif
from Bio import motifs
from Bio.Seq import Seq

instances = [Seq("TACAA"),
             Seq("TACGC"),
             Seq("TACAC"),
             Seq("TACCC"),
             Seq("AACCC"),
             Seq("AATGC"),
             Seq("AATGC")]
m = motifs.create(instances)

print(m.counts)
Motif.weblogo(m, 'test.png')

from Bio.Alphabet import IUPAC
# multiple sequence alignment 파일 읽기
alignment = AlignIO.read("C:/"
                         "Users/"
                         "SD NOH/"
                         "PycharmProjects/"
                         "First/"
                         "Bioinformatics_Biopython-master/"
                         "Bioinformatics_Biopython-master/"
                         "Section1/"
                         "Chap7/"
                         "HBA.aln", "clustal")

instance = []

# 서열 부분만 꺼내서 리스트에 넣기
for record in alignment:
    s = Seq(str(record.seq), IUPAC.protein)
    instance.append(s)

# Motif 객체 만들기기
m = motifs.create(instance)
Motif.weblogo(m, 'HBA_WebLogo.png')


# 계통수 그려보기
