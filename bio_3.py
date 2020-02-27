# 6장

# SeqIO.parse() 메서드로 파일 읽기
from Bio import SeqIO

# FASTA
seq = SeqIO.parse("C:/Users/SD NOH/PycharmProjects/First/Bioinformatics_Biopython-master"
                  "/Bioinformatics_Biopython-master/Section1/Chap6/sample_1.fasta", "fasta")
print(type(seq))
print(seq)

for s in seq:
    print(type(s))
    print(s)
    print("-"*150)


seq = SeqIO.parse("C:/"
                  "Users/"
                  "SD NOH/"
                  "PycharmProjects/"
                  "First/"
                  "Bioinformatics_Biopython-master/"
                  "Bioinformatics_Biopython-master/"
                  "Section1/"
                  "Chap6/"
                  "sample_2.fasta", "fasta")
print(type(seq))
for s in seq:
    print(type(s))
    print(s)
    print("-"*150)

# FASTQ

seq = SeqIO.parse("C:/"
                  "Users/"
                  "SD NOH/"
                  "PycharmProjects/"
                  "First/"
                  "Bioinformatics_Biopython-master/"
                  "Bioinformatics_Biopython-master/"
                  "Section1/"
                  "Chap6/"
                  "sample_1.fastq", "fastq")
for s in seq:
    print(s.seq)

# 압축된 FASTQ 파일 읽기

import gzip
handle = gzip.open("C:/"
                   "Users/"
                   "SD NOH/"
                   "PycharmProjects/"
                   "First/"
                   "Bioinformatics_Biopython-master/"
                   "Bioinformatics_Biopython-master/"
                   "Section1/"
                   "Chap6/"
                   "sample_1.fastq.gz", "rt")
seq = SeqIO.parse(handle, "fastq")
for s in seq:
    print(s.seq)

# GenBank 파일 읽기

gbk = SeqIO.read("C:/"
                 "Users/"
                 "SD NOH/"
                 "PycharmProjects/"
                 "First/"
                 "Bioinformatics_Biopython-master/"
                 "Bioinformatics_Biopython-master/"
                 "Section1/"
                 "Chap6/"
                 "kt225476.2.gbk", "genbank")

print(type(gbk))
print(gbk)

print(gbk.id)
print(gbk.description)
print(gbk.annotations['molecule_type'])
print(gbk.annotations['organism'])

# 인터넷을 통한 파일 읽기: CCR5 유전자와 HIV 저항성

from Bio import Entrez

Entrez.email = "lynx92@naver.com"
handle = Entrez.efetch(db="nucleotide", rettype = "gb", id = "AY463215", retmode="text")

for s in handle:
    print(s.strip())

# CCR5 유전자 정보 읽기

from Bio import Entrez
from Bio import SeqIO

Entrez.email = "lynx92@naver.com"

with Entrez.efetch(db="nucleotide", rettype="fasta", retmode = "text", id="42540826") as handle:
    seq = SeqIO.read(handle, "fasta")

print(seq)
print(len(seq))


