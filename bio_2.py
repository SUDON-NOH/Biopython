# 5장

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord # Bio.SeqRecord에서 SeqRecord를 import 함

seq = Seq("ACGT")
seqRecord = SeqRecord(seq)
print(seqRecord)

# 속성값 넣어서 출력
simple_seq = Seq("ACGT")
simple_seqRecord = SeqRecord(simple_seq, id = "NC_1111", name = "TEST")
print(simple_seqRecord)

simple_seqRecord.name = "Another name"
print(simple_seqRecord)

simple_seqRecord.description = "테스트 설명입니다"
print(simple_seqRecord)

# SeqRecord 객체 만들기
seq = Seq("ACGT")
seqRecord = SeqRecord(seq)
print(seqRecord)
print("-------------------")

# SeqRecord 객체에 설명을 넣어준다.
seqRecord.id = "NC_1111"
seqRecord.name = "GeneA"
seqRecord.description = "This is a description"
seqRecord.annotations["Annotation_Key1"] = "Annotation_Value1"
seqRecord.annotations["Annotation_key2"] = "Annotation_Value2"
print(seqRecord)

# FASTA 파일로부터 SeqRecord 만들기
from Bio import SeqIO

record = SeqIO.read("C:/Users/SD NOH/PycharmProjects/First/Bioinformatics_Biopython-master/"
                    "Bioinformatics_Biopython-master/Section1/Chap5/J01636.1.fasta", "fasta")

print(type(record))
print(record)

print(type(record.seq))

# GenBank 파일로부터 SeqRecord 객체 만들기
record = SeqIO.read("C:/Users/SD NOH/PycharmProjects/First/Bioinformatics_Biopython-master/"
                    "Bioinformatics_Biopython-master/Section1/Chap5/J01636.1.gbk", "genbank")

print(type(record))
print(record)

# SeqRecord 객체간 비교하기

seq1 = Seq("ACGT")
record1 = SeqRecord(seq1)
print(record1)

seq2 = Seq("ACGT")
record2 = SeqRecord(seq2)
print(record2)

print(record1.seq == record2.seq)


