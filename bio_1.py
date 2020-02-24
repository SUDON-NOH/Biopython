
from Bio.Seq import Seq

tatabox_seq = Seq("tataaaggcAATATGCAGTAG")

print(tatabox_seq)
# tataaaggcAATATGCAGTAG
print(type(tatabox_seq))
# <class 'Bio.Seq.Seq'>


print(tatabox_seq.transcribe())
# DNA 서열을 RNA로 전사해준다.
# uauaaaggcAAUAUGCAGUAG

print(tatabox_seq.translate())
# DNA 또는 RNA 서열을 단백질 서열로 번역해준다.
# YKGNMQ*

print(tatabox_seq.complement())
# 객체가 가진 서열의 상보적 서열을 만들어 반환한다.
# atatttccgTTATACGTCATC

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

tatabox_seq = Seq("tataaaggcAATATGCAGTAG", IUPAC.unambiguous_dna)
print(tatabox_seq)
print(type(tatabox_seq))

print(help(IUPAC))

# Sequence 객체의 메서드 중 count 메서드는 객체 안에 들어 있는 문자를 세어줌
exon_seq = Seq("ATGCAGTAG")
count_a = exon_seq.count("A")
print(count_a)

# GC-contents(%)
exon_seq = Seq("ATGCAGTAG")
g_count = exon_seq.count("G")
c_count = exon_seq.count("C")
gc_contents = (g_count + c_count) / len(exon_seq) * 100
print(gc_contents) # 44.44%


# 객체 서열 대소문자 변환하기
tatabox_seq = Seq("tataaaggcAATATGCAGTAG")
print(tatabox_seq.lower())
print(tatabox_seq.upper())

# Sequence 객체 전사, 번역하기
dna = Seq("ATGCAGTAG")
mrna = dna.transcribe()
ptn = dna.translate()
print(mrna)     # AUGCAGUAG
print(ptn)      # MQ*   *기호는 단백질 번역 과정이 끝나는 종결 코돈을 의미

# 첫 번째 종결 코돈에서 번역 종료하기
mRNA = Seq("AUGAACUAAGUUUAGAAU")
ptn = mRNA.translate()
print(ptn) # MN*V*N

ptn = mRNA.translate(to_stop=True)
print(ptn)

# 종결 코돈 기준으로 서열 나누기
mRNA = Seq("AUGAACUAAGUUUAGAAU")
ptn = mRNA.translate()
print(ptn) # MN*V*N

for seq in ptn.split("*"):
    print(seq)

# DNA Sequence 상보적, 역상보적 서열 만들기
seq = Seq("TATAAAGGCAATATGCAGTAG")
comp_seq = seq.complement()
rev_comp_seq = seq.reverse_complement()
print(seq)      # TATAAAGGCAATATGCAGTAG
print(comp_seq) # ATATTTCCGTTATACGTCATC
print(rev_comp_seq) # CTACTGCATATTGCCTTTATA

# 코돈 테이블
from Bio.Data import CodonTable

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
print(codon_table)

codon_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
print(codon_table)


# Sequence 객체에서 ORF 찾기
tatabox_seq = Seq("tataaaggcAATATGCAGTAG")
start_idx = tatabox_seq.find("ATG")
end_idx = tatabox_seq.find("TAG", start_idx)
orf = tatabox_seq[start_idx:end_idx]
print(orf)       # ATGCAG
print(start_idx) # 12

# Bio.SeqUtils 모듈 활용
from Bio.SeqUtils import GC

# GC-contents 계산하기
exon_seq = Seq("ATGCAGTAG")
gc_contents = GC(exon_seq)
print(gc_contents)

# 서열의 무게 계산하기
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import molecular_weight

seq1 = Seq("ATGCAGTAG")
seq2 = Seq("ATGCAGTAG", IUPAC.unambiguous_dna)
seq3 = Seq("ATGCAGTAG", IUPAC.protein)

print(molecular_weight(seq1))
print(molecular_weight(seq2))
print(molecular_weight(seq3))

# 가능한 모든 번역 구하기
from Bio.Seq import Seq
from Bio.SeqUtils import six_frame_translations

seq1 = Seq("ATGCCTTGAAATGTATAG")
print(six_frame_translations(seq1))

# Tm 계산
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

myseq = Seq("AGTCTGGGACGGCGCGGCAATCGCA")
print(mt.Tm_Wallace(myseq))

# 아미노산 서열의 약자와 기호간 변환하기
from Bio.SeqUtils import seq1

essential_amino_acid3 = "LeuLysMetValIleThrTrpPhe"
print(seq1(essential_amino_acid3)) # LKMVITWF

from Bio.SeqUtils import seq3

essential_amino_acid3 = "LKMVITWF"
print(seq3(essential_amino_acid3)) # LeuLysMetValIleThrTrpPhe

