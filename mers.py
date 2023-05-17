from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight

# 메르스 바이러스 시퀀스 데이터 다운로드
from Bio import Entrez

# Entrez 이메일 등록
Entrez.email = "jaemin@o.cnu.ac.kr"

# NCBI에서 메르스 바이러스 시퀀스 검색
term = "MERS coronavirus"
handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
record = Entrez.read(handle)
handle.close()

# 검색된 시퀀스의 accession 번호 가져오기
accession = record["IdList"][0]

# 메르스 바이러스 시퀀스 데이터 가져오기
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
sequence = SeqIO.read(handle, "fasta")
handle.close()

# 시퀀스 특성 분석
sequence_length = len(sequence)
gc_content = GC(sequence.seq)
mw = molecular_weight(sequence.seq)

# 결과 출력
print("Accession 번호: ", accession)
print("메르스 바이러스 시퀀스 길이: ", sequence_length)
print("메르스 바이러스 시퀀스 GC 함량: {:.2f}%".format(gc_content))
print("메르스 바이러스 시퀀스 분자량: {:.2f} kDa".format(mw / 1000))