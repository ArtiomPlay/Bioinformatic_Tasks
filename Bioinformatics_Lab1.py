import math
from textwrap import wrap
from Bio.Seq import Seq

StartCodon='ATG'
StopCodons=['TAA','TAG','TGA']
Directory="viruses/data/"
Files=[
    "bacterial1.fasta",
    "bacterial2.fasta",
    "bacterial3.fasta",
    "bacterial4.fasta",
    "mamalian1.fasta",
    "mamalian2.fasta",
    "mamalian3.fasta",
    "mamalian4.fasta",
]
AminoAcids="ACDEFGHIKLMNPQRSTVWY"

def getReadingFrame(sequence,frame):
    return wrap(sequence[frame:],width=3)

def cosineDistance(vecA,vecB):
    dot=sum(vecA[k]*vecB[k] for k in vecA.keys())
    magnA=math.sqrt(sum(v**2 for v in vecA.values()))
    magnB=math.sqrt(sum(v**2 for v in vecB.values()))
    return 1-(dot/(magnA*magnB))

allFrequencies=dict()

for file in Files:
    # Opens file and skips fasta header
    f=open(Directory+file,"r")
    f.readline()

    # Reads rest file content, removes all new lines and closes file
    dnaSequence=f.read().replace("\n","")
    f.close()

    # Makes reverse complement
    complement_table=str.maketrans("ATGC","TACG")
    complement=dnaSequence.translate(complement_table)
    reverseComplement=complement[::-1]

    # 3 Part of task
    filteredOrfRegions=[]

    # Makes all 6 reading frames
    for strand,seq in [("forward",dnaSequence),("reverse",reverseComplement)]:
        for frame in range(3):
            starts=[]
            stops=[]
            orfRegions=[]
            readingFrame=getReadingFrame(seq,frame)

            # Searches for start and stop codons in reading frame
            for codonIndex,codon in enumerate(readingFrame):
                if codon==StartCodon:
                    starts.append(codonIndex)
                if codon in StopCodons:
                    stops.append(codonIndex)

            # 1 Part of task
            for start in starts:
                for stop in stops:
                    if start<stop:
                        orfRegions.append((start,stop,"".join(readingFrame[start:stop+1])))
                        break

            # 2 Part of task
            lastStop=-1
            for stop in stops:
                for start in starts:
                    if lastStop<start<stop:
                        orfRegions.append((start,stop,"".join(readingFrame[start:stop+1])))
                        break
                lastStop=stop

            for orfRegion in orfRegions:
                if len(orfRegion[2])>=100:
                    filteredOrfRegions.append((strand,frame,orfRegion[2]))

    # 4 Part of task
    proteinSeq=[]
    for orfRegion in filteredOrfRegions:
        orf_seq=Seq(orfRegion[2])
        protein_seq=orf_seq.translate(to_stop=True)
        proteinSeq.append((orfRegion[0],orfRegion[1],str(protein_seq)))

    # 5 Part of task
    # Codons
    codons={a:0.0 for a in AminoAcids}
    codonCount=0
    for protein in proteinSeq:
        codonCount=codonCount+len(protein[2])
        for codon in protein[2]:
            codons[codon]+=1

    for codon in codons:
        codons[codon]=codons[codon]/codonCount

    # Dicodons
    dicodons={a+b:0.0 for a in AminoAcids for b in AminoAcids}
    dicodonCount=0
    for proteinRegion in proteinSeq:
        for i in range(0,len(proteinRegion[2])-1):
            dicodon=proteinRegion[2][i]+proteinRegion[2][i+1]
            dicodons[dicodon]+=1
            dicodonCount=dicodonCount+1

    for dicodon in dicodons:
        dicodons[dicodon]=dicodons[dicodon]/dicodonCount

    allFrequencies[file.split(".")[0]]={
        "codons":codons,
        "dicodons":dicodons
    }

for frequency in allFrequencies:
    print(allFrequencies[frequency])

with open("CodonDistanceMatrix.txt","w") as f:
    f.write(f"{len(allFrequencies)}")

    for file1 in allFrequencies:
        f.write(f"\n{file1}")

        for file2 in allFrequencies:
            dist=max(0.0,cosineDistance(allFrequencies[file1]["codons"],allFrequencies[file2]["codons"]))
            f.write(" %.3f" % round(dist,3))

    f.close()

with open("DicodonDistanceMatrix.txt","w") as f:
    f.write(f"{len(allFrequencies)}")

    for file1 in allFrequencies:
        f.write(f"\n{file1}")

        for file2 in allFrequencies:
            dist=max(0.0,cosineDistance(allFrequencies[file1]["dicodons"],allFrequencies[file2]["dicodons"]))
            f.write(" %.3f" % round(dist,3))

    f.close()