from textwrap import wrap

StartCodon='ATG'
StopCodons=['TAA','TAG','TGA']
Directory="viruses/data/"

def getReadingFrame(sequence,frame):
    return wrap(sequence[frame:],width=3)

# Read file name with type. Example bacterial1.fasta
print("Enter file name: ")
fileName=input()
#fileName="bacterial1.fasta"
# Opens file and skips fasta header
f=open(Directory+fileName,"r")
f.readline()

# Reads rest file content, removes all new lines and closes file
dnaSequence=f.read().replace("\n","")
f.close()

# Makes reverse complement
complement_table=str.maketrans("ATGC", "TACG")
complement=dnaSequence.translate(complement_table)
reverseComplement=complement[::-1]

# Makes all 6 reading frames
for strand, seq in [("forward", dnaSequence), ("reverse", reverseComplement)]:
    for frame in range(3):
        starts=[]
        stops=[]
        orfRegions1Part=[]
        orfRegions2Part=[]
        readingFrame=getReadingFrame(seq, frame)

        # Searches for start codons in reading frame
        for codonIndex,codon in enumerate(readingFrame):
            if codon==StartCodon:
                starts.append(codonIndex)
            if codon in StopCodons:
                stops.append(codonIndex)

        for start in starts:
            for stop in stops:
                if start<stop:
                    orfRegions1Part.append((start,stop,".".join(readingFrame[start:stop+1])))
                    break

        lastStop=-1
        for stop in stops:
            for start in starts:
                if lastStop<start<stop:
                    orfRegions2Part.append((start,stop,".".join(readingFrame[start:stop+1])))
                    break
            lastStop=stop

        print(f"{strand} frame {frame}: {readingFrame}")
        print(f"Start codons: {starts}")
        print(f"Stop codons: {stops}")
        print(f"ORF regions for 1 part: {orfRegions1Part}")
        print(f"ORF regions for 2 part: {orfRegions2Part}\n")