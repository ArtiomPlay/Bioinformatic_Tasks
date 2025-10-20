from textwrap import wrap

StartCodon='ATG'
StopCodons=['TAA','TAG','TGA']
Directory="viruses/data/"

def getReadingFrame(sequence,frame):
    return wrap(sequence[frame:],width=3)

# Read file name with type. Example bacterial1.fasta
print("Enter file name: ")
#fileName=input()
fileName="bacterial1.fasta"
# Opens file and skips fasta header
f=open(Directory+fileName,"r")
f.readline()

# Reads rest file content, removes all new lines and closes file
dnaSequence=f.read().replace("\n","")
f.close()

# Makes reverse complement
complement_table=str.maketrans("ATGC", "TACG")
compliment=dnaSequence.translate(complement_table)
reverseComplement=compliment[::-1]

# Makes all 6 reading frames
for strand, seq in [("forward", dnaSequence), ("reverse", reverseComplement)]:
    for frame in range(3):
        starts=[]
        stops=[]
        readingFrame=getReadingFrame(seq, frame)

        # Searches for start codons in reading frame
        for codonIndex,codon in enumerate(readingFrame):
            if codon==StartCodon:
                starts.append(codonIndex)
            if codon in StopCodons:
                stops.append(codonIndex)

        print(f"{strand} frame {frame}: {readingFrame}")
        print(starts)
        print(stops)
