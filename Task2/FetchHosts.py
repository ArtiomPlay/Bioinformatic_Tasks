from Bio import Entrez
import time

Entrez.email="artiomgarbul@gmail.com"

with open("seq_names.txt") as f:
    accessions=[line.strip() for line in f if line.strip()]

output=[]

for acc in accessions:
    try:
        handle=Entrez.esearch(db="nucleotide", term=acc)
        record=Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            seq_id=record["IdList"][0]
            handle=Entrez.efetch(db="nucleotide",id=seq_id,rettype="gb",retmode="text")
            data=handle.read()
            handle.close()

            host="Unknown"
            for line in data.splitlines():
                if "/host=" in line:
                    host=line.split("=")[1].strip().strip('"')
                    break

            print(f"{acc}: {host}")
            output.append(f"{acc}\t{host}")
        else:
            print(f"{acc}: Not found")
            output.append(f"{acc}\tNot found")

        time.sleep(0.4)

    except Exception as e:
        print(f"{acc}: Error ({e})")
        output.append(f"{acc}\tError")

with open("accessions_hosts.tsv","w") as f:
    f.write("\n".join(output))