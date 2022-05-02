import os
from Bio import SeqIO

def extract_single_fasta(fasta_address, singlefasta_folder):
    all_fasta = open(fasta_address)
    fasta_name = []
    file_name = []
    num = 0
    for record in SeqIO.parse(all_fasta, "fasta"):
        num += 1
        name = record.name
        description=record.description
        seq = str(record.seq)[:len(str(record.seq)) - 1]
        if '|' in description:
            description=description.replace('|','_')
        if '.' in description:
            description=description.replace('.','_')
        if ' ' in description:
            description=description.replace(' ','@')
        fasta_name.append(name)
        single_file_name = str(num) + ".fasta"
        file_name.append(single_file_name)
        path=os.path.join(singlefasta_folder,single_file_name)
        single_fasta = open(path, 'w')
        single_fasta.write('>' + description + '\n')
        single_fasta.write(seq)
        single_fasta.close()
    all_fasta.close()

    return fasta_name,file_name

def findSubStrIndex(substr, str, time):
    times = str.count(substr)
    if (times == 0) or (times < time):
        pass
    else:
        i = 0
        index = -1
        while i < time:
            index = str.find(substr, index+1)
            i+=1
        return int(index)