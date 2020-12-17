from Bio import SeqIO

out = open("../../Final_test/gisadiFull.filtered.fasta", "w")

for record in SeqIO.parse("../../Final_test/gisadiFull.fasta", "fasta"):
    if len(record.seq) >= 29000:
        seqString = record.seq.upper()
        if seqString.count('N') <= len(record.seq) * 0.01:
            SeqIO.write(record, out, "fasta")
        else:
            print(f'{record.id}\t{len(record.seq)}\t{seqString.count("N")/len(record.seq)}\tN')
    else:
        print(f'{record.id}\t{len(record.seq)}\t{seqString.count("N") / len(record.seq)}\tL')



