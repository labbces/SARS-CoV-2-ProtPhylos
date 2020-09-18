import re

# name = input('Enter file name:')
name = 'SeqBrasilJun.fasta'
try:
    hand = open(name, 'r')
except:
    print('File cannot be open')
    exit()

ids = list()

for line in hand:
    line = line.strip()
    if not line.startswith('>hCoV'):
        continue
    x = re.findall('EPI.+?P', line)
    y = x[0]
    id =y[:-1]
    ids.append(id)
print(ids)


# tabela = input('Nome do arquivo possuindo a Tabela de metadados do Gisaid:')
tabela = 'MetadadosBrasil.tsv'
try:
    handle = open(tabela, 'r')
except:
    print('File cannot be open')
    exit()

CompCoverageTable = open('nova_tabela.txt','w')
for line in handle:
    title = line
    break
CompCoverageTable.write(title)
count = 0


for epi in ids:
    print(epi)
    for line in handle:
        epis = line.split()
        if epi == epis[1]:
            print(epi, 'MATCH')
            count = count + 1
            print(count)
            CompCoverageTable.write(line)
            break

