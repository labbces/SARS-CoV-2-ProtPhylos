import re
name = input('Enter file name:')


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
    identificador = y[:-1]
    ids.append(identificador)


tabela = input('Nome do arquivo possuindo a Tabela de metadados do Gisaid:')

try:
    handle = open(tabela, 'r')
except:
    print('File cannot be open')
    exit()

CompCoverageTable = open('nova_tabela2.tsv', 'w')
for line in handle:
    title = line
    break
CompCoverageTable.write(title)


next_table = {}

for line in handle:
        string = line
        line = line.split()
        next_table[line[2]] = string # código para tabela nextmeta
#       next_table[line[1]] = string # código para a tabela do browser


matches = 0
chaves_out = list()
for epi in ids:
    if epi in next_table:
        string = next_table[epi]
        matches += 1
        CompCoverageTable.write(string)
    else:
        chaves_out.append(epi)
print(f'Número de matches: {matches}, IDs Totais: {len(ids)}')
print(f'Chaves que não batem:(opcional) ; número de chaves fora chaves_out: {len(chaves_out)}')


