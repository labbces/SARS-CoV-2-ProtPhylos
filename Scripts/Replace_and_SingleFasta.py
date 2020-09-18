# Lendo o arquivo
file_name = input('Name of the file:')
try:
    with open(file_name, 'r') as file:
        filedata = file.read()
except:
    print("File cannot be open")
    exit()

# Substituindo o "/" e o "|"
filedata = filedata.replace('/', 'BARRA')
filedata = filedata.replace('|', 'PIPE')

# Escrevendo o arquivo
with open(file_name, 'w') as file:
    file.write(filedata)
    file.close()

# Separando o arquivo multifasta em singlefasta
f = open(file_name, "r")
opened = False
for line in f:
    if line[0] == ">":
        if opened:
            of.close()
        opened = True
        of = open("%s.fa" % (line[1:].rstrip()), "w")  # Salvando como o identificador fasta
        print(line[1:].rstrip())
    of.write(line)  # Reescrevendo a sequência no identificador correspondente
of.close()
