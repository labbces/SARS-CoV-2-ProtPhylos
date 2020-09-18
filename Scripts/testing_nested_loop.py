
teste = open('teste_DO_nested.txt', 'w')
count = 0
for number in (1, 5, 26, 59, 48, 69):
    print(number)
    for i in range(0, 100):
        print(i)
        if number == i:
            print(number, 'MATCH')
            teste.write(f'{i}\n')
            count += 1
            print(count)
            break
