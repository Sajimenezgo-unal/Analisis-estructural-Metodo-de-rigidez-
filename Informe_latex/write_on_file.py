

count = 0
doc = []
with open('main.tex', 'w+',  encoding='utf8') as f:
    for line in f:
        print(line)
        doc.append(line)
for j in doc:
    count += 1
    if j == '%---Matrices Locales\n':
        doc.insert()
