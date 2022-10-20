data = {}
with open("train.csv","r") as fin:
    keys = fin.readline()
    for aline in fin:
        linelist = aline.strip().split(',')
        if linelist[1] != '':
            data[int(linelist[0])] = aline

with open("train_updates_20220929.csv","r") as fin:
    keys = fin.readline()
    for aline in fin:
        linelist = aline.strip().split(',')
        if linelist[1] != '':
            data[int(linelist[0])] = aline

with open("train_new.csv","w") as fout:
    fout.write(keys)
    for values in data.values():
        fout.write(values)
    
