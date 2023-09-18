import csv
from parse import parse
csvfile = open("/Users/jfeng1/kinase_data_all.csv","rb")
csv_reader = csv.reader(csvfile)
i = 0
colnameDict = {}
colnames = []
concConversionDict = {
    "0.100":"%cntrl@100nM",
    "1.000":"%cntrl@1000nM"
}

bionumbers = [" "]

concentrationDict = {}
data_dict = {}
kinase_dict = {}
conc_dict = {}
for row in csv_reader:
    if i == 0:
        for col_id,name in enumerate(row):
            p = parse("DISC RX KINASE ASSAY;Mean;%Binding (%);Protocol: {kinase};Concentration: {concentration} (uM);(Num)",name)
            if p is not None:
                concentration = p["concentration"]
                kinase = p["kinase"]
                if concentration in ["0.100","1.000"]:
                    kinase_dict[col_id] = kinase
                    conc_dict[col_id] = concentration
                    if not concentrationDict.has_key(kinase):
                        concentrationDict[kinase] =[]
                    if concentration not in concentrationDict[kinase]:
                        converted_conc = concConversionDict[concentration]
                        colnameDict["%s_%s"%(kinase,converted_conc)] = col_id
                        concentrationDict[kinase].append(converted_conc)
                        sorted(concentrationDict[kinase])
                        if len(concentrationDict[kinase])==2:
                            if kinase not in colnames:
                                colnames.append(kinase)
    else:
        for col_id,data in enumerate(row):
            bionumber = row[0]
            if col_id == 0:
                bionumbers.append(bionumber)
                bionumbers.append(bionumber)
            else:
                if kinase_dict.has_key(col_id) and conc_dict.has_key(col_id):
                    concentration = conc_dict[col_id]
                    converted_conc = concConversionDict[concentration]
                    kinase_name = kinase_dict[col_id]
                    key = "%s_%s_%s"%(bionumber,kinase_name,converted_conc)
                    try:
                        data_dict[key] = 100-float(data)
                    except:
                        data_dict[key] = 100

    i += 1
sorted(colnames)
print len(data_dict)
conc_row = ["DiscoveRx Gene Symbol"]
for id,bionumber in enumerate(bionumbers[1:]):
    if id%2==0:
        conc_row.append(concConversionDict['0.100'])
    else:
        conc_row.append(concConversionDict['1.000'])
# for key in colnames:
#     key1 = "%s_%s"%(key,concentrationDict[key][0])
#     key2 = "%s_%s"%(key,concentrationDict[key][1])
#     print key,colnameDict[key1],colnameDict[key2]

csv_output = open("/Users/jfeng1/kinase_output.csv","w")
csv_writer = csv.writer(csv_output)
csv_writer.writerow(bionumbers)
csv_writer.writerow(conc_row)
for kinase in colnames:
    row = [kinase]
    for id,converted_conc in enumerate(conc_row[1:]):
        bionumber = bionumbers[id+1]
        key = "%s_%s_%s"%(bionumber,kinase,converted_conc)
        print key
        if data_dict.has_key(key):
            row.append(data_dict[key])
        else:
            row.append(100)
    csv_writer.writerow(row)
csv_output.close()