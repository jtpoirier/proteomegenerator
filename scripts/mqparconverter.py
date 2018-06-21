import os
raw_files = []
for file in os.listdir(snakemake.input[2]):
    if file.endswith(".raw"):
        raw_files.append(os.path.join(snakemake.input[2], file))

with open(snakemake.input[1]) as oldMQPar, open(snakemake.output[0],"w") as newMQPar:
    for line in oldMQPar:
        newMQPar.write(line)
        if '<FastaFileInfo>' in line:
            newMQPar.write("<fastaFilePath>" + snakemake.input[0] + "</fastaFilePath>\n")
        if '<pluginFolder></pluginFolder>' in line:
            newMQPar.write("<numThreads>"+snakemake.params.t+"</numThreads>\n")
        if '<filePaths>' in line:
            for k in range(len(raw_files)):
                newMQPar.write("<string>" + raw_files[k] + "</string>\n")
        if '<experiments>' in line:
            for k in range(len(raw_files)):
                newMQPar.write("<string></string>\n")
        if '<fractions>' in line:
            for k in range(len(raw_files)):
                newMQPar.write("<short>32767</short>\n")
        if '<ptms>' in line:
            for k in range(len(raw_files)):
                newMQPar.write("<boolean>False</boolean>\n")
        if '<paramGroupIndices>' in line:
            for k in range(len(raw_files)):
                newMQPar.write("<int>0</int>\n")