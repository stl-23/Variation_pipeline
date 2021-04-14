import os
def modify(sample1,sample2,ref,outdir,species='human',omic_type='WGS',control='N'):
    outconfig = ''
    chrLenFile = ref+'.fai'
    chrFiles = os.path.dirname(ref)+'/split_fa/'
    if species == 'human':
        if omic_type == 'WGS':
            if control == 'Y':
                outconfig = """[general]
chrLenFile = {chrLenFile}
ploidy = 2
breakPointThreshold = .8
#coefficientOfVariation = 0.01
window = 50000
#step=10000
chrFiles = {chrFiles}
outputDir = {out}
[sample]
mateFile = {sample1}
#mateCopyNumberFile = test/sample.cpn
inputFormat = BAM
mateOrientation = 0
[control]
mateFile = {sample2}
#mateCopyNumberFile = test/sample.cpn
inputFormat = BAM
mateOrientation = 0
""".format(sample1=sample1,sample2=sample2,out=outdir,chrLenFile=chrLenFile,chrFiles=chrFiles)
            elif control == 'N':
                outconfig = """[general]
chrLenFile = {chrLenFile}
ploidy = 2
breakPointThreshold = .8
#coefficientOfVariation = 0.01
window = 50000
#step=10000
chrFiles = {chrFiles}
outputDir = {out}
[sample]
mateFile = {sample1}
#mateCopyNumberFile = test/sample.cpn
inputFormat = BAM
mateOrientation = 0
""".format(sample1=sample1,out=outdir,chrLenFile=chrLenFile,chrFiles=chrFiles)
        elif omic_type == 'WES':
            if control == 'Y':
                outconfig = """[general]
chrLenFile = {chrLenFile}
window = 0
ploidy = 2, 3
outputDir = {out}
breakPointType=4
chrFiles = {chrFiles}
maxThreads=6
breakPointThreshold=1.2
noisyData=TRUE
printNA=FALSE
readCountThreshold=50
[sample]
mateFile = {sample1}
inputFormat = BAM
mateOrientation = 0
[control]
mateFile = {sample2}
inputFormat = BAM
mateOrientation = 0
""".format(out=outdir,sample1=sample1,sample2=sample2,chrLenFile=chrLenFile,chrFiles=chrFiles)
            elif control == 'N':
                outconfig = """[general]
chrLenFile = {chrLenFile}
window = 0
ploidy = 2, 3
outputDir = {out}
breakPointType=4
chrFiles = {chrFiles}
maxThreads=6
breakPointThreshold=1.2
noisyData=TRUE
printNA=FALSE
readCountThreshold=50
[sample]
mateFile = {sample1}
inputFormat = BAM
mateOrientation = 0
""".format(out=outdir,sample1=sample1,chrLenFile=chrLenFile,chrFiles=chrFiles)

    elif species == 'non-human':
        if omic_type == 'WGS':
            if control == 'Y':
                outconfig = """[general]
chrLenFile = {chrLenFile}
ploidy = 2
breakPointThreshold = .8
#coefficientOfVariation = 0.01
window = 50000
#step=10000
chrFiles = {chrFiles}
outputDir = {out}
##if you are working with something non-human, we may need to modify these parameters:
#minExpectedGC = 0.35
#maxExpectedGC = 0.55

#readCountThreshold=10

#numberOfProcesses = 4
#outputDir = test
#contaminationAdjustment = TRUE
#contamination = 0.4
#minMappabilityPerWindow = 0.95
                
[sample]
mateFile = {sample1}
#mateCopyNumberFile = test/sample.cpn
inputFormat = BAM
mateOrientation = 0
[control]
mateFile = {sample2}
#mateCopyNumberFile = test/sample.cpn
inputFormat = BAM
mateOrientation = 0
""".format(sample1=sample1,sample2=sample2,out=outdir,chrLenFile=chrLenFile,chrFiles=chrFiles)
            elif control == 'N':
                outconfig = """[general]
chrLenFile = {chrLenFile}
ploidy = 2
breakPointThreshold = .8
#coefficientOfVariation = 0.01
window = 50000
#step=10000
chrFiles = {chrFiles}
outputDir = {out}
[sample]
mateFile = {sample1}
#mateCopyNumberFile = test/sample.cpn
inputFormat = BAM
mateOrientation = 0
""".format(sample1=sample1,out=outdir,chrLenFile=chrLenFile,chrFiles=chrFiles)
        elif omic_type == 'WES':
            if control == 'Y':
                outconfig = """[general]
chrLenFile = {chrLenFile}
window = 0
ploidy = 2, 3
outputDir = {out}
breakPointType=4
chrFiles = {chrFiles}
maxThreads=6
breakPointThreshold=1.2
noisyData=TRUE
printNA=FALSE
readCountThreshold=50
[sample]
mateFile = {sample1}
inputFormat = BAM
mateOrientation = 0
[control]
mateFile = {sample2}
inputFormat = BAM
mateOrientation = 0
""".format(out=outdir,sample1=sample1,sample2=sample2,chrLenFile=chrLenFile,chrFiles=chrFiles)
            elif control == 'N':
                outconfig = """[general]
chrLenFile = {chrLenFile}
window = 0
ploidy = 2, 3
outputDir = {out}
breakPointType=4
chrFiles = {chrFiles}
maxThreads=6
breakPointThreshold=1.2
noisyData=TRUE
printNA=FALSE
readCountThreshold=50
[sample]
mateFile = {sample1}
inputFormat = BAM
mateOrientation = 0
""".format(out=outdir,sample1=sample1,chrLenFile=chrLenFile,chrFiles=chrFiles)
    return outconfig

