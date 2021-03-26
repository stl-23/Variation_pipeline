
def modify(sample1,sample2='',species='human',omic_type='WGS',control='N',outdir):
    outconfig = ''
    if species == 'human':
        if omic_type == 'WGS':
            if control == 'Y':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                ploidy = 2
                breakPointThreshold = .8
                #coefficientOfVariation = 0.01
                window = 50000
                #step=10000
                chrFiles = index/split_fa/
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
                """.format(sample1=sample1,sample2=sample2,out=outdir)
            elif control == 'N':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                ploidy = 2
                breakPointThreshold = .8
                #coefficientOfVariation = 0.01
                window = 50000
                #step=10000
                chrFiles = index/split_fa/
                outputDir = {out}
                [sample]
                mateFile = {sample1}
                #mateCopyNumberFile = test/sample.cpn
                inputFormat = BAM
                mateOrientation = 0
                """.format(sample1=sample1,out=outdir)
        elif omic_type == 'WES':
            if control == 'Y':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                window = 0
                ploidy = 2, 3
                outputDir = {out}
                breakPointType=4
                chrFiles = index/split_fa/
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
                """.format(out=outdir,sample1=sample1,sample2=sample2)
            elif control == 'N':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                window = 0
                ploidy = 2, 3
                outputDir = {out}
                breakPointType=4
                chrFiles = index/split_fa/
                maxThreads=6
                breakPointThreshold=1.2
                noisyData=TRUE
                printNA=FALSE
                readCountThreshold=50
                [sample]
                mateFile = {sample1}
                inputFormat = BAM
                mateOrientation = 0
                """.format(out=outdir,sample1=sample1)

    elif species == 'non-human':
        if omic_type == 'WGS':
            if control == 'Y':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                ploidy = 2
                breakPointThreshold = .8
                #coefficientOfVariation = 0.01
                window = 50000
                #step=10000
                chrFiles = index/split_fa/
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
                """.format(sample1=sample1,sample2=sample2,out=outdir)
            elif control == 'N':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                ploidy = 2
                breakPointThreshold = .8
                #coefficientOfVariation = 0.01
                window = 50000
                #step=10000
                chrFiles = index/split_fa/
                outputDir = {out}
                [sample]
                mateFile = {sample1}
                #mateCopyNumberFile = test/sample.cpn
                inputFormat = BAM
                mateOrientation = 0
                """.format(sample1=sample1,out=outdir)
        elif omic_type == 'WES':
            if control == 'Y':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                window = 0
                ploidy = 2, 3
                outputDir = {out}
                breakPointType=4
                chrFiles = index/split_fa/
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
                """.format(out=outdir,sample1=sample1,sample2=sample2)
            elif control == 'N':
                outconfig = """[general]
                chrLenFile = index/ref.fa.fai
                window = 0
                ploidy = 2, 3
                outputDir = {out}
                breakPointType=4
                chrFiles = index/split_fa/
                maxThreads=6
                breakPointThreshold=1.2
                noisyData=TRUE
                printNA=FALSE
                readCountThreshold=50
                [sample]
                mateFile = {sample1}
                inputFormat = BAM
                mateOrientation = 0
                """.format(out=outdir,sample1=sample1)
    return outconfig

