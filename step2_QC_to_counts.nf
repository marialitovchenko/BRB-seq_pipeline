#!/usr/bin/env nextflow

def helpMessage() {

    log.info """
    -      \033[41m A L I T H E A  G E N O M I C S   P I P E L I N E v1.0\033[0m-
    ================================================================================
    Welcome to the Nextflow BRB-seq analysis command line pipeline!

    Usage:
    The \033[1;91mtypical\033[0m command for running the pipeline is as follows:
    nextflow step2_QC_to_counts.nf \033[1;91m--user\033[0m UserName \\
                        \033[1;91m--pi\033[0m PiName \\    
                        \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--FQdir\033[0m fqDir \\
                        \033[1;91m--genomeDir\033[0m allGenomes \\
                        \033[1;91m--outputDir\033[0m theResult

    or

    nextflow forTest.nf \033[1;91m--user\033[0m UserName \\
                        \033[1;91m--pi\033[0m PiName \\    
                        \033[1;91m--inputTab\033[0m table.csv \\

    to put the pipeline into \033[1;91mbackground\033[0m mode:
    nextflow step2_QC_to_counts.nf \033[1;91m--user\033[0m UserName \\
                        \033[1;91m--pi\033[0m PiName \\    
                        \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--FQdir\033[0m fqDir \\
                        \033[1;91m--genomeDir\033[0m allGenomes \\
                        \033[1;91m--outputDir\033[0m theResult \\
                        \033[1;91m-bg\033[0m \\
                        \033[1;91m-N\033[0m your.email@gmail.com

    \033[1;91mMandatory\033[0m arguments:
      \033[1;91m--user\033[0m        User name, should not contain spaces
      \033[1;91m--pi\033[0m        PI name, should not contain spaces
      \033[1;91m--inputTab\033[0m        Path to the table containing information about input
                        data. The table should have following columns: RunID,
                        (i.e. NXT0540), LibraryID (i.e. nxid12916), SampleID
                        (i.e. BRBseq_v3_plate_1_S25), R1len (length of read 1,
                        i.e. 21), BU_ptrn (Barcode-UMI pattern, one of BU or 
                        UM), Specie (i.e. Hsapiens), Genome (i.e. hg38), 
                        Specie and Genome indicate to which genome version of which specie sample should be aligned to.

                        If no FQdir is provided (see below), the system will
                        assume that input fastq files are located in
                        [current dir]/RunID/LibraryID.

                        If no genomeDir is provided (see below), the system
                        will that STAR indexed genome is located in
                        [current dir]/Specie/Genome

                        If no outputDir is provided (see below), the system
                        will output files in the current directory

    
    \033[1;91mOptional\033[0m arguments:
    This arguments are not going to be needed with use of graphical user
    interface
      \033[1;91m--FQdir\033[0m           Path to the directory containing folders (one per run)
                        with fastq files
      \033[1;91m--genomeDir\033[0m       Path to the directory containing all your genome
                        versions for all your species. For example, a valid
                        genome directory TestGenomeDir would contain two
                        folders names mus_musculus and homo_sapiens.
                        Consequently, homo_sapiens folder would contain
                        GRCh37.75 and GRCh38.99, and mus_musculus would contain
                        GRCm38.68 and GRCm38.98. \033[93m Please use then homo_sapiens
                        or mus_musculus in a Specie column of your input table,
                        and use GRCh37.75/GRCh38.99/GRCm38.68/GRCm38.98 in a
                        Genome column.\033[0m
      \033[1;91m--outputDir\033[0m       Path to the output directory
      \033[1;91m--help\033[0m            Displays this message
      \033[1;91m-bg\033[0m               Puts execution of the pipeline into background mode
      \033[1;91m-N\033[0m                email adress in order to get notified upon pipeline complition.
                                                    Do not use epfl email address, because emails can't pass firewall. Use gmail.
      \033[1;91m-resume\033[0m           Resumes execution of the pipeline from the moment it
                        was interrupted
      """.stripIndent()
}

// Show help message
params.help = ''
if (params.help) {
    helpMessage()
    exit 0
}

/* ----------------------------------------------------------------------------
* Input handling
*----------------------------------------------------------------------------*/
// user and PI
user = params.user
pi = params.pi

// path to the input table with samples
sampleTabPath = file(params.inputTab)

//user directory: directory which contains all the fastq files
params.FQdir = file('.')
userDir = file(params.FQdir)
// genomes directory: directory with all compiled STAR indexed genomes
params.genomeDir = file('.')
genomePath = file(params.genomeDir)
// output folder
params.outputDir = file('.')
outputDir = file(params.outputDir) 
mapStatsTab = outputDir + "/mapStatsTab.csv"
// subfolder for files with used selected barcodes
usedBarcodeDir = outputDir + "/barcodeTables"

// technical directory, contains all support files, like scripts, jars, etc
params.techDir = 'techDir'
params.brbseqTools = file(params.techDir + '/BRBseqTools.1.5.jar')
params.barcodefile = file(params.techDir + '/barcodes_v3.txt')
params.compile_report = file(params.techDir + '/compile_report.R')
params.markdown = file(params.techDir + '/Generate_UserReport.Rmd')

/* ----------------------------------------------------------------------------
* LOG: inform user about all the inputs
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
sampleTabInfoCh = Channel.fromPath( sampleTabPath )
sampleTabInfoCh
    .splitCsv(header: true, sep:'\t')
    .into { logFqFiles; logGenomesFiles }

logFqFiles
    .map{ row ->  userDir.toString() + '/' + row.RunID.toString() + '/' + row.LibraryID.toString() }
    .unique()
    .toList()
    .set{logFqFiles}

logGenomesFiles
    .map{ row ->  genomePath.toString() + '/' + row.Specie.toString() + '/' + row.Genome.toString() }
    .unique()
    .toList()
    .set{logGenomesFiles}

log.info """\
         -        \033[1;91m   B R B - s e q   N E X T F L O W   P I P E L I N E \033[0m-   
         ================================================================================
         \033[1;91mInput summary: \033[0m
         Submitted input table          : ${sampleTabPath}
         Expect to find fastq-s in      : ${logFqFiles.toString().replaceAll(/DataflowVariable.value../, '').replaceAll(/..$/, '')}
         Expect to find genomes in      : ${logGenomesFiles.toString().replaceAll(/DataflowVariable.value../, '').replaceAll(/..$/, '')}
         BRBseq tools in                : ${params.brbseqTools}
         Output folder                  : ${outputDir}

         \033[1;91mExpected output summary:\033[0m
         Upon completion, following folders are going to be created:
         ${outputDir}/trimmed	:	folder containing trimmed fastqs
         ${outputDir}/demultiplexed	:	folder containing demultiplexed fastqs
         ${outputDir}/mapped	:	folder containing mapped bam files
         ${outputDir}/mapStats	:	folder containing log files produced by STAR
         ${outputDir}/counts	:	folder containing counts for individual samples
         \033[1;93m${outputDir}/countTables\033[0m	:	folder containing \033[1;93mfinal count tables\033[0m
         \033[1;93m${mapStatsTab}\033[0m  : a \033[1;93mfile with mapping statistics\033[0m and further visualized with \033[1;93mplotMapStats.R\033[0m


         \033[1;91mImportant note:\033[0m: you may not see some of the samples in the final
         count tables due to 0 reads being accosiated to the give barcodes.

         ================================================================================

         \033[1;91m L E T' S   G O ! ! ! \033[0m    
         """
         .stripIndent()

/* ----------------------------------------------------------------------------
* Read input table
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
sampleTabCh = Channel.fromPath( sampleTabPath )
sampleTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.RunID, row.LibraryID, row.SampleID, row.R1len, 
                       row.BU_ptrn, row.SampleName, row.pos, row.Specie,
                       row.Genome) }
    .into{ sampleTab; fqForQCtrim; barcodesPerRun}

// sampleTab will have the full information about submitted samples, and 
// fqForQCtrim will only have info determining unique fastq files. This is done
// in order to restrict number of QC and trimming processes only to the 
// nessecary ones, without repetition. barcodesPerRun will be further used to 
// only select barcodes which were used for sample multiplexing.

// groupTuple assures that there is no occasional slipage between the different
// lines coming from the table, nextflow has it time to time.
fqForQCtrim.flatMap { item ->
        RunID = item[0];
        LibraryID = item[1];
        SampleID = item[2];
        R1len = item[3];
        BU_ptrn = item[4];
        collect { onefile -> return [ RunID, LibraryID, SampleID, R1len,
                            BU_ptrn ] }
    }
    .groupTuple(by : [0, 1, 2, 3, 4])
    .unique()
    .set{ uniq_fqForQCtrim }

// read in avaible barcodes
barcodeTabCh = Channel.fromPath( params.barcodefile  )
barcodeTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.Name, row.B1) }
    .set{ barcodeTab}
// create a map telling which barcodes were used per sample and output it to a
// file. But header is missing! groupTuple assures that there is no occasional 
// slipage between the different lines coming from the table, nextflow has it 
// time to time.
barcodesPerRun
  .flatMap { item ->
        RunID = item[0];
        LibraryID = item[1];
        SampleID = item[2];
        Name = item[6];
        SampleName = item[5];
        collect { onefile -> return [ Name, RunID, LibraryID, SampleID, 
                                      SampleName ] }
    }
    .groupTuple(by : [0, 1, 2, 3, 4])
    .combine(barcodeTab, by : 0)
    .collectFile(storeDir: usedBarcodeDir) { item ->
      [ "${item[1]}_${item[2]}_${item[3]}.txt", 
      item[4] + '\t' + item[5] +  '\n']
    }

/* ----------------------------------------------------------------------------
* Perform QC check
*----------------------------------------------------------------------------*/
process qcCheck {
    label 'low_memory'

    input:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn from uniq_fqForQCtrim

    output:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn into qcFiles

    shell:
    '''
    # perform quality check on original fastq
    # full paths for R1 and R2
    R1path=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
             grep !{SampleID} | grep !{params.R1code} | \
             grep "!{params.fastqExtens}")
    R2path=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
             grep !{SampleID} | grep !{params.R2code} | \
             grep "!{params.fastqExtens}")
    fastqc $R1path $R2path --threads 2
    # As fastqc produced files are located in the sam folder as fastq reads, we
    # have to retrieve fastqc results manually
    fastqcRes=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
                grep !{SampleID} | grep -E "zip")
    # target dir
    targetDir=!{outputDir}"/fastQC/"!{RunID}"/"!{LibraryID}"/"!{SampleID}
    mkdir -p $targetDir
    mv $fastqcRes $targetDir
    '''
}

/* ----------------------------------------------------------------------------
* Trim reads by quality and adapterss with trimgalore
*----------------------------------------------------------------------------*/
process trimReads {
    label 'low_memory'

    publishDir "${outputDir}/trimmed/${RunID}/${LibraryID}/${SampleID}",
                   mode: 'copy', pattern: '*_val_*.fq.gz', overwrite: true
    publishDir "${outputDir}/trimmed/${RunID}/${LibraryID}/${SampleID}",
                   mode: 'copy', pattern: '*.{txt,zip}', overwrite: true

    input:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn from qcFiles

    output:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, 
          "${SampleID}*_val_1.fq.gz",
          "${SampleID}*_val_2.fq.gz" into trimmedFiles
    tuple path("*.txt"), path("*.zip") into trimQCfiles

    shell:
    '''
    # full paths for R1 and R2
    R1path=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
             grep !{SampleID} | grep !{params.R1code} | \
             grep "!{params.fastqExtens}")
    R2path=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
             grep !{SampleID} | grep !{params.R2code} | \
             grep "!{params.fastqExtens}")

    # perform trimming with trim galore ONLY ON R2 because we lose a lot of 
    # reads if we trim R1
    R1result=$(basename $R1path | sed 's/[.].*//g')
    R1result=$(echo $R1result"_val_1.fq.gz")
    R2result=$(basename $R2path | sed 's/[.].*//g')
    R2result=$(echo $R2result"_val_2.fq.gz")
    reportFile=!{RunID}'_'!{LibraryID}'_'!{SampleID}'_trimming_report.txt'
    cutadapt !{params.cutadapt_allParams} --minimum-length=!{R1len} \
             -o $R1result -p $R2result $R1path $R2path 1>$reportFile

    # perform fastqc
    fastqc $R1result $R2result
    '''
}


/* ----------------------------------------------------------------------------
* !!! IMPORTANT NOTE: !!!
* With use of BRB-seq tools, it is not nessecary to demultiplex and map files
* one by one in order to get to the count table. Count table can be obtained
* directly from maped trimmed file. However, it's not possible then to derive
* percentage of unmapped reads, multiple mapping percentage, etc, from it. This
* is why we still need to map the individual demultiplexed files. But it 
* doesn't make sense to do two mapping runs: for just trimmed bam and for the
* demultiplexed one, it will take twice much time. For example, just trimming +
* mapping + counting takes 14 minutes for the test file, and this is without
* demultiplexing and mapping individual files. On the other hand, trimming +
* demultiplexing + mapping + counting takes 26 minutes in total. This is due to
* the higher degree of parallelization. So, I will demultiplex and
* map and count and then constract one count table per library.
*----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
* Demultiplex reads
*----------------------------------------------------------------------------*/
process demultiplex {
    label 'mid_memory'

    publishDir "${outputDir}/demultiplexed/${RunID}/${LibraryID}/${SampleID}",
                mode: 'copy', pattern: '*.{fastq.gz,txt}', overwrite: true

    input:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, trimmedR1, 
          trimmedR2 from trimmedFiles

    output:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, trimmedR1, trimmedR2,
          path('*.fastq.gz') into demultiplexBundle
    path('*.txt') into demultiplexStats

    shell:
    '''
    # calculate UMI length
    umiLen=$((!{R1len} - 10)) 
    # get run - specific barcode file (created above) and add header to it
    barcodeFile=!{usedBarcodeDir}/!{RunID}_!{LibraryID}_!{SampleID}.txt
    echo 'Name\tB1' | cat - $barcodeFile > temp && mv temp $barcodeFile

    java -jar !{params.brbseqTools} Demultiplex -r1 !{trimmedR1} \
              -r2 !{trimmedR2} -c $barcodeFile -o "." \
              -UMI $umiLen -p !{BU_ptrn}
    '''
}

// fork into #(of demultiplexed file) channels, add sample-specific metadata
demultiplexBundle
    .flatMap { item ->
        RunID = item[0];
        LibraryID = item[1];
        SampleID = item[2];
        R1len = item[3];
        BU_ptrn = item[4];
        trimmedR1 = item[5];
        trimmedR2 = item[6];
        files = item[7];
        files.collect { onefile -> return [ RunID, LibraryID, SampleID, R1len,
                                            BU_ptrn, trimmedR1, trimmedR2,
                                            onefile ] }
    }
    .flatMap { item ->
        RunID = item[0];
        LibraryID = item[1];
        SampleID = item[2];
        R1len = item[3];
        BU_ptrn = item[4];
        trimmedR1 = item[5];
        trimmedR2 = item[6];
        demultiplexed = item[7];
        SampleName = demultiplexed.toString().replaceAll(/.fastq.gz/, '').replaceAll(/.*\//, '');
        files.collect { onefile -> return [ RunID, LibraryID, SampleID, R1len,
                                            BU_ptrn, SampleName, trimmedR1,
                                            trimmedR2, demultiplexed ] }
    }
    .unique()
    .combine(sampleTab, by : [0, 1, 2, 3, 4, 5])
    .set{ demultiplexFiles }

/* ----------------------------------------------------------------------------
* Map demultiplexed reads to reference genome with STAR
*----------------------------------------------------------------------------*/
process mapWithStar {
    label 'high_memory'

    publishDir "${outputDir}/mapped/${RunID}/${LibraryID}/${SampleID}", 
               mode: 'copy', pattern: '*.sortedByCoord.out.bam',
               overwrite: true
    publishDir "${outputDir}/mapStats/${RunID}/${LibraryID}/${SampleID}", 
                mode: 'copy', pattern: '*_Log.final.out',
                overwrite: true

    input:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, SampleName, trimmedR1,
          trimmedR2, demultiplexfq, pos, Specie, Genome from demultiplexFiles

    output:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, pos, SampleName, Specie, 
          Genome, trimmedR1, trimmedR2, demultiplexfq, 
          path('*.sortedByCoord.out.bam'),
          path('*_Log.final.out') into mappedBundle

    shell:
    '''
    mapPrefName=`basename !{demultiplexfq} | sed 's/[.].*//g'`
    mapPrefName=$mapPrefName"_"
    STAR --runMode alignReads --readFilesIn !{demultiplexfq} \
         --genomeDir !{genomePath}'/'!{Specie}'/'!{Genome}'/STAR_Index' \
         --runThreadN !{task.cpus} \
         --outFileNamePrefix $mapPrefName \
         !{params.star_allParams}
    '''
}

/* ----------------------------------------------------------------------------
* I will split the channel here, one will go to the aggregation of mapping
* statistics, and the other one - into counting reads
*----------------------------------------------------------------------------*/
mappedBundle.into{ mappedForStats; mappedForCounts }

/* ----------------------------------------------------------------------------
* Aggregats mapping statistics
* I'll create a string with all sample info and will append to it mapping stats
* taken from STAR
*----------------------------------------------------------------------------*/
process aggregateMapStats {
    label 'low_memory'

    input:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, pos, SampleName, Specie, 
          Genome, trimmedR1, trimmedR2, demultiplexfq, mappedBam, 
          mappedLog from mappedForStats

    output:
    stdout mappingStatsAggr
 
    shell:
    '''
    # initial sample info
    statsAggr=(!{RunID})
    statsAggr+=(!{LibraryID})
    statsAggr+=(!{SampleID})
    statsAggr+=(!{pos})
    statsAggr+=(!{SampleName})
    statsAggr+=(!{Specie})
    statsAggr+=(!{Genome})

    # append mapping stats info
    statsAggr+=(`grep "Number of input reads" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Uniquely mapped reads number" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads mapped to multiple loci" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads mapped to too many loci" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads unmapped: too many mismatches" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads unmapped: too short" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads unmapped: other" !{mappedLog} | sed 's/.*|//'`)
    echo "${statsAggr[@]}"
    '''
}

mappingStatsAggr
    .collectFile(name: mapStatsTab, newLine: false)

/* ----------------------------------------------------------------------------
* Count reads in demultiplexed trimmed bams
*----------------------------------------------------------------------------*/
process countReads {
    label 'high_memory'

    publishDir "${outputDir}/counts/${RunID}/${LibraryID}/${SampleID}",
               mode: 'copy', pattern: '{*.detailed.txt}', overwrite: true

    input:
    tuple RunID, LibraryID, SampleID, R1len, BU_ptrn, pos, SampleName, Specie, 
          Genome, trimmedR1, trimmedR2, demultiplexfq, mappedBam, 
          mappedLog from mappedForCounts

    output:
    tuple RunID, LibraryID, SampleID, Specie, Genome, 
          path('*.reads.detailed.txt') into readBundle
    tuple RunID, LibraryID, SampleID, Specie, Genome,
          path('*.umis.detailed.txt') into umiBundle

    shell:
    '''
    umiLen=$((!{R1len} - 10)) 
    # get run - specific barcode file (created above)
    barcodeFile=!{usedBarcodeDir}/!{RunID}_!{LibraryID}_!{SampleID}.txt

    gtfPath=`find !{genomePath}'/'!{Specie}'/'!{Genome} | grep .gtf$`
    java -jar -Xmx2g !{params.brbseqTools} CreateDGEMatrix -f !{trimmedR1} \
         -b !{mappedBam} -c $barcodeFile -o "." \
         -gtf $gtfPath -UMI $umiLen -p !{BU_ptrn} 
   
    samplName=`basename !{mappedBam} | sed 's/_Aligned.sortedByCoord.out.bam/.count/g'`
    mv output.dge.reads.detailed.txt $samplName".dge.reads.detailed.txt"
    mv output.dge.reads.txt $samplName".dge.reads.txt"
    mv output.dge.umis.detailed.txt $samplName".dge.umis.detailed.txt"
    mv output.dge.umis.txt $samplName".dge.umis.txt"
    '''
}

/* ----------------------------------------------------------------------------
* Merge count tables per sample into 1 count table: reads
*----------------------------------------------------------------------------*/
readBundle
     .groupTuple(by: [0, 1, 2, 3])
     .set{readBundleMerged}

process mergeReadCounts {
   label 'mid_memory'
   publishDir "${outputDir}/countTables",  mode: 'copy',
               pattern: '{*readsCombined.csv}', overwrite: true

   input: 
   tuple RunID, LibraryID, SampleID, Specie, Genome, 
         Reads from readBundleMerged

   output:
   file '*readsCombined.csv' into finalReadsTabs

   shell:
   '''
   function getIndexOfSubSample {
       # first of all, extract sample name
       fileName=`echo $1 | sed 's@.*/@@g'`
       # I would simply replace everything after ".", but nextflow doesn't like
       # backslash, which is used as escape character in bash
       IFS='.' read -ra SubSample <<< "$fileName"
       SubSample=${SubSample[0]}
       if [[ "${SubSample}" == "undetermined" ]]
       then
           SubSample="Unknown_Barcode"
       fi
       
       # determine in which column counts for our subsample are
       fileHeader=`grep $SubSample $1 | tr '\t' ','`
       # split so it's an array
       IFS=',' read -ra fileHeader <<< "$fileHeader"
       # determine the index
       for i in "${!fileHeader[@]}"; do
          if [[ "${fileHeader[$i]}" == "${SubSample}" ]]; then
              subSampleIndex=`echo "${i}"`;
          fi
      done
      subSampleIndex=$((subSampleIndex + 1))
      echo $subSampleIndex
   }
   
   # output file
   outputCombinedTable=`echo !{RunID} "_" !{LibraryID} "_" !{Genome} "_" !{SampleID} "_readsCombined.csv"`
   outputCombinedTable=`echo $outputCombinedTable | sed 's@ @@g'`

   # loop over all individual count files, extract column corresponding to 
   # a sample and paste it to the output file
   for oneFile in !{Reads}
   do
      oneFile=`echo $oneFile | sed "s/[^t]$//" | sed "s@^[^/]@@g"`

      subSampleInd=`getIndexOfSubSample $oneFile`
      if test -f "$outputCombinedTable"
      then
        subSampleColumnFile=`echo $subSampleInd ".txt" | sed 's@ @@g' `
        cut -f $subSampleInd $oneFile > $subSampleColumnFile
        paste -d'\t' $outputCombinedTable $subSampleColumnFile > tmp.txt
        rm $subSampleColumnFile
        mv tmp.txt $outputCombinedTable
      else
        cut -f 1,2,$subSampleInd $oneFile > $outputCombinedTable
      fi
   done
   '''
}

/* ----------------------------------------------------------------------------
* Merge count tables per sample into 1 count table: UMIs
*----------------------------------------------------------------------------*/
umiBundle
     .groupTuple(by: [0, 1, 2, 3])
     .set{umiBundleMerged}

process mergeUMICounts {
   label 'mid_memory'

   publishDir "${outputDir}/countTables",  mode: 'copy',
               pattern: '{*umisCombined.csv}', overwrite: true

   input: 
   tuple RunID, LibraryID, SampleID, Specie, Genome, 
         UMIs from umiBundleMerged

   output:
   tuple RunID, LibraryID, SampleID, Specie, Genome into forUserReport
   file '*umisCombined.csv' into finalUMIsTabs

   shell:
   '''
   function getIndexOfSubSample {
       # first of all, extract SubSample name, i.e. A01, B12, C05, etc
       fileName=`echo $1 | sed 's@.*/@@g'`
       # I would simply replace everything after ".", but nextflow doesn't like
       # backslash, which is used as escape character in bash
       IFS='.' read -ra SubSample <<< "$fileName"
       SubSample=${SubSample[0]}
       if [[ "${SubSample}" == "undetermined" ]]
       then
           SubSample="Unknown_Barcode"
       fi
       
       # determine in which column counts for our subsample are
       fileHeader=`grep $SubSample $1 | tr '\t' ','`
       # split so it's an array
       IFS=',' read -ra fileHeader <<< "$fileHeader"
       # determine the index
       for i in "${!fileHeader[@]}"; do
          if [[ "${fileHeader[$i]}" == "${SubSample}" ]]; then
              subSampleIndex=`echo "${i}"`;
          fi
      done
      subSampleIndex=$((subSampleIndex + 1))
      echo $subSampleIndex
   }
   
   # output file
   outputCombinedTable=`echo !{RunID} "_" !{LibraryID} "_" !{Genome} "_" !{SampleID} "_umisCombined.csv"`
   outputCombinedTable=`echo $outputCombinedTable | sed 's@ @@g'`

   # loop over all individual count files, extract column corresponding to 
   # a sample and paste it to the output file
   for oneFile in !{UMIs}
   do
      oneFile=`echo $oneFile | sed "s/[^t]$//" | sed "s@^[^/]@@g"`

      subSampleInd=`getIndexOfSubSample $oneFile`
      if test -f "$outputCombinedTable"
      then
        subSampleColumnFile=`echo $subSampleInd ".txt" | sed 's@ @@g' `
        cut -f $subSampleInd $oneFile > $subSampleColumnFile
        paste -d'\t' $outputCombinedTable $subSampleColumnFile > tmp.txt
        rm $subSampleColumnFile
        mv tmp.txt $outputCombinedTable
      else
        cut -f 1,2,$subSampleInd $oneFile > $outputCombinedTable
      fi
   done
   '''
}

/* ----------------------------------------------------------------------------
* Generate Rmarkdown user reports
*----------------------------------------------------------------------------*/
process generateUserReport {
  label 'low_memory'

  publishDir "${outputDir}/user_reports",  mode: 'copy',
               pattern: '{*.html}', overwrite: true

  input: 
  tuple RunID, LibraryID, SampleID, Specie, Genome from forUserReport

  output:
  file '*.html' into userReport

  shell:
  '''
  Rscript !{params.compile_report} !{params.markdown} !{user} !{pi} \
          !{outputDir} !{RunID} !{LibraryID} !{SampleID} \
          !{Specie} !{Genome} !{sampleTabPath}
  '''
}

// clean up in case of successful completion
workflow.onComplete {
    if(workflow.success){
        file('work').deleteDir()
    }
}
