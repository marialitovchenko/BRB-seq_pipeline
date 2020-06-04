#!/usr/bin/env nextflow

def helpMessage() {

    log.info """
    -      \033[41m B R B - s e q   N E X T F L O W  I N D E X  G E N O M E v1.0\033[0m-
    ================================================================================
    Welcome to the Nextflow BRB-seq genome indexing pipeline. This pipeline can be 
    used in two modes: 1) creating and indexing a "usual" = from ENSEML genome and
    2) creating and indexing a custom genome, i.e. then you have plasmids included
    into the genome or specie-mixed genome. For the mode 1 you don't need to 
    download anything, the download will be done automatically. For the mode 2 you
    need to provide fasta file and gtf file. 

    Usage:
    The \033[1;91mtypical\033[0m command for running the pipeline is as follows:
    nextflow forTest.nf \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--genomeDir\033[0m folderWithGenomes 

    or

    to put the pipeline into \033[1;91mbackground\033[0m mode:
    nextflow forTest.nf \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--genomeDir\033[0m folderWithGenomes \\
                        \033[1;91m-bg\033[0m \\
                        \033[1;91m-N\033[0m your.email@gmail.com

    \033[1;91mMandatory\033[0m arguments:
      \033[1;91m--inputTab\033[0m        Path to the table containing information about input
                        data.
                        In \033[93m mode 1 \033[0m (creating and indexing a "usual" = from ENSEML)
                        table should have just two columns: Specie (i.e. homo_sapience),
                        and GenomeCode (i.e. hs19, hs38, ... , etc; check proper code 
                        for your specie of interest on
                        https://hgdownload.soe.ucsc.edu/downloads.html).
                        In \033[93m mode 2 \033[0m (creating and indexing a custom genome) table
                        should have 4 columns: Specie (i.e. homo_sapience), GenomeCode 
                        (desired name, i.e. MixedHsMM), Fasta (name of the fasta file for
                        your genome, i.e. "MixedHsMM.fa", it must be located in the 
                        genomeDir/specie_name/GenomeCode 
                        (i.e. folderWithGenomes/homo_sapience/myCustomGenome) 
                        and GTF (name of the GTF file for your genome, i.e. "MixedHsMM.gtf",
                        it must be located in the genomeDir/specie_name
                        (i.e. folderWithGenomes/homo_sapience/myCustomGenome)
                        \033[93m You can mix two modes in one table, i.e. \033[0m
                        Specie  GenomeCode  Fasta   GTF
                        homo_sapience   hg19
                        homo_sapience   myCustomGenome    customGenome.fa    customGenome.gtf
      \033[1;91m--genomeDir\033[0m       Directory
                        In \033[93m mode 1 \033[0m (creating and indexing a "usual" = from ENSEML):
                        just a directory name. It will be created.
                        In \033[93m mode 2 \033[0m (creating and indexing a custom genome):
                        a directory containing specie_name directory containing Fasta and 
                        GTF files of your custom genome (i.e. folderWithGenomes/homo_sapience) 

    \033[1;91mOptional\033[0m arguments:
    This arguments are not going to be needed with use of graphical user
    interface
      \033[1;91m--help\033[0m            Displays this message
      \033[1;91m-bg\033[0m               Puts execution of the pipeline into background mode
      \033[1;91m-N\033[0m                email adress in order to get notified upon pipeline complition.
                        Do not use epfl email address, because emails can't pass firewall.
                        Use gmail.
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
// path to the input table (with columns in dependence with mode )
genomeTabPath = file(params.inputTab)

// path to fasta file with sequences of markers, i.e. GFP
params.markerFasta = file('.')
markerFastaPath = file(params.markerFasta)

// genomes directory (basically output folder): directory where indexed genome
// will be put to
params.genomeDir = file('.')
genomePath = file(params.genomeDir)

// common inputs: path to UCSC genomes
goldenPath="https://hgdownload.soe.ucsc.edu/goldenPath/"
// path to picard jar
picardJar="/home/litovche/bin/picard.jar"

threadsNumb=6

/* ----------------------------------------------------------------------------
* Read input table
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
genomeTabCh = Channel.fromPath( genomeTabPath )
genomeTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.Specie, row.GenomeCode, row.Fasta, row.GTF) }
    .into{ genomeTab_download; genomeTab_custom }

// separate genomes which needed to be downloaded from ensemle 
// (genomeTab_download_flt) from the custom ones
genomeTab_download
    .filter{ it[2] == null }
    .filter{ it[3] == null}
    .set{ genomeTab_download_flt }

genomeTab_custom
    .filter{ it[2] != null }
    .filter{ it[3] != null}
    .set{genomeTab_custom_flt}

/* ----------------------------------------------------------------------------
* Download the genome if user wants new specie/version which we don't have
*----------------------------------------------------------------------------*/
process downloadGenome {
    publishDir "${genomePath}/${Specie}/${GenomeCode}",
                mode: 'copy', overwrite: true

    input:
    tuple Specie, GenomeCode, Fasta, GTF from genomeTab_download_flt

    output:
    tuple Specie, GenomeCode, "*.fa", 
          "*.gtf" optional true into genomes_ensembl

    shell:
    '''
    # if both are null: this is the case to download from ENSEMBL
    if [ !{Fasta} = "null" ] && [ !{GTF} = "null" ]; then
        # Download genome from UCSC
        wget !{goldenPath}!{GenomeCode}/bigZips/!{GenomeCode}.fa.gz
        gunzip !{GenomeCode}.fa.gz

        # Download gene annotation file
        urlBase=!{goldenPath}!{GenomeCode}/bigZips/genes/!{GenomeCode}
        # the best to use with UCSC genome is refGene, but some species don't
        # have it. In That case one needs to try other annotations
        annoTypes=(refGene ensGene ncbiRefSeq)
        # will try until something is downloaded
        wget "$urlBase".refGene.gtf.gz || continue
        if test -f !{GenomeCode}.refGene.gtf.gz; then
            gunzip !{GenomeCode}.refGene.gtf.gz
        else
            wget "$urlBase".ensGene.gtf.gz || continue
            if test -f !{GenomeCode}.ensGene.gtf.gz; then
                gunzip !{GenomeCode}.ensGene.gtf.gz
            else
                wget "$urlBase".ncbiRefSeq.gtf.gz || continue
                if test -f !{GenomeCode}.ncbiRefSeq.gtf.gz; then
                    gunzip !{GenomeCode}.ncbiRefSeq.gtf.gz
                fi
            fi
        fi
    fi
    '''
}

genomeTab_custom_flt
    .map { item ->
        Specie = item[0];
        GenomeCode = item[1];
        Fasta = genomePath.toString() + "/" + item[0] + "/" + 
                item[1] + "/" + item[2];
        GTF = genomePath.toString() +  "/" + item[0] + "/" + 
                item[1] + "/" + item[3];
        return [ Specie, GenomeCode, Fasta, GTF ] 
    }
    .mix(genomes_ensembl)
    .set{allGenomes}

/* ----------------------------------------------------------------------------
* [Optional] create custom GTF for markers fasta
*----------------------------------------------------------------------------*/
if( markerFastaPath != file('.') ) {
    // get the legth of each fasta sequence marker
    Channel
     .fromPath(markerFastaPath)
     .splitFasta( record: [id: true, seqString: true ])
     .map{ record -> tuple(record.id, record.seqString.length()) }
     .set{markerFasta}

     // create a GTF file corresponding to the markers
     process createMarkerGTF {
        echo true

        input:
        tuple markerID, markerSeqLen from markerFasta

        output:
        path('*.gtf') into markerGTF

        shell:
        '''
        geneID_exonID='gene_id "'!{markerID}'-gene"; transcript_id "'!{markerID}'-tr"; exon_number "1"; '
        geneName_biotype='gene_name "'!{markerID}'-gene"; gene_source "user"; gene_biotype "protein_coding";'
        transcriptid_Source=' transcript_name "'!{markerID}'-001"; transcript_source "user";'

        echo -e !{markerID}' \t 'protein_coding' \t 'gene' \t '1' \t '!{markerSeqLen}' \t '.' \t '+' \t '.' \t ''gene_id "'!{markerID}'-gene"; '$geneName_biotype > !{markerID}'.gtf'
        echo -e !{markerID}' \t 'protein_coding' \t 'transcript' \t '1' \t '!{markerSeqLen}' \t '.' \t '+' \t '.' \t ''gene_id "'!{markerID}'-gene"; transcript_id "'!{markerID}'-tr"; '$geneName_biotype$transcriptid_Source >> !{markerID}'.gtf'
        echo -e !{markerID}' \t 'protein_coding' \t 'exon' \t '1' \t '!{markerSeqLen}' \t '.' \t '+' \t '.' \t '$geneID_exonID$geneName_biotype$transcriptid_Source' exon_id "'!{markerID}'-exon";' >> !{markerID}'.gtf'
        echo -e !{markerID}' \t 'protein_coding' \t 'CDS' \t '1' \t '!{markerSeqLen}' \t '.' \t '+' \t '0' \t '$geneID_exonID$geneName_biotype$transcriptid_Source' protein_id "'!{markerID}'-protein";' >> !{markerID}'.gtf'
        echo -e !{markerID}' \t 'protein_coding' \t 'start_codon' \t '1' \t '3' \t '.' \t '+' \t '0' \t '$geneID_exonID$geneName_biotype$transcriptid_Source >> !{markerID}'.gtf'
        echo -e !{markerID}' \t 'protein_coding' \t 'stop_codon' \t '$((!{markerSeqLen} - 2))' \t '!{markerSeqLen}' \t '.' \t '+' \t '0' \t '$geneID_exonID$geneName_biotype$transcriptid_Source >> !{markerID}'.gtf'
        '''
     }

    process addMarkerToRefGen {
        input:
        path markerGTF from markerGTF
        tuple Specie, GenomeCode, Fasta, GTF from allGenomes

        output:
        tuple Specie, GenomeCode, path('*.fa'), path('*.gtf') into genomesToIndex1

        shell:
        '''
        refGenFaWithMarkers=$(basename !{Fasta} | sed 's/.*//g')
        refGenFaWithMarkers=$(echo $refGenFaWithMarkers'_withMarkers.fa')

        refGenGTFwithMarkers=$(basename !{Fasta} | sed 's/.*//g')
        refGenGTFwithMarkers=$(echo $refGenGTFwithMarkers'_withMarkers.gtf')

        cat !{Fasta} !{markerFastaPath} > $refGenFaWithMarkers
        cat !{GTF} !{markerGTF} > $refGenGTFwithMarkers
        '''
    }
}

/* ----------------------------------------------------------------------------
* Index genome for use with STAR
*----------------------------------------------------------------------------*/
process indexGenomeForSTAR {
    publishDir "${genomePath}/${Specie}/${GenomeCode}",
                mode: 'copy', overwrite: true

    input:
    tuple Specie, GenomeVersion, Fasta, GTF from genomesToIndex

    shell:
    '''
    # index with STAR
    STAR --runMode genomeGenerate --runThreadN "!{threadsNumb}" --genomeDir "!{genomeDir}" \
    --genomeFastaFiles "!{genomeFasta}" --sjdbGTFfile "!{genomeGTF}"
    # created files: chrLength.txt, chrStart.txt, exonGeTrInfo.tab,
    # genomeParameters.txt, SA, sjdbList.fromGTF.out.tab,
    # chrNameLength.txt, exonInfo.tab, geneInfo.tab, Log.out, SAindex,
    # sjdbList.out.tab, chrName.txt, Genome, sjdbInfo.txt, transcriptInfo.tab

    # index with samtools
    samtools faidx "!{genomeFasta}"
    # created files :  *.fai

    # index picard
    dictFile=$(echo "!{genomeFasta}" | sed "s/fa$/dict/g")
    java -jar "!{picardJar}" CreateSequenceDictionary \
              REFERENCE="!{genomeFasta}" OUTPUT=$dictFile
    # created files: *.dict
    '''
}
