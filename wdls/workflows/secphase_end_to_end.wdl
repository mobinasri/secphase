version 1.0


workflow runSecphaseEndToEnd {
    call secphaseEndToEnd
    output {
        File correctedBam = secphaseEndToEnd.correctedBam
        File correctedBamIndex = secphaseEndToEnd.correctedBamIndex
        File outLog = secphaseEndToEnd.outLog
        File variantBlocksBed = secphaseEndToEnd.variantBlocksBed
        File markerBlocksBed = secphaseEndToEnd.markerBlocksBed
        File modifiedReadBlocksVariantsBed = secphaseEndToEnd.modifiedReadBlocksVariantsBed
        File modifiedReadBlocksMarkersBed = secphaseEndToEnd.modifiedReadBlocksMarkersBed
        File initalVariantBlocksBed = secphaseEndToEnd.initalVariantBlocksBed
    }
}


task secphaseEndToEnd {
    input {
        File inputBam
        File? mapqTableText
        File diploidAssemblyFastaGz
        File? phasedVcf
        File? variantBed
        String secphaseOptions = "--hifi"
        String? correctBamOptions
        String version = "v0.4.3"
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=1024
        String dockerImage="mobinasri/secphase:v0.4.3"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        BAM_FILENAME=$(basename ~{inputBam})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output

        # sort by read name
        ln -s ~{inputBam} ${BAM_PREFIX}.bam
        samtools sort -n -@~{threadCount} ${BAM_PREFIX}.bam > output/${BAM_PREFIX}.sorted_by_qname.bam
        
        # index diploid fasta
        ln -s ~{diploidAssemblyFastaGz} asm.fa.gz
        gunzip -c asm.fa.gz > asm.fa
        samtools faidx asm.fa

        # run secphase
        if [[ -n "~{phasedVcf}" ]];then
            ln -s ~{phasedVcf} phased.vcf
            ln -s ~{variantBed} variant.bed
            echo "Running variant/marker dual mode"
            secphase ~{secphaseOptions} \
                -@~{threadCount}  \
                -i output/${BAM_PREFIX}.sorted_by_qname.bam \
                -f asm.fa \
                --outDir output \
                --prefix ${BAM_PREFIX} \
                -v phased.vcf \
                -B variant.bed
        else
            echo "Running marker mode"
            secphase ~{secphaseOptions} \
                -@~{threadCount}  \
                -i output/${BAM_PREFIX}.sorted_by_qname.bam \
                -f asm.fa \
                --outDir output \
                --prefix ${BAM_PREFIX}
        fi

        # we don't need the bam sorted by read name any more
        rm -rf  output/${BAM_PREFIX}.sorted_by_qname.bam

        # make an environment variable that contains all input args for correct_bam
        OPTIONS="~{correctBamOptions} --phasingLogText output/${BAM_PREFIX}.out.log" 

        if [ -n "~{mapqTableText}" ]
        then
            OPTIONS="${OPTIONS} --mapqTable ~{mapqTableText}"
        fi

        correct_bam ${OPTIONS} -i ~{inputBam} -o output/${BAM_PREFIX}.secphase_~{version}.bam -n~{threadCount}
        samtools index -@~{threadCount} output/${BAM_PREFIX}.secphase_~{version}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File correctedBam = glob("output/*.bam")[0]
        File correctedBamIndex = glob("output/*.bam.bai")[0]
        File outLog = glob("output/*.out.log")[0]
        File initalVariantBlocksBed = glob("output/*.initial_variant_blocks.bed")[0]
        File modifiedReadBlocksVariantsBed = glob("output/*.modified_read_blocks.variants.bed")[0]
        File modifiedReadBlocksMarkersBed = glob("output/*.modified_read_blocks.markers.bed")[0]
        File markerBlocksBed = glob("output/*.marker_blocks.bed")[0]
        File variantBlocksBed = glob("output/*.variant_blocks.bed")[0]
    }
}

