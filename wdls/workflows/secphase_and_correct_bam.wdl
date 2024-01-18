version 1.0

import "secphase.wdl" as secphase_t
import "correct_bam.wdl" as correct_bam_t

workflow runSecPhaseAndCorrectBam {
    input {
        File inputBam
        File diploidAssemblyFastaGz
        File? phasedVcf
        File? variantBed
        String secphaseOptions = "--hifi"
        String secphaseDockerImage = "mobinasri/secphase:v0.4.3"
        String version = "v0.4.3"
    }
    call secphase_t.sortByName {
        input:
            bamFile = inputBam,
            dockerImage = secphaseDockerImage,
            diskSize = 7 * ceil(size(inputBam, "GB")) + 64
    }

    call secphase_t.secphase {
        input:
            bam = sortByName.outputBam,
            diploidAssemblyFastaGz = diploidAssemblyFastaGz,
            phasedVcf = phasedVcf,
            variantBed = variantBed,
            options = secphaseOptions,
            dockerImage = secphaseDockerImage,
            prefix = "secphase_${version}",
            diskSize = ceil(size(sortByName.outputBam, "GB")) + 64
    }

    call correct_bam_t.correctBam{
        input:
            bam = inputBam,
            suffix = "secphase_${version}",
            phasingLogText = secphase.outLog,
            dockerImage = secphaseDockerImage
    }

    output {
        File correctedBam = correctBam.correctedBam
        File correctedBamIndex = correctBam.correctedBamIndex
        File excludedReadIdsText = correctBam.excludedReadIdsText
        File outLog = secphase.outLog
        File variantBlocksBed = secphase.variantBlocksBed
        File markerBlocksBed = secphase.markerBlocksBed
        File modifiedReadBlocksVariantsBed = secphase.modifiedReadBlocksVariantsBed
        File modifiedReadBlocksMarkersBed = secphase.modifiedReadBlocksMarkersBed
        File initalVariantBlocksBed = secphase.initalVariantBlocksBed
    }
}

