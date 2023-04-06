#include "ptAlignment.h"


ptAlignment* ptAlignment_construct(bam1_t* record, sam_hdr_t* sam_hdr){
        ptAlignment* alignment = (ptAlignment*) malloc(sizeof(ptAlignment));
        alignment->record = bam_init1();
        assert(bam_copy1(alignment->record, record) != NULL);
	strcpy(alignment->contig, sam_hdr_tid2name(sam_hdr, record->core.tid));
        alignment->score = 0.0;
        alignment->conf_blocks = NULL;
        alignment->flank_blocks = NULL;
        alignment->rfs = -1;
        alignment->rfe = -1;
        alignment->rds_f = -1;
        alignment->rde_f = -1;
        return alignment;
}

void ptAlignment_destruct(ptAlignment* alignment){
        if (alignment->record){
                bam_destroy1(alignment->record);
                alignment->record = NULL;
        }
        if (alignment->conf_blocks){
                stList_destruct(alignment->conf_blocks);
                alignment->conf_blocks = NULL;
        }
        if (alignment->flank_blocks){
                stList_destruct(alignment->flank_blocks);
                alignment->flank_blocks = NULL;
        }
        free(alignment);
}

bool contain_supp(ptAlignment** alignments, int alignments_len){
        for(int i=0;i < alignments_len; i++){
                if (alignments[i]->record->core.flag & BAM_FSUPPLEMENTARY) return 1;
        }
        return 0;
}


void print_contigs(ptAlignment** alignments, int alignments_len){
        DEBUG_PRINT("\nalignments:\n");
        DEBUG_PRINT("#alignment_idx\tcontig_name\t\n");
        for(int t =0 ; t < alignments_len; t++) {
                DEBUG_PRINT("%d\t%s\n", t, alignments[t]->contig);
        }
}

