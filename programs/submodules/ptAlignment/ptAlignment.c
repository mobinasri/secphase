#include "ptAlignment.h"


ptAlignment *ptAlignment_construct(bam1_t *record, sam_hdr_t *sam_hdr) {
    ptAlignment *alignment = (ptAlignment *) malloc(sizeof(ptAlignment));
    alignment->record = bam_init1();
    assert(bam_copy1(alignment->record, record) != NULL);
    strcpy(alignment->contig, sam_hdr_tid2name(sam_hdr, record->core.tid));
    alignment->score = 0.0;
    alignment->conf_blocks = NULL;
    alignment->flank_blocks = NULL;
    ptAlignment_init_coordinates(alignment);
    return alignment;
}

void ptAlignment_init_coordinates(ptAlignment *alignment) {
    // initialize to -1
    alignment->rfs = -1;
    alignment->rfe = -1;
    alignment->rde_f = -1;
    alignment->rds_f = -1;

    bam1_t *b = alignment->record;
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    uint8_t *quality = bam_get_qual(b);
    while (ptCigarIt_next(cigar_it)) {
        //set the start coordinates of the alignment
        if (alignment->rfs == -1 &&
            (cigar_it->op == BAM_CMATCH ||
             cigar_it->op == BAM_CEQUAL ||
             cigar_it->op == BAM_CDIFF)) {
            alignment->rfs = cigar_it->rfs;
            if (bam_is_rev(b)) {
                alignment->rde_f = cigar_it->rde_f;
            } else {
                alignment->rds_f = cigar_it->rds_f;
            }
        }
        //set the end coordinates of the alignment
        //if the alignment ends with hard or soft clipping
        //alignment->rfs != -1 is to make sure we have reached
        //the end of the alignment
        if (alignment->rfe == -1 &&
            alignment->rfs != -1 &&
            (cigar_it->op == BAM_CHARD_CLIP ||
             cigar_it->op == BAM_CSOFT_CLIP)) {
            alignment->rfe = cigar_it->rfe;
            if (bam_is_rev(b)) {
                alignment->rds_f = cigar_it->rde_f + 1;
            } else {
                alignment->rde_f = cigar_it->rds_f - 1;
            }
        }
    }
    //set the end coordinates of the alignment
    //if the alignment ends with mis/matches
    if (alignment->rfe == -1 &&
        (cigar_it->op == BAM_CMATCH ||
         cigar_it->op == BAM_CEQUAL ||
         cigar_it->op == BAM_CDIFF)) {
        alignment->rfe = cigar_it->rfe;
        if (bam_is_rev(b)) {
            alignment->rds_f = cigar_it->rds_f;
        } else {
            alignment->rde_f = cigar_it->rde_f;
        }
    }
}

void ptAlignment_destruct(ptAlignment *alignment) {
    if (alignment->record) {
        bam_destroy1(alignment->record);
        alignment->record = NULL;
    }
    if (alignment->conf_blocks) {
        stList_destruct(alignment->conf_blocks);
        alignment->conf_blocks = NULL;
    }
    if (alignment->flank_blocks) {
        stList_destruct(alignment->flank_blocks);
        alignment->flank_blocks = NULL;
    }
    free(alignment);
}

bool ptAlignment_contain_supp(ptAlignment **alignments, int alignments_len) {
    for (int i = 0; i < alignments_len; i++) {
        if (alignments[i]->record->core.flag & BAM_FSUPPLEMENTARY) return 1;
    }
    return 0;
}


void print_contigs(ptAlignment **alignments, int alignments_len) {
    DEBUG_PRINT("\nalignments:\n");
    DEBUG_PRINT("#alignment_idx\tcontig_name\t\n");
    for (int t = 0; t < alignments_len; t++) {
        DEBUG_PRINT("%d\t%s\n", t, alignments[t]->contig);
    }
}

int get_best_record_index(ptAlignment **alignments, int alignments_len, double prim_margin, double min_score,
                          double prim_margin_random) {
    assert(alignments_len > 0);
    if (alignments_len == 1) return 0;
    double max_score = -DBL_MAX;
    int max_idx = -1;
    double prim_score = -DBL_MAX;
    int prim_idx = -1;
    for (int i = 0; i < alignments_len; i++) {
        if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0) {
            prim_idx = i;
            prim_score = alignments[i]->score;
        } else if (max_score < alignments[i]->score) {
            max_idx = i;
            max_score = alignments[i]->score;
        }
    }
    int rnd = rand() % 2; // 50% chance for rnd=0 (same for rnd=1)
    if (abs(max_score - prim_score) < prim_margin_random) {
        return rnd == 0 ? prim_idx : max_idx; // if the scores were closer than prim_margin_random return one randomly
    }
    if (prim_idx == -1 || max_score < (prim_score + prim_margin) || max_score < min_score) return prim_idx;
    else return max_idx;
}

