#include "ptBlock.h"


ptBlock *ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f) {
    ptBlock *block = malloc(sizeof(ptBlock));
    block->rfs = rfs;
    block->rfe = rfe;
    block->sqs = sqs;
    block->sqe = sqe;
    block->rds_f = rds_f;
    block->rde_f = rde_f;
    block->data = NULL;
    block->destruct_data = free;
    block->copy_data = NULL;
    block->extend_data = NULL;
    return block;
}


void ptBlock_set_data(ptBlock *block, void *data, void (*destruct_data)(void *), void *(*copy_data)(void *),
                      void (*extend_data)(void *, void *)) {
    if (block->data != NULL) {
        block->free_data(block->data);
    }
    block->data = data;
    block->destruct_data = destruct_data;
    block->copy_data = copy_data;
    block->extend_data = extend_data;
}

void ptBlock_extend_data(ptBlock *block, void *data) {
    if (block->extend_data != NULL) {
        block->extend_data(block->data, data);
    }
}

void ptBlock_destruct_data(ptBlock *block) {
    if (block->destruct_data != NULL) {
        block->destruct_data(block->data);
    }
}

void *ptBlock_copy_data(ptBlock *block) {
    if (block->copy_data != NULL) {
        return block->copy_data(block->data);
    } else {
        return NULL;
    }
}

ptBlock *ptBlock_copy(ptBlock *block) {
    ptBlock *block_copy = ptBlock_construct(block->rfs,
                                            block->rfe,
                                            block->sqs,
                                            block->sqe,
                                            block->rds_f,
                                            block->rde_f);
    void *data_copy = ptBlock_copy_data(block);
    ptBlock_set_data(block_copy, data_copy, block->destruct_data, block->copy_data, block->extend_data);
    return block_copy;
}


void ptBlock_destruct(ptBlock *block) {
    ptBlock_destruct_data(block);
    free(block);
}


/// Get Functions ////
int ptBlock_get_rfs(ptBlock *block) {
    return block->rfs;
}


int ptBlock_get_rfe(ptBlock *block) {
    return block->rfe;
}

int ptBlock_get_sqs(ptBlock *block) {
    return block->sqs;
}

int ptBlock_get_sqe(ptBlock *block) {
    return block->sqe;
}

int ptBlock_get_rds_f(ptBlock *block) {
    return block->rds_f;
}

int ptBlock_get_rde_f(ptBlock *block) {
    return block->rde_f;
}


//// Set Functions ////

void ptBlock_set_rfs(ptBlock *block, int rfs) {
    block->rfs = rfs;
}

void ptBlock_set_rfe(ptBlock *block, int rfe) {
    block->rfe = rfe;
}

void ptBlock_set_sqs(ptBlock *block, int sqs) {
    block->sqs = sqs;
}

void ptBlock_set_sqe(ptBlock *block, int sqe) {
    block->sqe = sqe;
}

void ptBlock_set_rds_f(ptBlock *block, int rds_f) {
    block->rds_f = rds_f;
}

void ptBlock_set_rde_f(ptBlock *block, int rde_f) {
    block->rde_f = rde_f;
}

//// Comparing Functions ////
int ptBlock_cmp_rfs(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->rfs - b2->rfs;
}

int ptBlock_cmp_rds_f(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->rds_f - b2->rds_f;
}

int ptBlock_cmp_sqs(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->sqs - b2->sqs;
}


stHash *ptBlock_parse_bed(char *bed_path) {
    FILE *fp = fopen(bed_path, "r+");
    size_t read;
    size_t len;
    char *line = NULL;
    char *contig_name;
    char *token;
    ptBlock *block = NULL;
    stList *blocks = NULL;
    stHash *blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                  (void (*)(void *)) stList_destruct);
    while ((read = getline(&line, &len, fp)) != -1) {
        // replace '\n' by '\0'
        if (line[strlen(line) - 1] == '\n') {
            line[strlen(line) - 1] = '\0';
        }
        token = strtok(line, "\t");
        contig_name = malloc(strlen(token) + 1);
        strcpy(contig_name, token);
        token = strtok(NULL, "\t");
        start = atoi(token); // 0-based
        token = strtok(NULL, "\t");
        end = atoi(token) - 1; // 0-based
        block = ptBlock_construct(start, end, -1, -1, -1, -1);
        blocks = stHash_search(blocks_per_contig, contig_name);
        // if contig does not exist as a key in the table
        // add contig along with an empty list as its value
        if (blocks == NULL) {
            blocks = stList_construct3(0, ptBlock_destruct);
            stHash_insert(blocks_per_contig, contig_name, blocks);
        }
        stList_append(blocks, block);
    }
    // sort blocks per contig
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        stList_sort(blocks, ptBlock_cmp_rfs);
    }
    return blocks_per_contig;
}

stList *ptBlock_merge_blocks(stList *blocks,
                             int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                             void (*set_end)(ptBlock *, int)) {
    stList *blocks_merged = stList_construct3(0, ptBlock_destruct);
    if (stList_length(blocks) == 0) return blocks_merged;
    ptBlock *b = NULL;
    ptBlock *b_merged = NULL;
    for (int i = 0; i < stList_length(blocks); i++) {
        b = stList_get(blocks, i);
        //printf("%d\t%d\n", b->rds_f, b->rde_f);
        if (i == 0) { // Initiate b_merged for the first block
            b_merged = ptBlock_copy(b);
            continue;
        }
        if (get_end(b_merged) < get_start(b)) {// no overlap with previous merged block
            //save the merged block
            stList_append(blocks_merged, b_merged);
            // Initiate a new merged block
            b_merged = ptBlock_copy(b);
        } else { //there is overlap
            set_end(b_merged, get_end(b)); // extend end pos of the merged block
            if (b->data != NULL) {
                // add the data of the new block to the merged block
                ptBlock_extend_data(b_merged, b->data);
            }
        }
    }

    // Add the last merged block
    if (b_merged) {
        stList_append(blocks_merged, b_merged);
    }

    return blocks_merged;
}

stHash *ptBlock_merge_blocks_per_contig(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                        void (*set_end)(ptBlock *, int)) {
    char *contig_name;
    stList *blocks;
    stList *merged_blocks;
    stHash *merged_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                         (void (*)(void *)) stList_destruct);
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        // get blocks
        blocks = stHash_search(blocks_per_contig, contig_name);
        // merge blocks
        ptBlock_merge_blocks(blocks, get_start, get_end, set_end);
        // add merged blocks to the new table
        stHash_insert(merged_blocks_per_contig, contig_name, merged_blocks);
    }
    return merged_blocks_per_contig;
}

stHash *ptBlock_merge_blocks_per_contig_by_rf(stHash *blocks_per_contig){
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_rfs,
                                           ptBlock_get_rfe,
                                           ptBlock_set_rfe);
}

stHash *ptBlock_merge_blocks_per_contig_by_rd_f(stHash *blocks_per_contig){
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_rds_f,
                                           ptBlock_get_rde_f,
                                           ptBlock_set_rde_f);
}

stHash *ptBlock_merge_blocks_per_contig_by_sq(stHash *blocks_per_contig){
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_sqs,
                                           ptBlock_get_sqe,
                                           ptBlock_set_sqe);
}

int ptBlock_get_total_length(stHash *blocks_per_contig, int (*get_start)(ptBlock *), int (*get_end)(ptBlock *)){
    int total_len = 0;
    char *contig_name;
    stList *blocks;
    ptBlock* block;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        for (int i=0; i < stList_length(blocks); i++){
            block = stList_get(blocks, i);
            total_len += get_end(block) - get_start(block);
        }
    }
    return total_len;
}

int ptBlock_get_total_length_by_rf(stHash *blocks_per_contig){
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rfs, ptBlock_get_rfe);
}

int ptBlock_get_total_length_by_rd_f(stHash *blocks_per_contig){
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rds_f, ptBlock_get_rde_f);
}

int ptBlock_get_total_length_by_sq(stHash *blocks_per_contig){
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_sqs, ptBlock_get_sqe);
}