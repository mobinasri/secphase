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
    block->free_data = free;
    return block;
}


void ptBlock_add_data(ptBlock *block, void *data, void (*free_data)(void *)) {
    block->data = data;
    block->free_data = free_data;
}


ptBlock *ptBlock_copy(ptBlock *block, void *(*copy_data)(void *)) {
    ptBlock *block_copy = ptBlock_construct(block->rfs,
                                            block->rfe,
                                            block->sqs,
                                            block->sqe,
                                            block->rds_f,
                                            block->rde_f);
    void *data_copy = copy_data(block->data);
    ptBlock_add_data(block_copy, data_copy, block->free_data);
    return block_copy;
}


void ptBlock_destruct(ptBlock *block) {
    block->free_data(block->data);
    free(block);
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
