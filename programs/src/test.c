//
// Created by mobin on 7/16/23.
//

#include <getopt.h>
#include "sam.h"
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"
#include "vcf.h"
#include "edlib.h"
#include <time.h>
#include <string.h>
#include "ptBlock.h"
#include "ptVariant.h"
#include "ptAlignment.h"
#include "ptMarker.h"
#include "tpool.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

bool test_sortingBlocks(){
    stHash* blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                    (void (*)(void *)) stList_destruct);

    char ctg1_name[10] = "ctg1";
    char ctg2_name[10] = "ctg2";

    // create list of unsorted blocks

    stList* ctg1_blocks = stList_construct3(0, ptBlock_destruct);
    stList* ctg2_blocks = stList_construct3(0, ptBlock_destruct);

    stList_append(ctg1_blocks, ptBlock_construct(15, 50, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(5, 15, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(10, 20, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(0, 10, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(15, 60, -1, -1, -1, -1));

    stList_append(ctg1_blocks, ptBlock_construct(50, 60, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(5, 15, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(10, 10, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(8, 8, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(0, 10, -1, -1, -1, -1));


    // add blocks to the table
    stHash_insert(blocks_per_contig, copyString(ctg1_name), ctg1_blocks);
    stHash_insert(blocks_per_contig, copyString(ctg2_name), ctg2_blocks);

    // sort blocks by rfs
    ptBlock_sort_stHash_by_rfs(blocks_per_contig); // sort in place

    // get sorted blocks
    stList* ctg1_sorted_blocks = stHash_search(blocks_per_contig, ctg1_name);
    stList* ctg2_sorted_blocks = stHash_search(blocks_per_contig, ctg2_name);

    int ctg1_sorted_start[5] = {0, 5, 10, 15 ,15};
    int ctg2_sorted_start[5] = {0, 5, 8, 10 ,50};

    bool allPassed = true;
    for(int i =0; i < 5; i ++){
        ptBlock* ctg1_block = stList_get(ctg1_sorted_blocks, i);
        ptBlock* ctg2_block = stList_get(ctg2_sorted_blocks, i);

        allPassed &= ctg1_sorted_start[i] == ctg1_block->rfs;
        allPassed &= ctg2_sorted_start[i] == ctg2_block->rfs;
    }
    return allPassed;
}




int main(int argc, char *argv[]) {
    fprintf(stdout, "Start testing ....\n");
    fprintf(stdout, "Test sorting blocks:");
    fprintf(stdout, test_sortingBlocks() ? "\x1B[32m PASSED \x1B[0m\n" : "\x1B[31m FAILED \x1B[0m\n");
}
