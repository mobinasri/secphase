#ifndef PT_BLOCK_H
#define PT_BLOCK_H

#include <assert.h>
#include <math.h>
#include <float.h>
#include "common.h"

/*
 * @abstract Structure for saving a block
 * @field rfs           Start coordinate on reference (0-based inclusive)
 * @field rfe           End coordinate on reference (0-based inclusive)
 * @field sqs           Start coordinate on sequence (0-based inclusive)
 * @field sqe           End coordinate on sequence (0-based inclusive)
 * @field rds_f         Start coordinate on read's forward strand (0-based inclusive)
 * @field rde_f         End coordinate on read's forward strand (0-based inclusive)
 * @field data		Pointer to some data attached to the block
 * @field free_data	Function for freeing the memory allocated to data
 */
typedef struct {
        int rfs; // ref start
        int rfe; // ref end
        int sqs; // seq start
        int sqe; // seq end
        int rds_f;
        int rde_f;
        void* data;
        void (*free_data)(void*);
}ptBlock;


/* Construct a ptBlock structure
 * Note that this constructor only receives the coordinates
 * and related data should be added by another function; ptBlock_add_data()
 */
ptBlock* ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f);


/* Add data to a ptBlock structure. 
 * It also receives a function for freeing the memory allocated to the data.
 */
void ptBlock_add_data(ptBlock* block, void* data, void (*free_data)(void*));


/* Make a copy of a ptBlock structure
 * It takes a function for copying data
 */
ptBlock* ptBlock_copy(ptBlock* block, void* (*copy_data)(void*));


/* Destruct a ptBlock structure
 * It first frees the data augmented to the block
 * then frees the ptBlock struct
*/
void ptBlock_destruct(ptBlock* block);


// Compare two blocks by rds_f
int ptBlock_cmp_rds_f(const void *a, const void *b);

// Compare two blocks by sqs
int ptBlock_cmp_sqs(const void *a, const void *b);

#endif /* PT_BLOCK_H */

