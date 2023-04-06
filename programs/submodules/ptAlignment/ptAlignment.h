#ifndef PT_ALIGNMENT_H
#define PT_ALIGNMENT_H

#include "sam.h"
#include "faidx.h"
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"

/*! @typedef
 * @abstract Structure for an alignment record and its marker consistency score
 * @field record	Pointer to a bam1_t struct that contains the alignment record
 * @field score		marker consistency score (summation of the base quality values of the markers)
 * @field conf_blocks	confident blocks (blocks around markers. They do not contain any insertions longer than a threshold. 
 * 			They are useful for faster BAQ calculation)
 * @field flank_blocks	flanking blocks (windows of a specific length around markers. 
 * 			They are useful for faster BAQ calculation)
 * @field rfs		start coordinate on reference (0-based inclusive)
 * @field rfe		end coordinate on reference (0-based inclusive)
 * @field rds_f		start coordinate on read's forward strand (0-based inclusive)
 * @field rde_f		end coordinate on read's forward strand (0-based inclusive)
 */
typedef struct {
        bam1_t* record;
	char contig[50];
        double score;
        stList* conf_blocks;
        stList* flank_blocks;
        int rfs;
        int rfe;
        int rds_f;
        int rde_f;
}ptAlignment;

// Construct a ptAlignment struct
/*
 * @param record 	Alignment record
 * @param score		Initial score of the alignment (usually just 0.0)
 * @return alignment	Alignment saved as a ptAlignment struct
 *
 */
ptAlignment* ptAlignment_construct(bam1_t* record, double score);


// Destruct a ptAlignment struct
/*
 * @param alignment        Alignment saved as a ptAlignment struct 
 *
 */
void ptAlignment_destruct(ptAlignment* alignment);

// Check if there is a supplementary alignment among 
// the given list of alignments
bool ptAlignment_contain_supp(ptAlignment** alignments, int alignments_len);


void print_contigs(ptAlignment** alignments, int alignments_len, sam_hdr_t* h);


#endif /* PT_ALIGNMENT_H */