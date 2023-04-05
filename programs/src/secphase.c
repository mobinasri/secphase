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

#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))

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
ptAlignment* ptAlignment_construct(bam1_t* record, double score){
        ptAlignment* alignment = (ptAlignment*) malloc(sizeof(ptAlignment));
        bam1_t* record_cpy = bam_init1();
        assert(bam_copy1(record_cpy, record) != NULL);
        alignment->record = record_cpy;
        alignment->score = score;
        alignment->conf_blocks = NULL;
        alignment->flank_blocks = NULL;
        alignment->rfs = -1;
        alignment->rfe = -1;
        alignment->rds_f = -1;
        alignment->rde_f = -1;
        return alignment;
}


// Destruct a ptAlignment struct
/*
 * @param alignment        Alignment saved as a ptAlignment struct 
 *
 */
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


/*! @typedef
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
ptBlock* ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f){
        ptBlock* block = malloc(sizeof(ptBlock));
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


/* Add data to a ptBlock structure. 
 * It also receives a function for freeing the memory allocated to the data.
 */
void ptBlock_add_data(ptBlock* block, void* data, void (*free_data)(void*)){
        block->data = data;
        block->free_data = free_data;
}

/* Make a copy of a ptBlock structure
 * It takes a function for copying data
 */
ptBlock* ptBlock_copy(ptBlock* block, void* (*copy_data)(void*)){
	ptBlock* block_copy = ptBlock_construct(block->rfs, 
			                        block->rfe,
						block->sqs, 
						block->sqe, 
						block->rds_f, 
						block->rde_f);
	void* data_copy = copy_data(block->data);
	ptBlock_add_data(block_copy, data_copy, block->free_data);
	return block_copy;
}

/* Destruct a ptBlock structure
 * It first frees the data augmented to the block
 * then frees the ptBlock struct
*/
void ptBlock_destruct(ptBlock* block){
        block->free_data(block->data);
        free(block);
}

int ptBlock_cmp_rds_f(const void *a, const void *b){
        ptBlock* b1 = (ptBlock*) a;
        ptBlock* b2 = (ptBlock*) b;
        return b1->rds_f - b2->rds_f;
}

int ptBlock_cmp_sqs(const void *a, const void *b){
        ptBlock* b1 = (ptBlock*) a;
        ptBlock* b2 = (ptBlock*) b;
        return b1->sqs - b2->sqs;
}


//#define DEBUG

#ifdef DEBUG
//void DEBUG_PRINT(const char *, ...);
void DEBUG_PRINT(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

#else
static inline void DEBUG_PRINT(const char *fmt, ...) {};
#endif



/*! @typedef
 @abstract Structure for a variant record
 @field contig		Name of the contig where the variant is located
 @field pos		Position of the variant in contig
 @field type		Type of the variant; can be 
 			(Look at https://github.com/samtools/htslib/blob/develop/htslib/vcf.h)
 @field vaf		Variant allele frequency
 @field gq		Genotype quality
 @field gt              Genotype list
 @field gt_len		Length of gt
 @field alleles		allele sequences (starts with the REF allele)
 @field longest_allele_len	Length of the longest allele (could be insertion, deletion or snp)
 @field ps              Phase block number
 */
typedef struct {
        char contig[50];
        int32_t	pos;
	int8_t	type;
        float	vaf;
        int8_t	gq;
       	int8_t*	gt;
	int32_t gt_len;
	stList* alleles;
	int32_t longest_allele_len;
        int64_t ps;
}ptVariant;


/// Create a ptVariant structure
/**
 * Parameters are listed in the definition of the structure
 * @return 	Pointer to a new variant struct on success; NULL on failure
 * 		The ptMarker struct returned by a successful call should be freed
 * 		via ptMarker_destruct() when it is no longer needed.
 */
ptVariant* ptVariant_construct(char* contig, int32_t pos, int8_t  type, float vaf, int8_t  gq, int64_t ps){
        ptVariant* variant = (ptVariant*) malloc(sizeof(ptVariant));
        strcpy(variant->contig, contig);
        variant->pos = pos;
        variant->type = type;
        variant->vaf = vaf;
        variant->gq = gq;
        variant->ps = ps;
	// Initialize gt and alleles variables
	// They can be filled later with ptMarker_append_gt() and ptMarker_append_allele()
	variant->gt = NULL;
	variant->gt_len = 0;
	variant->alleles = stList_construct3(0, free);
	variant->longest_allele_len = 0;
        return variant;
}


/// Free a ptVariant structure
/*
 * @param variant       Pointer to the variant record
 */
void ptVariant_destruct(ptVariant* variant){
        free(variant->gt);
        stList_destruct(variant->alleles);
        free(variant);
}



int ptVariant_cmp(const void *a, const void *b){
        ptVariant* v1 = (ptVariant*) a;
        ptVariant* v2 = (ptVariant*) b;
	// first compare by contig then pos
	if (strcmp(v1->contig, v2->contig) !=0){
		return strcmp(v1->contig, v2->contig);
	}
	else {
        	return v1->pos - v2->pos;
	}
}

/// Add an allele sequence to a variant record
/**
   @param variant	Pointer to the variant record
   @param allele	Allele sequence to be added

   Note that alleles should be added in the correct order
   starting from the REF allele to be consistent with the
   genotype indices
 */
void ptVariant_append_allele(ptVariant* variant, char* allele){
	char* allele_copy = (char*) malloc((strlen(allele) + 1) * sizeof(char));
	strcpy(allele_copy, allele);
	stList_append(variant->alleles, allele_copy);
	
	// update longest allele len
	if (variant->longest_allele_len < strlen(allele)){
		variant->longest_allele_len = strlen(allele);
	}
}

/// Add a genotype to a variant record
/**
   @param variant       Pointer to the variant record
   @param gt        	genotype index to be added
 */
void ptVariant_append_gt(ptVariant* variant, int8_t gt){
	// Increase size of gt list
	if(variant->gt == NULL){
		variant->gt = malloc(1 * sizeof(int8_t));
	}else{
		variant->gt = realloc(variant->gt, (variant->gt_len + 1) * sizeof(int8_t));
	}
	variant->gt[variant->gt_len] = gt;
	variant->gt_len += 1;
}

ptVariant* ptVariant_copy(ptVariant* src){
        ptVariant* dest = (ptVariant*) malloc(sizeof(ptVariant));
        strcpy(dest->contig, src->contig);
        dest->pos = src->pos;
        dest->type = src->type;
        dest->vaf = src->vaf;
        dest->gq = src->gq;
        dest->ps = src->ps;
        dest->gt = (int8_t*) malloc(src->gt_len * sizeof(int8_t));
        dest->gt_len = src->gt_len;
	memcpy(dest->gt, src->gt, src->gt_len);
        dest->alleles = stList_construct3(0, free);
	for(int i=0; i < stList_length(src->alleles); i++){
		char* allele = stList_get(src->alleles, i);
		char* allele_copy = (char*) malloc((strlen(allele) + 1) * sizeof(char));
		strcpy(allele_copy, allele);
		stList_append(dest->alleles, allele_copy);
	}
	dest->longest_allele_len = src->longest_allele_len;
        return dest;
}


// make a copy of a list of variants
stList* ptVariant_stList_copy(stList* variants){
	stList* variants_copy = stList_construct3(0, ptVariant_destruct);
	for(int i=0; i < stList_length(variants); i++){
		stList_append(variants_copy, ptVariant_copy(stList_get(variants, i)));
	}
	return variants_copy;
}


// Returns true if two variants are at the same location otherwise false
bool ptVariant_is_equal(ptVariant* var1, ptVariant* var2){
        return (var1->pos == var2->pos) && (strcmp(var1->contig, var2->contig) == 0);
}

//Returns true if variant exists in the list otherwise false
bool ptVariant_exist_in_list(ptVariant* var, stList* var_list){
	for (int i=0; i < stList_length(var_list); i++){
		if (ptVariant_is_equal(var, stList_get(var_list,i))){
			return true;
		}
	}
	return false;
}

void ptVariant_print(ptVariant* variant){
	printf("\n## Variant##\n");
	switch(variant->type){
		case VCF_REF:
			printf("REF\n");
			break;
		case VCF_SNP:
			printf("SNP\n");
			break;
		case VCF_DEL:
                        printf("DEL\n");
                        break;
                case VCF_INS:
                        printf("INS\n");
                        break;
		case VCF_INDEL:
			printf("INDEL\n");
			break;
		case VCF_ANY:
			printf("ANY\n");
			break;
		default:
			printf("OTHER\n");
			break;
	}
	printf("contig\tpos\tvaf\tgq\tps\n");
	printf("%s\t%d\t%.2f\t%d\t%d\n", variant->contig, variant->pos, variant->vaf, variant->gq, variant->ps);
	for(int i=0; i < variant->gt_len; i++){
		printf("Genotype %d:%d %s\n", i, variant->gt[i], stList_get(variant->alleles, variant->gt[i]));
		//printf("Genotype %d\n",variant->gt[i]);
	}
}

/// Swap the genotypes of the variants in a phase block
// It is assumed that the variants are sorted by phase block
/**
   @param variants	A list of variant records
   @param start		The index of the first variant in the phase block
   @param end		The index of the last variant in the phase block
*/
void ptVariant_swap_gt(stList* variants, int start, int end){
	for(int i = start; i <= end; i++){
		ptVariant* variant = stList_get(variants, i);
		int8_t gt1 = variant->gt[1];
		variant->gt[1] = variant->gt[0];
		variant->gt[0] = gt1;
	}
}


/// Parse variants from a vcf file and save the phased variants
// in a list. Each variant is saved in a ptVariant structure.
// If consistent_gt is True this function will check if the 
// first genotypes of all the phased variants in a phase block 
// are more consistent with the reference allele than the second
// genotypes.If it was not consistent then all the first and 
// second genotypes will be swapped.
// It can be useful for polishing purposes when we assume that 
// the first genotypes show the changes that has to be made on 
// the assembly sequence.
// Only the variants with two genotypes (diploid genome) are included.
/**
   @param vcf_path	The path to a vcf file
*/
stList* read_phased_variants(char* vcf_path, bool consistent_gt){

	stList* variants = stList_construct3(0,  ptVariant_destruct);

	//open vcf file
	htsFile *fp    = hts_open(vcf_path,"rb");
	
	//read header
	bcf_hdr_t *hdr = bcf_hdr_read(fp);	
	bcf1_t *rec    = bcf_init();
	

	int64_t ps_pre = -1; // The id of the previous phase block
	int32_t n_ref_ps_variants = 0; // Number of variants with their first genotype as ref in the current phase block
	int32_t n_ps_variants = 0; // Number of total variants in the current phase block
	int32_t n_variants = 0;


	// Each phase block id can take variable number of bytes
	// in the vcf so we are defining ps variables with different
	// sizes for casting properly
	int8_t* ps_int8;
	int16_t* ps_int16;
	int32_t* ps_int32;
	int64_t* ps_int64;
	int64_t ps;

	printf("Reading phased variants\n\n##Phase_Block_ID\tN_Ref\tN_All\tRef_Ratio\tSwapped\n");
	//save each vcf record
	while ( bcf_read(fp, hdr, rec)>=0 )
	{
		//unpack for read REF,ALT,INFO,etc
		bcf_unpack(rec, BCF_UN_STR);
		bcf_unpack(rec, BCF_UN_FMT);

		int32_t *gt_arr = NULL, ngt_arr = 0;
		int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);

		// get genotypes
                bcf_fmt_t * fmt_gt = bcf_get_fmt(hdr, rec, "GT");

		// skip if the number of genotypes is not equal to 2
		if(fmt_gt->n != 2) continue;
                int8_t gt1 = bcf_gt_allele(fmt_gt->p[0]);
		int8_t gt2 = bcf_gt_allele(fmt_gt->p[1]);

		// check if the variant is phased or not
                // an issue was opened here
                // https://github.com/samtools/htslib/issues/1113
                // saying that the 1st genotype never has the phase bit
                // so we look at the 2nd genotype to check if it's phased or not
                int is_phased = bcf_gt_is_phased(fmt_gt->p[1]);

		// skip unphased variants
                if (! is_phased) continue;

		// get contig
		char* contig = bcf_hdr_id2name(hdr, rec->rid);
		
		// get position
		int32_t pos = rec->pos;

		// get allele frequency
		bcf_fmt_t * fmt_vaf = bcf_get_fmt(hdr, rec, "VAF");
		float* vaf_float_ptr = (float*) fmt_vaf->p;
		float vaf = vaf_float_ptr[0];

		// get genotype quality
		bcf_fmt_t * fmt_gq = bcf_get_fmt(hdr, rec, "GQ");
		assert(fmt_gq->type == BCF_BT_INT8);
		int8_t gq = fmt_gq->p[0];

		// get phase block id
                bcf_fmt_t * fmt_ps = bcf_get_fmt(hdr, rec, "PS");
		switch (fmt_ps->type){
			case BCF_BT_INT8:
                                ps_int8 = (int8_t*) fmt_ps->p;
                                ps = ps_int8[0];
                                break;
			case BCF_BT_INT16:
				ps_int16 = (int16_t*) fmt_ps->p;
				ps = ps_int16[0];
				break;
			case BCF_BT_INT32:
				ps_int32 = (int32_t*) fmt_ps->p;
				ps = ps_int32[0];
				break;
			case BCF_BT_INT64:
				ps_int64 = (int64_t*) fmt_ps->p;
				ps = ps_int64[0];
                                break;
		}

		// get var type
		int8_t var_type = bcf_get_variant_types(rec);

		if (ps_pre != ps && stList_length(variants) > 0){ // if phase block changed
			printf("%d\t%d\t%d\t%.2f", ps_pre, n_ref_ps_variants, n_ps_variants, (double)n_ref_ps_variants/ n_ps_variants);
			// if less than half of the variants in the previous 
			// phase block (ps_pre) have reference allele as their first 
			// genotype then iterate over the variants in that
			// phase block and swap the genotypes for each variant
			if (n_ref_ps_variants < (0.5 * n_ps_variants) && consistent_gt){
				n_variants = stList_length(variants);
				//Swap genotypes
				ptVariant_swap_gt(variants, n_variants - n_ps_variants, n_variants - 1);
				printf("\tYES\n");
			}else{
				printf("\tNO\n");
			}
			// reset variables related to phase block
			n_ps_variants = 0;
			n_ref_ps_variants = 0;
		}
		// Make a ptVariant structure
		ptVariant* variant = ptVariant_construct(contig, pos, var_type, vaf, gq, ps);

		// add alleles to the variant structure
		// rec->d.allele[0] is REF
                for (int i=0; i<rec->n_allele; ++i){
                        ptVariant_append_allele(variant, rec->d.allele[i]);
		}

		// add genotypes
		ptVariant_append_gt(variant, gt1);
		ptVariant_append_gt(variant, gt2);

		// add variant to the list
		stList_append(variants, variant);

		// update variables related to phase block
		ps_pre = ps;
		n_ps_variants += 1;
		n_ref_ps_variants += gt1 == 0 ? 1 : 0;
	}

	// check the last phase block
	if (n_ref_ps_variants < (0.5 * n_ps_variants) && consistent_gt){
		printf("%d\t%d\t%d\t%.2f", ps_pre, n_ref_ps_variants, n_ps_variants, (double)n_ref_ps_variants/ n_ps_variants);
		n_variants = stList_length(variants);
		//Swap genotypes
		ptVariant_swap_gt(variants, n_variants - n_ps_variants, n_variants - 1);
		printf("\tYES\n");
	}else{
		printf("\tNO\n");
	}

	// Free variant record and vcf header
	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	int ret;
	if ( (ret=hts_close(fp)) )
	{
		fprintf(stderr,"hts_close(%s): non-zero status %d\n",vcf_path,ret);
		exit(ret);
	}
	return variants;
}

// 
// Make a new list of variants after removing the variants
// whose first genotypes are 0, which means that their first 
// genotypes are same as reference
//
/*
 * @param variants	stList of variants (can be the output of read_phased_variants())
 *
 */
stList* filter_ref_variants(stList* variants){
	stList* selected_variants = stList_construct3(0, ptVariant_destruct);
	for(int i=0; i < stList_length(variants); i++){
		ptVariant* variant = stList_get(variants, i);
		if (variant->gt[0] == 0) continue;
		ptVariant* variant_copy = ptVariant_copy(variant);
		stList_append(selected_variants, variant_copy);
	}
	return selected_variants;
}


// 
// Make a table of blocks that surround the given variants.
// The keys of the table are contig names and the values
// are related lists of variant blocks.
// 
// This function first makes a symmetric window around each variant
// The length of window is (2 * min_margin + 1) initially
// If min_margin is smaller than the length of variant then the
// length of the window will increase to (2 * longest_allele_len + 1)
// This is to make sure that we take enough flanking sequence to
// have a reliable alignment between reference and read sequence
// (which will be performed later by edlib functions)
// 
// If the blocks of close variants have overlap they will be merged
// into one block. 
// The data attribute of each block contains a stList of variants 
// that it spans over. This information is neccessary for polishing the
// reference sequence later.
//
// Note that variants should be sorted by contig name and then by position
// ptVariant_cmp and stList_sort can be used for doing so
//
/*
 * @param variants      	stList of variants (can be the output of filter_ref_variants())
 * @param fai			Fasta index (needed for not exceeding the length of contigs)
 * @param min_margin		Minimum margin length from each side of a variant
 * @return variant_blocks	Table of variant blocks (keys=contig names)
 */
stHash* ptBlock_extract_variant_blocks(stList* variants, const faidx_t* fai, int min_margin){
	stHash* variant_blocks = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, (void (*)(void *)) stList_destruct);
	ptVariant* variant = NULL;
	ptVariant* pre_variant = NULL;
	stList* blocks;
	int rfs=-1;
	int rfe=-1;
	int contig_len;
	ptBlock* block = NULL;
	for(int i=0; i < stList_length(variants); i++){
		variant = stList_get(variants, i);
		//ptVariant_print(variant);
		// if this is the first variant 
		// or the contig has changed
		if (i == 0 || strcmp(variant->contig, pre_variant->contig) != 0){
			// make a list of blocks for the new contig
			blocks = stList_construct3(0, ptBlock_destruct);
			char* contig = malloc(50);
			strcpy(contig, variant->contig);
			contig_len = faidx_seq_len64(fai, contig);
			// add the blocks list to the hash table with the key equal to contig name
			stHash_insert(variant_blocks, contig, blocks);
		}
		char* allele_ref = stList_get(variant->alleles, 0);
		char* allele_alt = stList_get(variant->alleles, variant->gt[0]);
		// Adjust margin if variant is longer than min_margin
		int margin = max(min_margin, variant->longest_allele_len);
		if (strlen(allele_ref) <= strlen(allele_alt)){ // insertion or snp
			rfs = variant->pos - margin;
			rfe = variant->pos + margin;
		}
		else if (strlen(allele_alt) < strlen(allele_ref)){ // deletion
			rfs = variant->pos - margin;
			rfe = variant->pos + strlen(allele_ref) - strlen(allele_alt) + margin;
		}
		// correct rfs and rfe if exceeded feasible coordinates
		rfs = rfs < 0 ? 0 : rfs;
		rfe = contig_len < rfe ? contig_len - 1 : rfe;
		// make a ptBlock struct
                block = ptBlock_construct(rfs, rfe,
                                          -1, -1, // sqs and sqe
                                          -1, -1); // rfd_s and rfd_e
		// make a list of variants for the new block
		stList* vars = stList_construct3(0, ptVariant_destruct);
                stList_append(vars, ptVariant_copy(variant));
                ptBlock_add_data(block, vars, stList_destruct); // add variant records as the block data
		int i_pre=stList_length(blocks) - 1;
		// This loop is for expanding the currect block
		// if it had overlap with the previous blocks.
		// Once this loop reaches a previous block
		// with no overlap we can be sure that it
		// does not exist any overlapping blocks before
		// that too so the loop will break there.
		// After the loop finished we remove all the overlapping blocks
		// and add the new block that spans all of them.
		for(; 0 <= i_pre; i_pre--){
			ptBlock* pre_block = stList_get(blocks, i_pre);
			stList* pre_vars = (stList*) pre_block->data;
			ptVariant* pre_var0 = stList_get(pre_vars, 0);
			// check for overlap with previous block
            		if (pre_block != NULL &&
                            strcmp(variant->contig, pre_var0->contig) == 0 &&
                            pre_block->rfs <= block->rfe &&
                            block->rfs <= pre_block->rfe){
				block->rfs = min(pre_block->rfs, block->rfs); // expand the start coordinate
				block->rfe = max(pre_block->rfe, block->rfe); // expand the end coordinate
				assert(pre_block->data != NULL);
				// copy variants to the new expanded block
				for(int j=0; j < stList_length(pre_vars); j++){
					ptVariant* var = stList_get(pre_vars, j);
					stList_append((stList*) block->data, ptVariant_copy(var)); // add variant record
				}
			}
			else{
				break;
			}
		}
		/*if (stList_length(blocks) > i_pre + 1){
			printf("XXXXXXXXXXXXXXXx\n");
			stList* x=(stList*) block->data;
			for(int j=0; j < stList_length(x); j++){
				ptVariant_print(stList_get(x,j));
			}
		}*/
		// remove previous blocks that had overlap
		stList_removeInterval(blocks, i_pre + 1, stList_length(blocks) - i_pre - 1);
		stList_append(blocks, block);
		pre_variant = variant;
	}

	// sort variants in each block
	stHashIterator* it = stHash_getIterator(variant_blocks);
        char* key;
        while((key = stHash_getNext(it)) != NULL){
                blocks = stHash_search(variant_blocks, key);
		for(int i=0; i<stList_length(blocks);i++){
			block = stList_get(blocks, i);
			//printf("##%d %d\n", i, stList_length(block->data));
			stList_sort((stList*)block->data, ptVariant_cmp);
		}
	}
	return variant_blocks;
}


//
// Fetch the read sequence of a given block
// The attributes sqs and sqe should have been set in the block
//
/*
 * @param alignment	Alignment saved as a ptAlignment struct
 * @param block		Block saved as a ptBlock struct
 * @return block_seq	Read sequence of the given block
 */
char* fetch_read_seq(ptAlignment* alignment, ptBlock* block){
	uint8_t* seq = bam_get_seq(alignment->record);
	//printf("block->sqs - block->sqe= %d-%d\t len=%d\n", block->sqs , block->sqe, block->sqe - block->sqs + 1);
	assert(block->sqe < alignment->record->core.l_qseq);
	assert(0 <= block->sqs);
        int block_len = block->sqe - block->sqs + 1;
        //allocate read sequence
        char* block_seq = (uint8_t*) malloc(block_len + 1);
        for(int k=0; k < block_len; k++)
                block_seq[k] = seq_nt16_str[bam_seqi(seq, block->sqs + k)];
	block_seq[block_len] = '\0';
	return block_seq;
}

//
// Fetch the corrected reference sequence of a given block
// The variants in the block will be applied to the reference sequence.
// The attributes rfs and rfe should have been set in the block.
//
/*
 * @param fai			Fasta index
 * @param block         	Block saved as a ptBlock struct
 * @param conting_name		Contig name
 * @return seq_corrected	Corrected reference sequence of the block    
 */

char* fetch_corrected_ref_seq(const faidx_t* fai, ptBlock* block, char* contig_name){
	int offset_orig = 0;
	int offset_corrected = 0;
	stList* vars = (stList*) block->data;
	ptVariant* var = NULL;
	char* reg = malloc(200);
        memset(reg,'\0',200);
        sprintf(reg, "{%s}:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);
        int len;
	//printf("%s\n",reg);
        char* seq_orig = fai_fetch(fai, reg, &len);
	char* seq_corrected = malloc((strlen(seq_orig) * 2) * sizeof(char));
	int size = strlen(seq_orig) * 2;
	memset(seq_corrected,'\0',size);
	// copy the bases before the first variant
	//memcpy(seq_corrected, seq_orig, var->pos - block->rfs);
	//offset_orig += var->pos - block->rfs;
	//offset_corrected += var->pos - block->rfs;
	// apply variations and copy the bases between variants
	for(int i=0; i < stList_length(vars); i++){
		var = stList_get(vars, i);
		if (block->rfs + offset_orig < var->pos){
			if (size <= var->pos - block->rfs - offset_orig +  offset_corrected){
                        	seq_corrected = realloc(seq_corrected, size * 2);
                        	size *= 2;
                	}
			memcpy(seq_corrected + offset_corrected, seq_orig + offset_orig, var->pos - block->rfs - offset_orig);
			//printf("pos=%d, rfs=%d\n", var->pos , block->rfs);
        		offset_corrected += var->pos - block->rfs - offset_orig;
			offset_orig = var->pos - block->rfs;
			seq_corrected[offset_corrected] = '\0';
			//printf("orig=%d, corrected=%d\n", offset_orig, offset_corrected);
			//printf("@@%s\n", seq_corrected);
		}
		char* allele_ref = stList_get(var->alleles, 0);
                char* allele_alt = stList_get(var->alleles, var->gt[0]);
		if (size <= strlen(allele_alt) + strlen(seq_corrected) + offset_corrected){
			seq_corrected = realloc(seq_corrected, size * 2);
			size *= 2;
		}
		memcpy(seq_corrected + offset_corrected, allele_alt, strlen(allele_alt));
		offset_orig += strlen(allele_ref);
		offset_corrected += strlen(allele_alt);
		//seq_corrected[offset_corrected] = '\0';
                        //printf("@#%s\n", seq_corrected);
	}
	if (offset_orig < (block->rfe - block->rfs + 1)){ // last bases to be copied
		strcpy(seq_corrected + offset_corrected, seq_orig + offset_orig);
		offset_corrected += strlen(seq_orig + offset_orig);
	}
	seq_corrected[offset_corrected] = '\0';

	//printf("orig\t%s\n", seq_orig);
	//printf("corr\t%s\n", seq_corrected);
	return seq_corrected;
}


//
// Project blocks onto read coordinates
// Make a new list of blocks. For each block rds_f and rde_f attributes are set after 
// projecting the reference coordinates to read by iterating over cigar operations.
// The new blocks contain the variant information
// Note that the output blocks of this function may have overlaps on read coordinates.
// To merge the overlapping blocks merge_variant_read_blocks() can be called.
//
/*
 * @param alignment		Alignment saved as a ptAlignment struct
 * @param ref_variant_blocks	stList of blocks with rfs and rfe
 * @return read_blocks		stList of blocks with rds_f and rde_f
 */
// ref_variant_blocks should be sorted by rfs with no overlap
stList* project_blocks_to_read(ptAlignment* alignment, stList* ref_variant_blocks){
	bam1_t* b = alignment->record;
	// Find the block that starts after the given alignment
	// based on the reference coordinates
	int j = 0;
	ptBlock* ref_block = stList_get(ref_variant_blocks, j);
	while(ref_block->rfs < alignment->rfs){
		if (stList_length(ref_variant_blocks) - 1 == j){
			ref_block = NULL;
			break;
		}
		j += 1;
		ref_block = stList_get(ref_variant_blocks, j);
	}
	// If the alignment does not include any block
	if (ref_block && alignment->rfe < ref_block->rfe){
		ref_block = NULL;
	}
	// Initiate a list for keeping the blocks with read coordinates
	stList* read_blocks = stList_construct3(0, ptBlock_destruct);

	// If the alignment does not include any block return empty list
	if (ref_block == NULL){
		return read_blocks;
	}

	//printf("first block: %d\t%d\n", ref_block->rfs, ref_block->rfe);
	//printf("alignment interval: %d\t%d\n", alignment->rfs, alignment->rfe);
	ptBlock* read_block;
	int rds_f=-1; int rde_f=-1;
	int cigar_rds_f; int cigar_rde_f;
	// Iterate over cigar operations
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	// Since projection is from reference to read coordinates
	// there is no need to do anything about insertions
	// Insertions do not advance in reference coordinates
	int i = 0;
	while(ptCigarIt_next(cigar_it)){
		if (bam_is_rev(b)){
			cigar_rds_f = -1 * cigar_it->rde_f;
			cigar_rde_f = -1 * cigar_it->rds_f;
		}
		else {
			cigar_rds_f = cigar_it->rds_f;
			cigar_rde_f = cigar_it->rde_f;
		}
		if(cigar_it->op == BAM_CMATCH ||
                   cigar_it->op == BAM_CEQUAL ||
                   cigar_it->op == BAM_CDIFF ||
                   cigar_it->op == BAM_CDEL) {
			/*if(i == 0){
				printf("first cigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
				i=1;
			}*/
			//The loop below iterates the blocks' overlap with the current cigar operation
			//There may exist multiple blocks within the current op
			while(ref_block && ref_block->rfe <= cigar_it->rfe){
				//printf("block end cigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);

				// update start locations of the projected coordinates
				if(cigar_it->rfs <= ref_block->rfs){
					/*printf("rds_f updated\n");
					printf("\tcigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
					printf("\tblock: %d\t%d\n", ref_block->rfs, ref_block->rfe);*/
					rds_f = cigar_it->op == BAM_CDEL ? cigar_rds_f : cigar_rds_f + (ref_block->rfs - cigar_it->rfs);
					//printf("rds_f = %d\n",rds_f);
					//sqs = cigar_it->op == BAM_CDEL ? cigar_it->sqs : cigar_it->sqs + (ref_block->rfs - cigar_it->rfs);
                                }
				/*printf("rde_f updated\n");
                                printf("\tcigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
                                printf("\tblock: %d\t%d\n", ref_block->rfs, ref_block->rfe);*/
				// update end locations of the projected coordinates
				rde_f = cigar_it->op == BAM_CDEL ? cigar_rde_f : cigar_rds_f + (ref_block->rfe - cigar_it->rfs);
				//printf("rde_f = %d\n",rde_f);
				//sqe = cigar_it->op == BAM_CDEL ? cigar_it->sqe : cigar_it->sqs + (ref_block->rfe - cigar_it->rfs);
				
				// Construct block with read coordinates
				if(bam_is_rev(b)){ //positive strand
					read_block = ptBlock_construct(ref_block->rfs, 
								       ref_block->rfe, 
								        -1,
									-1,
									-1 * rde_f,
									-1 * rds_f);
				} else{ // negative strand
					read_block = ptBlock_construct(ref_block->rfs,
                                                                                ref_block->rfe,
                                                                                -1,
                                                                                -1,
                                                                                rds_f,
                                                                                rde_f);
                                                                            
				}
				
				// append new read block
				if (read_block->rds_f <= read_block->rde_f){
					// make a copy of all the variants in the ref block
					stList* variants_copy = ptVariant_stList_copy((stList*) ref_block->data);
					// add all copied variants to the read block
					ptBlock_add_data(read_block, variants_copy, stList_destruct);
					//printf("read block added: %d\t%d\n", read_block->rds_f, read_block->rde_f);
					// add read block
					stList_append(read_blocks, read_block);
				}
				else{ // delete block if the whole block is within a deletion
					ptBlock_destruct(read_block);
				}

                                // update block index; j
                                if (j < stList_length(ref_variant_blocks) - 1) {
					j++;
					ref_block = stList_get(ref_variant_blocks, j);
				}
				else if (j == stList_length(ref_variant_blocks) - 1){
					ref_block = NULL;
                                }
				if(ref_block == NULL) break;// no more block remaining so break the cigar iteration
			}//end while
			// if the start of the next block is within the current operation
			if(ref_block && 
			   cigar_it->rfs <= ref_block->rfs &&
			   ref_block->rfs <= cigar_it->rfe ){
				/*printf("rds_f updated\n");
				printf("\tcigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
				printf("\tblock: %d\t%d\n", ref_block->rfs, ref_block->rfe);*/
				rds_f = cigar_it->op == BAM_CDEL ? cigar_rds_f : cigar_rds_f + (ref_block->rfs - cigar_it->rfs);
				//printf("rds_f = %d\n",rds_f);
				//sqs = cigar_it->op == BAM_CDEL ? cigar_it->sqs : cigar_it->sqs + (ref_block->rfs - cigar_it->rfs);
                        }//end if
                }// end if
	}//end of iteration over cigar ops
	//printf("CIGAR DONE!\n\n\n");
        ptCigarIt_destruct(cigar_it);
	return read_blocks;
}

/**
 * Merge overlapping blocks in the given list of variant blocks
 * When two overlapping blocks are merged all the variants in the
 * two blocks will be added to the final merged block
   @param blocks     a list of blocks sorted by rds_f
   @return      A list of merged blocks
 */
stList* merge_variant_read_blocks(stList* blocks){
        stList* blocks_merged = stList_construct3(0, ptBlock_destruct);
        if (stList_length(blocks) == 0) return blocks_merged;
        ptBlock* b = NULL;
	ptBlock* b_merged = NULL;
        for (int i = 0; i < stList_length(blocks); i++){
		b = stList_get(blocks, i);
		//printf("%d\t%d\n", b->rds_f, b->rde_f);
		if(i == 0){ // Initiate b_merged for the first block
			b_merged = ptBlock_copy(b, ptVariant_stList_copy);
			continue;
		}
		if(b_merged->rde_f < b->rds_f){// no overlap with previous merged block
			//save the merged block
			stList_append(blocks_merged, b_merged);
			// Initiate a new merged block
			b_merged = ptBlock_copy(b, ptVariant_stList_copy);
		}
		else { //there is overlap
			b_merged->rde_f = b->rde_f; // extend end pos of the merged block
			// add new variants if there is any
			stList* new_vars = (stList*) b->data;
			stList* curr_vars = (stList*) b_merged->data;
			for(int j = 0; j < stList_length(new_vars); j++){
				ptVariant* new_var = stList_get(new_vars, j);
				// if the new variant does not exist make a copy and add
				// to the variants of the merged block
				if (! ptVariant_exist_in_list(new_var, curr_vars)){
					stList_append(curr_vars, ptVariant_copy(new_var));
				}
			}
		}
        }

	/*printf("##########################\n");
	for (int i = 0; i < stList_length(blocks_merged); i++){
		b = stList_get(blocks_merged, i);
                printf("%d\t%d\n", b->rde_f, b->rds_f);
	}*/

	// Add the last merged block
	if(b_merged){
		stList_append(blocks_merged, b_merged);
	}

	/*printf("##########################\n");
	printf("stList_length(blocks_merged)=%d\n",stList_length(blocks_merged));
        for (int i = 0; i < stList_length(blocks_merged); i++){
                b = stList_get(blocks_merged, i);
                printf("%d\t%d\n", b->rde_f, b->rds_f);
		for(int j=0; j <stList_length((stList*) b->data); j++){
                        ptVariant_print(stList_get((stList*) b->data, j));
                }
        }*/

        return blocks_merged;
}

//
// Extract read sequence and corrected reference sequence 
// for the given block and then compute the edit distance 
// between the two sequences.
//
/*
 * @param alignment     Alignment saved as a ptAlignment struct
 * @param fai		Fasta index
 * @param block         Block saved as a ptBlock struct
 * @param contig_name	Contig name
 * @return edit_distance	Edit distance between read and ref
 */
int get_edit_distance(ptAlignment* alignment, faidx_t* fai, ptBlock* block, char* contig_name){
	char* read_seq = fetch_read_seq(alignment, block);
	char* corrected_ref_seq = fetch_corrected_ref_seq(fai, block, contig_name);
	printf("read\t%s\nref\t%s\n",read_seq, corrected_ref_seq);
	int edit_distance = 0;
	EdlibAlignResult result = edlibAlign(read_seq, strlen(read_seq), corrected_ref_seq, strlen(corrected_ref_seq), edlibDefaultAlignConfig());
    	if (!(result.status == EDLIB_STATUS_OK)) {
		fprintf(stderr, "Edlib didn't work!\n");
        	exit(EXIT_FAILURE);
    	}
	edit_distance = result.editDistance;
	printf("edit=%d\n\n",edit_distance);
    	edlibFreeAlignResult(result);
	free(corrected_ref_seq);
	free(read_seq);
	return edit_distance;
}


//
// Project read blocks to reference by iterating over the alignment 
// cigar operations.
// Extract read sequence and corrected reference sequence for each 
// block and then compute the edit distance between the two sequences. 
// The output of this function would be the summation of edit distances 
// between read and reference for all blocks.
// Note that the input blocks should be in read coordinates (recommended to be 
// the output of merge_variant_read_blocks())
//
/*
 * @param alignment     	Alignment saved as a ptAlignment struct
 * @param fai           	Fasta index
 * @param variant_read_blocks	stList of blocks each of which saved as a ptBlock struct
 * @param contig_name   	Contig name
 * @return edit_distance	Summation of edit distances between read and reference for all blocks
 */
int get_edit_distance_all_blocks(ptAlignment* alignment, const faidx_t* fai, char* contig_name, stList* variant_read_blocks){
	bam1_t* b = alignment->record;
	int j = bam_is_rev(b) ? stList_length(variant_read_blocks) - 1 : 0;
	int edit_distance = 0;
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	int block_rds_f, block_rde_f, cigar_rds_f, cigar_rde_f, rfs, rfe, sqs, sqe;
	ptBlock* read_block = stList_get(variant_read_blocks, j);
	// reverse each block interval to make it work with the projecting algorithm below
	if (bam_is_rev(b)){
		block_rds_f = -1 * read_block->rde_f;
		block_rde_f = -1 * read_block->rds_f;
	}
	else {
		block_rds_f = read_block->rds_f;
		block_rde_f = read_block->rde_f;
	}
        //iterate over cigar operations from left to right (w.r.t reference)
	while(ptCigarIt_next(cigar_it)){
		// reverse each cigar operation interval to make it work with the projecting algorithm below
                if (bam_is_rev(b)){
                	cigar_rds_f = -1 * cigar_it->rde_f;
                        cigar_rde_f = -1 * cigar_it->rds_f;
		}
		else {
			cigar_rds_f = cigar_it->rds_f;
			cigar_rde_f = cigar_it->rde_f;
		}
		// match/ mismatch/ insertion
		// for M and I the main algorithm is the same except
		// some minor parts that are corrected by
		// the conditional statement "cigar_it->op == BAM_CINS ? : "
		if(cigar_it->op == BAM_CMATCH ||
		   cigar_it->op == BAM_CEQUAL ||
		   cigar_it->op == BAM_CDIFF ||
		   cigar_it->op == BAM_CINS) {
			//The loop below iterate the blocks overlap with the current cigar operation
			//There may exist multiple blocks within the current op
			while(read_block && block_rde_f <= cigar_rde_f){
				// update start locations of the projected coordinates
				if(cigar_rds_f <= block_rds_f){
					rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
					sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
				}
				// update end locations of the projected coordinates
				rfe = cigar_it->op == BAM_CINS ? cigar_it->rfe : cigar_it->rfs + (block_rde_f - cigar_rds_f);
				sqe = cigar_it->sqs + (block_rde_f - cigar_rds_f);
				/*if(1000 < (rfe-rfs)){
					printf("####%s\t%d\t%d\n",contig_name,rfs,rfe);
					printf("%d\t%d\n",read_block->rds_f, read_block->rde_f);
				}*/
				// get block
				ptBlock* block = ptBlock_construct(rfs,
						                   rfe,
							           sqs,
							           sqe,
							           read_block->rds_f,
								   read_block->rde_f);
				
				//printf("%d\t%d\n", rfs, rfe);
				// all_variants contains all variants from all regions in the genome 
				// that have overlap with alignment of this read block
				stList* all_variants = (stList*) read_block->data;
				//printf("All variants %d\n", stList_length(all_variants));
				//for(int k=0; k < stList_length(all_variants);k++){
                                //        ptVariant_print(stList_get(all_variants,k));
                               // }
				// kept_variants will contain only the variants from this contig
				// and located within the block
				stList* kept_variants = stList_construct3(0, ptVariant_destruct);
				for(int i=0; i < stList_length(all_variants); i++){
					ptVariant* variant = stList_get(all_variants, i);
					if(strcmp(variant->contig, contig_name) == 0 && 
					   block->rfs <= variant->pos &&
					   variant->pos <= block->rfe){
						stList_append(kept_variants, ptVariant_copy(variant));
					}
				}
				ptBlock_add_data(block, (void*) kept_variants, stList_destruct);
				//printf("Kept variants %d\n", stList_length(kept_variants));
				//for(int k=0; k < stList_length(kept_variants);k++){
				//	ptVariant_print(stList_get(kept_variants,k));
				//}
				edit_distance += get_edit_distance(alignment, fai, block, contig_name);
		                // update block index; j
				if (bam_is_rev(b) && j > 0){
					j--;
					read_block = stList_get(variant_read_blocks, j);
					block_rds_f = -1 * read_block->rde_f;
					block_rde_f = -1 * read_block->rds_f;
				}
				else if (!bam_is_rev(b) && j < stList_length(variant_read_blocks) - 1) {
					j++;
					read_block = stList_get(variant_read_blocks, j);
					block_rds_f = read_block->rds_f;
					block_rde_f = read_block->rde_f;
				}
				else if (j == 0 || j == stList_length(variant_read_blocks) - 1){
					read_block = NULL;
				}
			}
			if(read_block == NULL) break;// no more block remaining so break the cigar iteration
			// if the start of the next block is within the current operation
                        if(cigar_rds_f <= block_rds_f &&
                           block_rds_f <= cigar_rde_f ){
				rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
                                sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
			}
		}
	}	
	return edit_distance;
}



/*! @typedef
 * @abstract Structure for saving a marker's information.
 * @field alignment_idx		Index of the alignment (could be secondary or primary) from where this marker is taken
 * @field read_pos_f		Position of the marker in the read (All positions should be computed w.r.t
 * 				the forward strand for the inconsistency between negative and positive alignments)
 * @field base_idx		Index of the base in the query seq of the corresponding alignment. It may be different
 * 				from read_pos_f because of two reasons; negative orientation and hard clipping
 * @field base_q		Base quality of the marker (could be the raw quality or any processed form e.g BAQ)
 * @field is_match		True if the marker base is equal to the reference and False otherwise
 */
typedef struct {
	int32_t alignment_idx;
	int32_t read_pos_f;
	int32_t base_idx;
	int32_t base_q;
	bool is_match;
	int32_t ref_pos;
}ptMarker;


// Create a ptMarker struct
/**
 * @param alignment_idx		The index of the alignment (could be secondary or primary) from where this marker is taken
 * @param read_pos_f		The position of the marker in the read (All positions should be computed w.r.t to
 * 				the forward strand for the inconsistency between negative and positive alignments)
 * @param base_q    		The base quality of the marker (could be the raw quality or any processed form e.g BAQ)
 * @param is_match		True if the marker base is equal to the reference and False otherwise
 * @return 			Pointer to a new marker struct on success; NULL on failure
 * 				The ptMarker struct returned by a successful call should be freed
 * 				via free() when it is no longer needed.
 */
ptMarker* ptMarker_construct(int32_t alignment_idx, int32_t base_idx, int32_t read_pos_f, uint8_t base_q, bool is_match, int32_t ref_pos){
	ptMarker* marker = (ptMarker*) malloc(sizeof(ptMarker));
	marker->alignment_idx = alignment_idx;
        marker->is_match = is_match;
        marker->read_pos_f = read_pos_f;
	marker->base_idx = base_idx;
        marker->base_q = base_q;
	marker->ref_pos = ref_pos;
	return marker;
}


//
// Make a copy of a marker saved as a ptMarker struct
//
/*
 * @param src		Marker to be copied
 * @return dest		Copied marker
 *
 */
ptMarker* ptMarker_copy(ptMarker* src){
        ptMarker* dest = (ptMarker*) malloc(sizeof(ptMarker));
        dest->alignment_idx = src->alignment_idx;
        dest->is_match = src->is_match;
        dest->read_pos_f = src->read_pos_f;
        dest->base_idx = src->base_idx;
        dest->base_q = src->base_q;
	dest->ref_pos = src->ref_pos;
        return dest;
}

//
// Compare two ptMarker structs.
// The comparison is first based on their positions in the read and if both given markers
// are located at the same position then it looks at alignment_idx.
//
/*
 * @param a	Pointer to the first ptMarker structure
 * @param b	Pointer to the second ptMarker structure
 * @return 	Nagative when a should be before b and positive otherwise
 *
 */
int ptMarker_cmp(const void *a, const void *b){
        ptMarker* m1 = (ptMarker*) a;
        ptMarker* m2 = (ptMarker*) b;
	if (m1->read_pos_f == m2->read_pos_f) {
		return m1->alignment_idx - m2->alignment_idx;
	}
	else {
		return m1->read_pos_f - m2->read_pos_f;
	}
}

//
// This function will create a marker with its attribute is_match set to true
// It iterates over the cigar string to find the related seq and ref coordinates too
// and put this information into the created marker
//
/*
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignment_idx		Index of the alignment that should be used for making a match marker
 * @param read_pos_f		Read position (not seq position) for which a match marker should be created
 * @return match		Match marker created based on the given coordinate
 *
 */
ptMarker* ptMarker_construct_match(ptAlignment** alignments, int32_t alignment_idx, int32_t read_pos_f){
	bam1_t* record = alignments[alignment_idx]->record;
        uint8_t* q = bam_get_qual(record);
	ptMarker* match;
	uint32_t* cigar = bam_get_cigar(record);
	int lclip=0;
	int rclip=0;
	if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
                lclip = bam_cigar_oplen(cigar[0]);
        }
        if (bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP) {
                rclip = bam_cigar_oplen(cigar[record->core.n_cigar - 1]);
        }
	// construct the match marker
	if (bam_is_rev(record)){
		match = ptMarker_construct(alignment_idx,
					   record->core.l_qseq + rclip - read_pos_f - 1,
                                           read_pos_f,
                                           q[record->core.l_qseq + rclip - read_pos_f - 1],
                                           true,
					   -1);
        }
        else {
		match = ptMarker_construct(alignment_idx,
				           read_pos_f - lclip,
                                           read_pos_f,
                                           q[read_pos_f - lclip],
                                           true,
					   -1);
        }
	return match;
}

//
// If there exists a marker with a quality lower than the given threshold
// this function will remove all the markers in the same read position (read_pos_f)
//
// Note that a single read position may have as many markers as up to the number of alignments
// If BAQ is calculated for read bases a single position may obtain different quality 
// values for different alignments. This function removes markers from the same position 
// if in at least one alignment it obtained a BAQ lower than the threshold.
//
/*
 * @param markers_p             Pointer to a stList of markers each of which saved as ptMarker struct
 * @param threshold		Minimum quality value for a marker
 *
 */
void filter_lowq_markers(stList** markers_p, int threshold){
	stList* markers = *markers_p;
	int markers_len = stList_length(markers);
	if(markers_len == 0) return;
        stList* keep_markers = stList_construct3(0, free);
	ptMarker* pre_marker =NULL;
	ptMarker* marker;
	bool keep_flag = false;
	int idx_s = 0;
	int min_q = 100;
	ptMarker* marker_copy;
	for(int j=0; j < markers_len; j++){
		marker = stList_get(markers, j);
		if (pre_marker && marker->read_pos_f != pre_marker->read_pos_f){
			if (min_q > threshold){
				for(int k=idx_s; k<j; k++){
					// copy all markers with the same read_pos_f and add them
					marker_copy = ptMarker_copy(stList_get(markers, k));
					marker_copy->base_q = min_q;
					stList_append(keep_markers, marker_copy);
				}
			}
			idx_s = j;
			// reset min_q for the next read position
			min_q = 100;
		}
		if (min_q > marker->base_q) min_q = marker->base_q;
		/*
		if(marker->base_q > threshold) {
			keep_flag |= true;
			DEBUG_PRINT("q = %d ( > %d) marker->read_pos_f = %d\n", marker->base_q, threshold, marker->read_pos_f);
		}*/
		pre_marker = marker;
	}
	if (min_q > threshold){
		for(int k=idx_s; k<markers_len; k++){
			marker_copy = ptMarker_copy(stList_get(markers, k));
                        marker_copy->base_q = min_q;
                        stList_append(keep_markers, marker_copy);
                }
        }
	stList_destruct(markers);
	*markers_p = keep_markers;
}


//
// If there exist a read position which is located in an insertion in at least
// one alignment this function removes all the related markers in that position
//
/*
 * @param markers_p             Pointer to a stList of markers each of which saved as ptMarker struct
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 *
 */
void filter_ins_markers(stList** markers_p, ptAlignment** alignments, int alignments_len){
	stList* markers = *markers_p;
	int markers_len = stList_length(markers);
	if(markers_len == 0) return;
	bool* markers_keep_flags = (bool*) malloc(markers_len * sizeof(bool));
	for(int i=0; i< markers_len; i++) markers_keep_flags[i] =true;
	bam1_t* b;
	ptMarker* marker;
	ptCigarIt* cigar_it;
	int j; // marker index
        for(int i=0; i < alignments_len; i++){
		b = alignments[i]->record;
		j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
		marker = stList_get(markers, j);
		cigar_it = ptCigarIt_construct(b, true, true);
		//iterate over cigars
		while(ptCigarIt_next(cigar_it)){
			//check if marker is located within the cigar operation
			while(marker &&
			      (cigar_it->rds_f <= marker->read_pos_f) && 
			      (cigar_it->rde_f >= marker->read_pos_f)){
				//check if the marker is located within an insertion
				if(cigar_it->op == BAM_CINS ||
				   cigar_it->op == BAM_CSOFT_CLIP ||
				   cigar_it->op == BAM_CHARD_CLIP){
					markers_keep_flags[j] = false;
				}
				//update reference positions in the meanwhile!
				if(cigar_it->op == BAM_CEQUAL && marker->alignment_idx == i){
                                	marker->ref_pos = bam_is_rev(b) ? cigar_it->rfs + cigar_it->rde_f - marker->read_pos_f : cigar_it->rfs + marker->read_pos_f - cigar_it->rds_f;
                                }
				j += bam_is_rev(b) ? -1 : 1;
				marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
			}	
		}//end of iteration over cigar ops
		ptCigarIt_destruct(cigar_it);
	}// end of iteration over alignments
	// make a list for saving the remaining markers
	stList* keep_markers = stList_construct3(0, free);
	for(j=0; j < markers_len; j++){
		if(markers_keep_flags[j]){
			marker = stList_get(markers, j);
			stList_append(keep_markers, ptMarker_copy(marker));
		}
	}
	stList_destruct(markers);
        *markers_p = keep_markers;
}


//
// If there exist a read position which is mismatch in all alignments
// this function will remove those related markers since it is probably
// a read error. 
//
/*
 * @param markers_p             Pointer to a stList of markers each of which saved as ptMarker struct
 * @param alignments_len        Length of alignments array
 *
 */

void remove_all_mismatch_markers(stList** markers_p, int alignments_len){
	stList* markers = *markers_p;
	stList_sort(markers, ptMarker_cmp);
	int markers_len = stList_length(markers);
	if(markers_len == 0) return;
        bool* markers_keep_flags = (bool*) malloc(markers_len * sizeof(bool));
        memset(markers_keep_flags, 1, markers_len);
	ptMarker* marker;
	ptMarker* pre_marker = NULL;
	int occ = 0;
	for(int i=0; i < markers_len; i++){
		marker = stList_get(markers, i);
		if (pre_marker  &&
                    pre_marker->read_pos_f < marker->read_pos_f){
			if(occ == alignments_len){// if all markers are mismatch
				memset(markers_keep_flags + i - alignments_len, 0, alignments_len);
			}
			occ = 0;
		}
		pre_marker = marker;
		occ += 1;
	}
	if(occ == alignments_len){
		memset(markers_keep_flags + markers_len - alignments_len, 0, alignments_len);
        }
	// make a list for saving the remaining markers
        stList* keep_markers = stList_construct3(0, free);
        for(int j=0; j < markers_len; j++){
                if(markers_keep_flags[j]){
                        marker = stList_get(markers, j);
                        stList_append(keep_markers, ptMarker_copy(marker));
                }
        }
	// free previous list of markers
	stList_destruct(markers);
	// update the list
	*markers_p = keep_markers;
}

// 
// It is assumed that all markers in the given list are mismatch (is_match == false)
// A mismatch marker maybe a match in another alignment. 
// This function sorts markers first based on read_pos_f and then alignment index.
// It adds match markers in a way that after calling this function for each marker 
// position we will have exactly one match or mismatch marker per alignment. 
//
/*
 * @param markers_p       	Pointer to a stList of markers each of which saved as ptMarker struct
 *                      	All of these markers are assumed to be mismatch
 * @param alignments    	An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len	Length of alignments array
 *
 */
void sort_and_fill_markers(stList** markers_p, ptAlignment** alignments, int alignments_len){
	stList* markers = *markers_p;
	stList_sort(markers, ptMarker_cmp);
	int idx=0;
	ptMarker* match;
	ptMarker* marker;
	ptMarker* pre_marker = NULL;

	// markers_len is the number of initially given markers (all mismatch)
	// it does not get updated while adding match markers
	int markers_len = stList_length(markers);
	for(int64_t i=0; i < markers_len; i++){
		marker = stList_get(markers, i);
		// none of the initially given markers should be match
		assert(!marker->is_match);
		// if the previous location of read is passed
		// then we start adding matched markers for the
		// remaining alignments (if it existed any)
		if (pre_marker  &&
		    pre_marker->read_pos_f < marker->read_pos_f){
			for (int64_t j=idx; j < alignments_len; j++){
                        	match = ptMarker_construct_match(alignments, j, pre_marker->read_pos_f);
                        	stList_append(markers, match);
                	}
			idx = 0;
		}
		pre_marker = marker;
		// add match markers for the alignments whose indices start
		// from idx and end one before marker->alignment_idx
		for (int64_t j=idx; j < marker->alignment_idx; j++){
			match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
                	stList_append(markers, match);
        	}
		// update the first index of the next alignments that may lack this marker
		idx = marker->alignment_idx + 1;
	}
	if(markers_len > 0){
		// add the last matches if there still exist to be added
		for (int64_t j=idx; j < alignments_len; j++){
			match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
        		stList_append(markers, match);
        	}
	}
	stList_sort(markers, ptMarker_cmp); // sort again to locate the matches in their correct positions
}

// Quality value (q) is calculated by this formula as a function of probability (p)
// q(p) = -10 * log(p) 
// This function reveives q(p) and returns q(1-p)
// As edge cases:
//	 1. If q == 0 it reutrns  93
//	 2. If q >= 93 it returns 0
//
/*
 * @param q		The base quality q(p) (either raw or BAQ)
 * @return rev_q 	The reversed base quality q(1-p)
 */
double reverse_quality(uint8_t q){
	if (q >= 93) return 0;
	if (q == 0) return 93;
	double p = 1 - pow(10, (double)q / -10);
	double rev_q = -10 * log(p);
	return rev_q;
}


// Calculate the alignment scores based on the base qualities
// of the selected markers (could be raw base quality or BAQ)
// Each score is in logarithmic scale and represents the probability
// of all markers being match in the corresponding alignment. 
// We will have two cases:
//
// 1. If a marker (base) is match :
// 		The probability of being a real match is equal to 1 - P(error)
// 		Q = -10 * log(P(error)) is saved as base quality so 
// 		by calling reverse_quality we can calculate 10 * log(1 - P(error)) (with a negative sign) 
// 		is which a quality value of being a real match. It's usually close to 0.
// 2. If a marker (base) is mismatch:
//		The probability of being a real match is equal to P(error) / 3
// 		assuming that it is evenly probable to be a base other than what
// 		is now in the read sequence and one of those bases can be match to the
// 		reference sequence. Therefore the related quality value can be 
// 		calculated by 10 * log(P(error)) - 10 * log(3) = -Q - 10 * log(3)
//
// Note that the quality value does not have a negative sign in the formula so it is 
// always a negative value. The total score of each alignment is the summation of the 
// quality values of markers. One marker may obtain different quality values for 
// different alignments both because of different BAQ values and being match or 
// mismatch in different alignments.An alignment with higher score is more probable to 
// be to the correct haplotype.
//
/*
 * @param markers	stList of markers each of which saved as ptMarker struct
 * @param alignments    An array of alignments each of which saved as ptAlignment struct
 */
void calc_alignment_score(stList* markers, ptAlignment** alignments){
        ptMarker* marker;
        for(int64_t j=0; j < stList_length(markers); j++){
                marker =  stList_get(markers, j);
		if (marker->is_match){

			// p(all matches are real match) = (1-e_1) * (1-e_2) * ... (1-e_n1)
			// q(all matches are real match) = -rev(q_1) + -rev(q_2) + ... + -rev(q_n1)
			alignments[marker->alignment_idx]->score -= reverse_quality(marker->base_q) ;
		}
		else{
			// p(all mismtaches are real match) = (e1/3) * (e2/3) * ... (en/3)
                        // q(all mismatches are real match) = (-q1 - 10log(3)) + (-q2 - 10log(3)) + ... + (-qn - 10log(3))
			alignments[marker->alignment_idx]->score -= -1 * marker->base_q - 10 * log(3);
		}
		// p(all markers are real match for this alignment) = p(all matches are real match) * p(all mismtaches are real match)
		// p(all markers are real match for this alignment) = p(this alignment is correct)
		// final q = q(all matches are real match) + q(all mismatches are real match)
	}
}

// After setting alignment score values by calling calc_alignment_score
// this function can be called to get the index of the best alignment
// with highest score
//
/*
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 * @param prim_margin		The score of the best alignment should be higher than the 
 * 				score of the current primary alignment by this value
 * @param min_score		The minimum score of the best alignment. If the best alignment
 * 				had a score lower than this value the primary alignment will be
 * 				returned as the best alignment
 * @param prim_margin_random	If the scores of the best alignment and current primary alignment
 * 				was closer than this value then one of them will be returned randomly 
 *
 */

int get_best_record(ptAlignment** alignments, int alignments_len, double prim_margin, double min_score, double prim_margin_random){
	assert(alignments_len > 0);
	if (alignments_len == 1) return 0;
	double max_score = -DBL_MAX;
	int max_idx = -1;
	double prim_score = -DBL_MAX;
	int prim_idx = -1;
	for(int i=0; i < alignments_len; i++){
		if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
			prim_idx = i;
			prim_score = alignments[i]->score;
		}
		else if (max_score < alignments[i]->score){
			max_idx = i;
			max_score = alignments[i]->score;
		}
	}
	int rnd = rand() % 2; // 50% chance for rnd=0 (same for rnd=1)
	if (abs(max_score - prim_score) < prim_margin_random){
		return rnd == 0 ? prim_idx : max_idx; // if the scores were closer than prim_margin_random return one randomly
	}
	if (prim_idx == -1 || max_score < (prim_score + prim_margin) || max_score < min_score) return prim_idx;
	else return max_idx;
}


// Find the confident blocks for the given alignment
// Confident blocks are maximal blocks with no insertion longer
// than the given threshold
/*
 * @param alignment	Alignment record saved as a ptAlignment struct
 * @param threshold	Maximum length of an insertion that may exist in a confident block
 */
stList* find_confident_blocks(ptAlignment* alignment, int threshold){
	bam1_t* b = alignment->record;
	stList* conf_blocks = stList_construct3(0, ptBlock_destruct);
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	int conf_sqs = 0; // the seq start of the confident block
	int conf_rfs = b->core.pos; // the ref start of the confident block
	int conf_rd_f = bam_is_rev(b) ? cigar_it->rde_f : cigar_it->rds_f;
	ptBlock* block = NULL;
	while(ptCigarIt_next(cigar_it)){
		switch(cigar_it->op) {
                	case BAM_CINS:
			case BAM_CDEL:
                        	if(cigar_it->len > threshold && 
				   conf_sqs < cigar_it->sqs &&
				   conf_rfs < cigar_it->rfs){
					if (bam_is_rev(b)) {
						block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
								          conf_sqs, cigar_it->sqs - 1,
									  cigar_it->rde_f + 1, conf_rd_f);
					}
					else {
						block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
								          conf_sqs, cigar_it->sqs - 1,
									  conf_rd_f, cigar_it->rds_f - 1);
					}
					stList_append(conf_blocks, block);
				}
				if(cigar_it->len > threshold){
					conf_sqs = cigar_it->sqe + 1;
                                	conf_rfs = cigar_it->rfe + 1;
					conf_rd_f = bam_is_rev(b) ? cigar_it->rds_f - 1 : cigar_it->rde_f + 1;
				}
                        	break;
                	case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
                        	if(conf_sqs < cigar_it->sqs &&
				   conf_rfs < cigar_it->rfs){
					if (bam_is_rev(b)) {
                                        	block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                                  	  conf_sqs, cigar_it->sqs - 1,
									  cigar_it->rde_f + 1, conf_rd_f);
					}
					else {
						block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                                          conf_sqs, cigar_it->sqs - 1,
                                                                          conf_rd_f, cigar_it->rds_f - 1);
					}
					stList_append(conf_blocks, block);
				}
                                conf_sqs = cigar_it->sqe + 1;
                                conf_rfs = cigar_it->rfe + 1;
				conf_rd_f = bam_is_rev(b) ? cigar_it->rds_f - 1 : cigar_it->rde_f + 1;
                        	break;
        	}

	}
	// when the last block is not terminated with SOFT or HARD clipping
	if(conf_sqs <= cigar_it->sqe){
		if (bam_is_rev(b)) {
			block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                                  conf_sqs, cigar_it->sqe,
						  cigar_it->rds_f, conf_rd_f);
		}
		else {
			block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                                  conf_sqs, cigar_it->sqe,
						  conf_rd_f, cigar_it->rde_f);
		}
		stList_append(conf_blocks, block);
	}
	ptCigarIt_destruct(cigar_it);
	return conf_blocks;
}

/// Get the intersection of two arrays of blocks
/**
 * @param blocks1     1st set of blocks
 * @param blocks2     2nd set of blocks
 *
 * @return      stList of intersected blocks
 *
 * @note        Blocks should be sorted by rds_f and no overlap 
 * 		should exist between the blocks of each set
 */
stList* intersect_by_rd_f(stList* blocks1, stList* blocks2){
	stList* blocks_intersect = stList_construct3(0, ptBlock_destruct);
	if (stList_length(blocks1) == 0 || stList_length(blocks2) == 0) return blocks_intersect;
	ptBlock* b1;
	ptBlock* b2;
	ptBlock* b;
	int j = 0; 
	b2 = stList_get(blocks2, j);
	for (int i = 0; i < stList_length(blocks1); i++){
		b1 = stList_get(blocks1, i);
		// Find the first block (b2) in the 2nd set whose end point is after 
		// the start point of the current block (b1) in the 1st set
                while (b2 && b2->rde_f < b1->rds_f) {
			j++;
                        b2 = j == stList_length(blocks2) ? NULL : stList_get(blocks2, j);
                }
		// Because of the previous loop, b2's end point is now after b1's
		// start point. If b2's start point is also before the b1's end point
		// it means there exist an overlap so we save that overlap in 
		// blocks_intersect and go to the next block (update b2)
                while (b2 && b2->rds_f < b1->rde_f){
			b = ptBlock_construct(-1, 
					      -1, 
					      -1, 
					      -1, 
					      max(b1->rds_f, b2->rds_f), 
					      min(b1->rde_f, b2->rde_f));
                        stList_append(blocks_intersect, b);
			// Go to next b2 if current b2's end point is
			// before b1's end point
                        if (b2->rde_f <= b1->rde_f) {
                                j++;
                        	b2 = j == stList_length(blocks2) ? NULL : stList_get(blocks2, j);
                        }
                        else break;
                }
        }
	return blocks_intersect;
}


/// Given a list of alignments for a single read 
/// find the confident blocks for each alignment 
/// and update the corresponding attribute (alignment->conf_blocks)
/// for each alignment
/**
 * @param alignments		An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len	Length of alignments array
 * @param threshoold		Maximum length of an insertion in each confident block
 */
void set_confident_blocks(ptAlignment** alignments, int alignments_len, int threshold){
	for(int i=0; i < alignments_len; i++){
		alignments[i]->conf_blocks = find_confident_blocks(alignments[i], threshold);
	}
}

/// Given an alignment and the list of markers
/// find the blocks flanking the markers
/// and update the corresponding attribute (alignment->flank_blocks)
/**
 * @param alignments		An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len	Length of alignments array
 * @param margin		The flanking block around each marker
 * 				will extend from each side by this margin
 *
 * @note 			The markers should be sorted by read_pos_f
 */
stList* find_flanking_blocks(ptAlignment* alignment, stList* markers, int margin){
	bam1_t* b = alignment->record;
        stList* flank_blocks = stList_construct3(0, ptBlock_destruct);
        ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	// get the first merker for initializing start and end
	int i = 0;
	ptMarker* marker = stList_get(markers, i);
	int start = max(alignment->rds_f, marker->read_pos_f - margin);
	int end = min(alignment->rde_f, marker->read_pos_f + margin);
	i++;
	int curr_start;
	int curr_end;
	ptBlock* block = NULL;
        while(i < stList_length(markers)){
		marker = stList_get(markers, i);
		curr_start = max(alignment->rds_f, marker->read_pos_f - margin);
		curr_end =  min(alignment->rde_f, marker->read_pos_f + margin);
		if(curr_start < end) {
			end = curr_end;
		}
		else {
			// add the previous flanking block
			block = ptBlock_construct(-1, -1,
                                                  -1, -1,
                                                  start, end);
                        stList_append(flank_blocks, block);
			//update start and end for the next flanking block
			start = curr_start;
			end = curr_end;
		}
		i++;
	}
	// add the last flanking block
	block = ptBlock_construct(-1, -1,
                                  -1, -1,
                                  start, end);
        stList_append(flank_blocks, block);
	ptCigarIt_destruct(cigar_it);
	return flank_blocks;
}

/// Given a list of alignments for a single read
/// find the flanking blocks within each alignment
/// and update the corresponding attribute (alignment->flank_blocks)
/**
 * @param alignments		An array of ptAlignment structure
 * @param alignments_len	Length of alignments array
 * @param margin		For each marker we will dedicate a flanking block that spans
 * 				a window of length (2 * margin + 1) with the marker in its center
 */
void set_flanking_blocks(ptAlignment** alignments, int alignments_len, stList* markers, int margin){
        for(int i=0; i < alignments_len; i++){
                alignments[i]->flank_blocks = find_flanking_blocks(alignments[i], markers, margin);
        }
}


/// Given a list of alignments for a single read, 
/// find the consensus confident blocks using all previously 
/// found confident blocks and flanking blocks. In other words 
//  it takes the intersection of all sets of confident and 
//  flanking blocks for a single read. There should exist one 
//  set of confident blocks and one set of flanking blocks per 
//  alignment. They can be set by calling the functions set_flanking_blocks
//  and set_confident_blocks 
/**
 * @param alignments		An array of ptAlignment structure
 * @param alignments_len	Length of alignments array
 * @param threshold		The indels strictly longer than this threshold will cut the confident blocks
 *
 * @return 			The number of output confident blocks
 *
 * @note 			The conf_blocks attribute of all the alignments will be updated.
 * 				So all alignments will have the same conf_blocks after 
 * 				calling this function.
 */
int correct_conf_blocks(ptAlignment** alignments, int alignments_len, int threshold){
	assert(alignments_len > 0);
	stList_sort(alignments[0]->conf_blocks, ptBlock_cmp_rds_f);
        stList* blocks = stList_copy(alignments[0]->conf_blocks, NULL);
	ptBlock* block;
        stList* blocks_new;
	//intersect confident blocks
        for(int i=1; i < alignments_len; i++){
		stList_sort(alignments[i]->conf_blocks, ptBlock_cmp_rds_f);
                blocks_new = intersect_by_rd_f(blocks, alignments[i]->conf_blocks);
                stList_destruct(blocks);
                blocks = blocks_new;
        }
	//intersect flanking blocks
	for(int i=0; i < alignments_len; i++){
                stList_sort(alignments[i]->flank_blocks, ptBlock_cmp_rds_f);
                blocks_new = intersect_by_rd_f(blocks, alignments[i]->flank_blocks);
                stList_destruct(blocks);
                blocks = blocks_new;
        }
	if(stList_length(blocks) == 0){
		DEBUG_PRINT("\t\t\t### No Consensus Blocks!\n");
		for(int i=0; i < alignments_len; i++){
			stList_destruct(alignments[i]->conf_blocks);
			alignments[i]->conf_blocks = stList_construct3(0, ptBlock_destruct);
		}
		return 0;
	}
	stList* corrected_conf_blocks;
	// For each alignment project the consensus confident blocks to the ref and seq coordinates
	// The main function of the loop below is to fill the attributes rfs, rfe, sqs, and sqe for each block
	// based on the previously calculated rds_f and rde_f
        for(int i=0; i < alignments_len; i++){
		bam1_t* b = alignments[i]->record;
		int j = bam_is_rev(b) ? stList_length(blocks) - 1 : 0;
		corrected_conf_blocks = stList_construct3(0, ptBlock_destruct);
        	ptCigarIt* cigar_it = ptCigarIt_construct(alignments[i]->record, true, true);
        	int block_rds_f, block_rde_f, cigar_rds_f, cigar_rde_f, rfs, rfe, sqs, sqe;
		bool del_flag = false;
        	block = stList_get(blocks, j);
		// reverse each block interval to make it work with the projecting algorithm below
		if (bam_is_rev(b)){
			block_rds_f = -1 * block->rde_f;
                        block_rde_f = -1 * block->rds_f;
                }
                else {
                        block_rds_f = block->rds_f;
                	block_rde_f = block->rde_f;
               	}
		//iterate over cigar operations from left to right (w.r.t reference)
		// There is no need to check for BAM_CHARD_CLIP and BAM_CSOFT_CLIP since
		// we know that consensus blocks will have no overlap with those operations
		// Since these blocks are the outcomes of calling find_confident_blocks() and 
		// then they are intersected with intersect_by_rd_f()
        	while(ptCigarIt_next(cigar_it)){
			// reverse each cigar operation interval to make it work with the projecting algorithm below
			if (bam_is_rev(b)){
				cigar_rds_f = -1 * cigar_it->rde_f;
				cigar_rde_f = -1 * cigar_it->rds_f;
			}
			else {
				cigar_rds_f = cigar_it->rds_f;
                                cigar_rde_f = cigar_it->rde_f;
			}
			// match/ mismatch/ insertion
			// for M and I the main algorithm is the same except
			// some minor parts that are corrected by 
			// the conditional statement "cigar_it->op == BAM_CINS ? : "
			if(cigar_it->op == BAM_CMATCH ||
			   cigar_it->op == BAM_CEQUAL ||
			   cigar_it->op == BAM_CDIFF ||
			   cigar_it->op == BAM_CINS) {
				//The loop below iterate the blocks overlap with the current cigar operation
				//There may exist multiple blocks within the current op
				while(block && block_rde_f <= cigar_rde_f){
					// update start locations of the projected coordinates
					if(cigar_rds_f <= block_rds_f && !(del_flag && cigar_rds_f == block_rds_f)){
						rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
						sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
					}
					// update end locations of the projected coordinates
					rfe = cigar_it->op == BAM_CINS ? cigar_it->rfe : cigar_it->rfs + (block_rde_f - cigar_rds_f);
					sqe = cigar_it->sqs + (block_rde_f - cigar_rds_f);
					// add block
					stList_append(corrected_conf_blocks, ptBlock_construct(rfs, 
									                       rfe, 
									                       sqs, 
											       sqe,
											       block->rds_f,
											       block->rde_f));
					// update block index; j
					if (bam_is_rev(b) && j > 0){
                                		j--;
                                                block = stList_get(blocks, j);
                                                block_rds_f = -1 * block->rde_f;
                                                block_rde_f = -1 * block->rds_f;
                        		}
					else if (!bam_is_rev(b) && j < stList_length(blocks) - 1) {
						j++;
                                		block = stList_get(blocks, j);
                                		block_rds_f = block->rds_f;
                                		block_rde_f = block->rde_f;
                        		}
					else if (j == 0 || j == stList_length(blocks) - 1){
						block = NULL;
					}
				}
				if(block == NULL) break;// no more block remaining so break the cigar iteration
				// !(del_flag && cigar_rds_f == block_rds_f) means that,
				// del_flag and cigar_rds_f == block_rds_f cannot be true at the same time
				// if both were true it means that there is a deletion (shorter than threshold)
				// on the edge of the block (This is a rare event can happen when there is a deletion
				// right before an insertion)
				if(cigar_rds_f <= block_rds_f && 
				   block_rds_f <= cigar_rde_f && 
				   !(del_flag && cigar_rds_f == block_rds_f)) {
					rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
                                        sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
				}
				del_flag = false;
			}
			//deletion
			else if (cigar_it->op == BAM_CDEL){
				int prev_j = bam_is_rev(b) ? j + 1 : j - 1;
				ptBlock* prev_block =  (0 <= prev_j && prev_j <= stList_length(blocks) - 1) ? stList_get(blocks, prev_j) : NULL;
				// if this is a short deletion where the previous block ended right before it
				// ref end coordinate has to be corrected 
				if (prev_block && cigar_it->len <= threshold){
					if ((bam_is_rev(b) && prev_block->rds_f == cigar_it->rds_f) ||
					    (!bam_is_rev(b) && prev_block->rde_f == cigar_it->rde_f)){
						prev_block->rfe = cigar_it->rfe;
					}
				}
				// if this is a short deletion where the current block starts right after it
				// ref start coordinate should be updated
				if (block && block_rds_f == cigar_rds_f && cigar_it->len <= threshold){
                                        del_flag = true;
                                        rfs = cigar_it->rfs;
					sqs = cigar_it->sqs;
                                }
			}
			// assume that the blocks do not overlap with BAM_CSOFT and BAM_CHARD 
		}
		ptCigarIt_destruct(cigar_it);
		// delete previous blocks
		stList_destruct(alignments[i]->conf_blocks);
		//sort by seq (or ref) start coordinates
		stList_sort(corrected_conf_blocks, ptBlock_cmp_sqs);
		// update confident blocks for each alignment object
		alignments[i]->conf_blocks = corrected_conf_blocks;
	}
	return stList_length(blocks);
}

/// Given a list of alignments for a single read, 
/// Determine if there is any confident block whose length is exceeding the given
//  threshold. This length can be either in the reference or read coordinates.
//  This function can be used to shorten the blocks not to be 
//  longer than a threshold
/**
   @param alignments     Array of ptAlignment structures
   @param alignments_len     Length of alignments array
   @param threshold     Length threshold
   @param sam_hdr	SAM file header structure
 */
bool needs_to_find_blocks(ptAlignment** alignments, int alignments_len, int threshold, sam_hdr_t* sam_hdr){
	stList* blocks;
	ptBlock* block;
	bool flag = false;
	int tid; const char* contig_name;
	for(int j=0; j < alignments_len; j++){
		blocks = alignments[j]->conf_blocks;
		tid = alignments[j]->record->core.tid;
                contig_name = sam_hdr_tid2name(sam_hdr, tid);
		if (blocks == NULL) return true;
		if (stList_length(blocks) == 0) return true;
		for(int i=0; i < stList_length(blocks); i++){
                	block = stList_get(blocks, i);
                	DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
			if ((block->sqe - block->sqs) > threshold || (block->rfe - block->rfs)>threshold) flag = true;
		}
	}
	return flag;
}


//
// 
//
/*
 * @param fai			Fasta index
 * @param contig_name		Contig name
 * @param alignment 		Alignment saved as a ptAlignment struct
 * @param alignment_idx		Index of the alignment 
 * @param markers		stList of markers
 * @param conf_d		Gap open probability for realignment
 * @param conf_e		Gap extension probability for realignment
 * @param conf_bw		Bandwidth for realignment
 * @param set_q			Base qualities are set to this value prior to realignment
 *
 */
void calc_local_baq(const faidx_t* fai, const char* contig_name, ptAlignment* alignment, int alignment_idx, stList* markers, double conf_d, double conf_e, double conf_bw, int set_q){

	uint8_t* tseq; // translated seq A=>0,C=>1,G=>2,T=>3,other=>4
	uint8_t* tref; // translated ref
	uint8_t* bq;
	uint8_t* block_qual;
	int* state;
	uint8_t* q; // Probability of incorrect alignment from probaln_glocal()
	probaln_par_t conf = { conf_d, conf_e, conf_bw };
	char* ref;
	char* reg;
	bam1_t* b = alignment->record;
	uint8_t* seq = bam_get_seq(b);
	uint8_t* qual = bam_get_qual(b);
	stList* blocks = alignment->conf_blocks;
	int markers_len = stList_length(markers);
	int j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
	ptMarker* marker = stList_get(markers, j);
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	ptBlock* block;
	int ref_len; int seq_len;
	DEBUG_PRINT("\t\t\t## Number of Blocks: %ld\n", stList_length(blocks));
	// block margin determines how close a marker can be to the ends of a confident block
	// if it was closer than this margin the BAQ value will be set to 0
	int block_margin = 10;
	// Start iterating over confident blocks
	for(int i=0; i < stList_length(blocks); i++){
		block = stList_get(blocks, i);
		//DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
		// Find the first cigar operation that has overlap with the block
		while((cigar_it->sqe < block->sqs) || (cigar_it->rfe < block->rfs)){
			if(ptCigarIt_next(cigar_it) == 0) break;
		}// end of while cigar_it->sqs should be now less than or equal to block->sqs
		assert(cigar_it->sqs <= block->sqs);
		assert(cigar_it->rfs <= block->rfs);
		//printf("$%d-%d : %d\n", block->sqs, block->sqe, marker->base_idx);
		while(marker &&
                      ((marker->base_idx < block->sqs + block_margin) ||
		        (marker->alignment_idx != alignment_idx))){
			if(block->sqs <= marker->base_idx) qual[marker->base_idx] = 0;
			j += bam_is_rev(alignment->record) ? -1 : 1;
                        marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
		}// end of while the marker is now the first one located within the current block (or after)
		/*if(marker == NULL) printf("NULL\n");
		else printf("%d-%d : %d\n", block->sqs, block->sqe, marker->base_idx);*/
		if(marker &&
                   (block->sqe >= marker->base_idx + block_margin) &&
		   (block->sqs <= marker->base_idx - block_margin)){
			seq_len = block->sqe - block->sqs + 1;
			ref_len = block->rfe - block->rfs + 1;
			DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
			DEBUG_PRINT("\t\t\t###     seq_len:%d\tref_len:%d\n", seq_len, ref_len);
			if (seq_len < 5 || ref_len < 5) continue; // skip small blocks
			//allocate read sequence
			tseq = (uint8_t*) malloc(seq_len);

			for(int k=0; k < seq_len; k++)
				tseq[k] = seq_nt16_int[bam_seqi(seq, block->sqs + k)];
			//allocate reference sequence
			tref = (uint8_t*) malloc(ref_len);
			reg = malloc(200);
			memset(reg,'\0',200);
			sprintf(reg, "{%s}:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);
			int len;
			ref = fai_fetch(fai, reg, &len);
			assert(len == (block->rfe - block->rfs + 1));
			for(int k=0; k < ref_len; k++)
                                tref[k] = seq_nt16_int[seq_nt16_table[(unsigned char)ref[k]]];
			//allocate the quality of this confident block
			block_qual = (uint8_t*) malloc(seq_len);
                        for(int k=0; k < seq_len; k++)
				block_qual[k] = set_q;
			//allocate neccessary arrays for HMM BAQ
			state = (int*) malloc((block->sqe - block->sqs + 1) * sizeof(int));
			q = (uint8_t*) malloc(block->sqe - block->sqs + 1);
			//DEBUG_PRINT("Starting local BAQ :))\n");
			conf.bw = abs(ref_len - seq_len) + conf_bw;
			if (probaln_glocal(tref, ref_len, 
					   tseq, seq_len,
					   block_qual, &conf, state, q) == INT_MIN) {
            			fprintf(stderr, "probaln_glocal ERROR\n");
				fprintf(stderr, "%s:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);
			}
			/*printf("##%d\n",alignment_idx);
			for(int k=0; k < seq_len; k++)
                                printf("%d\t%d\t%c\n", block->sqs + k, q[k], q[k] < 20 ? '*' : ' ');*/
			//DEBUG_PRINT("local BAQ Finished:))\n");
			// the state and q are now updated if there is any marker located within the block
			//apply BAQ to the quality array
			bq = (uint8_t*) malloc(seq_len);
			memcpy(bq, block_qual, seq_len);
			int x;
			int y;
			while(cigar_it->sqs <= block->sqe || cigar_it->rfs <= block->rfe){
				x = cigar_it->rfs - block->rfs;
				x = x < 0 ? 0 : x;
				y = cigar_it->sqs - block->sqs;
				y = y < 0 ? 0 : y;
				//DEBUG_PRINT("rf:%d\t%d\n", cigar_it->rfe, block->rfe);
				//DEBUG_PRINT("sq:%d\t%d\n", cigar_it->sqe, block->sqe);
				if (cigar_it->op == BAM_CMATCH || 
			    	    cigar_it->op == BAM_CEQUAL || 
			            cigar_it->op == BAM_CDIFF) {
					int len = min(cigar_it->len, min(cigar_it->sqe, block->sqe) - max(cigar_it->sqs, block->sqs) + 1);
					//DEBUG_PRINT("\t\t\t#len:%d\tstart:%d\tend:%d\n", len, y, y +len -1);
					for (int t = y; t < (y + len); t++) {
						assert(t < seq_len);
                        			if (((state[t]&3) != 0) || (state[t]>>2 != x + (t - y))) bq[t] = 0;
						else bq[t] = qual[block->sqs + t] < q[t] ? qual[block->sqs + t] : q[t];
                    			}
				}
				if(cigar_it->sqe <= block->sqe || cigar_it->rfe <= block->rfe){
					if(ptCigarIt_next(cigar_it) == 0) break;
				}
				else break;
			}
			for(int k=0; k < seq_len; k++) qual[block->sqs + k] = bq[k] < 94 ? bq[k] : 93;
			free(tseq);
			free(tref);
			free(state);
			free(q);
			free(bq);
			free(block_qual);
			free(ref);
			free(reg);
		}
		else { // for the markers close to the borders of the confident blocks
			while (marker &&
                      	       marker->base_idx <= block->sqe){
				if (marker->alignment_idx == alignment_idx) qual[marker->base_idx] = 0;
                        	j += bam_is_rev(alignment->record) ? -1 : 1;
                        	marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
                	}
		}
	}
	ptCigarIt_destruct(cigar_it);
}

void calc_update_baq_all(const faidx_t* fai, 
		        ptAlignment** alignments, int alignments_len, 
			stList* markers, 
			const sam_hdr_t *h,
			double conf_d, double conf_e, double conf_b, int set_q){
	const char* contig_name;
	int tid;
	for(int i=0; i < alignments_len; i++){
		tid = alignments[i]->record->core.tid;
		contig_name = sam_hdr_tid2name(h, tid);
		calc_local_baq(fai, contig_name, alignments[i], i, markers, conf_d, conf_e, conf_b, set_q);
	}
	ptMarker* marker;
	uint8_t* q;
	//DEBUG_PRINT("\t\t## start updating marker base qualities\n");
	for(int i=0; i < stList_length(markers); i++){
		marker = stList_get(markers, i);
		q = bam_get_qual(alignments[marker->alignment_idx]->record);
		marker->base_q = q[marker->base_idx];
	}
}

bool contain_supp(ptAlignment** alignments, int alignments_len){
	for(int i=0;i < alignments_len; i++){
		if (alignments[i]->record->core.flag & BAM_FSUPPLEMENTARY) return 1;
	}
	return 0;
}

void print_markers(stList* markers){
	if(stList_length(markers) == 0) DEBUG_PRINT("NO MARKERS!\n");
	DEBUG_PRINT("#alignment_idx\tseq_pos\tread_pos_f\tq\tis_match\tref_pos\n");
	for(int t =0 ; t < stList_length(markers); t++) {
        	ptMarker* m = stList_get(markers, t);
                DEBUG_PRINT("%d\t%d\t%d\t%d\t%d\t%d\n", m->alignment_idx, m->base_idx, m->read_pos_f, m->base_q, m->is_match, m->ref_pos);
        }
}

void print_contigs(ptAlignment** alignments, int alignments_len, sam_hdr_t* h){
        DEBUG_PRINT("\nalignments:\n");
	DEBUG_PRINT("#alignment_idx\tcontig_name\t\n");
        for(int t =0 ; t < alignments_len; t++) {
                bam1_t* b = alignments[t]->record;
                DEBUG_PRINT("%d\t%s\n", t, sam_hdr_tid2name(h, b->core.tid));
        }
}

static struct option long_options[] =
{
    {"inputBam", required_argument, NULL, 'i'},
    {"inputFasta", required_argument, NULL, 'f'},
    {"baq", no_argument, NULL, 'q'},
    {"gapOpen", required_argument, NULL, 'd'},
    {"gapExt", required_argument, NULL, 'e'},
    {"bandwidth", required_argument, NULL, 'b'},
    {"consensus", no_argument, NULL, 'c'},
    {"indelThreshold", required_argument, NULL, 't'},
    {"initQ", required_argument, NULL, 's'},
    {"minQ", required_argument, NULL, 'm'},
    {"primMarginScore", required_argument, NULL, 'p'},
    {"primMarginRandom", required_argument, NULL, 'r'},
    {"minScore", required_argument, NULL, 'n'},
    {"hifi", no_argument, NULL, 'x'},
    {"ont", no_argument, NULL, 'y'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char *argv[]){
	int c;
	bool baq_flag = false;
	bool consensus = false;
	int threshold = 10;
	int min_q = 20;
	int min_score = -50;
	double prim_margin_score = 50;
	double prim_margin_random = 50;
	int set_q = 40;
	double conf_d=1e-4;
	double conf_e=0.1;
	double conf_b=20;
	char* inputPath;
	char* fastaPath;
	char* vcfPath;
   	char *program;
	bool preset_ont = false;
	bool preset_hifi = false;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   	while (~(c=getopt_long(argc, argv, "i:p:f:v:qd:e:b:n:r:m:ct:s:xyh", long_options, NULL))) {
		switch (c) {
                        case 'i':
                                inputPath = optarg;
                                break;
			case 'f':
				fastaPath = optarg;
				break;
			case 'v':
				vcfPath = optarg;
				break;
			case 'x':
				preset_hifi = true;
				baq_flag = true;
			       	consensus = true;
				threshold = 10; // indel size threshold
			       	conf_d = 1e-4;
			       	conf_e = 0.1;
				conf_b = 20;
			       	min_q = 10;
			        set_q = 40;
			       	prim_margin_score = 50;
				prim_margin_random = 50;
			        min_score = -50;
				break;
			case 'y':
				preset_ont = true;
                                baq_flag = true;
                                consensus = true;
                                threshold = 20; // indel size threshold
                                conf_d = 1e-3;
                                conf_e = 0.1;
                                conf_b = 20;
                                min_q = 10;
                                set_q = 20;
                                prim_margin_score = 10;
				prim_margin_random = 10;
                                min_score = -50;
				break;
			case 'q':
				baq_flag = true;
				break;
			case 'd':
                                conf_d = atof(optarg);
                                break;
			case 'e':
                                conf_e = atof(optarg);
                                break;
			case 'b':
                                conf_b = atof(optarg);
                                break;
			case 'c':
                                consensus = true;
                                break;
			case 't':
				threshold = atoi(optarg);
				break;
			case 's':
                                set_q = atoi(optarg);
                                break;
			case 'm':
                                min_q = atoi(optarg);
                                break;
			case 'p':
				prim_margin_score = atof(optarg);
				break;
			case 'r':
                                prim_margin_random = atof(optarg);
                                break;
			case 'n':
				min_score = atoi(optarg);
				break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -f <FASTA> \n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         --inputBam, -i         Input bam file\n");
				fprintf(stderr, "         --inputFasta, -f         Input fasta file\n");
				fprintf(stderr, "         --hifi, -x         hifi preset params [-q -c -t10 -d 1e-4 -e 0.1 -b20 -m10 -s40 -p50 -r50 -n -50] (Only one of --hifi or --ont should be enabled)\n");
				fprintf(stderr, "         --ont, -y        ont present params [-q -c -t20 -d 1e-3 -e 0.1 -b20 -m10 -s20 -p10 -r10 -n -50] (Only one of --hifi or --ont should be enabled) \n");
				fprintf(stderr, "         --baq, -q         Calculate BAQ [Disabled by default]\n");
				fprintf(stderr, "         --gapOpen, -d         Gap prob [Default: 1e-4, (for ONT use 1e-2)]\n");
				fprintf(stderr, "         --gapExt, -e         Gap extension [Default: 0.1]\n");
				fprintf(stderr, "         --bandwidth, -b         DP bandwidth [Default: 20]\n");
				fprintf(stderr, "         --consensus, -c         Use consensus confident blocks [Disabled by default]\n");
				fprintf(stderr, "         --indelThreshold, -t         Indel size threshold for confident blocks [Default: 10 (for ONT use 20)]\n");
				fprintf(stderr, "         --initQ, -s         Before calculating BAQ set all base qualities to this number [Default: 40 (for ONT use 20)]\n");
				fprintf(stderr, "         --minQ, -m         Minimum base quality (or BAQ if -q is set) to be considered as a marker  [Default: 20 (for ONT use 10)]\n");
				fprintf(stderr, "         --primMarginScore, -p         Minimum margin between the consistensy score of primary and secondary alignment to select the secondary alignment [Default: 50]\n");
				fprintf(stderr, "         --primMarginRandom, -r         Maximum margin between the consistensy score of primary and secondary alignment to select one randomly [Default: 50]\n");
				
				fprintf(stderr, "         --minScore, -n         Minimum consistency score of the selected secondary alignment [Default: -50]\n");
                                return 1;
		}
	}

	stList* phased_variants = read_phased_variants(vcfPath, true);
	printf("Number of parsed phased variants = %d\n", stList_length(phased_variants));
	stList* selected_variants = filter_ref_variants(phased_variants);
	//stList_destruct(phased_variants);
	printf("Number of selected variants = %d\n", stList_length(selected_variants));
	//for(int i=0; i < stList_length(selected_variants); i++){
	//	ptVariant_print(stList_get(selected_variants, i));
	//}
	faidx_t* fai = fai_load(fastaPath);
	//printf("Make blocks!\n");
	stHash* variant_blocks = ptBlock_extract_variant_blocks(selected_variants, fai, 30);
	//printf("Blocks made!\n");
	stHashIterator* it = stHash_getIterator(variant_blocks);
	char* key;
	stList* list;
	ptBlock* block;
	while((key = stHash_getNext(it)) != NULL){
		//printf("%s\n",key);
		list = stHash_search(variant_blocks, key);
		for (int i=0 ; i < stList_length(list); i++){
			block = stList_get(list, i);
			stList* vars = (stList*) block->data;
			ptVariant* var = stList_get(vars,0);
			printf("%s\t%d\t%d\n",var->contig, block->rfs, block->rfe + 1);
			//if (100 < (block->rfe - block->rfs) ){
                        //        printf("@@@%s\t%d\t%d\n", var->contig, block->rfs, block->rfe);
                        //}
			/*if(1 < stList_length(vars)){
				printf("\n\n########%s\t%d\t%d\n",key,block->rfs, block->rfe+1);
				for(int j=0; j< stList_length(vars); j++)
					ptVariant_print(stList_get(vars,j));
				ptVariant* var = stList_get(vars,0);
				char* seq = fetch_corrected_ref_seq(fai, block, var->contig);
				printf("%s\n",seq);
				//stList* vars = (stList*) block->data;
				//printf("%d\n",stList_length(vars));
			}*/
		}
	}
	//printf("Done!\n");
	//stHash_destruct(variant_blocks);
	//stList_destruct(selected_variants);
	//exit(EXIT_FAILURE);
	if(preset_ont && preset_hifi){
		fprintf(stderr, "Presets --hifi and --ont cannot be enabled at the same time. Select one of them!\n");
		exit(EXIT_FAILURE);
	}
	samFile* fp = sam_open(inputPath, "r");
	sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	bam1_t* b = bam_init1();
	char read_name[100];
	char read_name_new[100];
	memset(read_name, '\0', 100);
	memset(read_name_new, '\0', 100);
        uint8_t* quality;
	ptMarker* marker;
	stList* markers = stList_construct3(0, free);
	int alignments_len=0;
	int32_t read_pos;
	ptAlignment* alignments[10];
	int bytes_read;
	const char* contig_name;
	bool conf_blocks_length;
	while(true) {
		bytes_read = sam_read1(fp, sam_hdr, b);
		if (bytes_read > - 1){
			strcpy(read_name_new, bam_get_qname(b));
			if (read_name[0] == '\0') {
				strcpy(read_name, read_name_new);
			}
		}
		// If read name has changed or file is finished
		if ((strcmp(read_name_new, read_name) != 0) || (bytes_read <= -1)){
			printf("###%s\n", read_name);
			printf("N alignments: %d\n", alignments_len);
			if((alignments_len > 1) && (alignments_len <= 10) && !contain_supp(alignments, alignments_len)){
				stList** read_blocks = (stList*) malloc(alignments_len * sizeof(stList*));
				for(int i=0; i < alignments_len; i++){
					char* ctg = sam_hdr_tid2name(sam_hdr, alignments[i]->record->core.tid);
					stList* blocks_ctg = stHash_search(variant_blocks, ctg);
					if (blocks_ctg == NULL){
						read_blocks[i] = stList_construct3(0, ptBlock_destruct);
					}
					else {
						/*printf("###REF:\n");
						for(int u=0; u < stList_length(blocks_ctg); u++){
                                        	        ptBlock* t = stList_get(blocks_ctg, u);
							printf("%d\t%d\n",t->rfs,t->rfe);
                               		                if (1000 < (t->rfe - t->rfs)){
                        	                                printf("larger than 1000\n");
                	                                }
        	                                }*/

						//printf("start projecting to %d\n",i);
						read_blocks[i] = project_blocks_to_read(alignments[i], blocks_ctg);

						//printf("read_blocks[%d] #= %d\n", i, stList_length(read_blocks[i]));
						/*printf("###READ:\n");
                                                for(int u=0; u < stList_length(read_blocks[i]); u++){
                                                        ptBlock* t = stList_get(read_blocks[i], u);
							printf("read\t%d\t%d\t%d\n",t->rds_f, t->rde_f, t->rde_f - t->rds_f + 1);
							printf("ref\t%d\t%d\t%d\n",t->rfs, t->rfe, t->rfe - t->rfs + 1);
                                                        if (1000 < (t->rde_f - t->rds_f)){
                                                                printf(" read larger than 1000\n");
                                                        }
							if (1000 < (t->rfe - t->rfs)){
                                                                printf("ref larger than 1000\n");
                                                        }

                                                }*/

					}
				}
				//find the maximum start point and minimum end point for alignments
				int max_rds_f = alignments[0]->rds_f;
				int min_rde_f = alignments[0]->rde_f;
				for(int i=0; i < alignments_len; i++){
					max_rds_f = max_rds_f < alignments[i]->rds_f ? alignments[i]->rds_f : max_rds_f;
					min_rde_f = alignments[i]->rde_f < min_rde_f ? alignments[i]->rde_f : min_rde_f;
				}
				// put all read blocks in a single list
				stList* all_read_blocks = stList_construct3(0, ptBlock_destruct);
				for(int i=0; i < alignments_len; i++){
                                        int len = stList_length(read_blocks[i]);
					for(int j=0; j < len; j++){
						ptBlock* b = stList_get(read_blocks[i], j);
						// add block only if it is within the region which has alignments 
						// to all haplotypes
						if ((max_rds_f <= b->rds_f) && (b->rde_f <= min_rde_f)){
							stList_append(all_read_blocks, ptBlock_copy(b, ptVariant_stList_copy));
						}
					}
                                }
				if (0 < stList_length(all_read_blocks)){
					// sort all read blocks by start pos
					stList_sort(all_read_blocks, ptBlock_cmp_rds_f);
					/*for(int u=0; u < stList_length(all_read_blocks); u++){
						ptBlock* t = stList_get(all_read_blocks,u);
						if (1000 < (t->rde_f - t->rds_f)){
							printf("larger than 1000: %d\n", t->rde_f - t->rds_f + 1);
						}
					}*/
				 	printf("blocks are sorted\n");
					// merge blocks
					stList* read_blocks_merged = merge_variant_read_blocks(all_read_blocks);
					printf("blocks are merged\n");
					// get edit distances of the variant blocks for each alignment
					int* edit_distance_list = (int*) malloc(alignments_len * sizeof(int));
					int min_edit_distance = 1000000; // a random big number
					int best_idx = -1;
					int primary_idx = -1;
					printf("Get edit distance\n");
					printf("N blocks = %d\n",stList_length(read_blocks_merged));
					for(int i=0; i < alignments_len; i++){
                                        	char* ctg = sam_hdr_tid2name(sam_hdr, alignments[i]->record->core.tid);
						edit_distance_list[i] = get_edit_distance_all_blocks(alignments[i], fai, ctg, read_blocks_merged);
						printf("i=%d,edit=%d\n",i, edit_distance_list[i]);
						if (edit_distance_list[i] < min_edit_distance){
							min_edit_distance = edit_distance_list[i];
							best_idx = i;
						}
						if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){ // primary
							primary_idx = i;
						}
					}
					if (edit_distance_list[best_idx] == edit_distance_list[primary_idx]){
						best_idx = primary_idx;
					}
					bam1_t* best = 0 <= best_idx ? alignments[best_idx]->record : NULL;
					if(best &&
					   !contain_supp(alignments, alignments_len) &&
                                           (best->core.flag & BAM_FSECONDARY)){
                                        	printf("$\t%s\n", read_name);
                                        	for(int i=0; i<alignments_len; i++){
                                                	// change primary to secondary
                                                	if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
                                                        	printf("*\t");
                                                        	//alignments[i]->record->core.flag |= BAM_FSECONDARY;
                                                	}
                                                	//change secondary to primary for the best alignment
                                                	else if (i == best_idx){
                                                        	printf("@\t");
                                                        	//alignments[i]->record->core.flag &= ~BAM_FSECONDARY;
                                                	}
                                                	else printf("!\t");
                                                	contig_name = sam_hdr_tid2name(sam_hdr, alignments[i]->record->core.tid);
                                                	printf("%d\t%s\t%ld\t%ld\n", edit_distance_list[i], contig_name, alignments[i]->record->core.pos, alignments[i]->rfe);
                                        }
                                        printf("\n");
                                }
					stList_destruct(read_blocks_merged);
					free(edit_distance_list);
				}
				stList_destruct(all_read_blocks);
				for(int i=0; i < alignments_len; i++){
					stList_destruct(read_blocks[i]);
				}
				free(read_blocks);
				//printf("free blocks\n");
			}
			if(alignments_len > 0){
				stList_destruct(markers);
                                markers = stList_construct3(0, free);
                                for(int i = 0; i < alignments_len; i++){
                                        ptAlignment_destruct(alignments[i]);
                                        alignments[i] = NULL;
                                }
                                // initialize for new alignments
                                alignments_len = 0;
			}
			/*
			// If we have at least one marker and at least two alignments
			// then we can decide which one is the best alignment
			if ((stList_length(markers) > 0) && (alignments_len > 1) && (alignments_len <= 10) && !contain_supp(alignments, alignments_len)){

DEBUG_PRINT("\n@@ READ NAME: %s\n\t$ Number of alignments: %d\t Read l_qseq: %d\n", read_name, alignments_len, alignments[0]->record->core.l_qseq);
				print_contigs(alignments, alignments_len, sam_hdr);
				//DEBUG_PRINT("Initial markers:\n");
				//print_markers(markers);
				remove_all_mismatch_markers(&markers, alignments_len);
				//DEBUG_PRINT("No all-mismatch markers:\n");
                                //print_markers(markers);
				sort_and_fill_markers(markers, alignments, alignments_len);
				//DEBUG_PRINT("Match-added markers:\n");
                                //print_markers(markers);
				filter_ins_markers(&markers, alignments, alignments_len);
				DEBUG_PRINT("Insertion-removed markers\n");
				print_markers(markers);
				if(markers && stList_length(markers) > 0){
					int flank_margin = 625;
					DEBUG_PRINT("\t# Set confident & flanking blocks\n");
					while(needs_to_find_blocks(alignments, alignments_len, 5000, sam_hdr)){
						flank_margin *= 0.8;
						DEBUG_PRINT("\t# Margin = %d\n", flank_margin);
						set_confident_blocks(alignments, alignments_len, threshold);
						set_flanking_blocks(alignments, alignments_len, markers, flank_margin);
						if (consensus) conf_blocks_length = correct_conf_blocks(alignments, alignments_len, threshold);//TODO: remove if
						if (conf_blocks_length == 0) break;
					}
					if (conf_blocks_length > 0 || consensus == false){
						if(baq_flag){
							DEBUG_PRINT("\t# There are more than 0 markers and confident blocks are not empty \n\t# So calc baq and update quality\n");
							calc_update_baq_all(fai, 
									    alignments, alignments_len,
									    markers, sam_hdr,
								            conf_d, conf_e, conf_b, set_q);
						}
						filter_lowq_markers(&markers, min_q);
						DEBUG_PRINT("High-quality markers\n");
						print_markers(markers);
						DEBUG_PRINT("#Calculate likelihoods\n");
						calc_alignment_score(markers, alignments);
					}
				}
			}
			if (alignments_len > 0){ // maybe the previous alignment was unmapped
				// get the best alignment
				int best_idx = get_best_record(alignments, alignments_len, prim_margin_score, min_score, prim_margin_random);
				bam1_t* best = 0 <= best_idx ? alignments[best_idx]->record : NULL;
				// write all alignments without any change if they are either chimeric or best alignment is primary
				if(best &&
				   !contain_supp(alignments, alignments_len) && 
				   (best->core.flag & BAM_FSECONDARY)){
					printf("$\t%s\n", read_name);
					for(int i=0; i<alignments_len; i++){
						// change primary to secondary
						if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
							printf("*\t");
							alignments[i]->record->core.flag |= BAM_FSECONDARY;
						}
						//change secondary to primary for the best alignment
						else if (i == best_idx){
							printf("@\t");
							alignments[i]->record->core.flag &= ~BAM_FSECONDARY;
						}
						else printf("!\t");
						contig_name = sam_hdr_tid2name(sam_hdr, alignments[i]->record->core.tid);
						printf("%.2f\t%s\t%ld\t%ld\n", alignments[i]->score, contig_name, alignments[i]->record->core.pos, alignments[i]->rfe);
					}
					printf("\n");
				}
				stList_destruct(markers);
				markers = stList_construct3(0, free);
				for(int i = 0; i < alignments_len; i++){
					ptAlignment_destruct(alignments[i]);
                			alignments[i] = NULL;
        			}
				// initialize for new alignments
				alignments_len = 0;
			}*/
			strcpy(read_name, read_name_new);
		}
		if (bytes_read <= -1) break; // file is finished so break
		if(b->core.flag & BAM_FUNMAP) continue; // unmapped
		if(alignments_len > 9) continue;
		alignments[alignments_len] = ptAlignment_construct(b, 0.0);
		ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
		quality = bam_get_qual(b);
		while(ptCigarIt_next(cigar_it)){
			if (cigar_it->op == BAM_CDIFF) {
				for(int j=0; j < cigar_it->len; j++){
					if(quality[cigar_it->sqs + j] < min_q) continue;
                                        if (bam_is_rev(b)) {
						marker = ptMarker_construct(alignments_len,
                                                                            cigar_it->sqs + j,
                                                                            cigar_it->rde_f - j,
                                                                            quality[cigar_it->sqs + j],
                                                                            false,
									    cigar_it->rfs + j);
					}
					else{
						marker = ptMarker_construct(alignments_len,
                                                                            cigar_it->sqs + j,
                                                                            cigar_it->rds_f + j,
                                                                            quality[cigar_it->sqs + j],
                                                                            false,
									    cigar_it->rfs + j);
					}
					stList_append(markers, marker);
				}
			}
			//set the start coordinates of the alignment
			if (alignments[alignments_len]->rfs == -1 &&
			    (cigar_it->op == BAM_CMATCH ||
                             cigar_it->op == BAM_CEQUAL ||
                             cigar_it->op == BAM_CDIFF)) {
				alignments[alignments_len]->rfs = cigar_it->rfs;
				if(bam_is_rev(b)){
					alignments[alignments_len]->rde_f = cigar_it->rde_f;
				}
				else{ 
					alignments[alignments_len]->rds_f = cigar_it->rds_f;
				}
			}
			//set the end coordinates of the alignment
                        //if the alignment ends with hard or soft clipping
			//alignments[alignments_len]->rfs != -1 is to make sure we have reached
			//the end of the alignment
			if (alignments[alignments_len]->rfe == -1 &&
			    alignments[alignments_len]->rfs != -1 &&
                            (cigar_it->op == BAM_CHARD_CLIP ||
                             cigar_it->op == BAM_CSOFT_CLIP )) {
                                alignments[alignments_len]->rfe = cigar_it->rfe;
				if(bam_is_rev(b)){
					alignments[alignments_len]->rds_f = cigar_it->rde_f + 1;
				}else{
					alignments[alignments_len]->rde_f = cigar_it->rds_f - 1;
				}
                        }
		}
		//set the start coordinates of the alignment
                //if the alignment ends with mis/matches
                if (alignments[alignments_len]->rfe == -1 &&
		    (cigar_it->op == BAM_CMATCH ||
                     cigar_it->op == BAM_CEQUAL ||
                     cigar_it->op == BAM_CDIFF)) {
			alignments[alignments_len]->rfe = cigar_it->rfe;
			if(bam_is_rev(b)){
				alignments[alignments_len]->rds_f = cigar_it->rds_f;
			}
			else{
                		alignments[alignments_len]->rde_f = cigar_it->rde_f;
			}
		}
		alignments_len += 1;
		ptCigarIt_destruct(cigar_it);
	}
	// free memory
	stHash_destruct(variant_blocks);
        stList_destruct(selected_variants);
	stList_destruct(markers);
	fai_destroy(fai);
	sam_hdr_destroy(sam_hdr);
	sam_close(fp);
	bam_destroy1(b);
}

//main();
