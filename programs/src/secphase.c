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
#include "ptBlock.h"
#include "ptVariant.h"
#include "ptAlignment.h"
#include "ptMarker.h"

#define SCORE_TYPE_MARKER		0
#define SCORE_TYPE_EDIT_DISTANCE	1

void print_alignment_scores(ptAlignment** alignments, int alignments_len, int score_type, FILE* fout){
	for(int i=0; i< alignments_len; i++){
		// change primary to secondary
		if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
			fprintf(fout, "*\t");
			alignments[i]->record->core.flag |= BAM_FSECONDARY;
		}
		//change secondary to primary for the best alignment
		else if (i == best_idx){
			fprintf(fout, "@\t");
			alignments[i]->record->core.flag &= ~BAM_FSECONDARY;
		}
		else fprintf(fout, "!\t");
		if (score_type == SCORE_TYPE_MARKER){
			int edit_distance = -1 * alignments[i]->score;
			fprintf(fout, "%d\t%s\t%ld\t%ld\n", edit_distance, alignments[i]->contig, alignments[i]->record->core.pos, alignments[i]->rfe);
		}else if( score_type == SCORE_TYPE_EDIT_DISTANCE){
			fprintf(fout, "%.2f\t%s\t%ld\t%ld\n", alignments[i]->score, alignments[i]->contig, alignments[i]->record->core.pos, alignments[i]->rfe);
		}
	}
	fprintf(fout, "\n");
}

static struct option long_options[] =
{
    {"inputBam", required_argument, NULL, 'i'},
    {"inputFasta", required_argument, NULL, 'f'},
    {"inputVcf", required_argument, NULL, 'v'},
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
    {"minVariantMargin", required_argument, NULL, 'g'},
    {"prefix", required_argument, NULL, 'P'},
    {"debug", no_argument, NULL, 'D'},
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
	int min_var_margin = 50;
	double conf_d=1e-4;
	double conf_e=0.1;
	double conf_b=20;
	char* inputPath;
	char* fastaPath;
	char* vcfPath;
	char prefix[50] = "secphase";
	bool debug = false;
   	char *program;
	bool preset_ont = false;
	bool preset_hifi = false;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   	while (~(c=getopt_long(argc, argv, "i:p:P:Df:v:qd:e:b:n:r:m:ct:s:g:xyh", long_options, NULL))) {
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
			case 'P':
				prefix = optarg;
				break;
			case 'D':
				debug = true;
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
			case 'g':
				min_var_margin = atoi(optarg);
				break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -f <FASTA> \n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         --inputBam, -i         Input BAM file\n");
				fprintf(stderr, "         --inputFasta, -f         Input FASTA file\n");
				fprintf(stderr, "         --inputVcf, -v         Input phased VCF file\n");
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
				
				fprintf(stderr, "         --minScore, -n         Minimum marker score of the selected secondary alignment [Default: -50]\n");
				fprintf(stderr, "         --minVariantMargin, -g         Minimum margin for creating blocks around phased variants [Default: 50]\n");
                                return 1;
		}
	}


	faidx_t* fai = fai_load(fastaPath);
	stHash* variant_ref_blocks_per_contig = ptBlock_extract_variant_blocks(selected_variants, fai, min_var_margin);
	
	if(debug == true){
		char bed_path_ref_blocks[50];
		snprintf(bed_path_ref_blocks, "%s.variant_ref_blocks.bed", prefix);
		ptVariant_save_variant_ref_blocks(variant_ref_blocks_per_contig, bed_path_ref_blocks);
	}
	
	if(preset_ont && preset_hifi){
		fprintf(stderr, "[%s] Presets --hifi and --ont cannot be enabled at the same time. Select only one of them!\n", get_timestamp());
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

			if((alignments_len > 1) && 
				(alignments_len <= 10) && 
				!contain_supp(alignments, alignments_len) &&
				overlap_variant_ref_blocks(variant_ref_blocks_per_contig, alignments, alignments_len)){

				stList* merged_variant_read_blocks = ptVariant_get_merged_variant_read_blocks(variant_ref_blocks_per_contig, alignments, alignments_len);
				set_edit_distances_as_scores(merged_variant_read_blocks, alignments, alignments_len);
				int best_idx = get_best_record_index(alignments, alignments_len, 0, -100, 0);
				bam1_t* best = 0 <= best_idx ? alignments[best_idx]->record : NULL;
				if(best &&
					!contain_supp(alignments, alignments_len) &&
                                        (best->core.flag & BAM_FSECONDARY)){
					fprintf(fout, "$\t%s\n", read_name);
                                }
				stList_destruct(merged_variant_read_blocks);
				stHash_destruct(variant_ref_blocks_per_contig);
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
		alignments[alignments_len] = ptAlignment_construct(b, sam_hdr);
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
