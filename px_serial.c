#include "p7_config.h"

/* SIMD-vectorized acceleration filters, local only: */
//#include "dp_vector/msvfilter.h"          // MSV/SSV primary acceleration filter
//#include "dp_vector/vitfilter.h"          // Viterbi secondary acceleration filter
//#include "dp_vector/fwdfilter.h"          // Sparsification w/ checkpointed local Forward/Backward

/* Sparse DP, dual-mode glocal/local:    */
/*#include "dp_sparse/sparse_fwdback.h"     // sparse Forward/Backward
#include "dp_sparse/sparse_viterbi.h"     // sparse Viterbi
#include "dp_sparse/sparse_decoding.h"    // sparse Decoding
#include "dp_sparse/sparse_anchors.h"     // most probable anchor set (MPAS) 

/* Sparse anchor-set-constrained (ASC): 
#include "dp_sparse/sparse_asc_fwdback.h"   // ASC Forward/Backward
#include "dp_sparse/sparse_envelopes.h"     // Envelope inference
#include "dp_sparse/sparse_null2.h"         // Null2 score correction
#include "dp_sparse/sparse_aec_align.h"     // anchor/envelope constrained alignment
*/
#include "hmmer/src/dp_sparse/p7_engine.h"  // FIXME: we'll move the engine somewhere else, I think

#include "easel.h"
#include "esl_random.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
#include "easel.h"
#include "esl_dsqdata.h"

#include "hmmer/src/hmmer.h"
//added
#include <time.h>

uint64_t SSV_time;
uint64_t MSV_time;
uint64_t Viterbi_time;
uint64_t Forward_time;
uint64_t Backward_time;
uint64_t HMM_time;
uint64_t MSV_calls;
uint64_t Viterbi_calls;
uint64_t Forward_calls;
uint64_t Backward_calls;
uint64_t Main_calls;

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "-s",        eslARG_INT,     "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "px, the first parallel tests of H4";

int
p7_engine_Overthruster_timing(P7_ENGINE *eng, ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg)
{
  int   do_biasfilter    = (eng->params ? eng->params->do_biasfilter   : p7_ENGINE_DO_BIASFILTER);
  float sparsify_thresh  = (eng->params ? eng->params->sparsify_thresh : p7_SPARSIFY_THRESH);
  float seq_score;
  float P;
  int   status;
  //timing counters
 extern uint64_t SSV_time;
  extern uint64_t MSV_time;
  extern uint64_t Viterbi_time;
  extern uint64_t Forward_time;
  extern uint64_t Backward_time;
  extern uint64_t MSV_calls;
  extern uint64_t Viterbi_calls;
  extern uint64_t Forward_calls;
  extern uint64_t Backward_calls;
  extern uint64_t Main_calls;
  if (L == 0) return eslFAIL;

uint64_t filter_start_time, filter_end_time;

  if ((status = p7_bg_NullOne(bg, dsq, L, &(eng->nullsc))) != eslOK) return status; 
  
  filter_start_time = time(NULL);
  /*  First level: SSV and MSV filters */ 
  MSV_calls++;
  status = p7_MSVFilter(dsq, L, om, eng->fx, &(eng->mfsc));
  if (status != eslOK && status != eslERANGE) return status;

  seq_score = (eng->mfsc - eng->nullsc) / eslCONST_LOG2;          
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

  filter_end_time = time(NULL);
  MSV_time += (filter_end_time - filter_start_time);
  
  if (P > eng->F1) return eslFAIL;
  if (eng->stats) eng->stats->n_past_msv++;

  /* Biased composition HMM, ad hoc, acts as a modified null 
  if (do_biasfilter)
    {
 //     if ((status = p7_bg_FilterScore(bg, dsq, L, &(eng->biassc))) != eslOK) return status;
      seq_score = (eng->mfsc - eng->biassc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > eng->F1) return eslFAIL;
    }
  else eng->biassc = eng->nullsc;
  if (eng->stats) eng->stats->n_past_bias++;

  // TODO: in scan mode, you have to load the rest of the oprofile now,
  // configure its length model, and get GA/TC/NC thresholds.

  /* Second level: ViterbiFilter(), multihit with <om> 
  if (P > eng->F2)
    {
      if (eng->stats) eng->stats->n_ran_vit++;
      Viterbi_calls++;
      //printf("P = %.4f. Running Vit Filter\n", P);
      filter_start_time = __rdtsc();
//      status = p7_ViterbiFilter(dsq, L, om, eng->fx, &(eng->vfsc));  
      if (status != eslOK && status != eslERANGE) return status;

      seq_score = (eng->vfsc - eng->biassc) / eslCONST_LOG2;
      P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
 
      filter_end_time = __rdtsc();
      Viterbi_time += (filter_end_time - filter_start_time);   
      if (P > eng->F2) return eslFAIL;
    }
  if (eng->stats) eng->stats->n_past_vit++;


  /* Checkpointed vectorized Forward, local-only.
   
  filter_start_time = __rdtsc();
  Forward_calls++;
//  status = p7_ForwardFilter (dsq, L, om, eng->cx, &(eng->ffsc));
  if (status != eslOK) return status;

  seq_score = (eng->ffsc - eng->biassc) / eslCONST_LOG2;
  P  = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  filter_end_time = __rdtsc();
  Forward_time += (filter_end_time - filter_start_time);  
  if (P > eng->F3) return eslFAIL;
  if (eng->stats) eng->stats->n_past_fwd++;

  /* Sequence has passed all acceleration filters.
   * Calculate the sparse mask, by checkpointed vectorized decoding.
   
  Backward_calls++;
  filter_start_time = __rdtsc();
//  p7_BackwardFilter(dsq, L, om, eng->cx, eng->sm, sparsify_thresh);
  filter_end_time = __rdtsc();
  Backward_time += (filter_end_time - filter_start_time);  
*/
  return eslOK;
}



int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_BG          *bg      = NULL;
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_DSQDATA    *dd      = NULL;
  P7_ENGINE      *eng     = NULL;
  ESL_DSQDATA_CHUNK *chu = NULL;
  int             ncore   = 1;
  int  i;
  int             status;

  // init timing counters to 0
  extern uint64_t SSV_time;
  extern uint64_t MSV_time;
  extern uint64_t Viterbi_time;
  extern uint64_t Forward_time;
  extern uint64_t Backward_time;
  extern uint64_t HMM_time;
  extern uint64_t full_MSV_calls;
  extern uint64_t MSV_calls;
  extern uint64_t Viterbi_calls;
  extern uint64_t Forward_calls;
  extern uint64_t Backward_calls;
  extern uint64_t Main_calls;
/*
  SSV_time = 0;
  MSV_time = 0;
  Viterbi_time = 0;
  Forward_time = 0;
  Backward_time = 0;
  HMM_time = 0;
  full_MSV_calls = 0;
  MSV_calls=0;
  Viterbi_calls=0;
  Forward_calls =0;
  Backward_calls=0;
  Main_calls =0;
*/
  time_t program_start_time = time(NULL);
  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);

  p7_bg_SetFilter(bg, om->M, om->compo);

  /* Open s64_tequence database */
  status = esl_dsqdata_Open(&abc, seqfile, ncore, &dd);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   p7_Fail("Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error in opening dsqdata (code %d)", status);

  eng = p7_engine_Create(abc, NULL, NULL, gm->M, 400);

  uint64_t seqs = 0, chunks = 0;
   while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  
    {
	  seqs++;
      for (i = 0; i < chu->N; i++)
	{
	  chunks++;

	  p7_bg_SetLength(bg, (int) chu->L[i]);            // TODO: remove need for cast
	  p7_oprofile_ReconfigLength(om, (int) chu->L[i]); //         (ditto)
	  
	  //printf("seq %d %s\n", chu->i0+i, chu->name[i]);

	  status = p7_engine_Overthruster_timing(eng, chu->dsq[i], (int) chu->L[i], om, bg);  
	  if (status == eslFAIL) { 
	    p7_engine_Reuse(eng);
	    continue;
	  }

    Main_calls++;
    uint64_t HMM_start_time = time(NULL);
	  p7_profile_SetLength(gm, (int) chu->L[i]);
//	  status = p7_engine_Main(eng, chu->dsq[i], (int) chu->L[i], gm); 
    uint64_t HMM_end_time = time(NULL);
    HMM_time += (HMM_end_time- HMM_start_time);
	  p7_engine_Reuse(eng);

	}
      esl_dsqdata_Recycle(dd, chu);
    }
  
  esl_dsqdata_Close(dd);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  uint64_t program_end_time = time(NULL);
  /*
  printf("total_time, ssv_time, msv_time, viterbi_time, forward_time, backward_time, hmm_time\n");
  printf("%lu, %lu, %lu, %lu, %lu, %lu, %lu\n", program_end_time-program_start_time, SSV_time, MSV_time-SSV_time, Viterbi_time, Forward_time, Backward_time, HMM_time);
  printf("%lu calls to SSV\n", MSV_calls);
  printf("%lu calls to full MSV (passed SSV)\n;", full_MSV_calls);
  printf("%lu calls to Viterbi filter\n", Viterbi_calls);
  printf("%lu calls to Forward filter\n", Forward_calls);
  printf("%lu calls to Backward filter\n", Backward_calls);
  printf("%lu calls to Main comparison routine\n", Main_calls);
  */
  printf("Total time: %" PRId64 "\n", program_end_time-program_start_time);
  printf("Main_calls, sequences\n%" PRId64 ", %" PRId64 "\n", Main_calls, seqs);
  exit(0);
}





