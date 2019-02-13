#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include "htslib/sam.h"
#pragma clang diagnostic pop

#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


// no it doesn't work if you just modify these
#define GC_BINCOUNT 50
// #define GC_BINSIZE (1./GC_BINCOUNT)
#define GC_LOWERBOUND 20 // 0.4
#define GC_TELOMERIC_LOWERBOUND 0.48
#define GC_UPPERBOUND 29 // 0.6
#define GC_TELOMERIC_UPPERBOUND 0.52

#define TEL_MOTIF_CUTOFF 7

static double calcTelLength(const size_t *const telcounts, const size_t telcounts_len, const size_t gc_telomeric_count)
{
    double acc = 0;
    for (size_t i = TEL_MOTIF_CUTOFF; i < telcounts_len; ++i)
    {
        acc += telcounts[i];
    }

    if(gc_telomeric_count == 0){
        return (double)NAN;
    }
    // fprintf(stderr, "telacc: %f\n", acc);

    // * GENOME_LENGTH_AT_TEL_GC / LENGTH_UNIT / TELOMERE_ENDS
    return (acc / gc_telomeric_count) * 332720800 / 1000 / 46;
}

static size_t countMotif(const char *read)
{
    // this should use bm or better.
    const char* const motif1 = "TTAGGG";
    const char* const motif2 = "CCCTAA";
    size_t motif_length = strlen(motif1);

    size_t motifcount1 = 0;
    size_t motifcount2 = 0;

    const char *p = strstr(read, motif1);
    while(p)
    {
        motifcount1 += 1;
        p = strstr(p + motif_length, motif1);
    }

    p = strstr(read, motif2);
    while(p)
    {
        motifcount2 += 1;
        p = strstr(p + motif_length, motif2);
    }

    return motifcount1 > motifcount2 ? motifcount1 : motifcount2;
}

static double calcGC(const char *read, int32_t len)
{
    int32_t num_gc = 0;
    for(int32_t i = 0; i < len; ++i)
    {
        if(read[i] == 'C' || read[i] == 'G')
            ++num_gc;
    }
    return ((double)num_gc) / len;
}

int main(int argc, char const *argv[])
{
    if (argc != 2)
    {
        exit(EXIT_FAILURE);
    }

    samFile *f = sam_open(argv[1], "r");
    if (!f) exit(EXIT_FAILURE);

    bam_hdr_t *h = sam_hdr_read(f);
    if (!h) exit(EXIT_FAILURE);

    bam1_t *b1 = bam_init1();
    if (!b1) exit(EXIT_FAILURE);

    uint64_t ntotal = 0; // number of reads scanned in bam (we skip some reads, see below)
    uint64_t nmapped = 0;
    uint64_t nduplicates = 0;
    char *read = NULL;

    // even for long reads of many repeats memory usage is not too bad
    size_t telcounts_len = 2;
    size_t *telcounts = calloc(telcounts_len, sizeof(*telcounts));

    size_t gccounts_len = GC_BINCOUNT + 1;
    size_t *gccounts = calloc(gccounts_len, sizeof(*gccounts));
    size_t gc_telomeric_count = 0;

    int read_err;
    while ((read_err = sam_read1(f, h, b1)) >= 0)
    {
        ntotal++;

        if (!(b1->core.flag & BAM_FUNMAP)) {
            nmapped += 1;
        }

        if (b1->core.flag & BAM_FDUP) {
            nduplicates += 1;
        }

        read = realloc(read, b1->core.l_qseq + 1);
        if (!read) exit(EXIT_FAILURE);

        for (int i = 0; i < b1->core.l_qseq; ++i)
        {
            // TODO don't do the nibble -> char conversion, directly measure
            read[i] = seq_nt16_str[bam_seqi(bam_get_seq(b1), i)];
        }
        read[b1->core.l_qseq] = 0;

        double gc = calcGC(read, b1->core.l_qseq);
        // fprintf(stderr, "%f\n", gc);
        // get index for GC bin.
        int idx = (int)(gc * GC_BINCOUNT);
        // assert(idx >=0 && idx <= ScanParameters::GC_BIN_N-1);
        gccounts[idx] += 1;
        gc_telomeric_count += ((GC_TELOMERIC_LOWERBOUND <= gc) && (gc < GC_TELOMERIC_UPPERBOUND));

        size_t ptn_count = countMotif(read);
        if (ptn_count >= telcounts_len)
        {
            telcounts = realloc(telcounts, (ptn_count + 1) * sizeof(*telcounts));
            if (!telcounts) exit(EXIT_FAILURE);
            for (size_t i = telcounts_len; i < ptn_count + 1; ++i)
            {
                telcounts[i] = 0;
            }
            telcounts_len = ptn_count + 1;
        }
        telcounts[ptn_count] += 1;
    }
    bam_destroy1(b1);
    bam_hdr_destroy(h);
    sam_close(f);

    // header TODO: missing 'sample'
    printf("ReadGroup\tLibrary\tTotal\tMapped\tDuplicates\tLENGTH_ESTIMATE");
    for (int i = 0; i < 17; ++i)
    {
        printf("\tTEL%d", i);
    }
    for (int i = 0; i < 10; ++i)
    {
        printf("\tGC%d", i);
    }
    printf("\n");

    // LABEL_RG
    printf("UNKNOWN\t");
    // LABEL_LB
    printf("UNKNOWN\t");
    // LABEL_TOTAL
    printf("%lu\t", ntotal);
    // LABEL_MAPPED
    printf("%lu\t", nmapped);
    // LABEL_DUP
    printf("%lu\t", nduplicates);
    // LABEL_LEN
    printf("%f", calcTelLength(telcounts, telcounts_len, gc_telomeric_count));

    for(int i = 0; i < 17; i++){
        printf("\t%lu", telcounts[i]);
    }
    for(int i = GC_LOWERBOUND; i <= GC_UPPERBOUND; i++){
        printf("\t%lu", gccounts[i]);
    }
    printf("\n");

    free(gccounts);
    free(telcounts);
    free(read);
    return 0;
}
