// sw_antidiag_opt.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#ifndef MATCH
#define MATCH     2
#define MISMATCH -1
#define GAP      -2
#endif

// Allocate zeroed, 64-byte aligned int32 array of length 'len'
static inline int32_t *aligned_int32_alloc(size_t len) {
    void *p = NULL;
    if (posix_memalign(&p, 64, len * sizeof(int32_t)) != 0) return NULL;
    // do NOT memset here for large allocations; caller will clear as needed
    return (int32_t*)p;
}

// generate random DNA-like sequence (A/C/G/T)
void generate_sequence(char *seq, int n) {
    const char alphabet[] = "ACGT";
    for (int i = 0; i < n; ++i) seq[i] = alphabet[rand() & 3];
    seq[n] = '\0';
}

// Simple anti-diagonal Smith-Waterman: returns maximum score
int smith_waterman_antidiag_simple(const char *s1, const char *s2, int len1, int len2) {
    if (len1 <= 0 || len2 <= 0) return 0;

    // max cells on any diagonal
    int max_diag = (len1 < len2) ? len1 : len2;
    if (max_diag <= 0) return 0;

    // single contiguous block for three diagonals -> better locality
    int32_t *block = aligned_int32_alloc((size_t)max_diag * 3);
    if (!block) {
        fprintf(stderr, "Allocation failed\n");
        return 0;
    }
    int32_t *prev2 = block;
    int32_t *prev1 = block + max_diag;
    int32_t *curr  = block + 2 * max_diag;

    // initially zero all (cheap once)
    memset(block, 0, (size_t)max_diag * 3 * sizeof(int32_t));

    int32_t global_max = 0;

    for (int k = 2; k <= len1 + len2; ++k) {
        int i_start = (k - len2 > 1) ? (k - len2) : 1;
        int i_end   = (k - 1 < len1) ? (k - 1) : len1;
        int diag_len = (i_end >= i_start) ? (i_end - i_start + 1) : 0;

        // Precompute previous diagonals' ranges
        int i_start_k1 = ((k-1) - len2 > 1) ? ((k-1) - len2) : 1;
        int i_end_k1   = ((k-1) - 1 < len1) ? ((k-1) - 1) : len1;
        int diag_len_k1 = (i_end_k1 >= i_start_k1) ? (i_end_k1 - i_start_k1 + 1) : 0;

        int i_start_k2 = ((k-2) - len2 > 1) ? ((k-2) - len2) : 1;
        int i_end_k2   = ((k-2) - 1 < len1) ? ((k-2) - 1) : len1;
        int diag_len_k2 = (i_end_k2 >= i_start_k2) ? (i_end_k2 - i_start_k2 + 1) : 0;

        if (diag_len == 0) {
            // rotate pointers, but we don't need to zero whole buffer
            int32_t *tmp = prev2; prev2 = prev1; prev1 = curr; curr = tmp;
            // zero only the first few elements that might be read next time (safe)
            // but we can skip memset here; we'll zero needed region before writing
            continue;
        }

        // zero only the portion of curr we will write
        memset(curr, 0, (size_t)diag_len * sizeof(int32_t));

        // local copies to avoid repeated loads inside loop
        int i_start_local = i_start;
        int i_start_k1_local = i_start_k1;
        int i_start_k2_local = i_start_k2;
        int dk1 = diag_len_k1;
        int dk2 = diag_len_k2;

        if(len1 > 8000){
            
            omp_set_num_threads(omp_get_num_procs()/2);
        }else{
            omp_set_num_threads(1);
        }

        
        #pragma omp parallel for schedule(static) reduction(max: global_max)
      
        for (int t = 0; t < diag_len; ++t) {
            int i = i_start_local + t; // 1-based
            int j = k - i;             // 1-based

            int idx_up   = (i - 1) - i_start_k1_local; // H[i-1][j] on diag k-1
            int idx_left = (i)     - i_start_k1_local; // H[i][j-1] on diag k-1
            int idx_diag = (i - 1) - i_start_k2_local; // H[i-1][j-1] on diag k-2

            int32_t diag_val = 0, up_val = 0, left_val = 0;

            if (dk2 > 0 && (unsigned)idx_diag < (unsigned)dk2)
                diag_val = prev2[idx_diag];
            if (dk1 > 0 && (unsigned)idx_up < (unsigned)dk1)
                up_val = prev1[idx_up];
            if (dk1 > 0 && (unsigned)idx_left < (unsigned)dk1)
                left_val = prev1[idx_left];

            int score = (s1[i-1] == s2[j-1]) ? MATCH : MISMATCH;
            int v_diag = diag_val + score;
            int v_up   = up_val + GAP;
            int v_left = left_val + GAP;

            int best = v_diag;
            if (v_up > best) best = v_up;
            if (v_left > best) best = v_left;
            if (best < 0) best = 0;

            curr[t] = best;
            if (best > global_max) global_max = best;
        }

        // rotate buffers
        int32_t *tmp = prev2; prev2 = prev1; prev1 = curr; curr = tmp;
    }

    free(block);
    return global_max;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <sequence_length>\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    srand(42);
    printf("Using OpenMP with %d threads\n", omp_get_num_procs()/2);

    char *a = malloc((size_t)N + 1);
    char *b = malloc((size_t)N + 1);
    if (!a || !b) { fprintf(stderr, "alloc fail\n"); return 1; }
    generate_sequence(a, N);
    generate_sequence(b, N);

    double start = omp_get_wtime();
    int score = smith_waterman_antidiag_simple(a, b, N, N);
    double end = omp_get_wtime();

    double elapsed = end - start;
    printf("N=%d  score=%d  time=%.6f s\n", N, score, elapsed);
    double gflops = (double)(10.0 * N * N) / (elapsed * 1e9);
    printf("Performance: %.2f GFLOPS\n", gflops);

    free(a); free(b);
    return 0;
}

// ✔️ If you do NOT plan to use AVX2/AVX-512 vectorization

// Then yes, you may remove the aligned allocator and nothing breaks logically.