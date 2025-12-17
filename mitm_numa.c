// mpicc -O3 -march=native -Wall -o mitm_numa mitm_numa.c -fopenmp -lnuma
// mpiexec -n 1 ./mitm_numa --n 27 --online

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <numa.h>

#include "constants.h"

#ifdef __AVX512F__ 
#define VECTOR_SIZE 32
#else
#ifdef __AVX2__
#define VECTOR_SIZE 16
#else
#ifdef __ARM_NEON
#define VECTOR_SIZE 16
#else
#ifdef __AVX__
#define VECTOR_SIZE 8
#else
#define VECTOR_SIZE 1
#endif
#endif
#endif
#endif


typedef uint64_t u64;       /* portable 64-bit integer */
typedef uint32_t u32;       /* portable 32-bit integer */
struct __attribute__((packed)) entry { u32 k; u64 v; };  /* hash table entry */

/***************************** global variables ******************************/

u64 n = 0;         /* block size (in bits) */
u64 reduce = 0; /* number of bits to reduce memory by */
u64 mask;          /* this is 2**n - 1 */

/* (P, C) : two plaintext-ciphertext pairs */
u32 P[2][2] = { {0, 0}, {0xffffffff, 0xffffffff} };
u32 C[2][2];

/************************ tools and utility functions *************************/

int world_size, rank, num_threads, numa_nodes;

double wtime() {
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

// murmur64 hash functions, tailorized for 64-bit ints / Cf. Daniel Lemire
u64 murmur64(u64 x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return x;
}

void murmur64_vectorized(u64* restrict  x, u64* restrict  result, u64 dict_size) {
    for (int i = 0; i < VECTOR_SIZE; i++) {
        u64 val = x[i];
        val ^= val >> 33;
        val *= 0xff51afd7ed558ccdull;
        val ^= val >> 33;
        val *= 0xc4ceb9fe1a85ec53ull;
        val ^= val >> 33;
        val %= (dict_size * world_size * numa_nodes);
        result[i] = val;
    }
}

/* represent n in 4 bytes */
void human_format(u64 n, char* target) {
    if (n < 1000) {
        sprintf(target, "%" PRId64, n);
        return;
    }
    if (n < 1000000) {
        sprintf(target, "%.1fK", n / 1e3);
        return;
    }
    if (n < 1000000000) {
        sprintf(target, "%.1fM", n / 1e6);
        return;
    }
    if (n < 1000000000000ll) {
        sprintf(target, "%.1fG", n / 1e9);
        return;
    }
    if (n < 1000000000000000ll) {
        sprintf(target, "%.1fT", n / 1e12);
        return;
    }
}

/******************************** SPECK block cipher **************************/

#define ROTL32(x,r) (((x)<<(r)) | (x>>(32-(r))))
#define ROTR32(x,r) (((x)>>(r)) | ((x)<<(32-(r))))

#define ER32(x,y,k) (x=ROTR32(x,8), x+=y, x^=k, y=ROTL32(y,3), y^=x)
#define DR32(x,y,k) (y^=x, y=ROTR32(y,3), x^=k, x-=y, x=ROTL32(x,8))

void Speck64128KeySchedule(const u32 K[4], u32 rk[]) {
    u32 i, D = K[3], C = K[2], B = K[1], A = K[0];
    for (i = 0;i < 27;) {
        rk[i] = A; ER32(B, A, i++);
        rk[i] = A; ER32(C, A, i++);
        rk[i] = A; ER32(D, A, i++);
    }
}

void Speck64128KeySchedule_vectorized(u32(*restrict K)[VECTOR_SIZE], u32(*restrict rk)[VECTOR_SIZE]) {
    for (u32 i = 0;i < 27;) {
        for (int vi = 0; vi < VECTOR_SIZE; vi++) {
            rk[i][vi] = K[0][vi]; ER32(K[1][vi], K[0][vi], i);
        }
        i++;
        for (int vi = 0; vi < VECTOR_SIZE; vi++) {
            rk[i][vi] = K[0][vi]; ER32(K[2][vi], K[0][vi], i);
        }
        i++;
        for (int vi = 0; vi < VECTOR_SIZE; vi++) {
            rk[i][vi] = K[0][vi]; ER32(K[3][vi], K[0][vi], i);
        }
        i++;
    }
}

void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[]) {
    u32 i;
    Ct[0] = Pt[0]; Ct[1] = Pt[1];
    for (i = 0;i < 27;)
        ER32(Ct[1], Ct[0], rk[i++]);
}

void Speck64128Encrypt_vectorized(const u32 Pt[2], u32(*restrict Ct)[VECTOR_SIZE], const u32(*restrict rk)[VECTOR_SIZE]) {
    for (int vi = 0; vi < VECTOR_SIZE; vi++) {
        Ct[0][vi] = Pt[0];
        Ct[1][vi] = Pt[1];
    }

    for (u32 i = 0;i < 27;) {
        for (int vi = 0; vi < VECTOR_SIZE; vi++) {
            ER32(Ct[1][vi], Ct[0][vi], rk[i][vi]);
        }
        i++;
    }
}

void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[]) {
    int i;
    Pt[0] = Ct[0]; Pt[1] = Ct[1];
    for (i = 26;i >= 0;)
        DR32(Pt[1], Pt[0], rk[i--]);
}

void Speck64128Decrypt_vectorized(u32(*restrict Pt)[VECTOR_SIZE], const u32 Ct[2], u32 const rk[27][VECTOR_SIZE]) {
    for (int vi = 0; vi < VECTOR_SIZE; vi++) {
        Pt[0][vi] = Ct[0];
        Pt[1][vi] = Ct[1];
    }

    for (int i = 26; i >= 0;) {
        for (int vi = 0; vi < VECTOR_SIZE; vi++) {
            DR32(Pt[1][vi], Pt[0][vi], rk[i][vi]);
        }
        i--;
    }
}

/***************************** MITM problem ***********************************/

/* f : {0, 1}^n --> {0, 1}^n.  Speck64-128 encryption of P[0], using k */
u64 f(u64 k) {
    assert((k & mask) == k);
    u32 K[4] = { k & 0xffffffff, k >> 32, 0, 0 };
    u32 rk[27];
    Speck64128KeySchedule(K, rk);
    u32 Ct[2];
    Speck64128Encrypt(P[0], Ct, rk);
    return ((u64)Ct[0] ^ ((u64)Ct[1] << 32)) & mask;
}

/* f : {0, 1}^n --> {0, 1}^n.  Speck64-128 encryption of P[0], using k */
void f_vectorized(u64 k, u64 vector[VECTOR_SIZE]) {
    assert((k & mask) == k);

    u32 K[4][VECTOR_SIZE];
    for (int i = 0; i < VECTOR_SIZE; i++) {
        u64 ki = k + i;
        K[0][i] = ki & 0xffffffff;
        K[1][i] = ki >> 32;
        K[2][i] = 0;
        K[3][i] = 0;
    }
    u32 rk[27][VECTOR_SIZE];
    Speck64128KeySchedule_vectorized(K, rk);
    u32 Ct[2][VECTOR_SIZE];
    Speck64128Encrypt_vectorized(P[0], Ct, rk);
    for (int i = 0; i < VECTOR_SIZE; i++)
        vector[i] = ((u64)Ct[0][i] ^ ((u64)Ct[1][i] << 32)) & mask;
}

/* g : {0, 1}^n --> {0, 1}^n.  speck64-128 decryption of C[0], using k */
u64 g(u64 k) {
    assert((k & mask) == k);
    u32 K[4] = { k & 0xffffffff, k >> 32, 0, 0 };
    u32 rk[27];
    Speck64128KeySchedule(K, rk);
    u32 Pt[2];
    Speck64128Decrypt(Pt, C[0], rk);
    return ((u64)Pt[0] ^ ((u64)Pt[1] << 32)) & mask;
}

/* g : {0, 1}^n --> {0, 1}^n.  speck64-128 decryption of C[0], using k */
void g_vectorized(u64 k, u64 vector[VECTOR_SIZE]) {
    assert((k & mask) == k);

    u32 K[4][VECTOR_SIZE];
    for (int i = 0; i < VECTOR_SIZE; i++) {
        u64 ki = k + i;
        K[0][i] = ki & 0xffffffff;
        K[1][i] = ki >> 32;
        K[2][i] = 0;
        K[3][i] = 0;
    }
    u32 rk[27][VECTOR_SIZE];
    Speck64128KeySchedule_vectorized(K, rk);
    u32 Pt[2][VECTOR_SIZE];
    Speck64128Decrypt_vectorized(Pt, C[0], rk);
    for (int i = 0; i < VECTOR_SIZE; i++)
        vector[i] = ((u64)Pt[0][i] ^ ((u64)Pt[1][i] << 32)) & mask;
}

bool is_good_pair(u64 k1, u64 k2) {
    u32 Ka[4] = { k1 & 0xffffffff, k1 >> 32, 0, 0 };
    u32 Kb[4] = { k2 & 0xffffffff, k2 >> 32, 0, 0 };
    u32 rka[27];
    u32 rkb[27];
    Speck64128KeySchedule(Ka, rka);
    Speck64128KeySchedule(Kb, rkb);
    u32 mid[2];
    u32 Ct[2];
    Speck64128Encrypt(P[1], mid, rka);
    Speck64128Encrypt(mid, Ct, rkb);
    return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
}

void verify_good_pairs(u64 z, u64* x, int nx, int maxres, u64* k1, u64* k2, int* nres) {
    for (int i = 0; i < nx; i++) {
        if (is_good_pair(x[i], z)) {
#pragma omp critical
            {
                if (*nres < maxres) {
                    k1[*nres] = x[i];
                    k2[*nres] = z;
                    printf("SOLUTION FOUND!\n");
                }
                (*nres) += 1;
            }
        }
    }
}

/******************************** dictionary ********************************/

/*
 * "classic" hash table for 64-bit key-value pairs, with linear probing.
 * It operates under the assumption that the keys are somewhat random 64-bit integers.
 * The keys are only stored modulo 2**32 - 5 (a prime number), and this can lead
 * to some false positives.
 */
#define EMPTY 0xffffffff
#define PRIME 0xfffffffb

void reset_dict(struct entry* dict, u64 size, int local_thread_id) {
    for(u64 i = BUFFER_SIZE *  local_thread_id; i < size; i+= num_threads * BUFFER_SIZE) {
        for(u64 j = i; j < i + BUFFER_SIZE && j < size; j++) {
            dict[j].k = EMPTY;
        }
    }

    #pragma omp barrier
}

/* allocate a hash table with `size` slots (12*size bytes) */
struct entry* dict_setup(u64 size, int node, int local_thread_id) {
    static struct entry** dicts = NULL;

    if(node == 0 && local_thread_id == 0) {
        dicts = numa_alloc_local(sizeof(*dicts) * numa_nodes);
    }

    if (rank == 0 && node == 0 && local_thread_id == 0) {
        char hdsize[8];
        human_format(size * sizeof(**dicts), hdsize);
        printf("Dictionary size: %d x %sB\n", numa_nodes, hdsize);
    }

    #pragma omp barrier
    if(local_thread_id == 0) {
        dicts[node] = numa_alloc_local(sizeof(**dicts) * size);
    }

    
    #pragma omp barrier
    if (dicts[node] == NULL) {
        fprintf(stderr, "impossible to allocate the dictionnary\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    #pragma omp barrier
    struct entry* dict = dicts[node];
    reset_dict(dict, size, local_thread_id);

    return dict;
}

/* Insert the binding key |----> value in the dictionnary */
void dict_insert_hash(struct entry* dict, u64 hash, u64 key, u64 value, u64 dict_size) {
    u32 expected = EMPTY;
    while (!__atomic_compare_exchange_n(&dict[hash].k, &expected, key % PRIME, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
        expected = EMPTY;
        hash += 1;
        if (hash == dict_size)
            hash = 0;
    }
    dict[hash].v = value;
}

/* Insert the binding key |----> value in the dictionnary */
void dict_insert(struct entry* dict, u64 key, u64 value, u64 dict_size) {
    u32 h = murmur64(key) % dict_size;
    dict_insert_hash(dict, h, key, value, dict_size);
}

/* Query the dictionnary with this `key`.  Write values (potentially)
 *  matching the key in `values` and return their number. The `values`
 *  array must be preallocated of size (at least) `maxval`.
 *  The function returns -1 if there are more than `maxval` results.
 */
int dict_probe_hash(struct entry* dict, u64 hash, u64 key, int maxval, u64 values[], u64 dict_size) {
    u32 k = key % PRIME;
    int nval = 0;
    for (;;) {
        if (dict[hash].k == EMPTY)
            return nval;
        if (dict[hash].k == k) {
            if (nval == maxval)
                return -1;
            values[nval] = dict[hash].v;
            nval += 1;
        }
        hash += 1;
        if (hash == dict_size)
            hash = 0;
    }
}

/* Query the dictionnary with this `key`.  Write values (potentially)
 *  matching the key in `values` and return their number. The `values`
 *  array must be preallocated of size (at least) `maxval`.
 *  The function returns -1 if there are more than `maxval` results.
 */
int dict_probe(struct entry* dict, u64 key, int maxval, u64 values[], u64 dict_size) {
    u64 h = murmur64(key) % dict_size;
    return dict_probe_hash(dict, h, key, maxval, values, dict_size);
}

/******************************** NUMA queues ********************************/

struct queue_entry {
    u64 key;
    struct entry value;
};

struct SPSC_Ring {
    struct queue_entry* buffer;
    //padding to avoid false sharing
    u32 head __attribute__((aligned(64)));
    char pad1[64 - sizeof(u32)];
    u32 tail __attribute__((aligned(64)));
    char pad2[64 - sizeof(u32)];
} __attribute__((aligned(64)));

struct QueuePairs {
    struct SPSC_Ring ***queues; // queues[to_node][thread_id][from_node]
    struct queue_entry ****prefill_buffers; // prefill_buffers[node][thread_id][to_node][PREFILL_BUFFER_SIZE]
    u32 ***prefill_counts; // prefill_counts[node][thread_id][to_node]
};

struct SPSC_Ring spsc_create() {
    struct SPSC_Ring r;
    r.buffer = numa_alloc_local(BUFFER_SIZE * sizeof(struct queue_entry));
    // touch the memory to ensure allocation
    for (u32 i = 0; i < BUFFER_SIZE; i +=4096 / sizeof(struct queue_entry)) {
        r.buffer[i].key = 0;
    }
    __atomic_store_n(&r.head, 0, __ATOMIC_RELEASE);
    __atomic_store_n(&r.tail, 0, __ATOMIC_RELEASE);
    return r;
}

/* Consumer: attempt to pop up to 'count' elements into dst (contiguous) */
void spsc_pop_all(struct SPSC_Ring *r, struct entry* dict, u64 dict_size, int is_insert, int maxres, u64* k1, u64* k2, int* nres) {
    u32 tail = __atomic_load_n(&r->tail, __ATOMIC_ACQUIRE);
    u32 head = __atomic_load_n(&r->head, __ATOMIC_RELAXED);
    for (; head != tail; head = (head + 1) & (BUFFER_SIZE - 1)) {
        struct queue_entry entry = r->buffer[head];
        if(is_insert) {
            dict_insert_hash(dict, entry.key, entry.value.k, entry.value.v, dict_size);
        } else {
            u64 x[256];
            int nx = dict_probe_hash(dict, entry.key, entry.value.k, 256, x, dict_size);
            assert(nx >= 0);
            verify_good_pairs(entry.value.v, x, nx, maxres, k1, k2, nres);
        }
    }

    // publish new head
    __atomic_store_n(&r->head, head, __ATOMIC_RELEASE);
}

void pairring_pop_all(struct QueuePairs* pr, int node, int local_thread_id, struct entry* dict, u64 dict_size, int is_insert, int maxres, u64* k1, u64* k2, int* nres) {
    for(int from_node = 0; from_node < numa_nodes; from_node++) {
        if(node == from_node) continue;
        spsc_pop_all(&pr->queues[node][local_thread_id][from_node], dict, dict_size, is_insert, maxres, k1, k2, nres);
    }
}

void spsc_push_many(struct SPSC_Ring *r, const struct queue_entry* src, size_t count, struct QueuePairs* pr, int node, int local_thread_id, struct entry* dict, u64 dict_size, int is_insert, int maxres, u64* k1, u64* k2, int* nres) {
    u32 head, tail, free_entries;
    while(true) {
        head = __atomic_load_n(&r->head, __ATOMIC_ACQUIRE);
        tail = __atomic_load_n(&r->tail, __ATOMIC_RELAXED);
        free_entries = BUFFER_SIZE - (tail - head);
        if (free_entries >= count)
            break;
        // not enough space in remote queue, try to pop some entries locally to avoid deadlock
        pairring_pop_all(pr, node, local_thread_id, dict, dict_size, is_insert, maxres, k1, k2, nres);
    };

    // write elements
    for (size_t i = 0; i < count; ++i) {
        r->buffer[(tail + i) & (BUFFER_SIZE - 1)] = src[i];
    }
    // publish
    __atomic_store_n(&r->tail, (tail + count) & (BUFFER_SIZE - 1), __ATOMIC_RELEASE);
}

struct QueuePairs* pairring_create(int node, int local_thread_id) {
    static struct QueuePairs* pr;
    if(node == 0 && local_thread_id == 0) {
        pr = (struct QueuePairs*) malloc(sizeof(struct QueuePairs));
        pr->queues = (struct SPSC_Ring***) malloc(numa_nodes * sizeof(struct SPSC_Ring**));
        pr->prefill_buffers = (struct queue_entry****) malloc(numa_nodes * sizeof(struct queue_entry***));
        pr->prefill_counts = (u32 ***) malloc(numa_nodes * sizeof(u32 **));
    }
    
    #pragma omp barrier
    if(local_thread_id == 0) {
        pr->queues[node] = (struct SPSC_Ring**) numa_alloc_local(num_threads * sizeof(struct SPSC_Ring*));
        pr->prefill_buffers[node] = (struct queue_entry***) numa_alloc_local(num_threads * sizeof(struct queue_entry**));
        pr->prefill_counts[node] = (u32 **) numa_alloc_local(num_threads * sizeof(u32 *));
    }
    #pragma omp barrier
    pr->queues[node][local_thread_id] = (struct SPSC_Ring*) numa_alloc_local(numa_nodes * sizeof(struct SPSC_Ring));
    pr->prefill_buffers[node][local_thread_id] = (struct queue_entry**) numa_alloc_local(numa_nodes * sizeof(struct queue_entry*));
    pr->prefill_counts[node][local_thread_id] = (u32 *) numa_alloc_local(numa_nodes * sizeof(u32));
    for(int other_node = 0; other_node < numa_nodes; other_node++) {
        if(node == other_node) continue;
        pr->queues[node][local_thread_id][other_node] = spsc_create();
        pr->prefill_buffers[node][local_thread_id][other_node] = numa_alloc_local(PREFILL_BUFFER_SIZE * sizeof(struct queue_entry));
        pr->prefill_buffers[node][local_thread_id][other_node][0].key = 0; // touch memory
        pr->prefill_counts[node][local_thread_id][other_node] = 0;
    }
    #pragma omp barrier
    return pr;
}

void pairring_flush_buffer(struct QueuePairs* pr, int node, int local_thread_id, int target_node, struct entry* dict, u64 dict_size, int is_insert, int maxres, u64* k1, u64* k2, int* nres) {
    spsc_push_many(&pr->queues[target_node][local_thread_id][node], pr->prefill_buffers[node][local_thread_id][target_node], pr->prefill_counts[node][local_thread_id][target_node], pr, node, local_thread_id, dict, dict_size, is_insert, maxres, k1, k2, nres);
    pr->prefill_counts[node][local_thread_id][target_node] = 0;
}

void pairring_flush_all_buffers(struct QueuePairs* pr, int node, int local_thread_id, struct entry* dict, u64 dict_size, int is_insert, int maxres, u64* k1, u64* k2, int* nres) {
    for(int target_node = 0; target_node < numa_nodes; target_node++) {
        if(node == target_node) continue;
        pairring_flush_buffer(pr, node, local_thread_id, target_node, dict, dict_size, is_insert, maxres, k1, k2, nres);
    }
}

void pairring_push(struct QueuePairs* pr, int node, int local_thread_id, struct queue_entry entry, int target_node, struct entry* dict, u64 dict_size, int is_insert, int maxres, u64* k1, u64* k2, int* nres) {
    pr->prefill_buffers[node][local_thread_id][target_node][pr->prefill_counts[node][local_thread_id][target_node]] = entry;
    if(++pr->prefill_counts[node][local_thread_id][target_node] == PREFILL_BUFFER_SIZE) {
        pairring_flush_buffer(pr, node, local_thread_id, target_node, dict, dict_size, is_insert, maxres, k1, k2, nres);
    }
}

/******************************************************************************/

int done = 0;
/* search the "golden collision" */
int golden_claw_search(int maxres, u64 k1[], u64 k2[], int node, int local_thread_id, struct QueuePairs* pairrings, struct entry* dict, u64 dict_size) {
    __atomic_store_n(&done, 0, __ATOMIC_RELEASE);
    u64 N = 1ull << n;
    #pragma omp barrier

    /* chunking */
    u64 chunk_count = 1ull << reduce; // number of chunks
    int n_chunk = (int)n - reduce; // low bits per chunk
    assert((1ull << n_chunk) % BUFFER_SIZE == 0);
    assert(BUFFER_SIZE % VECTOR_SIZE == 0);

    int nres = 0;
    
    for (u64 chunk = 0; chunk < chunk_count; chunk++) {
        double start = wtime();
        if (rank == 0 && node == 0 && local_thread_id == 0) {
            printf("Processing chunk %" PRId64 "/%" PRId64 "...\n", chunk + 1, chunk_count);
        }
        
        for (u64 xo = (local_thread_id + node * num_threads) * BUFFER_SIZE + (chunk << n_chunk); xo < ((chunk + 1) << n_chunk); xo += BUFFER_SIZE * num_threads * numa_nodes) {
            for(u64 x = xo; x <  xo + BUFFER_SIZE; x += VECTOR_SIZE) {
                u64 vector[VECTOR_SIZE];
                f_vectorized(x, vector);
                u64 hash[VECTOR_SIZE];
                murmur64_vectorized(vector, hash, dict_size);
                for (int j = 0; j < VECTOR_SIZE; j++) {
                    u64 f_value = vector[j];
                    u64 z = x + j;
                    int target_rank = hash[j] % world_size;
                    int target_node = (hash[j] / world_size) % numa_nodes;
                    // throw away values that do not belong to this process
                    if (target_rank == rank) {
                        if (target_node == node) {
                            // insert directly
                            dict_insert_hash(dict, hash[j] / (world_size * numa_nodes), f_value, z, dict_size);
                        } else {
                            // enqueue for remote insertion
                            struct entry entry = { .k = (u32)(f_value % PRIME), .v = z };
                            struct queue_entry qentry = { .key = hash[j] / (world_size * numa_nodes), .value = entry };
                            pairring_push(pairrings, node, local_thread_id, qentry, target_node, dict, dict_size, 1, 0, NULL, NULL, NULL);
                        }
                    }
                }
            }
            pairring_pop_all(pairrings, node, local_thread_id, dict, dict_size, 1, 0, NULL, NULL, NULL);
        }
        pairring_flush_all_buffers(pairrings, node, local_thread_id, dict, dict_size, 1, 0, NULL, NULL, NULL);
        __atomic_fetch_add(&done, 1, __ATOMIC_RELEASE);
        int local_done = __atomic_load_n(&done, __ATOMIC_ACQUIRE);
        while(local_done < numa_nodes * num_threads) {
            pairring_pop_all(pairrings, node, local_thread_id, dict, dict_size, 1, 0, NULL, NULL, NULL);
            local_done = __atomic_load_n(&done, __ATOMIC_ACQUIRE);
        }
        #pragma omp barrier
        
        __atomic_store_n(&done, 0, __ATOMIC_RELEASE);
        pairring_pop_all(pairrings, node, local_thread_id, dict, dict_size, 1, 0, NULL, NULL, NULL);
        #pragma omp barrier
        
        double mid = wtime();
        if (rank == 0 && node == 0 && local_thread_id == 0) {
            printf("Fill: %.3fs\n", mid - start);
        }

        for (u64 yo = (local_thread_id + node * num_threads) * BUFFER_SIZE; yo < N; yo += BUFFER_SIZE * num_threads * numa_nodes) {
            for (u64 y = yo; y < yo + BUFFER_SIZE; y += VECTOR_SIZE) {
                u64 x[256];
                u64 vector[VECTOR_SIZE];
                g_vectorized(y, vector);
                u64 hash[VECTOR_SIZE];
                murmur64_vectorized(vector, hash, dict_size);
                for (int j = 0; j < VECTOR_SIZE; j++) {
                    u64 z = y + j;
                    u64 y = vector[j];
                    int target_rank = hash[j] % world_size;
                    int target_node = (hash[j] / world_size) % numa_nodes;
                    // throw away values that do not belong to this process
                    if (target_rank == rank) {
                        if(target_node == node) {
                            int nx = dict_probe_hash(dict, hash[j] / (world_size * numa_nodes), y, 256, x, dict_size);
                            assert(nx >= 0);
                            verify_good_pairs(z, x, nx, maxres, k1, k2, &nres);
                        } else {
                            // enqueue for remote probing
                            struct entry entry = { .k = (u32)(y % PRIME), .v = z };
                            struct queue_entry qentry = { .key = hash[j] / (world_size * numa_nodes), .value = entry };
                            pairring_push(pairrings, node, local_thread_id, qentry, target_node, dict, dict_size, 0, maxres, k1, k2, &nres);
                        }
                    }
                }
            }
            pairring_pop_all(pairrings, node, local_thread_id, dict, dict_size, 0, maxres, k1, k2, &nres);
        }
        pairring_flush_all_buffers(pairrings, node, local_thread_id, dict, dict_size, 0, maxres, k1, k2, &nres);
        __atomic_fetch_add(&done, 1, __ATOMIC_RELEASE);
        local_done = __atomic_load_n(&done, __ATOMIC_ACQUIRE);
        while(local_done < numa_nodes * num_threads) {
            pairring_pop_all(pairrings, node, local_thread_id, dict, dict_size, 0, maxres, k1, k2, &nres);
            local_done = __atomic_load_n(&done, __ATOMIC_ACQUIRE);
        }
        #pragma omp barrier
        
        __atomic_store_n(&done, 0, __ATOMIC_RELEASE);
        pairring_pop_all(pairrings, node, local_thread_id, dict, dict_size, 0, maxres, k1, k2, &nres);
        #pragma omp barrier

        reset_dict(dict, dict_size, local_thread_id);
        #pragma omp barrier
        
        if (rank == 0 && node == 0 && local_thread_id == 0) {
            printf("Probe: %.3fs\n", wtime() - mid);
        }
    }

    if (nres > maxres)
        return -1;
    return nres;
}

/************************** command-line options ****************************/

void usage(char** argv) {
    printf("%s [OPTIONS]\n\n", argv[0]);
    printf("Options:\n");
    printf("--n N                       block size [default 24]\n");
    printf("--C0 N                      1st ciphertext (in hex)\n");
    printf("--C1 N                      2nd ciphertext (in hex)\n");
    printf("--online                    get a problem online\n");
    printf("--reduce N                  reduce memory footprint by a factor 2^n, requiring 2^n times more computation\n");
    printf("\n");
    printf("n is required, either (CO and C1) or online is required\n");
    exit(0);
}

void process_command_line_options(int argc, char** argv) {
    struct option longopts[] = {
            {"n", required_argument, NULL, 'n'},
            {"C0", required_argument, NULL, '0'},
            {"C1", required_argument, NULL, '1'},
            {"online", no_argument, NULL, 'o'},
            {"reduce", required_argument, NULL, 'r'},
            {NULL, 0, NULL, 0}
    };
    char ch;
    int set = 0;
    int online = 0;
    u64 c0, c1;
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {
        case 'n':
            n = atoi(optarg);
            mask = (1ull << n) - 1;
            break;
        case '0':
            set |= 1;
            c0 = strtoull(optarg, NULL, 16);
            C[0][0] = c0 & 0xffffffff;
            C[0][1] = c0 >> 32;
            break;
        case '1':
            set |= 2;
            c1 = strtoull(optarg, NULL, 16);
            C[1][0] = c1 & 0xffffffff;
            C[1][1] = c1 >> 32;
            break;
        case 'o':
            online = 1;
            break;
        case 'r':
            reduce = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Unknown option\n");
            exit(1);
        }
    }
    if (n == 0 || (set != 3 && !online)) {
        usage(argv);
        exit(1);
    }
    if (online && rank == 0) {
        /* fetch problem from server */
        srand(time(NULL));
        rand();
        int version = rand() % 1000;
        char url[256];
        sprintf(url, "https://ppar.tme-crypto.fr/mathis.poppe.%d/%llu", version, (unsigned long long)n);
        printf("Fetching problem from %s\n", url);
        char filename[128];
        sprintf(filename, "%d_%llu.txt", (int)version, (unsigned long long)n);

        char command[512];
        sprintf(command, "curl -s %s > %s", url, filename);
        system(command);
        FILE* f = fopen(filename, "r");
        if (!f) {
            fprintf(stderr, "error opening file containing problem");
            exit(1);
        }
        // parse file
        /* example content:
        C0 = (668ad8ab, 7dc4d315)
        C1 = (17900ae3, c0eaa928)
        */
        char line[256];
        while (fgets(line, sizeof(line), f)) {
            if (line[0] == 'C' && line[1] == '0') {
                sscanf(line, "C0 = (%x, %x)", &C[0][0], &C[0][1]);
            } else if (line[0] == 'C' && line[1] == '1') {
                sscanf(line, "C1 = (%x, %x)", &C[1][0], &C[1][1]);
            }
        }
        fclose(f);
        // remove file
        remove(filename);
    }
    if (online) {
        // broadcast C0 and C1 to all processes
        MPI_Bcast(&C[0][0], 4, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    }
}

/******************************************************************************/

u64 k1_global[16], k2_global[16];
int nres_global = 0;

int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (provided < MPI_THREAD_MULTIPLE) {
        if (rank == 0) {
            fprintf(stderr, "Error: MPI does not provide needed threading level\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double start_time = wtime();

    numa_nodes = numa_num_configured_nodes();
    if (numa_nodes <= 1) {
        if (rank == 0) {
            fprintf(stderr, "Error: NUMA nodes not detected\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    numa_set_localalloc();

    process_command_line_options(argc, argv);
    num_threads = omp_get_max_threads() / numa_nodes;

    if (rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n",
            (int)n, C[0][0], C[0][1], C[1][0], C[1][1]);
        printf("Using %d processes, %d numa nodes, %d threads and vector size %d\n", world_size, numa_nodes, num_threads, VECTOR_SIZE);
    }
    
    #pragma omp parallel
    {
        // seperate threads per NUMA node
        int id = omp_get_thread_num();
        int node = id % numa_nodes;
        int local_thread_id = id / numa_nodes;

        // bind thread to node
        numa_run_on_node(node);
        numa_set_localalloc();

        // setup dictionary
        u64 dict_size = 1.125 * (1ull << (n - reduce)) / (world_size * numa_nodes);
        struct entry* dict = dict_setup(dict_size, node, local_thread_id);
        // setup queues
        struct QueuePairs* pairrings = pairring_create(node, local_thread_id);

        /* search */
        u64 k1[16], k2[16];
        int nkey = golden_claw_search(16, k1, k2, node, local_thread_id, pairrings, dict, dict_size);
        int index = __atomic_fetch_add(&nres_global, nkey, __ATOMIC_RELAXED);
        for (int i = 0; i < nkey; i++) {
            k1_global[index + i] = k1[i];
            k2_global[index + i] = k2[i];
        }
    }
    
    if (rank == 0) {
        int nres_per_process[world_size];
        MPI_Gather(&nres_global, 1, MPI_INT, nres_per_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int offsets[world_size];
        offsets[0] = 0;
        for (int i = 1; i < world_size; i++) {
            offsets[i] = offsets[i - 1] + nres_per_process[i - 1];
        }
        MPI_Gatherv(MPI_IN_PLACE, nres_global, MPI_UINT64_T, k1_global, nres_per_process, offsets, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, nres_global, MPI_UINT64_T, k2_global, nres_per_process, offsets, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        nres_global = offsets[world_size - 1] + nres_per_process[world_size - 1];
    } else {
        MPI_Gather(&nres_global, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(k1_global, nres_global, MPI_UINT64_T, NULL, NULL, NULL, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(k2_global, nres_global, MPI_UINT64_T, NULL, NULL, NULL, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    }

    double end_time = wtime();

    if (rank == 0) {
        printf("Total time: %.3fs\n", end_time - start_time);
        assert(nres_global > 0);
        /* validation */
        for (int i = 0; i < nres_global; i++) {
            assert(f(k1_global[i]) == g(k2_global[i]));
            assert(is_good_pair(k1_global[i], k2_global[i]));
            printf("Solution found: (%" PRIx64 ", %" PRIx64 ") [checked OK]\n", k1_global[i], k2_global[i]);
        }
    }
    MPI_Finalize();
    return 0;
}