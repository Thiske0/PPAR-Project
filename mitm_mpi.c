// mpicc -O3 -march=native -Wall -o mitm_mpi mitm_mpi.c -fopenmp
// mpiexec -n 1 ./mitm_mpi --n 27 --online

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


#define BUFFER_SIZE 8192  /* buffer size for each process before sending data over the network */

#ifdef __AVX512F__ 
#define VECTOR_SIZE 32
#else
#ifdef __AVX2__
#define VECTOR_SIZE 16
#else
#ifdef __AVX__
#define VECTOR_SIZE 8
#else
#define VECTOR_SIZE 1
#endif
#endif
#endif


typedef uint64_t u64;       /* portable 64-bit integer */
typedef uint32_t u32;       /* portable 32-bit integer */
struct __attribute__((packed)) entry { u32 k; u64 v; };  /* hash table entry */

/***************************** global variables ******************************/

u64 n = 0;         /* block size (in bits) */
u64 mask;          /* this is 2**n - 1 */

u64 dict_size;     /* number of slots in the hash table */
struct entry* A;   /* the hash table */

/* (P, C) : two plaintext-ciphertext pairs */
u32 P[2][2] = { {0, 0}, {0xffffffff, 0xffffffff} };
u32 C[2][2];

/************************ tools and utility functions *************************/

int world_size, rank;

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

void murmur64_vectorized(u64* restrict  x, u64* restrict  result) {
    for (int i = 0; i < VECTOR_SIZE; i++) {
        u64 val = x[i];
        val ^= val >> 33;
        val *= 0xff51afd7ed558ccdull;
        val ^= val >> 33;
        val *= 0xc4ceb9fe1a85ec53ull;
        val ^= val >> 33;
        val %= (dict_size * world_size);
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

/******************************** dictionary ********************************/

/*
 * "classic" hash table for 64-bit key-value pairs, with linear probing.
 * It operates under the assumption that the keys are somewhat random 64-bit integers.
 * The keys are only stored modulo 2**32 - 5 (a prime number), and this can lead
 * to some false positives.
 */
static const u32 EMPTY = 0xffffffff;
static const u64 PRIME = 0xfffffffb;

/* allocate a hash table with `size` slots (12*size bytes) */
void dict_setup(u64 size) {
    dict_size = size;
    if (rank == 0) {
        char hdsize[8];
        human_format(dict_size * sizeof(*A), hdsize);
        printf("Dictionary size: %sB\n", hdsize);
    }

    A = malloc(sizeof(*A) * dict_size);
    if (A == NULL) {
        fprintf(stderr, "impossible to allocate the dictionnary");
        exit(1);
    }
    for (u64 i = 0; i < dict_size; i++)
        A[i].k = EMPTY;
}

/* Insert the binding key |----> value in the dictionnary */
void dict_insert_hash(u32 hash, u64 key, u64 value) {
    u32 expected = EMPTY;
    while (!__atomic_compare_exchange_n(&A[hash].k, &expected, key % PRIME, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
        expected = EMPTY;
        hash += 1;
        if (hash == dict_size)
            hash = 0;
    }
    A[hash].v = value;
}

/* Insert the binding key |----> value in the dictionnary */
void dict_insert(u64 key, u64 value) {
    u32 h = murmur64(key) % dict_size;
    dict_insert_hash(h, key, value);
}

/* Query the dictionnary with this `key`.  Write values (potentially)
 *  matching the key in `values` and return their number. The `values`
 *  array must be preallocated of size (at least) `maxval`.
 *  The function returns -1 if there are more than `maxval` results.
 */
int dict_probe_hash(u32 hash, u64 key, int maxval, u64 values[]) {
    u32 k = key % PRIME;
    int nval = 0;
    for (;;) {
        if (A[hash].k == EMPTY)
            return nval;
        if (A[hash].k == k) {
            if (nval == maxval)
                return -1;
            values[nval] = A[hash].v;
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
int dict_probe(u64 key, int maxval, u64 values[]) {
    u64 h = murmur64(key) % dict_size;
    return dict_probe_hash(h, key, maxval, values);
}

/************************** buffer ****************************/

#define FULL_BUFFER_TAG 0x4
#define PARTIAL_BUFFER_TAG 0x8

struct buffer_entry {
    u32 key; // 32-bit key (modulo dict_size)
    struct entry value; // value to insert at this key
};

struct buffer_entry** buffers;
u32* buffer_indices;
u32* writers;
struct buffer_entry* recieve_buffer;

void init_buffers() {
    buffer_indices = malloc(world_size * sizeof(*buffer_indices));
    writers = malloc(world_size * sizeof(*writers));
    for (int i = 0; i < world_size; i++) {
        buffer_indices[i] = 0;
        writers[i] = 0;
    }

    buffers = malloc(world_size * sizeof(*buffers));
    for (int i = 0; i < world_size; i++) {
        buffers[i] = malloc(BUFFER_SIZE * sizeof(**buffers));
    }
    recieve_buffer = malloc(BUFFER_SIZE * sizeof(*recieve_buffer));

    int size = (BUFFER_SIZE * sizeof(*recieve_buffer) + MPI_BSEND_OVERHEAD) * 40; // probably overkill
    void* bsend_buf = malloc(size);
    MPI_Buffer_attach(bsend_buf, size);
}

void swap(int a, int b) {
    struct buffer_entry* temp = buffers[a];
    buffers[a] = buffers[b];
    buffers[b] = temp;
}

void buffer_add(int target, struct buffer_entry entry) {
    __atomic_fetch_add(&writers[target], 1, __ATOMIC_ACQUIRE);
    u32 index = __atomic_fetch_add(&buffer_indices[target], 1, __ATOMIC_ACQUIRE);
    if (index < BUFFER_SIZE - 1) {
        buffers[target][index] = entry;
        __atomic_fetch_sub(&writers[target], 1, __ATOMIC_RELEASE);
        return;
    }
    // we need to send the buffer
    if (index >= BUFFER_SIZE) {
        __atomic_fetch_sub(&writers[target], 1, __ATOMIC_RELEASE);
        // wait for the send to complete
        while (__atomic_load_n(&buffer_indices[target], __ATOMIC_ACQUIRE) > BUFFER_SIZE) {
            // busy wait
        }
        // try again
        buffer_add(target, entry);
        return;
    }
    buffers[target][index] = entry;
    // wait for other writers to complete
    while (__atomic_load_n(&writers[target], __ATOMIC_ACQUIRE) > 1) {
        // busy wait
    }
    // send the buffer

#pragma omp critical
    {
        // swap buffer_indices[target]
        swap(target, rank);
        __atomic_store_n(&buffer_indices[target], 0, __ATOMIC_RELEASE);
        MPI_Bsend(buffers[rank], sizeof(**buffers) * BUFFER_SIZE, MPI_BYTE, target, FULL_BUFFER_TAG, MPI_COMM_WORLD);
    }
    __atomic_fetch_sub(&writers[target], 1, __ATOMIC_RELEASE);
}

int buffer_probe(u64 hash, u64 key, int maxval, u64 values[]) {
    int target = hash % world_size;
    if (target == rank) {
        return dict_probe_hash(hash / world_size, key, maxval, values);
    }
    struct buffer_entry entry = { .key = (u32)(hash / world_size), .value = {.k = (u32)(key % PRIME), .v = 0 } };
    buffer_add(target, entry);
    return 0;
}

void buffer_insert(u64 hash, u64 key, u64 value) {
    int target = hash % world_size;
    if (target == rank) {
        dict_insert_hash(hash / world_size, key, value);
        return;
    }
    struct buffer_entry entry = { .key = (u32)(hash / world_size), .value = {.k = (u32)(key % PRIME), .v = value } };
    buffer_add(target, entry);
}

void try_recieve_buffers() {
    int flag;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, FULL_BUFFER_TAG, MPI_COMM_WORLD, &flag, &status);
    while (flag) {
        MPI_Recv(recieve_buffer, sizeof(**buffers) * BUFFER_SIZE, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // insert entries
        for (u32 i = 0; i < BUFFER_SIZE; i++) {
            dict_insert_hash(recieve_buffer[i].key,
                ((u64)recieve_buffer[i].value.k), recieve_buffer[i].value.v);
        }

        MPI_Iprobe(MPI_ANY_SOURCE, FULL_BUFFER_TAG, MPI_COMM_WORLD, &flag, &status);
    }
}

void send_recieve_remaining_buffers() {
    MPI_Request requests[world_size];
    for (u32 i = 0;i < world_size;i++) {
        if (i == rank) continue;
        int count = buffer_indices[i];
        MPI_Isend(buffers[i], sizeof(**buffers) * count, MPI_BYTE, i, PARTIAL_BUFFER_TAG, MPI_COMM_WORLD, &requests[i]);
    }
    for (u32 i = 0;i < world_size;i++) {
        if (i == rank) continue;
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, PARTIAL_BUFFER_TAG, MPI_COMM_WORLD, &status);
        int count;
        MPI_Get_count(&status, MPI_BYTE, &count);
        MPI_Recv(recieve_buffer, count, MPI_BYTE, status.MPI_SOURCE, PARTIAL_BUFFER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        u32 nentries = count / sizeof(**buffers);
        for (u32 j = 0;j < nentries;j++) {
            dict_insert_hash(recieve_buffer[j].key,
                ((u64)recieve_buffer[j].value.k), recieve_buffer[j].value.v);
        }

    }
    //only try to revieve once we know that all sends are posted
    //we know that all are posted becuase we recieved the partial buffers
    try_recieve_buffers();
    for (u32 i = 0;i < world_size;i++) {
        if (i == rank) continue;
        MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
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

/******************************************************************************/

/* search the "golden collision" */
int golden_claw_search(int maxres, u64 k1[], u64 k2[]) {
    double start = wtime();
    u64 N = 1ull << n;
    assert(N % BUFFER_SIZE == 0);
    assert(BUFFER_SIZE % VECTOR_SIZE == 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (u64 y = rank * BUFFER_SIZE; y < N; y += world_size * BUFFER_SIZE) {
        for (u64 x = 0; x < BUFFER_SIZE; x += VECTOR_SIZE) {
            u64 vector[VECTOR_SIZE];
            f_vectorized(y + x, vector);
            u64 hash[VECTOR_SIZE];
            murmur64_vectorized(vector, hash);
            for (int j = 0; j < VECTOR_SIZE; j++) {
                u64 z = y + x + j;
                u64 y = vector[j];
                buffer_insert(hash[j], y, z);
            }
        }
        try_recieve_buffers();
    }
    send_recieve_remaining_buffers();

    double mid = wtime();
    if (rank == 0) {
        printf("Fill: %.1fs\n", mid - start);
    }

    int nres = 0;
    u64 ncandidates = 0;
#pragma omp parallel for schedule(dynamic, 4096) reduction(+:ncandidates)
    for (u64 zv = 0; zv < N; zv += VECTOR_SIZE) {
        u64 x[256];
        u64 vector[VECTOR_SIZE];
        g_vectorized(zv, vector);
        u64 hash[VECTOR_SIZE];
        murmur64_vectorized(vector, hash);
        for (int j = 0; j < VECTOR_SIZE; j++) {
            u64 z = zv + j;
            u64 y = vector[j];
            if (hash[j] % world_size != rank)
                continue;
            int nx = dict_probe_hash(hash[j] / world_size, y, 256, x);
            assert(nx >= 0);
            ncandidates += nx;
            for (int i = 0; i < nx; i++)
                if (is_good_pair(x[i], z)) {
#pragma omp critical
                    {
                        if (nres < maxres) {
                            k1[nres] = x[i];
                            k2[nres] = z;
                            printf("SOLUTION FOUND!\n");
                        }
                        nres += 1;
                    }
                }
        }
    }
    if (rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, &ncandidates, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid, ncandidates);
    } else {
        MPI_Reduce(&ncandidates, NULL, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
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

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    process_command_line_options(argc, argv);
    if (rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n",
            (int)n, C[0][0], C[0][1], C[1][0], C[1][1]);
        int num_threads = omp_get_max_threads();
        printf("Using %d processes, %d threads and vector size %d\n", world_size, num_threads, VECTOR_SIZE);
    }

    dict_setup(1.125 * (1ull << n) / world_size);
    init_buffers();

    /* search */
    u64 k1[16], k2[16];
    int nkey = golden_claw_search(16, k1, k2);

    if (rank == 0) {
        //assert(nkey > 0);
        /* validation */
        for (int i = 0; i < nkey; i++) {
            assert(f(k1[i]) == g(k2[i]));
            assert(is_good_pair(k1[i], k2[i]));
            printf("Solution found: (%" PRIx64 ", %" PRIx64 ") [checked OK]\n", k1[i], k2[i]);
        }
    }
    MPI_Finalize();
    return 0;
}