#pragma once

/* number of MPI groups in fill part */
#define GROUPS_COUNT_FILL 2

/* number of MPI groups in probe part */
#define GROUPS_COUNT_PROBE 2

/* buffer size for each process before sending data over the network */
#define BUFFER_SIZE 8192

/* number of buffers that can be sent without waiting for completion */
#define BSEND_AMOUNT 1000
