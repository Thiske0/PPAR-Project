// gcc -Wall -o make_header make_header.c
// ./make_header --fill-groups 2 --probe-groups 2 --buffer-size 8192 --bsend-amount 1000 --prefill-buffer-size 128

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

int groups_count_fill = 0, groups_count_probe = 0, buffer_size = 0, bsend_amount = 0, prefill_buffer_size = 0, queue_size = 0, block_size = 0;

void usage(char** argv) {
    printf("%s [OPTIONS]\n\n", argv[0]);
    printf("Options:\n");
    printf("--fill-gropus N\n");
    printf("--probe-groups N\n");
    printf("--buffer-size N\n");
    printf("--prefill-buffer-size N\n");
    printf("--bsend-amount N\n");
    printf("\n");
    printf("every option is required\n");
    exit(0);
}

void process_command_line_options(int argc, char** argv) {
    struct option longopts[] = {
            {"fill-groups", required_argument, NULL, 'f'},
            {"probe-groups", required_argument, NULL, 'p'},
            {"block-size", required_argument, NULL, 'B'},
            {"buffer-size", required_argument, NULL, 'b'},
            {"queue-size", required_argument, NULL, 'q'},
            {"prefill-buffer-size", required_argument, NULL, 'r'},
            {"bsend-amount", required_argument, NULL, 's'},
            {NULL, 0, NULL, 0}
    };
    char ch;
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {
        case 'f':
            groups_count_fill = atoi(optarg);
            break;
        case 'p':
            groups_count_probe = atoi(optarg);
            break;
        case 'b':
            buffer_size = atoi(optarg);
            break;
        case 'B':
            block_size = atoi(optarg);
            break;
        case 'q':
            queue_size = atoi(optarg);
            break;
        case 'r':
            prefill_buffer_size = atoi(optarg);
            break;
        case 's':
            bsend_amount = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Unknown option\n");
            exit(1);
        }
    }
    if (groups_count_fill == 0 || groups_count_probe == 0 || buffer_size == 0 || prefill_buffer_size == 0 || bsend_amount == 0 || queue_size == 0 || block_size == 0) {
        usage(argv);
        exit(1);
    }
}

int main(int argc, char** argv) {
    process_command_line_options(argc, argv);
    FILE* constants = fopen("constants.h", "w");
    fprintf(constants, "#pragma once\n\n");
    fprintf(constants, "/* number of MPI groups in fill part */\n");
    fprintf(constants, "#define GROUPS_COUNT_FILL %d\n\n", groups_count_fill);
    fprintf(constants, "/* number of MPI groups in probe part */\n");
    fprintf(constants, "#define GROUPS_COUNT_PROBE %d\n\n", groups_count_probe);
    fprintf(constants, "/* size of block to process data */\n");
    fprintf(constants, "#define BLOCK_SIZE %d\n\n", block_size);
    fprintf(constants, "/* buffer size before insterting into the numa queue */\n");
    fprintf(constants, "#define BUFFER_SIZE %d\n\n", buffer_size);
    fprintf(constants, "/* buffer size before insterting into the numa queue */\n");
    fprintf(constants, "#define PREFILL_BUFFER_SIZE %d\n\n", prefill_buffer_size);
    fprintf(constants, "/* number of buffers that can be sent without waiting for completion */\n");
    fprintf(constants, "#define BSEND_AMOUNT %d\n", bsend_amount);
    fprintf(constants, "/* number of elements in the NUMA communication queues*/\n");
    fprintf(constants, "#define QUEUE_SIZE %d\n", queue_size);
    fclose(constants);
}