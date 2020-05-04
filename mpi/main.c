/*
 *  CSE 5441 Autumn 2019
 *  Lab #1 Solution
 *  Mike Gargasz.13
 *  Compile with the provided makefile
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#define NS_PER_US 1000.0
#define EPSILON 0.1
#define AFFECT_RATE 0.1
#define NUM_THREADS 28

struct box {
    int id;
    int upper_left_y;
    int upper_left_x;
    int height;
    int width;
    int top_neighbors;
    int top_neighbors_list[10];
    int bottom_neighbors;
    int bottom_neighbors_list[10];
    int left_neighbors;
    int left_neighbors_list[10];
    int right_neighbors;
    int right_neighbors_list[10];
};

int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

int min(int num1, int num2)
{
    return (num1 > num2 ) ? num2 : num1;
}

void populate_grid(struct box boxes[], double current_dsv[], char *line, size_t size, FILE *input) {
    int count = 0;
    char *n;
    int bytes_read = 0;
    bytes_read = getline(&line, &size, input);
    //filter out blanks
    while(bytes_read <= 1 || line[0] == '\t') {
        bytes_read = getline(&line, &size, input);
    }
    int next_id = atoi(strtok(line, " \t"));
    while(next_id != -1){
        boxes[count].id = next_id;

        //coords
        getline(&line, &size, input);
        boxes[count].upper_left_y = atoi(strtok(line, " \t"));
        boxes[count].upper_left_x = atoi(strtok(NULL, " \t"));
        boxes[count].height = atoi(strtok(NULL, " \t"));
        boxes[count].width = atoi(strtok(NULL, " \t"));

        //top
        getline(&line, &size, input);
        boxes[count].top_neighbors = atoi(strtok(line, " \t"));
        if(boxes[count].top_neighbors != 0) {
            int i = 0;
            n = strtok(NULL, " \t");
            do {
                boxes[count].top_neighbors_list[i] = atoi(n);
                i++;
                n = strtok(NULL, " \t");
            } while (n);
        }
        //bottom
        getline(&line, &size, input);
        boxes[count].bottom_neighbors = atoi(strtok(line, " \t"));
        if(boxes[count].bottom_neighbors != 0) {
            int i = 0;
            n = strtok(NULL, " \t");
            do {
                boxes[count].bottom_neighbors_list[i] = atoi(n);
                i++;
                n = strtok(NULL, " \t");
            } while (n);
        }
        //left
        getline(&line, &size, input);
        boxes[count].left_neighbors = atoi(strtok(line, " \t"));
        if(boxes[count].left_neighbors != 0) {
            int i = 0;
            n = strtok(NULL, " \t");
            do {
                boxes[count].left_neighbors_list[i] = atoi(n);
                i++;
                n = strtok(NULL, " \t");
            } while (n);
        }
        //right
        getline(&line, &size, input);
        boxes[count].right_neighbors = atoi(strtok(line, " \t"));
        if(boxes[count].right_neighbors != 0) {
            int i = 0;
            n = strtok(NULL, " \t");
            do {
                boxes[count].right_neighbors_list[i] = atoi(n);
                i++;
                n = strtok(NULL, " \t");
            } while (n);
        }
        //temp
        getline(&line, &size, input);
        current_dsv[count] = atoi(strtok(line, " \t"));

        count++;
        bytes_read = getline(&line, &size, input);
        //filter out blanks
        while(bytes_read <= 1 || line[0] == '\t') {
            bytes_read = getline(&line, &size, input);
        }
        next_id = atoi(strtok(line, " \t"));
    }
}

double process_top_neighbors(int num_containers, struct box boxes[], double current_dsv[], int rank, double partial_new_temp[], int partial_perimeter[]) {
        omp_set_dynamic(0);
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for
        for(int index = 0;index < num_containers; index++) {
            int begin_overlap = 0;
            int end_overlap = 0;
            int overlap = 0;
            int perimeter = 0;
            double new_temp = 0;
            if (boxes[index].top_neighbors > 0) {
                for (int i = 0; i < boxes[index].top_neighbors; i++) {
                    begin_overlap = max(boxes[index].upper_left_x, boxes[boxes[index].top_neighbors_list[i]].upper_left_x);
                    end_overlap = min(boxes[index].upper_left_x + boxes[index].width,
                                      boxes[boxes[index].top_neighbors_list[i]].upper_left_x +
                                      boxes[boxes[index].top_neighbors_list[i]].width);
                    overlap = end_overlap - begin_overlap;
                    new_temp += (overlap * current_dsv[boxes[index].top_neighbors_list[i]]);
                    perimeter += overlap;
                }
            } else {
                partial_new_temp[index] += (boxes[index].width * current_dsv[index]);
                partial_perimeter[index] += boxes[index].width;
            }
        }
}

double process_bottom_neighbors(int num_containers, struct box boxes[], double current_dsv[], int rank, double partial_new_temp[], int partial_perimeter[]) {
        omp_set_dynamic(0);
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for
        for(int index = 0;index < num_containers; index++) {
            int begin_overlap = 0;
            int end_overlap = 0;
            int overlap = 0;
            int perimeter = 0;
            double new_temp = 0;
            if(boxes[index].bottom_neighbors > 0) {
                for(int i = 0; i < boxes[index].bottom_neighbors; i++) {
                    begin_overlap = max(boxes[index].upper_left_x, boxes[boxes[index].bottom_neighbors_list[i]].upper_left_x);
                    end_overlap = min(boxes[index].upper_left_x + boxes[index].width, boxes[boxes[index].bottom_neighbors_list[i]].upper_left_x + boxes[boxes[index].bottom_neighbors_list[i]].width);
                    overlap = end_overlap - begin_overlap;
                    new_temp += (overlap * current_dsv[boxes[index].bottom_neighbors_list[i]]);
                    perimeter += overlap;
                }
            } else {
                partial_new_temp[index] += (boxes[index].width * current_dsv[index]);
                partial_perimeter[index] += boxes[index].width;
            }
        }
}

double process_left_neighbors(int num_containers, struct box boxes[], double current_dsv[], int rank, double partial_new_temp[], int partial_perimeter[]) {
        omp_set_dynamic(0);
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for
        for(int index = 0;index < num_containers; index++) {
            int begin_overlap = 0;
            int end_overlap = 0;
            int overlap = 0;
            int perimeter = 0;
            double new_temp = 0;
            if(boxes[index].left_neighbors > 0) {
                for(int i = 0; i < boxes[index].left_neighbors; i++) {
                    begin_overlap = max(boxes[index].upper_left_y, boxes[boxes[index].left_neighbors_list[i]].upper_left_y);
                    end_overlap = min(boxes[index].upper_left_y + boxes[index].height, boxes[boxes[index].left_neighbors_list[i]].upper_left_y + boxes[boxes[index].left_neighbors_list[i]].height);
                    overlap = end_overlap - begin_overlap;
                    new_temp += (overlap * current_dsv[boxes[index].left_neighbors_list[i]]);
                    perimeter += overlap;
                }
            } else {
                partial_new_temp[index] += (boxes[index].height * current_dsv[index]);
                partial_perimeter[index] += boxes[index].height;
            }
        }
}

double process_right_neighbors(int num_containers, struct box boxes[], double current_dsv[], int rank, double partial_new_temp[], int partial_perimeter[]) {
        omp_set_dynamic(0);
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for
        for(int index = 0;index < num_containers; index++) {
            int begin_overlap = 0;
            int end_overlap = 0;
            int overlap = 0;
            int perimeter = 0;
            double new_temp = 0;
            if(boxes[index].right_neighbors > 0) {
                for(int i = 0; i < boxes[index].right_neighbors; i++) {
                    begin_overlap = max(boxes[index].upper_left_y, boxes[boxes[index].right_neighbors_list[i]].upper_left_y);
                    end_overlap = min(boxes[index].upper_left_y + boxes[index].height, boxes[boxes[index].right_neighbors_list[i]].upper_left_y + boxes[boxes[index].right_neighbors_list[i]].height);
                    overlap = end_overlap - begin_overlap;
                    new_temp += (overlap * current_dsv[boxes[index].right_neighbors_list[i]]);
                    perimeter += overlap;
                }
            } else {
                partial_new_temp[index] += (boxes[index].height * current_dsv[index]);
                partial_perimeter[index] += boxes[index].height;
            }
        }
}

void amr_dissipation(int num_containers, struct box boxes[], double current_dsv[], double results[], double difference, double temp_dsv[], double epsilon, double affect_rate, int rank) {
    double partial_new_temp[num_containers];
    int partial_perimeter[num_containers];
    double new_temp[num_containers];
    int perimeter[num_containers];

    while (difference > epsilon) {
        if (rank == 0) {

            //send current DSVs
            MPI_Send(current_dsv, num_containers, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Send(current_dsv, num_containers, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
            MPI_Send(current_dsv, num_containers, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);
            MPI_Send(current_dsv, num_containers, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD);

        } else if (rank == 1) {

            MPI_Recv(current_dsv, num_containers, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            process_top_neighbors(num_containers, boxes, current_dsv, rank, partial_new_temp, partial_perimeter);

        } else if (rank == 2) {

            MPI_Recv(current_dsv, num_containers, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            process_top_neighbors(num_containers, boxes, current_dsv, rank, partial_new_temp, partial_perimeter);

        } else if (rank == 3) {

            MPI_Recv(current_dsv, num_containers, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            process_top_neighbors(num_containers, boxes, current_dsv, rank, partial_new_temp, partial_perimeter);

        } else if (rank == 4) {

            MPI_Recv(current_dsv, num_containers, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            process_top_neighbors(num_containers, boxes, current_dsv, rank, partial_new_temp, partial_perimeter);

        }
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(partial_new_temp, new_temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(partial_perimeter, perimeter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if(rank == 0) {
            omp_set_dynamic(0);
            omp_set_num_threads(NUM_THREADS);
            #pragma omp parallel for
            for (int index = 0; index < num_containers; index++) {
                //new dsv
                new_temp[index] /= perimeter[index];
                if (new_temp[index] > current_dsv[index]) {
                    temp_dsv[index] = current_dsv[index] + (new_temp[index] - current_dsv[index]) * affect_rate;
                } else {
                    temp_dsv[index] = current_dsv[index] - (current_dsv[index] - new_temp[index]) * affect_rate;
                }
            }
            //End Parallel Region
            results[1] = current_dsv[0];
            results[2] = current_dsv[0];
            for (int container = 0; container < num_containers; container++) {
                current_dsv[container] = temp_dsv[container];
                if (current_dsv[container] > results[1]) {
                    results[1] = current_dsv[container];
                }
                if (current_dsv[container] < results[2]) {
                    results[2] = current_dsv[container];
                }
            }
            results[0]++;
            difference = (results[1] - results[2]) / results[1];
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    char *line = NULL;
    size_t size;
    int num_containers;
    FILE * input;
    int rank, num_processes;

    //read in first line
    getline(&line, &size, input);
    char *n = strtok(line, " \t");
    num_containers = atoi(n);
    //declare array of boxes
    struct box boxes[num_containers];
    double current_dsv[num_containers];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) {
        input = fopen("../test/testgrid_50_78", "r");
        if (input == NULL) {
            printf("FILE DOES NOT EXIST \n");
            exit(EXIT_FAILURE);
        }
        populate_grid(boxes, current_dsv, line, size, input);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //perform AMR Dissipation
    double results[3] = {0, 0, 0};
    struct timespec start_realtime, end_realtime;
    double chrono;
    double epsilon = EPSILON;
    double affect_rate = AFFECT_RATE;
    double difference = epsilon + 1;
    double temp_dsv[num_containers];
    clock_t clock_mes;
    time_t start_timer, end_timer, time_mes;

    clock_mes = clock();
    time(&start_timer);
    clock_gettime(CLOCK_REALTIME, &start_realtime);
    {
        amr_dissipation(num_containers, boxes, current_dsv, results, difference, temp_dsv, epsilon, affect_rate, rank);
    }
    clock_gettime(CLOCK_REALTIME, &end_realtime);
    time(&end_timer);
    time_mes = end_timer - start_timer;
    clock_mes = clock() - clock_mes;
    chrono = (double)(((end_realtime.tv_sec - start_realtime.tv_sec) * CLOCKS_PER_SEC) + ((end_realtime.tv_nsec - start_realtime.tv_nsec) / NS_PER_US));

    //print results
    printf("*******************************************************************\n");
    printf("dissipation converged in %d iterations\n", (int)results[0]);
    printf("    with max DSV = %lf and min DSV = %lf\n", results[1], results[2]);
    printf("    affect rate = %0.2f; \t epsilon = %0.2f;\n", affect_rate, epsilon);
    printf("\n");
    printf("elapsed convergence loop time  (clock): %lu\n", clock_mes);
    printf("elapsed convergence loop time   (time): %lu\n", time_mes);
    printf("elapsed convergence loop time (chrono): %0.1lf\n", chrono);
    printf("*******************************************************************\n");

    fclose(input);
}