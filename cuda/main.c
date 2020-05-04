/*
 *  CSE 5441 Autumn 2019
 *  Lab #4 Solution
 *  Mike Gargasz.13
 *  Compile with the provided makefile
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define NS_PER_US 1000.0
#define MEM_PER_BLOCK 49152

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
    double dsv;
};

int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}
int min(int num1, int num2)
{
    return (num1 > num2 ) ? num2 : num1;
}
void populate_grid(struct box boxes[], char *line, size_t size) {
    int count = 0;
    char *n;
    int bytes_read = 0;
    bytes_read = getline(&line, &size, stdin);
    //filter out blanks
    while(bytes_read <= 1 || line[0] == '\t') {
        bytes_read = getline(&line, &size, stdin);
    }
    int next_id = atoi(strtok(line, " \t"));
    while(next_id != -1){
        boxes[count].id = next_id;

        //coords
        getline(&line, &size, stdin);
        boxes[count].upper_left_y = atoi(strtok(line, " \t"));
        boxes[count].upper_left_x = atoi(strtok(NULL, " \t"));
        boxes[count].height = atoi(strtok(NULL, " \t"));
        boxes[count].width = atoi(strtok(NULL, " \t"));

        //top
        getline(&line, &size, stdin);
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
        getline(&line, &size, stdin);
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
        getline(&line, &size, stdin);
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
        getline(&line, &size, stdin);
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
        getline(&line, &size, stdin);
        boxes[count].dsv = atoi(strtok(line, " \t"));

        count++;
        bytes_read = getline(&line, &size, stdin);
        //filter out blanks
        while(bytes_read <= 1 || line[0] == '\t') {
            bytes_read = getline(&line, &size, stdin);
        }
        next_id = atoi(strtok(line, " \t"));
    }
}

void amr_dissipation(int num_containers, struct box boxes[], double results[], double difference, float epsilon, float affect_rate) {
    double temp_dsv[num_containers];
    int tpb = MEM_PER_BLOCK / sizeof(boxes[0]);
    int nBlocks = num_containers / tpb + 1;
    int *gpu_input;
    int *gpu_output;
    int *gpu_affect_rate;
    size_t inputMemsize = nBlocks * tpb * sizeof(boxes[0]);
    size_t outputMemsize = num_containers * sizeof(double);
    cudaMaclloc((void**)&gpu_input, inputMemsize);
    cudaMalloc((void**)&gpu_output, outputMemsize);
    cudaMalloc((void**)&gpu_affect_rate, sizeof(float));
    cudaMemcpy(gpu_input, boxes, inputMemsize, cudaMemcpyHostToDevice);
    dim3 dimGrid(nBlocks);
    dim3 dimBlock(tpb);

    while(difference > epsilon) {

        processBoxes<<< dimGrid, dimBlock >>>(gpu_input, gpu_output, gpu_affect_rate);
        cudaMemcpy(temp_dsv, gpu_output, outputMemsize, cudaMemcpyDeviceToHost);

        results[1] = boxes[0].dsv;
        results[2] = boxes[0].dsv;
        for(int container = 0;container < num_containers; container++) {
            boxes[container].dsv = temp_dsv[container];
            if(boxes[container].dsv > results[1]) {
                results[1] = boxes[container].dsv;
            }
            if(boxes[container].dsv < results[2]) {
                results[2] = boxes[container].dsv;
            }
        }
        results[0]++;
        difference = (results[1] - results[2]) / results[1];
    }
}

__global__ void processBoxes( int  *boxes, int *temp_dsv, float *affect_rate){
    printf("I RAN");
    int index = blockDim.x * blockidx.x + threadIdx.x;
    int begin_overlap = 0;
    int end_overlap = 0;
    int overlap = 0;
    int perimeter = 0;
    double new_temp = 0;
    //top
    if(boxes[index].top_neighbors > 0) {
        for(int i = 0; i < boxes[index].top_neighbors; i++) {
            begin_overlap = max(boxes[index].upper_left_x, boxes[boxes[index].top_neighbors_list[i]].upper_left_x);
            end_overlap = min(boxes[index].upper_left_x + boxes[index].width, boxes[boxes[index].top_neighbors_list[i]].upper_left_x + boxes[boxes[index].top_neighbors_list[i]].width);
            overlap = end_overlap - begin_overlap;
            new_temp += (overlap * boxes[boxes[index].top_neighbors_list[i]].dsv);
            perimeter += overlap;
        }
    } else {
        new_temp += (boxes[index].width * boxes[index].dsv);
        perimeter += boxes[index].width;
    }
    //bottom
    if(boxes[index].bottom_neighbors > 0) {
        for(int i = 0; i < boxes[index].bottom_neighbors; i++) {
            begin_overlap = max(boxes[index].upper_left_x, boxes[boxes[index].bottom_neighbors_list[i]].upper_left_x);
            end_overlap = min(boxes[index].upper_left_x + boxes[index].width, boxes[boxes[index].bottom_neighbors_list[i]].upper_left_x + boxes[boxes[index].bottom_neighbors_list[i]].width);
            overlap = end_overlap - begin_overlap;
            new_temp += (overlap * boxes[boxes[index].bottom_neighbors_list[i]].dsv);
            perimeter += overlap;
        }
    } else {
        new_temp += (boxes[index].width * boxes[index].dsv);
        perimeter += boxes[index].width;
    }
    //left
    if(boxes[index].left_neighbors > 0) {
        for(int i = 0; i < boxes[index].left_neighbors; i++) {
            begin_overlap = max(boxes[index].upper_left_y, boxes[boxes[index].left_neighbors_list[i]].upper_left_y);
            end_overlap = min(boxes[index].upper_left_y + boxes[index].height, boxes[boxes[index].left_neighbors_list[i]].upper_left_y + boxes[boxes[index].left_neighbors_list[i]].height);
            overlap = end_overlap - begin_overlap;
            new_temp += (overlap * boxes[boxes[index].left_neighbors_list[i]].dsv);
            perimeter += overlap;
        }
    } else {
        new_temp += (boxes[index].height * boxes[index].dsv);
        perimeter += boxes[index].height;
    }
    //right
    if(boxes[index].right_neighbors > 0) {
        for(int i = 0; i < boxes[index].right_neighbors; i++) {
            begin_overlap = max(boxes[index].upper_left_y, boxes[boxes[index].right_neighbors_list[i]].upper_left_y);
            end_overlap = min(boxes[index].upper_left_y + boxes[index].height, boxes[boxes[index].right_neighbors_list[i]].upper_left_y + boxes[boxes[index].right_neighbors_list[i]].height);
            overlap = end_overlap - begin_overlap;
            new_temp += (overlap * boxes[boxes[index].right_neighbors_list[i]].dsv);
            perimeter += overlap;
        }
    } else {
        new_temp += (boxes[index].height * boxes[index].dsv);
        perimeter += boxes[index].height;
    }
    //new dsv
    new_temp /= perimeter;
    if(new_temp > boxes[index].dsv) {
        temp_dsv[index] = boxes[index].dsv + (new_temp - boxes[index].dsv) * &affect_rate;
    } else {
        temp_dsv[index] = boxes[index].dsv - (boxes[index].dsv - new_temp) * &affect_rate;
    }
}

int main(int argc, char *argv[]) {
    char *line = NULL;
    size_t size;
    int num_containers;

    //read in first line
    getline(&line, &size, stdin);
    char *n = strtok(line, " \t");
    num_containers = atoi(n);

    //declare array of boxes
    struct box boxes[num_containers];
    populate_grid(boxes, line, size);

    //perform AMR Dissipation
    double results[3] = {0};
    struct timespec start_realtime, end_realtime;
    double chrono;
    float epsilon = atof(argv[1]);
    float affect_rate = atof(argv[2]);
    double difference = epsilon + 1;
    clock_t clock_mes;
    time_t start_timer, end_timer, time_mes;

    clock_mes = clock();
    time(&start_timer);
    clock_gettime(CLOCK_REALTIME, &start_realtime);
    {
        amr_dissipation(num_containers, boxes, results, difference, epsilon, affect_rate);
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


}