#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define dT 0.2f
#define G 0.6f
#define BLOCK_SIZE 512

// Global variables
int num_planets;
int num_timesteps;

// Host arrays
float2* velocities;
float4* planets;

// Device arrays
float2* velocities_d;
float4* planets_d;


double walltime() {
    static struct timeval t;
    gettimeofday(&t, NULL);
    return (t.tv_sec + 1e-6 * t.tv_usec);
}


// Parse command line arguments
void parse_args(int argc, char** argv){
    if(argc != 2){
        printf("Useage: nbody num_timesteps\n");
        exit(-1);
    }

    num_timesteps = strtol(argv[1], 0, 10);
}

// Reads planets from planets.txt
void read_planets(){

  //char* a;
  FILE* file = fopen("planets256.txt", "r");
  if(file == NULL){
      printf("'planets.txt' not found. Exiting\n");
      exit(-1);
  }

  char line[200];
  fgets(line, 200, file);
  sscanf(line, "%d", &num_planets);

    planets = (float4*)malloc(sizeof(float4)*num_planets);
    velocities = (float2*)malloc(sizeof(float2)*num_planets);

    for(int p = 0; p < num_planets; p++){
        fgets(line, 200, file);
        sscanf(line, "%f %f %f %f %f",
                &planets[p].x,
                &planets[p].y,
                &velocities[p].x,
                &velocities[p].y,
                &planets[p].z);
    }

    fclose(file);
}

// Writes planets to file
void write_planets(int timestep){
    char name[20];
    int n = sprintf(name, "gpu_planets_out.txt");

    FILE* file = fopen(name, "wr+");

    for(int p = 0; p < num_planets; p++){
        fprintf(file, "%f %f %f %f %f\n",
                planets[p].x,
                planets[p].y,
                velocities[p].x,
                velocities[p].y,
                planets[p].z);
    }

    fclose(file);
}

// TODO 7. Calculate the change in velocity for p, caused by the interaction with q
__device__ float2 calculate_velocity_change_planet(float4 p, float4 q){

  float2 dist;
  float2 dvel;
  dist.x = q.x - p.x;
  dist.y = q.y - p.y;
  if(dist.x == 0 && dist.y == 0){
  float2 vel = {0.0f, 0.0f};
  return vel;
  }
  float abs_dist= sqrt(dist.x*dist.x + dist.y*dist.y);
  float dist_cubed = abs_dist*abs_dist*abs_dist;
  // Calculate change in velocity
  dvel.x = dT*G*q.z/dist_cubed * dist.x;
  dvel.y = dT*G*q.z/dist_cubed * dist.y;
  return dvel;
}

// TODO 5. Calculate the change in velocity for my_planet, caused by the interactions with a block of planets
__device__ float2 calculate_velocity_change_block(float4 my_planet, float4* shared_planets){
  float2 velocity;
  velocity.x = 0.0f;
  velocity.y = 0.0f;
  for(int i = 0; i < blockDim.x; i++){
    float2 temp_vel = calculate_velocity_change_planet(my_planet, shared_planets[i]);
    velocity.x += temp_vel.x;
    velocity.y += temp_vel.y;
  }
  return velocity;
}

// TODO 4. Update the velocities by calculating the planet interactions
__global__ void update_velocities(float4* planets, float2* velocities, int num_planets){

  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  float4 my_planet = planets[thread_id];
  // Shared memory for this block
  __shared__ float4 shared_planets[BLOCK_SIZE];

  // Compute the velocity change for planets for one block at a time
  for(int i = 0; i < num_planets; i+=blockDim.x){
    shared_planets[threadIdx.x] = planets[i + threadIdx.x];
    __syncthreads();
    float2 tempv = calculate_velocity_change_block(my_planet, shared_planets);
    velocities[thread_id].x += tempv.x;
    velocities[thread_id].y += tempv.y;
    __syncthreads();
  }
}

// TODO 7. Update the positions of the planets using the new velocities
__global__ void update_positions(float4* planets, float2* velocities, int num_planets){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  planets[tid].x += velocities[tid].x * dT;
  planets[tid].y += velocities[tid].y * dT;
}


int main(int argc, char** argv){

    parse_args(argc, argv);
    read_planets();

    double calculation_time = 0;
    double memcopy_time = 0;



    // TODO 1. Allocate device memory, and transfer data to device

    cudaMalloc(&planets_d,    sizeof(float4)*num_planets);
    cudaMalloc(&velocities_d, sizeof(float2)*num_planets);

    double mem_start = walltime();
    cudaMemcpy(planets_d,    planets,     sizeof(float4)*num_planets,
                                          cudaMemcpyHostToDevice);
    cudaMemcpy(velocities_d, velocities,  sizeof(float2)*num_planets,
                                          cudaMemcpyHostToDevice);
    memcopy_time += walltime() - mem_start;


    // Calculating the number of blocks
    int num_blocks = num_planets/BLOCK_SIZE + ((num_planets%BLOCK_SIZE == 0) ? 0 : 1);


    // Main loop

    double calc_start = walltime();
    for(int t = 0; t < num_timesteps; t++){
        // TODO 2. Call kernels
        //Update velocities
        update_velocities<<<num_blocks, BLOCK_SIZE>>>(planets_d, velocities_d,
                                                      num_planets);
        //Update positions
        update_positions<<<num_blocks, BLOCK_SIZE>>>(planets_d, velocities_d,
                                                    num_planets);
    }
    calculation_time = walltime() - calc_start;



    double mem_start = walltime();
    // TODO 3. Transfer data back to host
    cudaMemcpy(planets,     planets_d,    sizeof(float4)*num_planets,
                                          cudaMemcpyDeviceToHost);
    cudaMemcpy(velocities,  velocities_d, sizeof(float2)*num_planets,
                                          cudaMemcpyDeviceToHost);
    memcopy_time += walltime() - mem_start;


    cudaFree(planets_d);
    cudaFree(velocities_d);
    // Output
    write_planets(num_timesteps);
    printf("%7.7f ms\n", calculation_time);
    printf("%7.7f ms\n", memcopy_time);



}
