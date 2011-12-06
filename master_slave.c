/*
 * master_slave.c
 *
 * Parallel Spatial Simulated Annealing using OpenMPI, and a master-slave
 * parallelism approach.
 *
 * Author(s): Guy Cobb <guy.cobb@gmail.com>
 *            Caleb Phillips <cphillips@smallwhitecube.com>
 *
 * License: Academic/Beerware/MIT or something equivalent
 */

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>

#define N 10          // number of coordinates in a solution
#define M 50          // number of solutions to maintain in pool
#define MAX_RUNS 100  // maximum number of slave-runs to do
#define VERBOSE 1     // debugging output (on = 1, off = 0)
#define COLOR 1       // colorize output (on = 1, off = 0)
#define GDB_ATTACH 0  // sleep in master to allow a gdb attachment

#define TAG_COORDS 0
#define TAG_SHUTDOWN 1
#define TAG_FITBEFORE 2
#define TAG_FITAFTER 3

/* some colorization stuff stolen from http://linuxgazette.net/issue65/padala.html */
#define RESET		0
#define BRIGHT 		1
#define DIM		2
#define UNDERLINE 	3
#define BLINK		4
#define REVERSE		7
#define HIDDEN		8

#define BLACK 		0
#define RED		1
#define GREEN		2
#define YELLOW		3
#define BLUE		4
#define MAGENTA		5
#define CYAN		6
#define	WHITE		7

struct coord{
	double n;  // northing
	double e;  // easting (EPSG:32160)
	int i;     // index into candidates array
}; 

struct job{
	struct coord *coords;
	int checked_out;
	int run_count;
	int owner;
	int last_owner;
	double fitness;
	double last_gain;
};

// zOMG GLOBALS!
int rank;                    // rank of process
int np;                      // number of processes (1 master, np-1 slaves)
MPI_Datatype coord_type; 
struct coord coords[N];
struct job pool[M];
struct coord *points;
double fitness_before;
double fitness_after;

/* colorize output text */
void textcolor(int attr, int fg, int bg) {
  if(!COLOR) return;
  char command[64];
  /* Command is the control command to the terminal */
  sprintf(command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
  printf("%s", command);
}

/* debugging output */ 
void debug(const char* fmt, ...) {
  if(!VERBOSE) return;
  va_list args;
  char buf[BUFSIZ];
  static int inc  = 0;
  static char hostname[BUFSIZ] = "";

  /* makes the output really long... */
  /*
  if (strlen(hostname) == 0) {
    gethostname(hostname, BUFSIZ);	
  }
  */

  sprintf(buf, "DEBUG: [ %d / %d  %d ]:\t%s", rank, np, inc++, fmt);

  va_start (args, fmt);
  textcolor(DIM,YELLOW,BLACK);
  vprintf (buf, args);
  textcolor(RESET, WHITE, BLACK);
  va_end (args);
  fflush(stdout);
}

/* initialize the list of all possible points (master) */
int initialize_points(void){
  const int buflen = 512;
  char buffer[buflen];
  char *strptr, *str, *tok;
  int j,i;
  FILE *fh;
  int n = 0;

  fh = fopen("candidates.txt","r");
  while(fgets(buffer,buflen,fh) != NULL) n++;
  fclose(fh);

  debug("%d candidates, allocating memory...\n",n);

  points = (struct coord *)malloc(sizeof(struct coord)*n);
  if(points == NULL){
    debug("Failed to allocate memory to hold candidate points. Exiting.");
    exit(1);
  }

  // 6 273 1 2017953.8167517 4581541.48079023 FALSE
  i = 0;
  fh = fopen("candidates.txt","r");
  while(fgets(buffer,buflen,fh) != NULL){
    for(j = 1, str = buffer; ;j++, str = NULL){
      tok = strtok_r(str," ",&strptr);
      if(tok == NULL) break;
      if(j == 4) points[i].e = atof(tok);
      else if(j == 5) points[i].n = atof(tok);
      else if(j == 1){ 
        points[i].i = atoi(tok);
      }
    }
    i++;
  }
  fclose(fh);
  return n;
}

/* fill pool with M randomly chosen (but not repeating) candidate solutions (master) */
void initialize_pool(int npoints){
  int i, j, k, r;
  for(i = 0; i < M; i++){
    pool[i].checked_out = 0;
    pool[i].run_count = 0;
    pool[i].owner = -1;
    pool[i].last_owner = -1;
    pool[i].fitness = -1.0;
    pool[i].last_gain = 0.0;
    pool[i].coords = (struct coord *)malloc(sizeof(struct coord)*N);
    if(pool[i].coords == NULL){
      debug("Failed to allocate memory for pool slot %d\n",i);
      exit(1);
    }
    for(j = 0; j < N; j++){
      r = -1;
      while(r < 0){
        r = floor((((double)rand())/((double) RAND_MAX)) * (npoints-1));
        for(k = 0; k < j; k++){
          if(pool[i].coords[k].i == r){
            r = -1;
            break;
          }
        }
      }
      pool[i].coords[j].e = points[r].e;
      pool[i].coords[j].n = points[r].n;
      pool[i].coords[j].i = r;
    }
  }
}

void print_coords(struct coord *c, int n){
  int i;
  for(i = 0; i < N; i++){ debug("%d: (%f,%f)\n",c[i].i,c[i].e,c[i].n); }
}

/* Prepare a candidate to send to a slave (master) */
void check_out_candidate(int id){
  int r,j;
  int found = -1;
  // if the pool has no non-taken candidates, this will loop forever (or until a client checks one in)
  // because of this: M should be >> np
  while(found < 0){
    r = floor((((double)rand())/((double) RAND_MAX)) * M);
    if(pool[r].checked_out == 0){
      found = r;
    }
  }    
  pool[found].checked_out = 1;
  pool[found].run_count++;
  pool[found].owner = id;
 
  for(j = 0; j < N; j++){
    coords[j].n = pool[found].coords[j].n;
    coords[j].e = pool[found].coords[j].e;
    coords[j].i = pool[found].coords[j].i;
  }

  debug("Checked out Candidate to %d\n",id);
  print_coords(coords,N);
}

/* Process a candidate returned by a slave (master) */
void check_in_candidate(int id){
  int i,j;
  int found = 0;
  for(i = 0; i < M; i++){
    if(pool[i].owner == id){
      pool[i].last_owner = id;
      pool[i].owner = -1;
      pool[i].checked_out = 0;
      if(fitness_after > 0){ // fitness_after is <0 if something went wrong on the slave
        pool[i].fitness = fitness_after;
        pool[i].last_gain = fitness_before-fitness_after;
        for(j = 0; j < N; j++){
          pool[i].coords[j].i = coords[j].i;
          pool[i].coords[j].n = points[coords[j].i].n;
          pool[i].coords[j].e = points[coords[j].i].e;
        }
      }
      found = 1;
      break;
    }
  }
  debug("Checked in candidate from %d with fitness %f (gain %f)\n",id,pool[i].fitness,pool[i].last_gain);
  print_coords(coords,N);

  if(found == 0){
    debug("Failed to find item to check into pool for slave %d\n",id);
  }
}

/* Do something with a candidate (slave) */
int process_candidate(double *fitness_before, double *fitness_after, int nruns){
  int i, len, j, maxlen;
  const char *prog = "/home/esl/caleb/R-2.13.0/bin/R --vanilla --args ";
  const char *progpost = " < slave.R 2>&1";
  char *cmd;
  char *pos;
  char *str;
  char *tmp;
  FILE *ph;
  char *tok;
  char *strptr;
  const unsigned int buflen = 512;
  char buffer[buflen];
  double wpe_before, wpe_after, akv_before, akv_after;
  int indices[N];

  // construct the run id which is <rank>.<#>
  int runidlen = 8;
  if(nruns <= 1) runidlen += 1;
  else runidlen += ceil(log((double)nruns)/log(10.0));
  char *runid = (char *)malloc(sizeof(char)*runidlen);
  snprintf(runid,runidlen,"%06d.%d",rank,nruns);

  debug("Processing Candidate: \n"); 
  print_coords(coords,N);

  // determine the maximum printed length of any given id
  maxlen = 0;
  for(i = 0; i < N; i++){
    len = ceil(log((double)(coords[i].i))/log(10.0));
    if(len > maxlen) maxlen = len;
  }
  // allocate enough space to whole entire command
  // <prog> <num 1>.<num 2>.[...].<num n> <runid>
  cmd = (char *)malloc(sizeof(char)*(strlen(prog) + strlen(progpost) + strlen(runid) + 1 + N*(maxlen+1) + 1));
  tmp = (char *)malloc(sizeof(char)*(maxlen+2));

  // construct the command string
  strcpy(cmd,prog);
  pos = &(cmd[strlen(cmd)]); // points to the null byte at the end of the cmd string
  for(i = 0; i < N; i++){
    sprintf(tmp,"%d.",coords[i].i);
    strcpy(pos,tmp);
    pos = &(cmd[strlen(cmd)]);
  }
  cmd[strlen(cmd)-1] = '\0'; // remove trailing period
  strcat(cmd," ");
  strcat(cmd,runid);
  strcat(cmd,progpost);
  debug("%s\n",cmd);

  // FITNESS 3.959667 3.958261 5.579055 5.574111 NULL
  // > print(cat("SAMPLE",e,""))
  // SAMPLE 7 1 2 3 9 4 6 8 10 37725 NULL
  wpe_before = -1.0;
  wpe_after = -1.0;
  indices[N-1] = -1;
  ph = popen(cmd,"r");
  while(fgets(buffer,buflen,ph) != NULL){
    debug(buffer);
    if(strncmp(buffer,"FITNESS",7) == 0){
      for(j = 1, str = buffer; ;j++, str = NULL){
        tok = strtok_r(str," ",&strptr);
        if(tok == NULL) break;

        if(j == 2) wpe_before = atof(tok);
        else if(j == 3) wpe_after = atof(tok);
        else if(j == 4) akv_before = atof(tok);
        else if(j == 5) akv_after = atof(tok);
      }
    }
    else if(strncmp(buffer,"SAMPLE",6) == 0){
      for(j = 1, str = buffer; ;j++, str = NULL){
        tok = strtok_r(str," ",&strptr);
        if(tok == NULL) break;
        if((j > 1) && (j <= (N+1))) indices[j-2] = atoi(tok);
      }
    }
  }
  pclose(ph);

  // sanity check that we actually read in something that looks like a result
  if((wpe_before < 0) || (wpe_after < 0) || (indices[N-1] < 0)){
    debug("Failed to parse output from R. wpe_before = %f, wpe_after = %f, last index = %d\n",wpe_before,wpe_after,indices[N-1]);
    return 1;
  }

  *fitness_before = wpe_before;
  *fitness_after = wpe_after;
  for(i = 0; i < N; i++){
    coords[i].i = indices[i];
  }

  return 0;
}

int main(int argc, char** argv) {
  int nruns = 0;
  int npoints = 0;
  int master = 0;
  int dummy_val = 0;
  int blocklens[4] = {1, 1, 1, 1}; 
  int i;
  int done = 0;
  int work_available = 0;
  int base;

  srand(42);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  MPI_Status status;
  MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_UB};
  MPI_Aint     disps[4]; 

  /* Calculate the displacements of the structure components and the structure overall*/
  MPI_Get_address(coords, disps);
  MPI_Get_address(&coords[0].e, disps+1);
  MPI_Get_address(&coords[0].i, disps+2);
  MPI_Get_address(&coords[1], disps+3);
  base = disps[0];
  for (i=0;i<4; i++) disps[i] -= base;

  MPI_Type_struct(4, blocklens, disps, types, &coord_type); 
  MPI_Type_commit(&coord_type); 

  debug("Node %d of %d reporting in\n", rank, np);

  if (rank == 0){ 	/* Master process creates work and sends to slaves */

    // gdb attachment hack ala http://www.open-mpi.org/faq/?category=debugging#serial-debuggers
    if(GDB_ATTACH){
      i = 0;
      debug("Sleeping until a GDB attach on PID %d\n",getpid());
      fflush(stdout);
      while(i == 0) sleep(5);
    }

    debug("Master is initializing points...\n");
    npoints = initialize_points();
    debug("Okay, points are initialized.\n");

    debug("Master is initializing pool...\n");
    initialize_pool(npoints);
    debug("Pool initialized!\n");

    /* Send initial candidates out to slaves */
    for (i=1; i<np; i++){
      check_out_candidate(i);
      MPI_Send(coords, N, coord_type, i, TAG_COORDS, MPI_COMM_WORLD);
    }

    /* Go through the rest of the candidate pool and send to slaves */
    for (;;){
      /* See if any of the slaves have a possible solution for us */
      MPI_Recv(coords, N, coord_type, MPI_ANY_SOURCE, TAG_COORDS, MPI_COMM_WORLD, &status);
      MPI_Recv(&fitness_before, 1, MPI_DOUBLE, status.MPI_SOURCE, TAG_FITBEFORE, MPI_COMM_WORLD, &status);
      MPI_Recv(&fitness_after, 1, MPI_DOUBLE, status.MPI_SOURCE, TAG_FITAFTER, MPI_COMM_WORLD, &status);

      nruns++; // count runs by all slaves

      // fitness_after will be negative if something went wrong
      check_in_candidate(status.MPI_SOURCE);

      if (nruns > MAX_RUNS) break;

      // FIXME: at some point, prune or evolve the pool or something

      /* Send out next candidate to the slave that just reported in */
      check_out_candidate(status.MPI_SOURCE);
      MPI_Send(coords, N, coord_type, status.MPI_SOURCE, TAG_COORDS, MPI_COMM_WORLD);
    }

    /* send shutdown signal to all slaves */
    for (i=1; i<np; i++) MPI_Send(&dummy_val, 1, MPI_INT, i, TAG_SHUTDOWN, MPI_COMM_WORLD);

    // look through pool and pick out winner
    double min_fitness = -1.0;
    double max_fitness = -1.0;
    int min_fitness_index = 0;
    for (i = 0; i < M; i++){
      if((min_fitness < 0) || ((pool[i].run_count > 0) && (pool[i].fitness < min_fitness))){
        min_fitness = pool[i].fitness;
        min_fitness_index = i;
      }
      if((max_fitness < 0) || ((pool[i].run_count > 0) && (pool[i].fitness > max_fitness))){
        max_fitness = pool[i].fitness;
      }
    }
    debug("Minimum fitness obtained = %f, worst (optimized) fitness is %f worse\n",min_fitness,max_fitness-min_fitness);
    debug("Best result is in sa_slave_%06d.*.RData, although I'm not sure which run id...\n",pool[min_fitness_index].last_owner);
    debug("Here's the sample:\n");
    print_coords(pool[min_fitness_index].coords,N);

  } else {

    /* Receive work from the master process until it says we're done */
    for (;;){

      /* Check for a shutdown signal from the master */
      MPI_Iprobe(master, TAG_SHUTDOWN, MPI_COMM_WORLD, &done, &status);
      if (done){ 
        MPI_Recv(&dummy_val, 1, MPI_INT, master, TAG_SHUTDOWN, MPI_COMM_WORLD, &status);
	break;
      }

      /* Check if work is available from the master, and if so receive and process it */
      MPI_Iprobe(master, TAG_COORDS, MPI_COMM_WORLD, &work_available, &status);
      if (work_available){
        /* Receive candidate from master to process */ 
	MPI_Recv(coords, N, coord_type, master, TAG_COORDS, MPI_COMM_WORLD, &status);

	/* Perform simulated annealing on candidate */
	if(process_candidate(&fitness_before,&fitness_after,nruns++) > 0){
          debug("Problems in process candidate\n");
          // make these negative so the master know's something went wrong.
          fitness_after = -1.0;
          fitness_before = -1.0;
        }

	/* Send our results back to the master */
        MPI_Send(&fitness_after, 1, MPI_DOUBLE, master, TAG_FITAFTER, MPI_COMM_WORLD);
        MPI_Send(&fitness_before, 1, MPI_DOUBLE, master, TAG_FITBEFORE, MPI_COMM_WORLD);
        MPI_Send(coords, N, coord_type, master, TAG_COORDS, MPI_COMM_WORLD);
      }
    } 
  }

  MPI_Finalize();
  return 0;
}
