/*
	Basic framework for master / slave workflow to be used for simulated annealing.
*/

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

int         rank;          /* rank of process      */
int         np;            /* number of processes  */
int ndebug; 

MPI_Datatype coord_type; 

#define N 10 // number of coordinates in a solution
#define M 10 // number of solutions to maintain in pool
#define MAX_RUNS 10

/* some colorization stuff stolen from http://linuxgazette.net/issue65/padala.html */
#include <stdio.h>
int enable_color = 1;

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

/* colorize output text */
void textcolor(int attr, int fg, int bg) {
    char command[64];
    if(enable_color) {
        /* Command is the control command to the terminal */
        sprintf(command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
        printf("%s", command);
    }
}

/* debugging output */ 
void debug(const char* fmt, ...) {
    if (ndebug) {
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
}

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
};

struct coord coords[N];
struct coord possible_solution[N];
struct job pool[M];
struct coord *points;
int npoints;

// initialize the list of all possible points (master)
void initialize_points(void){
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

  npoints = n; // global!

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
}

// fill pool with M randomly chosen (but not repeating) candidate solutions (master)
void initialize_pool(void){
  int i, j, k, r;
  for(i = 0; i < M; i++){
    pool[i].coords = (struct coord *)malloc(sizeof(struct coord)*N);
    pool[i].checked_out = 0;
    pool[i].run_count = 0;
    pool[i].owner = -1;
    // FIXME: check that malloc succeeded
    for(j = 0; j < N; j++){
      r = -1;
      while(r < 0){
        r = floor((((double)rand())/((double) RAND_MAX)) * npoints);
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

// master
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

void check_in_candidate(int id){
  int i,j;
  int found = 0;
  for(i = 0; i < M; i++){
    if(pool[i].owner == id){
      pool[i].owner = -1;
      pool[i].checked_out = 0;
      for(j = 0; j < N; j++){
        pool[i].coords[j].e = possible_solution[j].e;
        pool[i].coords[j].i = possible_solution[j].i;
        pool[i].coords[j].n = possible_solution[j].n;
      }
      found = 1;
      break;
    }
  }
  debug("Checked in candidate from %d\n",id);
  print_coords(possible_solution,N);

  if(found == 0){
    debug("Failed to find item to check into pool for slave %d\n",id);
  }
}

/* Do something with a candidate (slave) */
void process_candidate(){
	int i, len, maxlen;
        const char *prog = "R --vanilla < client.R --args ";
        char *cmd;
        char *pos;
        char *tmp;

	debug("Processing Candidate: \n"); 
        print_coords(coords,N);

        // determine the maximum printed length of any given id
        maxlen = 0;
        for(i = 0; i < N; i++){
           len = ceil(log((double)(coords[i].i))/log(10.0));
           if(len > maxlen) maxlen = len;
        }
        // allocate enough space to whole entire command
        cmd = (char *)malloc(sizeof(char)*(strlen(prog) + N*(maxlen+1)));
        tmp = (char *)malloc(sizeof(char)*(maxlen+1));

        // construct the command string
        strcpy(cmd,prog);
        pos = &(cmd[strlen(cmd)]); // points to the null byte at the end of the cmd string
        for(i = 0; i < N; i++){
            sprintf(tmp,"%d.",coords[i].i);
            strcpy(pos,tmp);
            pos = &(cmd[strlen(cmd)]);
        }
        cmd[strlen(cmd)-1] = '\0'; // remove trailing period
	debug("%s\n",cmd);

	// FIXME: actually run and process output

        for(i = 0; i < N; i++){ 
          possible_solution[i].e = coords[i].e;
          possible_solution[i].i = coords[i].i;
          possible_solution[i].n = coords[i].n;
        }
}

int main(int argc, char** argv) {
    int nruns = 0;
    ndebug=1;

    srand(42);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    MPI_Status status;
    int master = 0;
    int tag = 0;
    int tag_shutdown=1;

    debug("Node %d of %d reporting in\n", rank, np);

    /* Create MPI derived datatype for coordinates */
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_UB};
    int          blocklens[4] = {1, 1, 1, 1}; 
    MPI_Aint     disps[4]; 
    int i;

	/* Calculate the displacements of the structure components and the structure overall*/
	MPI_Get_address(coords, disps);
	MPI_Get_address(&coords[0].e, disps+1);
	MPI_Get_address(&coords[0].i, disps+2);
	MPI_Get_address(&coords[1], disps+3);
	int base = disps[0];
	for (i=0;i<4; i++){
		disps[i] -= base;
	}

    MPI_Type_struct(4, blocklens, disps, types, &coord_type); 
    MPI_Type_commit(&coord_type); 

    int candidate = 0;

    if (rank == 0){ 	/* Master process creates work and sends to slaves */

        debug("Master is initializing points...\n");
        initialize_points();
        debug("Okay, points are initialized.\n");

        debug("Master is initializing pool...\n");
        initialize_pool();
        debug("Pool initialized!\n");

	/* Send initial candidates out to slaves */
	for (i=1; i<np; i++){
	    check_out_candidate(i);
	    MPI_Send(coords, N, coord_type, i, tag, MPI_COMM_WORLD);
	}

	/* Go through the rest of the candidate pool and send to slaves */
	for (;;){
	    /* See if any of the slaves have a possible solution for us */
	    MPI_Recv(possible_solution, N, coord_type, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            nruns++;
	    check_in_candidate(status.MPI_SOURCE);

            // FIXME: do something when the possible solution, like put it
            // back in the pool

            // FIXME: at some point, prune or evolve the pool or something

	    /* Send out next candidate to the slave that just reported in */
	    check_out_candidate(status.MPI_SOURCE);
	    MPI_Send(coords, N, coord_type, status.MPI_SOURCE, tag, MPI_COMM_WORLD);

	    /* TODO: decide on stopping criteria */
	    if (nruns > MAX_RUNS){
	        break;
	    }
        }

	/* send shutdown signal to all slaves */
	int dummy_val = 0;
	for (i=1; i<np; i++){
	    MPI_Send(&dummy_val, 1, MPI_INT, i, tag_shutdown, MPI_COMM_WORLD);
	}

    } else { 		/* Slaves receive work from master and process it */		

	int done = 0;
	int work_available = 0;

	/* Receive work from the master process until it says we're done */
	for (;;){

	    /* Check for a shutdown signal from the master */
	    MPI_Iprobe(master, tag_shutdown, MPI_COMM_WORLD, &done, &status);
	    if (done){ 
		int dummy_val = 0;
		MPI_Recv(&dummy_val, 1, MPI_INT, master, tag_shutdown, MPI_COMM_WORLD, &status);
		break;
	    }

	    /* Check if work is available from the master, and if so receive and process it */
	    MPI_Iprobe(master, tag, MPI_COMM_WORLD, &work_available, &status);
	    if (work_available){
	        /* Receive candidate from master to process */ 
		MPI_Recv(coords, N, coord_type, master, tag, MPI_COMM_WORLD, &status);

		/* Perform simulated annealing on candidate */
		process_candidate();

		/* Send our results back to the master */
		MPI_Send(coords, N, coord_type, master, tag, MPI_COMM_WORLD);
	    }
	} 
    }

    MPI_Finalize();
    return 0;
}
