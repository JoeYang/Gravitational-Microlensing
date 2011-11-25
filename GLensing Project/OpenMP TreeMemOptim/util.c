#include "util.h"
#include "tree_struct.h"

/* error • Reports an error and then exits with a failed status */
void error(const char *msg) {
  fprintf(stderr, "Error: %s\n", msg);
  exit(1);
}

/* salloc • Safe allocation of memory that exits cleanly in case of error */
void *salloc(size_t bytes) {
  void *ptr = malloc(bytes);
  /* Returning NULL for size 0 is valid */
  if (ptr == NULL && bytes != 0) {
    error("Cannot allocate object");
  }
  return ptr;
}

/* nextline • Loads next line of text, errors if no input */
void nextline(FILE *input, char **buf, size_t *len) {
  if (getline(buf, len, input) == -1) {
    error("getline has returned an error");
  }
}


/* read_lenses • Loads a lens file of format {x, y, (optional)mass} and allocate the correct sized array for the attributes */
void read_lenses(const char *filename, cell ** root) {
  size_t i;
  char c, *line = NULL;
  FILE *fp;
	lens *temp_lens;

  fprintf(stderr, "Reading in lenses...\n");
  if (!(fp = fopen(filename, "r"))) error("Can't open lens file...");

  nobjects = 0;
  // Count the number of lenses we must allocate for (one per line)
  while ((c = getc(fp)) != EOF) {
    if (c == '\n') ++nobjects;
  }
  fprintf(stderr, "Total lenses found: %d\n", (int)nobjects);
  // Seek to the start of the file for actual reading
  fseek(fp, 0, SEEK_SET);

  /* creating the tree*/
  for(i=0; i < nobjects; i++){
    temp_lens = (lens *)salloc(sizeof(lens));                /*allocating memory to temp lens*/
    if(fscanf(fp, "%f %f", &(temp_lens->x), &(temp_lens->y))!=2)
      error("invalid input!");
    temp_lens->m = 1;
    make_tree(root, temp_lens);
    free(temp_lens);
  }

  if (fclose(fp) != 0) error("Can't close lens file...");
  // Deallocate memory used by line
  free(line);
}


/* write_pgm • Output the results as a PGM (portable gray map) image for review */
void write_pgm(unsigned int *results, int pixel_x, int pixel_y, int highest) {
  FILE *fout;
  fprintf(stderr, "Writing resulting image...\n");
  
  time_t clock;
  time(&clock);
  struct tm *timeinfo;
  timeinfo = localtime(&clock);
  char filename[40];
  strftime(filename,40, "img-%Y-%b-%d-%H:%M:%S", timeinfo);
  
  char pgmName[100];
  sprintf(pgmName, "results/%s.pgm", filename);

  if (!(fout = fopen(pgmName, "w"))) system("mkdir results/");
  if (!(fout = fopen(pgmName, "w"))) error("Can't open results file...");
  // Writing the PGM format which starts with P2
  fprintf(fout, "P2\n");
  // Followed by pixel width, height and the value considered white
  fprintf(fout, "%d %d\n", pixel_x, pixel_y);
  fprintf(fout, "%d\n", highest);
  // Print each value in a row of WIDTH length
  int px, py;
  for(py = 0; py < pixel_y; ++py) {
    for(px = 0; px < pixel_x; ++px) {
      fprintf(fout, "%d\n", results[py * pixel_x + px]);
    }
  }
  if (fclose(fout) != 0) error("Can't close results file...");
  
  //fprintf(stderr, "Converting the observation lensing image into png format.\n"); 
  //char pngName[100], command[100];
  //sprintf(pngName, "results/%s.png", filename);
  
  //sprintf(command, "convert %s %s", pgmName, pngName);
  //system(command);
  fprintf(stderr, "The image has been saved in the results folder.\n");
}


int total_r(unsigned int *results, unsigned int size){
  unsigned int i, total = 0;
  for(i = 0; i < size; ++i){
        total += results[i];
  }
  return total;
}

int highest(unsigned int *results, unsigned int size) {
  unsigned int i, highest_count = 0;
  for(i = 0; i < size; ++i){
    if (results[i] > highest_count)
      highest_count = results[i];
  }
  return highest_count;
}
