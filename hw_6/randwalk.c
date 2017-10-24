#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>




// Custom types //

typedef struct particle_t {
	long x;
	long y;
	long origin_sector;
	long steps_left;
} PARTICLE;




// Global variables //

const double DOUBLE_ERROR = 1e-7;

/// Sector side
long l;
/// Sector number at X axis
long a;
/// Sector number at Y axis
long b;
/// Total number of steps a particle should do
long n;
/// Total number of particles
long N;

double p_left;
double p_right;
double p_up;
double p_down;

/// Sector offset length
long offset;

long field_width;
long field_height;




// Service & processing functions //

void print_help(void) {
	printf("\x1B[4mUsage\x1B[0m: randwalk l a b n N p_left p_right p_up p_down\n");
	printf("Monte-Carlo random walk modelling MPI-based program\n");
	
	printf("\n\n");
	
	printf("\t\x1B[1ml\x1B[0m\n");
	printf("\tSquare side length\n");
	printf("\n");
	
	printf("\t\x1B[1ma\x1B[0m\n");
	printf("\tNumber of squares at X axis (horizontal)\n");
	printf("\n");
	
	printf("\t\x1B[1mb\x1B[0m\n");
	printf("\tNumber of squares at Y axis (vertical)\n");
	printf("\n");
	
	printf("\t\x1B[1mn\x1B[0m\n");
	printf("\tNumber of steps (for each particle)\n");
	printf("\n");
	
	printf("\t\x1B[1mN\x1B[0m\n");
	printf("\tNumber of particles\n");
	printf("\n");
	
	printf("\t\x1B[1mp_left\x1B[0m\n");
	printf("\tProbability of going to the left\n");
	printf("\n");
	
	printf("\t\x1B[1mp_right\x1B[0m\n");
	printf("\tProbability of going to the right\n");
	printf("\n");
	
	printf("\t\x1B[1mp_up\x1B[0m\n");
	printf("\tProbability of going up\n");
	printf("\n");
	
	printf("\t\x1B[1mp_down\x1B[0m\n");
	printf("\tProbability of going down\n");
	printf("\n");
	
	printf("\n");
}



/**
 * @return Sector the 'unit' actually belongs to
 */
inline long get_sector(PARTICLE* unit) {
	return (unit ->y / l) * a + (unit ->x / l);
}



/**
 * @param current_sector ID of the sector which holds 'unit' at the moment
 * @return Sector that should process 'unit'
 */
inline long get_sector_with_offset(PARTICLE* unit, long current_sector) {
	long ruled_sector = get_sector(unit);
	if (ruled_sector == current_sector) {
		return current_sector;
	}
	
	long x_min = l * (current_sector % a);
	long y_min = l * (current_sector / a);
	long x_max = x_min + l - 1;
	long y_max = y_min + l - 1;
	x_min -= offset;
	y_min -= offset;
	x_max += offset;
	y_max += offset;
	
	if (x_min < 0) {
		x_min += field_width;
	}
	else if (x_max > field_width) {
		x_max -= field_width;
	}
	if ((unit ->x < x_min) && (unit ->x > x_max)) {
		return ruled_sector;
	}
	
	if (y_min < 0) {
		y_min += field_height;
	}
	else if (y_max > field_height) {
		y_max -= field_height;
	}
	if ((unit ->y < y_min) && (unit ->y > y_max)) {
		return ruled_sector;
	}
	
	return current_sector;
}



/**
 * @return Position of 'unit' in 1-dimension space
 */
inline long get_through_position(PARTICLE* unit) {
	return unit ->y * field_width + unit ->x;
}











int main(int argc, char** argv) {
	// Read arguments //
	{
		// Check argc
		if (argc != 10) {
			print_help();
			return -1;
		}
		
		// Convert arguments
		l = strtol(argv[1], NULL, 10);
		a = strtol(argv[2], NULL, 10);
		b = strtol(argv[3], NULL, 10);
		n = strtol(argv[4], NULL, 10);
		N = strtol(argv[5], NULL, 10);
		p_left = strtod(argv[6], NULL);
		p_right = strtod(argv[7], NULL);
		p_up = strtod(argv[8], NULL);
		p_down = strtod(argv[9], NULL);
		
		// Check conversions
		if (errno != 0) {
			print_help();
			return -1;
		}
	}
	
	
	
	
	
	return 0;
}




