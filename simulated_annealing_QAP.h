// Qap4sa.cpp : Defines the entry point for the console application.
//http://www.adaptivebox.net/CILib/code/qapcodes_link.html

#include <iostream>
using namespace std;
#include <fstream>
#include <math.h>
#include <time.h>

const long n_max = 10;                       //matrix size (problem size)
const long infini = 1399999999;            
const long num_temp = 1000;                   // Number of temprature stages

typedef long  type_vecteur[n_max];           //Solution vector
typedef long type_matrice[n_max][n_max];     // for distance and flow matric

long max(long a, long b) { if (a > b) return(a); else return(b); };
double max(double a, double b) { if (a > b) return(a); else return(b); }
long min(long a, long b) { if (a < b) return(a); else return(b); }
double min(double a, double b) { if (a < b) return(a); else return(b); }
void swap(long &a, long &b) { long temp = a; a = b; b = temp; }
															

/*Innitialize the marices*/                                            
type_matrice dist_input = { {0, 1, 2, 3, 4, 5, 6, 7, 8 },
	                        {1, 0, 1, 2, 3, 1, 2, 3, 4 },
							{2, 1, 0, 1, 2, 2, 1, 2, 3 },
							{3, 2, 1, 0, 1, 3, 2, 1, 2 },
							{4, 3, 2, 1, 0, 4, 3, 2, 1 },
							{5, 1, 2, 3, 4, 0, 1, 2, 3 },
							{6, 2, 1, 2, 3, 1, 0, 1, 2 },
							{7, 3, 2, 1, 2, 2, 1, 0, 1 },
							{8, 4, 3, 2, 1, 3, 2, 1, 0 }
};

type_matrice flow_input = { {0, 1, 2, 3, 4, 5, 6, 7, 8 },
	                        {1, 0, 5, 2, 4, 1, 0, 0, 6 },
							{2, 5, 0, 3, 0, 2, 2, 2, 0 },
							{3, 2, 3, 0, 0, 0, 0, 0, 5 },
							{4, 4, 0, 0, 0, 5, 2, 2, 10 },
							{5, 1, 2, 0, 5, 0, 10, 0, 0 },
							{6, 0, 2, 0, 2, 10, 0, 5, 1 },
							{7, 0, 2, 0, 2, 0, 5, 0, 10 },
							{8, 6, 0, 5, 10, 0, 1, 10, 0 }
						};

long  n = 8;                   // size of the problem 
long num_iterations;           // number of iterations
long num_run;                  // How mant times to restart and run the program
long no_run;                  // number of run counter iteration counters
long Cost;                    // cost of current best solution       
type_matrice dist, flow;      // flow and distance matrices
type_vecteur p;               // Solution vector 

/* Constant variables for the random number generator */

const long m = 2147483647; const long m2 = 2145483479;
const long a12 = 63308; const long q12 = 33921; const long r12 = 12979;
const long a13 = -183326; const long q13 = 11714; const long r13 = 2883;
const long a21 = 86098; const long q21 = 24919; const long r21 = 7417;
const long a23 = -539608; const long q23 = 3976; const long r23 = 2071;
const double invm = 4.656612873077393e-10;
long x10 = 12345, x11 = 67890, x12 = 13579,
x20 = 24680, x21 = 98765, x22 = 43210;
 
/*Generate random number*/
double gen_random()
{
	long h, p12, p13, p21, p23;
	h = x10 / q13; p13 = -a13*(x10 - h*q13) - h*r13;
	h = x11 / q12; p12 = a12*(x11 - h*q12) - h*r12;
	if (p13 < 0) p13 = p13 + m; if (p12 < 0) p12 = p12 + m;
	x10 = x11; x11 = x12; x12 = p12 - p13; if (x12 < 0) x12 = x12 + m;
	h = x20 / q23; p23 = -a23*(x20 - h*q23) - h*r23;
	h = x22 / q21; p21 = a21*(x22 - h*q21) - h*r21;
	if (p23 < 0) p23 = p23 + m2; if (p21 < 0) p21 = p21 + m2;
	x20 = x21; x21 = x22; x22 = p21 - p23; if (x22 < 0) x22 = x22 + m2;
	if (x12 < x22) h = x12 - x22 + m; else h = x12 - x22;
	if (h == 0) return(1.0); else return(h*invm);
}

long unif(long low, long high)
{
	return(low + long(double(high - low + 1) *  gen_random()));
}

/*Calculate the change on the cost of flow for the given assignment*/
long calc_delta(long n, type_matrice & dist, type_matrice & flow, type_vecteur & p, long r, long s)
{
	/* calculate the value of move (r, s) on solution p*/
	long d;
	d = (dist[r][r] - dist[s][s])*(flow[p[s]][p[s]] - flow[p[r]][p[r]]) +
		(dist[r][s] - dist[s][r])*(flow[p[s]][p[r]] - flow[p[r]][p[s]]);
	for (long k = 1; k <= n; k = k + 1) if (k != r && k != s)
		d = d + (dist[k][r] - dist[k][s])*(flow[p[k]][p[s]] - flow[p[k]][p[r]]) +
		(dist[r][k] - dist[s][k])*(flow[p[s]][p[k]] - flow[p[r]][p[k]]);
	return(d);
}

void get_dist_flow_input(long &n, type_matrice &dist, type_matrice &flow) //lire ----read file...open--read--close
{

	long i, j;

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
			dist[i][j] = dist_input[i][j];                        //read distance
	}
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
			flow[i][j] = flow_input[i][j];                        //read flow
	}

}

/*Display the input distance and flow matrix */

void display_input_matrix() {

	long i, j;

	get_dist_flow_input(n, dist, flow);
	cout << "Qap4: Nugent et al. Eight departments  : \n";
	cout << "Distance matrix: \n";
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++) {
			cout << " " << dist[i][j];
		}

		cout << '\n';
	}
	cout << "Flow matrix: \n";
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++) {
			cout << " " << flow[i][j];
		}
		cout << '\n';
	}

}

/*Calculate the cost of the assignment (solution) p*/
long calc_cost(long n, type_matrice & dist, type_matrice & flow, type_vecteur & p)
{
	long i, j;
	long c = 0;
	for (i = 1; i <= n; i = i + 1) for (j = 1; j <= n; j = j + 1)
		c = c + dist[i][j] * flow[p[i]][p[j]];
	return(c);
}

/*Generate a random permutation of solution p: initial solution p[i] = {1, 2, 3, 4, 5, 6, 7, 8 } */
void generate_starting_sol(long n, type_vecteur  & p)
{
	long i;
	for (i = 8; i >= 1; i--) p[i] = i;                             // initial solution p[i] = {1, 2, 3, 4, 5, 6, 7, 8}
	for (i = 2; i <= n; i++) swap(p[i - 1], p[unif(i - 1, n)]);     // Randomly swap two solutions two swap
}                                                                              

/*Simulated Annealing Implementation*/
void simu_annealing(long n, type_matrice &dist, type_matrice &flow, type_vecteur & best_sol, long & best_cost, long nb_iterations)
{
	type_vecteur p;
	long i, r, s;                    //r and s are indexes for chossing random value from the solution vector p
	long delta;
	long k = n*(n - 1) / 2, mxfail = k, num_fail, no_iteration;
	long dmin = infini, dmax = 0;
	double t0, tf, beta, tfound, temperature;

	for (i = 1; i <= n; i = i + 1) p[i] = best_sol[i];  // store the best solution to the global p vector
	long Cost = calc_cost(n, dist, flow, p);            // calculate the cost of this solution
	best_cost = Cost;                                   // store cost of p in to to the best cost 

	for (no_iteration = 1; no_iteration <= num_temp; no_iteration++)
	{
		r = unif(1, n);                  //generate random value for r, between 1 and n
		s = unif(1, n - 1);
		if (s >= r) s = s + 1;

		delta = calc_delta(n, dist, flow, p, r, s);  // calculate the random change in cost for solution (r, s)
		if (delta > 0)
		{
			dmin = min(dmin, delta); dmax = max(dmax, delta); 
		};
		Cost = Cost + delta;                                   //the cost is changed by delta when we insert (r, s)
		swap(p[r], p[s]);                                      //swap r & s in the solution vector
															   
	};
	t0 = dmin + (dmax - dmin) / 10.0;
	tf = dmin;
	beta = (t0 - tf) / (nb_iterations*t0*tf);    /// change in temprature/cooling rate
	num_fail = 0;
	tfound = t0;
	temperature = t0;                    //the new temprature t0 after cooled...
	r = 1; s = 2;                        // reset them after the first loop (cost loop), indexes of of p
										 //cout << "simu_annealing step 3 : \n";
	for (no_iteration = 1; no_iteration <= num_iterations - num_temp; no_iteration++)
	{
		temperature = temperature / (1.0 + beta*temperature);   //update the temprature value
																//cout << "temprature value:" << temperature <<endl;
		s = s + 1;                    //increament s
		if (s > n)                    // if s exceded problem size, increment r, if r can still be increamented, increament r and assign in to s
		{
			r = r + 1;
			if (r > n - 1) r = 1;      
			s = r + 1;
		};

		delta = calc_delta(n, dist, flow, p, r, s);

		if ((delta < 0) || (gen_random() < exp(-double(delta) / temperature)) || mxfail == num_fail)  /*accept worst solution if
																								   1. delta negative
																								   2. p(accept_worst_sol)> random
																								   3. max number of fail passed */
		{
			Cost = Cost + delta;
			swap(p[r], p[s]);           /*move ---interchange r and s*/
			num_fail = 0;
		}
		else num_fail = num_fail + 1;

		if (mxfail == num_fail) { beta = 0; temperature = tfound; };

		if (Cost < best_cost)                                                      
		{
			best_cost = Cost;      //keep the current cost as the best cost
			for (i = 1; i <= n; i = i + 1) best_sol[i] = p[i];
			tfound = temperature;
			cout << "Iteration : " << no_iteration;
			cout << "   Solution: [";
			for (int i = 1; i <= n; i++)
			{		
				cout << " " << p[i]  ;
				
			}cout <<"] "<< "  Cost = " << best_cost << '\n';
		};

	};  
};


int main()
{
	get_dist_flow_input(n, dist, flow);             // Get the input matrix 
	display_input_matrix();                          // display the input matrix                      

	cout << "Number of iterations, Number of run : \n";
	cin >> num_iterations >> num_run;                      // input number of iteration and number of run
	for (no_run = 1; no_run <= num_run; no_run++)
	{
		generate_starting_sol(n, p);                               // generate random starting solution
		simu_annealing(n, dist, flow, p, Cost, num_iterations);      // search

	};

	//cout << "final best cost : " << Cost << endl;
	return 0;
}

