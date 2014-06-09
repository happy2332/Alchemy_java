/* updated 5 March 2013 */
/* costwalksat modified 10/24 by H. Kautz */
/* costwalksat by Jiang October 94 based on the following */
/* walksat by Bram Cohen 7/93 */
/* Modifications by Henry Kautz 9/93 */
/* This program continues to accept the old data except that*/
/* it treats each clause with a default cost of unit 1*/
/* Modified Wed May  8 20:47:37 PDT 2002 for -hard option */

#define MAXATOM 10000000		/* maximum possible number of atoms */

#ifdef Huge
#define STOREBLOCK 2000000	/* size of block to malloc each time */
#else
#define STOREBLOCK 1000000	/* size of block to malloc each time */
#define MAXCLAUSE 5000000	/* maximum possible number of clauses */
#endif

#define TRUE 1
#define FALSE 0

#define BIG 100000000

#define MAXLENGTH 1000          /* maximum number of literals which can be in any clause */

#define RANDOM 0
#define PRODUCTSUM 1
#define RECIPROCAL 2
#define ADDITIVE 3
#define BEST 4
#define EXPONENTIAL 5
#define TABU 6

#define NOVALUE -1

#define INIT_PARTIAL 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <signal.h>
#include <limits.h>

/* NOTE: if the -O3 option is not available, then comment out this macro! */
static  long int  scratch;
#define abs(x) ((scratch=(x))>0?scratch: -scratch)


	/* Atoms start at 1 */
	/* Not a is recorded as -1 * a */
	/* One dimensional arrays are statically allocated. */
	/* Two dimensional arrays are dynamically allocated in */
	/* the second dimension only.  */
	/* (Arrays are not entirely dynamically allocated, because */
	/* doing so slows execution by 25% on SGI workstations.) */


 long int  numatom;
 long int  numclause;
 long int  numliterals;

#ifdef Huge
 long int  ** clause;			/* clauses to be satisfied */
				/* indexed as clause[clause_num][literal_num] */
 long int  * cost;			/* cost of each clause */
 long int  * size;			/* length of each clause */
 long int  * false;			/* clauses which are false */
 long int  * wherefalse;		/* where each clause is listed in false */
 long int  * numtruelit;		/* number of true literals in each clause */
#else
 long int  * clause[MAXCLAUSE];	/* clauses to be satisfied */
				/* indexed as clause[clause_num][literal_num] */
 long int  cost[MAXCLAUSE];		/* cost of each clause */
 long int  size[MAXCLAUSE];		/* length of each clause */
 long int  false[MAXCLAUSE];		/* clauses which are false */
 long int  lowfalse[MAXCLAUSE];
 long int  wherefalse[MAXCLAUSE];	/* where each clause is listed in false */
 long int  numtruelit[MAXCLAUSE];	/* number of true literals in each clause */
#endif

 long int  *occurence[2*MAXATOM+1];	/* where each literal occurs */
				/* indexed as occurence[literal+MAXATOM][occurence_num] */

 long int  numoccurence[2*MAXATOM+1];	/* number of times each literal occurs */


 long int  atom[MAXATOM+1];		/* value of each atom */ 
 long int  lowatom[MAXATOM+1];

 long int  changed[MAXATOM+1];		/* step at which atom was last flipped */

 long int  breakcost[MAXATOM+1];	/* the cost of clauses that introduces if var is flipped */

 long int  numfalse;			/* number of false clauses */
 long int  costofnumfalse;		/* cost associated with the number of false clauses */


 long int  costexpected = FALSE;   /*indicate whether cost is expected from the input*/
 long int  abort_flag;

 long int  heuristic;			/* heuristic to be used */
 long int  numerator;			/* make random flip with numerator/denominator frequency */
 long int  denominator;		
 long int  tabu_length;		/* length of tabu list */

long int  numflip;		/* number of changes so far */
long int  numnullflip;		/*  number of times a clause was picked, but no  */
				/*  variable from it was flipped  */

 long int  highestcost=1;               /*The highest cost of a clause violation*/
 long int  eqhighest=0;	/*Is there a clause violated with the highest cost*/
 long int  numhighest=0;      /*Number of violated clauses with the highest cost*/

 long int  hard=0; /* if true never break a highest cost clause */

 long int  selecthigh( long int );   /*find the  long int 'th clause with the highest cost*/

double elapsed_seconds(void);
 long int  countunsat(void);
void countlowunsatcost(  long int  *,  long int  *);
void scantwo( long int  argc, char *argv[],  long int  i, double *varptr);
void scanone( long int  argc, char *argv[],  long int  i,  long int  *varptr);
void init(char initfile[],  long int  initoptions);
void initprob(void);                 /* create a new problem */
void fix( long int  tofix);                 /* changes an atom to make the given clause true */
void flipatom( long int  toflip);           /* changes the assignment of the given literal */
 long int  pickzero( long int  *numbreak, long int  clausesize);	/* return a randomly chosen zero */
                                                /* if none exists, return NOVALUE */
 long int  pickweight( long int  *weight, long int  clausesize); 	/* picks a number based */
						/* on the given weights */
void print_false_clauses_cost(long int  lowbad);
void print_low_assign(long int  lowbad);
void save_low_assign(void);

 long int  pickrandom( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  pickproductsum( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  pickreciprocal( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  pickadditive( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  pickbest( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  pickexponential( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  pickfair( long int  *numbreak, long int  clausesize,  long int  tofix);
 long int  picktabu( long int  *numbreak, long int  clausesize,  long int  tofix);

void handle_interrupt( long int  sig);

long super( long int  i);

 long int  main( long int  argc,char *argv[])
{
     long int  i;			/* loop counter */
     long int  j;			/* another loop counter */
     long int  k;			/* yet another loop counter */
    char initfile[100] = { 0 };
     long int  numrun = 10;
     long int  cutoff = 100000000000L;
     long int  base_cutoff = 100000000000L;
     long int  printonlysol = FALSE;
     long int  printsolcnf = FALSE;
     long int  printfalse = FALSE;
     long int  printlow = FALSE;
     long int  initoptions = FALSE;
     long int  superlinear = FALSE;
     long int  printtrace = FALSE;
    long int  totalflip = 0;	/* total number of flips in all tries so far */
    long int  totalsuccessflip = 0; /* total number of flips in all tries which succeeded so far */
     long int  numtry = 0;		/* total attempts at solutions */
     long int  numsuccesstry = 0;	/* total found solutions */
     long int  numsol = 1;	        /* stop after this many tries succeeds */
     long int  tochange;
     long int  seed;			/* seed for random */
    struct timeval tv;
    struct timezone tzp;
    double expertime;
    long flips_this_solution;
    long int  lowbad;		/* lowest number of bad clauses during try */
    long int  lowcost;		/* lowest cost of bad clauses during try */
     long int  targetcost = 0;		/* the cost at which the program terminates*/
    long x;
    long integer_sum_x = 0;
    double sum_x = 0.0;
    double sum_x_squared = 0.0;
    double mean_x;
    double second_moment_x;
    double variance_x;
    double std_dev_x;
    double std_error_mean_x;
    double seconds_per_flip;
     long int  r;
     long int  sum_r = 0;
    double sum_r_squared = 0.0;
    double mean_r;
    double variance_r;
    double std_dev_r;
    double std_error_mean_r;
     long int  worst_cost, computed_cost;
    float maxtime = 10000.0; // var added by Happy

    gettimeofday(&tv,&tzp);
    seed = (( tv.tv_sec & 0177 ) * 1000000) + tv.tv_usec;
    heuristic = BEST;
    numerator = NOVALUE;
    denominator = 100;

    for (i=1;i < argc;i++)
      {
          
 	  if (strcmp(argv[i],"-withcost") == 0)
	    costexpected = TRUE;
	  else if (strcmp(argv[i],"-seed") == 0)
	    scanone(argc,argv,++i,&seed);
	  else if (strcmp(argv[i],"-targetcost") == 0)
	    scanone(argc,argv,++i,&targetcost);
	  else if (strcmp(argv[i],"-cutoff") == 0)
	    scanone(argc,argv,++i,&cutoff);
	  else if (strcmp(argv[i],"-maxtime") == 0)
	    scantwo(argc,argv,++i,&maxtime);
	  else if (strcmp(argv[i],"-random") == 0)
	    heuristic = RANDOM;
	  else if (strcmp(argv[i],"-productsum") == 0)
	    heuristic = PRODUCTSUM;
	  else if (strcmp(argv[i],"-reciprocal") == 0)
	    heuristic = RECIPROCAL;
	  else if (strcmp(argv[i],"-additive") == 0)
	    heuristic = ADDITIVE;
	  else if (strcmp(argv[i],"-exponential") == 0)
	    heuristic = EXPONENTIAL;
	  else if (strcmp(argv[i],"-best") == 0)
	    heuristic = BEST;
	  else if (strcmp(argv[i],"-noise") == 0)
	    {
		scanone(argc,argv,++i,&numerator);
		scanone(argc,argv,++i,&denominator);
	    }
	  else if (strcmp(argv[i],"-init") == 0  && i < argc-1)
	    sscanf(argv[++i], " %s", initfile);
	  else if (strcmp(argv[i],"-partial") == 0)
	    initoptions = INIT_PARTIAL;
	  else if (strcmp(argv[i],"-super") == 0)
	    superlinear = TRUE;
	  else if (strcmp(argv[i],"-tries") == 0)
	    scanone(argc,argv,++i,&numrun);
	  else if (strcmp(argv[i],"-tabu") == 0)
	    {
		scanone(argc,argv,++i,&tabu_length);
		heuristic = TABU;
	    }
	  else if (strcmp(argv[i],"-low") == 0)
	    printlow = TRUE;
	  else if (strcmp(argv[i],"-sol") == 0)
	    {
		printonlysol = TRUE;
		printlow = TRUE;
	    }
	  else if (strcmp(argv[i],"-bad") == 0)
	    printfalse = TRUE;
	  else if (strcmp(argv[i],"-hard") == 0)
	    hard = TRUE;
	  else if (strcmp(argv[i],"-numsol") == 0)
	    scanone(argc,argv,++i,&numsol);
	  else if (strcmp(argv[i],"-trace") == 0)
	    scanone(argc,argv,++i,&printtrace);
	  else 
	    {
		fprintf(stderr, "Bad argument %s\n", argv[i]);
		fprintf(stderr, "General parameters:\n");
		fprintf(stderr, "  -seed N -cutoff N -tries N\n");
		fprintf(stderr, "  -numsol N = stop after finding N solutions\n");
		fprintf(stderr, "  -init FILE = set vars not included in FILE to false\n");
		fprintf(stderr, "  -partial FILE = set vars not included in FILE randomly\n");
		fprintf(stderr, "  -withcost = input is a set of weighted clauses\n");
		fprintf(stderr, "  -targetcost N = find assignments of cost <= N (MAXSAT)\n");
		fprintf(stderr, "  -hard = never break a highest-cost clause\n");
		fprintf(stderr, "Heuristics:\n");
		fprintf(stderr, "  -noise N M -best -super -tabu N\n");
		fprintf(stderr, "  -productsum -reciprocal -additive -exponential\n");
		fprintf(stderr, "Printing:\n");
		fprintf(stderr, "  -trace N = print statistics every N flips\n");
		fprintf(stderr, "  -sol = print assignments where cost < target\n");
		fprintf(stderr, "  -low = print lowest assignment each try\n");
		fprintf(stderr, "  -bad = print unsat clauses each try\n");
		fprintf(stderr, "  -solcnf = print sat assign in cnf format, and exit\n");
		exit(-1);
	    }
      }
    base_cutoff = cutoff;
    if (numerator==NOVALUE){
	if (heuristic==BEST)
	  numerator = 50;
	else
	  numerator = 0;
    }

    srandom(seed);
#ifdef Huge
    printf("maxwalksat version 20 (Huge)\n");
#else
    printf("maxwalksat version 20\n");
#endif
    printf("seed = %ld\n",seed);
    printf("cutoff = %ld\n",cutoff);
    printf("tries = %ld\n",numrun);
    printf("numsol = %ld\n",numsol);
    printf("targetcost = %ld\n",targetcost);

    printf("heuristic = ");

    switch(heuristic)
      {
	case RANDOM:
	  printf("random");
	  break;
	case PRODUCTSUM:
	  printf("productsum");
	  break;
	case RECIPROCAL:
	  printf("reciprocal");
	  break;
	case ADDITIVE:
	  printf("additive");
	  break;
	case BEST:
	  printf("best");
	  break;
	case EXPONENTIAL:
	  printf("exponential");
	  break;
	case TABU:
	  printf("tabu %ld", tabu_length);
	  break;
      }
    if (numerator>0){
	printf(", noise %ld / %ld", numerator, denominator);
    }
    if (superlinear)
      printf(", super");
    if (printtrace)
      printf(", trace %ld", printtrace);
    if (initfile[0]){
	printf(", init %s", initfile);
	if (initoptions == INIT_PARTIAL)
	  printf(", partial");
    }
    printf("\n");
    
    initprob();

    if (costexpected) printf("clauses contain explicit costs\n");
    else printf("clauses all assigned default cost of 1\n");

    printf("numatom = %ld, numclause = %ld, numliterals = %ld\n",numatom,numclause,numliterals);
    printf("wff read in\n");
    printf("                                           average             average       mean              standard\n");
    printf("    lowest     worst    number                when                over      flips                 error\n");
    printf("      cost    clause    #unsat    #flips     model   success       all      until        std         of\n");
    printf("  this try  this try  this try  this try     found      rate     tries     assign        dev       mean\n");

    signal(SIGINT, (void *) handle_interrupt);
    abort_flag = FALSE;
    numnullflip = 0;
    x = 0; r = 0;
    lowcost = BIG;
    printf("numrun = %ld\n",numrun);
    double total_elapsed_time = elapsed_seconds();
    printf("maxtime = %f\n",maxtime);
     long int  min_false_weight = INT_MAX;
     long int min_false_bad = INT_MAX;
    for(k = 0;k < numrun;k++)
      {
	  
	  init(initfile, initoptions);
	  lowbad = numfalse; 
	  lowcost = costofnumfalse;
	  save_low_assign();
	  numflip = 0;

	  if (superlinear) cutoff = base_cutoff * super(r+1);

	double nextTime = 2.0;
	total_elapsed_time = 0.0;

	  while((numflip < cutoff) && (costofnumfalse > targetcost))
	    {
		if (numflip%10000 == 0)
		{
			total_elapsed_time += elapsed_seconds();
			  if(total_elapsed_time >= nextTime)
			  {
				printf ("***\t%f\t%ld\n",nextTime, lowcost);
				//printf("num flips = %ld\n",numflip);
				nextTime += 2;
			    	if(nextTime > maxtime)
			    		break;
			  }
		}
		if (printtrace && (numflip % printtrace == 0)){
		    printf("%10ld          %10ld%10ld\n", costofnumfalse,numfalse,numflip);
		    fflush(stdout);
		}
		numflip++;
                if ((eqhighest) && (highestcost!=1))
                {/* fprintf(stderr, "number of highest %ld\n", numhighest);*/ 	
			fix(selecthigh(1+random()%numhighest));
                }
		else fix(false[random()%numfalse]);
	        if (costofnumfalse < lowcost)
	        {
		    lowcost = costofnumfalse;
		    lowbad = numfalse;
		    save_low_assign();
		}
	    }
	  numtry++;
	  totalflip += numflip;
	  x += numflip;
	  r ++;
	  if(costofnumfalse<=targetcost)
	    {
		numsuccesstry++;
		totalsuccessflip += numflip;
		integer_sum_x += x;
		sum_x = (double) integer_sum_x;
		sum_x_squared += ((double)x)*((double)x);
		mean_x = sum_x / numsuccesstry;
		if (numsuccesstry > 1){
		    second_moment_x = sum_x_squared / numsuccesstry;
		    variance_x = second_moment_x - (mean_x * mean_x);
		    /* Adjustment for small small sample size */
		    variance_x = (variance_x * numsuccesstry)/(numsuccesstry - 1);
		    std_dev_x = sqrt(variance_x);
		    std_error_mean_x = std_dev_x / sqrt((double)numsuccesstry);
		}
		sum_r += r;
		mean_r = ((double)sum_r)/numsuccesstry;
		sum_r_squared += ((double)r)*((double)r);
		x = 0;
		r = 0;
	    }

	  countlowunsatcost(&computed_cost, &worst_cost);
	  if(lowcost != computed_cost)
	    {
		fprintf(stderr, "Program error, verification of assignment cost fails!\n");
		exit(-1);
	    }

	  if(numsuccesstry == 0)
	    printf("%10ld%10ld%10ld%10ld         *         0         *          *          *          *\n",
		   lowcost,worst_cost,lowbad,numflip);
	  else if (numsuccesstry == 1)
	    printf("%10ld%10ld%10ld%10ld%10ld%10ld%10ld %10.1f          *          *\n",
		   lowcost,worst_cost,lowbad,numflip,totalsuccessflip/numsuccesstry,
		   (numsuccesstry*100)/numtry,totalflip/numsuccesstry,
		   mean_x);
	  else
	    printf("%10ld%10ld%10ld%10ld%10ld%10ld%10ld %10.1f %10.1f %10.1f\n",
		   lowcost,worst_cost,lowbad,numflip,totalsuccessflip/numsuccesstry,
		   (numsuccesstry*100)/numtry,totalflip/numsuccesstry,
		   mean_x, std_dev_x, std_error_mean_x);
	  if (numfalse>0 && printfalse)
	    print_false_clauses_cost(lowbad);
	  if (printlow && (!printonlysol || costofnumfalse<=targetcost))
	    print_low_assign(lowcost);
	  if (lowcost < min_false_weight)
	  {
		min_false_weight = lowcost;
		//min_false_bad = lowbad;
	  }
	  if (numsuccesstry >= numsol) break;
	  if (abort_flag) break;
	  fflush(stdout);
      }
    expertime = elapsed_seconds();
    seconds_per_flip = expertime / totalflip;
    //printf("\ntotal elapsed seconds = %f\n", expertime);
    printf("\ntotal elapsed seconds = %f\n", total_elapsed_time);
    printf("average flips per second = %ld\n", (long)(totalflip/expertime));
    printf("false weight = %ld\n",min_false_weight);
    //printf("false clauses : \n");
    //print_false_clauses_cost(min_false_bad);
    if (heuristic == TABU)
      printf("proportion null flips = %f\n", ((double)numnullflip)/totalflip);
    printf("number of solutions found = %ld\n", numsuccesstry);
    if (numsuccesstry > 0)
      {
	  printf("mean flips until assign = %f\n", mean_x);
	  if (numsuccesstry>1){
	      printf("  variance = %f\n", variance_x);
	      printf("  standard deviation = %f\n", std_dev_x);
	      printf("  standard error of mean = %f\n", std_error_mean_x);
	  }
	  printf("mean seconds until assign = %f\n", mean_x * seconds_per_flip);
	  if (numsuccesstry>1){
	      printf("  variance = %f\n", variance_x * seconds_per_flip * seconds_per_flip);
	      printf("  standard deviation = %f\n", std_dev_x * seconds_per_flip);
	      printf("  standard error of mean = %f\n", std_error_mean_x * seconds_per_flip);
	  }
	  printf("mean restarts until assign = %f\n", mean_r);
	  if (numsuccesstry>1){
	      variance_r = (sum_r_squared / numsuccesstry) - (mean_r * mean_r);
	      variance_r = (variance_r * numsuccesstry)/(numsuccesstry - 1);	   
	      std_dev_r = sqrt(variance_r);
	      std_error_mean_r = std_dev_r / sqrt((double)numsuccesstry);
	      printf("  variance = %f\n", variance_r);
	      printf("  standard deviation = %f\n", std_dev_r);
	      printf("  standard error of mean = %f\n", std_error_mean_r);
	  }
      }

    if (numsuccesstry > 0)
      {
	  printf("ASSIGNMENT ACHIEVING TARGET %ld FOUND\n", targetcost);
	  printf("printsolcnf = %d\n",printsolcnf);
	  if(printsolcnf == TRUE)
	  {
	    for(i = 1;i < numatom+1;i++)
	      printf("v %ld\n", atom[i] == 1 ? i : -i);
	  }
      }
    else
      printf("ASSIGNMENT NOT FOUND\n");
    return 0;
}

long super( long int  i)
{
    long power;
     long int  k;

    if (i<=0){
	fprintf(stderr, "bad argument super(%ld)\n", i);
	exit(1);
    }
    /* let 2^k be the least power of 2 >= (i+1) */
    k = 1;
    power = 2;
    while (power < (i+1)){
	k += 1;
	power *= 2;
    }
    if (power == (i+1)) return (power/2);
    return (super(i - (power/2) + 1));
}

void handle_interrupt( long int  sig)
{
    if (abort_flag) exit(-1);
    abort_flag = TRUE;
}

void scantwo( long int  argc, char *argv[],  long int  i, double *varptr)
{
    if (i>=argc || sscanf(argv[i],"%f",varptr)!=1){
	fprintf(stderr, "Bad argument %s\n", i<argc ? argv[i] : argv[argc-1]);
	exit(-1);
    }
}
void scanone( long int  argc, char *argv[],  long int  i,  long int  *varptr)
{
    if (i>=argc || sscanf(argv[i],"%ld",varptr)!=1){
	fprintf(stderr, "Bad argument %s\n", i<argc ? argv[i] : argv[argc-1]);
	exit(-1);
    }
}

void init(char initfile[],  long int  initoptions)
{
     long int  i;			/* loop counter */
     long int  j;			/* another loop counter */
     long int  thetruelit;
    FILE * infile;
     long int  lit;

    for(i = 0;i < numclause;i++)
      numtruelit[i] = 0;
    numfalse = 0;
    costofnumfalse = 0;
    eqhighest = 0;
    numhighest = 0;

    for(i = 1;i < numatom+1;i++)
      {
	  changed[i] = -BIG;
	  breakcost[i] = 0;
      }

    if (initfile[0] && initoptions!=INIT_PARTIAL){
	for(i = 1;i < numatom+1;i++)
	  atom[i] = 0;
    }
    else {
	for(i = 1;i < numatom+1;i++)
	  atom[i] = random()%2;
    }

    if (initfile[0]){
	if ((infile = fopen(initfile, "r")) == NULL){
	    fprintf(stderr, "Cannot open %s\n", initfile);
	    exit(1);
	}
	i=0;
	while (fscanf(infile, " %ld", &lit)==1){
	    i++;
	    if (abs(lit)>numatom){
		fprintf(stderr, "Bad init file %s\n", initfile);
		exit(1);
	    }
	    if (lit<0) atom[-lit]=0;
	    else atom[lit]=1;
	}
	if (i==0){
	    fprintf(stderr, "Bad init file %s\n", initfile);
	    exit(1);
	}
	close(infile);
	/* printf("read %ld values\n", i); */
    }

    /* Initialize breakcost in the following: */
    for(i = 0;i < numclause;i++)
      {
	  for(j = 0;j < size[i];j++)
	    {
		if((clause[i][j] > 0) == atom[abs(clause[i][j])])
		  {
		      numtruelit[i]++;
		      thetruelit = clause[i][j];
		  }
	    }
	  if(numtruelit[i] == 0)
	    {
		wherefalse[i] = numfalse;
		false[numfalse] = i;
		numfalse++;
		costofnumfalse += cost[i];
                if (highestcost == cost[i])
                {eqhighest = 1; numhighest++;}
	    }
	  else if (numtruelit[i] == 1)
	    {
		breakcost[abs(thetruelit)] += cost[i];
	    }
      }
}

void
print_false_clauses_cost(long int  lowbad)
{
     long int  i, j;
     long int  bad;
     long int  lit, sign;

    printf("Unsatisfied clauses:\n");
    for (i=0; i<numclause; i++){
	bad = TRUE;
	for (j=0; j < size[i]; j++) {
	    lit = clause[i][j];
	    sign = lit > 0 ? 1 : 0;
	    if ( lowatom[abs(lit)] == sign ){
		bad = FALSE;
		break;
	    }
	}
	if (bad){
	    printf("Cost %ld ", cost[i]);
	    for (j=0; j<size[i]; j++){
		printf("%ld ", clause[i][j]);
	    }
	    printf("0\n");
	}
    }
    printf("End unsatisfied clauses\n");
}


void initprob(void)
{
     long int  i;			/* loop counter */
     long int  j;			/* another loop counter */
     long int  lastc;
     long int  nextc;
     long int  *storeptr;
     long int  freestore;
     long int  lit;
    char buf[512];

    while ((lastc = getchar()) == 'c')
      {
	  while ((nextc = getchar()) != EOF && nextc != '\n');
      }
    ungetc(lastc,stdin);
    fgets(buf, 512, stdin);
    if (sscanf(buf, "p wcnf %ld %ld",&numatom,&numclause) == 2){
	costexpected = 1; }
    else if (sscanf(buf, "p cnf %ld %ld",&numatom,&numclause) != 2){
	  fprintf(stderr,"Bad input file\n");
	  exit(-1);}
    if(numatom > MAXATOM)
      {
	  fprintf(stderr,"ERROR - too many atoms\n");
	  exit(-1);
      }

#ifdef Huge
    clause = ( long int  **) malloc(sizeof( long int  *)*(numclause+1));
    cost = ( long int  *) malloc(sizeof( long int )*(numclause+1));
    size = ( long int  *) malloc(sizeof( long int )*(numclause+1));
    false = ( long int  *) malloc(sizeof( long int )*(numclause+1));
    lowfalse = ( long int  *) malloc(sizeof( long int )*(numclause+1));
    wherefalse = ( long int  *) malloc(sizeof( long int )*(numclause+1));
    numtruelit = ( long int  *) malloc(sizeof( long int )*(numclause+1));
#else
    if(numclause > MAXCLAUSE)                     
      {                                      
	  fprintf(stderr,"ERROR - too many clauses\n"); 
	  exit(-1);                              
      }                                        
#endif
    freestore = 0;
    numliterals = 0;
    for(i = 0;i < 2*MAXATOM+1;i++)
      numoccurence[i] = 0;
    for(i = 0;i < numclause;i++)
      {
	  if (costexpected) 
	    {
	     if (scanf("%ld ",&cost[i]) != 1)
		  {
		      fprintf(stderr, "Bad input file\n");
		      exit(-1);
		  }
             else if (cost[i]>highestcost) highestcost = cost[i];
	    }
	  else cost[i] = 1; /*the default cost of a clause violation is unit 1*/

	  size[i] = -1;
	  if (freestore < MAXLENGTH)
	    {
		storeptr = ( long int  *) malloc( sizeof( long int ) * STOREBLOCK );
		freestore = STOREBLOCK;
		//fprintf(stderr,"allocating memory...\n");
	    }
	  clause[i] = storeptr;
	  do
	    {
		size[i]++;
		if(size[i] > MAXLENGTH)
		  {
		      printf("ERROR - clause too long\n");
		      exit(-1);
		  }
		if (scanf("%ld ",&lit) != 1)
		  {
		      fprintf(stderr, "Bad input file\n");
		      exit(-1);
		  }
		if(lit != 0)
		  {
		      *(storeptr++) = lit; /* clause[i][size[i]] = j; */
		      freestore--;
		      numliterals++;
		      numoccurence[lit+MAXATOM]++;
		  }
	    }
	  while(lit != 0);
      }
    if(size[0] == 0)
      {
	  fprintf(stderr,"ERROR - incorrect problem format or extraneous characters\n");
	  exit(-1);
      }

    for(i = 0;i < 2*MAXATOM+1;i++)
      {
	  if (freestore < numoccurence[i])
	    {
		storeptr = ( long int  *) malloc( sizeof( long int ) * STOREBLOCK );
		freestore = STOREBLOCK;
		//fprintf(stderr,"allocating memory...\n");
	    }
	  occurence[i] = storeptr;
	  freestore -= numoccurence[i];
	  storeptr += numoccurence[i];
	  numoccurence[i] = 0;
      }

    for(i = 0;i < numclause;i++)
      {
	  for(j = 0;j < size[i];j++)
	    {
		occurence[clause[i][j]+MAXATOM]
		  [numoccurence[clause[i][j]+MAXATOM]] = i;
		numoccurence[clause[i][j]+MAXATOM]++;
	    }
      }
}


void 
fix( long int  tofix)
{
     long int  numbreak[MAXLENGTH];	/* number of clauses changing */
    /* each atoms would make false */
     long int  i;			/* loop counter */
     long int  j;			/* another loop counter */
     long int  choice;
    static  long int  (*pickcode[])( long int  *numbreak, long int  clausesize,  long int  tofix) = 
      {pickrandom,pickproductsum,pickreciprocal,pickadditive,
	 pickbest,pickexponential,picktabu};

    for(i = 0;i < size[tofix];i++)
      numbreak[i] = breakcost[abs(clause[tofix][i])];

    choice = (pickcode[heuristic])(numbreak,size[tofix],tofix);
    if (choice == NOVALUE)
      numnullflip++;
    else
      flipatom(abs(clause[tofix][choice]));
}


void flipatom( long int  toflip)
{
     long int  i, j;			/* loop counter */
     long int  toenforce;		/* literal to enforce */
    register  long int  cli;
    register  long int  lit;
     long int  numocc;
    register  long int  sz;
    register  long int  * litptr;
     long int  * occptr;

    changed[toflip] = numflip;
    if(atom[toflip] > 0)
      toenforce = -toflip;
    else
      toenforce = toflip;
    atom[toflip] = 1-atom[toflip];
    
    numocc = numoccurence[MAXATOM-toenforce];
    occptr = occurence[MAXATOM-toenforce];
    for(i = 0; i < numocc ;i++)
      {
	  /* cli = occurence[MAXATOM-toenforce][i]; */
	  cli = *(occptr++);

	  if (--numtruelit[cli] == 0){
	      false[numfalse] = cli;
	      wherefalse[cli] = numfalse;
	      numfalse++;
    	      costofnumfalse += cost[cli];
	      /* Decrement toflip's breakcost */
	      breakcost[toflip] -= cost[cli];
             if (cost[cli] == highestcost)
              { eqhighest = 1; numhighest++; }
	  }
	  else if (numtruelit[cli] == 1){
	      /* Find the lit in this clause that makes it true, and inc its breakcount */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  if((lit > 0) == atom[abs(lit)]){
		      breakcost[abs(lit)] += cost[cli];
		      break;
		  }
	      }
	  }
      }
    
    numocc = numoccurence[MAXATOM+toenforce];
    occptr = occurence[MAXATOM+toenforce];
    for(i = 0; i < numocc; i++)
      {
	  /* cli = occurence[MAXATOM+toenforce][i]; */
	  cli = *(occptr++);

	  if (++numtruelit[cli] == 1){
	      numfalse--;
	      costofnumfalse -= cost[cli];
           	      false[wherefalse[cli]] =
		false[numfalse];
	      wherefalse[false[numfalse]] =
		wherefalse[cli];
	      /* Increment toflip's breakcount */
	      breakcost[toflip] += cost[cli];
              if (cost[cli] == highestcost)
              { numhighest--; if (numhighest==0) eqhighest = 0; }
	  }
	  else if (numtruelit[cli] == 2){
	      /* Find the lit in this clause other than toflip that makes it true,
		 and decrement its breakcount */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  if( ((lit > 0) == atom[abs(lit)]) &&
		     (toflip != abs(lit)) ){
		      breakcost[abs(lit)] -= cost[cli];
		      break;
		  }
	      }
	  }
      }
}

 long int  pickrandom( long int  *numbreak, long int  clausesize,  long int  tofix)
	/* returns a random number */
	{
	return(random()%clausesize);
	}

 long int  pickproductsum( long int  *numbreak, long int  clausesize,  long int  tofix)
	/* weights each possibility by the */
	/* product of the product and sum of everything else */
	{
	 long int  i;                             /* a loop counter */
	 long int  weight[MAXLENGTH];             /* weights of each possibility */
	 long int  tochange;                      /* value to return */
	 long int  totalproduct = 1;              /* product of all numbreaks */
	 long int  totalsum = 0;                  /* sum of all numbreaks */

	if(clausesize == 1)
		return(0);
	if((tochange = pickzero(numbreak,clausesize)) != NOVALUE)
		return(tochange);
	for(i = 0;i < clausesize;i++)
		{
		totalproduct *= numbreak[i];
		totalsum += numbreak[i];
		}
	for(i = 0;i < clausesize;i++)
		{
		weight[i] = (totalproduct/numbreak[i])*
			(totalsum-numbreak[i]);
		}
	return(pickweight(weight,clausesize));
	}

 long int  pickreciprocal( long int  *numbreak, long int  clausesize,  long int  tofix)
	/* weights each possibility by its reciprocal*/
	{
	 long int  i;                             /* a loop counter */
	 long int  weight[MAXLENGTH];             /* weights of each possibility */
	 long int  tochange;                      /* value to return */
	 long int  totalproduct = 1;              /* product of all numbreaks */

	if(clausesize == 1)
		return(0);
	if((tochange = pickzero(numbreak,clausesize)) != NOVALUE)
		return(tochange);
	for(i = 0;i < clausesize;i++)
		totalproduct *= numbreak[i];
	for(i = 0;i < clausesize;i++)
		weight[i] = (totalproduct/numbreak[i]);
	return(pickweight(weight,clausesize));
	}

 long int  pickadditive( long int  *numbreak, long int  clausesize,  long int  tofix)
	/* weights each possibility by the sum of the others */
	{
	 long int  i;                             /* a loop counter */
	 long int  weight[MAXLENGTH];             /* weights of each possibility */
	 long int  tochange;                      /* value to return */
	 long int  totalsum = 0;                  /* sum of all numbreaks */

	if(clausesize == 1)
		return(0);
	if((tochange = pickzero(numbreak,clausesize)) != NOVALUE)
		return(tochange);
	for(i = 0;i < clausesize;i++)
		totalsum += numbreak[i];
	for(i = 0;i < clausesize;i++)
		weight[i] = (totalsum-numbreak[i]);
	return(pickweight(weight,clausesize));
	}

 long int  pickbest( long int  *numbreak, long int  clausesize,  long int  tofix)
{
     long int  i;			/* a loop counter */
     long int  best[MAXLENGTH];	/* best possibility so far */
     long int  numbest;		/* how many are tied for best */
     long int  bestvalue;		/* best value so far */

    numbest = 0;
    bestvalue = BIG;

    for (i=0; i< clausesize; i++){
	if (numbreak[i]<=bestvalue){
	    if (numbreak[i]<bestvalue) numbest=0;
	    bestvalue = numbreak[i];
	    best[numbest++] = i;
	}
    }

    if (bestvalue>0 && (random()%denominator < numerator)){

	if (!hard || cost[tofix]>=highestcost) return(random()%clausesize); /* allow random breaks of hard clauses */

	/* only do a random walk breaking non-hard clauses */

	numbest = 0;
	for (i=0; i< clausesize; i++){
	    if (numbreak[i]<highestcost){
		best[numbest++] = i;
	    }
	}
	/* if (numbest==0) { fprintf(stderr, "Wff is not feasible!\n"); exit(1); } */
	if (numbest==0) { return(random()%clausesize); }
    }
    if (numbest == 1) return best[0];
    return(best[random()%numbest]); 



}


 long int  
picktabu( long int  *numbreak, long int  clausesize,  long int  tofix)
{
     long int  i;			/* a loop counter */
     long int  best[MAXLENGTH];	/* best possibility so far */
     long int  numbest;		/* how many are tied for best */
     long int  bestvalue;		/* best value so far */
     long int  tabu_level;
     long int  val;

    numbest = 0;
    bestvalue = BIG;

    if (numerator>0 && (random()%denominator < numerator)){
	for (i=0; i< clausesize; i++){
	    if ((tabu_length < numflip - changed[abs(clause[tofix][i])]) ||
		numbreak[i]==0){
		if (numbreak[i]==0) numbest=0;
		best[numbest++] = i;
	    }
	}
    }
    else {
	for (i=0; i< clausesize; i++){
	    if (numbreak[i]<=bestvalue && 
		((tabu_length < numflip - changed[abs(clause[tofix][i])]) ||
		 numbreak[i]==0)){
		if (numbreak[i]<bestvalue) numbest=0;
		bestvalue = numbreak[i];
		best[numbest++] = i;
	    }
	}
    }
    
    if (numbest == 0) return NOVALUE;
    if (numbest == 1) return best[0];
    return (best[random()%numbest]);
}

 long int  pickexponential( long int  *numbreak, long int  clausesize,  long int  tofix)
	{
	 long int  i;                             /* a loop counter */
	 long int  best[MAXLENGTH];               /* best possibility so far */
	 long int  numbest;                       /* how many are tied for best */
	 long int  bestvalue;                     /* best value so far */
	 long int  weight[MAXLENGTH];             /* weights of each possibility */
	 long int  tochange;                      /* value to return */
	 long int  totalproduct = 1;              /* product of all numbreaks */
	 long int  totalsum = 0;                  /* sum of all numbreaks */

	if(clausesize == 1)
		return(0);
	if((tochange = pickzero(numbreak,clausesize)) != NOVALUE)
		return(tochange);
	for(i = 0;i < clausesize;i++)
		weight[i] = 2*2*2*2*2*2*2*2*2*2*2*2*2*2;
	for(i = 0;i < clausesize;i++)
		{
		weight[i] >>= numbreak[i]; /* weight[i] = weight[i]/(2^numbreak[i]) */
		}

	return(pickweight(weight,clausesize));
	}

 long int  pickzero( long int  *numbreak, long int  clausesize)
	{
	 long int  i;                /* loop counter */
	 long int  numzero = 0;      /* number of zeros so far */
	 long int  select;           /* random number */

	for(i = 0;i < clausesize;i++)
		{
		if(numbreak[i] == 0)
			numzero++;
		}
	if(numzero == 0)
		return(NOVALUE);
	select = random()%numzero;
	for(i = 0;select >= 0;i++)
		{
		if(numbreak[i] == 0)
			select--;
		}

	return(i-1);
	}

 long int  pickweight( long int  *weight, long int  clausesize)
	{
	 long int  i;                /* loop counter */
	 long int  total = 0;        /* sum of weights */
	 long int  select;           /* random number less than total */
	 long int  subtotal = 0;

	for(i = 0;i < clausesize;i++)
		total += weight[i];
	if(total == 0)
		return(random()%clausesize);
	select = random()%total;
	for(i = 0;subtotal <= select;i++)
		subtotal += weight[i];
	return(i-1);
	}

 long int  countunsat(void)
	{
	 long int  i, j, unsat, bad, lit, sign;

	unsat = 0;
	for (i=0;i < numclause;i++)
		{
		bad = TRUE;
		for (j=0; j < size[i]; j++)
			{
			lit = clause[i][j];
			sign = lit > 0 ? 1 : 0;
			if ( atom[abs(lit)] == sign )
				{
				bad = FALSE;
				break;
				}
			}
		if (bad)
			unsat++;
		}
	return unsat;
	}

void countlowunsatcost( long int  * unsatcostptr,  long int  * worstcostptr)
	{
	 long int  i, j, bad, lit, sign, unsatcost, worstcost;

	unsatcost = 0;
	worstcost = 0;
	for (i=0;i < numclause;i++)
		{
		bad = TRUE;
		for (j=0; j < size[i]; j++)
			{
			lit = clause[i][j];
			sign = lit > 0 ? 1 : 0;
			if ( lowatom[abs(lit)] == sign )
				{
				bad = FALSE;
				break;
				}
			}
		if (bad){
			unsatcost += cost[i] ;
			if (cost[i] > worstcost) worstcost = cost[i];
		    }
		}
	* unsatcostptr = unsatcost;
	* worstcostptr = worstcost;
	return;
	}


double 
elapsed_seconds(void)
	{
	double answer;

	static struct rusage prog_rusage;
	static long prev_rusage_seconds = 0;
	static long prev_rusage_micro_seconds = 0;

	getrusage(0, &prog_rusage);
  	answer = (double)(prog_rusage.ru_utime.tv_sec - prev_rusage_seconds)
    		+ ((double)(prog_rusage.ru_utime.tv_usec - prev_rusage_micro_seconds)) / 
		1000000.0 ;
	prev_rusage_seconds = prog_rusage.ru_utime.tv_sec ;
	prev_rusage_micro_seconds = prog_rusage.ru_utime.tv_usec ;
	return answer;
	}


void
print_low_assign(long int  lowbad)
{
     long int  i, j;

    printf("Begin assign with lowest # bad = %ld (atoms assigned true)\n", lowbad);
    j=1;
    for (i=1; i<=numatom; i++){
	if (lowatom[i]>0){
	    printf(" %ld", i);
	    if (j++ % 10 == 0) printf("\n");
	}
    }
    if ((j-1) % 10 != 0) printf("\n");
    printf("End assign\n");
}

void
save_low_assign(void)
{
     long int  i;

    for (i=1; i<=numatom; i++)
      lowatom[i] = atom[i];
}

 long int  selecthigh( long int  high)
{
   long int  counter=0;
   long int  i=0;
  while ((i<numfalse) && (counter<high))
  {
    if (cost[false[i]] == highestcost)
       counter++;
    i++;
  }
 /* fprintf(stderr, "The cost of the selected clause %ld\n", cost[false[i-1]]);*/
  return(false[i-1]);
}
   

