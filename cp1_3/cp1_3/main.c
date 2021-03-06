#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//float Bisection_compute_C(float a, float b);
//float Bisection_compute_f(float a, float N, float *coeffs);
void Bisection_Method(float N, float *coeffs, float a, float b);
void Newton_method(float N, float *coeffs, float a, float b);
double Horner_Algorithm(int N, float *coeffs, float x, int m);
float a, b;/* a and b are the intervals for evaluation of the
		   different methods */
double x_0; /*Newton's method initial guess */


int main(int argc, char *argv[]) /*argc is automatically calculated
								 which is the number of arguments passed to main including the name of
								 the program */
								 /* argv[] is a sequencial list of character array */

{	/* Memory allocation for coefficients*/
	float* coeffs = (float*)malloc(argc * sizeof(float));

	/* This section of the code is to store the coefficients in
	the memory allocated above by malloc */
	int n = argc - 3;
	int i;
	for (i = 0; i < argc - 1; ++i)
	{
		if (n != 0)
		{
			*(coeffs + i) = (float)atof(argv[n]);
			--n;
		}
		/* this section is to store a and b in the
		allocated memory*/
		else
		{
			*(coeffs + i) = (float)atof(argv[i]);
			*(coeffs + i + 1) = (float)atof(argv[i + 1]);
			a = *(coeffs + argc - 2);
			b = *(coeffs + argc - 1);
		}
	}

	
	Bisection_Method(argc, coeffs, a, b); /*calling for the Bisection solution*/
//	float x = 3;
	//printf("agrc = %d", argc);
	int Nx=0;
	if (argc < 3)
		Nx = 5;
	else
	{
	 Nx = argc;
	}
//	printf("\nargc = %d, Nx = %d\n", argc, Nx);
//	int m = 2;
//	Horner_Algorithm(Nx, coeffs, x, m);
	free(coeffs);
}

/* The Bisection Method*/
void Bisection_Method(float N, float *coeffs, float a, float b)
{	int Bisection_iterations_cont = 0;

	//printf("\na = %f; b = %f\n", a, b);
	double machine_ep_check = (fabs(b - a) / 2); /* this checks if the code reached the machine epsilon*/
	float c = 0;/* the midpoint of the interval*/
	float max_iterations_allowed = fabs((log10(b - a) - log10(0.0000012)) / log10(2));
	//printf("max interation = %f", max_iterations_allowed);


	while ((machine_ep_check > 0.0000012) && N != 1)/*machine epsilon for 32 bit system
													and checking for input from command line by making sure the argc (stored in variable
													N is not 1 i.e. not just the name of the program*/
	{	++Bisection_iterations_cont;
		/* this part is to check for available solution within interval given*/
		if (Bisection_iterations_cont>max_iterations_allowed)
		{
			printf(" error if no root is found within the interval a and b");
			goto lable1;
		}

		c = (a + b) / 2;
	
		float f_a = 0, f_c = 0;

		int n = 0;
		f_a = Horner_Algorithm(N, coeffs, a, 1);/* Horner one time to get f_a*/
		f_c = Horner_Algorithm(N, coeffs, c, 1);/* Horner one time to get f_c*/
////////////////////////////////////////////////////////////////////////////////*
//		while (n < (N - 3))
//		{
//			f_a = f_a + (*(coeffs + n))*(pow(a, n));
			/*printf("\ncoeffs[%d]= %f\n", n, *(coeffs + n));*/
			/*printf("\nf(n=%d)= %f; a = %f\n", n, f, a);*/
//			++n;
//		}
//////////////////////////////////////////////////////////////////////////////
//		n = 0;
//		while (n < (N - 3))
//		{
//			f_c = f_c + (*(coeffs + n))*(pow(c, n));
			//printf("\ncoeffs[%d]= %f\n", n, *(coeffs + n));
			//printf("\nf(n=%d)= %f; a = %f\n", n, f, a);
//			++n;
//		}
///////////////////////////////////////////////////////////////////////////////
		
		if ((f_a*f_c) <= 0)
		{
			b = c;
			//	printf("\nif b = %f\n", b);
			//	printf("\nif a = %f\n", a);
		}
		else {
			a = c;
			//	printf("\nelse a = %f\n", a);
			//	printf("\nelse b = %f\n", b);
		}
		machine_ep_check = (fabs(a - b) / 2); /*to check if the progrm is within the machine epsilon*/


		//	printf("\nmachine eps = %.16f\n", machine_ep_check);

	}
	printf("\nx = r_0 = %.5f\n", c);
	printf("\nBisection_iterations_counter = %d\n\n", Bisection_iterations_cont);
lable1:;
}


void Newton_method(float N, float *coeffs, float a, float b)
{
	x_0 = (a + b) / 2; /* find the first guess */


	


}

double Horner_Algorithm(int Nx, float *coeffs, float x, int m) /* if m = 1 find f(x) if m=2 the */
{
//	printf("\n1 N_argc = %d\n", Nx);
	int n = Nx - 4; /* n is the number of coefficients */
//	printf("\n1 n = %d\n", n);
///	printf("\n1 N = %d\n", Nx);
	float* coeffs_prime = (float*)malloc( (n-1) * sizeof(float));
	int i; /* for loop counter*/
	float p = *(coeffs+n); /*x is the initial value variable; p is the derivative*/
	*(coeffs_prime) = p;
	
//	printf("\np(3) = %f\n", p);
	int g = 1;
//	printf("\nhorners loop n= %d\n",n);
//	printf("\nhere\n");
	
	for (i = n-1; i >= 0; --i)
	{
		
	//	p = p + *(coeffs + i)*(pow(x, (i - n)));
		p = p*x + *(coeffs + i );
		*(coeffs_prime+g) = p;
		++g;
//		printf("\nhorners for loop\n");
//		printf("\np(%d) = %f\n",i+1, p);

	}
	n = Nx - 5;
//	printf("\n2 n = %d\n", n);
//	printf("\n2 N = %d\n", Nx);

	float p_prime = *(coeffs_prime);
	////////////////////////
//	for (i = 0; i < n; ++i)
//	{
//		printf("\nccheck p_prime(3) = %f\n", *(coeffs_prime+i));
//	}
	g = 1;
//	printf("\ninitial p_prime(3) = %f\n", p_prime);
	if (m = 2)
	{
		for (i = n - 1; i >= 0; --i)
		{
			p_prime = p_prime*x + *(coeffs_prime + g);
		//	*(coeffs_prime + g) = p;
			++g;
//			printf("\n2horners for loop\n");
//			printf("\np_prime(%d) = %f\n", i + 1, p_prime);

		}
	}
//	printf("\nEnd of horners for loop\n");
//	printf("\np(3) = %f\n\np_prime(3)= %f\n", p, p_prime);

	if (m = 1)
	{
		return p;/* to return the value of f(x)*/
	}
	else
	{
		return p_prime; /* to return the value of f'(x)*/
	}
}







/*Bisection's computing c
float Bisection_compute_C(float a, float b)
{
float z = (a + b) / 2;
//printf("\nc = %f\n", z);
return z;
}

/*Bisection's computing f(a) or f(c) or f(b)
float Bisection_compute_f(float a, float N, float *coeffs)
{
float f = 0;
int n = 0;
while (n < (N - 3))
{
f = f + (*(coeffs + n))*(pow(a, n));
//printf("\ncoeffs[%d]= %f\n", n, *(coeffs + n));
//printf("\nf(n=%d)= %f; a = %f\n", n, f, a);
++n;
}

//printf("\nf = %f\n", f);
return f;

}

/*
/* Newton's Method Function*/
//	Newton_Method();

/*Horner's Method Function*/
//	Horner_Method();



/*
float* coeffs = (float*)malloc(N * sizeof(float));
for (int i = 1; i < (argc - 1); i++)
{
coef[i-1]=
}
*/
