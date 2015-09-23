#include "qute.h"
#include "omp.h"

void do_test1(int num_threads,int N,int M) {

	omp_set_num_threads(num_threads);
	double local_vals[num_threads];
	#pragma omp parallel
	{
		//printf("%d threads: ",omp_get_num_threads());
		double local_val=0;
		int i1=omp_get_thread_num()*N/omp_get_num_threads();
		int i2=(omp_get_thread_num()+1)*N/omp_get_num_threads();
		if (omp_get_thread_num()==omp_get_num_threads()-1) {
			i2=N;
		}
		for (long j=i1; j<i2; j++) {
			int sgn=1;
			for (long i=0; i<M; i++) {
				for (long k=0; k<M; k++) {
					local_val+=sgn*(j+i)*(j+k);
					//sgn*=-1;
				}
			}
		}
		local_vals[omp_get_thread_num()]=local_val;
	}
	double val=0;
	for (int ii=0; ii<num_threads; ii++) val+=local_vals[ii];
	printf("Result = %.0f\n",val);
}

void do_test2(int num_threads,int N,int M) {
	omp_set_num_threads(num_threads);
	double local_vals[num_threads];
	#pragma omp parallel
	{
		//printf("%d threads: ",omp_get_num_threads());
		double local_val=0;
		#pragma omp for
		for (long j=0; j<N; j++) {
			int sgn=1;
			for (long i=0; i<M; i++) {
				for (long k=0; k<M; k++) {
					local_val+=sgn*(j+i)*(j+k);
					//sgn*=-1;
				}
			}
		}
		local_vals[omp_get_thread_num()]=local_val;
	}
	double val=0;
	for (int ii=0; ii<num_threads; ii++) val+=local_vals[ii];
	printf("Result = %.0f\n",val);
}

void do_test3(int num_threads,int N,int M,double *data) {
	omp_set_num_threads(num_threads);
	double local_vals[num_threads];
	#pragma omp parallel
	{
		//printf("%d threads: ",omp_get_num_threads());
		double local_val=0;
		int i1=omp_get_thread_num()*N/omp_get_num_threads();
		int i2=(omp_get_thread_num()+1)*N/omp_get_num_threads();
		if (omp_get_thread_num()==omp_get_num_threads()-1) {
			i2=N;
		}
		for (long j=i1; j<i2; j++) {
			int sgn=1;
			double *data0=&data[j*M];
			for (long i=0; i<M; i++) {
				for (long k=0; k<M; k++) {
					local_val+=sgn*data0[i]*data0[k];
					//sgn*=-1;
				}
			}
		}
		local_vals[omp_get_thread_num()]=local_val;
	}
	double val=0;
	for (int ii=0; ii<num_threads; ii++) val+=local_vals[ii];
	printf("Result = %.0f\n",val);
}

void do_test4(int num_threads,int N,int M,double *data) {
	omp_set_num_threads(num_threads);
	double local_vals[num_threads];
	#pragma omp parallel
	{
		//printf("%d threads: ",omp_get_num_threads());
		double local_val=0;
		#pragma omp for
		for (long j=0; j<N; j++) {
			int sgn=1;
			double *data0=&data[j*M];
			for (long i=0; i<M; i++) {
				for (long k=0; k<M; k++) {
					local_val+=sgn*data0[i]*data0[k];
					//sgn*=-1;
				}
			}
		}
		local_vals[omp_get_thread_num()]=local_val;
	}
	double val=0;
	for (int ii=0; ii<num_threads; ii++) val+=local_vals[ii];
	printf("Result = %.0f\n",val);
}

int main(int argc, char *argv[])
{
	int N=2e2;
	int M=4e3;
	int max_threads=20;

	double *data=(double *)malloc(sizeof(double)*N*M);
	for (long i=0; i<N; i++) {
		for (long j=0; j<M; j++) {
			data[i*M+j]=i+j;
		}
	}

	printf("\n");
    if (0) {
		double elapsed1=0;
		for (int num_threads=1; num_threads<=max_threads; num_threads++) {
			QTime timer; timer.start();
			do_test1(num_threads,N,M);
			double elapsed=timer.elapsed()*1.0/1000;
			if (num_threads==1) elapsed1=elapsed;
			printf("%d threads: %.3f s, rel rate: %.3f, norm rate: %.3f\n",num_threads,elapsed,elapsed1/elapsed,elapsed1/elapsed/num_threads);
		}
	}
	printf("\n");
    if (0) {
		double elapsed1=0;
		for (int num_threads=1; num_threads<=max_threads; num_threads++) {
			QTime timer; timer.start();
			do_test2(num_threads,N,M);
			double elapsed=timer.elapsed()*1.0/1000;
			if (num_threads==1) elapsed1=elapsed;
			printf("%d threads: %.3f s, rel rate: %.3f, norm rate: %.3f\n",num_threads,elapsed,elapsed1/elapsed,elapsed1/elapsed/num_threads);
		}
	}
	printf("\n");
    if (1) {
		double elapsed1=0;
		for (int num_threads=1; num_threads<=max_threads; num_threads++) {
			QTime timer; timer.start();
			do_test3(num_threads,N,M,data);
			double elapsed=timer.elapsed()*1.0/1000;
			if (num_threads==1) elapsed1=elapsed;
			printf("%d threads: %.3f s, rel rate: %.3f, norm rate: %.3f\n",num_threads,elapsed,elapsed1/elapsed,elapsed1/elapsed/num_threads);
		}
	}
	printf("\n");
    if (0) {
		double elapsed1=0;
		for (int num_threads=1; num_threads<=max_threads; num_threads++) {
			QTime timer; timer.start();
			do_test4(num_threads,N,M,data);
			double elapsed=timer.elapsed()*1.0/1000;
			if (num_threads==1) elapsed1=elapsed;
			printf("%d threads: %.3f s, rel rate: %.3f, norm rate: %.3f\n",num_threads,elapsed,elapsed1/elapsed,elapsed1/elapsed/num_threads);
		}
	}

	free(data);

	return 0;
}
