#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define epsilon 0.01
#define probability 0.95

#define d_dimension 10
#define d_low_dimension 6

using namespace std;

void matrixmult(float *anchor_pts, float *JLT_matrix, float *transformed_pts)
{

	int x=d_dimension+1;
    int y=d_dimension;
    int m=d_dimension;
    int n=d_low_dimension;


    if (y == m)
    {
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < n; j++)
            {
                transformed_pts[i*n+j] = 0;
                for (int k = 0; k < m; k++)
                {
                    transformed_pts[i*n+j] = transformed_pts[i*n+j] + anchor_pts[i*y+k] * JLT_matrix[k*n+j];
                }
            }
        }
        cout
                << "\n-----------------------------------------------------------\n";
        cout << "\n\nMultiplication of Matrix A and Matrix B :\n\n";
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << "\t" << transformed_pts[i*n+j];
            }
            cout << "\n\n";
        }
    }
    else
    {
        cout << "\n\nMultiplication is not possible";
    }

}

int main(){

	// --- parameter initialization
	double beta = -(log((float)(1.0 - probability)) / log(double(d_dimension + 1)));
	//int d_low_dimension = (((4.0 + 2.0 * beta) * log(double(dimension + 1))) / ((pow(epsilon, 2))/2.0 - ((pow(epsilon, 3))/3.0)) + 1;

	int epsilon_prime; // --- calculating epsilon prime
	if ((1 - sqrt(1 - epsilon)) < (sqrt(1 + epsilon) - 1)) { epsilon_prime = 1 - sqrt(1 - epsilon); }
	else { epsilon_prime = (sqrt(1 + epsilon) + 1); }


	float *anchor_pts= new float[(d_dimension+1)*(d_dimension)];
	float *JLT_matrix = new float[(d_dimension)*(d_low_dimension)];
	float *distances_vector= new float[(d_dimension+1)];
	float *transformed_pts = new float[(d_dimension+1)*(d_low_dimension)];

	// --- defining anchor_pts
	for(int i = 0;i<(d_dimension+1);i++) {
		for(int j=0; j< d_dimension;j++) {

				if(j==i) {
					anchor_pts[i*(d_dimension)+j] = (float)1;
				}
				else {
					anchor_pts[i*(d_dimension)+j] = (float)0;
				}

				if(i==d_dimension){
					if(j==0){ anchor_pts[(d_dimension)*(d_dimension)+j]=(float)1;}
					if(j==1){anchor_pts[(d_dimension)*(d_dimension)+j]=(float)1;}
				}

			}
	}





	//--- initializing a distance vector storing distances from
	//--- each transformed anchor point to the embedded solution
	for (int i =0; i<d_dimension + 1;i++) {
		distances_vector[i]=0.0f;
	}



	// --- defining the JLT transformation matrix

	float jlt_bound1=(1.0/6.0);
	float jlt_bound2=(5.0/6.0);
	float sqrt_3=(float)sqrt(3.0);
	for (int i=0; i< d_dimension; i++) {
		for (int j=0; j<d_low_dimension; j++) {

				float r=(float)rand() / (float)RAND_MAX;
				if (r < jlt_bound1) { JLT_matrix[i*d_low_dimension+j] = -1.0f*sqrt_3;}

				else if((r>=jlt_bound1) && (r<=jlt_bound2)) {JLT_matrix[i*d_low_dimension+j] = 0.0;}

				else if(r>jlt_bound2) { JLT_matrix[i*d_low_dimension+j] = sqrt_3; }

		}

	}

	 cout
	                << "\n-----------------------------------------------------------\n";
	        cout << "\n\nJLT matrix: :\n\n";

	for (int i=0; i<d_dimension; i++){
		for (int j=0; j<d_low_dimension; j++) {

			cout <<"\t"<< JLT_matrix[i*d_low_dimension+j];
		}
		cout << "\n\n";
	}

	 cout
	                << "\n-----------------------------------------------------------\n";
	        cout << "\n\n anchor points:\n\n";

	for (int i=0; i<d_dimension+1; i++){
			for (int j=0; j<d_dimension; j++) {

				cout << "\t" << anchor_pts[i*(d_dimension)+j];
			}
			cout << "\n\n";
		}

	matrixmult(anchor_pts, JLT_matrix, transformed_pts);

	delete [] anchor_pts;
	delete [] JLT_matrix;
	delete [] transformed_pts;
	delete [] distances_vector;

}
