
#include<stdio.h>
#include<math.h>
 
#include<time.h>
#include<stdlib.h>
#include<lapacke.h>
#include<complex.h>
#include"encoder_functions.h"


//----------------------------------------------------------------------------------------------
//SHOR'S [[9,1,3]]
//----------------------------------------------------------------------------------------------

void shor_encoder( Qureg qureg, int bit_value)
{
        if (bit_value < 0 || bit_value > 1){
	       	printf("Invalid\n");
		exit(1);
	}
	//Putting reg into |0>_L state
	//
	//// CN^0_2 CN^0_1 H_0 |000>|000>|000>
	 applyHadamard(qureg,0);
         applyControlledPauliX(qureg,0,1);
	 applyControlledPauliX(qureg,0,2);
	 // CN^3_5 CN^3_4 H_3 1/n(|000> + |111>)|000>|000>
         applyHadamard(qureg,3);
         applyControlledPauliX(qureg,3,4);
	 applyControlledPauliX(qureg,3,5);

         // CN^3_5 CN^3_4 H_3 1/n(|000> + |111>)(|000> + |111>) |000>
         applyHadamard(qureg,6);
	 applyControlledPauliX(qureg,6,7);
	 applyControlledPauliX(qureg,6,8);

	 if (bit_value == 1){
		 applyPauliZ(qureg,0);
		 applyPauliZ(qureg,3);
		 applyPauliZ(qureg,6);
	 }
}

void shor_syndrome(Qureg qureg, int syndrome[8])
{
       	

        //Z_1Z_2
        int anc = 9;
        applyControlledPauliX(qureg,0,anc);
        applyControlledPauliX(qureg,1,anc);

        syndrome[0] = applyQubitMeasurement(qureg,anc);


        //Z_2Z_3
        anc ++;
        applyControlledPauliX(qureg,1,anc);
        applyControlledPauliX(qureg,2,anc);

        syndrome[1] = applyQubitMeasurement(qureg,anc);

        //Z_4Z_5
        anc ++;
        applyControlledPauliX(qureg,3,anc);
        applyControlledPauliX(qureg,4,anc);

        syndrome[2] = applyQubitMeasurement(qureg,anc);

        //Z_4Z_5
        anc ++;
        applyControlledPauliX(qureg,4,anc);
        applyControlledPauliX(qureg,5,anc);

        syndrome[3] = applyQubitMeasurement(qureg,anc);

        //Z_2Z_3
        anc ++;
        applyControlledPauliX(qureg,6,anc);
        applyControlledPauliX(qureg,7,anc);

        syndrome[4] = applyQubitMeasurement(qureg,anc);
        //Z_2Z_3
        anc ++;
        applyControlledPauliX(qureg,7,anc);
        applyControlledPauliX(qureg,8,anc);

        syndrome[5] = applyQubitMeasurement(qureg,anc);

        //X_1X_2X_3X_4X_5X_6
        anc++;
        applyHadamard(qureg,anc);
        applyControlledPauliX(qureg,anc,0);
        applyControlledPauliX(qureg,anc,1);
        applyControlledPauliX(qureg,anc,2);
        applyControlledPauliX(qureg,anc,3);
        applyControlledPauliX(qureg,anc,4);
        applyControlledPauliX(qureg,anc,5);

        applyHadamard(qureg,anc);
        syndrome[6] = applyQubitMeasurement(qureg,anc);
        //X_4X_5X_6X_7X_8X_9
        anc++;
        applyHadamard(qureg,anc);
        applyControlledPauliX(qureg,anc,3);
        applyControlledPauliX(qureg,anc,4);
        applyControlledPauliX(qureg,anc,5);
        applyControlledPauliX(qureg,anc, 6);
        applyControlledPauliX(qureg,anc, 7);
        applyControlledPauliX(qureg,anc, 8);

        applyHadamard(qureg,anc);
        syndrome[7] = applyQubitMeasurement(qureg,anc);
}

void shor_error_correction(Qureg qureg, int syndrome[8])
{
       if (syndrome[0] == 1 && syndrome[1] == 0){
                applyPauliX(qureg,0);
        }
        else if (syndrome[0] == 1 && syndrome[1] == 1){
                applyPauliX(qureg,1);
        }

        else if (syndrome[1] == 1 && syndrome[2] == 0){
                applyPauliX(qureg,2);
        }

        else if (syndrome[2] == 1 && syndrome[3] == 0){
                applyPauliX(qureg,3);
        }
        else if (syndrome[2] == 1 && syndrome[3] == 1){
                applyPauliX(qureg,4);
        }
        else if (syndrome[3] == 1 && syndrome[4] == 0){
                applyPauliX(qureg,5);
        }
        else if (syndrome[4] == 1 && syndrome[5] == 0){
                applyPauliX(qureg,6);
        }
        else if (syndrome[4] == 1 && syndrome[5] == 1){
                applyPauliX(qureg,7);
        }
        else if (syndrome[5] == 1 && syndrome[6] == 0){
                applyPauliX(qureg,8);
        }

        else if (syndrome[6] == 1 && syndrome[7] == 0){
                applyPauliZ(qureg,0);
        }
        else if (syndrome[6] == 1 && syndrome[7] == 1){
                applyPauliZ(qureg,3);
        }
        else if (syndrome[6] == 0 && syndrome[7] == 1){
                applyPauliZ(qureg,6);
        }
}



void shor_decoder(Qureg qureg)
{

        applyControlledPauliX(qureg,0,1);
        applyControlledPauliX(qureg,0,2);


        applyControlledPauliX(qureg,3,4);
        applyControlledPauliX(qureg,3,5);

        applyControlledPauliX(qureg,6,7);
        applyControlledPauliX(qureg,6,8);


        applyHadamard(qureg,0);
        applyHadamard(qureg,3);
        applyHadamard(qureg,6);
}

/*
void shor_print(Qureg qureg, int syndrome[6])
{
        int final[9];
        int count = 0;

        for(int i =0; i<9 ; i++){
                final[i] = applyQubitMeasurement(qureg,i);

        }
        for(int i = 0; i< 3; i++){
                if (final[3*i] == 0){
                        count++;
                }
        }

        printf("The final value in the reg is :");
        for(int i =0; i < 9; i++){
                printf("%d",final[i]);
        }
        printf("\n");

        printf("count is %d\n", count);

        if (count >= 2){
                printf("Hence the encode bit is %d.\n", 0);
        } else {
                printf("Hence the encoded bit is %d.\n",1);
        }

        printf("The syndrome measurement is: ");

        for (int i = 0; i < 6; i++){
                printf("%d", syndrome[i]);
        }
        printf("\n");

        printf("The qubit registor is: ");
        for (int i =0; i< 9; i++){
                printf("%d", applyQubitMeasurement(qureg,i));
        }
        printf("\n");
}
*/


//-------------------------------------------------------------------------------------------------
//STEANE [[7,1,3]]
//-------------------------------------------------------------------------------------------------

void steane_encoder(Qureg qureg, int bit) //Encodeds bit into the first 7 qubits of registor
{
        if( bit != 1 &&  bit != 0){
                printf("Invalid");
                exit(1);
        }


        //The [7,4] Hamming code is used to generate |0>_L and |1>_L

        //|0>_L = 1/8^1/2(|0000000> + |1010101> + |0110011> + |1100110> + |0001111> + |1011010> + |011110> + |1101001>)
        if (bit == 0){

                applyHadamard(qureg, 6);
                applyHadamard(qureg, 5);
                applyHadamard(qureg, 4);

                applyControlledPauliX(qureg,6,3);
                applyControlledPauliX(qureg,5,3);
                applyControlledPauliX(qureg,6,1);
                applyControlledPauliX(qureg,0,2);
                applyControlledPauliX(qureg,0,1);
                applyControlledPauliX(qureg,4,3);
                applyControlledPauliX(qureg,5,2);
                applyControlledPauliX(qureg,6,0);
                applyControlledPauliX(qureg,4,2);
                applyControlledPauliX(qureg,5,0);
                applyControlledPauliX(qureg,4,1);
        } else if (bit == 1){

        // |1>_L = 1/8^1/2(|1111111> + |0101010> + |1001100> + |0011001> + |1110000> + |0100101> + |1000011> + |0010110> )
                applyPauliX(qureg,0);

                applyHadamard(qureg, 6);
                applyHadamard(qureg, 5);
                applyHadamard(qureg, 4);

                applyControlledPauliX(qureg,6,3);
                applyControlledPauliX(qureg,5,3);
                applyControlledPauliX(qureg,6,1);
                applyControlledPauliX(qureg,0,2);
                applyControlledPauliX(qureg,0,1);
                applyControlledPauliX(qureg,4,3);
                applyControlledPauliX(qureg,5,2);
                applyControlledPauliX(qureg,6,0);
                applyControlledPauliX(qureg,4,2);
                applyControlledPauliX(qureg,5,0);
                applyControlledPauliX(qureg,4,1);

        }
}

void steane_decoder(Qureg qureg, int start)
{
        int n = start;
        applyControlledPauliX(qureg,4 + n,1 + n);
        applyControlledPauliX(qureg,5 + n,0 + n);
        applyControlledPauliX(qureg,4 + n,2 + n);
        applyControlledPauliX(qureg,6 + n,0 + n);
        applyControlledPauliX(qureg,5 + n,2 + n);
        applyControlledPauliX(qureg,4 + n,3 + n);
        applyControlledPauliX(qureg,0 + n,1 + n);
        applyControlledPauliX(qureg,0 + n,2 + n);
        applyControlledPauliX(qureg,6 + n,1 + n);
        applyControlledPauliX(qureg,5 + n,3 + n);
        applyControlledPauliX(qureg,6 + n,3 + n);

        applyHadamard(qureg, 4 + n);
        applyHadamard(qureg, 5 + n);
        applyHadamard(qureg, 6 + n);
}


void stablizer_measure_steane(Qureg qureg, int start_code, int start_aux, int x_synd[3], int z_synd[3]){
        //x_synd[] = {0,0,0};
        //z_synd[] = {0,0,0};

        int s = start_code;
        int anc = start_aux;



        //IIIXXXX


        applyHadamard(qureg, anc);

        applyControlledPauliX(qureg,anc,3 + s);
        applyControlledPauliX(qureg,anc,4 + s);
        applyControlledPauliX(qureg,anc,5 + s);
        applyControlledPauliX(qureg,anc,6 + s);

        applyHadamard(qureg,anc);

        x_synd[2] = applyQubitMeasurement(qureg,anc);


        //IXXIIXX

        anc ++;

        applyHadamard(qureg, anc);

        applyControlledPauliX(qureg,anc,1 + s);
        applyControlledPauliX(qureg,anc,2 + s);
        applyControlledPauliX(qureg,anc,5 + s);
        applyControlledPauliX(qureg,anc,6 + s);

        applyHadamard(qureg,anc);
        x_synd[1] = applyQubitMeasurement(qureg,anc);

        //XIXIXIX
        anc ++;
        applyHadamard(qureg, anc);

        applyControlledPauliX(qureg,anc,0 + s);
        applyControlledPauliX(qureg,anc,2 + s);
        applyControlledPauliX(qureg,anc,4 + s);
        applyControlledPauliX(qureg,anc,6 + s);

        applyHadamard(qureg,anc);
        x_synd[0] = applyQubitMeasurement(qureg,anc);
        // IIIZZZZ

        anc++;

        applyControlledPauliX(qureg,3 + s,anc);
        applyControlledPauliX(qureg,4 + s,anc);
        applyControlledPauliX(qureg,5 + s,anc);
        applyControlledPauliX(qureg,6 + s,anc);

        z_synd[2] = applyQubitMeasurement(qureg,anc);

        // IZZIIZZ
        anc ++;
        applyControlledPauliX(qureg,1 + s,anc);
        applyControlledPauliX(qureg,2 + s,anc);
        applyControlledPauliX(qureg,5 + s,anc);
        applyControlledPauliX(qureg,6 + s,anc);

        z_synd[1] = applyQubitMeasurement(qureg,anc);

        // ZIZIZIZ
        anc++;
        applyControlledPauliX(qureg,0 + s,anc);
        applyControlledPauliX(qureg,2 + s,anc);
        applyControlledPauliX(qureg,4 + s,anc);
        applyControlledPauliX(qureg,6 + s,anc);

        z_synd[0] = applyQubitMeasurement(qureg,anc);
}




void error_correction_steane(Qureg qureg, int start, int z_synd[3], int x_synd[3] /*, char *name */)
{
        int Xsynd = 0;
        int Zsynd = 0;
        for (int i = 0; i< 3 ; i++){
                Xsynd += x_synd[i] * pow(2,i);
                Zsynd += z_synd[i] * pow(2,i);


        }



        //printf("%s qubit.\n", name);

        //printf("\n");

        if (Xsynd > 0) {
                applyPauliZ(qureg,start + Xsynd - 1);
                //printf("Acted on p[%d] with Z gate.\n", Xsynd - 1);
        }

        //printf("\n");
        if (Zsynd > 0){
                applyPauliX(qureg, start +  Zsynd - 1);
                //printf("Acted on p[%d] with X gate.\n", Zsynd - 1);
        }
        //printf("\n");
        //printf("\n");
}






void entangle(Qureg qureg, int start_of_original_registor, int length, int start_of_entangled_reg){

        for(int i = 0; i< length; i++){
                applyControlledPauliX(qureg, start_of_original_registor + i, start_of_entangled_reg + i);
        }
}

void error(Qureg qureg, int start, int length)
{
        int error_loc = abs(rand() % length);
        int error_type = abs(rand() % 3);
        int error_loc_2 = abs((3*rand()) % length);

        if (error_type == 0){
                applyPauliX(qureg, start + error_loc);
                //printf("X gate was applied to %d qubit.\n", error_loc + 1);
        }else if (error_type == 1){
                applyPauliZ(qureg, start + error_loc);
                //printf("Z gate was applied to %d qubit.\n", error_loc + 1);
        }else if(error_type == 2){
                applyPauliZ(qureg,start + error_loc);
                applyPauliX(qureg, start + error_loc_2);
                //printf("A Z gate was applied to %d qubit and X gate to %d qubit.\n", error_loc, error_loc_2);
        }


}
//------------------------------------------------------------------------------------------------
//PRINTING FUNCTIONS
//------------------------------------------------------------------------------------------------


void prac_print(int num_bits, int Alice[num_bits], int Bob[num_bits])
{
        printf("Alice's message is calculated to be: [ ");

        for (int i = 0; i < num_bits; i++){
                printf(" %d ", Alice[i]);
        }
        printf("].\n");

        printf("Bob's message is calculated to be: [ ");

        for (int i = 0; i < num_bits; i++){
                printf(" %d ", Bob[i]);
        }
        printf("].\n");
}

//------------------------------------------------------------------------------------------------
//ENTROPY FUNCTIONS
//------------------------------------------------------------------------------------------------

double VonNeumannEntropy( Qureg qureg, int num){

        //create tempory matrix and tempory vector for the eigenvalues of said matrix.
        int n = pow(2,num);

        lapack_complex_double *vector = malloc(sizeof(lapack_complex_double)*n*n);
        double *eigvals = malloc(sizeof(double)*n);

        for(int col = 0; col < n; col++){
                for (int row = 0; row < n ; row++){
                        qcomp z = getDensityQuregAmp(qureg,row,col);
                        vector[col*n + row] = lapack_make_complex_double(creal(z), cimag(z));
                       // printf("(%f,%f) ", creal(z), cimag(z));
                }
                //printf("\n");
        }

        //Calculate eigenvalues
        lapack_int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', n, vector,n,eigvals);

        if (info !=0){
                fprintf(stderr, "zheev failed (%d)\n", info);
        }


        double entropy = 0.0;

        for( int i = 0; i< n; i++){
                if (eigvals[i] > 1e-15){

                        entropy -=eigvals[i] * log2(eigvals[i]);
                }
        }

        free(vector);
        free(eigvals);

        return entropy;
}



double TestVonNeumannEntropy(Qureg qureg, int num){ //Assumes a real,matrix with non zero on the main diagonal

        int n = pow(2,num);
        double *vector = malloc(sizeof(double)*n);
        for (int i = 0; i< n; i++){
                vector[i] = getDensityQuregAmp(qureg,i,i);
        }

        double entropy = 0.0;

        for (int i = 0; i < n; i++){
                if (vector[i] > 1e-15){
                        entropy -= vector[i] * log2(vector[i]);
                }
        }

        free(vector);

        return entropy;
}


void error_channel(Qureg qureg, int num, double error){

        for (int i = 0; i < num; i++){
                int random = rand();
                int type = random & 1;
                double unit_random = (double)random / (double)RAND_MAX ;

                if (unit_random <= error ){
                        if (type == 0){
                                applyPauliX(qureg, i);
                        } else{
                                applyPauliZ(qureg,i);
                        }
                }
        }

}

