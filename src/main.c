#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "quest.h"
#include<time.h>
#include<lapacke.h>
#include<complex.h>
#include"encoder_functions.h"

int main(){

	srand(time(NULL));


	printf("\n");
	printf("--------------------------------------------\n");
	printf("----------QUANTUM ERROR CORRECTION----------\n");
	printf("--------------------------------------------\n");
	printf("\n");

	printf("The intention of this program is to illustrate the effectiveness of different Quatum error correction protocols at protecting a message under different error rates.\n");
	printf("At the moment this program compares Shor's [[9,1,3]] code to Steane's [[7,1,3]] code.\n");

	printf("---------------------------------------------\n");
	printf("\n");
	printf("Please enter a message (no more then 100 characters):\n");
	char message[100];
	fgets(message, sizeof(message),stdin);

	printf("Please enter the error rate per qubit: \n");
	double error;
	scanf("%lf", &error);
	if (error < 0 && error > 1){
		printf("Invalid error rate.\n");
		printf("Error rate must be between 0 and 1.\n");
		return -1;
	}

	printf("\n");
	printf("----------------------------------------------------\n");
	printf("---------RUNNING ERROR CORRECTION PROTOCOL----------\n");
	printf("----------------------------------------------------\n");

	initQuESTEnv();


	//First we consider the Steane code
	unsigned char steane[100];
        int steane_end = 0;	

	Qureg qureg = createQureg(13);
	int *x_synd = malloc(sizeof(int)*3);
	int *z_synd = malloc(sizeof(int)*3);
	for (int i = 0; message[i] != '\0' && message[i] != '\n'; i++){ 

		steane[i] = 0;
		steane_end ++;


		for (int j = 7; j >=0; j--){
			int bit = (message[i] >> j) & 1;

			// The encoding process
			//initQuESTEnv();


			initZeroState(qureg);


			steane_encoder(qureg,bit);

			// The error rate is applied to each qubit
			error_channel(qureg, 13, error);
			// Syndrome measurement
			stablizer_measure_steane( qureg, 0, 7,  x_synd, z_synd);

			//error correction

			error_correction_steane( qureg,0, z_synd, x_synd);

			//Decoding

			steane_decoder(qureg, 0);

			bit = applyQubitMeasurement(qureg,0);

			//build the new byte

			steane[i] = (steane[i]<<1) | bit;


			//free memory
			//
			//finalizeQuESTEnv();
		}
	}
	free(x_synd);
	free(z_synd);
	destroyQureg(qureg);

        //Now we consider the Shor code
        unsigned char shor[100];
	int shor_end = 0;

        Qureg reg = createQureg(17);
        int *syndrome = malloc(sizeof(int)*8);
        for (int i = 0; message[i] != '\0' && message[i] != '\n'; i++){

                shor[i] = 0;
		shor_end ++;


                for (int j = 7; j >=0; j--){
                        int bit = (message[i] >> j) & 1;

                        // The encoding process
                        //initQuESTEnv();


                        initZeroState(reg);


                        shor_encoder(reg,bit);

                        // The error rate is applied to each qubit

                        error_channel(reg, 17, error);

                        // Syndrome measurement

                        shor_syndrome(reg, syndrome);

                        //error correction

                        shor_error_correction( reg, syndrome);

                        //Decoding

                        shor_decoder(reg);

                        bit = applyQubitMeasurement(reg,0);

                        //build the new byte

                        shor[i] = (shor[i]<<1) | bit;



                        //free memory
                        //
                        //finalizeQuESTEnv();
                }
        }
        free(syndrome);
        destroyQureg(reg);



	printf("Using the Shor error correction protocol we recover the following message: \n");
	for (int i = 0; i < shor_end ; i++){
		printf("%c", shor[i]);
	}
	printf("\n");


	printf("Using the Steane error correction protocol we recover the following message: \n");
	for (int i = 0; i < steane_end ; i++){
		printf("%c", steane[i]);
	}
	printf("\n");



        printf("\n");

        printf("-----------------------------------------------\n");
        printf("---------------ENTROPY FUNCTIONS---------------\n");
        printf("-----------------------------------------------\n");

        printf("\n");
        printf("Enter a single bit of information.\n");
        int bit;
        scanf("%d", &bit);
        if (bit != 0 && bit != 1){
                printf("Invalid.\n");
                return -1;
        }
        printf("We are going to calculate the entropy of a certain number of qubits in both the Shor and Steane code\n");
        printf("Enter the number of qubits: \n");
        int num;
        scanf("%d", &num);
        while (num <= 0 && num >= 7){
                printf("Invalid number of qubits. The number should be between 1 and 6\n");
                printf("Try again\n");

                scanf("%d", &num);
        }

        printf("\n");
        printf("----------------------------------------------------\n");
        printf("----------------CALCULATING ENTROPY-----------------\n");
        printf("----------------------------------------------------\n");
        printf("\n");
        Qureg shor_reg = createDensityQureg(9);
        Qureg steane_reg = createDensityQureg(7);

        steane_encoder( steane_reg, bit);
        shor_encoder(shor_reg, bit);

        int traceout_shor[9 - num];
        int traceout_steane[7 - num];

        for (int i = 0; i < 9-num ; i++){
                traceout_shor[i] = i;
        }

        for (int j = 0; j < 7-num ; j++){
                traceout_steane[j] = j;
        }

        Qureg reduced_shor = calcPartialTrace(shor_reg, traceout_shor, 9-num);
        Qureg reduced_steane = calcPartialTrace(steane_reg, traceout_steane , 7-num);


        double steane_entropy = VonNeumannEntropy( reduced_steane, num);
        double shor_entropy = VonNeumannEntropy(reduced_shor, num);

        destroyQureg(shor_reg);
        destroyQureg(steane_reg);
        destroyQureg(reduced_shor);
        destroyQureg(reduced_steane);

        finalizeQuESTEnv();

        printf("The entropy of %d qubits of the reduced steane registor is : %lf\n", num,  steane_entropy);
        printf("The entropy of %d qubits of the reduced shor registor is : %lf\n", num, shor_entropy);

	return 0;
}

