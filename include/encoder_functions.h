#ifndef ENCODER_FUNCTIONS_H
#define ENCODER_FUNCTIONS_H

#include "quest.h"

//--------------------------------------------------------------------------------------------------
//SHOR [[9,1,3]]
//-------------------------------------------------------------------------------------------------

void shor_encoder( Qureg qureg, int bit_value);

void shor_syndrome(Qureg qureg, int syndrome[8]);

void shor_error_correction(Qureg qureg, int syndrome[8]);

void shor_decoder(Qureg qureg);

//void shor_print(Qureg qureg);

//-------------------------------------------------------------------------------------------------
//STEANE [[7,1,3]]
//-------------------------------------------------------------------------------------------------

void steane_encoder(Qureg qureg, int bit); //Encodeds bit into the first 7 qubits of registor

void steane_decoder(Qureg qureg, int start);

void stablizer_measure_steane(Qureg qureg, int start_code, int start_aux, int x_synd[3], int z_synd[3]);

void error_correction_steane(Qureg qureg, int start, int z_synd[3], int x_synd[3]);

void steane_decoder(Qureg qureg, int start);

void entangle(Qureg qureg, int start_of_original_registor, int length, int start_of_entangled_reg);

void error(Qureg qureg, int start, int length);

//------------------------------------------------------------------------------------------------
//ENTROPY FUCNTIONS
//------------------------------------------------------------------------------------------------

void prac_print(int num_bits, int Alice[num_bits], int Bob[num_bits]);

double VonNeumannEntropy( Qureg qureg, int num);

double TestVonNeumannEntropy(Qureg qureg, int num);

void error_channel(Qureg qureg, int num, double error);











#endif
