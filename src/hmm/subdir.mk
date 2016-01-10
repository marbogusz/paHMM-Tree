# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hmm/BackwardPairHMM.cpp \
../src/hmm/DpMatrixFull.cpp \
../src/hmm/DpMatrixLoMem.cpp \
../src/hmm/EvolutionaryPairHMM.cpp \
../src/hmm/ForwardPairHMM.cpp \
../src/hmm/PairwiseHmmDeleteState.cpp \
../src/hmm/PairwiseHmmInsertState.cpp \
../src/hmm/PairwiseHmmMatchState.cpp \
../src/hmm/ViterbiPairHMM.cpp 

OBJS += \
./src/hmm/BackwardPairHMM.o \
./src/hmm/DpMatrixFull.o \
./src/hmm/DpMatrixLoMem.o \
./src/hmm/EvolutionaryPairHMM.o \
./src/hmm/ForwardPairHMM.o \
./src/hmm/PairwiseHmmDeleteState.o \
./src/hmm/PairwiseHmmInsertState.o \
./src/hmm/PairwiseHmmMatchState.o \
./src/hmm/ViterbiPairHMM.o 

CPP_DEPS += \
./src/hmm/BackwardPairHMM.d \
./src/hmm/DpMatrixFull.d \
./src/hmm/DpMatrixLoMem.d \
./src/hmm/EvolutionaryPairHMM.d \
./src/hmm/ForwardPairHMM.d \
./src/hmm/PairwiseHmmDeleteState.d \
./src/hmm/PairwiseHmmInsertState.d \
./src/hmm/PairwiseHmmMatchState.d \
./src/hmm/ViterbiPairHMM.d 


# Each subdirectory must supply rules for building sources it contributes
src/hmm/%.o: ../src/hmm/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


