# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/core/BandingEstimator.cpp \
../src/core/BioNJ.cpp \
../src/core/BrentOptimizer.cpp \
../src/core/CommandReader.cpp \
../src/core/Definitions.cpp \
../src/core/Dictionary.cpp \
../src/core/DistanceMatrix.cpp \
../src/core/FileLogger.cpp \
../src/core/FileParser.cpp \
../src/core/Maths.cpp \
../src/core/OptimizedModelParameters.cpp \
../src/core/Optimizer.cpp \
../src/core/PMatrix.cpp \
../src/core/PMatrixDouble.cpp \
../src/core/PMatrixTriple.cpp \
../src/core/PairHmmCalculationWrapper.cpp \
../src/core/PhylogeneticTree.cpp \
../src/core/HmmException.cpp \
../src/core/SequenceElement.cpp \
../src/core/Sequences.cpp \
../src/core/TransitionProbabilities.cpp 

OBJS += \
./src/core/BandingEstimator.o \
./src/core/BioNJ.o \
./src/core/BrentOptimizer.o \
./src/core/CommandReader.o \
./src/core/Definitions.o \
./src/core/Dictionary.o \
./src/core/DistanceMatrix.o \
./src/core/FileLogger.o \
./src/core/FileParser.o \
./src/core/Maths.o \
./src/core/OptimizedModelParameters.o \
./src/core/Optimizer.o \
./src/core/PMatrix.o \
./src/core/PMatrixDouble.o \
./src/core/PMatrixTriple.o \
./src/core/PairHmmCalculationWrapper.o \
./src/core/PhylogeneticTree.o \
./src/core/HmmException.o \
./src/core/SequenceElement.o \
./src/core/Sequences.o \
./src/core/TransitionProbabilities.o 

CPP_DEPS += \
./src/core/BandingEstimator.d \
./src/core/BioNJ.d \
./src/core/BrentOptimizer.d \
./src/core/CommandReader.d \
./src/core/Definitions.d \
./src/core/Dictionary.d \
./src/core/DistanceMatrix.d \
./src/core/FileLogger.d \
./src/core/FileParser.d \
./src/core/Maths.d \
./src/core/OptimizedModelParameters.d \
./src/core/Optimizer.d \
./src/core/PMatrix.d \
./src/core/PMatrixDouble.d \
./src/core/PMatrixTriple.d \
./src/core/PairHmmCalculationWrapper.d \
./src/core/PhylogeneticTree.d \
./src/core/HmmException.d \
./src/core/SequenceElement.d \
./src/core/Sequences.d \
./src/core/TransitionProbabilities.d 


# Each subdirectory must supply rules for building sources it contributes
src/core/%.o: ../src/core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


