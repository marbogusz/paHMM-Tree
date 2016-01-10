# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/heuristics/Band.cpp \
../src/heuristics/BandCalculator.cpp \
../src/heuristics/GuideTree.cpp \
../src/heuristics/ModelEstimator.cpp \
../src/heuristics/Node.cpp \
../src/heuristics/StateTransitionEstimator.cpp \
../src/heuristics/StateTransitionML.cpp \
../src/heuristics/SubstitutionModelEstimator.cpp \
../src/heuristics/TripletAligner.cpp \
../src/heuristics/TripletSamplingTree.cpp 

OBJS += \
./src/heuristics/Band.o \
./src/heuristics/BandCalculator.o \
./src/heuristics/GuideTree.o \
./src/heuristics/ModelEstimator.o \
./src/heuristics/Node.o \
./src/heuristics/StateTransitionEstimator.o \
./src/heuristics/StateTransitionML.o \
./src/heuristics/SubstitutionModelEstimator.o \
./src/heuristics/TripletAligner.o \
./src/heuristics/TripletSamplingTree.o 

CPP_DEPS += \
./src/heuristics/Band.d \
./src/heuristics/BandCalculator.d \
./src/heuristics/GuideTree.d \
./src/heuristics/ModelEstimator.d \
./src/heuristics/Node.d \
./src/heuristics/StateTransitionEstimator.d \
./src/heuristics/StateTransitionML.d \
./src/heuristics/SubstitutionModelEstimator.d \
./src/heuristics/TripletAligner.d \
./src/heuristics/TripletSamplingTree.d 


# Each subdirectory must supply rules for building sources it contributes
src/heuristics/%.o: ../src/heuristics/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


