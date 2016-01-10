# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/models/AminoacidSubstitutionModel.cpp \
../src/models/GTRModel.cpp \
../src/models/HKY85Model.cpp \
../src/models/IndelModel.cpp \
../src/models/NegativeBinomialGapModel.cpp \
../src/models/NucleotideSubstitutionModel.cpp \
../src/models/SubstitutionModelBase.cpp 

OBJS += \
./src/models/AminoacidSubstitutionModel.o \
./src/models/GTRModel.o \
./src/models/HKY85Model.o \
./src/models/IndelModel.o \
./src/models/NegativeBinomialGapModel.o \
./src/models/NucleotideSubstitutionModel.o \
./src/models/SubstitutionModelBase.o 

CPP_DEPS += \
./src/models/AminoacidSubstitutionModel.d \
./src/models/GTRModel.d \
./src/models/HKY85Model.d \
./src/models/IndelModel.d \
./src/models/NegativeBinomialGapModel.d \
./src/models/NucleotideSubstitutionModel.d \
./src/models/SubstitutionModelBase.d 


# Each subdirectory must supply rules for building sources it contributes
src/models/%.o: ../src/models/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


