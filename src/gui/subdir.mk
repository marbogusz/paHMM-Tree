# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/gui/Gui.cpp \

GUIOBJS += \
./src/gui/Gui.o \

CPP_DEPS += \
./src/gui/Gui.d \


# Each subdirectory must supply rules for building sources it contributes
src/gui/%.o: ../src/heuristics/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


