all: bin/qemu.bin

MPS2_DATA_IN_FLASH = 1
LDSCRIPT = obj/ldscript.ld
ARCH_FLAGS += -mcpu=cortex-m4 -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16

CFLAGS += -O3 \
	-Wno-error -fPIC \
	-Wall -Wextra -Wimplicit-function-declaration \
	-Wredundant-decls -Wmissing-prototypes -Wstrict-prototypes \
	-Wundef -Wshadow \
	$(ARCH_FLAGS) \
	-I../common -I../common/mps2 \
	-fno-common -MD \
	-DMPS2_AN386

LDFLAGS += \
	--specs=nosys.specs \
	-Wl,--wrap=_sbrk \
	-Wl,--wrap=_close \
	-Wl,--wrap=_isatty \
	-Wl,--wrap=_kill \
	-Wl,--wrap=_lseek \
	-Wl,--wrap=_read \
	-Wl,--wrap=_write \
	-Wl,--wrap=_fstat \
	-Wl,--wrap=_getpid \
	-ffreestanding \
	-Lobj \
	-T$(LDSCRIPT) \
	$(ARCH_FLAGS)

LIBHAL_OBJS := \
	obj/hal-mps2.c.o \
	obj/randombytes.c.o \
	obj/startup_MPS2.S.o

obj/libpqm4hal.a: $(LIBHAL_OBJS)

obj/hal-mps2.c.o: ../common/hal-mps2.c
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/randombytes.c.o: ../common/randombytes.c
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/startup_MPS2.S.o: ../common/mps2/startup_MPS2.S
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $(if $(MPS2_DATA_IN_FLASH),-DDATA_IN_FLASH) $<

LDLIBS += -lpqm4hal$(if $(NO_RANDOMBYTES),-nornd)
LIBDEPS += obj/libpqm4hal.a

$(LDSCRIPT): ../common/mps2/MPS2.ld
	@printf "  GENLNK  $@\n"; \
	[ -d $(@D) ] || $(Q)mkdir -p $(@D); \
	arm-none-eabi-gcc -x assembler-with-cpp -E -Wp,-P $(CPPFLAGS) $< -o $@

$(LDSCRIPT): CPPFLAGS += -I../common/mps2 $(if $(MPS2_DATA_IN_FLASH),-DDATA_IN_FLASH)

ENABLE_QEMU_TESTS = 1
QEMU = qemu-system-arm
QEMUFLAGS = -M mps2-an386 -nographic -semihosting

elf/qemu.elf: $(PROJECT_OBJS) $(LIBDEPS) $(LDSCRIPT) obj/_ELFNAME_qemu.elf.o

run-qemu: bin/qemu.bin
	$(QEMU) $(QEMUFLAGS) -kernel bin/qemu.bin