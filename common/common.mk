# Common Makefile infrastructure for ce2025 projects
# This file contains the common build rules and is included by project Makefiles

# Check if PLATFORM is provided (except for clean target)
ifneq ($(MAKECMDGOALS),clean)
ifndef PLATFORM
$(error PLATFORM not specified. Usage: make PLATFORM=qemu)
endif
endif

PREFIX	?= arm-none-eabi
CC		= $(PREFIX)-gcc
LD		= $(PREFIX)-gcc
AR		= $(PREFIX)-ar
OBJCOPY	= $(PREFIX)-objcopy
OBJDUMP	= $(PREFIX)-objdump
GDB		= $(PREFIX)-gdb

# Silent make
Q ?= @

# Include platform-specific configuration (except for clean target)
ifneq ($(MAKECMDGOALS),clean)
ifeq ($(PLATFORM),stm32)
include ../common/stm32f407.mk
else ifeq ($(PLATFORM),qemu)
include ../common/qemu.mk
else
$(error Invalid PLATFORM '$(PLATFORM)'. Valid platforms: qemu)
endif
endif

# Common build rules
obj/_ELFNAME_%.o:
	@echo "  GEN     $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)echo "const char _elf_name[] = \"$*\";" | \
		$(CC) -x c -c -o $@ $(filter-out -g3,$(CFLAGS)) -

elf/%.elf: obj/_ELFNAME_%.elf.o $(LINKDEPS)
	@echo "  LD      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(LD) $(LDFLAGS) -o $@ $(filter %.o,$^) -Wl,--start-group $(LDLIBS) -Wl,--end-group

obj/%.a:
	@echo "  AR      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(AR) rcs $@ $(filter %.o,$^)

bin/%.bin: elf/%.elf
	@echo "  OBJCOPY $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(OBJCOPY) -Obinary $< $@

bin/%.hex: elf/%.elf
	@echo "  OBJCOPY $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(OBJCOPY) -Oihex $< $@

obj/%.c.o: %.c
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/%.S.o: %.S
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/%.s.o: %.s
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

size: $(PROJECT_OBJS)
	@echo "=== Code Size Analysis ==="
	@arm-none-eabi-size -t $(PROJECT_OBJS)

clean:
	find . -name \*.o -type f -exec rm -f {} \;
	find . -name \*.d -type f -exec rm -f {} \;
	rm -rf obj/ bin/ elf/