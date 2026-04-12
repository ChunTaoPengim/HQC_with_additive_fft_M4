
CC		= gcc
LD		= gcc


# Silent make
Q ?=


# Common build rules

obj/%.c.o: %.c
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<


portabletest: $(PROJECT_OBJS)
	@mkdir -p bin
	$(LD) -o bin/portabletest $(PROJECT_OBJS)


clean:
	find . -name \*.o -type f -exec rm -f {} \;
	find . -name \*.d -type f -exec rm -f {} \;
	rm -rf obj/ bin/ elf/