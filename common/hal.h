#include <stdint.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef HAL_H
#define HAL_H

enum clock_mode {
    CLOCK_FAST,
    CLOCK_BENCHMARK
};

void hal_setup(const enum clock_mode clock);
void hal_send_str(const char* in);
uint64_t hal_get_time(void);
size_t hal_get_stack_size(void);
void hal_spraystack(void);
size_t hal_checkstack(void);

// System call wrapper prototypes
int __wrap__close(int fd);
int __wrap__fstat(int fd, struct stat *buf);
pid_t __wrap__getpid(void);
int __wrap__isatty(int fd);
int __wrap__kill(pid_t pid, int sig);
off_t __wrap__lseek(int fd, off_t offset, int whence);
ssize_t __wrap__read(int fd, void *buf, size_t count);
ssize_t __wrap__write(int fd, const void *buf, size_t count);
void* __wrap__sbrk(int incr);

#endif
