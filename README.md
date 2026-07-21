The TCHES paper is available at https://eprint.iacr.org/2026/739 .
# HQC on Arm Cortex-M4

This repository contains implementations of HQC for Arm Cortex-M4 microcontrollers. Code can be tested using QEMU emulation or deployed on a STM32L4R5 Nucleo-144 board for performance benchmarking.

## benchmark results
### multiplications

### New Multiplication
| hqc-1 | hqc-3 | hqc-5|
| :---| :--- | :--- | 
| 1753860| 4836147 | 9290004 |

The table is part of table 4 in the paper.
### New codec 
| operation | hqc-1 | hqc-3 | hqc-5|
| :--- |:---| :--- | :---|
| Encode | 7945| 14878 |26752|
| Decode | 909832| 1348110 | 2337350|

The table is part of table 5 in the paper.
### HQC component
| operation | hqc-1 | hqc-3 | hqc-5 |
| :--- | :--- | :--- | :--- |
| Keygen | 2928092 | 8054559 | 12701974 |
| Encap  | 5187762 | 14427435| 21890093 |
| Decap  | 8748583 | 22876873| 34771718 |

The table is part of table 6 in the paper.
## Prerequisites

Linux (recommended: Ubuntu) or MacOS is required to use this project.
If you are using Windows, we highly recommend using an Ubuntu virtual machine, e.g., using [Virtual Box](https://www.virtualbox.org/).
WSL2 also works for building and QEMU; flashing a board from WSL2 additionally requires USB passthrough via usbipd-win, which is why we suggest VirtualBox for the full hardware flow.

**Arm GCC Toolchain**:
- **Preferred**: Download from [Arm GNU Toolchain Downloads](https://developer.arm.com/downloads/-/arm-gnu-toolchain-downloads)(In our paper we use version 10.3.1, which can be found in https://developer.arm.com/downloads/-/gnu-rm)
- Ubuntu/Debian: `sudo apt install gcc-arm-none-eabi`
- macOS: `brew install gcc-arm-embedded`

**QEMU with Arm support**:
- Ubuntu/Debian: `sudo apt install qemu-system-arm`
- macOS: `brew install qemu`

**For  stlink**:
- Ubuntu/Debian: `sudo apt install stlink-tools`
- macOS: `brew install stlink`

**For serial output**:
- Ubuntu/Debian: `sudo apt install python3-serial`
- Cross-platform: `pip install pyserial`

## Quick Start
We use qemu and pqm4(https://github.com/mupq/pqm4) to test and benchmark our projects:  

### Building and Running

Navigate to any project directory and use these commands:

**QEMU (Emulation)**:
```bash
cd hqc-1/                    # or any project directory
make PLATFORM=qemu            # Build
make run-qemu PLATFORM=qemu   # Build and run
```

### Test multiplications

```bash
cd hqc-1/                    # or any project directory
make PLATFORM=qemu            # Build
make MODE=MUL run-qemu PLATFORM=qemu   # Build and run
```

### Test All Projects
```bash
./run-all-tests.sh            # Run tests for all projects
```


## Output Files

Each build generates platform-specific binaries:
- `bin/qemu.bin` - QEMU binary
- `elf/qemu.elf` - QEMU ELF (for debugging)

## Hardware Setup

### Connecting the STM32L4R5 Nucleo-144 board
- Connect the Micro-USB cable to the top of development board.

### Follow instructions in pqm4
Use the script to clone pqm4 and copy files to the directory:

```bash
./pqm4.sh
cd pqm4
```

### Config and serial output

To view output from the board, make file and connect via serial terminal, use HQC-1 as an example:

```bash
# make files
make -j4 PLATFORM=nucleo-l4r5zi IMPLEMENTATION_PATH=crypto_kem/hqc-1/m4f

# flash the board
st-flash write bin/crypto_kem_hqc-1_m4f_speed.bin 0x8000000
# connect to serial terminal:
pyserial-miniterm /dev/tty.usbmodem* 38400
# or for linux
pyserial-miniterm /dev/ttyACM* 38400
```
**Note on `benchmarks.py`.** pqm4's `speed` test benchmarks the KEM API only
— `keypair`, `encaps`, `decaps` — and pqm4 offers no facility for timing
individual components, so Tables 4 and 5 of the paper cannot be obtained
through it. `pqm4.sh` therefore replaces `mupq/crypto_kem/speed.c` with a
version that additionally times the polynomial multiplication and the
encoder/decoder.

This replacement targets HQC only: it references symbols that exist in our
implementations, so other schemes in the pqm4 tree will no longer build
against it. Please use the commands above rather than `benchmarks.py`, and do
not build other schemes in this working copy. To restore the original
behaviour, run `git -C pqm4/mupq checkout crypto_kem/speed.c`.

The stack usage in Table 6 is obtained by flashing
`bin/crypto_kem_hqc-1_m4f_stack.bin` instead. Repeat with
`IMPLEMENTATION_PATH=crypto_kem/hqc-3/m4f` and `hqc-5/m4f` for the remaining
columns.

## Development

### Common Infrastructure
- `common/`: Shared build system, HAL, and utilities
- `common/common.mk`: Common Makefile infrastructure
- `common/qemu.mk`: QEMU-specific build configuration  
