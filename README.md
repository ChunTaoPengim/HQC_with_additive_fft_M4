# HQC on cortex m4

This repository contains implementations of various cryptographic algorithms cross-compiled for an Arm Cortex-M4 microcontrollers. Code can be tested using QEMU emulation or deployed on a STM32L4R5 Nucleo-144 board for performance benchmarking.

## benchmark results
### multiplications

### New Multiplication
| hqc-1 | hqc-3 | 
| :---| :--- | 
| 1753860| 4836147 |

### New codec 
| operation | hqc-1 | hqc-3 | hqc-5|
| :--- |:---| :--- | :---|
| Encode | 7945| 14878 |26752|
| Decode | 909832| 1348110 | 2337350|

### HQC component
| operation | hqc-1 | hqc-3 | hqc-5 |
| :--- | :--- | :--- | :--- |
| Keygen | 2928092 | 8054559 | 12701974 |
| Encap  | 5187762 | 14427435| 21890093 |
| Decap  | 8748583 | 22876873| 34771718 |

## Prerequisites

Linux (recommended: Ubuntu) or MacOS is required to use this project.
If you are using Windows, we highly recommend using an Ubuntu virtual machine, e.g., using [Virtual Box](https://www.virtualbox.org/).

**Arm GCC Toolchain**:
- **Preferred**: Download from [Arm GNU Toolchain Downloads](https://developer.arm.com/downloads/-/arm-gnu-toolchain-downloads)(In our paper we use version 10.3.1)
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

## Project Structure

Each project follows a consistent structure:
```
project-name/
├── lib/    # Library for fips202
├── src/          # Source codes                
├── tests/        # Test files
├── Makefile                # rules for building and running
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
```

## Development

### Common Infrastructure
- `common/`: Shared build system, HAL, and utilities
- `common/common.mk`: Common Makefile infrastructure
- `common/qemu.mk`: QEMU-specific build configuration  
