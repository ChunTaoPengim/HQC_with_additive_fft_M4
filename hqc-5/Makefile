COMM_DIR = src/common/
LIB_DIR = lib
TARGET_DIR := tests/

#SRC_DIR := src/portable_c/
SRC_DIR := src/m4/

CFLAGS += -I$(LIB_DIR) -I$(COMM_DIR) -I$(SRC_DIR) -I$(TARGET_DIR)

MODE ?= KEM

# ----------------------------------------------------
# 變數定義
# ----------------------------------------------------
# 根據 MODE 來定義要編譯的檔案
ifeq ($(MODE), KEM)
    TARGET := tests/test_kem.c
else ifeq ($(MODE), MUL)
    TARGET := tests/test_mul.c
else ifeq ($(MODE), C)
    SRC_DIR = src/portable_c/
    #TARGET := tests/test_0x11dgf256.c
    #TARGET := tests/test_rs_encoding.c
    #TARGET := tests/test_rs_syndrome.c
    TARGET := tests/test_portableCkem.c
else
    # 預設或錯誤處理
    $(error Unknown MODE: $(MODE). Must be KEM or MUL)
endif


#HQC_SRC := src/code.c src/fft.c src/fips202.c src/gf.c src/hqc.c src/kem.c src/parsing.c src/reduce.c src/reed_muller.c src/reed_solomon.c src/shake_ds.c src/shake_prng.c src/vector.c
#PMUL_SRC:= src/common/basisconversion.c src/common/btfy_consts_gft.c src/common/btfy_ffft.c src/common/encode.c src/common/gf2x.c src/common/gfmul_fft.c src/common/radix_mul.c

HQC_SRC := $(wildcard src/common/*.c) $(wildcard $(SRC_DIR)/*.c) $(wildcard $(LIB_DIR)/*.c)

# M4_ASM_SRC := src/gft_mul_vi_asm.S src/gfv_tower_asm.S src/encode_asm.S src/radix_16_basemul.S
M4_ASM_SRC := $(wildcard $(SRC_DIR)/*.S)

# src = all c files in the project folder
# Separate C and ASM sources
PROJECT_C_SOURCES = $(HQC_SRC) $(TARGET)
PROJECT_ASM_SOURCES = $(M4_ASM_SRC)
PROJECT_C_OBJS = $(addprefix obj/,$(PROJECT_C_SOURCES:.c=.c.o))
PROJECT_ASM_OBJS = $(addprefix obj/,$(PROJECT_ASM_SOURCES:.S=.S.o))
PROJECT_OBJS = $(PROJECT_C_OBJS) $(PROJECT_ASM_OBJS)


ifeq ($(SRC_DIR), src/portable_c/)
    CFLAGS += -I./portable_c/
    include ./portable_c/portable_c.mk
else
    include ../common/common.mk
endif
