# GRlite core build — two targets:
#   make native       build the native test binary (host gcc/clang)
#   make wasm         build grlite.{js,wasm} via Emscripten, output to web/public/grlite/
#   make test         build + run native tests
#   make clean
#
# Invocation on Windows: from git-bash or PowerShell, run `mingw32-make ...`
# (MinGW's make installs as that name to avoid clashing with cmd's built-in).
#
# Tooling (per memory grlite-emscripten):
#   - Emscripten must be activated in the shell first. PowerShell:
#       & 'C:\Users\dmarks\src\emsdk\emsdk_env.ps1' | Out-Null
#   - Native compiler: gcc (MinGW-w64 on Windows). Override via `make CC=clang`.

OS_NAME := $(shell uname -s 2>/dev/null || echo Windows_NT)

# --- Native compiler -----------------------------------------------------
CC      ?= gcc
# -Wpedantic intentionally omitted: we rely on the well-established GCC
# ##__VA_ARGS__ extension in the test assertion macro (TEST_ASSERT) for
# zero-format-arg assertions. Everything else is strict C11.
CFLAGS  ?= -O2 -g -std=c11 -Wall -Wextra -Wno-unused-parameter
INCLUDES = -Iinclude -Isrc

ifdef SANITIZE
CFLAGS  += -fsanitize=address,undefined -fno-omit-frame-pointer
LDFLAGS += -fsanitize=address,undefined
endif

# --- Emscripten ----------------------------------------------------------
ifeq ($(OS),Windows_NT)
EMCC ?= emcc.bat
else
EMCC ?= emcc
endif

EMCC_FLAGS = -O3 -std=c11 -Wall -Wextra -msimd128 \
             -s ENVIRONMENT=web \
             -s MODULARIZE=1 \
             -s EXPORT_ES6=1 \
             -s ALLOW_MEMORY_GROWTH=1 \
             -s EXPORTED_FUNCTIONS='["_gr_sim_create","_gr_sim_destroy","_gr_sim_step","_gr_sim_step_n","_gr_sim_step_count","_gr_sim_time","_gr_sim_dt","_gr_sim_dx","_gr_sim_width","_gr_sim_height","_gr_sim_field_ptr","_gr_sim_load_scenario","_malloc","_free"]' \
             -s EXPORTED_RUNTIME_METHODS='["HEAPF32","HEAPU8","ccall","cwrap","UTF8ToString","stringToUTF8"]'

# --- Sources -------------------------------------------------------------
CORE_SRCS = src/sim.c src/field.c scenarios/registry.c scenarios/wave_pulse.c
TEST_SRCS = tests/stage01_wave.c

BUILD       = build
NATIVE_DIR  = $(BUILD)/native

WASM_PUBLIC_DIR = ../web/public/grlite
WASM_OUT        = $(WASM_PUBLIC_DIR)/grlite.js

# Native object files mirror the source layout under build/native/.
NATIVE_OBJS = $(patsubst %.c,$(NATIVE_DIR)/%.o,$(CORE_SRCS))

# Test binary names (no extension; MinGW-w64 will silently add .exe — Make tracks both).
TEST_BINS = $(patsubst tests/%.c,$(NATIVE_DIR)/%,$(TEST_SRCS))

# --- Default target ------------------------------------------------------
.PHONY: all native wasm test clean

all: native wasm

# --- Native compile + link -----------------------------------------------
$(NATIVE_DIR)/%.o: %.c include/grlite.h src/sim_internal.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Each test under tests/ links against all core objects.
$(NATIVE_DIR)/%: tests/%.c $(NATIVE_OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(NATIVE_OBJS) -lm $(LDFLAGS) -o $@

native: $(TEST_BINS)

test: native
	@echo "=== running native tests ==="
	@fail=0; for t in $(TEST_BINS); do \
	  echo ""; echo "--- $$t ---"; \
	  "$$t" || fail=1; \
	done; \
	echo ""; \
	if [ $$fail -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "TESTS FAILED"; exit 1; fi

# --- WASM build ----------------------------------------------------------
wasm: $(WASM_OUT)

$(WASM_OUT): $(CORE_SRCS) include/grlite.h src/sim_internal.h
	@mkdir -p $(WASM_PUBLIC_DIR)
	$(EMCC) $(EMCC_FLAGS) $(INCLUDES) $(CORE_SRCS) -o $(WASM_OUT)

# --- Clean ---------------------------------------------------------------
clean:
	rm -rf $(BUILD) $(WASM_PUBLIC_DIR)
