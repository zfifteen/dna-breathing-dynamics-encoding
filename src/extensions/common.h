/**
 * Common definitions and utilities for C/C++ extensions.
 *
 * Provides platform detection, optimization macros, and shared utilities
 * for high-performance DNA breathing dynamics computations.
 */

#ifndef DNA_BREATHING_COMMON_H
#define DNA_BREATHING_COMMON_H

#include <stdint.h>
#include <stdbool.h>

// Platform detection
#ifdef __APPLE__
    #include <TargetConditionals.h>
    #if defined(__arm64__) || defined(__aarch64__)
        #define APPLE_SILICON 1
    #endif
#endif

// Compiler optimization hints
#ifdef __GNUC__
    #define INLINE inline __attribute__((always_inline))
    #define RESTRICT __restrict__
    #define LIKELY(x) __builtin_expect(!!(x), 1)
    #define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
    #define INLINE inline
    #define RESTRICT
    #define LIKELY(x) (x)
    #define UNLIKELY(x) (x)
#endif

// Memory alignment for SIMD operations
#define CACHE_LINE_SIZE 64
#define SIMD_ALIGNMENT 64

#ifdef __GNUC__
    #define ALIGNED(n) __attribute__((aligned(n)))
#else
    #define ALIGNED(n)
#endif

// High-precision arithmetic includes
#ifdef HAVE_MPFR
    #include <mpfr.h>
    #include <gmp.h>
#endif

// Apple Silicon AMX optimization support
#ifdef APPLE_SILICON
    // AMX tile dimensions
    #define AMX_TILE_ROWS 16
    #define AMX_TILE_COLS 64

    // Metal compute shader support (future)
    // #include <Metal/Metal.h>
#endif

/**
 * Error codes for extension functions.
 */
typedef enum {
    DNA_SUCCESS = 0,
    DNA_ERROR_INVALID_SEQUENCE = -1,
    DNA_ERROR_INVALID_PARAMETER = -2,
    DNA_ERROR_MEMORY_ALLOCATION = -3,
    DNA_ERROR_COMPUTATION = -4,
} dna_error_t;

/**
 * Nucleotide encoding for efficient processing.
 */
typedef enum {
    NUC_A = 0,
    NUC_C = 1,
    NUC_G = 2,
    NUC_T = 3,
    NUC_N = 4,  // Any nucleotide
    NUC_INVALID = -1
} nucleotide_t;

/**
 * Convert ASCII nucleotide character to encoded value.
 */
INLINE nucleotide_t encode_nucleotide(char c) {
    switch (c) {
        case 'A': case 'a': return NUC_A;
        case 'C': case 'c': return NUC_C;
        case 'G': case 'g': return NUC_G;
        case 'T': case 't': return NUC_T;
        case 'N': case 'n': return NUC_N;
        default: return NUC_INVALID;
    }
}

/**
 * Check if nucleotide is valid DNA base.
 */
INLINE bool is_valid_dna_nucleotide(char c) {
    return encode_nucleotide(c) != NUC_INVALID;
}

/**
 * Check if nucleotide is GC (strong hydrogen bonding).
 */
INLINE bool is_gc(nucleotide_t nuc) {
    return nuc == NUC_G || nuc == NUC_C;
}

/**
 * Check if nucleotide is AT (weak hydrogen bonding).
 */
INLINE bool is_at(nucleotide_t nuc) {
    return nuc == NUC_A || nuc == NUC_T;
}

#endif  // DNA_BREATHING_COMMON_H
