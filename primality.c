/*
 * primality.c  --  Primality-test cost benchmarking tool
 *
 * Six methods: trial division, Fermat, Solovay-Strassen, Miller-Rabin,
 * AKS, Lucas-Lehmer.  Self-contained big integers (no GMP).  Cost is
 * reported as deterministic, machine-independent operation counts at
 * three levels (word multiplications / modular multiplications / method
 * specific counts) so the data is reproducible for research.
 *
 * Build:  gcc primality.c -O2 -o primality.exe
 *
 * Input: a decimal integer, or 0x-prefixed hex, up to 2048 bits.
 */

/* Use MinGW's C99-compliant stdio so %llu works (msvcrt printf does not). */
#define __USE_MINGW_ANSI_STDIO 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <io.h>      /* _dup, _fileno, _fdopen, _close (verbose redirect) */

/* ===================================================================== */
/*  Big integer core                                                     */
/* ===================================================================== */

#define LIMB_BITS 32
#define MAXLIMBS  160              /* ~5120 bits of headroom            */
#define MAXBITS   (MAXLIMBS * LIMB_BITS)

typedef uint32_t limb_t;
typedef uint64_t dlimb_t;

typedef struct {
    limb_t d[MAXLIMBS];            /* little-endian: d[0] = least sig.  */
    int    n;                      /* number of used limbs (0 == zero)  */
} bn;

/* ---- cost counters (reset before each method run) ------------------- */
typedef struct {
    uint64_t word_muls;            /* single limb*limb multiplications   */
    uint64_t word_div_steps;       /* single-limb division steps         */
    uint64_t mod_muls;             /* modular multiplications            */
    uint64_t mod_sqrs;             /* modular squarings                  */
    uint64_t mod_reductions;       /* big divisions used for reduction   */
} costs;

static costs G;
static void costs_reset(void) { memset(&G, 0, sizeof(G)); }

/* format a uint64_t into buf (snprintf %llu is unreliable on old MinGW) */
static char *u64str(uint64_t v, char *buf) {
    char tmp[24]; int p = 0;
    if (v == 0) { buf[0] = '0'; buf[1] = '\0'; return buf; }
    while (v) { tmp[p++] = (char)('0' + v % 10); v /= 10; }
    int j = 0; while (p) buf[j++] = tmp[--p]; buf[j] = '\0';
    return buf;
}

/* ---- basic helpers -------------------------------------------------- */
static void bn_zero(bn *a) { a->n = 0; }

static void bn_norm(bn *a) {
    while (a->n > 0 && a->d[a->n - 1] == 0) a->n--;
}

static int bn_is_zero(const bn *a) { return a->n == 0; }

static void bn_set_u32(bn *a, uint32_t v) {
    if (v == 0) { a->n = 0; } else { a->d[0] = v; a->n = 1; }
}

static void bn_set_u64(bn *a, uint64_t v) {
    a->n = 0;
    if (v & 0xffffffffULL) { a->d[0] = (limb_t)v; a->n = 1; }
    if (v >> 32)           { a->d[1] = (limb_t)(v >> 32); a->n = 2; }
}

static void bn_copy(bn *r, const bn *a) {
    r->n = a->n;
    if (a->n) memcpy(r->d, a->d, (size_t)a->n * sizeof(limb_t));
}

static int bn_cmp(const bn *a, const bn *b) {
    if (a->n != b->n) return a->n < b->n ? -1 : 1;
    for (int i = a->n - 1; i >= 0; i--)
        if (a->d[i] != b->d[i]) return a->d[i] < b->d[i] ? -1 : 1;
    return 0;
}

static int bn_is_one(const bn *a) { return a->n == 1 && a->d[0] == 1; }
static int bn_is_even(const bn *a) { return a->n == 0 || (a->d[0] & 1u) == 0; }

static int bn_bitlen(const bn *a) {
    if (a->n == 0) return 0;
    int b = (a->n - 1) * LIMB_BITS;
    limb_t hi = a->d[a->n - 1];
    while (hi) { b++; hi >>= 1; }
    return b;
}

static int bn_bit(const bn *a, int i) {
    int limb = i / LIMB_BITS, off = i % LIMB_BITS;
    if (limb >= a->n) return 0;
    return (a->d[limb] >> off) & 1u;
}

/* r = a + b */
static void bn_add(bn *r, const bn *a, const bn *b) {
    const bn *x = a->n >= b->n ? a : b;
    const bn *y = a->n >= b->n ? b : a;
    dlimb_t carry = 0;
    int i;
    for (i = 0; i < y->n; i++) {
        dlimb_t s = (dlimb_t)x->d[i] + y->d[i] + carry;
        r->d[i] = (limb_t)s; carry = s >> 32;
    }
    for (; i < x->n; i++) {
        dlimb_t s = (dlimb_t)x->d[i] + carry;
        r->d[i] = (limb_t)s; carry = s >> 32;
    }
    r->n = x->n;
    if (carry) r->d[r->n++] = (limb_t)carry;
}

/* r = a - b, requires a >= b */
static void bn_sub(bn *r, const bn *a, const bn *b) {
    int64_t borrow = 0;
    int i;
    for (i = 0; i < b->n; i++) {
        int64_t s = (int64_t)a->d[i] - b->d[i] - borrow;
        if (s < 0) { s += (int64_t)1 << 32; borrow = 1; } else borrow = 0;
        r->d[i] = (limb_t)s;
    }
    for (; i < a->n; i++) {
        int64_t s = (int64_t)a->d[i] - borrow;
        if (s < 0) { s += (int64_t)1 << 32; borrow = 1; } else borrow = 0;
        r->d[i] = (limb_t)s;
    }
    r->n = a->n;
    bn_norm(r);
}

static void bn_add_small(bn *r, const bn *a, uint32_t v) {
    bn t; bn_set_u32(&t, v); bn_add(r, a, &t);
}
static void bn_sub_small(bn *r, const bn *a, uint32_t v) {
    bn t; bn_set_u32(&t, v); bn_sub(r, a, &t);
}

/* shift left by k bits */
static void bn_shl(bn *r, const bn *a, int k) {
    if (a->n == 0 || k == 0) { bn_copy(r, a); return; }
    int limbs = k / LIMB_BITS, bits = k % LIMB_BITS;
    bn tmp; bn_zero(&tmp);
    if (bits == 0) {
        for (int i = a->n - 1; i >= 0; i--) tmp.d[i + limbs] = a->d[i];
    } else {
        for (int i = a->n; i >= 0; i--) {
            limb_t lo = (i - 1 >= 0) ? a->d[i - 1] : 0;
            limb_t hi = (i < a->n) ? a->d[i] : 0;
            tmp.d[i + limbs] = (limb_t)(((dlimb_t)hi << bits) | (lo >> (LIMB_BITS - bits)));
        }
    }
    tmp.n = a->n + limbs + 1;
    bn_norm(&tmp);
    bn_copy(r, &tmp);
}

/* shift right by k bits */
static void bn_shr(bn *r, const bn *a, int k) {
    int limbs = k / LIMB_BITS, bits = k % LIMB_BITS;
    if (limbs >= a->n) { bn_zero(r); return; }
    bn tmp; bn_zero(&tmp);
    int m = a->n - limbs;
    if (bits == 0) {
        for (int i = 0; i < m; i++) tmp.d[i] = a->d[i + limbs];
    } else {
        for (int i = 0; i < m; i++) {
            limb_t lo = a->d[i + limbs];
            limb_t hi = (i + limbs + 1 < a->n) ? a->d[i + limbs + 1] : 0;
            tmp.d[i] = (limb_t)((lo >> bits) | ((dlimb_t)hi << (LIMB_BITS - bits)));
        }
    }
    tmp.n = m;
    bn_norm(&tmp);
    bn_copy(r, &tmp);
}

/* r = a * v (small) */
static void bn_mul_small(bn *r, const bn *a, uint32_t v) {
    dlimb_t carry = 0;
    bn tmp; bn_zero(&tmp);
    for (int i = 0; i < a->n; i++) {
        dlimb_t p = (dlimb_t)a->d[i] * v + carry;
        G.word_muls++;
        tmp.d[i] = (limb_t)p; carry = p >> 32;
    }
    tmp.n = a->n;
    if (carry) tmp.d[tmp.n++] = (limb_t)carry;
    bn_norm(&tmp);
    bn_copy(r, &tmp);
}

/* r = a * b (schoolbook) */
static void bn_mul(bn *r, const bn *a, const bn *b) {
    if (a->n == 0 || b->n == 0) { bn_zero(r); return; }
    bn tmp; bn_zero(&tmp);
    for (int i = 0; i < a->n + b->n; i++) tmp.d[i] = 0;
    for (int i = 0; i < a->n; i++) {
        dlimb_t carry = 0;
        dlimb_t ai = a->d[i];
        for (int j = 0; j < b->n; j++) {
            dlimb_t p = ai * b->d[j] + tmp.d[i + j] + carry;
            G.word_muls++;
            tmp.d[i + j] = (limb_t)p; carry = p >> 32;
        }
        tmp.d[i + b->n] = (limb_t)carry;
    }
    tmp.n = a->n + b->n;
    bn_norm(&tmp);
    bn_copy(r, &tmp);
}

static void bn_sqr(bn *r, const bn *a) { bn_mul(r, a, a); }

static int clz32(uint32_t x) {
    int c = 0;
    for (int i = 31; i >= 0; i--) { if ((x >> i) & 1u) break; c++; }
    return c;
}

/* q = a / b, r = a % b.  Either q or r may be NULL.  b != 0. */
static void bn_divmod(bn *q, bn *r, const bn *a, const bn *b) {
    int cmp = bn_cmp(a, b);
    if (cmp < 0) { if (q) bn_zero(q); if (r) bn_copy(r, a); return; }
    if (cmp == 0) { if (q) bn_set_u32(q, 1); if (r) bn_zero(r); return; }

    G.mod_reductions++;

    if (b->n == 1) {                      /* single-limb divisor */
        uint32_t bb = b->d[0];
        bn qq; bn_zero(&qq);
        dlimb_t rem = 0;
        for (int i = a->n - 1; i >= 0; i--) {
            dlimb_t cur = (rem << 32) | a->d[i];
            qq.d[i] = (limb_t)(cur / bb);
            rem = cur % bb;
            G.word_div_steps++;
        }
        qq.n = a->n; bn_norm(&qq);
        if (q) bn_copy(q, &qq);
        if (r) { bn_zero(r); if (rem) { r->d[0] = (limb_t)rem; r->n = 1; } }
        return;
    }

    /* Knuth Algorithm D (Hacker's Delight divmnu, base 2^32) */
    int n = b->n;
    int m = a->n - n;
    int s = clz32(b->d[n - 1]);

    limb_t vn[MAXLIMBS];
    if (s) {
        for (int i = n - 1; i > 0; i--)
            vn[i] = (b->d[i] << s) | (b->d[i - 1] >> (LIMB_BITS - s));
        vn[0] = b->d[0] << s;
    } else {
        for (int i = 0; i < n; i++) vn[i] = b->d[i];
    }

    limb_t un[MAXLIMBS + 1];
    un[a->n] = s ? (a->d[a->n - 1] >> (LIMB_BITS - s)) : 0;
    if (s) {
        for (int i = a->n - 1; i > 0; i--)
            un[i] = (a->d[i] << s) | (a->d[i - 1] >> (LIMB_BITS - s));
        un[0] = a->d[0] << s;
    } else {
        for (int i = 0; i < a->n; i++) un[i] = a->d[i];
    }

    bn qq; bn_zero(&qq);
    for (int j = m; j >= 0; j--) {
        dlimb_t num = ((dlimb_t)un[j + n] << 32) | un[j + n - 1];
        dlimb_t qhat = num / vn[n - 1];
        dlimb_t rhat = num % vn[n - 1];
        while (qhat >= ((dlimb_t)1 << 32) ||
               qhat * vn[n - 2] > ((rhat << 32) | un[j + n - 2])) {
            qhat--; rhat += vn[n - 1];
            if (rhat >= ((dlimb_t)1 << 32)) break;
        }
        /* multiply and subtract */
        int64_t k = 0, t;
        for (int i = 0; i < n; i++) {
            dlimb_t p = qhat * vn[i];
            G.word_muls++;
            t = (int64_t)un[j + i] - k - (int64_t)(p & 0xffffffffULL);
            un[j + i] = (limb_t)t;
            k = (int64_t)(p >> 32) - (t >> 32);
        }
        t = (int64_t)un[j + n] - k;
        un[j + n] = (limb_t)t;
        qq.d[j] = (limb_t)qhat;
        if (t < 0) {                      /* add back */
            qq.d[j]--;
            dlimb_t c = 0;
            for (int i = 0; i < n; i++) {
                dlimb_t sum = (dlimb_t)un[j + i] + vn[i] + c;
                un[j + i] = (limb_t)sum; c = sum >> 32;
            }
            un[j + n] = (limb_t)(un[j + n] + c);
        }
    }
    qq.n = m + 1; bn_norm(&qq);
    if (q) bn_copy(q, &qq);
    if (r) {
        bn_zero(r);
        if (s) {
            for (int i = 0; i < n; i++)
                r->d[i] = (un[i] >> s) | ((dlimb_t)un[i + 1] << (LIMB_BITS - s));
        } else {
            for (int i = 0; i < n; i++) r->d[i] = un[i];
        }
        r->n = n; bn_norm(r);
    }
}

static void bn_mod(bn *r, const bn *a, const bn *m) { bn_divmod(NULL, r, a, m); }

/* r = (a * b) mod m */
static void bn_modmul(bn *r, const bn *a, const bn *b, const bn *m) {
    bn p; bn_mul(&p, a, b); bn_mod(r, &p, m);
    G.mod_muls++;
}
/* r = (a * a) mod m */
static void bn_modsqr(bn *r, const bn *a, const bn *m) {
    bn p; bn_sqr(&p, a); bn_mod(r, &p, m);
    G.mod_sqrs++;
}

/* r = base^exp mod m  (right-to-left square-and-multiply) */
static void bn_modpow(bn *r, const bn *base, const bn *exp, const bn *m) {
    bn result, b;
    bn_set_u32(&result, 1);
    bn_mod(&b, base, m);
    int bits = bn_bitlen(exp);
    for (int i = 0; i < bits; i++) {
        if (bn_bit(exp, i)) bn_modmul(&result, &result, &b, m);
        if (i + 1 < bits)   bn_modsqr(&b, &b, m);
    }
    bn_copy(r, &result);
}

/* gcd via Euclid */
static void bn_gcd(bn *r, const bn *a, const bn *b) {
    bn x, y, t;
    bn_copy(&x, a); bn_copy(&y, b);
    while (!bn_is_zero(&y)) {
        bn_mod(&t, &x, &y);
        bn_copy(&x, &y); bn_copy(&y, &t);
    }
    bn_copy(r, &x);
}

/* Jacobi symbol (a/n), n odd > 0; returns -1, 0, or 1.  steps counted. */
static int bn_jacobi(const bn *a_in, const bn *n_in, uint64_t *steps) {
    bn a, n, t;
    bn_mod(&a, a_in, n_in);
    bn_copy(&n, n_in);
    int result = 1;
    while (!bn_is_zero(&a)) {
        while (bn_is_even(&a)) {
            bn_shr(&a, &a, 1);
            uint32_t r8 = (n.n ? n.d[0] : 0) & 7u;
            if (r8 == 3 || r8 == 5) result = -result;
            if (steps) (*steps)++;
        }
        /* swap a, n */
        bn_copy(&t, &a); bn_copy(&a, &n); bn_copy(&n, &t);
        uint32_t a4 = (a.n ? a.d[0] : 0) & 3u;
        uint32_t n4 = (n.n ? n.d[0] : 0) & 3u;
        if (a4 == 3 && n4 == 3) result = -result;
        bn_mod(&a, &a, &n);
        if (steps) (*steps)++;
    }
    return bn_is_one(&n) ? result : 0;
}

/* integer b-th root, floored: returns r = floor(a^(1/b)); sets *exact */
static void bn_iroot(bn *root, const bn *a, int b, int *exact) {
    if (bn_is_zero(a)) { bn_zero(root); if (exact) *exact = 1; return; }
    int bits = bn_bitlen(a);
    int hi_bits = bits / b + 1;
    bn lo, hi, mid, p, one;
    bn_set_u32(&one, 1);
    bn_set_u32(&lo, 1);
    bn_shl(&hi, &one, hi_bits);            /* hi = 2^hi_bits  (upper bound) */
    /* binary search for largest mid with mid^b <= a */
    while (bn_cmp(&lo, &hi) < 0) {
        bn sum; bn_add(&sum, &lo, &hi);
        bn_add_small(&sum, &sum, 1);
        bn_shr(&mid, &sum, 1);             /* mid = ceil((lo+hi+1)/2) */
        /* p = mid^b */
        bn_copy(&p, &mid);
        int overflow = 0;
        for (int e = 1; e < b; e++) {
            bn_mul(&p, &p, &mid);
            if (bn_bitlen(&p) > MAXBITS - 64) { overflow = 1; break; }
        }
        if (!overflow && bn_cmp(&p, a) <= 0) bn_copy(&lo, &mid);
        else { bn_sub_small(&hi, &mid, 1); }
    }
    bn_copy(root, &lo);
    /* check exactness */
    bn_copy(&p, &lo);
    for (int e = 1; e < b; e++) bn_mul(&p, &p, &lo);
    if (exact) *exact = (bn_cmp(&p, a) == 0);
}

/* ===================================================================== */
/*  Parsing and printing                                                 */
/* ===================================================================== */

static int bn_from_dec(bn *a, const char *s) {
    bn_zero(a);
    int any = 0;
    for (; *s; s++) {
        if (isspace((unsigned char)*s)) continue;
        if (!isdigit((unsigned char)*s)) return 0;
        bn_mul_small(a, a, 10);
        bn_add_small(a, a, (uint32_t)(*s - '0'));
        any = 1;
        if (a->n > MAXLIMBS - 4) return 0;
    }
    return any;
}

static int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static int bn_from_hex(bn *a, const char *s) {
    bn_zero(a);
    int any = 0;
    for (; *s; s++) {
        if (isspace((unsigned char)*s)) continue;
        int v = hexval(*s);
        if (v < 0) return 0;
        bn_shl(a, a, 4);
        bn_add_small(a, a, (uint32_t)v);
        any = 1;
        if (a->n > MAXLIMBS - 4) return 0;
    }
    return any;
}

/* parse decimal or 0xHEX */
static int bn_parse(bn *a, const char *s) {
    while (isspace((unsigned char)*s)) s++;
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))
        return bn_from_hex(a, s + 2);
    return bn_from_dec(a, s);
}

/* decimal string into caller buffer; returns buf */
static char *bn_to_dec(const bn *a, char *buf, size_t bufsz) {
    if (bn_is_zero(a)) { snprintf(buf, bufsz, "0"); return buf; }
    char tmp[MAXBITS / 3 + 16];
    int pos = 0;
    bn cur; bn_copy(&cur, a);
    bn q; bn billion; bn_set_u32(&billion, 1000000000u);
    while (!bn_is_zero(&cur)) {
        bn rem;
        bn_divmod(&q, &rem, &cur, &billion);
        uint32_t chunk = rem.n ? rem.d[0] : 0;
        bn_copy(&cur, &q);
        for (int i = 0; i < 9; i++) {
            tmp[pos++] = (char)('0' + chunk % 10);
            chunk /= 10;
            if (bn_is_zero(&cur) && chunk == 0) break;
        }
    }
    int j = 0;
    for (int i = pos - 1; i >= 0 && j < (int)bufsz - 1; i--) buf[j++] = tmp[i];
    buf[j] = '\0';
    return buf;
}

/* ===================================================================== */
/*  PRNG (splitmix64) for witness bases                                  */
/* ===================================================================== */

static uint64_t prng_state;
static void prng_seed(uint64_t s) { prng_state = s; }
static uint64_t prng_next(void) {
    uint64_t z = (prng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* random a in [2, n-2], requires n > 4 */
static void bn_rand_base(bn *a, const bn *n) {
    bn nm3, r;
    bn_sub_small(&nm3, n, 3);              /* n-3 */
    int bits = bn_bitlen(n);
    int limbs = (bits + 31) / 32;
    bn_zero(&r);
    for (int i = 0; i < limbs; i++) {
        if ((i & 1) == 0) {
            uint64_t v = prng_next();
            r.d[i] = (limb_t)v;
            if (i + 1 < limbs) r.d[i + 1] = (limb_t)(v >> 32);
        }
    }
    r.n = limbs; bn_norm(&r);
    bn rem; bn_mod(&rem, &r, &nm3);        /* in [0, n-4] */
    bn_add_small(a, &rem, 2);              /* in [2, n-2] */
}

/* ===================================================================== */
/*  Settings & verdict types                                             */
/* ===================================================================== */

typedef struct {
    int      k;                /* rounds for probabilistic tests        */
    uint64_t trial_bound;      /* trial division ceiling                */
    int      aks_cap_bits;     /* skip AKS above this bit length        */
    int      aks_hard_cap;     /* never allow AKS above this            */
    uint64_t seed;             /* PRNG seed                             */
} settings;

typedef enum {
    V_PRIME, V_COMPOSITE, V_PROBABLY_PRIME, V_INCONCLUSIVE, V_NA, V_SKIPPED
} verdict;

static const char *verdict_str(verdict v) {
    switch (v) {
        case V_PRIME:           return "PRIME (proven)";
        case V_COMPOSITE:       return "COMPOSITE";
        case V_PROBABLY_PRIME:  return "PROBABLY PRIME";
        case V_INCONCLUSIVE:    return "INCONCLUSIVE";
        case V_NA:              return "N/A";
        case V_SKIPPED:         return "SKIPPED";
    }
    return "?";
}

/* short, CSV-safe verdict token */
static const char *verdict_token(verdict v) {
    switch (v) {
        case V_PRIME:          return "PRIME";
        case V_COMPOSITE:      return "COMPOSITE";
        case V_PROBABLY_PRIME: return "PROBABLY_PRIME";
        case V_INCONCLUSIVE:   return "INCONCLUSIVE";
        case V_NA:             return "NA";
        case V_SKIPPED:        return "SKIPPED";
    }
    return "UNKNOWN";
}

typedef struct {
    const char *method;
    verdict     v;
    costs       cost;
    char        hilabel[40];   /* high-level cost label                 */
    uint64_t    hival;         /* high-level cost value                 */
    int         applicable;    /* included in run-all table             */
} result_row;

static void print_costs(const costs *c) {
    printf("    cost (low)  : word_muls=%llu  word_div_steps=%llu\n",
           (unsigned long long)c->word_muls,
           (unsigned long long)c->word_div_steps);
    printf("    cost (mid)  : mod_muls=%llu  mod_sqrs=%llu  reductions=%llu\n",
           (unsigned long long)c->mod_muls,
           (unsigned long long)c->mod_sqrs,
           (unsigned long long)c->mod_reductions);
}

/* ===================================================================== */
/*  Small deterministic prime check (for AKS exponent / sanity)          */
/* ===================================================================== */

static int small_is_prime(uint64_t x) {
    if (x < 2) return 0;
    if (x % 2 == 0) return x == 2;
    for (uint64_t d = 3; d * d <= x; d += 2)
        if (x % d == 0) return 0;
    return 1;
}

/* ===================================================================== */
/*  1. Trial division                                                    */
/* ===================================================================== */

static verdict test_trial(const bn *n, const settings *cfg, result_row *row) {
    costs_reset();
    uint64_t divisions = 0;
    verdict v;
    char factbuf[64] = "";

    bn two; bn_set_u32(&two, 2);
    if (bn_cmp(n, &two) < 0) { v = V_COMPOSITE; goto done; }
    if (bn_cmp(n, &two) == 0) { v = V_PRIME; goto done; }
    if (bn_is_even(n)) { v = V_COMPOSITE; strcpy(factbuf, "2"); goto done; }

    {
        int proven = 0;
        uint64_t d = 3;
        for (; d <= cfg->trial_bound; d += 2) {
            /* d*d > n  ->  proven prime */
            bn dd; bn_set_u64(&dd, d * d);
            if (bn_cmp(&dd, n) > 0) { proven = 1; break; }
            bn rem; bn small; bn_set_u64(&small, d);
            bn_mod(&rem, n, &small);
            divisions++;
            if (bn_is_zero(&rem)) {
                u64str(d, factbuf);
                v = V_COMPOSITE; goto done;
            }
        }
        v = proven ? V_PRIME : V_INCONCLUSIVE;
    }

done:
    printf("[1] Trial division\n");
    printf("    verdict     : %s", verdict_str(v));
    if (factbuf[0]) printf("  (factor found: %s)", factbuf);
    else if (v == V_INCONCLUSIVE)
        printf("  (no factor <= %llu; sqrt(n) exceeds bound)",
               (unsigned long long)cfg->trial_bound);
    printf("\n");
    printf("    high level  : trial_divisions=%llu\n", (unsigned long long)divisions);
    print_costs(&G);
    if (row) {
        row->method = "Trial division"; row->v = v; row->cost = G;
        snprintf(row->hilabel, sizeof row->hilabel, "divisions");
        row->hival = divisions; row->applicable = 1;
    }
    return v;
}

/* ===================================================================== */
/*  2. Fermat                                                            */
/* ===================================================================== */

static verdict test_fermat(const bn *n, const settings *cfg, result_row *row) {
    costs_reset();
    prng_seed(cfg->seed);
    verdict v = V_PROBABLY_PRIME;
    int rounds = 0;
    char note[96] = "";

    bn two; bn_set_u32(&two, 2);
    bn four; bn_set_u32(&four, 4);
    if (bn_cmp(n, &two) < 0)  { v = V_COMPOSITE; goto done; }
    if (bn_cmp(n, &two) == 0) { v = V_PRIME; goto done; }
    if (bn_is_even(n))        { v = V_COMPOSITE; strcpy(note, "even"); goto done; }
    if (bn_cmp(n, &four) <= 0){ v = V_PRIME; goto done; }   /* n == 3 */

    {
        bn nm1; bn_sub_small(&nm1, n, 1);
        printf("[2] Fermat   (k = %d rounds, seed = %llu)\n",
               cfg->k, (unsigned long long)cfg->seed);
        for (int i = 0; i < cfg->k; i++) {
            bn a; bn_rand_base(&a, n);
            rounds++;
            bn g; bn_gcd(&g, &a, n);
            char abuf[64]; bn_to_dec(&a, abuf, sizeof abuf);
            if (!bn_is_one(&g)) {
                char gb[64]; bn_to_dec(&g, gb, sizeof gb);
                printf("      base %-20s -> gcd(a,n)=%s  COMPOSITE\n", abuf, gb);
                snprintf(note, sizeof note, "gcd witness");
                v = V_COMPOSITE; goto done;
            }
            bn r; bn_modpow(&r, &a, &nm1, n);
            if (!bn_is_one(&r)) {
                printf("      base %-20s -> a^(n-1) != 1   COMPOSITE\n", abuf);
                snprintf(note, sizeof note, "Fermat witness");
                v = V_COMPOSITE; goto done;
            }
            printf("      base %-20s -> a^(n-1) = 1     pass\n", abuf);
        }
        v = V_PROBABLY_PRIME;
        goto report;
    }

done:
    printf("[2] Fermat   (k = %d rounds, seed = %llu)\n",
           cfg->k, (unsigned long long)cfg->seed);
report:
    printf("    verdict     : %s", verdict_str(v));
    if (v == V_PROBABLY_PRIME) printf("  (passed all %d bases; Carmichael caveat)", rounds);
    if (note[0]) printf("  [%s]", note);
    printf("\n");
    printf("    high level  : rounds=%d\n", rounds);
    print_costs(&G);
    if (row) {
        row->method = "Fermat"; row->v = v; row->cost = G;
        snprintf(row->hilabel, sizeof row->hilabel, "rounds");
        row->hival = (uint64_t)rounds; row->applicable = 1;
    }
    return v;
}

/* ===================================================================== */
/*  3. Solovay-Strassen                                                  */
/* ===================================================================== */

static verdict test_solovay(const bn *n, const settings *cfg, result_row *row) {
    costs_reset();
    prng_seed(cfg->seed);
    verdict v = V_PROBABLY_PRIME;
    int rounds = 0;
    uint64_t jsteps = 0;
    char note[96] = "";

    bn two; bn_set_u32(&two, 2);
    bn four; bn_set_u32(&four, 4);
    if (bn_cmp(n, &two) < 0)  { v = V_COMPOSITE; goto done; }
    if (bn_cmp(n, &two) == 0) { v = V_PRIME; goto done; }
    if (bn_is_even(n))        { v = V_COMPOSITE; strcpy(note, "even"); goto done; }
    if (bn_cmp(n, &four) <= 0){ v = V_PRIME; goto done; }

    {
        bn nm1; bn_sub_small(&nm1, n, 1);
        bn e;   bn_shr(&e, &nm1, 1);       /* (n-1)/2 */
        printf("[3] Solovay-Strassen   (k = %d rounds, seed = %llu)\n",
               cfg->k, (unsigned long long)cfg->seed);
        for (int i = 0; i < cfg->k; i++) {
            bn a; bn_rand_base(&a, n);
            rounds++;
            char abuf[64]; bn_to_dec(&a, abuf, sizeof abuf);
            bn g; bn_gcd(&g, &a, n);
            if (!bn_is_one(&g)) {
                printf("      base %-20s -> gcd(a,n)>1      COMPOSITE\n", abuf);
                snprintf(note, sizeof note, "gcd witness");
                v = V_COMPOSITE; goto done;
            }
            int j = bn_jacobi(&a, n, &jsteps);        /* -1 or 1 (0 ruled out) */
            bn x; bn_modpow(&x, &a, &e, n);           /* a^((n-1)/2) mod n */
            /* expected: x == j (mod n): 1 -> 1, -1 -> n-1 */
            int ok;
            if (j == 1) ok = bn_is_one(&x);
            else /* j == -1 */ ok = (bn_cmp(&x, &nm1) == 0);
            if (!ok) {
                printf("      base %-20s -> Euler mismatch   COMPOSITE\n", abuf);
                snprintf(note, sizeof note, "Euler witness");
                v = V_COMPOSITE; goto done;
            }
            printf("      base %-20s -> (a/n)=%+d, Euler ok  pass\n", abuf, j);
        }
        v = V_PROBABLY_PRIME;
        goto report;
    }

done:
    printf("[3] Solovay-Strassen   (k = %d rounds, seed = %llu)\n",
           cfg->k, (unsigned long long)cfg->seed);
report:
    printf("    verdict     : %s", verdict_str(v));
    if (v == V_PROBABLY_PRIME) printf("  (error <= 2^-%d)", rounds);
    if (note[0]) printf("  [%s]", note);
    printf("\n");
    printf("    high level  : rounds=%d  jacobi_steps=%llu\n",
           rounds, (unsigned long long)jsteps);
    print_costs(&G);
    if (row) {
        row->method = "Solovay-Strassen"; row->v = v; row->cost = G;
        snprintf(row->hilabel, sizeof row->hilabel, "rounds");
        row->hival = (uint64_t)rounds; row->applicable = 1;
    }
    return v;
}

/* ===================================================================== */
/*  4. Miller-Rabin                                                      */
/* ===================================================================== */

static verdict test_miller(const bn *n, const settings *cfg, result_row *row) {
    costs_reset();
    prng_seed(cfg->seed);
    verdict v = V_PROBABLY_PRIME;
    int rounds = 0, s = 0;
    uint64_t squarings = 0;
    char note[96] = "";

    bn two; bn_set_u32(&two, 2);
    bn four; bn_set_u32(&four, 4);
    if (bn_cmp(n, &two) < 0)  { v = V_COMPOSITE; goto done; }
    if (bn_cmp(n, &two) == 0) { v = V_PRIME; goto done; }
    if (bn_is_even(n))        { v = V_COMPOSITE; strcpy(note, "even"); goto done; }
    if (bn_cmp(n, &four) <= 0){ v = V_PRIME; goto done; }

    {
        bn nm1; bn_sub_small(&nm1, n, 1);
        bn d; bn_copy(&d, &nm1);
        while (bn_is_even(&d)) { bn_shr(&d, &d, 1); s++; }
        printf("[4] Miller-Rabin   (k = %d rounds, seed = %llu)\n",
               cfg->k, (unsigned long long)cfg->seed);
        printf("    decomposition: n-1 = 2^%d * d\n", s);
        for (int i = 0; i < cfg->k; i++) {
            bn a; bn_rand_base(&a, n);
            rounds++;
            char abuf[64]; bn_to_dec(&a, abuf, sizeof abuf);
            bn x; bn_modpow(&x, &a, &d, n);
            if (bn_is_one(&x) || bn_cmp(&x, &nm1) == 0) {
                printf("      base %-20s -> pass\n", abuf);
                continue;
            }
            int found = 0;
            for (int r = 1; r < s; r++) {
                bn_modsqr(&x, &x, n); squarings++;
                if (bn_cmp(&x, &nm1) == 0) { found = 1; break; }
            }
            if (!found) {
                printf("      base %-20s -> no -1 in chain  COMPOSITE\n", abuf);
                snprintf(note, sizeof note, "MR witness");
                v = V_COMPOSITE; goto done;
            }
            printf("      base %-20s -> pass\n", abuf);
        }
        v = V_PROBABLY_PRIME;
        goto report;
    }

done:
    printf("[4] Miller-Rabin   (k = %d rounds, seed = %llu)\n",
           cfg->k, (unsigned long long)cfg->seed);
report:
    printf("    verdict     : %s", verdict_str(v));
    if (v == V_PROBABLY_PRIME) printf("  (error <= 4^-%d)", rounds);
    if (note[0]) printf("  [%s]", note);
    printf("\n");
    printf("    high level  : rounds=%d  extra_squarings=%llu\n",
           rounds, (unsigned long long)squarings);
    print_costs(&G);
    if (row) {
        row->method = "Miller-Rabin"; row->v = v; row->cost = G;
        snprintf(row->hilabel, sizeof row->hilabel, "rounds");
        row->hival = (uint64_t)rounds; row->applicable = 1;
    }
    return v;
}

/* ===================================================================== */
/*  5. AKS                                                               */
/* ===================================================================== */

/* Euler totient of small r */
static uint64_t totient(uint64_t r) {
    uint64_t result = r, m = r;
    for (uint64_t p = 2; p * p <= m; p++) {
        if (m % p == 0) {
            while (m % p == 0) m /= p;
            result -= result / p;
        }
    }
    if (m > 1) result -= result / m;
    return result;
}

/* multiplicative order of (n mod r) modulo r, or 0 if gcd(n,r) != 1 */
static uint64_t mult_order(const bn *n, uint64_t r) {
    bn rr; bn_set_u64(&rr, r);
    bn nm; bn_mod(&nm, n, &rr);
    uint64_t a = nm.n ? ((uint64_t)nm.d[0] | ((nm.n > 1) ? ((uint64_t)nm.d[1] << 32) : 0)) : 0;
    a %= r;
    /* gcd */
    uint64_t x = a, y = r;
    while (y) { uint64_t t = x % y; x = y; y = t; }
    if (x != 1) return 0;
    uint64_t k = 1, cur = a % r;
    while (cur != 1) { cur = (cur * (a % r)) % r; k++; if (k > r) return 0; }
    return k;
}

/* poly arithmetic over Z_n[x]/(x^r - 1); coeff arrays length r */
static void poly_mulmod(bn *res, const bn *A, const bn *B, int r, const bn *n,
                        uint64_t *polymuls) {
    bn *acc = (bn *)malloc(sizeof(bn) * r);
    for (int i = 0; i < r; i++) bn_zero(&acc[i]);
    for (int i = 0; i < r; i++) {
        if (bn_is_zero(&A[i])) continue;
        for (int j = 0; j < r; j++) {
            if (bn_is_zero(&B[j])) continue;
            int kk = (i + j) % r;
            bn prod; bn_modmul(&prod, &A[i], &B[j], n);
            bn_add(&acc[kk], &acc[kk], &prod);
            bn_mod(&acc[kk], &acc[kk], n);
            (*polymuls)++;
        }
    }
    for (int i = 0; i < r; i++) bn_copy(&res[i], &acc[i]);
    free(acc);
}

static verdict test_aks(const bn *n, const settings *cfg, result_row *row) {
    costs_reset();
    uint64_t polymuls = 0;
    int bits = bn_bitlen(n);
    verdict v;
    char note[160] = "";
    uint64_t r_used = 0;

    printf("[5] AKS\n");
    if (bits > cfg->aks_cap_bits) {
        v = V_SKIPPED;
        snprintf(note, sizeof note,
                 "n is %d bits > AKS cap (%d). AKS is impractical here.",
                 bits, cfg->aks_cap_bits);
        goto report;
    }

    bn two; bn_set_u32(&two, 2);
    if (bn_cmp(n, &two) < 0)  { v = V_COMPOSITE; goto report; }
    if (bn_cmp(n, &two) == 0) { v = V_PRIME; goto report; }

    /* Step 1: perfect power? */
    for (int b = 2; b <= bits; b++) {
        int exact; bn root; bn_iroot(&root, n, b, &exact);
        if (exact && bn_cmp(&root, &two) >= 0) {
            char rb[64]; bn_to_dec(&root, rb, sizeof rb);
            snprintf(note, sizeof note, "perfect power: %s^%d", rb, b);
            v = V_COMPOSITE; goto report;
        }
    }

    /* Step 2: find r with ord_r(n) > (log2 n)^2 */
    {
        uint64_t target = (uint64_t)bits * (uint64_t)bits;   /* >= ceil(log2 n)^2 */
        uint64_t r = 2;
        for (;; r++) {
            bn rr; bn_set_u64(&rr, r);
            if (bn_cmp(n, &rr) == 0) continue;
            uint64_t ord = mult_order(n, r);
            if (ord == 0) continue;        /* gcd(n,r) != 1: handled in step 3 */
            if (ord > target) break;
        }
        r_used = r;

        /* Step 3: gcd check for a = 2..min(r, n-1) */
        bn nm1; bn_sub_small(&nm1, n, 1);
        uint64_t amax = r;
        {
            bn rb; bn_set_u64(&rb, r);
            if (bn_cmp(&nm1, &rb) < 0) amax = nm1.n ? nm1.d[0] : 0;
        }
        for (uint64_t a = 2; a <= amax; a++) {
            bn ab; bn_set_u64(&ab, a);
            bn g; bn_gcd(&g, &ab, n);
            if (!bn_is_one(&g) && bn_cmp(&g, n) != 0) {
                char ab2[24];
                snprintf(note, sizeof note, "small factor %s (gcd)", u64str(a, ab2));
                v = V_COMPOSITE; goto report;
            }
        }

        /* Step 4: n <= r  ->  prime */
        bn rb; bn_set_u64(&rb, r);
        if (bn_cmp(n, &rb) <= 0) { v = V_PRIME; goto report; }

        /* Step 5: polynomial congruences */
        uint64_t phir = totient(r);
        /* limit = floor( sqrt(phi(r)) * log2(n) ); use bits as log2 upper bound */
        uint64_t sq = (uint64_t)(/*floor sqrt*/ 0);
        { uint64_t x = phir; while ((sq + 1) * (sq + 1) <= x) sq++; }
        uint64_t limit = sq * (uint64_t)bits;
        if (limit < 1) limit = 1;

        bn nmodr; { bn t; bn_set_u64(&t, r); bn_mod(&t, n, &t);
                    bn_copy(&nmodr, &t); }
        uint64_t nmr = nmodr.n ? nmodr.d[0] : 0;

        printf("    r=%llu  phi(r)=%llu  poly_checks(a)=1..%llu\n",
               (unsigned long long)r, (unsigned long long)phir,
               (unsigned long long)limit);

        bn *base = (bn *)malloc(sizeof(bn) * r);
        bn *res  = (bn *)malloc(sizeof(bn) * r);
        int composite = 0;
        for (uint64_t a = 1; a <= limit && !composite; a++) {
            /* base poly = x + a   (mod n) */
            for (int i = 0; i < (int)r; i++) bn_zero(&base[i]);
            bn amod; bn_set_u64(&amod, a); bn_mod(&amod, &amod, n);
            bn_copy(&base[0], &amod);
            if (r > 1) bn_set_u32(&base[1 % (int)r], 1);
            /* res = base^n  mod (x^r - 1, n)  via square-and-multiply */
            for (int i = 0; i < (int)r; i++) bn_zero(&res[i]);
            bn_set_u32(&res[0], 1);
            int nb = bn_bitlen(n);
            bn cur; (void)cur;
            bn *acc = (bn *)malloc(sizeof(bn) * r);
            for (int i = 0; i < (int)r; i++) bn_copy(&acc[i], &base[i]);
            for (int bit = 0; bit < nb; bit++) {
                if (bn_bit(n, bit)) poly_mulmod(res, res, acc, r, n, &polymuls);
                if (bit + 1 < nb)   poly_mulmod(acc, acc, acc, r, n, &polymuls);
            }
            free(acc);
            /* expected = x^(n mod r) + a  (mod n) */
            bn *exp = (bn *)malloc(sizeof(bn) * r);
            for (int i = 0; i < (int)r; i++) bn_zero(&exp[i]);
            bn_set_u32(&exp[nmr % (int)r], 1);
            bn av; bn_set_u64(&av, a);
            bn_add(&exp[0], &exp[0], &av); bn_mod(&exp[0], &exp[0], n);
            for (int i = 0; i < (int)r; i++)
                if (bn_cmp(&res[i], &exp[i]) != 0) { composite = 1; break; }
            free(exp);
        }
        free(base); free(res);
        v = composite ? V_COMPOSITE : V_PRIME;
        if (composite) snprintf(note, sizeof note, "poly congruence failed");
    }

report:
    printf("    verdict     : %s", verdict_str(v));
    if (note[0]) printf("  [%s]", note);
    printf("\n");
    if (v != V_SKIPPED) {
        printf("    high level  : r=%llu  poly_mults=%llu\n",
               (unsigned long long)r_used, (unsigned long long)polymuls);
        print_costs(&G);
    }
    if (row) {
        row->method = "AKS"; row->v = v; row->cost = G;
        snprintf(row->hilabel, sizeof row->hilabel, "poly_mults");
        row->hival = polymuls; row->applicable = (v != V_SKIPPED);
    }
    return v;
}

/* ===================================================================== */
/*  6. Lucas-Lehmer  (Mersenne numbers only)                             */
/* ===================================================================== */

/* if n == 2^p - 1 return p (>0), else 0 */
static int mersenne_exponent(const bn *n) {
    bn np1; bn_add_small(&np1, n, 1);      /* n+1 must be a power of two */
    if (np1.n == 0) return 0;
    /* exactly one bit set? */
    int ones = 0, pos = -1;
    for (int i = 0; i < bn_bitlen(&np1); i++)
        if (bn_bit(&np1, i)) { ones++; pos = i; }
    if (ones != 1) return 0;
    return pos;                            /* n = 2^pos - 1 */
}

/* reduce x mod (2^p - 1), x may be up to ~2p bits */
static void mersenne_reduce(bn *x, int p, const bn *M) {
    while (bn_bitlen(x) > p) {
        bn hi, lo, mask;
        bn_shr(&hi, x, p);                 /* hi = x >> p */
        /* lo = x & (2^p - 1) */
        bn_copy(&mask, M);
        /* low p bits: copy then clear above */
        bn_copy(&lo, x);
        {
            int limbs = p / LIMB_BITS, bitsr = p % LIMB_BITS;
            if (bitsr == 0) { lo.n = limbs; }
            else {
                lo.n = limbs + 1;
                if (lo.n > x->n) lo.n = x->n;
                if (limbs < lo.n) lo.d[limbs] &= ((limb_t)1 << bitsr) - 1;
            }
            bn_norm(&lo);
        }
        bn_add(x, &hi, &lo);
    }
    if (bn_cmp(x, M) >= 0) bn_sub(x, x, M);
}

static verdict test_lucas(const bn *n, const settings *cfg, result_row *row) {
    (void)cfg;
    costs_reset();
    verdict v;
    char note[120] = "";
    uint64_t iters = 0;
    int p = mersenne_exponent(n);

    printf("[6] Lucas-Lehmer\n");
    if (p == 0) {
        v = V_NA;
        snprintf(note, sizeof note, "n is not of the form 2^p - 1 (not Mersenne)");
        goto report;
    }
    printf("    n = 2^%d - 1  (Mersenne number M_%d)\n", p, p);
    if (p == 2) { v = V_PRIME; snprintf(note, sizeof note, "M_2 = 3"); goto report; }
    if (!small_is_prime((uint64_t)p)) {
        /* M_p with composite p is always composite */
        v = V_COMPOSITE;
        snprintf(note, sizeof note, "exponent p=%d is composite => M_p composite", p);
        goto report;
    }

    {
        bn M; bn_copy(&M, n);
        bn s; bn_set_u32(&s, 4);
        for (int i = 0; i < p - 2; i++) {
            bn sq; bn_sqr(&sq, &s);        /* s^2 */
            mersenne_reduce(&sq, p, &M);
            /* s = sq - 2 (mod M): if sq < 2 borrow */
            if (bn_cmp(&sq, &s /*dummy*/) >= 0) { /* always, sq>=0 */ }
            bn two; bn_set_u32(&two, 2);
            if (bn_cmp(&sq, &two) < 0) { bn_add(&sq, &sq, &M); }
            bn_sub(&s, &sq, &two);
            iters++;
        }
        v = bn_is_zero(&s) ? V_PRIME : V_COMPOSITE;
    }

report:
    printf("    verdict     : %s", verdict_str(v));
    if (note[0]) printf("  [%s]", note);
    printf("\n");
    if (v == V_PRIME || v == V_COMPOSITE) {
        if (iters) printf("    high level  : LL_iterations=%llu (= p-2)\n",
                          (unsigned long long)iters);
        print_costs(&G);
    }
    if (row) {
        row->method = "Lucas-Lehmer"; row->v = v; row->cost = G;
        snprintf(row->hilabel, sizeof row->hilabel, "LL_iters");
        row->hival = iters; row->applicable = (v == V_PRIME || v == V_COMPOSITE);
    }
    return v;
}

/* ===================================================================== */
/*  Run all + comparison table                                           */
/* ===================================================================== */

static void run_all(const bn *n, const settings *cfg) {
    result_row rows[6];
    memset(rows, 0, sizeof rows);
    printf("\n================ RUN ALL ================\n\n");
    test_trial  (n, cfg, &rows[0]); printf("\n");
    test_fermat (n, cfg, &rows[1]); printf("\n");
    test_solovay(n, cfg, &rows[2]); printf("\n");
    test_miller (n, cfg, &rows[3]); printf("\n");
    test_aks    (n, cfg, &rows[4]); printf("\n");
    test_lucas  (n, cfg, &rows[5]); printf("\n");

    printf("================ COMPARISON ================\n");
    printf("%-18s %-16s %14s %12s  %s\n",
           "method", "verdict", "word_muls", "mod_mul+sqr", "high-level");
    printf("--------------------------------------------------------------------------------\n");
    for (int i = 0; i < 6; i++) {
        result_row *r = &rows[i];
        char hl[64];
        if (r->applicable) {
            char vb[24];
            snprintf(hl, sizeof hl, "%s %s", u64str(r->hival, vb), r->hilabel);
        } else {
            snprintf(hl, sizeof hl, "-");
        }
        printf("%-18s %-16s %14llu %12llu  %s\n",
               r->method ? r->method : "?",
               verdict_str(r->v),
               (unsigned long long)r->cost.word_muls,
               (unsigned long long)(r->cost.mod_muls + r->cost.mod_sqrs),
               hl);
    }
    printf("--------------------------------------------------------------------------------\n");
}

/* ===================================================================== */
/*  Menu / main                                                          */
/* ===================================================================== */

static void read_line(char *buf, size_t sz) {
    if (!fgets(buf, (int)sz, stdin)) { buf[0] = '\0'; return; }
    size_t L = strlen(buf);
    while (L && (buf[L - 1] == '\n' || buf[L - 1] == '\r')) buf[--L] = '\0';
}

static void settings_menu(settings *cfg) {
    char line[128];
    for (;;) {
        printf("\n--- Settings ---\n");
        printf("  1) rounds k          = %d\n", cfg->k);
        printf("  2) trial_bound       = %llu\n", (unsigned long long)cfg->trial_bound);
        printf("  3) aks_cap_bits      = %d   (hard max %d)\n",
               cfg->aks_cap_bits, cfg->aks_hard_cap);
        printf("  4) seed              = %llu\n", (unsigned long long)cfg->seed);
        printf("  0) back\n");
        printf("  choice> ");
        read_line(line, sizeof line);
        if (line[0] == '0' || line[0] == '\0') break;
        printf("  new value> ");
        char val[128]; read_line(val, sizeof val);
        switch (line[0]) {
            case '1': { int v = atoi(val); if (v > 0) cfg->k = v; } break;
            case '2': { unsigned long long v = strtoull(val, NULL, 10);
                        if (v >= 2) cfg->trial_bound = v; } break;
            case '3': { int v = atoi(val);
                        if (v >= 1 && v <= cfg->aks_hard_cap) cfg->aks_cap_bits = v;
                        else printf("  (must be 1..%d)\n", cfg->aks_hard_cap); } break;
            case '4': { cfg->seed = strtoull(val, NULL, 10); } break;
            default: break;
        }
    }
}

/* ===================================================================== */
/*  Batch / CSV mode (non-interactive, for generating article datasets)  */
/* ===================================================================== */

enum { M_TRIAL=1, M_FERMAT=2, M_SOLOVAY=4, M_MILLER=8, M_AKS=16, M_LUCAS=32,
       M_ALL=63 };

/* CSV row -> the real stdout stream (passed in as csv). */
static void emit_csv(FILE *csv, const char *label, int bits,
                     const result_row *r, const settings *cfg) {
    char wm[24], wd[24], mm[24], ms[24], md[24], hv[24], sd[24];
    fprintf(csv, "%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%s\n",
            label[0] ? label : "-", bits,
            r->method ? r->method : "?", verdict_token(r->v),
            u64str(r->cost.word_muls, wm),
            u64str(r->cost.word_div_steps, wd),
            u64str(r->cost.mod_muls, mm),
            u64str(r->cost.mod_sqrs, ms),
            u64str(r->cost.mod_reductions, md),
            r->hilabel[0] ? r->hilabel : "-",
            u64str(r->hival, hv),
            cfg->k, u64str(cfg->seed, sd));
}

static int starts_with(const char *s, const char *p) {
    return strncmp(s, p, strlen(p)) == 0;
}

/*
 * Batch mode: reads one entry per line of stdin: "<number>[,<label>]"
 * (number is decimal or 0xHEX). Emits one CSV row per selected method to
 * stdout. The verbose per-method report is sent to stderr (redirect/ignore).
 *
 * Flags: --csv  --method=all|trial|fermat|solovay|miller|aks|lucas
 *        --k=N  --seed=N  --trial-bound=N  --aks-cap=N
 */
static int run_batch(int argc, char **argv) {
    settings cfg;
    cfg.k = 20; cfg.trial_bound = 1000000ULL;
    cfg.aks_cap_bits = 20; cfg.aks_hard_cap = 40;
    cfg.seed = 88172645463325252ULL;
    int method = M_ALL;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if      (!strcmp(a, "--csv")) { /* mode flag */ }
        else if (starts_with(a, "--method=")) {
            const char *m = a + 9;
            if      (!strcmp(m, "all"))     method = M_ALL;
            else if (!strcmp(m, "trial"))   method = M_TRIAL;
            else if (!strcmp(m, "fermat"))  method = M_FERMAT;
            else if (!strcmp(m, "solovay")) method = M_SOLOVAY;
            else if (!strcmp(m, "miller"))  method = M_MILLER;
            else if (!strcmp(m, "aks"))     method = M_AKS;
            else if (!strcmp(m, "lucas"))   method = M_LUCAS;
        }
        else if (starts_with(a, "--k="))           cfg.k = atoi(a + 4);
        else if (starts_with(a, "--seed="))        cfg.seed = strtoull(a + 7, NULL, 10);
        else if (starts_with(a, "--trial-bound=")) cfg.trial_bound = strtoull(a + 14, NULL, 10);
        else if (starts_with(a, "--aks-cap="))     cfg.aks_cap_bits = atoi(a + 10);
    }

    /*
     * The six test functions print a verbose report to stdout. In batch
     * mode we want only CSV on the real stdout. So: dup the real stdout to
     * a fresh stream `csv`, then point stdout at the null device for the
     * verbose reports. CSV is written to `csv`.
     */
    fflush(stdout);
    int   saved_fd = _dup(_fileno(stdout));
    FILE *csv      = (saved_fd >= 0) ? _fdopen(saved_fd, "w") : NULL;
    if (!csv) csv = stdout;                 /* fallback: CSV+verbose mixed */
    else      freopen("nul", "w", stdout);  /* discard verbose reports     */

    fprintf(csv, "label,bits,method,verdict,word_muls,word_div_steps,"
                 "mod_muls,mod_sqrs,mod_reductions,high_label,high_value,k,seed\n");

    char line[8192];
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line);
        while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = '\0';
        if (line[0] == '\0' || line[0] == '#') continue;

        char label[256] = "";
        char *comma = strchr(line, ',');
        if (comma) { *comma = '\0'; strncpy(label, comma + 1, sizeof label - 1);
                     label[sizeof label - 1] = '\0'; }

        bn n;
        if (!bn_parse(&n, line)) { fprintf(csv, "%s,,PARSE_ERROR\n", label); continue; }
        int bits = bn_bitlen(&n);

        result_row row;
        if (method & M_TRIAL)   { test_trial  (&n, &cfg, &row); emit_csv(csv, label, bits, &row, &cfg); }
        if (method & M_FERMAT)  { test_fermat (&n, &cfg, &row); emit_csv(csv, label, bits, &row, &cfg); }
        if (method & M_SOLOVAY) { test_solovay(&n, &cfg, &row); emit_csv(csv, label, bits, &row, &cfg); }
        if (method & M_MILLER)  { test_miller (&n, &cfg, &row); emit_csv(csv, label, bits, &row, &cfg); }
        if (method & M_AKS)     { test_aks    (&n, &cfg, &row); emit_csv(csv, label, bits, &row, &cfg); }
        if (method & M_LUCAS)   { test_lucas  (&n, &cfg, &row); emit_csv(csv, label, bits, &row, &cfg); }
        fflush(csv);
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc > 1) return run_batch(argc, argv);

    settings cfg;
    cfg.k = 20;
    cfg.trial_bound = 1000000ULL;
    cfg.aks_cap_bits = 20;
    cfg.aks_hard_cap = 40;
    cfg.seed = 88172645463325252ULL;

    printf("=====================================================\n");
    printf(" Primality-test cost benchmark  (self-contained, C)\n");
    printf(" Methods: trial division, Fermat, Solovay-Strassen,\n");
    printf("          Miller-Rabin, AKS, Lucas-Lehmer\n");
    printf(" Cost = deterministic operation counts (machine-independent)\n");
    printf("=====================================================\n");

    char line[4096];
    bn n; int have = 0;

    for (;;) {
        if (!have) {
            printf("\nEnter number (decimal or 0xHEX, up to 2048 bits), or 'q' to quit:\n> ");
            read_line(line, sizeof line);
            if (line[0] == 'q' || line[0] == 'Q') break;
            if (line[0] == '\0') continue;
            if (!bn_parse(&n, line)) { printf("  ! could not parse a number.\n"); continue; }
            int bits = bn_bitlen(&n);
            if (bits > 2048) {
                printf("  ! %d bits exceeds the 2048-bit limit.\n", bits);
                continue;
            }
            char dbuf[800];
            printf("  parsed: %s\n", bn_to_dec(&n, dbuf, sizeof dbuf));
            printf("  bit length: %d\n", bits);
            have = 1;
        }

        printf("\n  Method:\n");
        printf("    1) Trial division     2) Fermat\n");
        printf("    3) Solovay-Strassen   4) Miller-Rabin\n");
        printf("    5) AKS                6) Lucas-Lehmer\n");
        printf("    7) RUN ALL            8) Settings\n");
        printf("    9) New number         0) Quit\n");
        printf("  choice> ");
        read_line(line, sizeof line);

        printf("\n");
        switch (line[0]) {
            case '1': test_trial  (&n, &cfg, NULL); break;
            case '2': test_fermat (&n, &cfg, NULL); break;
            case '3': test_solovay(&n, &cfg, NULL); break;
            case '4': test_miller (&n, &cfg, NULL); break;
            case '5': test_aks    (&n, &cfg, NULL); break;
            case '6': test_lucas  (&n, &cfg, NULL); break;
            case '7': run_all     (&n, &cfg);       break;
            case '8': settings_menu(&cfg);          break;
            case '9': have = 0;                     break;
            case '0': return 0;
            default:  printf("  ? unknown choice\n");
        }
    }
    return 0;
}
