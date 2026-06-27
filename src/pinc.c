/**
 * Pinc - a simple probabilistic AlphaFold interaction score
 *
 * Computes Pinc score from an AlphaFold PAE matrix (JSON)
 * and the corresponding structure file (PDB/CIF).
 *
 * git: https://git.mpi-cbg.de/tothpetroczylab/Pinc
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <io.h>
#define F_OK 0
#define W_OK 2
#define access _access
#define strcasecmp _stricmp
#else
#include <unistd.h>
#endif

/* -------------------------------------------------------------------------- */
/* Constants                                                                  */
/* -------------------------------------------------------------------------- */
#define CONTACT_RADIUS 12.0
#define MAX_LINE 2048
#define MAX_CHAIN_LEN 16
#define MAX_SEQ_LEN 16
#define MAX_TOKENS 50000
#define MAX_CIF_COLUMNS 128

#ifndef M_PI
#define M_PI 3.141593
#endif

// packed 1D index for a symmetric NxN upper-triangular matrix
// caller must guarantee i <= j; debug builds assert this, release (NDEBUG) builds trust it
#ifdef NDEBUG
#define PACKED_IDX(i, j, n) \
    ((size_t)(i) * (size_t)(n) - ((size_t)(i) * ((size_t)(i) + 1)) / 2 + (size_t)(j))
#else
#define PACKED_IDX(i, j, n) \
    (assert((size_t)(i) <= (size_t)(j)), \
     (size_t)(i) * (size_t)(n) - ((size_t)(i) * ((size_t)(i) + 1)) / 2 + (size_t)(j))
#endif

// safe for i/j in either order; swaps them before delegating to PACKED_IDX
#define SYM_IDX(i, j, n) ((i) <= (j) ? PACKED_IDX(i, j, n) : PACKED_IDX(j, i, n))

/* -------------------------------------------------------------------------- */
/* Data structures                                                            */
/* -------------------------------------------------------------------------- */

// string and identifier metadata for one residue; kept separate from coordinates
// (Structure-of-Arrays layout) so the hot O(n^2) distance loop reads only dense
// coordinate data and stays cache-friendly and auto-vectorizable
typedef struct {
    char chain[MAX_CHAIN_LEN];
    char seq_id[MAX_SEQ_LEN];
    char token[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
    int chain_id;
} ResidueMeta;

typedef char ChainName[MAX_CHAIN_LEN];

// per-atom coordinates tagged with their owning residue index; collected while
// parsing but never consumed elsewhere beyond being freed
typedef struct { int res_idx; double x, y, z; } RawAtom;

// parsed residues in Structure-of-Arrays form: meta and the x/y/z coordinates live
// in separate parallel flat arrays (see ResidueMeta above for the rationale)
typedef struct {
    ResidueMeta *meta;
    double *x;
    double *y;
    double *z;
    int res_count;
    ChainName *uchains;
    int num_uchains;
    RawAtom *atoms;
    int atom_count;
} ParsedStruct;

// aggregated score for one unique pair of chains; pinc1/pinc2 are the mean contact
// probability viewed from each chain's own (asymmetric) PAE column, pinc_score is
// their average and the symmetric Pinc score reported by default.
// order records the pair's generation index so the descending sort breaks ties stably
typedef struct {
    char chain1[MAX_CHAIN_LEN];
    char chain2[MAX_CHAIN_LEN];
    double pinc1;
    double pinc2;
    double pinc_score;
    int order;
} ChainPair;

// one row of the --pairlist output: a single (already symmetrized, see cont_sym)
// residue-residue contact together with its centroid-to-centroid distance.
// lo/hi are the residue indices (lo < hi); they let the descending sort break ties in
// the same column-major upper-triangle order the R reference emits for equal probabilities
typedef struct {
    char token1[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
    char token2[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
    double contact_p;
    double distance;
    int lo;
    int hi;
} TokenPair;

// per-chain running counter used to invent sequential seq_ids for mmCIF residues
// whose label_seq_id is "." (typically heteroatoms/waters with no canonical numbering)
typedef struct {
    char name[MAX_CHAIN_LEN];
    int count;
} ChainDotCounter;

// CSR-style grouping of residue indices by chain: indices belonging to chain c live
// in flat[offsets[c] .. offsets[c]+counts[c]-1]
typedef struct {
    int *flat;
    int *offsets;
    int *counts;
    int num_chains;
} ChainResIndex;

// open-addressing hash table (linear probing) mapping a residue token string to its
// index in the metadata array; capacity is always a power of two so '& (capacity-1)'
// replaces the modulo operation, and empty slots are marked with -1
typedef struct {
    int *res_indices;
    int capacity;
    int size;
} ResHash;

/* -------------------------------------------------------------------------- */
/* Utils                                                                      */
/* -------------------------------------------------------------------------- */

// fail with an error message and exit
void fail(const char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(EXIT_FAILURE);
}

// trim whitespace from the start and end of a string in-place
char *trim_in_place(char *str) {
    if (!str || !*str) return str;

    char *start = str;
    while (isspace((unsigned char)*start)) start++;

    if (start != str) memmove(str, start, strlen(start) + 1);

    if (*str == 0) return str;

    char *end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';

    return str;
}

// thread-safe tokenization that collapses runs of delimiters, so consecutive
// whitespace never yields empty tokens; fine for mmCIF's whitespace-separated
// columns, but would swallow a genuinely empty field
char *safe_strtok_r(char *str, const char *delim, char **saveptr) {
    if (!str) str = *saveptr;
    if (!str) return NULL;

    str += strspn(str, delim);
    if (!*str) {
        *saveptr = NULL;
        return NULL;
    }

    char *token = str;
    str += strcspn(str, delim);
    if (*str) {
        *str++ = '\0';
    }
    *saveptr = str;

    return token;
}

// robust double parsing wrapper around strtod
// returns false (rather than 0.0) when no characters were consumed, e.g. blank fields
bool safe_parse_double(const char *str, double *val) {
    if (!str || !*str) return false;
    char *end;
    *val = strtod(str, &end);
    return (end != str);
}

bool has_extension(const char *file, const char *ext1, const char *ext2) {
    const char *dot = strrchr(file, '.');
    if (!dot) return false;
    dot++;
    if (strcasecmp(dot, ext1) == 0) return true;
    if (ext2 && strcasecmp(dot, ext2) == 0) return true;
    return false;
}

// atomic mass by element symbol; anything not S/P/O/N (hydrogen included)
// falls back to carbon's mass
double get_mass(const char *elem) {
    if (strcasecmp(elem, "S") == 0) return 32.0650;
    if (strcasecmp(elem, "P") == 0) return 30.9738;
    if (strcasecmp(elem, "O") == 0) return 15.9994;
    if (strcasecmp(elem, "N") == 0) return 14.0067;
    return 12.0107;
}

double vol_sphere(double r) {
    return (4.0 / 3.0) * M_PI * (r * r * r);
}

// sphere-sphere intersection volume
double vol_intersect(double Ru, double D) {
    double rc = CONTACT_RADIUS;
    if (!isfinite(Ru) || !isfinite(D) || Ru <= 0.0) return 0.0;
    // spheres are too far apart to overlap at all
    if (D >= (rc + Ru)) return 0.0;
    // one sphere fully contains the other; intersection is just the smaller sphere's volume
    if (D <= fabs(rc - Ru)) return vol_sphere(fmin(rc, Ru));

    // standard closed-form volume of the lens formed by two partially overlapping spheres
    double t1 = rc + Ru - D;
    double t2 = Ru - rc;
    double num = M_PI * (t1 * t1) * (D * D + 2.0 * D * (Ru + rc) - 3.0 * (t2 * t2));
    return fmax(0.0, num / (12.0 * D));
}

// contact probability from a single PAE value and centroid distance, clamped to [0,1];
// factored out to avoid repeating the same guard/clamp pattern at every call site
double contact_prob(double pae, double d) {
    if (pae <= 0.0) return 0.0;
    double v_unc = vol_sphere(pae);
    if (v_unc <= 0.0) return 0.0;
    double p = vol_intersect(pae, d) / v_unc;
    if (!isfinite(p) || p < 0.0) return 0.0;
    return (p > 1.0) ? 1.0 : p;
}

// FNV-1a hash for strings
uint32_t hash_str(const char *str) {
    uint32_t hash = 2166136261u;
    while (*str) {
        hash ^= (unsigned char)(*str++);
        hash *= 16777619u;
    }
    return hash;
}

// allocate a hash table; capacity is rounded up to the next power of two so lookups
// can mask instead of modulo
ResHash *rhash_create(int target_capacity) {
    ResHash *h = malloc(sizeof(ResHash));
    if (!h) fail("Hash table allocation failed.");

    int p2_cap = 1;
    while (p2_cap < target_capacity) p2_cap <<= 1;

    h->capacity = p2_cap;
    h->size = 0;
    h->res_indices = malloc(h->capacity * sizeof(int));
    if (!h->res_indices) fail("Hash indices allocation failed.");

    for (int i = 0; i < h->capacity; i++) h->res_indices[i] = -1;
    return h;
}

void rhash_resize(ResHash *h, ResidueMeta *meta) {
    int old_cap = h->capacity;
    int *old_idx = h->res_indices;

    h->capacity *= 2;
    h->res_indices = malloc(h->capacity * sizeof(int));
    if (!h->res_indices) fail("Hash resize allocation failed.");

    for (int i = 0; i < h->capacity; i++) h->res_indices[i] = -1;

    for (int i = 0; i < old_cap; i++) {
        if (old_idx[i] != -1) {
            int ri = old_idx[i];
            uint32_t idx = hash_str(meta[ri].token) & (h->capacity - 1);
            while (h->res_indices[idx] != -1) {
                idx = (idx + 1) & (h->capacity - 1);
            }
            h->res_indices[idx] = ri;
        }
    }
    free(old_idx);
}

// residue index for a token, or -1 if not found; linear probing walks forward
// until the token matches or an empty (-1) slot is hit
int rhash_get(ResHash *h, ResidueMeta *meta, const char *token) {
    if (h->size == 0) return -1;
    uint32_t idx = hash_str(token) & (h->capacity - 1);
    while (h->res_indices[idx] != -1) {
        int ri = h->res_indices[idx];
        if (strcmp(meta[ri].token, token) == 0) return ri;
        idx = (idx + 1) & (h->capacity - 1);
    }
    return -1;
}

// insert a residue index, resizing first if the load factor would exceed 50% so the
// new entry is never written into a table that's about to be reallocated
void rhash_insert(ResHash *h, ResidueMeta *meta, int ri) {
    if (h->size >= h->capacity / 2) {
        rhash_resize(h, meta);
    }
    uint32_t idx = hash_str(meta[ri].token) & (h->capacity - 1);
    while (h->res_indices[idx] != -1) {
        idx = (idx + 1) & (h->capacity - 1);
    }
    h->res_indices[idx] = ri;
    h->size++;
}

void rhash_free(ResHash *h) {
    free(h->res_indices);
    free(h);
}

// build chain index: two-pass (count then fill) with a single flat array
ChainResIndex build_chain_index(const ResidueMeta *meta, int n_res, int num_uchains) {
    ChainResIndex ci;
    ci.num_chains = num_uchains;
    ci.counts = calloc(num_uchains, sizeof(int));
    ci.offsets = malloc(num_uchains * sizeof(int));
    if (!ci.counts || !ci.offsets) fail("Chain index allocation failed.");

    // first pass: count residues per chain
    for (int i = 0; i < n_res; i++) {
        ci.counts[meta[i].chain_id]++;
    }

    // compute offsets from counts
    ci.offsets[0] = 0;
    for (int c = 1; c < num_uchains; c++) {
        ci.offsets[c] = ci.offsets[c - 1] + ci.counts[c - 1];
    }

    // second pass: fill flat array using a temporary write-head per chain
    ci.flat = malloc(n_res * sizeof(int));
    int *write_pos = malloc(num_uchains * sizeof(int));
    if (!ci.flat || !write_pos) fail("Chain index flat array allocation failed.");
    memcpy(write_pos, ci.offsets, num_uchains * sizeof(int));

    for (int i = 0; i < n_res; i++) {
        int cid = meta[i].chain_id;
        ci.flat[write_pos[cid]++] = i;
    }
    free(write_pos);

    return ci;
}

void free_chain_index(ChainResIndex *ci) {
    free(ci->flat);
    free(ci->offsets);
    free(ci->counts);
}

// sorting helpers (for consistency with R version)
int cmp_chainpair(const void *a, const void *b) {
    const ChainPair *pa = (const ChainPair *)a;
    const ChainPair *pb = (const ChainPair *)b;
    char sa[32], sb[32];
    snprintf(sa, sizeof(sa), "%.4f", pa->pinc_score);
    snprintf(sb, sizeof(sb), "%.4f", pb->pinc_score);
    int c = strcmp(sb, sa);
    if (c != 0) return c;
    return (pa->order > pb->order) - (pa->order < pb->order);
}

int cmp_tokenpair(const void *a, const void *b) {
    const TokenPair *pa = (const TokenPair *)a;
    const TokenPair *pb = (const TokenPair *)b;
    char sa[32], sb[32];
    snprintf(sa, sizeof(sa), "%.4f", pa->contact_p);
    snprintf(sb, sizeof(sb), "%.4f", pb->contact_p);
    int c = strcmp(sb, sa);
    if (c != 0) return c;
    if (pa->hi != pb->hi) return (pa->hi > pb->hi) - (pa->hi < pb->hi);
    return (pa->lo > pb->lo) - (pa->lo < pb->lo);
}

/* -------------------------------------------------------------------------- */
/* File parsers                                                               */
/* -------------------------------------------------------------------------- */

// minimal hand-rolled scanner for AlphaFold's PAE output; it does not parse JSON
// structurally, it just locates the PAE array and reads every numeric token inside,
// using bracket-depth tracking only to find where that array ends
double *parse_pae_json(const char *json_file, int *n_res) {
    FILE *f = fopen(json_file, "rb");
    if (!f) fail("could not open JSON file.");

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *buffer = malloc(fsize + 1);
    if (!buffer) fail("memory allocation failed for JSON file.");
    if (fread(buffer, 1, fsize, f) != (size_t)fsize) fail("failed to read JSON.");
    buffer[fsize] = '\0';
    fclose(f);

    // AlphaFold2/Multimer and AlphaFold3 use different key names for the same matrix
    char *ptr = strstr(buffer, "\"predicted_aligned_error\"");
    if (!ptr) ptr = strstr(buffer, "\"pae\"");
    if (!ptr) fail("no PAE matrix found in JSON.");

    while (*ptr && *ptr != '[') ptr++;
    if (!*ptr) fail("malformed JSON: expected '[' after PAE key.");

    // depth starts at 1 for the '[' just consumed and only returns to 0 at that same
    // array's matching ']'; inner row brackets are tracked but otherwise ignored, so
    // every number inside (regardless of nesting) lands in one flat buffer
    int depth = 1;
    ptr++;

    // generous initial allocation, doubled via realloc below if it is exceeded
    size_t max_alloc = 1024 * 1024;
    double *pae = malloc(max_alloc * sizeof(double));
    if (!pae) {
        free(buffer);
        fail("memory allocation failed for PAE array.");
    }

    size_t count = 0;
    while (*ptr && depth > 0) {
        if (*ptr == '[') {
            depth++;
            ptr++;
        } else if (*ptr == ']') {
            depth--;
            ptr++;
        } else if (isdigit((unsigned char)*ptr) || *ptr == '-' || *ptr == '.') {
            char *end;
            double val = strtod(ptr, &end);
            pae[count++] = val;
            if (count >= max_alloc) {
                max_alloc *= 2;
                double *tmp = realloc(pae, max_alloc * sizeof(double));
                if (!tmp) {
                    free(pae);
                    free(buffer);
                    fail("OOM for PAE array.");
                }
                pae = tmp;
            }
            ptr = end;
        } else {
            ptr++;
        }
    }
    free(buffer);

    // the PAE matrix is always n_res x n_res, so n_res can be recovered from the total
    // element count; a non-square count means parsing went wrong or the file is malformed
    *n_res = (int)sqrt((double)count);
    if ((size_t)(*n_res) * (size_t)(*n_res) != count) fail("parsed PAE matrix is not perfectly square.");
    return pae;
}

// parses a PDB (fixed-column) or mmCIF (whitespace-delimited, columns taken from the
// atom_site loop header) file into one row per residue, each residue's position set to
// the mass-weighted centroid of its heavy (non-H) atoms.
// the whole file is read into one bulk buffer to avoid per-line syscall overhead;
// coordinates are stored in separate SoA flat arrays for cache efficiency
ParsedStruct parse_structure(const char *str_file) {
    FILE *f = fopen(str_file, "rb");
    if (!f) fail("could not open structure file.");

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *buffer = malloc(fsize + 1);
    if (!buffer) fail("memory allocation failed for structure file buffer.");
    if (fread(buffer, 1, fsize, f) != (size_t)fsize) fail("failed to read structure file.");
    buffer[fsize] = '\0';
    fclose(f);

    int cap = 1000;
    // SoA coordinate arrays, parallel to meta (see ResidueMeta for the rationale)
    ResidueMeta *meta   = malloc(cap * sizeof(ResidueMeta));
    double *res_x       = malloc(cap * sizeof(double));
    double *res_y       = malloc(cap * sizeof(double));
    double *res_z       = malloc(cap * sizeof(double));
    double *mass_sum    = malloc(cap * sizeof(double));
    if (!meta || !res_x || !res_y || !res_z || !mass_sum) fail("OOM for SoA structure arrays.");

    int count = 0;

    int uc_cap = 16;
    ChainName *uchains = malloc(uc_cap * sizeof(ChainName));
    if (!uchains) fail("memory allocation failed for unique chains.");
    int num_uchains = 0;

    ResHash *rhash = rhash_create(4096);

    int at_cap = 8192, at_cnt = 0;
    RawAtom *atoms = malloc(at_cap * sizeof(RawAtom));
    if (!atoms) fail("OOM atoms.");

    bool is_cif = false;

    // c_* hold the column ordinal of each required mmCIF field, discovered from the
    // atom_site loop header below; -1 means "not found yet"
    int c_chain = -1, c_seq = -1, c_elem = -1, c_x = -1, c_y = -1, c_z = -1;
    bool in_atom_site = false;
    bool cif_columns_validated = false;
    int cif_col_count = 0;

    // tracks, per chain, how many '.' seq_ids have already been substituted (see below)
    int dc_cap = 16;
    ChainDotCounter *dot_counters = malloc(dc_cap * sizeof(ChainDotCounter));
    if (!dot_counters) fail("memory allocation failed for dot counters.");
    int num_chains_tracked = 0;

    char *ptr = buffer;
    while (*ptr) {
        // locate the end of the current line within the bulk buffer and null-terminate
        // it in-place so downstream string functions see a normal C string; DOS-style
        // CRLF is handled by also nulling the '\r' that precedes the '\n'
        char *line_end = strchr(ptr, '\n');
        if (line_end) {
            *line_end = '\0';
            if (line_end > ptr && *(line_end - 1) == '\r') {
                *(line_end - 1) = '\0';
            }
        }

        char *line = ptr;

        // presence of any _atom_site.* tag identifies the file as mmCIF and marks the
        // start of its column-definition header
        if (strstr(line, "_atom_site.group_PDB")) {
            is_cif = true;
            in_atom_site = true;
        }

        // each header line declares the next column in the data rows that follow, in
        // order; record the ordinal of every field this program needs and let
        // cif_col_count keep advancing for fields it does not care about.
        // the line is copied to a stack buffer before trimming because line points into
        // the bulk buffer and trim_in_place would corrupt subsequent lines if applied
        // directly to it
        if (is_cif && in_atom_site && line[0] == '_') {
            char key_buf[MAX_LINE];
            strncpy(key_buf, line, sizeof(key_buf) - 1);
            key_buf[sizeof(key_buf) - 1] = '\0';
            char *key = trim_in_place(key_buf);

            if (strcmp(key, "_atom_site.label_asym_id") == 0) c_chain = cif_col_count;
            else if (strcmp(key, "_atom_site.label_seq_id") == 0) c_seq = cif_col_count;
            else if (strcmp(key, "_atom_site.type_symbol") == 0) c_elem = cif_col_count;
            else if (strcmp(key, "_atom_site.Cartn_x") == 0) c_x = cif_col_count;
            else if (strcmp(key, "_atom_site.Cartn_y") == 0) c_y = cif_col_count;
            else if (strcmp(key, "_atom_site.Cartn_z") == 0) c_z = cif_col_count;
            cif_col_count++;
        }
        else if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {

            // validate CIF columns once on the first ATOM/HETATM line;
            // fails early if any required _atom_site field is missing
            if (is_cif && !cif_columns_validated) {
                if (c_chain < 0 || c_seq < 0 || c_elem < 0 ||
                    c_x < 0 || c_y < 0 || c_z < 0) {
                    free(buffer);
                    free(meta); free(res_x); free(res_y); free(res_z); free(mass_sum);
                    free(dot_counters);
                    rhash_free(rhash);
                    free(uchains);
                    free(atoms);
                    fail("CIF file missing required _atom_site columns "
                         "(need label_asym_id, label_seq_id, type_symbol, "
                         "Cartn_x, Cartn_y, Cartn_z).");
                }
                cif_columns_validated = true;
            }

            char chain[MAX_CHAIN_LEN] = {0};
            char seq_id[MAX_SEQ_LEN] = {0};
            char elem[4] = "C";
            double x = 0, y = 0, z = 0;

            if (is_cif) {
                // safe_strtok_r mutates its input in place, so copy the line before
                // tokenizing to avoid corrupting the bulk buffer for subsequent lines
                char line_copy[MAX_LINE];
                strncpy(line_copy, line, sizeof(line_copy) - 1);
                line_copy[sizeof(line_copy) - 1] = '\0';

                char *tokens[MAX_CIF_COLUMNS];
                int t_count = 0;
                char *saveptr = NULL;
                char *tok = safe_strtok_r(line_copy, " \t\n\r", &saveptr);

                while (tok && t_count < MAX_CIF_COLUMNS) {
                    tokens[t_count++] = tok;
                    tok = safe_strtok_r(NULL, " \t\n\r", &saveptr);
                }

                if (c_chain < t_count) {
                    if ((size_t)snprintf(chain, sizeof(chain), "%s", tokens[c_chain]) >= sizeof(chain)) goto next_line;
                }
                if (c_seq < t_count) {
                    if ((size_t)snprintf(seq_id, sizeof(seq_id), "%s", tokens[c_seq]) >= sizeof(seq_id)) goto next_line;
                }
                if (c_elem < t_count) strncpy(elem, tokens[c_elem], 3);

                if (c_x < t_count && !safe_parse_double(tokens[c_x], &x)) goto next_line;
                if (c_y < t_count && !safe_parse_double(tokens[c_y], &y)) goto next_line;
                if (c_z < t_count && !safe_parse_double(tokens[c_z], &z)) goto next_line;

            } else {
                // fixed-width PDB columns (0-indexed here): chain at 21, seq number at
                // 22-25, x/y/z at 30-37/38-45/46-53, element symbol at 76-77; a line
                // shorter than 54 chars is skipped rather than read out of bounds
                size_t len = strlen(line);
                if (len < 54) goto next_line;

                chain[0] = line[21];
                chain[1] = '\0';

                memcpy(seq_id, line + 22, 4);
                seq_id[4] = '\0';
                trim_in_place(seq_id);

                char coord_buf[16] = {0};

                memcpy(coord_buf, line + 30, 8);
                coord_buf[8] = '\0';
                if (!safe_parse_double(coord_buf, &x)) goto next_line;

                memcpy(coord_buf, line + 38, 8);
                coord_buf[8] = '\0';
                if (!safe_parse_double(coord_buf, &y)) goto next_line;

                memcpy(coord_buf, line + 46, 8);
                coord_buf[8] = '\0';
                if (!safe_parse_double(coord_buf, &z)) goto next_line;

                if (len >= 78) {
                    memcpy(elem, line + 76, 2);
                    elem[2] = '\0';
                    trim_in_place(elem);
                }
            }

            trim_in_place(chain);
            trim_in_place(elem);
            // hydrogens are excluded here because get_mass() has no entry for H and
            // would otherwise silently weight it as carbon in the centroid below
            if (strcasecmp(elem, "H") == 0) goto next_line;

            // mmCIF marks residues with no canonical numbering (often heteroatoms or
            // waters) with label_seq_id "."; substitute a synthetic, per-chain
            // sequential number so each one still gets a distinct token below
            if (strcmp(seq_id, ".") == 0) {
                int dc = 1;
                for (int i = 0; i < num_chains_tracked; i++) {
                    if (strcmp(dot_counters[i].name, chain) == 0) {
                        dc = ++dot_counters[i].count;
                        break;
                    }
                }
                if (dc == 1) {
                    if (num_chains_tracked >= dc_cap) {
                        dc_cap *= 2;
                        ChainDotCounter *tmp = realloc(dot_counters, dc_cap * sizeof(ChainDotCounter));
                        if (!tmp) fail("OOM dot counters.");
                        dot_counters = tmp;
                    }
                    strcpy(dot_counters[num_chains_tracked].name, chain);
                    dot_counters[num_chains_tracked].count = 1;
                    num_chains_tracked++;
                }
                if ((size_t)snprintf(seq_id, sizeof(seq_id), "%d", dc) >= sizeof(seq_id)) goto next_line;
            }

            double mass = get_mass(elem);
            // token uniquely identifies a residue across the whole structure; it is
            // the key used by both the fast-path check and the hash table below
            char token[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
            if ((size_t)snprintf(token, sizeof(token), "%s_%s", chain, seq_id) >= sizeof(token)) goto next_line;

            // register chain globally; linear scan is fine since the number of
            // distinct chains stays small
            int cid = -1;
            for (int i = 0; i < num_uchains; i++) {
                if (strcmp(uchains[i], chain) == 0) {
                    cid = i;
                    break;
                }
            }
            if (cid == -1) {
                if (num_uchains >= uc_cap) {
                    uc_cap *= 2;
                    uchains = realloc(uchains, uc_cap * sizeof(ChainName));
                    if (!uchains) fail("OOM chains.");
                }
                strcpy(uchains[num_uchains], chain);
                cid = num_uchains++;
            }

            // check the last entry before falling through to the hash table;
            // atoms of one residue are almost always contiguous in the file, so this
            // fast path resolves the vast majority of atoms without a hash lookup
            int ri = -1;
            if (count > 0 && strcmp(meta[count - 1].token, token) == 0) {
                ri = count - 1;
            } else {
                ri = rhash_get(rhash, meta, token);
            }

            if (ri == -1) {
                // guard before any allocation so the error fires before leaving a
                // partially-constructed state that would complicate cleanup
                if (count >= MAX_TOKENS) {
                    free(buffer);
                    free(meta); free(res_x); free(res_y); free(res_z); free(mass_sum);
                    free(dot_counters);
                    rhash_free(rhash);
                    free(uchains);
                    free(atoms);
                    fprintf(stderr, "Structure exceeds %d residues.\n", MAX_TOKENS);
                    fail("Structure too large.");
                }

                // resize all parallel SoA arrays together so they always stay the
                // same length; a failure on any one is caught before the others change
                if (count >= cap) {
                    cap *= 2;
                    ResidueMeta *tm  = realloc(meta,     cap * sizeof(ResidueMeta));
                    double      *tx  = realloc(res_x,    cap * sizeof(double));
                    double      *ty  = realloc(res_y,    cap * sizeof(double));
                    double      *tz  = realloc(res_z,    cap * sizeof(double));
                    double      *tms = realloc(mass_sum, cap * sizeof(double));
                    if (!tm || !tx || !ty || !tz || !tms) fail("OOM structure arrays.");
                    meta = tm; res_x = tx; res_y = ty; res_z = tz; mass_sum = tms;
                }

                snprintf(meta[count].chain,  sizeof(meta[count].chain),  "%s", chain);
                snprintf(meta[count].seq_id, sizeof(meta[count].seq_id), "%s", seq_id);
                snprintf(meta[count].token,  sizeof(meta[count].token),  "%s", token);
                meta[count].chain_id = cid;

                res_x[count]    = 0.0;
                res_y[count]    = 0.0;
                res_z[count]    = 0.0;
                mass_sum[count] = 0.0;
                ri = count;

                rhash_insert(rhash, meta, ri);
                count++;
            }

            // accumulate mass-weighted sums; divided by mass_sum once parsing is
            // complete (see below) to get each residue's final centroid coordinate
            res_x[ri] += x * mass;
            res_y[ri] += y * mass;
            res_z[ri] += z * mass;
            mass_sum[ri] += mass;

            if (at_cnt >= at_cap) {
                at_cap *= 2;
                RawAtom *tmp = realloc(atoms, at_cap * sizeof(RawAtom));
                if (!tmp) fail("OOM atoms array.");
                atoms = tmp;
            }
            atoms[at_cnt++] = (RawAtom){ ri, x, y, z };
        }

next_line:
        // advance the read pointer to the character after the '\n' that was nulled
        // above; if no '\n' was found we have consumed the last (unterminated) line
        if (line_end) ptr = line_end + 1;
        else break;
    }

    free(buffer);
    free(dot_counters);
    rhash_free(rhash);

    // crude sanity check to catch a wrong/non-structure file being passed in, rather
    // than silently proceeding with a near-empty structure
    if (count < 10) fail("non-structure file or insufficient ATOM records.");

    // divide accumulated mass-weighted sums by total mass to get each residue's
    // centroid; a residue with mass_sum == 0 (e.g. all-hydrogen) is left at the
    // origin rather than flagged — worth inspecting if contacts look wrong
    for (int i = 0; i < count; i++) {
        if (mass_sum[i] > 0.0) {
            res_x[i] /= mass_sum[i];
            res_y[i] /= mass_sum[i];
            res_z[i] /= mass_sum[i];
        }
    }
    // internal accumulator no longer needed after centroid coords are finalized
    free(mass_sum);

    ParsedStruct ps = { meta, res_x, res_y, res_z, count, uchains, num_uchains, atoms, at_cnt };
    return ps;
}

/* -------------------------------------------------------------------------- */
/* Main program                                                               */
/* -------------------------------------------------------------------------- */

// CLI entry point: validates inputs, computes per-residue distances and PAE-derived
// contact probabilities, aggregates them into a Pinc score per chain pair, and writes
// the default CSV plus any of the optional --all / --pairlist outputs
int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Description:\n"
               "   A program to compute the Pinc score from an AlphaFold PAE\n"
               "   matrix (JSON) and the corresponding structure file (PDB/CIF).\n\n"
               "Usage:\n"
               "   Pinc <json_file> <structure_file> [options]\n\n"
               "Options:\n"
               "   None        Pinc score for all unique chain pairs (default)\n"
               "   --all       Full contact probability matrix\n"
               "   --pairlist  Residue pair list of non-zero contact probabilities\n"
               "   --results   Optional path to output folder\n\n"
               "Output format:\n"
               "   - default:  <structure_file_name>_Pinc.csv\n"
               "   - all:      <structure_file_name>_contact_probability.json\n"
               "   - pairlist: <structure_file_name>_pairlist.csv\n\n");
        return 1;
    }

    const char *jsn_file = argv[1];
    const char *str_file = argv[2];
    bool opt_all = false;
    bool opt_pairlist = false;
    char results_dir[512] = ".";

    // access() with F_OK only checks existence, not readability; a permissions
    // problem will instead surface later as a fopen() failure inside the parsers
    if (!has_extension(jsn_file, "json", NULL)) fail("JSON file must have .json extension.");
    if (!has_extension(str_file, "pdb", "cif")) fail("structure file must have .pdb/.cif extension.");
    if (access(jsn_file, F_OK) != 0 || access(str_file, F_OK) != 0) fail("input file(s) not found.");

    // unrecognized arguments are silently ignored rather than rejected
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--all") == 0) opt_all = true;
        else if (strcmp(argv[i], "--pairlist") == 0) opt_pairlist = true;
        else if (strcmp(argv[i], "--results") == 0) {
            if (i + 1 < argc) {
                if ((size_t)snprintf(results_dir, sizeof(results_dir), "%s", argv[++i]) >= sizeof(results_dir))
                    fail("Results directory path too long.");
            } else fail("option --results requires a folder path.");
        }
    }

    if (access(results_dir, W_OK) != 0) fail("results directory does not exist or is not writable.");

    // strip trailing slash(es) so the path joins below never produce a double slash
    size_t rlen = strlen(results_dir);
    while (rlen > 0 && (results_dir[rlen - 1] == '/' || results_dir[rlen - 1] == '\\')) {
        results_dir[rlen - 1] = '\0';
        rlen--;
    }

    // derive a bare stem from the structure filename (drop directory and extension)
    // so all output files share one consistent base name under results_dir
    const char *base_name = str_file;
    for (const char *p = str_file; *p; p++) {
        if (*p == '/' || *p == '\\') base_name = p + 1;
    }
    char base_no_ext[256];
    if ((size_t)snprintf(base_no_ext, sizeof(base_no_ext), "%s", base_name) >= sizeof(base_no_ext))
        fail("Base filename too long.");
    char *dot = strrchr(base_no_ext, '.');
    if (dot) *dot = '\0';

    char out_path[1024];
    if ((size_t)snprintf(out_path, sizeof(out_path), "%s/%s", results_dir, base_no_ext) >= sizeof(out_path))
        fail("Output path too long.");

    int n_pae = 0;
    double *pae_mat = parse_pae_json(jsn_file, &n_pae);

    ParsedStruct ps = parse_structure(str_file);
    int n_res       = ps.res_count;
    ResidueMeta *meta = ps.meta;
    double *res_x   = ps.x;
    double *res_y   = ps.y;
    double *res_z   = ps.z;
    int num_uchains = ps.num_uchains;
    ChainName *uchains = ps.uchains;
    RawAtom *atoms  = ps.atoms;

    // the single most important correctness check in the program: everything below
    // assumes pae_mat row/column i and meta[i] refer to the same i-th residue, in the
    // same order, so an undetected size mismatch here would silently corrupt every score
    if (n_pae != n_res) {
        fprintf(stderr, "Dim mismatch. PAE: %dx%d. Struct: %dx%d.\n", n_pae, n_pae, n_res, n_res);
        fail("Check that structure and PAE correspond to the same model.");
    }

    // allocate the single packed symmetric upper-triangular contact matrix;
    // dist_mat is intentionally omitted: storing it would cost O(n^2/2) doubles and
    // profiling showed the resulting cache pressure outweighed the savings from
    // avoiding distance recomputation in the downstream chain-pair and pairlist loops,
    // which only execute for the small in-contact fraction of pairs anyway
    size_t sym_size = ((size_t)n_res * (size_t)(n_res + 1)) / 2;
    double *cont_sym = malloc(sym_size * sizeof(double));
    if (!cont_sym) fail("Symmetric matrix allocation failed.");

    // compute symmetric contact probabilities for all residue pairs;
    // SoA coordinate arrays (res_x/y/z) keep cache lines dense for this hot loop
    for (size_t i = 0; i < (size_t)n_res; i++) {
        cont_sym[PACKED_IDX(i, i, n_res)] = 1.0;

        for (size_t j = i + 1; j < (size_t)n_res; j++) {
            double dx = res_x[i] - res_x[j];
            double dy = res_y[i] - res_y[j];
            double dz = res_z[i] - res_z[j];
            double d  = sqrt(dx * dx + dy * dy + dz * dz);

            size_t p_idx = PACKED_IDX(i, j, n_res);

            // the PAE matrix is NOT symmetric (pae[i][j] is the error in j's position
            // when aligned on i, and vice versa), so both directions are evaluated
            // separately and only averaged into cont_sym at the end
            double p_ij = contact_prob(pae_mat[i * n_res + j], d);
            double p_ji = contact_prob(pae_mat[j * n_res + i], d);
            cont_sym[p_idx] = (p_ij + p_ji) / 2.0;
        }
    }

    // chain index: single flat array with offset/count per chain
    ChainResIndex ci = build_chain_index(meta, n_res, num_uchains);

    size_t total_pairs = ((size_t)num_uchains * ((size_t)num_uchains - 1)) / 2;
    // every field defaults to zero even for a pair with no in-range contacts
    ChainPair *cpairs = calloc(total_pairs, sizeof(ChainPair));
    if (!cpairs && total_pairs > 0) fail("OOM chain pairs.");

    size_t cp_idx = 0;

    for (int i = 0; i < num_uchains; i++) {
        for (int j = i + 1; j < num_uchains; j++) {
            snprintf(cpairs[cp_idx].chain1, sizeof(cpairs[cp_idx].chain1), "%s", uchains[i]);
            snprintf(cpairs[cp_idx].chain2, sizeof(cpairs[cp_idx].chain2), "%s", uchains[j]);

            double sum1 = 0.0, sum2 = 0.0;
            int cnt = 0;

            for (int u = 0; u < ci.counts[i]; u++) {
                int ri = ci.flat[ci.offsets[i] + u];
                for (int v = 0; v < ci.counts[j]; v++) {
                    int rj = ci.flat[ci.offsets[j] + v];

                    double dx = res_x[ri] - res_x[rj];
                    double dy = res_y[ri] - res_y[rj];
                    double dz = res_z[ri] - res_z[rj];
                    double d  = sqrt(dx * dx + dy * dy + dz * dz);

                    if (isfinite(d) && d < CONTACT_RADIUS) {
                        // ri belongs to chain i (aligned), rj belongs to chain j (scored)
                        // Pinc1 = chain i residues are scored → use pae[rj][ri]
                        // Pinc2 = chain j residues are scored → use pae[ri][rj]
                        sum1 += contact_prob(pae_mat[(size_t)rj * n_res + ri], d);
                        sum2 += contact_prob(pae_mat[(size_t)ri * n_res + rj], d);

                        cnt++;
                    }
                }
            }
            // pinc1/pinc2: mean contact probability over this chain pair's contacts,
            // viewed from each chain's own PAE column; pinc_score: their symmetric average
            cpairs[cp_idx].pinc1      = (cnt > 0) ? (sum1 / cnt) : 0.0;
            cpairs[cp_idx].pinc2      = (cnt > 0) ? (sum2 / cnt) : 0.0;
            cpairs[cp_idx].pinc_score = (cnt > 0) ? ((sum1 + sum2) / (2.0 * cnt)) : 0.0;
            cpairs[cp_idx].order      = (int)cp_idx;
            cp_idx++;
        }
    }

    // highest-scoring chain pairs first (see cmp_chainpair)
    qsort(cpairs, total_pairs, sizeof(ChainPair), cmp_chainpair);

    // default Pinc.csv output
    char pinc_csv[1024];
    if ((size_t)snprintf(pinc_csv, sizeof(pinc_csv), "%s_Pinc.csv", out_path) >= sizeof(pinc_csv)) fail("Path long.");

    FILE *fp = fopen(pinc_csv, "w");
    if (fp) {
        fprintf(fp, "chain1,chain2,Pinc1,Pinc2,Pinc\n");
        for (size_t i = 0; i < total_pairs; i++) {
            fprintf(fp, "%s,%s,%.4f,%.4f,%.4f\n",
                    cpairs[i].chain1, cpairs[i].chain2,
                    cpairs[i].pinc1, cpairs[i].pinc2, cpairs[i].pinc_score);
        }
        fclose(fp);
    } else fprintf(stderr, "warning: could not write %s\n", pinc_csv);

    // optional --all output
    if (opt_all) {
        // unlike cont_sym this is the full unpacked, asymmetric n_res x n_res matrix
        // (memory cost O(n^2) rather than O(n^2/2)), since --all exposes both directions;
        // distances are recomputed here rather than read from a stored dist_mat
        double *cont_asym = malloc((size_t)n_res * (size_t)n_res * sizeof(double));
        if (!cont_asym) fail("Full contact matrix allocation failed.");

        for (size_t i = 0; i < (size_t)n_res; i++) {
            cont_asym[i * n_res + i] = 1.0;
            for (size_t j = i + 1; j < (size_t)n_res; j++) {
                double dx = res_x[i] - res_x[j];
                double dy = res_y[i] - res_y[j];
                double dz = res_z[i] - res_z[j];
                double d    = sqrt(dx * dx + dy * dy + dz * dz);
                double p_ij = contact_prob(pae_mat[i * n_res + j], d);
                double p_ji = contact_prob(pae_mat[j * n_res + i], d);
                cont_asym[i * n_res + j] = p_ij;
                cont_asym[j * n_res + i] = p_ji;
            }
        }

        char jsn_out[1024];
        if ((size_t)snprintf(jsn_out, sizeof(jsn_out), "%s_contact_probability.json", out_path) >= sizeof(jsn_out))
            fail("JSON output path too long.");

        // written by hand with fprintf rather than a library, mirroring the
        // hand-rolled scanning in parse_pae_json(); not minified but valid JSON
        FILE *fj = fopen(jsn_out, "w");
        if (fj) {
            fprintf(fj, "{\n  \"token_chain_ids\":[\n");
            for (int i = 0; i < n_res; i++) {
                fprintf(fj, "    \"%s\"%s\n", meta[i].chain, (i < n_res - 1) ? "," : "");
            }
            fprintf(fj, "  ],\n  \"contact_probability\":[\n");
            for (int i = 0; i < n_res; i++) {
                fprintf(fj, "[\n");
                for (int j = 0; j < n_res; j++) {
                    fprintf(fj, "      %.4f%s\n", cont_asym[(size_t)i * n_res + j], (j < n_res - 1) ? "," : "");
                }
                fprintf(fj, "    ]%s\n", (i < n_res - 1) ? "," : "");
            }
            fprintf(fj, "  ]\n}\n");
            fclose(fj);
        }

        free(cont_asym);
    }

    // optional --pairlist output
    // reuses cont_sym computed in the main loop; distance is recomputed on demand
    // for the small fraction of pairs with non-zero contact probability
    if (opt_pairlist) {
        int pl_cap = 1000, pcnt = 0;
        TokenPair *tpairs = malloc(pl_cap * sizeof(TokenPair));
        if (!tpairs) fail("OOM pairlist.");

        for (int ca = 0; ca < num_uchains; ca++) {
            for (int cb = ca + 1; cb < num_uchains; cb++) {
                for (int u = 0; u < ci.counts[ca]; u++) {
                    int ri = ci.flat[ci.offsets[ca] + u];
                    for (int v = 0; v < ci.counts[cb]; v++) {
                        int rj = ci.flat[ci.offsets[cb] + v];

                        // enforce consistent ordering (lo < hi) for
                        // upper-triangle packed index and output semantics
                        int lo = (ri < rj) ? ri : rj;
                        int hi = (ri < rj) ? rj : ri;

                        // reuse already-averaged cont_sym; distance is recomputed below
                        // so the pair can be filtered the same way the R reference does:
                        // keep iff it is a non-zero contact AND within CONTACT_RADIUS
                        size_t p_idx = PACKED_IDX(lo, hi, n_res);
                        double cp = cont_sym[p_idx];

                        double dx = res_x[lo] - res_x[hi];
                        double dy = res_y[lo] - res_y[hi];
                        double dz = res_z[lo] - res_z[hi];
                        double dd = sqrt(dx * dx + dy * dy + dz * dz);

                        if (cp > 0.0 && dd < CONTACT_RADIUS) {
                            if (pcnt >= pl_cap) {
                                pl_cap *= 2;
                                TokenPair *tmp = realloc(tpairs, pl_cap * sizeof(TokenPair));
                                if (!tmp) fail("OOM pairlist array.");
                                tpairs = tmp;
                            }
                            snprintf(tpairs[pcnt].token1, sizeof(tpairs[pcnt].token1), "%s", meta[lo].token);
                            snprintf(tpairs[pcnt].token2, sizeof(tpairs[pcnt].token2), "%s", meta[hi].token);
                            tpairs[pcnt].contact_p = cp;
                            tpairs[pcnt].distance  = dd;
                            tpairs[pcnt].lo        = lo;
                            tpairs[pcnt].hi        = hi;
                            pcnt++;
                        }
                    }
                }
            }
        }

        // highest contact probability first (see cmp_tokenpair)
        qsort(tpairs, pcnt, sizeof(TokenPair), cmp_tokenpair);

        char pair_csv[1024];
        if ((size_t)snprintf(pair_csv, sizeof(pair_csv), "%s_pairlist.csv", out_path) >= sizeof(pair_csv)) fail("Path.");

        FILE *ft = fopen(pair_csv, "w");
        if (ft) {
            fprintf(ft, "token1,token2,contact_p,distance\n");
            for (int i = 0; i < pcnt; i++) {
                fprintf(ft, "%s,%s,%.4f,%.4f\n",
                        tpairs[i].token1, tpairs[i].token2,
                        tpairs[i].contact_p, tpairs[i].distance);
            }
            fclose(ft);
        }
        free(tpairs);
    }

    // cleanup
    free_chain_index(&ci);
    free(meta);
    free(res_x);
    free(res_y);
    free(res_z);
    free(uchains);
    free(atoms);
    free(pae_mat);
    free(cont_sym);
    free(cpairs);

    return 0;
}