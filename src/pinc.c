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
#ifdef NDEBUG
#define PACKED_IDX(i, j, n) \
    ((size_t)(i) * (size_t)(n) - ((size_t)(i) * ((size_t)(i) + 1)) / 2 + (size_t)(j))
#else
#define PACKED_IDX(i, j, n) \
    (assert((size_t)(i) <= (size_t)(j)), \
     (size_t)(i) * (size_t)(n) - ((size_t)(i) * ((size_t)(i) + 1)) / 2 + (size_t)(j))
#endif

#define SYM_IDX(i, j, n) ((i) <= (j) ? PACKED_IDX(i, j, n) : PACKED_IDX(j, i, n))

/* -------------------------------------------------------------------------- */
/* Data structures                                                            */
/* -------------------------------------------------------------------------- */
typedef struct {
    char chain[MAX_CHAIN_LEN];
    char seq_id[MAX_SEQ_LEN];
    char token[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
    int chain_id;
    double x, y, z, mass_sum; 
} Residue;

typedef char ChainName[MAX_CHAIN_LEN];

typedef struct {
    Residue *res_arr;
    int res_count;
    ChainName *uchains;
    int num_uchains;
} ParsedStruct;

typedef struct {
    char chain1[MAX_CHAIN_LEN];
    char chain2[MAX_CHAIN_LEN];
    double pinc_score;
} ChainPair;

typedef struct {
    char token1[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
    char token2[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
    double contact_p;
    double distance;
} TokenPair;

typedef struct {
    char name[MAX_CHAIN_LEN];
    int count;
} ChainDotCounter;

typedef struct {
    int *flat;
    int *offsets;
    int *counts;
    int num_chains;
} ChainResIndex;

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

// thread-safe tokenization
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
bool safe_parse_double(const char *str, double *val) {
    if (!str || !*str) return false;
    char *end;
    *val = strtod(str, &end);
    return (end != str);
}

// check extension
bool has_extension(const char *file, const char *ext1, const char *ext2) {
    const char *dot = strrchr(file, '.');
    if (!dot) return false;
    dot++;
    if (strcasecmp(dot, ext1) == 0) return true;
    if (ext2 && strcasecmp(dot, ext2) == 0) return true;
    return false;
}

// assign atomic masses based on element symbol
double get_mass(const char *elem) {
    if (strcasecmp(elem, "S") == 0) return 32.0650;
    if (strcasecmp(elem, "P") == 0) return 30.9738;
    if (strcasecmp(elem, "O") == 0) return 15.9994;
    if (strcasecmp(elem, "N") == 0) return 14.0067;
    return 12.0107;
}

// sphere volume helper
double vol_sphere(double r) {
    return (4.0 / 3.0) * M_PI * (r * r * r);
}

// sphere-sphere intersection volume
double vol_intersect(double Ru, double D) {
    double rc = CONTACT_RADIUS;
    if (!isfinite(Ru) || !isfinite(D) || Ru <= 0.0) return 0.0;
    if (D >= (rc + Ru)) return 0.0;
    if (D <= fabs(rc - Ru)) return vol_sphere(fmin(rc, Ru));

    double t1 = rc + Ru - D;
    double t2 = Ru - rc;
    double num = M_PI * (t1 * t1) * (D * D + 2.0 * D * (Ru + rc) - 3.0 * (t2 * t2));
    return fmax(0.0, num / (12.0 * D));
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

// initialize hash table
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

// resize and rehash
void rhash_resize(ResHash *h, Residue *res_arr) {
    int old_cap = h->capacity;
    int *old_idx = h->res_indices;

    h->capacity *= 2;
    h->res_indices = malloc(h->capacity * sizeof(int));
    if (!h->res_indices) fail("Hash resize allocation failed.");
    
    for (int i = 0; i < h->capacity; i++) h->res_indices[i] = -1;

    for (int i = 0; i < old_cap; i++) {
        if (old_idx[i] != -1) {
            int ri = old_idx[i];
            uint32_t idx = hash_str(res_arr[ri].token) & (h->capacity - 1);
            while (h->res_indices[idx] != -1) {
                idx = (idx + 1) & (h->capacity - 1);
            }
            h->res_indices[idx] = ri;
        }
    }
    free(old_idx);
}

// retrieve residue index from hash table (-1 if not found)
int rhash_get(ResHash *h, Residue *res_arr, const char *token) {
    if (h->size == 0) return -1;
    uint32_t idx = hash_str(token) & (h->capacity - 1);
    while (h->res_indices[idx] != -1) {
        int ri = h->res_indices[idx];
        if (strcmp(res_arr[ri].token, token) == 0) return ri;
        idx = (idx + 1) & (h->capacity - 1);
    }
    return -1;
}

// insert residue index into hash table
void rhash_insert(ResHash *h, Residue *res_arr, int ri) {
    if (h->size >= h->capacity / 2) {
        rhash_resize(h, res_arr);
    }
    uint32_t idx = hash_str(res_arr[ri].token) & (h->capacity - 1);
    while (h->res_indices[idx] != -1) {
        idx = (idx + 1) & (h->capacity - 1);
    }
    h->res_indices[idx] = ri;
    h->size++;
}

// free hash table
void rhash_free(ResHash *h) {
    free(h->res_indices);
    free(h);
}

// build chain index: two-pass (count then fill) with a single flat array
ChainResIndex build_chain_index(const Residue *res, int n_res, int num_uchains) {
    ChainResIndex ci;
    ci.num_chains = num_uchains;
    ci.counts = calloc(num_uchains, sizeof(int));
    ci.offsets = malloc(num_uchains * sizeof(int));
    if (!ci.counts || !ci.offsets) fail("Chain index allocation failed.");

    // first pass: count residues per chain
    for (int i = 0; i < n_res; i++) {
        ci.counts[res[i].chain_id]++;
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
        int cid = res[i].chain_id;
        ci.flat[write_pos[cid]++] = i;
    }
    free(write_pos);

    return ci;
}

// free chain index
void free_chain_index(ChainResIndex *ci) {
    free(ci->flat);
    free(ci->offsets);
    free(ci->counts);
}

// sorting helpers
int cmp_chainpair(const void *a, const void *b) {
    double diff = ((ChainPair *)b)->pinc_score - ((ChainPair *)a)->pinc_score;
    return (diff > 0) - (diff < 0);
}

int cmp_tokenpair(const void *a, const void *b) {
    double diff = ((TokenPair *)b)->contact_p - ((TokenPair *)a)->contact_p;
    return (diff > 0) - (diff < 0);
}

/* -------------------------------------------------------------------------- */
/* File parsers                                                               */
/* -------------------------------------------------------------------------- */

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

    char *ptr = strstr(buffer, "\"predicted_aligned_error\"");
    if (!ptr) ptr = strstr(buffer, "\"pae\"");
    if (!ptr) fail("no PAE matrix found in JSON.");

    while (*ptr && *ptr != '[') ptr++;
    if (!*ptr) fail("malformed JSON: expected '[' after PAE key.");

    int depth = 1;
    ptr++;

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

    *n_res = (int)sqrt((double)count);
    if ((size_t)(*n_res) * (size_t)(*n_res) != count) fail("parsed PAE matrix is not perfectly square.");
    return pae;
}

ParsedStruct parse_structure(const char *str_file) {
    FILE *f = fopen(str_file, "r");
    if (!f) fail("could not open structure file.");

    int cap = 1000;
    Residue *res_arr = malloc(cap * sizeof(Residue));
    if (!res_arr) fail("memory allocation failed for structure array.");
    int count = 0;

    int uc_cap = 16;
    ChainName *uchains = malloc(uc_cap * sizeof(ChainName));
    if (!uchains) fail("memory allocation failed for unique chains.");
    int num_uchains = 0;

    ResHash *rhash = rhash_create(4096);

    char line[MAX_LINE];
    bool is_cif = false;

    int c_chain = -1, c_seq = -1, c_elem = -1, c_x = -1, c_y = -1, c_z = -1;
    bool in_atom_site = false;
    bool cif_columns_validated = false;
    int cif_col_count = 0;

    int dc_cap = 16;
    ChainDotCounter *dot_counters = malloc(dc_cap * sizeof(ChainDotCounter));
    if (!dot_counters) fail("memory allocation failed for dot counters.");
    int num_chains_tracked = 0;

    while (fgets(line, sizeof(line), f)) {
        if (strchr(line, '\n') == NULL && !feof(f)) {
            int ch;
            while ((ch = fgetc(f)) != '\n' && ch != EOF);
        }

        if (strstr(line, "_atom_site.group_PDB")) {
            is_cif = true;
            in_atom_site = true;
        }

        if (is_cif && in_atom_site && line[0] == '_') {
            char *key = trim_in_place(line);
            if (strcmp(key, "_atom_site.label_asym_id") == 0) c_chain = cif_col_count;
            else if (strcmp(key, "_atom_site.label_seq_id") == 0) c_seq = cif_col_count;
            else if (strcmp(key, "_atom_site.type_symbol") == 0) c_elem = cif_col_count;
            else if (strcmp(key, "_atom_site.Cartn_x") == 0) c_x = cif_col_count;
            else if (strcmp(key, "_atom_site.Cartn_y") == 0) c_y = cif_col_count;
            else if (strcmp(key, "_atom_site.Cartn_z") == 0) c_z = cif_col_count;
            cif_col_count++;
            continue;
        }

        if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {

            // validate CIF columns once on the first ATOM/HETATM line;
            // fails early if any required _atom_site field is missing
            if (is_cif && !cif_columns_validated) {
                if (c_chain < 0 || c_seq < 0 || c_elem < 0 ||
                    c_x < 0 || c_y < 0 || c_z < 0) {
                    fclose(f);
                    free(dot_counters);
                    rhash_free(rhash);
                    free(res_arr);
                    free(uchains);
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
                    if ((size_t)snprintf(chain, sizeof(chain), "%s", tokens[c_chain]) >= sizeof(chain)) continue;
                }
                if (c_seq < t_count) {
                    if ((size_t)snprintf(seq_id, sizeof(seq_id), "%s", tokens[c_seq]) >= sizeof(seq_id)) continue;
                }
                if (c_elem < t_count) strncpy(elem, tokens[c_elem], 3);

                if (c_x < t_count && !safe_parse_double(tokens[c_x], &x)) continue;
                if (c_y < t_count && !safe_parse_double(tokens[c_y], &y)) continue;
                if (c_z < t_count && !safe_parse_double(tokens[c_z], &z)) continue;

            } else {
                size_t len = strlen(line);
                if (len < 54) continue;

                chain[0] = line[21];
                chain[1] = '\0';

                memcpy(seq_id, line + 22, 4);
                seq_id[4] = '\0';
                trim_in_place(seq_id);

                char coord_buf[16] = {0};

                memcpy(coord_buf, line + 30, 8);
                coord_buf[8] = '\0';
                if (!safe_parse_double(coord_buf, &x)) continue;

                memcpy(coord_buf, line + 38, 8);
                coord_buf[8] = '\0';
                if (!safe_parse_double(coord_buf, &y)) continue;

                memcpy(coord_buf, line + 46, 8);
                coord_buf[8] = '\0';
                if (!safe_parse_double(coord_buf, &z)) continue;

                if (len >= 78) {
                    memcpy(elem, line + 76, 2);
                    elem[2] = '\0';
                    trim_in_place(elem);
                }
            }

            trim_in_place(chain);
            trim_in_place(elem);
            if (strcasecmp(elem, "H") == 0) continue;

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
                if ((size_t)snprintf(seq_id, sizeof(seq_id), "%d", dc) >= sizeof(seq_id)) continue;
            }

            double mass = get_mass(elem);
            char token[MAX_CHAIN_LEN + MAX_SEQ_LEN + 1];
            if ((size_t)snprintf(token, sizeof(token), "%s_%s", chain, seq_id) >= sizeof(token)) continue;

            // register chain globally
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

            // check the last entry before falling through to the hash table
            int ri = -1;
            if (count > 0 && strcmp(res_arr[count - 1].token, token) == 0) {
                ri = count - 1;
            } else {
                ri = rhash_get(rhash, res_arr, token);
            }

            if (ri == -1) {
                // fail before allocating downstream matrices
                if (count >= MAX_TOKENS) {
                    fclose(f);
                    free(dot_counters);
                    rhash_free(rhash);
                    free(res_arr);
                    free(uchains);
                    fprintf(stderr, "Structure exceeds %d residues.\n", MAX_TOKENS);
                    fail("Structure too large.");
                }

                if (count >= cap) {
                    cap *= 2;
                    Residue *tmp = realloc(res_arr, cap * sizeof(Residue));
                    if (!tmp) fail("OOM structure array.");
                    res_arr = tmp;
                }
                snprintf(res_arr[count].chain, sizeof(res_arr[count].chain), "%s", chain);
                snprintf(res_arr[count].seq_id, sizeof(res_arr[count].seq_id), "%s", seq_id);
                snprintf(res_arr[count].token, sizeof(res_arr[count].token), "%s", token);
                res_arr[count].chain_id = cid;
                res_arr[count].x = 0;
                res_arr[count].y = 0;
                res_arr[count].z = 0;
                res_arr[count].mass_sum = 0;
                ri = count;

                rhash_insert(rhash, res_arr, ri);
                count++;
            }

            res_arr[ri].x += x * mass;
            res_arr[ri].y += y * mass;
            res_arr[ri].z += z * mass;
            res_arr[ri].mass_sum += mass;
        }
    }
    fclose(f);
    free(dot_counters);
    rhash_free(rhash);

    if (count < 10) fail("non-structure file or insufficient ATOM records.");

    // precompute COM to spare nested loop divisions
    for (int i = 0; i < count; i++) {
        if (res_arr[i].mass_sum > 0) {
            res_arr[i].x /= res_arr[i].mass_sum;
            res_arr[i].y /= res_arr[i].mass_sum;
            res_arr[i].z /= res_arr[i].mass_sum;
        }
    }

    ParsedStruct ps = {
        res_arr,
        count,
        uchains,
        num_uchains
    };
    return ps;
}

/* -------------------------------------------------------------------------- */
/* Main program                                                               */
/* -------------------------------------------------------------------------- */

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

    // validate extensions
    if (!has_extension(jsn_file, "json", NULL)) fail("JSON file must have .json extension.");
    if (!has_extension(str_file, "pdb", "cif")) fail("structure file must have .pdb/.cif extension.");
    if (access(jsn_file, F_OK) != 0 || access(str_file, F_OK) != 0) fail("input file(s) not found.");

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

    size_t rlen = strlen(results_dir);
    while (rlen > 0 && (results_dir[rlen - 1] == '/' || results_dir[rlen - 1] == '\\')) {
        results_dir[rlen - 1] = '\0';
        rlen--;
    }

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
    int n_res = ps.res_count;
    Residue *res = ps.res_arr;
    int num_uchains = ps.num_uchains;
    ChainName *uchains = ps.uchains;

    if (n_pae != n_res) {
        fprintf(stderr, "Dim mismatch. PAE: %dx%d. Struct: %dx%d.\n", n_pae, n_pae, n_res, n_res);
        fail("Check that structure and PAE correspond to the same model.");
    }

    // allocate packed symmetric upper-triangular matrices
    size_t sym_size = ((size_t)n_res * (size_t)(n_res + 1)) / 2;
    double *dist_mat = malloc(sym_size * sizeof(double));
    double *cont_sym = malloc(sym_size * sizeof(double));
    if (!dist_mat || !cont_sym) fail("Symmetric matrix allocation failed.");

    // compute distances and symmetric contact probabilities
    for (size_t i = 0; i < (size_t)n_res; i++) {
        size_t diag_idx = PACKED_IDX(i, i, n_res);
        dist_mat[diag_idx] = 0.0;
        cont_sym[diag_idx] = 1.0;

        for (size_t j = i + 1; j < (size_t)n_res; j++) {
            double dx = res[i].x - res[j].x;
            double dy = res[i].y - res[j].y;
            double dz = res[i].z - res[j].z;
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            
            size_t p_idx = PACKED_IDX(i, j, n_res);
            dist_mat[p_idx] = d;

            double pae_ij = pae_mat[i * n_res + j];
            double p_ij = 0.0;
            if (pae_ij > 0) {
                double v_int = vol_intersect(pae_ij, d);
                double v_unc = vol_sphere(pae_ij);
                p_ij = (v_unc > 0) ? (v_int / v_unc) : 0.0;
                if (!isfinite(p_ij)) p_ij = 0.0;
                if (p_ij > 1.0) p_ij = 1.0;
            }

            double pae_ji = pae_mat[j * n_res + i];
            double p_ji = 0.0;
            if (pae_ji > 0) {
                double v_int = vol_intersect(pae_ji, d);
                double v_unc = vol_sphere(pae_ji);
                p_ji = (v_unc > 0) ? (v_int / v_unc) : 0.0;
                if (!isfinite(p_ji)) p_ji = 0.0;
                if (p_ji > 1.0) p_ji = 1.0;
            }

            cont_sym[p_idx] = (p_ij + p_ji) / 2.0;
        }
    }

    // chain index: single flat array with offset/count per chain
    ChainResIndex ci = build_chain_index(res, n_res, num_uchains);

    size_t total_pairs = ((size_t)num_uchains * ((size_t)num_uchains - 1)) / 2;
    ChainPair *cpairs = calloc(total_pairs, sizeof(ChainPair));
    if (!cpairs && total_pairs > 0) fail("OOM chain pairs.");

    size_t cp_idx = 0;

    for (int i = 0; i < num_uchains; i++) {
        for (int j = i + 1; j < num_uchains; j++) {
            snprintf(cpairs[cp_idx].chain1, sizeof(cpairs[cp_idx].chain1), "%s", uchains[i]);
            snprintf(cpairs[cp_idx].chain2, sizeof(cpairs[cp_idx].chain2), "%s", uchains[j]);

            double sum = 0.0;
            int count = 0;

            // iterate strictly over matching inter-chain indices
            for (int u = 0; u < ci.counts[i]; u++) {
                int ri = ci.flat[ci.offsets[i] + u];
                for (int v = 0; v < ci.counts[j]; v++) {
                    int rj = ci.flat[ci.offsets[j] + v];
                    
                    size_t s_idx = SYM_IDX(ri, rj, n_res);
                    if (dist_mat[s_idx] < CONTACT_RADIUS) {
                        sum += cont_sym[s_idx];
                        count++;
                    }
                }
            }
            cpairs[cp_idx].pinc_score = (count > 0) ? (sum / count) : 0.0;
            cp_idx++;
        }
    }

    qsort(cpairs, total_pairs, sizeof(ChainPair), cmp_chainpair);

    // output default
    char pinc_csv[1024];
    if ((size_t)snprintf(pinc_csv, sizeof(pinc_csv), "%s_Pinc.csv", out_path) >= sizeof(pinc_csv)) fail("Path long.");

    FILE *fp = fopen(pinc_csv, "w");
    if (fp) {
        fprintf(fp, "chain1,chain2,Pinc\n");
        for (size_t i = 0; i < total_pairs; i++) {
            fprintf(fp, "%s,%s,%.4f\n", cpairs[i].chain1, cpairs[i].chain2, cpairs[i].pinc_score);
        }
        fclose(fp);
    } else fprintf(stderr, "warning: could not write %s\n", pinc_csv);

    // optional --all output
    if (opt_all) {
        double *cont_asym = malloc((size_t)n_res * (size_t)n_res * sizeof(double));
        if (!cont_asym) fail("Full contact matrix allocation failed.");

        for (size_t i = 0; i < (size_t)n_res; i++) {
            cont_asym[i * n_res + i] = 1.0;
            for (size_t j = i + 1; j < (size_t)n_res; j++) {
                double d = dist_mat[PACKED_IDX(i, j, n_res)];

                double pae_ij = pae_mat[i * n_res + j];
                double p_ij = 0.0;
                if (pae_ij > 0) {
                    double v_int = vol_intersect(pae_ij, d);
                    double v_unc = vol_sphere(pae_ij);
                    p_ij = (v_unc > 0) ? (v_int / v_unc) : 0.0;
                    if (!isfinite(p_ij)) p_ij = 0.0;
                    if (p_ij > 1.0) p_ij = 1.0;
                }
                cont_asym[i * n_res + j] = p_ij;

                double pae_ji = pae_mat[j * n_res + i];
                double p_ji = 0.0;
                if (pae_ji > 0) {
                    double v_int = vol_intersect(pae_ji, d);
                    double v_unc = vol_sphere(pae_ji);
                    p_ji = (v_unc > 0) ? (v_int / v_unc) : 0.0;
                    if (!isfinite(p_ji)) p_ji = 0.0;
                    if (p_ji > 1.0) p_ji = 1.0;
                }
                cont_asym[j * n_res + i] = p_ji;
            }
        }

        char jsn_out[1024];
        if ((size_t)snprintf(jsn_out, sizeof(jsn_out), "%s_contact_probability.json", out_path) >= sizeof(jsn_out))
            fail("JSON output path too long.");

        FILE *fj = fopen(jsn_out, "w");
        if (fj) {
            fprintf(fj, "{\n  \"token_chain_ids\":[\n");
            for (int i = 0; i < n_res; i++) {
                fprintf(fj, "    \"%s\"%s\n", res[i].chain, (i < n_res - 1) ? "," : "");
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
    if (opt_pairlist) {
        int cap = 1000, pcnt = 0;
        TokenPair *tpairs = malloc(cap * sizeof(TokenPair));
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

                        size_t p_idx = PACKED_IDX(lo, hi, n_res);
                        double cp = cont_sym[p_idx];
                        double dd = dist_mat[p_idx];

                        if (cp > 0 && dd < CONTACT_RADIUS) {
                            if (pcnt >= cap) {
                                cap *= 2;
                                TokenPair *tmp = realloc(tpairs, cap * sizeof(TokenPair));
                                if (!tmp) fail("OOM pairlist array.");
                                tpairs = tmp;
                            }
                            snprintf(tpairs[pcnt].token1, sizeof(tpairs[pcnt].token1), "%s", res[lo].token);
                            snprintf(tpairs[pcnt].token2, sizeof(tpairs[pcnt].token2), "%s", res[hi].token);
                            tpairs[pcnt].contact_p = cp;
                            tpairs[pcnt].distance = dd;
                            pcnt++;
                        }
                    }
                }
            }
        }

        qsort(tpairs, pcnt, sizeof(TokenPair), cmp_tokenpair);

        char pair_csv[1024];
        if ((size_t)snprintf(pair_csv, sizeof(pair_csv), "%s_pairlist.csv", out_path) >= sizeof(pair_csv)) fail("Path.");

        FILE *ft = fopen(pair_csv, "w");
        if (ft) {
            fprintf(ft, "token1,token2,contact_p,distance\n");
            for (int i = 0; i < pcnt; i++) {
                fprintf(ft, "%s,%s,%.4f,%.4f\n", tpairs[i].token1, tpairs[i].token2, tpairs[i].contact_p, tpairs[i].distance);
            }
            fclose(ft);
        }
        free(tpairs);
    }

    // cleanup
    free_chain_index(&ci);
    free(res);
    free(uchains);
    free(pae_mat);
    free(dist_mat);
    free(cont_sym);
    free(cpairs);

    return 0;
}