// Host run of the CUDA port: reads truth CSVs (id,cat,r1(3),r2(3),tof,dir,nrev,...),
// solves the nrev = 0 rows and writes "id status v1(3) v2(3)" with %a doubles.
#include "../excali_zero_rev.cuh"
#include <cstdio>
#include <cstdlib>
#include <cstring>
int main(int argc, char **argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: host_driver in.csv out.txt\n"); return 1; }
    FILE *in = std::fopen(argv[1], "r"), *out = std::fopen(argv[2], "w");
    if (!in || !out) return 1;
    char line[4096];
    if (!std::fgets(line, sizeof line, in)) return 1;      // header
    while (std::fgets(line, sizeof line, in)) {
        double v[13]; int n = 0; char *p = line;
        while (n < 13) { char *e; v[n++] = std::strtod(p, &e); p = std::strchr(e, ','); if (!p) break; ++p; }
        if (n < 13 || (int)v[10] != 0) continue;
        double r1[3] = {v[2], v[3], v[4]}, r2[3] = {v[5], v[6], v[7]}, v1[3], v2[3];
        int st, half;
        exc::excali_zero_rev(r1, r2, v[8], (int)v[9], v1, v2, st, half);
        std::fprintf(out, "%ld %d %a %a %a %a %a %a\n", (long)v[0], st, v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
    }
    return 0;
}
