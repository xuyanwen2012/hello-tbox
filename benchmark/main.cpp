#include <benchmark/benchmark.h>
#include <cglm/cglm.h>
#include <omp.h>
#include <tbox/tbox.h>

#include <cstdlib>

#include "morton.h"

namespace bm = benchmark;

constexpr auto kN = 10'000'000;

static void BM_Morton32(bm::State& st) {
  const auto num_threads = st.range(0);

  omp_set_num_threads(num_threads);

  //   // allocate kN vec4 elements
  //   vec4* data = tb_nalloc_type(kN, vec4);
  //   tb_assert(data);

  //   // allocate kN morton code
  //   tb_uint32_t* morton_keys = tb_nalloc_type(kN, tb_uint32_t);
  //   tb_assert(morton_keys);

  // use malloc to allocate memory
  vec4* data = new vec4[kN];
  tb_uint32_t* morton_keys = new tb_uint32_t[kN];

  tb_srand(114514);

  const tb_float_t min_coord = 0.0f;
  const tb_float_t max_coord = 1024.0f;
  const tb_float_t range = max_coord - min_coord;

  // // #pragma omp parallel for schedule(static)
  for (tb_int_t i = 0; i < kN; i++) {
    data[i][0] =
        (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min_coord;
    data[i][1] =
        (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min_coord;
    data[i][2] =
        (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min_coord;
    data[i][3] = 1.0f;
  }

  for (auto _ : st) {
    bm::DoNotOptimize(data[0]);

    convert_xyz_to_morton_code(data, morton_keys, kN, min_coord, range);
  }

  delete[] data;
  delete[] morton_keys;

  //   tb_free(data);
  //   tb_free(morton_keys);
}

BENCHMARK(BM_Morton32)
    ->RangeMultiplier(2)
    ->Range(1, 48)
    ->Unit(bm::kMillisecond);

BENCHMARK_MAIN();
