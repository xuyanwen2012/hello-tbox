#include <benchmark/benchmark.h>
#include <cglm/cglm.h>
#include <omp.h>
#include <tbox/tbox.h>

#include <cstdlib>

#include "morton.h"

namespace bm = benchmark;

constexpr auto kN = 10'000'000;
constexpr auto kMin = 0.0f;
constexpr auto kMax = 1024.0f;
constexpr auto kRange = kMax - kMin;

static void GenerateInput(vec4* data, tb_float_t min, tb_float_t max) {
  tb_srand(114514);
  const tb_float_t range = max - min;

#pragma omp parallel for schedule(static)
  for (tb_int_t i = 0; i < kN; i++) {
    data[i][0] = (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min;
    data[i][1] = (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min;
    data[i][2] = (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min;
    data[i][3] = 1.0f;
  }
}

static void BM_Morton32(bm::State& st) {
  const auto num_threads = st.range(0);

  // use malloc to allocate memory
  vec4* data = new vec4[kN];
  tb_uint32_t* morton_keys = new tb_uint32_t[kN];

  GenerateInput(data, kMin, kMax);

  omp_set_num_threads(num_threads);
  for (auto _ : st) {
    convert_xyz_to_morton_code(data, morton_keys, kN, kMin, kRange);
    bm::DoNotOptimize(morton_keys[0]);
  }

  delete[] data;
  delete[] morton_keys;
}

static void BM_QSort(bm::State& st) {
  //   const auto num_threads = st.range(0);

  // use malloc to allocate memory
  tb_uint32_t* morton_keys = new tb_uint32_t[kN];

  std::generate_n(morton_keys, kN, [n = kN]() mutable { return n--; });

  //   omp_set_num_threads(num_threads);
  for (auto _ : st) {
    qsort(morton_keys, kN, sizeof(tb_uint32_t), compare_uint32_t);
  }

  delete[] morton_keys;
}

BENCHMARK(BM_Morton32)
    ->RangeMultiplier(2)
    ->Range(1, 48)
    ->Unit(bm::kMillisecond);

BENCHMARK(BM_QSort)->Unit(bm::kMillisecond);

BENCHMARK_MAIN();
