add_requires("benchmark")

target("bm_all")
    set_kind("binary")
    add_files("main.cpp")
    add_files(
        "../src/morton.c",
        "../src/octree.c",
        "../src/radix_tree.c"
    )
    add_includedirs("../src")
    add_defines("__tb_prefix__=\"hello-tbox\"")
    add_packages("tbox", "openmp", "cglm")
    add_packages("benchmark")
