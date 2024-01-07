add_requires("tbox >1.6.3", "openmp", "cglm")

target("hello-tbox")
    set_kind("binary")
    add_defines("__tb_prefix__=\"hello-tbox\"")
    add_packages("tbox", "openmp", "cglm")
    add_headerfiles("*.h")
    add_files("*.c")

