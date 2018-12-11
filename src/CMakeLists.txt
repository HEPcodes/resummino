set(resummino_SRCS
    main.cc
    hxs.cc
    hxs_dpt2.cc
    hxs_dlnm2.cc
    prm.cc
    pxs.cc
    fi.cc
    mbn.cc
    mtr.cc
    mre.cc
    mbo.cc
    mbu.cc
    mct.cc
    npf.cc
    vir.cc
    pdf.cc
    res.cc
    mth.cc
)

add_executable(resummino ${resummino_SRCS})
target_link_libraries(resummino ${GSL_LIB} ${CBLAS_LIB} qcdloop ff ${LHAPDF_LIB})
install(TARGETS resummino RUNTIME DESTINATION bin)
