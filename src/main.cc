// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Command-line interface.

#include <unistd.h>
#include <errno.h>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <string>
#include <complex>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "utils.h"
#include <map>
#include "gsl_all.h"
#include "prm.h"
#include "npf.h"
#include "hxs.h"
#include "fi.h"
#include "pdf.h"
#include "getopt.h"
#include "globals.h"

using namespace std;

#define UNITS (result_type == total ? "pb" : "pb/GeV")
#define get_option(name) (arguments[name] == "" ? config[name] : arguments[name])

static char * my_get_current_dir_name(void)
{
    // This is necessary for systems without GNU extensions.
    size_t output_size = 100;
    char *output = (char*)malloc(output_size);
    output = getcwd(output, output_size);
    while (output == NULL && errno == ERANGE) {
        output_size += 100;
        output = getcwd(output, output_size);
    }
    return output;
}

static int pdg_to_internal_id(int pdg_code)
{
    switch (pdg_code) {
    case 1000012: // ~nu_eL
        return 10;
        break;
    case 1000014: // ~nu_muL
        return 11;
        break;
    case 1000016: // ~nu_tau
        return 12;
        break;
    case 1000011: // ~e_L
        return 13;
        break;
    case 2000011: // ~e_R
        return 14;
        break;
    case 1000013: // ~mu_L
        return 15;
        break;
    case 2000013: // ~mu_R
        return 16;
        break;
    case 1000015: // ~tau_1
        return 17;
        break;
    case 2000015: // ~tau_2
        return 18;
        break;
    case 1000022: // ~chi_10
        return 0;
        break;
    case 1000023: // ~chi_20
        return 1;
        break;
    case 1000025: // ~chi_30
        return 2;
        break;
    case 1000035: // ~chi_40
        return 3;
        break;
    case 1000024: // ~chi_1+
        return 4;
        break;
    case 1000037: // ~chi_2+
        return 5;
        break;
    case 12: // nu_eL
        return 20;
        break;
    case 14: // nu_muL
        return 21;
        break;
    case 16: // nu_tauL
        return 22;
        break;
    case 11: // eL
        return 23;
        break;
    case 13: // nu_muL
        return 24;
        break;
    case 15: // nu_tauL
        return 25;
        break;
    case 1000001: // ~d_L
        return 31;
        break;
    case 2000001: // ~d_R
        return 32;
        break;
    case 1000002: // ~u_L
        return 33;
        break;
    case 2000002: // ~u_R
        return 34;
        break;
    case 1000003: // ~s_L
        return 35;
        break;
    case 2000003: // ~s_R
        return 36;
        break;
    case 1000004: // ~c_L
        return 37;
        break;
    case 2000004: // ~c_R
        return 38;
        break;
    case 1000005: // ~b_1
        return 39;
        break;
    case 2000005: // ~b_2
        return 40;
        break;
    case 1000006: // ~t_1
        return 41;
        break;
    case 2000006: // ~t_2
        return 42;
        break;
    case 1000021: // ~g
        return 30;
        break;
    default:
        fprintf(stderr, "error: Wrong PDG code.\n");
        exit(1);
    }
}

static void print_banner()
{
    printf("\n"                                                  \
           "  Welcome to Resummino " RESUMMINO_VERSION "\n"             \
           "  Copyright 2008-2010 Jonathan Debove.\n"                   \
           "  Copyright 2011-2013 David R. Lamprea and Marcel Rothering.\n" \
           "\n"                                                         \
           "  Licensed under the EUPL 1.1 or later.\n"                  \
           "  See http://www.resummino.org/ for more information.\n"    \
           "\n"                                                         \
           "  If using Resummino for slepton pair production, please cite.\n" \
           "  - G. Bozzi, B. Fuks, M. Klasen, PRD 74, 015001 (2006); NPB 777, 157 (2007); NPB 794, 46 (2007).\n" \
           "  - B. Fuks, M. Klasen, D. R. Lamprea, M. Rothering, arXiv:1304.0790 [hep-ph].\n" \
           "\n"                                                         \
           "  For gaugino pair production, please cite:\n"              \
           "  - J. Debove, B. Fuks, M. Klasen, PLB 688, 208 (2010); NPB 842, 51 (2011); NPB 849, 64 (2011).\n" \
           "  - B. Fuks, M. Klasen, D. R. Lamprea, M. Rothering, JHEP 1210, 081 (2012); arXiv:1304.0790 [hep-ph].\n" \
           "\n");
}

int main(int argc, char *argv[])
{
    string input_file("resummino.in");
    string log_file("");
    int stop_after_lo = 0;
    int stop_after_nlo = 0;
    enum {total, pt, m, ptj} result_type;
    Parameters *jp = new Parameters();

    // A value different than "" for an argument means that the user has set
    // that argment through the command-line.
    map<string, string> arguments;
    arguments["key"] = "";
    arguments["m"] = "";
    arguments["pt"] = "";
    arguments["pdfset_lo"] = "";
    arguments["pdfset_nlo"] = "";
    arguments["mu_f"] = "";
    arguments["mu_r"] = "";
    arguments["particle1"] = "";
    arguments["particle2"] = "";
    arguments["output_file"] = "";
    arguments["paramfile"] = "";
    arguments["center_of_mass_energy"] = "";

    for (;;) {
        static struct option long_options[] = {
            {"version", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {"parameter-log", required_argument, 0, 'p'},
            {"lo", no_argument, &stop_after_lo, 'l'},
            {"nlo", no_argument, &stop_after_nlo, 'n'},
            {"invariant-mass", required_argument, 0, 'm'},
            {"key", required_argument, 0, 'k'},
            {"transverse-momentum", required_argument, 0, 't'},
            {"pdfset_lo", required_argument, 0, 'a'},
            {"pdfset_nlo", required_argument, 0, 'b'},
            {"mu_f", required_argument, 0, 'f'},
            {"mu_r", required_argument, 0, 'r'},
            {"particle1", required_argument, 0, 'c'},
            {"particle2", required_argument, 0, 'd'},
            {"output_file", required_argument, 0, 'o'},
            {"slha", required_argument, 0, 'i'},
            {"center_of_mass_energy", required_argument, 0, 'e'},
            {0, 0, 0, 0}
        };
        int option_index;
        int c;
        c = getopt_long(argc, argv, "vhp:lnm:t:a:b:f:r:k:",
                        long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'v':
            fprintf(stdout, "resummino " RESUMMINO_VERSION "\n");
            exit(0);
            break;
        case 'h':
            fprintf(stdout,
                    "\n  Please see http://www.resummino.org/ for instructions.\n\n");
            exit(0);
            break;
        case 'p':
            log_file = optarg;
            break;
        case 'm':
            arguments["m"] = optarg;
            break;
        case 't':
            arguments["pt"] = optarg;
            break;
        case 'a':
            arguments["pdfset_lo"] = optarg;
            break;
        case 'b':
            arguments["pdfset_nlo"] = optarg;
            break;
        case 'k':
            arguments["key"] = optarg;
            break;
        case 'f':
            arguments["mu_f"] = optarg;
            break;
        case 'r':
            arguments["mu_r"] = optarg;
            break;
        case 'c':
            arguments["particle1"] = optarg;
            break;
        case 'd':
            arguments["particle2"] = optarg;
            break;
        case 'o':
            arguments["output_file"] = optarg;
            break;
        case 'i':
            arguments["slha"] = optarg;
            break;
        case 'e':
            arguments["center_of_mass_energy"] = optarg;
            break;
        case '?':
            fprintf(stderr, "error: Unknown option given.\n");
            exit(1);
            break;
        default:
            break;
        }
    }

    print_banner();

    // There can be either 0 or 1 left arguments (input file name).
    if (optind > argc) {
        fprintf(stderr, "error: too many arguments.\n");
        exit(1);
    } else if (optind == argc - 1) {
        input_file = argv[optind];
    } else {
        fprintf(stderr, "warning: you did not specify an input file: "
                "using default 'resummino.in'.\n");
    }

    map<string, string> config = jp->read_input_file(input_file);

    // Changes working directory to file folder for relative paths.
    string old_wd(my_get_current_dir_name());
    //string real_name(canonicalize_file_name(input_file.c_str()));
    string real_name(realpath(input_file.c_str(), NULL));
    chdir(real_name.substr(0, real_name.rfind("/")).c_str());

    // Sets collider type.
    if (config["collider_type"] == "proton-proton") {
        jp->ic = 0;
    } else if (config["collider_type"] == "proton-antiproton"
               || config["collider_type"] == "antiproton-proton") {
        jp->ic = 1;
    } else {
        fprintf(stderr, "error: collider type not recognized.\n");
        exit(1);
    }

    // Sets computation type.
    if (config["result"] == "total") {
        result_type = total;
    } else if (config["result"] == "pt") {
        result_type = pt;
    } else if (config["result"] == "m") {
        result_type = m;
    } else if (config["result"] == "ptj") {
        result_type = ptj;
    } else {
        fprintf(stderr, "error: computation '%s' not recognized. "
                "'total', 'pt', 'ptj' or 'm' expected.\n",
                config["result"].c_str());
        exit(1);
    }

    // Sets process.
    jp->in1 = 0;
    jp->in2 = 0;

    int pdg_p1 = atoi(get_option("particle1").c_str());
    int pdg_p2 = atoi(get_option("particle2").c_str());

    // For the chargino-neutralino case. Anti-chargino is negatively charged.
    if (pdg_p1 < 0 && pdg_p2 > 0 && pdg_to_internal_id(abs(pdg_p1)) <= 5
            && pdg_to_internal_id(abs(pdg_p2)) <= 5) {
        jp->out1 = pdg_to_internal_id(abs(pdg_p2));
        jp->out2 = pdg_to_internal_id(abs(pdg_p1));
    }
    // Same as before. Just p1 and p2 exchanged.
    else if (pdg_p2 < 0 && pdg_p1 > 0 && pdg_to_internal_id(abs(pdg_p1)) <= 5
             && pdg_to_internal_id(abs(pdg_p2)) <= 5) {
        jp->out1 = pdg_to_internal_id(abs(pdg_p1));
        jp->out2 = pdg_to_internal_id(abs(pdg_p2));
    }
    // For the sleptons.
    else if (pdg_p1 < 0 && pdg_p2 > 0) {
        jp->out1 = pdg_to_internal_id(abs(pdg_p1));
        jp->out2 = pdg_to_internal_id(abs(pdg_p2));
    }
    // Same as above. Just p1 and p2 changed.
    else if (pdg_p2 < 0 && pdg_p1 > 0) {
        jp->out1 = pdg_to_internal_id(abs(pdg_p2));
        jp->out2 = pdg_to_internal_id(abs(pdg_p1));
    }
    // For the chargino-neutralino case. Chargino is positively charged.
    else if (pdg_p1 > 0 && pdg_p2 > 0 && pdg_to_internal_id(abs(pdg_p1)) <  4
             && pdg_to_internal_id(abs(pdg_p2)) >= 4
             && pdg_to_internal_id(abs(pdg_p2)) <= 5) {
        jp->out1 = pdg_to_internal_id(abs(pdg_p2));
        jp->out2 = pdg_to_internal_id(abs(pdg_p1));
    }
    // For the chargino-neutralino case. Chargino is positively charged.
    else if (pdg_p2 > 0 && pdg_p1 > 0 && pdg_to_internal_id(abs(pdg_p2)) <  4
             && pdg_to_internal_id(abs(pdg_p1)) >= 4
             && pdg_to_internal_id(abs(pdg_p1)) <= 5) {
        jp->out1 = pdg_to_internal_id(abs(pdg_p1));
        jp->out2 = pdg_to_internal_id(abs(pdg_p2));
    } else {
        fprintf(stderr, "error: Process not recognized.\n");
        exit(1);
    }

    jp->sh = pow2(atof(config["center_of_mass_energy"].c_str()));

    // Diagonal CKM matrix.
    for (int i0 = 0; i0 < 3; i0++) {
        for (int i1 = 0; i1 < 3; i1++) {
            jp->ckm[i0][i1] = (i0 == i1 ? 1.0 : 0.0);
        }
    }

    // Sets Yukawa couplings to 0 for (d,s,b,u,c).
    for (int i0 = 0; i0 < 6; i0++) {
        jp->yq[i0] = 0.;
    }

    // Reads parameters from SLHA.
    jp->read_slha(get_option("slha").c_str());

    // Sets SM and SUSY couplings.
    jp->set_couplings();

    // Sets the transverse momentum.
    if (get_option("pt") == "auto" || get_option("pt") == "") {
        jp->pts = -1;
    } else {
        jp->pts = pow2(atof(get_option("pt").c_str()));
    }

    // Sets the invariant mass.
    if (get_option("M") == "auto" || get_option("M") == "") {
        // Sets automatic invariant mass.
        if (jp->out1 < 10) {
            jp->mis = pow2(jp->mCH[jp->out1] + jp->mCH[jp->out2]);
        } else if (jp->out1 >= 10 && jp->out1 < 20) {
            jp->mis = pow2(jp->mSL[jp->out1 - 10] + jp->mSL[jp->out2 - 10]);
        }
    } else {
        // Sets user-defined invariant mass.
        jp->mis = pow2(atof(get_option("M").c_str()));
    }

    // Sets scales.
    jp->mufs = pow2(atof(get_option("mu_f").c_str())) * jp->mis;
    jp->murs = pow2(atof(get_option("mu_r").c_str())) * jp->mis;

    // Integration parameters.
    jp->precision = atof(config["precision"].c_str());
    jp->max_iters = atoi(config["max_iters"].c_str());
    jp->fout = stderr;

    // Restores working directory and saves the parameter log.
    chdir(old_wd.c_str());
    if (log_file != "") {
        jp->write_log(log_file.c_str());
    }

    // Computes the processes at LO, NLO und NLO+NLL.

    double res0 = 0.0;
    double err0 = 0.0;
    double res1 = 0.0;
    double err1 = 0.0;
    double res2 = 0.0;
    double err2 = 0.0;
    double res3 = 0.0;
    double err3 = 0.0;
    double res4 = 0.0;
    double err4 = 0.0;
    double chisq = 0.0;

    // Sets the PDF at LO.
    if (config["pdf_format"] == "lhpdf") {
        LHAPDF::initPDFSet(config["pdf_lo"], LHAPDF::LHPDF,
                           atoi(get_option("pdfset_lo").c_str()));
    } else {
        LHAPDF::initPDFSet(config["pdf_lo"], LHAPDF::LHGRID,
        atoi(get_option("pdfset_lo").c_str()));
    }

    // LO calculation
    switch (result_type) {
    case total:
        hadronic_xs(res0, err0, chisq, 0, -0, jp);
        break;
    case pt:
        hadronic_xs_dPT2(res0, err0, chisq, 0, -0, jp);
        res0 *= 2.0 * sqrt(jp->pts);
        err0 *= 2.0 * sqrt(jp->pts);
        break;
    case ptj:
        hadronic_xs_dPT2(res0, err0, chisq, 0, -0, jp);
        res0 *= 2.0 * sqrt(jp->pts);
        err0 *= 2.0 * sqrt(jp->pts);
        break;
    case m:
        hadronic_xs_dlnM2(res0, err0, chisq, 0, -0, jp);
        res0 *= 2.0 / sqrt(jp->mis);
        err0 *= 2.0 / sqrt(jp->mis);
        break;
    }
    res1 = res0;
    err1 = err0;

    printf("LO = (%.7e +- %.7e) %s [chi**2 = %.7e]\n",
           res1, sqrt(err1), UNITS, chisq);

    /*
    if (stop_after_lo) {
        // Prints again the results.
        printf("\nResults:\n"
               "LO = (%.7e +- %.7e) %s\n",
               res1, sqrt(err1), UNITS);

        exit(0);
    }
    */
    if (!stop_after_lo) {

        // Sets the PDF at NLO.
        if (config["pdf_format"] == "lhpdf") {
            LHAPDF::initPDFSet(config["pdf_nlo"], LHAPDF::LHPDF,
                               atoi(get_option("pdfset_nlo").c_str()));
        } else {
            LHAPDF::initPDFSet(config["pdf_nlo"], LHAPDF::LHGRID,
            atoi(get_option("pdfset_nlo").c_str()));
        }
        // NLO
        for (int i0 = 0; i0 < 5; i0++) {
            switch (result_type) {
            case total:
                hadronic_xs(res0, err0, chisq, i0, -0, jp);
                break;
            case pt:
                hadronic_xs_dPT2(res0, err0, chisq, i0, -0, jp);
                res0 *= 2.0 * sqrt(jp->pts);
                err0 *= 2.0 * sqrt(jp->pts);
                break;
            case ptj:
                hadronic_xs_dPT2(res0, err0, chisq, i0, -0, jp);
                res0 *= 2.0 * sqrt(jp->pts);
                err0 *= 2.0 * sqrt(jp->pts);
                break;
            case m:
                hadronic_xs_dlnM2(res0, err0, chisq, i0, -0, jp);
                res0 *= 2.0 / sqrt(jp->mis);
                err0 *= 2.0 / sqrt(jp->mis);
                break;
            }
            if (i0 == 0)    {
                res2 = res0;
                err2 = err0;
            } else if (i0 < 5) {
                res2 += res0;
                err2 += err0;
            }
        }

        printf("NLO = (%.7e +- %.7e) %s [chi**2 = %.7e]\n",
               res2, sqrt(err2), UNITS, chisq);
        /*
        if (stop_after_nlo) {
            // Prints again the results.
            printf("\nResults:\n"
                   "LO = (%.7e +- %.7e) %s\n"
                   "NLO = (%.7e +- %.7e) %s\n",
                   res1, sqrt(err1), UNITS,
                   res2, sqrt(err2), UNITS);
            exit(0);
        }
        */
        if (!stop_after_nlo) {
            // NLO+NLL
            pdfFit(jp->a1min, jp->afit, jp->mis / jp->sh, jp->mufs);
            for (int i0 = 5; i0 < 6; i0++) {
                switch (result_type) {
                case total:
                    hadronic_xs(res0, err0, chisq, i0, -0, jp);
                    break;
                case pt:
                    hadronic_xs_dPT2(res0, err0, chisq, i0, -0, jp);
                    res0 *= 2.0 * sqrt(jp->pts);
                    err0 *= 2.0 * sqrt(jp->pts);
                    break;
                case ptj:
                    hadronic_xs_dPT2(res0, err0, chisq, i0, -0, jp);
                    res0 *= 2.0 * sqrt(jp->pts);
                    err0 *= 2.0 * sqrt(jp->pts);
                    break;
                case m:
                    hadronic_xs_dlnM2(res0, err0, chisq, i0, -0, jp);
                    res0 *= 2.0 / sqrt(jp->mis);
                    err0 *= 2.0 / sqrt(jp->mis);
                    break;
                }
                res3 = res2 + res0;
                err3 = err2 + err0;
            }

            printf("NLO+NLL = (%.7e +- %.7e) %s [chi**2 = %.7e]\n",
                   res3, sqrt(err3), UNITS, chisq);

            // Joint resummation.
            if (result_type == ptj) {
                hadronic_xs_dPT2(res0, err0, chisq, 6, -0, jp);
                res0 *= 2.0 * sqrt(jp->pts);
                err0 *= 2.0 * sqrt(jp->pts);
                res4 = res2 + res0;
                err4 = err2 + err0;
                printf("NLO+NLLj = (%.7e +- %.7e) %s [chi**2 = %.7e]\n",
                       res4, sqrt(err4), UNITS, chisq);
            }

        }

    }

    // Prints again the results.
    if (result_type == ptj) {
        printf("\nResults:\n"
               "LO = (%.7e +- %.7e) %s\n"
               "NLO = (%.7e +- %.7e) %s\n"
               "NLO+NLL = (%.7e +- %.7e) %s\n"
               "NLO+NLLj = (%.7e +- %.7e) %s\n",
               res1, sqrt(err1), UNITS,
               res2, sqrt(err2), UNITS,
               res3, sqrt(err3), UNITS,
               res4, sqrt(err4), UNITS);
    } else {
        printf("\nResults:\n"
               "LO = (%.7e +- %.7e) %s\n"
               "NLO = (%.7e +- %.7e) %s\n"
               "NLO+NLL = (%.7e +- %.7e) %s\n",
               res1, sqrt(err1), UNITS,
               res2, sqrt(err2), UNITS,
               res3, sqrt(err3), UNITS);
    }

    // Prints results to a file for a subsequent analysis.
    // The current format is JSON.
    if (get_option("output_file") != "") {
        FILE *fout = fopen(get_option("output_file").c_str(), "w");
        if (!fout) {
            fprintf(stderr,
                    "error: Could not open output file %s.\n",
                    get_option("output_file").c_str());
            return 1;
        }
        fprintf(fout, "{\n"                                        \
                "\"key\": \"%s\",\n"                               \
                "\"pt\": %s,\n"                                    \
                "\"m\": %s,\n"                                     \
                "\"pdflo\": \"%s\",\n"                             \
                "\"pdfsetlo\": %s,\n"                              \
                "\"pdfnlo\": \"%s\",\n"                            \
                "\"pdfsetnlo\": %s,\n"                             \
                "\"muf\": %s,\n"                                   \
                "\"mur\": %s,\n"                                   \
                "\"lo\": %.7e,\n"                                  \
                "\"nlo\": %.7e,\n"                                 \
                "\"nll\": %.7e,\n"                                 \
                "\"nllj\": %.7e,\n"                                \
                "\"units\": \"%s\"\n"                              \
                "}",
                get_option("key").c_str(),
                (get_option("pt") == "auto" ? "-1" : get_option("pt").c_str()),
                (get_option("M") == "auto" ? "-1" : get_option("M").c_str()),
                get_option("pdf_lo").c_str(),
                get_option("pdfset_lo").c_str(),
                get_option("pdf_nlo").c_str(),
                get_option("pdfset_nlo").c_str(),
                get_option("mu_f").c_str(),
                get_option("mu_r").c_str(),
                res1,
                res2,
                res3,
                res4,
                UNITS
            );
    }

    return 0;
}
