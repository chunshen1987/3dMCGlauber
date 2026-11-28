// Copyright @ Chun Shen 2018

#include "EventGenerator.h"

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>

namespace MCGlb {

namespace {

constexpr char minBiasTableVersion[] = "# 3dMCGlauber minimum-bias Nstrings v1";

std::string minBiasTableParameters(
    const Parameters& parameters, real collisionEnergy) {
    // Exclude production event selection, output settings, and the random seed.
    const std::set<std::string> excluded {
        "cenMin",
        "cenMax",
        "b_min",
        "b_max",
        "NpartMin",
        "NpartMax",
        "seed",
        "only_event_statistics",
        "batch_density_output",
        "batch_2Ddensity_output",
        "batch_eccentricity_output",
        "outputInitialEst",
        "minBiasTableFile"};
    std::ostringstream header;
    header << "# collision_energy_GeV="
           << std::setprecision(std::numeric_limits<real>::max_digits10)
           << collisionEnergy << " parameters";
    for (const auto& [name, value] : parameters.get_parameters()) {
        if (excluded.count(name) == 0) {
            header << ' ' << name << '=' << std::quoted(value);
        }
    }
    return header.str();
}

bool readMinBiasTable(
    std::istream& input, const std::string& parameters,
    std::vector<float>& estimators, unsigned int& seed) {
    std::string version, storedParameters;
    if (!std::getline(input, version) || version != minBiasTableVersion
        || !std::getline(input, storedParameters)
        || storedParameters != parameters || !(input >> seed)) {
        return false;
    }
    for (auto& estimator : estimators) {
        if (!(input >> estimator) || !std::isfinite(estimator) || estimator < 0
            || estimator != std::floor(estimator)) {
            return false;
        }
    }
    input >> std::ws;
    return input.eof() && !input.bad()
           && std::is_sorted(
               estimators.begin(), estimators.end(), std::greater<float>());
}

bool writeMinBiasTable(
    const std::string& filename, const std::string& parameters,
    const std::vector<float>& estimators, unsigned int seed) {
    // Publish the complete table in one rename, even when jobs share a path.
    std::string temporary = filename + ".XXXXXX";
    int fd = mkstemp(temporary.data());
    if (fd == -1) return false;
    close(fd);

    std::ofstream output(temporary);
    output << minBiasTableVersion << '\n'
           << parameters << '\n'
           << seed << '\n'
           << std::setprecision(std::numeric_limits<float>::max_digits10);
    for (const auto estimator : estimators) {
        output << estimator << '\n';
    }
    output.close();
    if (output && std::rename(temporary.c_str(), filename.c_str()) == 0) {
        return true;
    }
    std::remove(temporary.c_str());
    return false;
}

}  // namespace

std::vector<CollisionEvent> EventGenerator::get_CollisionEventvector() {
    return mc_glauber_ptr_->get_collision_information();
}

EventGenerator::EventGenerator(
    std::string input_filename, int argc, char* argv[], int seed) {
    parameter_list_.read_in_parameters_from_file(input_filename);
    parameter_list_.read_in_parameters_from_arguments(argc, argv, "=", 4);
    parameter_list_.print_parameter_list();
    int ran_seed = parameter_list_.get_seed();
    if (seed != 0) ran_seed = seed;
    auto gamma_beta = parameter_list_.get_tau_form_fluct_gamma_beta();
    ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(
        new RandomUtil::Random(ran_seed, 0.0, 1.0, gamma_beta));
    mc_glauber_ptr_ =
        std::unique_ptr<Glauber>(new Glauber(parameter_list_, ran_gen_ptr_));

    statistics_only_ = parameter_list_.get_only_event_statistics();
    batchDensityOutput_ = parameter_list_.get_batch_density_output();
    batchEccOutput_ = parameter_list_.get_batch_eccentricity_output();
    initialEstOutput_ = parameter_list_.get_initialEst_output();
    density_maker_ptr_ = std::unique_ptr<MakeDensity>(new MakeDensity());
    // set grid information
    density_maker_ptr_->set_1D_grid_info_eta(72, 0.2);
    density_maker_ptr_->set_2D_grid_info_etax(72, 0.2, 100, 0.2);
    density_maker_ptr_->set_3D_grid_info(100, 0.2, 100, 0.2, 72, 0.2);

    // set Gaussian widths for sigma_x, sigma_eta
    density_maker_ptr_->setGaussianWidths(0.4, 0.5);
    density_maker_ptr_->setStringTransShiftFrac(0.0);
    if (batchDensityOutput_) {
        batchDensity2DOutput_ = parameter_list_.get_batch_2Ddensity_output();
    } else {
        batchDensity2DOutput_ = false;
    }
    cenEstMax_ = 1e16;
    cenEstMin_ = 0;
    NpartMin_ = parameter_list_.getParam("NpartMin", 1);
    NpartMax_ = parameter_list_.getParam("NpartMax", 1000);
}

float EventGenerator::computeCenEstimator(
    const int Npart, const int Ncoll, const int Nstrings) const {
    return (static_cast<float>(Nstrings));
}

void EventGenerator::generateMinBiasEventList() {
    real cenMin = parameter_list_.getParam("cenMin", 0.);
    real cenMax = parameter_list_.getParam("cenMax", 100.);
    if (cenMin < 1e-8 && 100 - cenMax < 1e-8) return;

    const int nev = 50000;
    const auto filename =
        parameter_list_.checkParamIsDefined("minBiasTableFile")
            ? parameter_list_.get_param_val("minBiasTableFile")
            : "minimum_bias_collisions.dat";
    const auto parameters = minBiasTableParameters(
        parameter_list_, mc_glauber_ptr_->get_collision_energy());
    const auto bMin = parameter_list_.get_b_min();
    const auto bMax = parameter_list_.get_b_max();
    std::vector<float> estimators(nev);
    unsigned int tableSeed = 0;
    std::ifstream input(filename);
    if (readMinBiasTable(input, parameters, estimators, tableSeed)) {
        messager << "Loaded minimum-bias table from " << filename
                 << " (table seed = " << tableSeed << ")";
        messager.flush("info");
    } else {
        if (input.is_open()) {
            messager << "Invalid or incompatible minimum-bias table in "
                     << filename << "; regenerating";
            messager.flush("warning");
            input.close();
        }
        estimators = sampleMinBiasEventList(nev);
        if (writeMinBiasTable(
                filename, parameters, estimators, ran_gen_ptr_->get_seed())) {
            messager << "Saved minimum-bias table to " << filename;
            messager.flush("info");
        } else {
            messager << "Could not save minimum-bias table to " << filename;
            messager.flush("warning");
        }
    }

    // Apply the same parameter normalization on cache hits as after sampling.
    parameter_list_.set_b_min(bMin);
    parameter_list_.set_b_max(bMax);
    int idx = std::min(nev - 1, static_cast<int>(nev * cenMax / 100.));
    cenEstMin_ = estimators[idx];
    idx = std::max(0, static_cast<int>(nev * cenMin / 100.));
    cenEstMax_ = estimators[idx];
    messager << "centrality cut [" << cenMin << ", " << cenMax
             << "]: cenEstMin = " << cenEstMin_
             << ", cenEstMax = " << cenEstMax_;
    messager.flush("info");
}

std::vector<float> EventGenerator::sampleMinBiasEventList(int nev) {
    messager << "Random seed = " << ran_gen_ptr_->get_seed();
    messager.flush("info");
    messager << "Generating minimum bias with " << nev << " events ... ";
    messager.flush("info");
    parameter_list_.set_b_min(0.);
    parameter_list_.set_b_max(25.);
    auto b_max_local = 0.;
    std::vector<float> estimators;
    estimators.reserve(nev);

    int iev = 0;
    int nev_progress = std::max(1, nev / 10);
    while (iev < nev) {
        mc_glauber_ptr_->make_nuclei();
        auto Ncoll = mc_glauber_ptr_->make_collision_schedule();
        auto Npart = mc_glauber_ptr_->get_Npart();

        if (Npart < 2) continue;

        auto b = mc_glauber_ptr_->get_impact_parameter();
        if (b > b_max_local) b_max_local = b;

        if (iev % nev_progress == 0) {
            messager << "Progress: " << iev << " out of " << nev << " is done.";
            messager.flush("info");
            if (iev > 0) {
                parameter_list_.set_b_max(b_max_local + 1.);
            }
        }

        auto Nstrings = mc_glauber_ptr_->decide_QCD_strings_production();
        estimators.push_back(computeCenEstimator(Npart, Ncoll, Nstrings));
        iev++;
    }
    std::sort(estimators.begin(), estimators.end(), std::greater<float>());

    return estimators;
}

void EventGenerator::generate_events(int nev, int event_id_offset) {
    messager << "Random seed = " << ran_gen_ptr_->get_seed();
    messager.flush("info");
    messager << "Generating " << nev << " events ... ";
    messager.flush("info");
    // this file records all the essential information for the generated events
    std::ofstream record_file("events_summary.dat", std::ios::out);
    if (statistics_only_) {
        record_file << "# event_id  Npart  Ncoll  Nstrings  b(fm)  "
                    << "PartProtons  PartNeutrons" << std::endl;
    } else {
        record_file << "# event_id  Npart  Ncoll  Nstrings  b(fm)  "
                    << "PartProtons  PartNeutrons  trans_area(fm^2)"
                    << std::endl;
    }

    int iev = 0;
    long long int icollisions = 0;
    int nev_progress = std::max(1, nev / 10);
    int mean_Npart = 0;
    while (iev < nev) {
        mc_glauber_ptr_->make_nuclei();
        auto Ncoll = mc_glauber_ptr_->make_collision_schedule();
        auto Npart = mc_glauber_ptr_->get_Npart();
        auto NpartNeut = mc_glauber_ptr_->get_Neutronspart();
        auto NpartProt = mc_glauber_ptr_->get_Protonspart();
        auto Nstrings = mc_glauber_ptr_->decide_QCD_strings_production();
        icollisions++;
        if (event_of_interest_trigger(Npart, Ncoll, Nstrings)) {
            int event_id = iev + event_id_offset;
            mean_Npart += Npart;
            if (iev % nev_progress == 0) {
                messager << "Progress: " << iev << " out of " << nev
                         << " is done.";
                messager.flush("info");
            }
            if (!batchDensityOutput_) {
                std::ostringstream binaryCollFilename;
                binaryCollFilename << "binaryCollisions_event_" << event_id
                                   << ".dat";
                mc_glauber_ptr_->outputBinaryCollisions(
                    binaryCollFilename.str());
            }

            Ncoll = mc_glauber_ptr_->perform_string_production();
            auto b = mc_glauber_ptr_->get_impact_parameter();

            // write event information to the record file
            record_file << event_id << "  " << Npart << "  " << Ncoll << "  "
                        << Nstrings << "  " << b << "  " << NpartProt << "  "
                        << NpartNeut;
            iev++;
            if (statistics_only_) {
                record_file << std::endl;
                continue;
            }

            density_maker_ptr_->set_QCD_string_output_arr(
                mc_glauber_ptr_->get_QCD_strings_output_list());

            double s_perp = density_maker_ptr_->compute_transverse_area();
            record_file << "  " << s_perp << std::endl;

            if (initialEstOutput_) {
                density_maker_ptr_->output_eccentricity(
                    "ecc_ed_n", event_id, 1);
                density_maker_ptr_->output_netBaryon_eta_distribution(
                    "nB_etas_distribution", event_id, 1);
                density_maker_ptr_->output_netElectricCharges_eta_distribution(
                    "nQ_etas_distribution", event_id, 1);
                density_maker_ptr_->output_energyDensity_eta_distribution(
                    "ed_etas_distribution", event_id, 1);
            }

            if (batchDensityOutput_) {
                density_maker_ptr_->output_netBaryon_eta_distribution(
                    "nB_etas_distribution", event_id);
                density_maker_ptr_->output_netElectricCharges_eta_distribution(
                    "nQ_etas_distribution", event_id);
                density_maker_ptr_->output_energyDensity_eta_distribution(
                    "ed_etas_distribution", event_id);
                if (batchDensity2DOutput_) {
                    density_maker_ptr_->output_energyDensity_xeta_distribution(
                        "ed2D_xetas_distribution", event_id);
                    density_maker_ptr_->output_energyDensity_3d(
                        "ed3D", event_id);
                }
                if (batchEccOutput_) {
                    density_maker_ptr_->output_eccentricity(
                        "ecc_ed_n", event_id);
                    density_maker_ptr_->setParticipantList(
                        mc_glauber_ptr_->getParticipantList());
                    density_maker_ptr_->outputTATBEccentricity(
                        "ecc_ed_n", event_id);
                }
            } else {
                std::ostringstream filename;
                filename << "strings_event_" << event_id << ".dat";
                mc_glauber_ptr_->output_QCD_strings(
                    filename.str(), Npart, Ncoll, Nstrings, b,
                    ran_gen_ptr_->get_seed());
                std::ostringstream specFileName;
                specFileName << "spectators_event_" << event_id << ".dat";
                mc_glauber_ptr_->output_spectators(specFileName.str());

                std::ostringstream partFileName;
                partFileName << "participants_event_" << event_id << ".dat";
                mc_glauber_ptr_->outputParticipants(partFileName.str());
            }
        }
    }
    record_file.close();
    mean_Npart = static_cast<real>(mean_Npart) / static_cast<real>(nev);
    messager << "Completed. <Npart> = " << mean_Npart;
    messager.flush("info");

    if (nev > 1000) {
        auto b_max = parameter_list_.get_b_max();
        auto b_min = parameter_list_.get_b_min();
        auto total_cross_section =
            (M_PI * (b_max * b_max - b_min * b_min) * static_cast<real>(nev)
             / static_cast<real>(icollisions) / 100.);
        messager << "Total cross section sig_tot = " << total_cross_section
                 << " b";
        messager.flush("info");
    }
}

bool EventGenerator::event_of_interest_trigger(
    const int Npart, const int Ncoll, const int Nstrings) const {
    bool pick = false;
    float cenEst = computeCenEstimator(Npart, Ncoll, Nstrings);
    pick =
        ((Npart >= NpartMin_) && (Npart <= NpartMax_) && (cenEst >= cenEstMin_)
         && (cenEst <= cenEstMax_));
    return (pick);
}

};  // namespace MCGlb
