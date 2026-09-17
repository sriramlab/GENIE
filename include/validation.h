#ifndef VALIDATION_H
#define VALIDATION_H

// UPG-007: Centralized input validation for GENIE
// This module provides validation functions to catch common user errors early
// with clear error messages.
//
// Note: This header must be included AFTER the options struct is defined
// (see arguments.h for include location).

#include <iostream>
#include <cstring>

// command_line_opts is extern declared in arguments.h, defined in ge_flexible.cpp

inline void validate_memory_mode_model_compatibility() {
    // UPG-001: Memory-efficient modes only work with Model G
    // GxE models crash (SIGSEGV) with memeff enabled
    const char* model_arg = command_line_opts.model.c_str();
    if (command_line_opts.memeff && strcmp(model_arg, "G") != 0) {
        std::cerr << "ERROR: Memory-efficient mode (memeff=1) is only supported for Model G." << std::endl;
        std::cerr << "  Model requested: " << command_line_opts.model << std::endl;
        std::cerr << "  Hint: Use model=G, or disable memory-efficient mode (memeff=0 or memory_mode=0)." << std::endl;
        exit(1);
    }
}

inline void validate_parameter_ranges() {
    // Warn on unusual num_vec values
    if (command_line_opts.num_of_vec < 5) {
        std::cerr << "Warning: num_vec=" << command_line_opts.num_of_vec
                  << " is low. Recommended: 10 for stable estimates." << std::endl;
    } else if (command_line_opts.num_of_vec > 100) {
        std::cerr << "Warning: num_vec=" << command_line_opts.num_of_vec
                  << " is unusually high. Recommended: 10-20." << std::endl;
    }

    // Warn on unusual num_jack values
    if (command_line_opts.jack_number < 10) {
        std::cerr << "Warning: num_jack=" << command_line_opts.jack_number
                  << " is low. Recommended: 100 for stable standard errors." << std::endl;
    } else if (command_line_opts.jack_number > 1000) {
        std::cerr << "Warning: num_jack=" << command_line_opts.jack_number
                  << " is unusually high. Recommended: 100." << std::endl;
    }

    // Warn on unusual thread count
    if (command_line_opts.nthreads > 64) {
        std::cerr << "Warning: nthreads=" << command_line_opts.nthreads
                  << " is very high. Consider reducing if performance degrades." << std::endl;
    }
}

inline void validate_params() {
    // Run all validation checks
    validate_memory_mode_model_compatibility();
    validate_parameter_ranges();
}

#endif // VALIDATION_H
