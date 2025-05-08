#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <stdbool.h>
#include <time.h>

// Log levels
typedef enum {
    LOG_DEBUG,
    LOG_INFO,
    LOG_SUCCESS,
    LOG_WARNING,
    LOG_ERROR,
    LOG_FATAL
} LogLevel;

// ANSI color codes (same as your defines)
#define LOG_RESET   "\033[0m"
#define LOG_BLUE    "\033[34m"
#define LOG_GREEN   "\033[32m"
#define LOG_YELLOW  "\033[33m"
#define LOG_RED     "\033[31m"

// Initialize logger
void log_init(FILE* terminal_output, const char* log_file_path, bool use_colors);

// Close logger (call before program exit)
void log_close();

// Main logging function
void log_message(LogLevel level, const char* file, int line, const char* format, ...);

// Log file content
void log_file_content(const char* file_path, const char* description);

// Convenience macros with your style
#define LOG_DEBUG(...)    log_message(LOG_DEBUG,   __FILE__, __LINE__, __VA_ARGS__)
#define LOG_INFO(...)     log_message(LOG_INFO,    __FILE__, __LINE__, __VA_ARGS__)
#define LOG_SUCCESS(...)  log_message(LOG_SUCCESS, __FILE__, __LINE__, __VA_ARGS__)
#define LOG_WARN(...)     log_message(LOG_WARNING, __FILE__, __LINE__, __VA_ARGS__)
#define LOG_ERROR(...)    log_message(LOG_ERROR,   __FILE__, __LINE__, __VA_ARGS__)
#define LOG_FATAL(...)    log_message(LOG_FATAL,   __FILE__, __LINE__, __VA_ARGS__)

// Special line message (similar to your LINEMSG)
#define LOG_LINE()        log_message(LOG_DEBUG, __FILE__, __LINE__, "Line execution")

#endif // LOGGER_H