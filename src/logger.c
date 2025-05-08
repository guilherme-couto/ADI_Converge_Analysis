#include "../include/logger.h"
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

static FILE* terminal_out = NULL;
static FILE* file_out = NULL;
static bool colors_enabled = false;

void log_init(FILE* terminal_output, const char* log_file_path, bool use_colors) {
    terminal_out = terminal_output ? terminal_output : stderr;
    colors_enabled = use_colors;
    
    if (log_file_path) {
        file_out = fopen(log_file_path, "a");
        if (!file_out) {
            fprintf(terminal_out, "Failed to open log file: %s\n", log_file_path);
            file_out = NULL;
        }
    }
}

void log_close() {
    if (file_out) {
        fclose(file_out);
        file_out = NULL;
    }
}

static const char* get_level_string(LogLevel level) {
    switch (level) {
        case LOG_DEBUG:   return "DEBUG";
        case LOG_INFO:    return "INFO";
        case LOG_SUCCESS: return "SUCCESS";
        case LOG_WARNING: return "WARNING";
        case LOG_ERROR:   return "ERROR";
        case LOG_FATAL:   return "FATAL";
        default:          return "UNKNOWN";
    }
}

static const char* get_level_color(LogLevel level) {
    if (!colors_enabled) return "";
    
    switch (level) {
        case LOG_DEBUG:   return LOG_YELLOW;
        case LOG_INFO:    return LOG_BLUE;
        case LOG_SUCCESS: return LOG_GREEN;
        case LOG_WARNING: return LOG_YELLOW;
        case LOG_ERROR:   return LOG_RED;
        case LOG_FATAL:   return LOG_RED;
        default:          return "";
    }
}

static const char* get_level_prefix(LogLevel level) {
    switch (level) {
        case LOG_DEBUG:   return "[.]";
        case LOG_INFO:    return "[i]";
        case LOG_SUCCESS: return "[+]";
        case LOG_WARNING: return "[!]";
        case LOG_ERROR:   return "[x]";
        case LOG_FATAL:   return "[X]";
        default:          return "[?]";
    }
}

void log_message(LogLevel level, const char* file, int line, const char* format, ...) {
    time_t now;
    time(&now);
    char time_buf[20];
    strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", localtime(&now));
    
    // Format the message
    va_list args;
    va_start(args, format);
    
    // Terminal output (with colors)
    if (terminal_out) {
        fprintf(terminal_out, "%s%s%s %s %s:%d: ",
                get_level_color(level), get_level_prefix(level), 
                colors_enabled ? LOG_RESET : "",
                time_buf, file, line);
        
        vfprintf(terminal_out, format, args);
        fputc('\n', terminal_out);
        fflush(terminal_out);
    }
    
    // File output (no colors)
    if (file_out) {
        fprintf(file_out, "%s %-7s %s:%d: ",
                time_buf, get_level_string(level), file, line);
        
        vfprintf(file_out, format, args);
        fputc('\n', file_out);
        fflush(file_out);
    }
    
    va_end(args);
}

void log_file_content(const char* file_path, const char* description) {
    FILE* config_file = fopen(file_path, "r");
    if (!config_file) {
        LOG_ERROR("Failed to open file for logging: %s", file_path);
        return;
    }

    // Log header
    LOG_INFO("%s (file content below):", description);
    LOG_INFO("===== START OF FILE: %s =====", file_path);

    // Read and log line by line
    char line[1024];
    while (fgets(line, sizeof(line), config_file)) {
        // Remove trailing newline
        line[strcspn(line, "\n")] = '\0';
        
        // Log the line (using INFO level but without metadata)
        if (terminal_out) {
            fprintf(terminal_out, "    %s\n", line);
        }
        if (file_out) {
            fprintf(file_out, "    %s\n", line);
        }
    }

    LOG_INFO("===== END OF FILE: %s =====", file_path);
    fclose(config_file);
}