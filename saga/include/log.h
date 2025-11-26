#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <mutex>
#include <ostream>

/**
* @file
* @author Matthew F. McEneaney
* @date 25/Nov./2025
* @version 0.0.0
* @brief Basic logging utility for SAGA.
*/

namespace saga {

namespace log {

/**
* @class Logger
* @brief Singleton logger class for SAGA.
*/
class Logger {
public:
    enum class Level {
        Debug = 0,
        Info,
        Warning,
        Error
    };

    /**
    * @brief Get the singleton instance of the Logger.
    * @return Reference to the Logger instance.
    */
    static Logger& instance() {
        static Logger inst;
        return inst;
    }

    /**
    * @brief Enable or disable colored output.
    * @param enable True to enable color, false to disable.
    */
    void enableColor(bool enable) {
        std::lock_guard<std::mutex> lock(mtx_);
        colorEnabled_ = enable;
    }

    /**
    * @brief Set the minimum log level.
    * @param level The minimum log level.
    */
    void setLogLevel(Level level) {
        std::lock_guard<std::mutex> lock(mtx_);
        minLevel_ = level;
    }

    /**
    * @brief Set the log level from a string.
    * @param levelStr The log level as a string.
    */
    void setLogLevelFromString(std::string levelStr) {
        std::transform(levelStr.begin(), levelStr.end(),
                       levelStr.begin(), ::tolower);

        if (levelStr == "debuging" ||
                 levelStr == "debug")  setLogLevel(Level::Debug);
        else if (levelStr == "info")   setLogLevel(Level::Info);
        else if (levelStr == "warn" ||
                 levelStr == "warning") setLogLevel(Level::Warning);
        else if (levelStr == "error")  setLogLevel(Level::Error);
    }

    /**
    * @brief Set the output stream.
    * @param os The output stream to use.
    */
    void setOutputStream(std::ostream& os) {
        std::lock_guard<std::mutex> lock(mtx_);
        outStream_ = &os;
    }

    /**
    * @brief Set the log file.
    * @param filename The name of the log file.
    */
    void setLogFile(const std::string& filename) {
        std::lock_guard<std::mutex> lock(mtx_);
        file_.open(filename, std::ios::app);
    }

    /**
    * @brief Log a message.
    * @param level The log level.
    * @param msg The log message.
    */
    void log(Level level, const std::string& msg) {
        std::lock_guard<std::mutex> lock(mtx_);

        // Filter out messages below the current level
        if (level < minLevel_) {
            return;
        }

        std::string out = timestamp() + " [" + levelToString(level) + "] " + msg;

        if (outStream_) {
            // Apply color (if enabled and printing to a terminal)
            if (colorEnabled_ && outStream_ == &std::cout) {  // NEW
                (*outStream_) << colorCode(level) << out << resetCode();
            } else {
                (*outStream_) << out;
            }
            (*outStream_) << std::endl;
        }
        if (file_.is_open()) {
            file_ << out << std::endl;
        }
    }

private:

    /**
    * @brief Constructor for Logger.
    */
    Logger() 
        : outStream_(&std::cout),
          minLevel_(Level::Debug),
          colorEnabled_(true) {}
    /**
    * @brief Destructor for Logger.
    */
    ~Logger() {
        if (file_.is_open()) file_.close();
    }

    std::ofstream file_;
    std::ostream* outStream_;  
    Level minLevel_;
    bool colorEnabled_;
    std::mutex mtx_;

    /**
    * @brief Get the current timestamp.
    * @return The current timestamp as a string.
    */
    std::string timestamp() {
        using namespace std::chrono;
        auto now = system_clock::now();
        std::time_t t = system_clock::to_time_t(now);

        char buffer[32];
        std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S",
                      std::localtime(&t));
        return buffer;
    }

    /**
    * @brief Convert log level to string.
    * @param lvl The log level.
    * @return The log level as a string.
    */
    std::string levelToString(Level lvl) {
        switch (lvl) {
            case Level::Debug:   return "DEBUG";
            case Level::Info:    return "INFO";
            case Level::Warning: return "WARN";
            case Level::Error:   return "ERROR";
        }
        return "UNKNOWN";
    }

    // ---------- Color codes (ANSI) ----------

    /**
    * @brief Get the color code for a log level.
    * @param lvl The log level.
    * @return The color code as a string.
    */
    const char* colorCode(Level lvl) {
        switch (lvl) {
            case Level::Debug:   return "\033[34m";      // Blue
            case Level::Info:    return "\033[32m";      // Green
            case Level::Warning: return "\033[33m";      // Yellow
            case Level::Error:   return "\033[1;31m";    // Bold Red
        }
        return "";
    }

    /**
    * @brief Get the color reset code.
    * @return The color reset code as a string.
    */
    const char* resetCode() {
        return "\033[0m";
    }
};

} // namespace log

} // namespace saga

// Macros
#define LOG_DEBUG(msg) saga::log::Logger::instance().log(saga::log::Logger::Level::Debug,   msg)
#define LOG_INFO(msg)  saga::log::Logger::instance().log(saga::log::Logger::Level::Info,    msg)
#define LOG_WARN(msg)  saga::log::Logger::instance().log(saga::log::Logger::Level::Warning, msg)
#define LOG_ERROR(msg) saga::log::Logger::instance().log(saga::log::Logger::Level::Error,   msg)
