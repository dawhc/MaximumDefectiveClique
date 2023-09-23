#ifndef LOG_HPP
#define LOG_HPP

#pragma once

#include <cstring>
#include <cstdio>
#include <chrono>
#include <string>
#include <iomanip>

namespace logging {
	static std::string formatCurrentTime() {
		auto timepoint = std::chrono::system_clock::now();
		auto seconds = std::chrono::duration_cast<std::chrono::seconds>(timepoint.time_since_epoch());
		auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(timepoint.time_since_epoch() - seconds);
		auto sec = static_cast<time_t>(seconds.count());
		std::tm tm;
		gmtime_r(&sec, &tm);
		char buf[32] = {0};
		auto size = std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm);
		return std::string(buf, size) + "." + std::to_string(milliseconds.count());
	}

};

#define log(fmt, args...) do {														\
	fprintf(stderr, "\033[2;33m[%s]\033[0m[\033[34m%s\033[0m:\033[32m%d\033[0m in \033[35m%s\033[0m] ",											\
		logging::formatCurrentTime().c_str(), 										\
		strrchr(__FILE__, '/')+1, __LINE__, __FUNCTION__);							\
	fprintf(stderr, fmt, ##args);													\
	fprintf(stderr, "\n");															\
} while (0);																	    \

#endif