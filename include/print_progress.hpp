
#pragma once

#include <cstdio>
#include <chrono>

template <typename T>
inline void print_progress(T percentage) {
	static const auto PBSTR = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
	static const auto PBWIDTH = 60;

	static std::chrono::steady_clock::time_point start;

	if (percentage==.0f)
		start = std::chrono::steady_clock::now();

	// Format duration
	std::string suffix = "";
	const auto duration = std::chrono::steady_clock::now() - start;
	const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
	const auto minutes = static_cast<int>(seconds)/60;
	const auto hours   = minutes/60;
	suffix +=         std::to_string(static_cast<int>(hours)) + "h";
	suffix += " : " + std::to_string(static_cast<int>(minutes%60)) + "m";
	suffix += " : " + std::to_string(static_cast<int>(seconds%60)) + "s";

	// Remaining time
	if (percentage < 1.f && percentage > 0.01f && seconds>3) {
		const auto remaining = static_cast<int>(seconds / percentage - seconds);
		const auto minutes = static_cast<int>(remaining)/60;
		const auto hours   = minutes/60;
		suffix += "        Remaining: " + std::to_string(static_cast<int>(hours)) + "h";
		suffix += " : " + std::to_string(static_cast<int>(minutes%60)) + "m";
		suffix += " : " + std::to_string(static_cast<int>(remaining%60)) + "s";
	}

	suffix += "  ";

	const auto val  = static_cast<int>(percentage * T(100));
	const auto lpad = static_cast<int>(percentage * PBWIDTH);
	const auto rpad = PBWIDTH - lpad;
	std::printf("\r%3d%% [%.*s%*s]      %s", val, lpad, PBSTR, rpad, "", suffix.c_str());
	std::fflush(stdout);
}
