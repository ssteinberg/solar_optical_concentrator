
#pragma once

#include <cstddef>
#include <random>
#include <chrono>
#include <thread>
#include <functional>
#include <climits>

// Thread-safe random generator with strong initial seed
class random_generator {
public:
	struct engine {
		// Generates a seeded random engine, usually a Mersenne Twister engine.
		template<class T = std::mt19937, std::size_t N = T::state_size * T::word_size>
		static T generate_seeded_random_engine() noexcept {
			std::array<unsigned,N> data;
			std::random_device source("rdrand");
			std::generate(std::begin(data), std::end(data), std::ref(source));

			const auto tse = std::chrono::high_resolution_clock::now().time_since_epoch();
			const auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(tse).count();
			const auto tid = std::hash<std::thread::id>{}(std::this_thread::get_id());

			reinterpret_cast<std::uint64_t*>(data.data())[0] ^= std::hash<std::uint64_t>{}(static_cast<std::uint64_t>(time));
			reinterpret_cast<std::uint64_t*>(data.data())[1] ^= std::hash<std::uint64_t>{}(static_cast<std::uint64_t>(tid));

			std::seed_seq seq(std::begin(data), std::end(data));

			// Seed the engine
			return T{ seq };
		}
		std::mt19937 random_engine;
		engine() noexcept : random_engine(generate_seeded_random_engine()) {}
	};

private:
	static thread_local engine generator;

public:
	auto& operator()() noexcept { return generator.random_engine; }
};
