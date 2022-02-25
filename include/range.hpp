
#pragma once

#include <string>
#include <type_traits>
#include <iostream>

template <typename T>
struct range {
    struct range_iterator {
        const range<T> *rng;
        std::size_t idx;

        auto  operator*() const noexcept { return rng->s + rng->step*T(idx); }
        auto  operator+(std::size_t n) noexcept { return range_iterator{ rng,idx+n }; }
        auto& operator++() noexcept { ++idx; return *this; }
        auto  operator++(int) noexcept {
            const auto ret = *this;
            ++(*this);
            return ret;
        }
        bool operator==(const range_iterator &it) const noexcept { return rng==it.rng && idx==it.idx; }
        bool operator!=(const range_iterator &it) const noexcept { return !((*this)==it); }
    };
    
    static_assert(std::is_floating_point_v<T>);
    T s,e,step{};
    std::size_t n;
    
    range(T val) noexcept : s(val),e(val),step(0),n(1) {}
    range(T start, T end, std::size_t steps) noexcept : s(start),e(end),n(steps) {
        step = (e-s)/T(n);
    }
    range(const std::string &str) noexcept {
        s = (T)std::stod(str);
        const auto endidx = str.find('-');
        const auto stepidx = str.find(',');
        if (endidx==std::string::npos ||
            stepidx==std::string::npos ||
            stepidx<endidx || stepidx+1>=str.length()) {
            e=s;
            n=1;
            step=0;
            return;
        }
        e = (T)std::stod(str.substr(endidx+1));
        step = (T)std::stod(str.substr(stepidx+1));
        n = (std::size_t)((e-s)/step)+1;
    }

    range(range&&) noexcept = default;
    range(const range&) noexcept = default;
    range& operator=(range&&) noexcept = default;
    range& operator=(const range&) noexcept = default;

    auto begin() const { return range_iterator{ this,0 }; }
    auto end() const { return range_iterator{ this,n }; }
    auto operator[](std::size_t idx) const noexcept { return *(begin()+idx); }
    auto size() const noexcept { return n; }
};
