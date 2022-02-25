
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <optional>

inline std::ostream& operator<<(std::ostream& os, std::chrono::nanoseconds ns) {
    using days = std::chrono::duration<int, std::ratio<86400>>;
    
    const auto& fill = os.fill();
    os.fill('0');

    const auto& d = std::chrono::duration_cast<days>(ns);
    ns -= d;
    const auto& h = std::chrono::duration_cast<std::chrono::hours>(ns);
    ns -= h;
    const auto& m = std::chrono::duration_cast<std::chrono::minutes>(ns);
    ns -= m;
    const auto& s = std::chrono::duration_cast<std::chrono::seconds>(ns);
    ns -= s;

    if (d.count())
        os << std::setw(2) << d.count() << " days ";
    if (h.count())
       os << std::setw(2) << h.count() << "h:";
    os << std::setw(2) << m.count() << "m:"
       << std::setw(2) << s.count() << 's';

    std::optional<int> fs_count;
    switch (os.precision()) {
    case 9: fs_count = ns.count();
        break;
    case 6: fs_count = std::chrono::duration_cast<std::chrono::microseconds>(ns).count();
        break;
    case 3: fs_count = std::chrono::duration_cast<std::chrono::milliseconds>(ns).count();
        break;
    }
    if (fs_count)
        os << "." << std::setw(os.precision()) << fs_count.value();

    os.fill(fill);

    return os;
};
