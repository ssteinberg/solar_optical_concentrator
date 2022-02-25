
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once


#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <vector>
#include <queue>
#include <memory>

class thread_pool {
private:
    std::mutex m;
    std::condition_variable cv;
    std::vector<std::thread> threads;
    std::queue<std::function<void(void)>> queue;

    bool alive{ true };

public:
    inline thread_pool(const std::size_t thread_count = std::thread::hardware_concurrency()) {
        threads.reserve(thread_count);
        for (std::size_t t=0; t<thread_count; ++t) {
            threads.emplace_back([this]() {
                for(;this->alive;) {
                    typename decltype(queue)::value_type f;

                    {
                        std::unique_lock<std::mutex> l(m);
                        cv.wait(l, [this]() { return !this->alive || !queue.empty(); });

                        if (!queue.empty()) {
                            f = std::move(queue.front());
                            queue.pop();
                        }
                    }

                    for(;f;) {
                        f();

                        std::unique_lock<std::mutex> l(m);
                        if (queue.empty()) break;

                        f = std::move(queue.front());
                        queue.pop();
                    }
                }
            });
        }
    }
    inline ~thread_pool() noexcept {
        {
            std::unique_lock<std::mutex> l(m);
            alive=false;
            cv.notify_all();
        }
        for (auto &t : threads)
            if (t.joinable()) t.join();
    }

    template <typename Func>
    auto enqueue(Func &&f) noexcept {
        using R = typename std::invoke_result_t<Func>;

        auto pt = std::make_shared<std::packaged_task<R(void)>>(std::move(f));
        auto future = std::move(pt->get_future());

        {
            std::unique_lock<std::mutex> l(m);
            queue.emplace([pt = std::move(pt)]() { (*pt)(); });
            cv.notify_one();
        }

        return future;
    }
};
