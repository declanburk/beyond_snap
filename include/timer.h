#include <chrono>

struct Timer {
    std::chrono::time_point<std::chrono::system_clock> start_;
    std::chrono::time_point<std::chrono::system_clock> end_;

    Timer() {
        start_ = std::chrono::system_clock::now();
    }

    void start() {
        start_ = std::chrono::system_clock::now();
    }

    double stop() {
        end_ = std::chrono::system_clock::now();
        return getTime();
    }

    double getTime() const {
        if (end_ < start_) {
            return 0.0;
        }
        std::chrono::duration<double> duration = end_ - start_;
        return duration.count();
    }
};

struct LapTimer : Timer {
    int lps = 0;
    double avg = 0.0;

    LapTimer() : Timer() {}

    double stop() {
        end_ = std::chrono::system_clock::now();

        double dur = getTime();

        addTime(dur);

        return dur;
    }

    void addTime(double dur) {
        if (lps == 0) {
            avg = dur;
        } else {
            avg = (avg * (lps / (lps + 1.0))) + (dur / (lps + 1.0));
        }
        lps++;
    }
};