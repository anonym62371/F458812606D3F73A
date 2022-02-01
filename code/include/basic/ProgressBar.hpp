#ifndef PROGRESSBAR_PROGRESSBAR_HPP
#define PROGRESSBAR_PROGRESSBAR_HPP

#include <chrono>
#include <iostream>
#include <ostream>

void output_time(std::ostream& os, float time_in_sec) {
    if (time_in_sec < 60.0) {
        os << time_in_sec << "s";
    } else if (time_in_sec < 3600.0) {
        float remaining_time_in_sec = time_in_sec;

        int time_in_min = 0;
        while (remaining_time_in_sec >= 60.0) {
            remaining_time_in_sec -= 60;
            time_in_min++;
        }

        os << time_in_min << "m";

        if (remaining_time_in_sec > 0.0) {
            os << " " << remaining_time_in_sec << "s";
        }
    } else {
        float remaining_time_in_sec = time_in_sec;

        int time_in_hr = 0;
        while (remaining_time_in_sec >= 3600) {
            remaining_time_in_sec -= 3600;
            time_in_hr++;
        }

        os << time_in_hr << "h";

        int time_in_min = 0;
        while (remaining_time_in_sec >= 60) {
            remaining_time_in_sec -= 60;
            time_in_min++;
        }

        if (time_in_min > 0) {
            os << " " << time_in_min << "m";
        }

        if (remaining_time_in_sec > 0.0) {
            os << " " << remaining_time_in_sec << "s";
        }
    }
}

class ProgressBar {
private:
    unsigned int ticks = 0;

    const unsigned int total_ticks;
    const unsigned int bar_width;
    const char complete_char = '=';
    const char incomplete_char = ' ';
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

public:
    ProgressBar(unsigned int total, unsigned int width, char complete, char incomplete) :
            total_ticks {total}, bar_width {width}, complete_char {complete}, incomplete_char {incomplete} {}

    ProgressBar(unsigned int total, unsigned int width) : total_ticks {total}, bar_width {width} {}

    unsigned int operator++() { return ++ticks; }

    void display() const
    {
        float progress = (float) ticks / total_ticks;
        int pos = (int) (bar_width * progress);

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now-start_time).count();

        std::cout << "[";

        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << complete_char;
            else if (i == pos) std::cout << ">";
            else std::cout << incomplete_char;
        }
        std::cout << "] " << int(progress * 100.0) << "% ";
        output_time(std::cout, float(time_elapsed) / 1000.0);
        std::cout << "\r";
                  // << float(time_elapsed) / 1000.0 << "s\r";
        std::cout.flush();
    }

    void done() const
    {
        display();
        std::cout << std::endl;
    }
};

#endif //PROGRESSBAR_PROGRESSBAR_HPP
