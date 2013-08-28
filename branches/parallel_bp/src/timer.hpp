
#ifndef TIMER
#define TIMER

#ifndef WINDOWS
#  include <sys/time.h>
#  include <sys/times.h>
#endif // WINDOWS
#include <ctime>

struct time_data {
  std::time_t seconds;
  std::time_t microsecs;
};

class timer {
public:
  timer(const bool use_timer);
  void reset();
  void mark();
  void accumulate();
  void output(const char message[]);
    
private:
  std::clock_t start_cpu;
  time_data cpu;
  time_data start_wall;
  time_data wall;
  bool use;
};

#endif // TIMER
