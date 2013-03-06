
#include <iostream>
#include <unistd.h>
#include "timer.hpp"

static std::clock_t get_current_cpu_time();
static void get_elapsed_cpu_time(time_data &elapsed, std::clock_t &start_time,
				 const bool reset_start = false);

static void get_current_wall_time(time_data &current);
static void get_elapsed_wall_time(time_data &elapsed, time_data &start,
				  const bool reset_start = false);

inline std::clock_t get_current_cpu_time()
{
  struct tms data;
  (void) times(&data);
  return data.tms_utime + data.tms_stime;
}

inline void get_current_wall_time(time_data &current)
{
  timeval t;
  gettimeofday(&t, 0);
  current.seconds = t.tv_sec;
  current.microsecs = t.tv_usec;
}

void get_elapsed_cpu_time(time_data &elapsed,
			  std::clock_t &start_time,
			  const bool reset_start)
{
  static long ticks = 0;
  if (ticks == 0)
    ticks = sysconf(_SC_CLK_TCK);
  std::clock_t current = get_current_cpu_time();
  std::clock_t diff = current - start_time;
  if (diff < 0) {
    // clock wrapped around - Todo, needs to see sizeof(clock_t)
    diff += 0x7fffffff;
  }
  elapsed.seconds = diff / ticks;
  elapsed.microsecs = diff - (elapsed.seconds * ticks);
  elapsed.microsecs *= (1000000 / ticks);
  if (reset_start)
    start_time = current;
}

void get_elapsed_wall_time(time_data &elapsed, time_data &start,
			   const bool reset_start)
{
  timeval t;
  gettimeofday(&t, 0);
  elapsed.seconds = t.tv_sec - start.seconds;
  elapsed.microsecs = t.tv_usec - start.microsecs;
  if (elapsed.microsecs > 1000000) {
    elapsed.microsecs -= 1000000;
    elapsed.seconds++;
  } else if (elapsed.microsecs < 0) {
    elapsed.microsecs += 1000000;
    elapsed.seconds--;
  }
  if (reset_start) {
    start.seconds = t.tv_sec;
    start.microsecs = t.tv_usec;
  }
}

timer::timer(const bool use_timer)
{
  if (use_timer) {
    start_cpu = get_current_cpu_time();
    get_current_wall_time(start_wall);
  }
  cpu.seconds = 0;
  cpu.microsecs = 0;
  wall.seconds = 0;
  wall.microsecs = 0;
  use = use_timer;
}

void timer::reset()
{
  if (use) {
    start_cpu = get_current_cpu_time();
    get_current_wall_time(start_wall);
    cpu.seconds = 0;
    cpu.microsecs = 0;
    wall.seconds = 0;
    wall.microsecs = 0;
  }
}

void timer::mark()
{
  if (use) {
    start_cpu = get_current_cpu_time();
    get_current_wall_time(start_wall);
  }
}

void timer::accumulate()
{
  if (use) {
    time_data temp_cpu;
    time_data temp_wall;
    get_elapsed_cpu_time(temp_cpu, start_cpu, true);
    get_elapsed_wall_time(temp_wall, start_wall, true);
    cpu.seconds += temp_cpu.seconds;
    cpu.microsecs += temp_cpu.microsecs;
    if (cpu.microsecs > 999999) {
      cpu.seconds++;
      cpu.microsecs -= 1000000;
    }
    wall.seconds += temp_wall.seconds;
    wall.microsecs += temp_wall.microsecs;
    if (wall.microsecs > 999999) {
      wall.seconds++;
      wall.microsecs -= 1000000;
    }
  }
}

void timer::output(const char message[])
{
  if (use) {
    std::cout << message << ": ";
    std::cout << cpu.seconds << '.';
    std::cout.width(6);
    std::cout.fill('0');
    std::cout << cpu.microsecs << " cpu time, ";
    std::cout << wall.seconds << '.';
    std::cout.width(6);
    std::cout.fill('0');
    std::cout << wall.microsecs << " wall time\n";
  }
}
