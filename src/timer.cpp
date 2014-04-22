
#ifndef WIN32
#  include <unistd.h>
#endif // WIN32
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // types
#include "timer.hpp"
#include "ui_calls.hpp"

static std::clock_t get_current_cpu_time();
static void get_elapsed_cpu_time(time_data &elapsed, std::clock_t &start_time,
				 const bool reset_start = false);

static void get_current_wall_time(time_data &current);
static void get_elapsed_wall_time(time_data &elapsed, time_data &start,
				  const bool reset_start = false);

#ifdef WIN32

#  ifdef MATLAB_MEX_FILE

inline std::clock_t get_current_cpu_time()
{
	return 0;
}

inline void get_current_wall_time(time_data &current)
{
	current.seconds = 0;
	current.microsecs = 0;
}

#  else

#include <Windows.h>
#include <time.h>

inline std::clock_t get_current_cpu_time()
{
  FILETIME a, b, kernel, user;
  GetProcessTimes(GetCurrentProcess(), &a, &b, &kernel, &user);
  __int64 t1 = (((__int64)kernel.dwHighDateTime) << 32) + kernel.dwLowDateTime;
  __int64 t2 = (((__int64)user.dwHighDateTime) << 32) + user.dwLowDateTime;
  // 100 ns units
  __int64 t = (t1 + t2) / (10000000LL / CLOCKS_PER_SEC);
  return t;
}

inline void get_current_wall_time(time_data &current)
{
  long ticks = CLOCKS_PER_SEC;
  long timer = clock();
  current.seconds = timer / ticks;
  current.microsecs = (timer % ticks) * (1000000 / ticks);
}

#  endif // MEX_FILE

#else

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

#endif // WIN32

void get_elapsed_cpu_time(time_data &elapsed,
			  std::clock_t &start_time,
			  const bool reset_start)
{
#ifdef WIN32
  static long ticks = CLOCKS_PER_SEC;
#else
  static long ticks = 0;
  if (ticks == 0)
    ticks = sysconf(_SC_CLK_TCK);
#endif // WIN32
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
  time_data cur;
  get_current_wall_time(cur);
  elapsed.seconds = cur.seconds - start.seconds;
  elapsed.microsecs = cur.microsecs - start.microsecs;
  if (elapsed.microsecs > 1000000) {
    elapsed.microsecs -= 1000000;
    elapsed.seconds++;
  } else if (elapsed.microsecs < 0) {
    elapsed.microsecs += 1000000;
    elapsed.seconds--;
  }
  if (reset_start) {
    start.seconds = cur.seconds;
    start.microsecs = cur.microsecs;
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
    std::string m = message;
    add_output(m);
    add_output(": ");
    add_output((int)cpu.seconds);
    add_output('.');
    add_output(cpu.microsecs, 6, true);
    add_output(" cpu time, ");
    add_output((int)wall.seconds);
    add_output('.');
    add_output(wall.microsecs, 6, true);
    add_output(" wall time");
    send_output();
  }
}
