
#include <QApplication>
#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>

#include "base_types.hpp"
#include "ui_messages.hpp"

HxMessage *ccpi_recon::messages = 0;
HxWorkArea *ccpi_recon::progress = theWorkArea;
bool ccpi_recon::do_progress = false;

void report_error(const std::string message)
{
  ccpi_recon::messages->error(QApplication::translate("CCPiReconstruction",
						      message.c_str()));
}

void report_error(const std::string message, const std::string arg)
{
  std::string sum = message + arg;
  report_error(sum);
}

void report_error(const std::string message, const std::string arg1,
		  const std::string arg2)
{
  std::string sum = message + arg1;
  sum += arg2;
  report_error(sum);
}

void initialise_progress(const int length, const char label[])
{
  if (ccpi_recon::do_progress) {
    ccpi_recon::progress->stopWorking();
    ccpi_recon::progress->undivide();
  }
  // Turn into busy state, don't activate the Stop button.
  ccpi_recon::progress->subdivide(length);
  ccpi_recon::progress->startWorkingNoStop(QApplication::translate("CCPiReconstruction",
						       label));
  ccpi_recon::do_progress = true;
}

void update_progress(const int value)
{
  //progress->progressStep();
  ccpi_recon::progress->setProgressValue(value);
}

void add_output(const std::string str)
{
  ccpi_recon::messages->printf("%s", str.c_str());
}

void add_output(const char c)
{
  ccpi_recon::messages->printf("%c", c);
}

void add_output(const int i)
{
  ccpi_recon::messages->printf("%1d", i);
}

void add_output(const sl_int i)
{
  ccpi_recon::messages->printf("%1ld", i);
}

void add_output(const int i, const int w, const bool fill)
{
  if (fill)
    ccpi_recon::messages->printf("%0*d", w, i);
  else
    ccpi_recon::messages->printf("%*d", w, i);
}

void add_output(const real r)
{
  ccpi_recon::messages->printf("%f", r);
}

void send_output()
{
  ccpi_recon::messages->printf("\n");
}
