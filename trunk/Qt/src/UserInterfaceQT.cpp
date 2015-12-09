#include "UserInterfaceQT.h"


UserInterfaceQT::UserInterfaceQT(void)
{
}


UserInterfaceQT::~UserInterfaceQT(void)
{
}

void UserInterfaceQT::set_initialise_progress(const int length, const char label[])
{
	emit initialise_progress(length,QString(label));
}

void UserInterfaceQT::set_update_progress(const int value)
{
	emit update_progress(value);
}

void UserInterfaceQT::set_report_error(const std::string message)
{
	emit report_error(message);
}

void UserInterfaceQT::set_report_error(const std::string message, const std::string arg)
{
	emit report_error(message, arg);
}

void UserInterfaceQT::set_report_error(const std::string message, const std::string arg1,
	const std::string arg2)
{
	emit report_error(message, arg1,arg2);
}

void UserInterfaceQT::set_add_output(const std::string str)
{
	emit add_output(QString(str.c_str()));
}

void UserInterfaceQT::set_add_output(const char c)
{
	emit add_output(QString(c));
}

void UserInterfaceQT::set_send_output(const QString& str)
{
	emit send_output(str);
}