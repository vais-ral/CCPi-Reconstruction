#ifndef CCPI_USERINTERFACE_QT_H
#define CCPI_USERINTERFACE_QT_H
#include <QObject>
#include "base_types.hpp"
class UserInterfaceQT :
	public QObject
{
	Q_OBJECT

public:
	UserInterfaceQT(void);
	~UserInterfaceQT(void);
	void set_initialise_progress(const int length, const char label[]);
	void set_update_progress(const int value);
	void set_report_error(const std::string message);
	void set_report_error(const std::string message, const std::string arg);
	void set_report_error(const std::string message, const std::string arg1,
			const std::string arg2);
	void set_add_output(const std::string str);
	void set_add_output(const char c);
	void set_add_output(const int i);
	void set_add_output(const int i, const int w, const bool fill);
	void set_add_output(const sl_int i);
	void set_add_output(const real r);
	void set_send_output(const QString &str);

signals:
	void initialise_progress(const int length, const QString& label);
	void update_progress(const int value);
	void report_error(const std::string message);
	void report_error(const std::string message, const std::string arg);
	void report_error(const std::string message, const std::string arg1,
			const std::string arg2);
	void add_output(const QString &str);
	void send_output(const QString &str);
};
#endif
