#ifndef IO_H

#define IO_H

#include "std.h"
#include "printable.h"

class io {
	public:
		static int debug ;
		static void print (const string& s, int level, string file) ; 
		static void println (const string& s, int level) ;
		static void print (const string& s, int level )  ;
		static void print (printable &p, int level) ;
		static void println(const printable &p, int level) ;
		static void print (printable &p, int level, string file) ;
};



inline
void io::print (const string &s, int level )  {
	if (level <= debug)
		cout << s << flush;
}

inline
void io::println (const string &s, int level) {
	print (s+"\n",level);
}

inline
void io::print (const string &s, int level, string file) { 
	if (level <= debug) { 
		ofstream ofs (file.c_str());
		ofs << s;
		ofs.close();
	}
}

inline
void io::print (printable &p, int level) { 
	print (p.to_string(), level);
}

inline
void io::println (const printable &p, int level) {
	io::print (p.to_string()+"\n", level);
}

inline
void io::print (printable &p, int level, string file) { 
	print (p.to_string(), level, file);
}
#endif
